///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyInviscidBurger.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Unsteady inviscid Burger solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/UnsteadyInviscidBurger.h>

using namespace std;

namespace Nektar
{
    string UnsteadyInviscidBurger::className
        = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
                "UnsteadyInviscidBurger",
                UnsteadyInviscidBurger::create,
                "Inviscid Burger equation");
    
    UnsteadyInviscidBurger::UnsteadyInviscidBurger(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
    }
    
    /**
     * @brief Initialisation object for the inviscid Burger equation.
     */
    void UnsteadyInviscidBurger::v_InitObject()
    {
        // Call to the initialisation object of UnsteadySystem
        AdvectionSystem::v_InitObject();
                
        // Define the normal velocity fields
        if (m_fields[0]->GetTrace())
        {
            m_traceVn  = Array<OneD, NekDouble>(GetTraceNpoints());
        }
        
        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Continuous field 
            case MultiRegions::eGalerkin:
            {
                string advName;
                m_session->LoadSolverInfo("AdvectionType", advName, "NonConservative");
                m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
                m_advObject->SetFluxVector   (&UnsteadyInviscidBurger::GetFluxVector, this);
                break;
            }
            // Discontinuous field 
            case MultiRegions::eDiscontinuous:
            {
                string advName;
                string riemName;
                
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
                m_advObject->SetFluxVector   (&UnsteadyInviscidBurger::GetFluxVector, this);
                
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(riemName, m_session);
                m_riemannSolver->SetScalar("Vn", &UnsteadyInviscidBurger::GetNormalVelocity, this);
                
                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject      (m_session, m_fields);
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }
                
        // If explicit it computes RHS and PROJECTION for the time integration
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&UnsteadyInviscidBurger::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyInviscidBurger::DoOdeProjection, this);
        }
        // Otherwise it gives an error because (no implicit integration)
        else
        {
            m_ode.DefineOdeRhs     (&UnsteadyInviscidBurger::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyInviscidBurger::DoOdeProjection, this);
            m_ode.DefineImplicitSolve(&UnsteadyInviscidBurger::DoImplicitSolve, this);
            // ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }
    
    /**
     * @brief Inviscid Burger equation destructor.
     */
    UnsteadyInviscidBurger::~UnsteadyInviscidBurger()
    {
    }
    
    /**
     * @brief Get the normal velocity for the inviscid Burger equation.
     */
    Array<OneD, NekDouble> &UnsteadyInviscidBurger::GetNormalVelocity()
    {
        // Counter variable
        int i;
        
        // Number of trace (interface) points
        int nTracePts       = GetTraceNpoints();
        
        // Number of solution points
        int nSolutionPts    = GetNpoints();
        
        // Number of fields (variables of the problem)
        int nVariables      = m_fields.num_elements();
        
        // Auxiliary variables to compute the normal velocity
        Array<OneD, NekDouble>               Fwd        (nTracePts);
        Array<OneD, Array<OneD, NekDouble> > physfield  (nVariables);
        
        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        // The TimeIntegration Class does not update the physical values of the 
        // solution. It is thus necessary to transform back the coefficient into
        // the physical space and save them in physfield to compute the normal 
        // advection velocity properly. However it remains a critical point. 
        for(i = 0; i < nVariables; ++i)
        {
            physfield[i]    = Array<OneD, NekDouble>(nSolutionPts);
            m_fields[i]->BwdTrans_IterPerExp(m_fields[i]->GetCoeffs(), 
                                             physfield[i]);
        }

        /// Extract the physical values at the trace space
        m_fields[0]->ExtractTracePhys(physfield[0], Fwd);
        
        /// Compute the normal velocity
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nTracePts, 
                         m_traceNormals[i], 1, 
                         Fwd, 1, 
                         m_traceVn, 1, 
                         m_traceVn, 1);
            
            Vmath::Smul(nTracePts, 0.5, m_traceVn, 1, m_traceVn, 1);
        }
        return m_traceVn;
    }
    
    /**
     * @brief Compute the right-hand side for the inviscid Burger equation.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyInviscidBurger::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;
        
        // Number of fields (variables of the problem)
        int nVariables      = inarray.num_elements();
        
        // Number of solution points
        int nSolutionPts    = GetNpoints();
        
        // !Useless variable for WeakDG and FR!
        Array<OneD, Array<OneD, NekDouble> >    advVel;    
        
        // RHS computation using the new advection base class
        m_advObject->Advect(nVariables, m_fields, advVel, inarray,
                            outarray, time);
        
        // Negate the RHS
        for (i = 0; i < nVariables; ++i)
        {
            Vmath::Neg(nSolutionPts, outarray[i], 1);
        }
    }
    
    /**
     * @brief Compute the projection for the inviscid Burger equation.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyInviscidBurger::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;
        
        // Number of variables of the problem
        int nVariables = inarray.num_elements();
        
        // Set the boundary conditions
        SetBoundaryConditions(time);
        
        // Switch on the projection type (Discontinuous or Continuous)
        switch(m_projectionType)
        {
            // Discontinuous projection
            case MultiRegions::eDiscontinuous:
            {
                // Number of quadrature points
                int nQuadraturePts = GetNpoints();
                
                // Just copy over array                
                for(i = 0; i < nVariables; ++i)
                {
                    Vmath::Vcopy(nQuadraturePts, inarray[i], 1, outarray[i], 1);
                }
                break;
            }
                
            // Continuous projection
            case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
                
                for(i = 0; i < nVariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }

#ifdef DEMO_IMPLICITSOLVER_JFNK

    /* @brief Compute the diffusion term implicitly. 
     *      Solve the whole system implicitly 
     *      incontract with only implicitly solve the 
     *      linear(diffusion) part in the origin DoImplicitSolve.
     * @param forc       The forcing term of the equation.
     * @param inoutarray input initial guess/output Calculated solution.
     * @param inrhs      input initial guess rhs.
     * @param time       Time.
     * @param lambda     Diffusion coefficient.
     */
    
    void UnsteadyInviscidBurger::DoImplicitSolve(
                                                 const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                                       Array<OneD,       Array<OneD, NekDouble> >&out,
                                                 const NekDouble time,
                                                 const NekDouble lambda)
    {
        m_TimeIntegtSol_n   = inarray;
        m_TimeIntegtSol_k   = out;
        m_TimeIntegLambda   = lambda;
        m_BndEvaluateTime   = time;

        bool l_verbose      = m_session->DefinesCmdLineArgument("verbose");

        const unsigned int MaxNonlinIte =   500;
        unsigned int nvariables  = inarray.num_elements();
        unsigned int npoints     = inarray[0].num_elements();
        unsigned int ntotal      = nvariables*npoints;

        bool converged;
        // unsigned int nGlobal,nDir;
        // NekDouble dsolnorm;
        // NekDouble resfactor = 1.0;
        NekDouble resnorm;
        NekDouble tolrnc    = 1.0E-8;
        NekDouble tol2      = tolrnc*tolrnc;
        NekDouble LinSysTol = 0.0;


        m_PrecMatVars = Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            m_PrecMatVars[i] =  Array<OneD, DNekBlkMatSharedPtr> (nvariables);
        }

        AllocatePrecondBlkDiag(m_PrecMatVars);

        Array<OneD, NekDouble> NonlinSysRes_1D(ntotal,0.0),sol_k_1D(ntotal,0.0),dsol_1D(ntotal,0.0);
        //Array<OneD,       Array<OneD, NekDouble> > sol_k;
        Array<OneD,       Array<OneD, NekDouble> > NonlinSysRes(nvariables),dsol(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            int offset = i*npoints;
            NonlinSysRes[i] =  NonlinSysRes_1D + offset;
            dsol[i]         =  dsol_1D + offset;
        }
        m_SysEquatResid_k = NonlinSysRes;

        //sol_k = out;
        for(int i = 0; i < nvariables; i++)
        {
            Vmath::Vcopy(npoints,inarray[i],1,m_TimeIntegtSol_k[i],1);
        }

        //ASSERTL0((1==nvariables),"only 1==nvariables")


        // TODO: 
        // v_Comm is based on the explist used may be different(m_tracemap or m_locToGloMap) for diffrent
        // should generate GlobalLinSys first and get the v_Comm from GlobalLinSys. here just give it a m_Comm no parallel support yet!!
        //const std::weak_ptr<ExpList> 
        LibUtilities::CommSharedPtr v_Comm
             = m_fields[0]->GetComm()->GetRowComm();

        //LibUtilities::CommSharedPtr                 v_Comm;
        NekLinSysIterative linsol(m_session,v_Comm);
        m_LinSysOprtors.DefineMatrixMultiply(&UnsteadyInviscidBurger::MatrixMultiply, this);
        m_LinSysOprtors.DefinePrecond(&UnsteadyInviscidBurger::preconditioner_BlkDiag, this);
        linsol.setLinSysOperators(m_LinSysOprtors);
        // NonlinSysEvaluator(m_TimeIntegtSol_k,m_SysEquatResid_k);
        NonlinSysEvaluator(m_TimeIntegtSol_k,m_SysEquatResid_k);
        DebugNumCalJac(m_PrecMatVars);
        // Cout2DArrayBlkMat(m_PrecMatVars);
        ElmtVarInvMtrx(m_PrecMatVars);
        // Cout2DArrayBlkMat(m_PrecMatVars);
        converged = false;
        for (int k = 0; k < MaxNonlinIte; k++)
        {
            
            NonlinSysEvaluator(m_TimeIntegtSol_k,m_SysEquatResid_k);
            // GetpreconditionerNSBlkDiag(m_TimeIntegtSol_k,m_PrecMatVars);

            // DebugNumCalElmtJac(0);

            // NonlinSysRes_1D and m_SysEquatResid_k share the same storage
            resnorm = Vmath::Dot(ntotal,NonlinSysRes_1D,NonlinSysRes_1D);
            if (resnorm<tol2)
            {
                /// TODO: m_root
                if(l_verbose)
                {
                    cout <<"Newton iteration converged with residual:   " << sqrt(resnorm)<<"  using   "<<k<<"   iterations"<<endl;
                }
                converged = true;
                break;
            }
            

            //TODO: currently  NonlinSysRes is 2D array and SolveLinearSystem needs 1D array
            LinSysTol = 0.01*sqrt(resnorm);
            linsol.SolveLinearSystem(ntotal,NonlinSysRes_1D,dsol_1D,0,LinSysTol);

            for(int i = 0; i < nvariables; i++)
            {
                Vmath::Vsub(npoints,m_TimeIntegtSol_k[i],1,dsol[i],1,m_TimeIntegtSol_k[i],1);
            }

        }

        /// TODO: disconnect these from other arrays to avoid memory cannot release.
        // m_TimeIntegtSol_k   =   nullptr;
        // m_TimeIntegtSol_n   =   nullptr;
        // m_SysEquatResid_k   =   nullptr;
        ASSERTL0((converged),"Nonlinear system solver not converge in CompressibleFlowSystem::DoImplicitSolve ")
        return;
    }

    

    void UnsteadyInviscidBurger::NonlinSysEvaluator(
                                                       Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                       Array<OneD, Array<OneD, NekDouble> > &out)
    {
        Array<OneD, Array<OneD, NekDouble> > sol_n;
        sol_n                  = m_TimeIntegtSol_n;
        //inforc = m_TimeIntegForce;
        unsigned int nvariable = inarray.num_elements();
        unsigned int npoints   = inarray[0].num_elements();
        
        
        DoOdeProjection(inarray,inarray,m_BndEvaluateTime);
        DoOdeRhs(inarray,out,m_BndEvaluateTime);
        for (int i = 0; i < nvariable; i++)
        {
            Vmath::Svtvp(npoints,m_TimeIntegLambda,out[i],1,sol_n[i],1,out[i],1);
            Vmath::Vsub(npoints,inarray[i],1,out[i],1,out[i],1);
        }
        return;
    }

    
    void UnsteadyInviscidBurger::MatrixMultiply(
                                                 const Array<OneD, NekDouble> &inarray,
                                                 Array<OneD, NekDouble >&out)
    {
        MatrixMultiply_MatrixFree(inarray,out);
        return;
    }


    void UnsteadyInviscidBurger::MatrixMultiply_MatrixFree(
                                                 const  Array<OneD, NekDouble> &inarray,
                                                        Array<OneD, NekDouble >&out)
    {
        NekDouble eps = 1.0E-6;
        NekDouble oeps = 1.0/eps;
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        Array<OneD, NekDouble > tmp;
        Array<OneD,       Array<OneD, NekDouble> > solplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resplus(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            solplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
        }

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,eps,tmp,1,m_TimeIntegtSol_k[i],1,solplus[i],1);
        }
        NonlinSysEvaluator(solplus,resplus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = out + i*npoints;
            // cout << "resplus[i][0]           =" << resplus[i][0]           <<endl;
            // cout << "m_SysEquatResid_k[i][0] =" << m_SysEquatResid_k[i][0] <<endl;
            // cout << "tmp[0]                  =" << tmp[0]                  <<endl;
            Vmath::Vsub(npoints,&resplus[i][0],1,&m_SysEquatResid_k[i][0],1,&tmp[0],1);
            Vmath::Smul(npoints, oeps ,&tmp[0],1,&tmp[0],1);
        }
        return;
    }
    
    
    void UnsteadyInviscidBurger::preconditioner(
                                                 const Array<OneD, NekDouble> &inarray,
                                                 Array<OneD, NekDouble >&out)
    {
        int ntotal     = inarray.num_elements();
        Vmath::Vcopy(ntotal,inarray,1,out,1);
        return;
    }


    void UnsteadyInviscidBurger::DebugNumCalJac(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
            std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int nElmts    = (*pexp).size();
        int nvariables= gmtxarray.num_elements();

        Array<OneD, Array<OneD, DNekMatSharedPtr> >  ElmtPrecMatVars;
        for(int i = 0; i < nElmts; i++)
        {
            DebugNumCalElmtJac(ElmtPrecMatVars,i);
            for(int m = 0; m < nvariables; m++)
            {
                for(int n = 0; n < nvariables; n++)
                {
                    gmtxarray[m][n]->SetBlock(i,i,ElmtPrecMatVars[m][n]);
                }
            }
        }
    }

    
    void UnsteadyInviscidBurger::DebugNumCalElmtJac(Array<OneD, Array<OneD, DNekMatSharedPtr> > &ElmtPrecMatVars ,const int nelmt)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
            std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int nElmtPnt    = (*pexp)[nelmt]->GetTotPoints();
        int nElmtOffset = explist->GetPhys_Offset(nelmt);
        
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        unsigned int ntotpnt    = nvariables*npoints;
        Array<OneD, NekDouble > tmpinn_1d(ntotpnt,0.0);
        Array<OneD, NekDouble > tmpout_1d(ntotpnt,0.0);
        Array<OneD,       Array<OneD, NekDouble> > tmpinn(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > tmpout(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            int noffset = i*npoints;
            tmpinn[i] = tmpinn_1d+noffset;
            tmpout[i] = tmpout_1d+noffset;
        }

        DNekMatSharedPtr    tmpStdMat;

        ElmtPrecMatVars = Array<OneD, Array<OneD, DNekMatSharedPtr> >  (nvariables);

        for(int i = 0; i < nvariables; i++)
        {
            ElmtPrecMatVars[i] =  Array<OneD, DNekMatSharedPtr> (nvariables);
            for(int j = 0; j < nvariables; j++)
            {
                ElmtPrecMatVars[i][j] =  MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nElmtPnt, nElmtPnt, 0.0);
            }
        }

        for (int i = 0; i < nvariables; i++)
        {
            for (int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                tmpinn[i][nElmtOffset+npnt] = 1.0;
                MatrixMultiply(tmpinn_1d,tmpout_1d);
                for (int j = 0; j < nvariables; j++)
                {
                    for (int npntf = 0; npntf < nElmtPnt; npntf++)
                    {
                        tmpStdMat = ElmtPrecMatVars[j][i];
                        (*tmpStdMat)(npntf,npnt) = tmpout[j][nElmtOffset+npntf];
                    }
                }
                tmpinn[i][nElmtOffset+npnt] = 0.0;
            }
        }
        // Cout2DArrayBlkMat(ElmtPrecMatVars);
        return ElmtPrecMatVars;
    }

    void UnsteadyInviscidBurger::preconditioner_BlkDiag(
                                                 const Array<OneD, NekDouble> &inarray,
                                                 Array<OneD, NekDouble >&outarray)
    {
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        Array<OneD, NekVector<NekDouble> > tmparray(nvariables);


        PointerWrapper pwrapp = eWrapper;
        if(inarray.get() == outarray.get())
        {
            pwrapp = eCopy;
        }

        for(int m = 0; m < nvariables; m++)
        {
            int moffset = m*npoints;
            tmparray[m] =  NekVector<NekDouble> (npoints,inarray+moffset,pwrapp);
        }
        
        Vmath::Fill(outarray.num_elements(),0.0,outarray,1);
        for(int m = 0; m < nvariables; m++)
        {
            int moffset = m*npoints;
            NekVector<NekDouble> out(npoints,outarray+moffset,eWrapper);
            for(int n = 0; n < nvariables; n++)
            {
                out += (*m_PrecMatVars[m][n])*tmparray[m];
            }
        }
        
    }
    void UnsteadyInviscidBurger::AllocatePrecondBlkDiag(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {

        int nvars = m_fields.num_elements();
        int nelmts  = m_fields[0]->GetNumElmts();
        int nrowsVars,ncolsVars;
        int nelmtcoef,nelmtpnts;
        DNekMatSharedPtr loc_matNvar;
        Array<OneD, unsigned int > nelmtmatdim(nelmts);
        for(int i = 0; i < nelmts; i++)
        {
            nelmtcoef   =   m_fields[0]->GetExp(i)->GetNcoeffs();
            nelmtpnts   =   m_fields[0]->GetExp(i)->GetTotPoints();
            ASSERTL0(nelmtpnts>=nelmtcoef,"in AllocatePrecondBlkDiag nelmtpnts>=nelmtcoef must hold");
            nelmtmatdim[i]  =   nelmtpnts;
        }

        for(int i = 0; i < nvars; i++)
        {
            for(int j = 0; j < nvars; j++)
            {
                gmtxarray[i][j] = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(nelmtmatdim, nelmtmatdim, eDIAGONAL);
                for(int nelm = 0; nelm < nelmts; ++nelm)
                {
                    nrowsVars = nelmtmatdim[nelm];
                    ncolsVars = nrowsVars;
                    
                    loc_matNvar = MemoryManager<DNekMat>::AllocateSharedPtr(nrowsVars,ncolsVars,0.0);
                    gmtxarray[i][j]->SetBlock(nelm,nelm,loc_matNvar);
                }
                
            }
        }
    }

    void UnsteadyInviscidBurger::ElmtVarInvMtrx(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
            std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int ntotElmt            = (*pexp).size();
        int nConvectiveFields = m_fields.num_elements();

        DNekMatSharedPtr        ElmtMat;

        Array<OneD, NekDouble> coefarray(explist->GetNcoeffs(),0.0);
        Array<OneD, NekDouble> physarray(explist->GetTotPoints(),0.0);

        Array<OneD, int > elmtpnts(ntotElmt);
        Array<OneD, int > elmtcoef(ntotElmt);
        int nElmtPnt   =    (*pexp)[0]->GetTotPoints();
        int nElmtPnt0  =    nElmtPnt;

        int nElmtPntVAr = nElmtPnt0*nConvectiveFields;
        DNekMatSharedPtr    tmpGmtx = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nElmtPntVAr, nElmtPntVAr,0.0);
        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();
            
            if (nElmtPnt0!=nElmtPnt) 
            {
                nElmtPnt0 = nElmtPnt;
                nElmtPntVAr = nElmtPnt0*nConvectiveFields;
                tmpGmtx = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nElmtPntVAr, nElmtPntVAr,0.0);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    ElmtMat =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                    for(int nrw = 0; nrw < nElmtPnt; nrw++)
                    {
                        int nrwvar = m*nElmtPnt+nrw;
                        for(int ncl = 0; ncl < nElmtPnt; ncl++)
                        {
                            int nclvar = n*nElmtPnt+ncl;
                            (*tmpGmtx)(nrwvar,nclvar)=(*ElmtMat)(nrw,ncl);
                        }
                    }
                }
            }

            tmpGmtx->Invert();

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    ElmtMat =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                    for(int nrw = 0; nrw < nElmtPnt; nrw++)
                    {
                        int nrwvar = m*nElmtPnt+nrw;
                        for(int ncl = 0; ncl < nElmtPnt; ncl++)
                        {
                            int nclvar = n*nElmtPnt+ncl;
                            (*ElmtMat)(nrw,ncl) =   (*tmpGmtx)(nrwvar,nclvar);
                        }
                    }
                }
            }
        }
        return;
    }

    void UnsteadyInviscidBurger::GetpreconditionerNSBlkDiag(
                                            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {
        // DoOdeProjection(inarray,inarray,m_BndEvaluateTime);

        Fill2DArrayOfBlkDiagonalMat(gmtxarray,0.0);
        // Cout2DArrayBlkMat(gmtxarray);
        // AddMatNSBlkDiag_volume(inarray,gmtxarray);
        // Cout2DArrayBlkMat(gmtxarray);

        AddMatNSBlkDiag_boundary(inarray,gmtxarray);
        // Cout2DArrayBlkMat(gmtxarray);
        MultiplyElmtBwdInvMassFwd(gmtxarray,m_TimeIntegLambda);
        // Cout2DArrayBlkMat(gmtxarray);
        ElmtVarInvMtrx(gmtxarray);
    }

    void UnsteadyInviscidBurger::Fill2DArrayOfBlkDiagonalMat(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const NekDouble valu)
    {
        
        int n1d = gmtxarray.num_elements();
        int n2d = gmtxarray[0].num_elements();

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        DNekMatSharedPtr    loc_matNvar;

        for(int n1 = 0; n1 < n1d; ++n1)
        {
            for(int n2 = 0; n2 < n2d; ++n2)
            {
                
                gmtxarray[n1][n2]->GetBlockSizes(rowSizes,colSizes);
                int nelmts  = rowSizes.num_elements();

                for(int i = 0; i < nelmts; ++i)
                {
                    loc_matNvar =   gmtxarray[n1][n2]->GetBlock(i,i);

                    int nrows = loc_matNvar->GetRows();
                    int ncols = loc_matNvar->GetColumns();

                    for(int j = 0; j < nrows; j++)
                    {
                        for(int k = 0; k < ncols; k++)
                        {
                            (*loc_matNvar)(j,k)=valu;
                        }
                    }
                }
            }
        }

    }


    void UnsteadyInviscidBurger::MultiplyElmtBwdInvMassFwd(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const NekDouble dtlamda)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
            std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int ntotElmt            = (*pexp).size();
        int nElmtPnt,nElmtCoef;
        int nConvectiveFields = m_fields.num_elements();

        NekDouble Negdtlamda    =   -dtlamda;
        DNekMatSharedPtr        tmpGmtx,ElmtMat;
        DNekMatSharedPtr        Fwdmat;

        Array<OneD, NekDouble> coefarray(explist->GetNcoeffs(),0.0);
        Array<OneD, NekDouble> physarray(explist->GetTotPoints(),0.0);

        Array<OneD, int > elmtpnts(ntotElmt);
        Array<OneD, int > elmtcoef(ntotElmt);
        int nElmtCoef0 =    (*pexp)[0]->GetNcoeffs();
        int nElmtPnt0  =    (*pexp)[0]->GetTotPoints();
        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
            nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();
            ASSERTL0(nElmtCoef==nElmtCoef0,"nElmtCoef==nElmtCoef0");
            ASSERTL0(nElmtPnt==nElmtPnt0,"nElmtPnt==nElmtPnt0");
            elmtpnts[nelmt]     =   nElmtPnt;
            elmtcoef[nelmt]     =   nElmtCoef;
        }

        for(int m = 0; m < nConvectiveFields; m++)
        {
            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                {
                    StdRegions::StdMatrixKey  key(StdRegions::eFwdTrans,
                                          (*pexp)[nelmt]->DetShapeType(), *(*pexp)[nelmt]);
                    Fwdmat = (*pexp)[nelmt]->GetStdMatrix(key);

                    nElmtPnt            = elmtpnts[nelmt];
                    nElmtCoef           = elmtcoef[nelmt];
                    tmpGmtx =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);

                    DNekMatSharedPtr coefmat = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(1, nElmtCoef,0.0);
                    DNekMatSharedPtr pntsmat = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(1, nElmtPnt,0.0);

                    for(int nrw = 0; nrw < nElmtCoef; nrw++)
                    {
                        for(int ncl = 0; ncl < nElmtCoef; ncl++)
                        {
                            (*coefmat)(0,ncl)   =   (*tmpGmtx)(nrw,ncl);
                        }
                        (*pntsmat) = (*coefmat)*(*Fwdmat);
                        for(int ncl = 0; ncl < nElmtPnt; ncl++)
                        {
                            (*tmpGmtx)(nrw,ncl) =   (*pntsmat)(0,ncl);
                        }
                    }
                }
            }
        }

        for(int m = 0; m < nConvectiveFields; m++)
        {
            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int ncl = 0; ncl < nElmtPnt0; ncl++)
                {
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtPnt            = elmtpnts[nelmt];
                        nElmtCoef           = elmtcoef[nelmt];
                        tmpGmtx =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);

                        int n_offset    =   explist->GetCoeff_Offset(nelmt);
                        
                        for(int nrw = 0; nrw < nElmtCoef; nrw++)
                        {
                            coefarray[n_offset+nrw]   =   (*tmpGmtx)(nrw,ncl);
                        }
                    }
                    explist->MultiplyByElmtInvMass(coefarray, coefarray);
                    explist->BwdTrans             (coefarray, physarray);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtPnt            = elmtpnts[nelmt];
                        tmpGmtx =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        int n_offset    =   explist->GetPhys_Offset(nelmt);
                        for(int nrw = 0; nrw < nElmtPnt; nrw++)
                        {
                            (*tmpGmtx)(nrw,ncl)   = Negdtlamda*physarray[n_offset+nrw];
                        }
                    }
                }
            }
        }

        for(int m = 0; m < nConvectiveFields; m++)
        {
            for(int ncl = 0; ncl < nElmtPnt0; ncl++)
            {
                for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                {
                    tmpGmtx =   gmtxarray[m][m]->GetBlock(nelmt,nelmt);
                    (*tmpGmtx)(ncl,ncl)   += 1.0;
                }
            }
        }
        return;
    }

    void UnsteadyInviscidBurger::AddMatNSBlkDiag_volume(
                                        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {
        int nSpaceDim = m_graph->GetSpaceDimension();
        int nvariable = inarray.num_elements();

        Array<OneD, Array<OneD, DNekMatSharedPtr> > ElmtJac;
        
        for(int nDirctn = 0; nDirctn < nSpaceDim; nDirctn++)
        {
            ElmtJac = GetFluxVectorJacDirctn(nDirctn,inarray);
            m_advObject->AddVolumJac2Mat(nvariable,m_fields,ElmtJac,nDirctn,gmtxarray);
        }
    }


    void UnsteadyInviscidBurger::AddMatNSBlkDiag_boundary(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {
        int nvariables = inarray.num_elements();
        Array<OneD, DNekBlkMatSharedPtr > TraceJac;
        TraceJac    =   GetTraceJac(inarray);
        m_advObject->AddTraceJac2Mat(nvariables,m_fields, TraceJac,gmtxarray);
    }

    Array<OneD, DNekBlkMatSharedPtr> UnsteadyInviscidBurger::GetTraceJac(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray)
    {
        int nvariables = inarray.num_elements();
        // int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nvariables);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nvariables);

        Array<OneD, unsigned int> n_blks(nTracePts);
        for(int i=0;i<nTracePts;i++)
        {
            n_blks[i]    = nvariables;
        }
        DNekBlkMatSharedPtr FJac = MemoryManager<DNekBlkMat>
            ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);
        DNekBlkMatSharedPtr BJac = MemoryManager<DNekBlkMat>
            ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);

        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for(int i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }

        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        m_advObject->CalcTraceJac(nvariables, m_fields, advVel, inarray,
                            Fwd, Bwd, FJac, BJac);

        Array<OneD, DNekBlkMatSharedPtr>    TraceJac(2);
        TraceJac[0] = FJac;
        TraceJac[1] = BJac;

        return TraceJac;
    }

    Array<OneD, Array<OneD, DNekMatSharedPtr> > UnsteadyInviscidBurger::GetFluxVectorJacDirctn(const int nDirctn,
                                                        const Array<OneD, const Array<OneD, NekDouble> >&inarray)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();
        LocalRegions::ExpansionSharedPtr pexp = (*expvect)[0];
        // int nElmtPnt            = pexp->GetTotPoints();
        // int nElmtCoef           = pexp->GetNcoeffs();
        // int nSpaceDim           = m_graph->GetSpaceDimension();  
        Array<OneD, NekDouble> normals;
        Array<OneD, Array<OneD, NekDouble> > normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        normals =   normal3D[nDirctn];

        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);

        Array<OneD, Array<OneD, DNekMatSharedPtr> > ElmtJac(ntotElmt);

        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);


        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            // nElmtCoef           = (*expvect)[nelmt]->GetNcoeffs();
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
    
            for(int j = 0; j < nConvectiveFields; j++)
            {   
                locVars[j] = inarray[j]+GetPhys_Offset(nelmt);
            }

            ElmtJac[nelmt] =   Array<OneD, DNekMatSharedPtr>(nElmtPnt);
            for(int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                ElmtJac[nelmt][npnt] = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nConvectiveFields, nConvectiveFields);
                for(int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }

                GetFluxVectorJacPoint(pointVar,normals,ElmtJac[nelmt][npnt]);
            }
        }
        return ElmtJac;
    }




    void UnsteadyInviscidBurger::GetFluxVectorJacPoint(
            const Array<OneD, NekDouble>                &conservVar, 
            const Array<OneD, NekDouble>                &normals, 
                 DNekMatSharedPtr                       &fluxJac)
    {
        int nvariables      = conservVar.num_elements();
        int expDim          = m_spacedim;

        if (nvariables > expDim+2)
        {
            ASSERTL0(false,"nvariables > expDim+2 case not coded")
        }

        DNekMatSharedPtr PointFJac3D = MemoryManager<DNekMat>
            ::AllocateSharedPtr(nvariables, nvariables);
        
        Array<OneD, NekDouble> PointFwd(nvariables,0.0);


        int nj=0;
        int nk=0;
        for(int j=0; j< nvariables; j++)
        {
            nj = j;
            PointFwd[nj] = conservVar[j];
        }
        NekDouble fsw,efix_StegerWarming;
        efix_StegerWarming = 0.0;

        fsw = 0.0; // exact flux Jacobian if fsw=0.0
        PointFluxJacobian_pn(PointFwd,normals,PointFJac3D,efix_StegerWarming,fsw);

        for(int j=0; j< nvariables; j++)
        {
            nj = j;
            for(int k=0; k< nvariables; k++)
            {
                nk = k;
                (*fluxJac)(j,k) = (*PointFJac3D)(nj,nk); 
            }
        }
        
    }

    // Currently duplacate in compressibleFlowSys
    // if fsw=+-1 calculate the steger-Warming flux vector splitting flux Jacobian
    // if fsw=0   calculate the Jacobian of the exact flux 
    // efix is the numerical flux entropy fix parameter
    void UnsteadyInviscidBurger::PointFluxJacobian_pn(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw)
    {
        for(int i = 0; i < Fwd.num_elements(); i++)
        {
            (*FJac)(i,i)   =   Fwd[i];
        }
    }

    



#endif


    /**
     * @brief Return the flux vector for the inviscid Burger equation.
     * 
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyInviscidBurger::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        const int nq = GetNpoints();

        for (int i = 0; i < flux.num_elements(); ++i)
        {
            for (int j = 0; j < flux[0].num_elements(); ++j)
            {
                Vmath::Vmul(nq, physfield[i], 1, physfield[i], 1, 
                            flux[i][j], 1);
                Vmath::Smul(nq, 0.5, flux[i][j], 1, flux[i][j], 1);
            }
        }
    }
}

