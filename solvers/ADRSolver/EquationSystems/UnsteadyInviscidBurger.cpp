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
                
        m_ode.DefineImplicitSolve(&UnsteadyInviscidBurger::DoImplicitSolve, this);
        // If explicit it computes RHS and PROJECTION for the time integration
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&UnsteadyInviscidBurger::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyInviscidBurger::DoOdeProjection, this);
        }
        // Otherwise it gives an error because (no implicit integration)
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
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
                                                 const Array<OneD, const Array<OneD, NekDouble> >&forc,
                                                       Array<OneD,       Array<OneD, NekDouble> >&sol,
                                                 const NekDouble time,
                                                 const NekDouble lambda)
    {
        m_TimeIntegSoltn  = sol;
        m_TimeIntegForce  = forc;
        m_TimeIntegLambda = lambda;
        m_BndEvaluateTime = time;

        

        const unsigned int MaxNonlinIte =   500;
        unsigned int nvariables  = forc.num_elements();
        unsigned int npoints     = forc[0].num_elements();
        unsigned int ntotal      = nvariables*npoints;

        bool converged;
        unsigned int nGlobal,nDir;
        NekDouble resnorm,dsolnorm;
        NekDouble resfactor = 1.0;
        NekDouble tolrnc    = 1.0E-6;
        NekDouble tol2      = tolrnc*tolrnc;


        for (int i = 0; i < npoints; i++)
        {
            cout <<"forc["<<i<<"]= "<<forc[0][i]<<endl;
        }
        for (int i = 0; i < npoints; i++)
        {
            cout <<"sol["<<i<<"]= "<<sol[0][i]<<endl;
        }

        Array<OneD, NekDouble> NonlinSysRes_1D(ntotal,0.0),dsol_1D(ntotal,0.0);
        Array<OneD,       Array<OneD, NekDouble> > NonlinSysRes(nvariables),dsol(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            int offset = i*npoints;
            NonlinSysRes[i] =  NonlinSysRes_1D + offset;
            dsol[i]         =  dsol_1D + offset;
        }
        m_SysEquatResidu = NonlinSysRes;

        ASSERTL0((1==nvariables),"only 1==nvariables")

        NonlinSysEvaluator(sol,NonlinSysRes);

        resnorm = Vmath::Dot(ntotal,NonlinSysRes_1D,NonlinSysRes_1D);
            

        if (resnorm<tol2)
        {
            return;
        }

        // TODO: 
        // v_Comm is based on the explist used may be different(m_tracemap or m_locToGloMap) for diffrent
        // should generate GlobalLinSys first and get the v_Comm from GlobalLinSys. here just give it a m_Comm no parallel support yet!!
        //const std::weak_ptr<ExpList>       explist= *m_fields[0];
        LibUtilities::CommSharedPtr v_Comm
             = m_fields[0]->GetComm()->GetRowComm();

        //LibUtilities::CommSharedPtr                 v_Comm;
        NekLinSysIterative linsol(m_session,v_Comm);
        m_LinSysOprtors.DefineMatrixMultiply(&UnsteadyInviscidBurger::MatrixMultiply, this);
        m_LinSysOprtors.DefinePrecond(&UnsteadyInviscidBurger::preconditioner, this);
        linsol.setLinSysOperators(m_LinSysOprtors);

        converged = false;
        for (int iNonl = 0; iNonl < MaxNonlinIte; iNonl++)
        {
            //TODO: currently  NonlinSysRes is 2D array and SolveLinearSystem needs 1D array
            linsol.SolveLinearSystem(ntotal,NonlinSysRes_1D,dsol_1D,0);

            for(int i = 0; i < nvariables; i++)
            {
                Vmath::Vadd(npoints,dsol[i],1,sol[i],1,sol[i],1);
            }
            
            // dsolnorm = Vmath::Dot(ntotal,dsol_1D,dsol_1D);
            
            // // the resfactor between L2norm of nonlinear risidual and dsol;
            // if (0==iNonl)
            // {
            //     resfactor = resnorm/dsolnorm;
            // }
            // resfactor = 1.0;
            NonlinSysEvaluator(sol,NonlinSysRes);
            resnorm = Vmath::Dot(ntotal,dsol_1D,dsol_1D);
            cout << "   iNonl = "<<iNonl<<" resnorm = "<<resnorm<<"  resfactor = "<<resfactor<<endl;
            if (resnorm<tol2)
            {
                converged = true;
                break;
            }

        }

        ASSERTL0((converged),"Nonlinear system solver not converge in UnsteadyInviscidBurger::DoImplicitSolve ")
        return;
    }
    
    

    void UnsteadyInviscidBurger::NonlinSysEvaluator(
                                                       Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                       Array<OneD, Array<OneD, NekDouble> > &out)
    {
        Array<OneD, Array<OneD, NekDouble> > inforc;
        inforc = m_TimeIntegForce;
        unsigned int nvariable = inforc.num_elements();
        unsigned int npoints   = inforc[0].num_elements();
        

        DoOdeProjection(inarray,inarray,m_BndEvaluateTime);
        DoOdeRhs(inarray,out,m_BndEvaluateTime);
        for (int i = 0; i < nvariable; i++)
        {
            Vmath::Svtvp(npoints,-m_TimeIntegLambda,out[i],1,inarray[i],1,out[i],1);
            Vmath::Vsub(npoints,out[i],1,inforc[i],1,out[i],1);
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
        unsigned int nvariables = m_TimeIntegSoltn.num_elements();
        unsigned int npoints    = m_TimeIntegSoltn[0].num_elements();
        Array<OneD, NekDouble > tmp;
        Array<OneD,       Array<OneD, NekDouble> > solplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resplus(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            solplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
        }

        tmp = inarray;
        for (int i = 0; i < nvariables; i++)
        {
            tmp = tmp + i*npoints;
            Vmath::Svtvp(npoints,eps,tmp,1,m_TimeIntegSoltn[i],1,solplus[i],1);
        }
        NonlinSysEvaluator(solplus,resplus);

        tmp = out;
        for (int i = 0; i < nvariables; i++)
        {
            tmp = tmp + i*npoints;
            Vmath::Vvpts(npoints,&solplus[i][0],1,&m_TimeIntegSoltn[i][0],1,oeps,&tmp[0],1);
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

