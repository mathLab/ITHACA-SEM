///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyViscousBurgers.cpp
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
// Description: Unsteady advection-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyViscousBurgers.h>
#include <StdRegions/StdQuadExp.h>

using namespace std;

namespace Nektar
{
    string UnsteadyViscousBurgers::className
    = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
                "UnsteadyViscousBurgers",
                UnsteadyViscousBurgers::create);
    
    UnsteadyViscousBurgers::UnsteadyViscousBurgers(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph),
          m_varCoeffLap(StdRegions::NullVarCoeffMap)
    {
        m_planeNumber = 0;
    }
    
    /**
     * @brief Initialisation object for the unsteady linear advection 
     * diffusion equation.
     */
    void UnsteadyViscousBurgers::v_InitObject()
    {
        AdvectionSystem::v_InitObject();
        
        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        m_session->LoadParameter("epsilon",    m_epsilon,  0.0);
        
        m_session->MatchSolverInfo(
            "SpectralVanishingViscosity", "True", m_useSpecVanVisc, false);
        m_session->MatchSolverInfo(
            "SpectralVanishingViscosity", "VarDiff", m_useSpecVanViscVarDiff, false);
        if(m_useSpecVanViscVarDiff)
        {
            m_useSpecVanVisc = true;
        }
        
        if(m_useSpecVanVisc)
        {
            m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
            m_session->LoadParameter("SVVDiffCoeff",m_sVVDiffCoeff,0.1);
        }        

        // Type of advection and diffusion classes to be used
        switch(m_projectionType)
        {
            // Discontinuous field 
            case MultiRegions::eDiscontinuous:
            {
                ASSERTL0(false,"Need to implement for DG");
                // Do not forwards transform initial condition
                m_homoInitialFwd = false;

                // Advection term
                string advName;
                string riemName; 
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::GetAdvectionFactory().
                    CreateInstance(advName, advName);
                m_advObject->SetFluxVector(&UnsteadyViscousBurgers::
                                           GetFluxVectorAdv, this);
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
                    CreateInstance(riemName, m_session);
                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject      (m_session, m_fields);
                
                // Diffusion term
                std::string diffName;
                m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                m_diffusion = SolverUtils::GetDiffusionFactory().
                    CreateInstance(diffName, diffName);
                m_diffusion->SetFluxVector(&UnsteadyViscousBurgers::
                                           GetFluxVectorDiff, this);
                m_diffusion->InitObject(m_session, m_fields);
                break;
            }
            // Continuous field 
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                // Advection term
                std::string advName;
                m_session->LoadSolverInfo("AdvectionType", advName, 
                                          "NonConservative");
                m_advObject = SolverUtils::GetAdvectionFactory().
                    CreateInstance(advName, advName);
                m_advObject->SetFluxVector(&UnsteadyViscousBurgers::
                                           GetFluxVectorAdv, this);
                
                if(m_useSpecVanViscVarDiff)
                {
                    Array<OneD, Array<OneD, NekDouble> > vel(m_fields.num_elements());
                    for(int i = 0; i < m_fields.num_elements(); ++i)
                    {
                        vel[i] = m_fields[i]->UpdatePhys();
                    }
                    SVVVarDiffCoeff(vel,m_varCoeffLap);
                }

                // In case of Galerkin explicit diffusion gives an error
                if (m_explicitDiffusion)
                {
                    ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
                }
                // In case of Galerkin implicit diffusion: do nothing
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }
        
        // Forcing terms
        m_forcing = SolverUtils::Forcing::Load(m_session, m_fields,
                                               m_fields.num_elements());
        
        m_ode.DefineImplicitSolve (&UnsteadyViscousBurgers::DoImplicitSolve, this);
        m_ode.DefineOdeRhs        (&UnsteadyViscousBurgers::DoOdeRhs,        this);
        
        if (m_projectionType == MultiRegions::eDiscontinuous &&
            m_explicitDiffusion == 1)
        {    
            m_ode.DefineProjection(&UnsteadyViscousBurgers::DoOdeProjection, this);
        }
    }
	
    /**
     * @brief Unsteady linear advection diffusion equation destructor.
     */
    UnsteadyViscousBurgers::~UnsteadyViscousBurgers()
    {
    }
    
    /**
     * @brief Get the normal velocity for the unsteady linear advection 
     * diffusion equation.
     */
    Array<OneD, NekDouble> &UnsteadyViscousBurgers::GetNormalVelocity(
                   Array<OneD, Array<OneD, NekDouble> >&inarray)
    {   
        // Number of trace (interface) points
        int i;
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);
        m_traceVn = Array<OneD, NekDouble>(nTracePts, 0.0);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (i = 0; i < inarray.num_elements(); ++i)
        {
            m_fields[0]->ExtractTracePhys(inarray[i], tmp);

            Vmath::Vvtvp(nTracePts,
                         m_traceNormals[i], 1,
                         tmp, 1,
                         m_traceVn, 1,
                         m_traceVn, 1);
        }
        
        return m_traceVn;
    }
    
    /**
     * @brief Compute the right-hand side for the unsteady linear advection 
     * diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyViscousBurgers::DoOdeRhs(
                                          const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                                          const NekDouble time)
    {  
        // Number of fields (variables of the problem)
        int nVariables = inarray.num_elements();
        
        // Number of solution points
        int nSolutionPts = GetNpoints();
        
        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nVariables);
        
        for (int i = 0; i < nVariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
        }
        
        // RHS computation using the new advection base class
        m_advObject->Advect(nVariables, m_fields, inarray,
                            inarray, outarray, time);
        
        // Negate the RHS
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Neg(nSolutionPts, outarray[i], 1);
        }
        
        // No explicit diffusion for CG
        if (m_projectionType == MultiRegions::eDiscontinuous)
        {
            m_diffusion->Diffuse(nVariables, m_fields, inarray, outarrayDiff);

            for (int i = 0; i < nVariables; ++i)
            {
                Vmath::Vadd(nSolutionPts, &outarray[i][0], 1, 
                            &outarrayDiff[i][0], 1, &outarray[i][0], 1);
            }
        }
        
        // Add forcing terms
        for (auto &x : m_forcing)
        {
            // set up non-linear terms
            x->Apply(m_fields, inarray, outarray, time);
        }
    }
    
    /**
     * @brief Compute the projection for the unsteady advection 
     * diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyViscousBurgers::DoOdeProjection(
                                                 const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                                 Array<OneD,       Array<OneD, NekDouble> > &outarray,
                                                 const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        SetBoundaryConditions(time);
        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();
                
                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
                break;
            }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
            {
                // Do nothing for the moment. 
            }
        default:
            {
                ASSERTL0(false, "Unknown projection scheme");
                break;
            }
        }
    }
    
    
#ifndef DEMO_IMPLICITSOLVER_JFNK_VIS
    /* @brief Compute the diffusion term implicitly. 
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     * @param lambda     Diffusion coefficient.
     */
    void UnsteadyViscousBurgers::DoImplicitSolve(
                                                 const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                                 Array<OneD,       Array<OneD, NekDouble> >&outarray,
                                                 const NekDouble time,
                                                 const NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int nq = m_fields[0]->GetNpoints();
		
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = 1.0/lambda/m_epsilon;
        
        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff/m_epsilon;
        }

        Array<OneD, Array< OneD, NekDouble> > F(nvariables);
        F[0] = Array<OneD, NekDouble> (nq*nvariables);
        
        for (int n = 1; n < nvariables; ++n)
        {
            F[n] = F[n-1] + nq;
        }
        
        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], 
                        inarray[i], 1, F[i], 1);
        }
        
        //Setting boundary conditions
        SetBoundaryConditions(time);
        
        if(m_useSpecVanViscVarDiff)
        {
            static int cnt = 0;

            if(cnt %10 == 0)
            {
                Array<OneD, Array<OneD, NekDouble> > vel(m_fields.num_elements());
                for(int i = 0; i < m_fields.num_elements(); ++i)
                {
                    m_fields[i]->ClearGlobalLinSysManager();
                    vel[i] = m_fields[i]->UpdatePhys();
                }
                SVVVarDiffCoeff(vel,m_varCoeffLap);
            }
            ++cnt;
        }
        for (int i = 0; i < nvariables; ++i)
        {
            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), 
                                   NullFlagList, factors, m_varCoeffLap);
            
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
        }
    }
#endif

#ifdef DEMO_IMPLICITSOLVER_JFNK_VIS

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
    
    void UnsteadyViscousBurgers::DoImplicitSolve(
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
        NekDouble tolrnc    = 1.0E-10;
        NekDouble tol2      = tolrnc*tolrnc;

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
        LibUtilities::CommSharedPtr v_Comm
             = m_fields[0].lock()->GetComm()->GetRowComm();

        NekLinSysIterative linsol(m_session,v_Comm);
        m_LinSysOprtors.DefineMatrixMultiply(&UnsteadyViscousBurgers::MatrixMultiply, this);
        m_LinSysOprtors.DefinePrecond(&UnsteadyViscousBurgers::preconditioner, this);
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
            
            dsolnorm = Vmath::Dot(ntotal,dsol_1D,dsol_1D);
            
            // the resfactor between L2norm of nonlinear risidual and dsol;
            if (0==iNonl)
            {
                resfactor = resnorm/dsolnorm;
            }


            if (dsolnorm*resfactor<tol2)
            {
                converged = true;
                break;
            }
            else
            {
                NonlinSysEvaluator(sol,NonlinSysRes);
            }
        }

        ASSERTL0((converged),"Nonlinear system solver not converge in UnsteadyViscousBurgers::DoImplicitSolve ")
        return;
    }
    
    

    void UnsteadyViscousBurgers::NonlinSysEvaluator(
                                                 const Array<OneD, Array<OneD, NekDouble> > &inarray,
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
        
    }

    
    void UnsteadyViscousBurgers::MatrixMultiply(
                                                 const Array<OneD, NekDouble> &inarray,
                                                 Array<OneD, NekDouble >&out)
    {
        MatrixMultiply_MatrixFree(inarray,out);
        return;
    }


    void UnsteadyViscousBurgers::MatrixMultiply_MatrixFree(
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
    }
    
    
    void UnsteadyViscousBurgers::preconditioner(
                                                 const Array<OneD, NekDouble> &inarray,
                                                 Array<OneD, NekDouble >&out)
    {
        int ntotal     = inarray.num_elements();
        Vmath::Vcopy(ntotal,inarray,1,out,1);
    }

#endif
    
    /**
     * @brief Return the flux vector for the advection part.
     * 
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyViscousBurgers::GetFluxVectorAdv(
                                                  const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                                                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {

        const int nq = m_fields[0]->GetNpoints();

        for (int i = 0; i < flux.num_elements(); ++i)
        {
            for (int j = 0; j < flux[0].num_elements(); ++j)
            {
                Vmath::Vmul(nq, physfield[i], 1, physfield[j], 1,
                            flux[i][j], 1);
            }
        }
    }

    /**
     * @brief Return the flux vector for the diffusion part.
     *      
     * @param i           Equation number.
     * @param j           Spatial direction.
     * @param physfield   Fields.
     * @param derivatives First order derivatives.
     * @param flux        Resulting flux.
     */
    void UnsteadyViscousBurgers::GetFluxVectorDiff(
        const int i, 
        const int j,
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &derivatives,
              Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        for (int k = 0; k < flux.num_elements(); ++k)
        {
            Vmath::Zero(GetNpoints(), flux[k], 1);
        }
        Vmath::Vcopy(GetNpoints(), physfield[i], 1, flux[j], 1);
    }
    
    void UnsteadyViscousBurgers::v_GenerateSummary(
            SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
        if(m_useSpecVanVisc)
        {
            stringstream ss;
            ss << "SVV (cut off = " << m_sVVCutoffRatio
               << ", coeff = "      << m_sVVDiffCoeff << ")";
            AddSummaryItem(s, "Smoothing", ss.str());
        }
    }
}
