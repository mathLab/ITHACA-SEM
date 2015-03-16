///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvectionDiffusion.cpp
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

#include <ADRSolver/EquationSystems/UnsteadyAdvectionDiffusion.h>

namespace Nektar
{
    string UnsteadyAdvectionDiffusion::className
        = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
                "UnsteadyAdvectionDiffusion",
                UnsteadyAdvectionDiffusion::create);
    
    UnsteadyAdvectionDiffusion::UnsteadyAdvectionDiffusion(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession),
          AdvectionSystem(pSession)
    {
        m_planeNumber = 0;
    }
    
    /**
     * @brief Initialisation object for the unsteady linear advection 
     * diffusion equation.
     */
    void UnsteadyAdvectionDiffusion::v_InitObject()
    {
        AdvectionSystem::v_InitObject();
        
        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        m_session->LoadParameter("epsilon",    m_epsilon,  0.0);
        
        // Define Velocity fields
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");
        vel.resize(m_spacedim);
        
        EvaluateFunction(vel, m_velocity, "AdvectionVelocity");
        
        m_session->MatchSolverInfo(
            "SpectralVanishingViscosity", "True", m_useSpecVanVisc, false);
        
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
                // Do not forwards transform initial condition
                m_homoInitialFwd = false;

                // Advection term
                string advName;
                string riemName; 
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::GetAdvectionFactory().
                    CreateInstance(advName, advName);
                m_advObject->SetFluxVector(&UnsteadyAdvectionDiffusion::
                                           GetFluxVectorAdv, this);
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
                    CreateInstance(riemName);
                m_riemannSolver->SetScalar("Vn", &UnsteadyAdvectionDiffusion::
                                           GetNormalVelocity, this);
                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject      (m_session, m_fields);
                
                // Diffusion term
                std::string diffName;
                m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                m_diffusion = SolverUtils::GetDiffusionFactory().
                    CreateInstance(diffName, diffName);
                m_diffusion->SetFluxVector(&UnsteadyAdvectionDiffusion::
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
                m_advObject->SetFluxVector(&UnsteadyAdvectionDiffusion::
                                           GetFluxVectorAdv, this);
                
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
        
        m_ode.DefineImplicitSolve (&UnsteadyAdvectionDiffusion::DoImplicitSolve, this);
        m_ode.DefineOdeRhs        (&UnsteadyAdvectionDiffusion::DoOdeRhs,        this);
        
        if (m_projectionType == MultiRegions::eDiscontinuous &&
            m_explicitDiffusion == 1)
        {    
            m_ode.DefineProjection(&UnsteadyAdvectionDiffusion::DoOdeProjection, this);
        }
    }
	
    /**
     * @brief Unsteady linear advection diffusion equation destructor.
     */
    UnsteadyAdvectionDiffusion::~UnsteadyAdvectionDiffusion()
    {
    }
    
    /**
     * @brief Get the normal velocity for the unsteady linear advection 
     * diffusion equation.
     */
    Array<OneD, NekDouble> &UnsteadyAdvectionDiffusion::GetNormalVelocity()
    {   
        // Number of trace (interface) points
        int i;
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);
        m_traceVn = Array<OneD, NekDouble>(nTracePts, 0.0);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (i = 0; i < m_velocity.num_elements(); ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], tmp);

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
    void UnsteadyAdvectionDiffusion::DoOdeRhs(
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
        m_advObject->Advect(nVariables, m_fields, m_velocity,
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
        
    }
    
    /**
     * @brief Compute the projection for the unsteady advection 
     * diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvectionDiffusion::DoOdeProjection(
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
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
                
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
                break;
            }
            default:
            {
                ASSERTL0(false, "Unknown projection scheme");
                break;
            }
        }
    }
    
    
    /* @brief Compute the diffusion term implicitly. 
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     * @param lambda     Diffusion coefficient.
     */
    void UnsteadyAdvectionDiffusion::DoImplicitSolve(
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
        
        for (int i = 0; i < nvariables; ++i)
        {
            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), 
                                   NullFlagList, factors);
            
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
        }
    }
    
    /**
     * @brief Return the flux vector for the advection part.
     * 
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyAdvectionDiffusion::GetFluxVectorAdv(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].num_elements() == m_velocity.num_elements(),
                 "Dimension of flux array and velocity array do not match");

        const int nq = m_fields[0]->GetNpoints();

        for (int i = 0; i < flux.num_elements(); ++i)
        {
            for (int j = 0; j < flux[0].num_elements(); ++j)
            {
                Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1,
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
    void UnsteadyAdvectionDiffusion::GetFluxVectorDiff(
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
    
    void UnsteadyAdvectionDiffusion::v_GenerateSummary(
            SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
    }
}
