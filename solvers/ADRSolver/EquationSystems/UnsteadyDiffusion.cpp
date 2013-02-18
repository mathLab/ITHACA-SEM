///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyDiffusion.cpp
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
// Description: Unsteady diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/UnsteadyDiffusion.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    string UnsteadyDiffusion::className = GetEquationSystemFactory().RegisterCreatorFunction("UnsteadyDiffusion", UnsteadyDiffusion::create);

    UnsteadyDiffusion::UnsteadyDiffusion(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    /**
     * @brief Initialisation object for the unsteady diffusion problem.
     */
    void UnsteadyDiffusion::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        m_session->LoadParameter("epsilon",    m_epsilon,  0.0);

        int npoints = m_fields[0]->GetNpoints();
        if(m_session->DefinesParameter("d00"))
        {
            m_varcoeff[StdRegions::eVarCoeffD00]
            = Array<OneD, NekDouble>(npoints, m_session->GetParameter("d00"));
        }
        if(m_session->DefinesParameter("d11"))
        {
            m_varcoeff[StdRegions::eVarCoeffD11]
            = Array<OneD, NekDouble>(npoints, m_session->GetParameter("d11"));
        }
        if(m_session->DefinesParameter("d22"))
        {
            m_varcoeff[StdRegions::eVarCoeffD22]
            = Array<OneD, NekDouble>(npoints, m_session->GetParameter("d22"));
        }

        switch (m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                string diffName;
        
                m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
                m_diffusion->SetFluxVector(&UnsteadyDiffusion::GetFluxVector, this);
                m_diffusion->InitObject(m_session);
                break;
            }
        
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                // In case of Galerkin explicit diffusion gives an error
                if (m_explicitDiffusion)
                {
                    ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
                }
                // In case of Galerkin implicit diffusion: do nothing
            }
        }
        
        
        if (m_explicitDiffusion)
        {
            m_ode.DefineOdeRhs    (&UnsteadyDiffusion::DoOdeRhs,        this);
            m_ode.DefineProjection(&UnsteadyDiffusion::DoOdeProjection, this);
        }
        else
        {
            m_ode.DefineImplicitSolve(
                                &UnsteadyDiffusion::DoImplicitSolve, this);
        }
    }

    /**
     * @brief Unsteady diffusion problem destructor.
     */
    UnsteadyDiffusion::~UnsteadyDiffusion()
    {
    }

    /* @brief Compute the right-hand side for the unsteady diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyDiffusion::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> > &inarray,
              Array<OneD,        Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;
        
        // Number of fields (variables of the problem)
        int nVariables = inarray.num_elements();
        
        // Number of solution points
        int nSolutionPts = GetNpoints();
        
        // RHS computation using the new advection base class
        m_diffusion->Diffuse(nVariables, 
                             m_fields, 
                             inarray, 
                             outarray);
    }

    /**
     * @brief Compute the projection for the unsteady diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyDiffusion::DoOdeProjection(
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
    
    /** 
     * @brief Implicit solution of the unsteady diffusion problem.
     */
    void UnsteadyDiffusion::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble time,
        const NekDouble lambda)
    {
        StdRegions::ConstFactorMap factors;
        int nvariables = inarray.num_elements();
        int npoints    = m_fields[0]->GetNpoints();
        factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsilon;
        factors[StdRegions::eFactorTau]    = 1.0;
        
        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(npoints, 
                        -factors[StdRegions::eFactorLambda], 
                        inarray[i], 1, 
                        m_fields[i]->UpdatePhys(), 1);
            
            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(), 
                                   NullFlagList, 
                                   factors, 
                                   m_varcoeff);
            
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), 
                                  m_fields[i]->UpdatePhys());
            
            m_fields[i]->SetPhysState(false);
            
            // The solution is Y[i]
            outarray[i] = m_fields[i]->GetPhys();
        }
    }
    
    /** 
     * @brief Return the flux vector for the unsteady diffusion problem.
     */
    void UnsteadyDiffusion::GetFluxVector(
        const int i, 
        const int j,
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &derivatives,
              Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        for(int k = 0; k < flux.num_elements(); ++k)
        {
            Vmath::Zero(GetNpoints(), flux[k], 1);
        }
        Vmath::Vcopy(GetNpoints(), physfield[i], 1, flux[j], 1);
    }
}
