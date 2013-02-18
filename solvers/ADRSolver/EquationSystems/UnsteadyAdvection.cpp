/////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvection.cpp
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
// Description: Unsteady linear advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <ADRSolver/EquationSystems/UnsteadyAdvection.h>

namespace Nektar
{
    string UnsteadyAdvection::className = GetEquationSystemFactory().RegisterCreatorFunction("UnsteadyAdvection", UnsteadyAdvection::create, "Unsteady Advection equation.");

    UnsteadyAdvection::UnsteadyAdvection(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    /**
     * @brief Initialisation object for the unsteady linear advection equation.
     */
    void UnsteadyAdvection::v_InitObject()
    {
        // Call to the initialisation object of UnsteadySystem
        UnsteadySystem::v_InitObject();
                
        // Define the normal velocity fields
        if (m_fields[0]->GetTrace())
        {
            m_traceVn  = Array<OneD, NekDouble>(GetTraceNpoints());
        }
        
        // Read the advection velocities from session file 
        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");
        
        // Resize the advection velocities vector to dimension of the problem
        vel.resize(m_spacedim);

        // Store in the global variable m_velocity the advection velocities
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        EvaluateFunction(vel, m_velocity, "AdvectionVelocity");

        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Continuous field 
            case MultiRegions::eGalerkin:
            {
                string advName;
                m_session->LoadSolverInfo("AdvectionType", advName, "NonConservative");
                m_advection = SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
                m_advection->SetFluxVector   (&UnsteadyAdvection::GetFluxVector, this);
                break;
            }
            // Discontinuous field 
            case MultiRegions::eDiscontinuous:
            {
                string advName;
                string riemName;
                
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advection = SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
                m_advection->SetFluxVector   (&UnsteadyAdvection::GetFluxVector, this);
                
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(riemName);
                m_riemannSolver->AddScalar("Vn", &UnsteadyAdvection::GetNormalVelocity, this);
                
                m_advection->SetRiemannSolver(m_riemannSolver);
                m_advection->InitObject      (m_session, m_fields);
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
            m_ode.DefineOdeRhs     (&UnsteadyAdvection::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyAdvection::DoOdeProjection, this);
        }
        // Otherwise it gives an error (no implicit integration)
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    /**
     * @brief Unsteady linear advection equation destructor.
     */
    UnsteadyAdvection::~UnsteadyAdvection()
    {
    }

    /**
     * @brief Get the normal velocity for the linear advection equation.
     */
    Array<OneD, NekDouble> &UnsteadyAdvection::GetNormalVelocity()
    {
        // Number of trace (interface) points
        int nTracePts = GetTraceNpoints();
        
        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);
        
        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);
        
        // Compute the normal velocity
        for (int i = 0; i < m_velocity.num_elements(); ++i)
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
     * @brief Compute the right-hand side for the linear advection equation.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvection::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;
        
        // Number of fields (variables of the problem)
        int nVariables     = inarray.num_elements();
        
        // Number of solution points
        int nSolutionPts    = GetNpoints();
        
        // RHS computation using the new advection base class
        m_advection->Advect(nVariables, 
                            m_fields, 
                            m_velocity, 
                            inarray, 
                            outarray);
        
        // Negate the RHS
        for (i = 0; i < nVariables; ++i)
        {
            Vmath::Neg(nSolutionPts, outarray[i], 1);
        }
    }

    /**
     * @brief Compute the projection for the linear advection equation.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvection::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;
        
        // Number of fields (variables of the problem)
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
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs(),0.0);
                for(i = 0; i < nVariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
                break;
            }
                
            default:
                ASSERTL0(false,"Unknown projection scheme");
                break;
        }
    }
    
    /**
     * @brief Return the flux vector for the linear advection equation.
     * 
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyAdvection::GetFluxVector(
        const int i, 
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),
                 "Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(), 
                        physfield[i], 1, 
                        m_velocity[j], 1, 
                        flux[j], 1);
        }
    }

    void UnsteadyAdvection::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
    }
}
