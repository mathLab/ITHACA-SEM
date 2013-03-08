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

namespace Nektar
{
    string UnsteadyInviscidBurger::className = GetEquationSystemFactory().RegisterCreatorFunction("UnsteadyInviscidBurger", UnsteadyInviscidBurger::create, "Inviscid Burger equation");
    
    UnsteadyInviscidBurger::UnsteadyInviscidBurger(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }
    
    /**
     * @brief Initialisation object for the inviscid Burger equation.
     */
    void UnsteadyInviscidBurger::v_InitObject()
    {
        // Call to the initialisation object of UnsteadySystem
        UnsteadySystem::v_InitObject();
                
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
                m_advection = SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
                m_advection->SetFluxVector   (&UnsteadyInviscidBurger::GetFluxVector, this);
                break;
            }
            // Discontinuous field 
            case MultiRegions::eDiscontinuous:
            {
                string advName;
                string riemName;
                
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advection = SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
                m_advection->SetFluxVector   (&UnsteadyInviscidBurger::GetFluxVector, this);
                
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(riemName);
                m_riemannSolver->AddScalar("Vn", &UnsteadyInviscidBurger::GetNormalVelocity, this);
                
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
        m_advection->Advect(nVariables, m_fields, advVel, inarray, outarray);
        
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

