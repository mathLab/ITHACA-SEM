///////////////////////////////////////////////////////////////////////////////
//
// File PulseWaveSystem.cpp
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
// Description: Generic timestepping for APE solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>

namespace Nektar
{
    /**
     * @class PulseWaveSystem
     *
     * Provides the underlying timestepping framework for pulse wave solvers
     * including the general timestepping routines. This class is not intended
     * to be directly instantiated, but rather is a base class on which to
     * define pulse wave solvers
     *

    /**
     * Processes SolverInfo parameters from the session file and sets up
     * timestepping-specific code.
     * @param   m_Session        Session object to read parameters from.
     */
    PulseWaveSystem::PulseWaveSystem(
									 const LibUtilities::SessionReaderSharedPtr& m_session)
									: UnsteadySystem(m_session)
    {
    }

    void PulseWaveSystem::v_InitObject()
    {
        EquationSystem::v_InitObject();

        // Load SolverInfo parameters
        m_session->MatchSolverInfo("DIFFUSIONADVANCEMENT","Explicit",
								   m_explicitDiffusion,true);
        m_session->MatchSolverInfo("ADVECTIONADVANCEMENT","Explicit",
								   m_explicitAdvection,true);
		
        // Determine TimeIntegrationMethod to use.
        ASSERTL0(m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"),
				 "No TIMEINTEGRATIONMETHOD defined in session.");
        int i;
        for (i = 0; i < (int)LibUtilities::SIZE_TimeIntegrationMethod; ++i)
        {
            bool match;
            m_session->MatchSolverInfo("TIMEINTEGRATIONMETHOD",
									   LibUtilities::TimeIntegrationMethodMap[i], match, false);
            if (match)
            {
                m_timeIntMethod = (LibUtilities::TimeIntegrationMethod) i;
                break;
            }
        }
        ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod,
				 "Invalid time integration type.");
		
		
		// if discontinuous Galerkin determine numerical flux to use
		 if (m_projectionType == MultiRegions::eDiscontinuousGalerkin)
		 {
			 //ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPEPULSE"),
			 //		 "No UPWINDTYPEPULSE defined in session.");
		 
			 int i;
			 for (i = 0; i < (int)SIZE_UpwindTypePulse; ++i)
			 {
				 bool match;
				 m_session->MatchSolverInfo("UPWINDTYPEPULSE",
											UpwindTypeMapPulse[i], match, false);
		 
				 if (match)
				 {
					 m_upwindTypePulse = (UpwindTypePulse) i;
					 break;
				 }
			 }
			 //ASSERTL0(i != (int) SIZE_UpwindTypePulse,
			 //		 "Invalid upwind type Pulse.");
		 }
		 else
		 {
			 m_upwindTypePulse = (UpwindTypePulse) 0;
		 }
		
		
		m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
	    
		// Load constant artery wall thickness
		m_session->LoadParameter("h0", m_h0, 1.0);
		
		// Load Poission Ratio
		m_session->LoadParameter("nue", m_nue, 0.5);
		
		// Load blood density
		m_session->LoadParameter("rho", m_rho, 0.5);
		
		// Load blood density
		m_session->LoadParameter("pext", m_pext, 0.0);    }
	

    /**
     *
     */
	/// Destructor
	PulseWaveSystem::~PulseWaveSystem()
    {

    }


    
}
