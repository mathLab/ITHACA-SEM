///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterSystem.cpp
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
// Description: Generic timestepping for shallow water solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>

namespace Nektar
{
    /**
     * @class ShallowWaterSystem
     *
     * Provides the underlying timestepping framework for shallow water flow solvers
     * including the general timestepping routines. This class is not intended
     * to be directly instantiated, but rather is a base class on which to
     * define shallow water solvers, e.g. SWE, Boussinesq, linear and nonlinear versions.
     *
     * For details on implementing unsteady solvers see
     * \ref sectionADRSolverModuleImplementation 
     */

    /**
     * Processes SolverInfo parameters from the session file and sets up
     * timestepping-specific code.
     * @param   pSession        Session object to read parameters from.
     */
    ShallowWaterSystem::ShallowWaterSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : EquationSystem(pSession)
    {
    }

    void ShallowWaterSystem::v_InitObject()
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
	if (m_projectionType == MultiRegions::eDiscontinuous)
	  {
	    ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
		     "No UPWINDTYPE defined in session.");
	    
	    int i;
	    for (i = 0; i < (int)SIZE_UpwindType; ++i)
	      {
		bool match;
		m_session->MatchSolverInfo("UPWINDTYPE",
					  UpwindTypeMap[i], match, false);
		if (match)
		  {
		    m_upwindType = (UpwindType) i;
		    break;
		  }
	      }
	    ASSERTL0(i != (int) SIZE_UpwindType,
		     "Invalid upwind type.");
	  }
	else
	  {
	    m_upwindType = (UpwindType) 0;
	  }
	
	    
	// Load generic input parameters
        m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);

	// Load acceleration of gravity
	m_session->LoadParameter("Gravity", m_g, 9.81);

	// input/output in primitive variables
	m_primitive = true;
    }


    /**
     *
     */
    ShallowWaterSystem::~ShallowWaterSystem()
    {

    }


    /**
     * Initialises the time integration scheme (as specified in the session
     * file), and perform the time integration.
     */
  void ShallowWaterSystem::v_DoSolve()
  {
        int i,n,nchk = 0;
        int nq = m_fields[0]->GetTotPoints();
        int nvariables = m_fields.num_elements();

        // Set up wrapper to fields data storage.
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);

        for(i = 0; i < nvariables; ++i)
        {
            m_fields[i]->SetPhysState(false);
            fields[i]  = m_fields[i]->UpdatePhys();
        }
	
	
        // Declare an array of TimeIntegrationSchemes For multi-stage
        // methods, this array will have just one entry containing the
        // actual multi-stage method...
        // For multi-steps method, this can have multiple entries
        //  - the first scheme will used for the first timestep (this
        //    is an initialization scheme)
        //  - the second scheme will used for the second timestep
        //    (this is an initialization scheme)
        //  - ...
        //  - the last scheme will be used for all other time-steps
        //    (this will be the actual scheme)

        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
        int numMultiSteps;

        switch(m_timeIntMethod)
        {
	case LibUtilities::eForwardEuler:
        case LibUtilities::eClassicalRungeKutta4:
	  {
                numMultiSteps = 1;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                LibUtilities::TimeIntegrationSchemeKey IntKey(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

                u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
        case LibUtilities::eAdamsBashforthOrder2:
            {
                numMultiSteps = 2;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eClassicalRungeKutta4);

                // Used for all other time steps
                LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
	default:
            {
                ASSERTL0(false,"populate switch statement for integration scheme");
            }
        }

        // Perform integration in time.
        for(n = 0; n < m_steps; ++n)
        {
            // Integrate over timestep.
            if( n < numMultiSteps-1)
            {
                // Use initialisation schemes if time step is less than the
                // number of steps in the scheme.
                fields = IntScheme[n]->TimeIntegrate(m_timestep,u,m_ode);
            }
            else
            {
                fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,m_ode);
            }

            // Increment time.
            m_time += m_timestep;

            // Write out status information.
            if(!((n+1)%m_infosteps))
            {
                cout << "Steps: " << n+1 << "\t Time: " << m_time << "\t " << endl;
            }

            // Write out checkpoint files.
            if(n&&(!((n+1)%m_checksteps)))
            {
	      
	      // update m_fields
	      for(i = 0; i < nvariables; ++i)
		{
		  Vmath::Vcopy(nq,fields[i],1,m_fields[i]->UpdatePhys(),1);
		}
	      
	      // go to primitive variables 
	      v_ConservativeToPrimitive();
	      
	      for(i = 0; i < nvariables; ++i)
                {
		  m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
		}
	      Checkpoint_Output(nchk++);

	    v_PrimitiveToConservative();
	    }
        }

        // At the end of the time integration, store final solution.
	// update m_fields
	for(i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(nq,fields[i],1,m_fields[i]->UpdatePhys(),1);
	  }
	
	// to to primitive variables
	v_ConservativeToPrimitive();

	for(i = 0; i < nvariables; ++i)
	  {
	    m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
	  }

    }

    /**
     *
     */
    void ShallowWaterSystem::v_DoInitialise()
    {
        SetInitialConditions();
	EvaluateWaterDepth();
	EvaluateCoriolis();
	v_PrimitiveToConservative();
    }

    /**
     *
     */
    void ShallowWaterSystem::v_PrintSummary(std::ostream &out)
    {
        EquationSystem::v_PrintSummary(out);
	out << "\tUpwind Type     : " << UpwindTypeMap[m_upwindType] << endl;
        out << "\tAdvection       : " << (m_explicitAdvection ? "explicit" : "implicit") << endl;
	out << "\tIntegration Type: " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
        out << "\tTime Step       : " << m_timestep << endl;
        out << "\tNo. of Steps    : " << m_steps << endl;
        out << "\tCheckpoints     : " << m_checksteps << " steps" << endl;
	out << "\tVariables       : eta should be in field[0]" <<endl;
	out << "\t                  u   should be in field[1]" <<endl;
	out << "\t                  v   should be in field[2]" <<endl;
    }


    /**
     *
     */
    void ShallowWaterSystem::v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        ASSERTL0(false, "This function is not implemented for this equation.");
    }


    /**
     *
     */
    void ShallowWaterSystem::v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                Array<OneD, Array<OneD, NekDouble> > &numfluxY )
    {
        ASSERTL0(false, "This function is not implemented for this equation.");
    }


    /**
     *
     */
    void ShallowWaterSystem::v_GetFluxVector(const int i, const int j,
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        for(int k = 0; k < flux.num_elements(); ++k)
        {
            Vmath::Zero(GetNpoints(),flux[k],1);
        }
        Vmath::Vcopy(GetNpoints(),physfield[i],1,flux[j],1);
    }

  
  void ShallowWaterSystem::v_PrimitiveToConservative()
  {
    ASSERTL0(false, "This function is not implemented for this equation.");
  }
  
  void ShallowWaterSystem::v_ConservativeToPrimitive()
  {
    ASSERTL0(false, "This function is not implemented for this equation.");
  }

  void ShallowWaterSystem::EvaluateWaterDepth(void)
  {
    EvaluateFunction("d",m_depth,"WaterDepth");
  }
  
  
  void ShallowWaterSystem::EvaluateCoriolis(void)
  {
    EvaluateFunction("f",m_coriolis,"Coriolis");
  }

 

}
