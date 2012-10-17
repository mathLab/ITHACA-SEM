///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadySystem.cpp
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
// Description: Generic timestepping for unsteady compressible flow solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include <CompressibleFlowSolver/EquationSystems/UnsteadySystem.h>

namespace Nektar
{
  /**
   * @class UnsteadySystem
   *
   * Provides the underlying timestepping framework for unsteady compressible
   * flow solvers including the general timestepping routines. This class is not
   * intended to be directly instantiated, but rather is a base class on which
   * to define compressible flow solvers, e.g. Euler, Euler with artificial
   * diffusion, Navier-Stokes.
   *
   * For details on implementing unsteady solvers see
   * \ref sectionADRSolverModuleImplementation 
   */
  
  /**
   * Processes SolverInfo parameters from the session file and sets up
   * timestepping-specific code.
   * @param   pSession        Session object to read parameters from.
   */
  UnsteadySystem::UnsteadySystem(
          const LibUtilities::SessionReaderSharedPtr& pSession)
    : EquationSystem(pSession)
  {
  }

  void UnsteadySystem::v_InitObject()
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
	ASSERTL0(false,"No Continuous Galerkin implemented for compressible flow solver.");
      }
    
    m_session->LoadParameter("Gamma",m_gamma,1.4);
    m_session->LoadParameter("GasConstant",m_GasConstant,287.058);
    m_session->LoadParameter("CFLParameter",m_cfl,0.0);
    
    // Load generic input parameters
    m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
  }
  
  /**
   *
   */
  UnsteadySystem::~UnsteadySystem()
  {
    
  }

  /**
   * Initialises the time integration scheme (as specified in the session
   * file), and perform the time integration.
   */
  void UnsteadySystem::v_DoSolve()
  {
    int i,n,nchk = 0;
    int nq = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nvariables = m_fields.num_elements();
    int n_elements = m_fields[0]->GetExpSize();

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
    Array<OneD, Array<OneD, NekDouble> >   fieldsOld(nvariables);
    Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);

    for(i = 0; i < nvariables; ++i)
      {
	m_fields[i]->SetPhysState(false);
	fields[i]  = m_fields[i]->UpdatePhys();
	fieldsOld[i] = Array<OneD, NekDouble>(nq);
      }

    // Saving initial condition and starting residuals file
    NekDouble Rsd;
    std::string outname = m_sessionName + ".rsd";
    ofstream outfile(outname.c_str());
    outfile << "# it   rho   rhou   rhov   E" << endl;
    WriteTecplotFile(-1,"tec",false);
    
    // CFL parameters definition
    // Declaring the DG-CFL limit
    Array<OneD,NekDouble> CFL(n_elements,0.0);
    const Array<OneD,int> ExpOrder = GetNumExpModesPerExp();
    Array<OneD,NekDouble> DGStability = GetStabilityLimitVector(ExpOrder);
    
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
    NekDouble RKStability;

    switch(m_timeIntMethod)
      {
      case LibUtilities::eForwardEuler:
      case LibUtilities::eClassicalRungeKutta4:
	{
	  numMultiSteps = 1;

 	  // save it somewhere in the integration class
 	  RKStability = 2.784*m_cfl;
	  Vmath::Sdiv(n_elements,RKStability,DGStability,1,CFL,1);

	  IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

	  LibUtilities::TimeIntegrationSchemeKey IntKey(m_timeIntMethod);
	  IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
	  
	  u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,m_ode);
	  break;
	}
      case LibUtilities::eAdamsBashforthOrder2:
	{
	  numMultiSteps = 2;

 	  // save it somewhere in the integration class
 	  RKStability = m_cfl;
 	  Vmath::Sdiv(n_elements,RKStability,DGStability,1,CFL,1);
	  
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

    outname = m_session->GetFilename() + ".his";
    std::ofstream hisFile (outname.c_str());
    
    n = 0;

    // Check final condition
    if(m_fintime>0.0 && m_steps>0)
      ASSERTL0(false,"Final condition not unique: fintime>0.0 & Nsteps>0");
    // Check Timestep condition
    if(m_timestep>0.0 && m_cfl>0.0)
      ASSERTL0(false,"Timestep not unique: timestep>0.0 & CFL>0.0");

    // Perform integration in time.
    while(n < m_steps || m_time<m_fintime)
      {

	// Save old solution
        for(i = 0; i < nvariables; ++i)
	  Vmath::Vcopy(nq,fields[i],1,fieldsOld[i],1);
 
	if(m_cfl>0.0)
	  m_timestep = GetTimeStep(fields, ExpOrder, CFL, RKStability);
	
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
	if(m_time+m_timestep>m_fintime && m_fintime>0.0)
	  m_timestep = m_fintime - m_time;
	m_time += m_timestep;

	// Write out status information.
	if((!((n+1)%m_infosteps) || n==m_steps || m_time==m_fintime) && m_comm->GetRank() == 0)
        {
	    cout << "Steps: " << n+1 << "\t Time: " << m_time << "\t TimeStep: " << m_timestep << endl;
        }
	
	// Write out checkpoint files.
	if(n&&(!((n+1)%m_checksteps)) || (n==m_steps && m_steps!=0) || (m_time>=m_fintime && m_fintime>0.0))
	  {
	    
	    // update m_fields
	    for(i = 0; i < nvariables; ++i)
	      {
		Vmath::Vcopy(nq,fields[i],1,m_fields[i]->UpdatePhys(),1);
	      }
	    
	    for(i = 0; i < nvariables; ++i)
	      {
		m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
	      }
	    Checkpoint_Output(nchk++);
	  }

	// Calculate and save residuals
	outfile << n << " ";
	for(i = 0; i < nvariables; ++i)
	  { 
	    Vmath::Vsub(nq,fields[i],1,fieldsOld[i],1,fieldsOld[i],1);
	    Vmath::Vmul(nq,fieldsOld[i],1,fieldsOld[i],1,fieldsOld[i],1);
	    Rsd = sqrt(Vmath::Vsum(nq,fieldsOld[i],1));
	    outfile << Rsd << " ";
	  }
	outfile << endl;

	n++;
      }
  }
  
  /**
   *
   */
  void UnsteadySystem::v_DoInitialise()
  {
    SetInitialConditions();
  }

  /**
   *
   */
  void UnsteadySystem::v_PrintSummary(std::ostream &out)
  {
    EquationSystem::v_PrintSummary(out);
    out << "\tUpwind Type     : " << UpwindTypeMap[m_upwindType] << endl;
    out << "\tAdvection       : " << (m_explicitAdvection ? "explicit" : "implicit") << endl;
    out << "\tIntegration Type: " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
    out << "\tTime Step       : " << m_timestep << endl;
    out << "\tNo. of Steps    : " << m_steps << endl;
    out << "\tFinal Time      : " << m_fintime << endl;
    out << "\tCFL parameter   : " << m_cfl << endl;
    out << "\tCheckpoints     : " << m_checksteps << " steps" << endl;
    out << "\tVariables       : rho      should be in field[0]" <<endl;
    out << "\t                  rhou     should be in field[1]" <<endl;
    out << "\t                  rhov     should be in field[2]" <<endl;
    out << "\t                  E (rhoe) should be in field[3]" <<endl;
  }

  /**
   *
   */
  void UnsteadySystem::v_NumericalFlux(
					       Array<OneD, Array<OneD, NekDouble> > &physfield,
					       Array<OneD, Array<OneD, NekDouble> > &numflux)
  {
      ASSERTL0(false, "This function is not implemented for this equation.");
  }
  
  
  /**
   *
   */
  void UnsteadySystem::v_NumericalFlux(
					       Array<OneD, Array<OneD, NekDouble> > &physfield,
					       Array<OneD, Array<OneD, NekDouble> > &numfluxX,
					       Array<OneD, Array<OneD, NekDouble> > &numfluxY )
  {
      ASSERTL0(false, "This function is not implemented for this equation.");
  }
  
  
  /**
   *
   */
  void UnsteadySystem::v_GetFluxVector(const int i, const int j,
					       Array<OneD, Array<OneD, NekDouble> > &physfield,
					       Array<OneD, Array<OneD, NekDouble> > &flux)
  {
      for(int k = 0; k < flux.num_elements(); ++k)
      {
          Vmath::Zero(GetNpoints(), flux[k], 1);
      }
      Vmath::Vcopy(GetNpoints(), physfield[i], 1, flux[j], 1);
  }
  
}
