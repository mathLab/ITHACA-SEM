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
// Description: Generic timestepping for PulseWave solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include <MultiRegions/ContField1D.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>
#include <LibUtilities/BasicUtils/Timer.h>


namespace Nektar
{

    /**
     * @class PulseWaveSystem
     *
     * Provides the underlying timestepping framework for pulse wave solvers
     * including the general timestepping routines. 
     */

    /**
     * Processes SolverInfo parameters from the session file and sets up
     * timestepping-specific code.
     * @param   m_Session        Session object to read parameters from.
     */
    PulseWaveSystem::PulseWaveSystem(const LibUtilities::SessionReaderSharedPtr& m_session)
									: UnsteadySystem(m_session)
    {
    }

	/**
     * Destructor
	 */
	PulseWaveSystem::~PulseWaveSystem()
    {
    }
	
	/** 
	 * Initialisation routine for multidomain solver
	 */
    void PulseWaveSystem::v_InitObject()
    {       
		cout << "\n==== Start PulseWaveSystem::v_InitObject() ===="<<endl;

		m_filename = m_session->GetFilename();
		
        // Save the basename of input file name for output details.
        m_sessionName = m_filename;
        m_sessionName = m_sessionName.substr(0, m_sessionName.find_last_of("."));
        
        // Read the geometry and the expansion information
        m_graph = SpatialDomains::MeshGraph::Read(m_filename);
		m_domainsize = m_graph->GetDomain().size();
		
        m_UseContCoeff = false;
		
        // Also read and store the boundary conditions
        m_boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
		::AllocateSharedPtr(m_session, m_graph);
		
        // Read and store history point data
        m_historyPoints = MemoryManager<SpatialDomains::History>
		::AllocateSharedPtr(m_graph);
		
        m_historyPoints->Read(m_filename);
        if (m_historyPoints->GetNumHistoryPoints() > 0
			&& !m_session->DefinesParameter("IO_HistorySteps"))
        {
            cout << "Warning: Must set IO_HistorySteps parameter to enable "
					"output of history points." << endl;
        }
		
        // Set space dimension for use in class
        m_spacedim = m_graph->GetSpaceDimension();
        
        // Setting parameteres for homogenous problems
        m_HomoDirec       = 0;
        m_useFFT          = false;
		m_dealiasing      = false;
        m_HomogeneousType = eNotHomogeneous;
		
        		
        // Determine projectiontype
        if(m_session->DefinesSolverInfo("PROJECTION"))
        {
            std::string ProjectStr
			= m_session->GetSolverInfo("PROJECTION");
			
            if((ProjectStr == "Continuous")||(ProjectStr == "Galerkin")||
               (ProjectStr == "CONTINUOUS")||(ProjectStr == "GALERKIN"))
            {
                m_projectionType = MultiRegions::eGalerkin;
            }
            else if(ProjectStr == "DisContinuous")
            {
                m_projectionType = MultiRegions::eDiscontinuousGalerkin;
            }
            else
            {
                ASSERTL0(false,"PROJECTION value not recognised");
            }
        }
		
        int i;
        int nvariables = m_session->GetVariables().size();
        bool DeclareCoeffPhysArrays = true;
		
        m_fields   = Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);
        m_spacedim = m_graph->GetSpaceDimension()+m_HomoDirec;
        m_expdim   = m_graph->GetMeshDimension();
		
        // Continuous Galerkin projection
        if(m_projectionType == MultiRegions::eGalerkin)
        {
            switch(m_expdim)
            {
                case 1:
                {
					for(i = 0 ; i < m_fields.num_elements(); i++)
					{
						m_fields[i] = MemoryManager<MultiRegions::ContField1D>
						::AllocateSharedPtr(m_session,m_graph,m_session->GetVariable(i));
					}
				break;
                }
				case 2:
                {
                    ASSERTL0(false,"1D Solver; Expansion dimension 2 not allowed");					
                    break;
                }
                case 3:
				{
					ASSERTL0(false,"1D Solver; Expansion dimension 3 not allowed");
					break;
				}
				default:
					ASSERTL0(false,"Expansion dimension not recognised");
					break;
            }
        }
        else // Discontinuous Galerkin
        {
            switch(m_expdim)
            {
                case 1:
                {
					/* In case of a PulseWavePropagation Problem with multiple subdomains use this specialized constructor 
					 * to set up m_vessels and m_traces. i is the variable for the currently processed subdomain
					 */
					if(m_graph->GetDomain().size() > 1)
					{
						cout << "\n-- setting up m_vessels --";
						m_vessels = Array<OneD, MultiRegions::ExpListSharedPtr> (nvariables*m_domainsize);
						const SpatialDomains::CompositeMap domain = (m_graph->GetDomain());
						
						for(i = 0 ; i < m_vessels.num_elements(); i++)
						{
							m_vessels[i] = MemoryManager<MultiRegions::DisContField1D>
								::AllocateSharedPtr(m_session,domain,m_graph,m_session->GetVariable(i%nvariables),(i/nvariables));
						}
						cout << "-- m_vessels are set up --\n"<<endl;
							
							
						cout << "-- checking m_vessels --"<<endl;
						cout << "Number of m_vessels:  "<<m_vessels.num_elements()<<endl;
						for(i = 0 ; i < m_vessels.num_elements(); i++)
						{
							cout << "Elements in m_vessels["<<i<<"]:  "<<m_vessels[i]->GetExpSize()<<endl;
							for(int j = 0 ; j < m_vessels[i]->GetExpSize(); j++)
							{
								cout << "Number of modes in m_vessel["<<i<<"]["<<j<<"]:  "<<m_vessels[i]->GetExp(j)->GetCoeffs().num_elements()<<endl;
							}
							cout << endl;
						}
						for (i = 0; i<m_domainsize; i++)
						{
							cout << "Elements in m_traces["<<i<<"]:  "<<m_vessels[i*nvariables]->GetTrace1D(i)->GetExpSize()<<endl;
						}
						cout << "-- m_vessels are checked --\n"<<endl;
							
						
						// Only needed for output: whole field
						m_outfields = Array<OneD, MultiRegions::ExpListSharedPtr> (nvariables);
						for(i = 0 ; i < m_outfields.num_elements(); i++)
						{
							m_outfields[i] = MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(m_session,m_graph,
											m_session->GetVariable(i));
						}
					}
					else // If only one domain
					{
						for(i = 0 ; i < m_fields.num_elements(); i++)
						{
							m_fields[i] = MemoryManager<MultiRegions
                               ::DisContField1D>::AllocateSharedPtr(m_session,m_graph,m_session->GetVariable(i));
						}
					}					
                    break;
                }
				case 2:
				case 3:
                {
					ASSERTL0(false,"1D Solver; Expansion dimension 2 or 3 not allowed");										
                    break;
                }
				default:
					ASSERTL0(false,"Expansion dimension not recognised");
					break;
            }
			
            // Set up Normals.
            switch(m_expdim)
            {
				case 1:
				{
					if (m_graph->GetDomain().size() > 1)
					{
						m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_spacedim*m_domainsize);
						
						for(int j = 0; j < m_domainsize; ++j)
						{
							for(i = 0; i < m_spacedim; ++i)
							{						
								m_traceNormals[j*m_spacedim+i] = Array<OneD, NekDouble> (m_vessels[j*nvariables]->GetTrace1D(j)->GetExpSize());
							}
							//m_vessels[j*nvariables]->GetTrace1D(j)->GetNormals(m_traceNormals);
						}
					}
					else
					{
						m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
						
						for(i = 0; i < m_spacedim; ++i)
						{
							m_traceNormals[i] = Array<OneD, NekDouble> (GetTraceNpoints());
						}
						m_fields[0]->GetTrace1D()->GetNormals(m_traceNormals);
					}

					break;
				}
				case 2:
				case 3:
                {
					ASSERTL0(false,"1D Solver; Expansion dimension 2 or 3 not allowed");					
                    break;
                }
				default:
                    ASSERTL0(false,"Expansion dimension not recognised");
                    break;
            }
        }
		
        // Set Default Parameter
        m_session->LoadParameter("Time", m_time, 0.0);
        m_session->LoadParameter("TimeStep", m_timestep, 0.01);
        m_session->LoadParameter("NumSteps", m_steps, 0);
        m_session->LoadParameter("IO_CheckSteps", m_checksteps, m_steps);
        m_session->LoadParameter("FinTime", m_fintime, 0);
        m_session->LoadParameter("NumQuadPointsError", m_NumQuadPointsError, 0);
		
		if (m_graph->GetDomain().size() > 1)
		{
			m_fields[0] = m_vessels[0];
			m_fields[1] = m_vessels[1];
		}
		
        // Read in spatial data
        int nq = m_fields[0]->GetNpoints();
        m_spatialParameters = MemoryManager<SpatialDomains::SpatialParameters>::AllocateSharedPtr(nq);
        m_spatialParameters->Read(m_filename);
		
        Array<OneD, NekDouble> x(nq), y(nq), z(nq);
        m_fields[0]->GetCoords(x,y,z);
        m_spatialParameters->EvaluateParameters(x,y,z);
        ScanForHistoryPoints();
		
        // Zero all physical fields initially.
        ZeroPhysFields();
		
        // Load SolverInfo parameters
        m_session->MatchSolverInfo("DIFFUSIONADVANCEMENT","Explicit", m_explicitDiffusion,true);
        m_session->MatchSolverInfo("ADVECTIONADVANCEMENT","Explicit", m_explicitAdvection,true);
		
        // Determine TimeIntegrationMethod to use.
        ASSERTL0(m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"), "No TIMEINTEGRATIONMETHOD defined in session.");
        for (i = 0; i < (int)LibUtilities::SIZE_TimeIntegrationMethod; ++i)
        {
            bool match;
            m_session->MatchSolverInfo("TIMEINTEGRATIONMETHOD", LibUtilities::TimeIntegrationMethodMap[i], match, false);
            if (match)
            {
                m_timeIntMethod = (LibUtilities::TimeIntegrationMethod) i;
                break;
            }
        }
        ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
		
		// If Discontinuous Galerkin determine upwinding method to use
		if (m_projectionType == MultiRegions::eDiscontinuousGalerkin)
		{		 
			for (int i = 0; i < (int)SIZE_UpwindTypePulse; ++i)
			{
				bool match;
				m_session->MatchSolverInfo("UPWINDTYPEPULSE", UpwindTypeMapPulse[i], match, false);
				if (match)
				{
					m_upwindTypePulse = (UpwindTypePulse) i;
					break;
				}
			}
		}
		else
		{
			m_upwindTypePulse = (UpwindTypePulse) 0;
		}
		
		// Load solver- specific parameters
		m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
		// Load constant artery wall thickness
		m_session->LoadParameter("h0", m_h0, 1.0);
		// Load Poission Ratio
		m_session->LoadParameter("nue", m_nue, 0.5);
		// Load blood density
		m_session->LoadParameter("rho", m_rho, 0.5);
		// Load external pressure
		m_session->LoadParameter("pext", m_pext, 0.0);
	
		cout << "==== Finished PulseWaveSystem::v_InitObject() ====\n"<<endl;
	}
	
	
	
	/**
	 * Initialisation routine for multiple subdomain case;
	 * Sets the initial conditions for all fields
	 * in case of only subdomain call SetInitialCondition
	 * form EquationSystem.cpp
	 */
	void PulseWaveSystem::v_DoInitialise()
	{
		// In multidomain case
		if (m_graph->GetDomain().size() > 1)
		{
			NekDouble initialtime = 0.0;
			bool dumpInitialConditions = true;
		
			if (m_session->GetComm()->GetRank() == 0)
			{
				cout << "Initial Conditions:" << endl;
			}
	
			std::string root("InitialConditions[");
			std::string ziel;
			std::string close("]");

			/* Loop over all subdomains to initialize all with the Initial Conditions
			 * read from the inputfile*/
			for (int omega = 0; omega < m_domainsize; omega++)
			{
				m_fields[0] = m_vessels[0+2*omega];
				m_fields[1] = m_vessels[1+2*omega];
				
				if (m_session->GetComm()->GetRank() == 0)
				{
					cout << "Subdomain = " <<omega<<endl;
				}
				
				std:stringstream os;
				std::string omega_str;
				os << omega;
				omega_str = os.str();			
				ziel = root + omega_str + close;
			
				if (m_session->DefinesFunction(ziel))
				{
					int nq = m_fields[0]->GetNpoints();
					Array<OneD,NekDouble> x0(nq);
					Array<OneD,NekDouble> x1(nq);
					Array<OneD,NekDouble> x2(nq);
					m_fields[0]->GetCoords(x0,x1,x2);
				
					for(unsigned int i = 0 ; i < m_fields.num_elements(); i++)
					{
						LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction(ziel, i);
						for(int j = 0; j < nq; j++)
						{
							(m_fields[i]->UpdatePhys())[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],initialtime);
						}
						m_fields[i]->SetPhysState(true);
						m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),
														 m_fields[i]->UpdateCoeffs());
					
						if (m_session->GetComm()->GetRank() == 0)
						{
							cout << "\tField "<< m_session->GetVariable(i)<<": " << ifunc->GetExpression() << endl;
						}
					}									
				}
				else
				{
					ASSERTL0(false,"Initialisation failed");
				}
				
				m_vessels[0+2*omega] = m_fields[0];
				m_vessels[1+2*omega] = m_fields[1];
			}	
			
			/* Also Initialise the bifurcation and junction conditions at t = 0.0
			 * with the Initial Conditions*/
			for (int omega=0; omega<m_domainsize; omega++)
			{
				for (int l=0; l<m_vessels[2*omega]->GetBndConditions().num_elements(); l++)
				{
					if ((m_vessels[2*omega]->GetBndConditions()[l]->GetBoundaryConditionType() == SpatialDomains::eBifurcation)
						|| (m_vessels[2*omega]->GetBndConditions()[l]->GetBoundaryConditionType() == SpatialDomains::eJunction))
					{
						m_vessels[2*omega]->UpdateBndCondExpansion(l)->SetCoeff(m_vessels[2*omega]
																				->GetPhys()[m_vessels[2*omega]->GetPhys().num_elements()-1]);
						m_vessels[2*omega+1]->UpdateBndCondExpansion(l)->SetCoeff(m_vessels[2*omega+1]
																				->GetPhys()[m_vessels[2*omega+1]->GetPhys().num_elements()-1]);
					}
				}
			}
		}
		
		// In single domain case
		else
		{
			SetInitialConditions();
		}
		
		/* Also initialise the StaticArea from the inputfile
		 * and the material properties*/
		StaticArea();
		m_A_0 = m_A_0global[0];
		MaterialProperties();
		m_beta = m_betaglobal[0];
		
    }
	
	
	/**
	 * DoSolve routine for PulseWavePropagation with multiple subdomains
	 * taken from UnsteadySystem
     * Initialises the time integration scheme (as specified in the session
     * file), and perform the time integration.
     */
    void PulseWaveSystem::v_DoSolve()
    {
		if (m_graph->GetDomain().size() > 1)
		{
			cout << "==== Start PulseWaveSystem::v_DoSolve() ===="<<endl;
			NekDouble IntegrationTime = 0.0;
			int i,n,nchk = 1;
			int ncoeffs = 0;
			int npoints = 0;
			int nvariables = 0;
			
			Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  fields(m_domainsize);			
			Array <OneD, Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> > IntScheme(m_domainsize);
			Array <OneD, LibUtilities::TimeIntegrationSolutionSharedPtr > u(m_domainsize);
			int numMultiSteps;
			
			for (int omega = 0; omega < m_domainsize; omega++)
			{
				m_fields[0] = m_vessels[0+2*omega];
				m_fields[1] = m_vessels[1+2*omega];
				
				ncoeffs = m_fields[0]->GetNcoeffs();
				npoints = m_fields[0]->GetNpoints();
				nvariables = m_fields.num_elements();
				
				fields[omega] = Array<OneD, Array<OneD, NekDouble> >(nvariables);
				for(i = 0; i < nvariables; ++i)
				{
					fields[omega][i]  = m_fields[i]->UpdatePhys();
					m_fields[i]->SetPhysState(false);
				}
				
				/* Declare an array of TimeIntegrationSchemes For multi-stage
				 * methods, this array will have just one entry containing the
				 * actual multi-stage method...
				 * For multi-steps method, this can have multiple entries
				 *  - the first scheme will used for the first timestep (this
				 *    is an initialization scheme)
				 *  - the second scheme will used for the second timestep
				 *    (this is an initialization scheme)
				 *  - ...
				 *  - the last scheme will be used for all other time-steps
				 *   (this will be the actual scheme)
				 */
				switch(m_timeIntMethod)
				{
					case LibUtilities::eIMEXdirk_2_3_2:
					case LibUtilities::eIMEXdirk_3_4_3:
					case LibUtilities::eDIRKOrder2:
					case LibUtilities::eDIRKOrder3:
					case LibUtilities::eBackwardEuler:
					case LibUtilities::eForwardEuler:
					case LibUtilities::eClassicalRungeKutta4:
					case LibUtilities::eRungeKutta2_ModifiedEuler:
					case LibUtilities::eRungeKutta2_ImprovedEuler:
					case LibUtilities::eAdamsBashforthOrder1:
					case LibUtilities::eIMEXOrder1:
					{
						numMultiSteps = 1;
						
						IntScheme[omega] = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
						
						LibUtilities::TimeIntegrationSchemeKey IntKey(m_timeIntMethod);
						IntScheme[omega][0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
						
						m_omega = 0;
						u[omega] = IntScheme[omega][0]->InitializeScheme(m_timestep,fields[omega],m_time,m_ode);
						break;
					}
					case LibUtilities::eAdamsBashforthOrder2:
					{
						numMultiSteps = 2;
						
						IntScheme[omega] = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
						
						// Used in the first time step to initalize the scheme
						LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eForwardEuler);
						
						// Used for all other time steps
						LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
						IntScheme[omega][0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
						IntScheme[omega][1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
						
						// Initialise the scheme for the actual time integration scheme
						u[omega] = IntScheme[omega][1]->InitializeScheme(m_timestep,fields[omega],m_time,m_ode);
						break;
					}
					case LibUtilities::eIMEXOrder2:
					{
						numMultiSteps = 2;
						
						IntScheme[omega] = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
						
						// Used in the first time step to initalize the scheme
						LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eIMEXOrder1);
						
						// Used for all other time steps
						LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
						IntScheme[omega][0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
						IntScheme[omega][1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
						
						// Initialise the scheme for the actual time integration scheme
						u[omega] = IntScheme[omega][1]->InitializeScheme(m_timestep,fields[omega],m_time,m_ode);
						break;
					}
						
					case LibUtilities::eIMEXOrder3:
					{
						numMultiSteps = 3;
						
						IntScheme[omega] = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
						
						// Used in the first time step to initalize the scheme
						LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eIMEXdirk_3_4_3);
						LibUtilities::TimeIntegrationSchemeKey IntKey1(LibUtilities::eIMEXdirk_3_4_3);
						
						// Used for all other time steps
						LibUtilities::TimeIntegrationSchemeKey IntKey2(m_timeIntMethod);
						IntScheme[omega][0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
						IntScheme[omega][1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
						IntScheme[omega][2] = LibUtilities::TimeIntegrationSchemeManager()[IntKey2];
						
						// Initialise the scheme for the actual time integration scheme
						u[omega] = IntScheme[omega][2]->InitializeScheme(m_timestep,fields[omega],m_time,m_ode);
						break;
					}
					default:
					{
						ASSERTL0(false,"populate switch statement for integration scheme");
					}
				}
				
				m_vessels[0+2*omega] = m_fields[0];
				m_vessels[1+2*omega] = m_fields[1];
				
			}//end of subdomain loop
			
			
			//Write the history file
			std::string outname = m_session->GetSessionName() + ".his";
			std::ofstream hisFile (outname.c_str());
			
        
			// Time loop
			for(n = 0; n < m_steps; ++n)
			{
				cout << "\nNew TimeStep begins: timestep = "<<n<<"\t";
				cout << " m_time = "<<m_time<<endl;
				
				Timer timer;
				timer.Start();
				
				/* Link the domains:
				 * 1. Set up boundary conditions for all subdomains
				 * then modify the BCs in case of Junction or Bifurcation*/
				LinkSubdomains(fields);				
				
				/* Domains should now be linked successfully
				 * Calculate the domains seperately
				 * this already works*/
				for (int omega = 0; omega < m_domainsize; omega++)
				{
					m_omega = omega;
					m_A_0 = m_A_0global[omega];
					m_beta = m_betaglobal[omega];
					
					m_fields[0] = m_vessels[0+2*omega];
					m_fields[1] = m_vessels[1+2*omega];
					m_fields[0]->GetTrace1D() = m_vessels[2*omega]->GetTrace1D(omega);
					m_fields[1]->GetTrace1D() = m_vessels[2*omega+1]->GetTrace1D(omega);
					
					// Integrate over timestep.
					if( n < numMultiSteps-1)
					{
						// Use initialisation schemes if time step is less than the
						// number of steps in the scheme.
						fields[omega] = IntScheme[omega][n]->TimeIntegrate(m_timestep,u[omega],m_ode);
					}
					else
					{
						cout << "\n Timeintegrating subdomain "<<omega<<"..."<<endl;
						fields[omega] = IntScheme[omega][numMultiSteps-1]->TimeIntegrate(m_timestep,u[omega],m_ode);
						cout << " Subdomain "<<omega<<" done"<<endl;
					}
					
					m_vessels[0+2*omega] = m_fields[0];
					m_vessels[1+2*omega] = m_fields[1];
					
				}//end of subdomain loop
				
				m_time += m_timestep;
				timer.Stop();
				IntegrationTime += timer.TimePerTest(1);
				
				
				// Write out status information.
				if(m_session->GetComm()->GetRank() == 0 && !((n+1)%m_infosteps))
				{
					cout << "Steps: " << n+1
					<< "\t Time: " << m_time
					<< "\t Time-step: " << m_timestep << "\t" << endl;
				}
				
				// Transform data if needed
				if((m_historysteps && !((n+1)%m_historysteps)) || (n&&(!((n+1)%m_checksteps))))
				{
					for (int omega = 0; omega < m_domainsize; omega++)
					{
						m_fields[0] = m_vessels[0+2*omega];
						m_fields[1] = m_vessels[1+2*omega];
						m_fields[0]->GetTrace1D() = m_vessels[2*omega]->GetTrace1D(omega);
						m_fields[1]->GetTrace1D() = m_vessels[2*omega]->GetTrace1D(omega);
			 
						for (i = 0; i < nvariables; ++i)
						{
							m_fields[i]->FwdTrans(fields[omega][i],m_fields[i]->UpdateCoeffs());
							m_fields[i]->SetPhysState(false);
						}
					
						m_vessels[0+2*omega] = m_fields[0];
						m_vessels[1+2*omega] = m_fields[1];
					}
				}
				
				/*/ Write out history data to file
				 * Not needed now, do later...
				 if(m_historysteps && !((n+1)%m_historysteps))
				 {
					WriteHistoryData(hisFile);
				 }*/
				
				
				// Write out checkpoint files; first append all data to m_fields then write m_fields
				if(n&&(!((n+1)%m_checksteps)))
				{
					PrepareMultidomainOutput();
					Checkpoint_Output(nchk++);
				}
				
			}//end of timeintegration
			
			cout <<"\nCFL number              : " << m_cfl << endl;
			cout <<"Time-integration timing : " << IntegrationTime << " s" << endl << endl;
			
			
			// At the end of the time integration, store final solution.
			for (int omega = 0; omega < m_domainsize; omega++)
			{
				m_fields[0] = m_vessels[0+2*omega];
				m_fields[1] = m_vessels[1+2*omega];
				m_fields[0]->GetTrace1D() = m_vessels[2*omega]->GetTrace1D(omega);
				m_fields[1]->GetTrace1D() = m_vessels[2*omega]->GetTrace1D(omega);
		 
				for(i = 0; i < nvariables; ++i)
				{
					if(m_fields[i]->GetPhysState() == false)
					{
						m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
					}
					m_fields[i]->UpdatePhys() = fields[omega][i];
				}
			
				m_vessels[0+2*omega] = m_fields[0];
				m_vessels[1+2*omega] = m_fields[1];
			}
		}//end of case with multiple domains
		else
		{
			UnsteadySystem::v_DoSolve();
		}
		cout << "==== finished PulseWaveSystem::v_DoSolve() ====\n"<<endl;
    }
	
	
	
	/**
	 * Link the subdomains
	 */
	void PulseWaveSystem::LinkSubdomains(Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &fields)
    {
		Array<OneD, NekDouble> Au(3);
		Array<OneD, NekDouble> uu(3);
		Array<OneD, NekDouble> beta(3);
		Array<OneD, NekDouble> A_0(3);
		int in_BCExp = 0;
		NekDouble Q, A_r, u_r;
		
		int nel_p = 0;
		int p = 0;
		int d1 = 0;
		int d2 = 0;
		int p_BCExp = 0;
		int d1_BCExp = 0;
		int d2_BCExp = 0;
		
		/*cout << "Checking for Junctions..."<<endl;
		for (int l=0; l<m_fields[0]->GetBndConditions().num_elements(); l++)
		{
			if (m_fields[0]->GetBndConditions()[l]->GetBoundaryConditionType() == SpatialDomains::eJunction)
			{
				p = m_fields[0]->GetBndConditions()[l]->GetParent();
				d1 = m_fields[0]->GetBndConditions()[l]->GetDaughter1(); 
				cout << "    Junction found:\t VAR=A  parent domain = "<<p<<"\t daughter domain = "<<d1<<""<<endl;
				
				p = m_fields[1]->GetBndConditions()[l]->GetParent();
				d1 = m_fields[1]->GetBndConditions()[l]->GetDaughter1(); 
				cout << "    Junction found:\t VAR=u  parent domain = "<<p<<"\t daughter domain = "<<d1<<""<<endl;
			}
		}*/
		/*cout << "Checking for Bifurcations..."<<endl;
		for (int l=0; l<m_fields[0]->GetBndConditions().num_elements(); l++)
		{
			if (m_fields[0]->GetBndConditions()[l]->GetBoundaryConditionType() == SpatialDomains::eBifurcation)
			{
				p = m_fields[0]->GetBndConditions()[l]->GetParent();
				d1 = m_fields[0]->GetBndConditions()[l]->GetDaughter1(); 
				d2 = m_fields[0]->GetBndConditions()[l]->GetDaughter2();
				cout << "    Bifurcation found:\t VAR=A  parent domain = "<<p<<"\t daughter domain 1 = "<<d1<<"\t daughter domain 2 = "<<d2<<endl;
				
				p = m_fields[1]->GetBndConditions()[l]->GetParent();
				d1 = m_fields[1]->GetBndConditions()[l]->GetDaughter1();
				d2 = m_fields[1]->GetBndConditions()[l]->GetDaughter2();
				cout << "    Bifurcation found:\t VAR=u  parent domain = "<<p<<"\t daughter domain 1 = "<<d1<<"\t daughter domain 2 = "<<d2<<endl;
			}
		}
		cout << endl;*/
		
		
		
		// Set the values of all boundary conditions
		for (int k = 0; k<m_vessels.num_elements(); k++)
		{
			m_vessels[k]->EvaluateBoundaryConditions(m_time);
		}
		
		/*for (int k=0; k<m_domainsize; k++)
		{
			cout << "m_bndCondExpansions before linking "<<k<<endl;
			for (int l=0; l<m_fields[0]->GetBndConditions().num_elements(); l++)
			{
				cout <<" A bndCondExp("<<l<<") = "<<m_vessels[2*k]->UpdateBndCondExpansion(l)->GetCoeff(0)<<"\t";
				cout <<" u bndCondExp("<<l<<") = "<<m_vessels[2*k+1]->UpdateBndCondExpansion(l)->GetCoeff(0)<<endl;
			}
		}*/
		
		
		// Detect special network boundary conditions
		for (int omega = 0; omega < m_domainsize; omega++)
		{
			
			// Check if the endpoint of the domain is a "Bifurcation"
			if (m_vessels[2*omega]->GetBndConditions()[2*omega+1]->GetBoundaryConditionType() == SpatialDomains::eBifurcation)
			{
				cout << " Found the Bifurcation at m_bndcondexp["<<2*omega+1<<"]"<<" at the end of domain "<<omega<<"";
				cout << "\tQuadrature point = "<< fields[omega][0].num_elements()-1 << endl;
				
				// Parent vessel
				nel_p = fields[omega][0].num_elements()-1;
				Au[0] = fields[omega][0][nel_p];
				uu[0] = fields[omega][1][nel_p];
				beta[0] = m_betaglobal[omega][nel_p];
				A_0[0] = m_A_0global[omega][nel_p];
				
				// Daughter vessel 1
				d1 = m_vessels[2*omega]->GetBndConditions()[2*omega+1]->GetDaughter1();
				Au[1] = fields[d1][0][0];
				uu[1] = fields[d1][1][0];
				beta[1] = m_betaglobal[d1][0];
				A_0[1] = m_A_0global[d1][0];
				
				// Daughter vessel 2
				d2 = m_vessels[2*omega]->GetBndConditions()[2*omega+1]->GetDaughter2();
				Au[2] = fields[d2][0][0];
				uu[2] = fields[d2][1][0];
				beta[2] = m_betaglobal[d2][0];
				A_0[2] = m_A_0global[d2][0];
				
				cout << "Before the Riemann solve"<<endl;
				cout << " parent:    Au = "<<Au[0]<<"\tuu = "<<uu[0]<<endl;
				cout << " daughter1: Au = "<<Au[1]<<"\tuu = "<<uu[1]<<endl;
				cout << " daughter2: Au = "<<Au[2]<<"\tuu = "<<uu[2]<<endl;
				
				// Solve the Riemann problem for a bifurcation
				BifurcationRiemann1_to_2(Au, uu, beta, A_0);
				
				cout << "After the Riemann solve"<<endl;
				cout << " parent:    Au = "<<Au[0]<<"\tuu = "<<uu[0]<<endl;
				cout << " daughter1: Au = "<<Au[1]<<"\tuu = "<<uu[1]<<endl;
				cout << " daughter2: Au = "<<Au[2]<<"\tuu = "<<uu[2]<<endl;
				
				// Store the values into the right positions:
				// Parent vessel
				p_BCExp = 2*(omega+1)-1;
				fields[omega][0][nel_p] = Au[0];
				fields[omega][1][nel_p] = uu[0];
				//m_vessels[2*omega]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = Au[0];
				//m_vessels[2*omega+1]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = uu[0];
				
				// Daughter vessel 1
				d1_BCExp = 2*d1;
				m_vessels[2*d1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = Au[1];
				m_vessels[2*d1+1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = uu[1];
				
				// Daughter vessel 2
				d2_BCExp = 2*d2;
				m_vessels[2*d2]->UpdateBndCondExpansion(d2_BCExp)->UpdateCoeffs()[0] = Au[2];
				m_vessels[2*d2+1]->UpdateBndCondExpansion(d2_BCExp)->UpdateCoeffs()[0] = uu[2];
				
			}
			
			// Check if the endpoint of the domain is a "Junction"
			if (m_vessels[2*omega]->GetBndConditions()[2*omega+1]->GetBoundaryConditionType() == SpatialDomains::eJunction)
			{
				cout << " Found the junction at m_bndcondexp["<<2*omega+1<<"]"<<" at the end of domain "<<omega<<"";
				cout << "\tQuadrature point = "<< fields[omega][0].num_elements()-1 << endl;
				
				//Parent vessel
				nel_p = fields[omega][0].num_elements()-1;
				Au[0] = fields[omega][0][nel_p];
				uu[0] = fields[omega][1][nel_p];
				beta[0] = m_betaglobal[omega][nel_p];
				A_0[0] = m_A_0global[omega][nel_p];
				
				//Daughter vessel 1
				d1 = m_vessels[2*omega]->GetBndConditions()[2*omega+1]->GetDaughter1();
				Au[1] = fields[d1][0][0];
				uu[1] = fields[d1][1][0];
				beta[1] = m_betaglobal[d1][0];
				A_0[1] = m_A_0global[d1][0];
				
				cout << "Before the Riemann solve"<<endl;
				cout << " parent: domain = "<<omega<<"   Au = "<<Au[0]<<"\tuu = "<<uu[0]<<endl;
				cout << " daughter1: domain = "<<d1<<"   Au = "<<Au[1]<<"\tuu = "<<uu[1]<<endl;
				
				// Solve the Riemann problem for a junction
				JunctionRiemann(Au, uu, beta, A_0);
				
				cout << "After the Riemann solve"<<endl;
				cout << " parent:    Au = "<<Au[0]<<"\tuu = "<<uu[0]<<endl;
				cout << " daughter1: Au = "<<Au[1]<<"\tuu = "<<uu[1]<<endl;
				
				// Store the values into the right positions:
				// Parent vessel
				p_BCExp = 2*(omega+1)-1;
				fields[omega][0][nel_p] = Au[0];
				fields[omega][1][nel_p] = uu[0];
				//m_vessels[2*omega]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = Au[0];
				//m_vessels[2*omega+1]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = uu[0];
				
				// Daughter vessel 1
				d1_BCExp = 2*d1;
				m_vessels[2*d1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = Au[1];
				m_vessels[2*d1+1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = uu[1];
				
			}
			
			// Check if the endpoint of the domain is a "Merging" flow conditon
			if (m_vessels[2*omega]->GetBndConditions()[2*omega+1]->GetBoundaryConditionType() == SpatialDomains::eMerging)
			{
				cout << " Found the Merging at m_bndcondexp["<<2*omega+1<<"]"<<" at the end of domain "<<omega<<"";
			}
			
		}
		
		for (int k=0; k<m_domainsize; k++)
		{
			cout << "m_bndCondExpansions after linking in domain "<<k<<endl;
			for (int l=0; l<m_fields[0]->GetBndConditions().num_elements(); l++)
			{
				cout <<" A bndCondExp("<<l<<") = "<<m_vessels[2*k]->UpdateBndCondExpansion(l)->GetCoeff(0)<<"\t";
				cout <<" u bndCondExp("<<l<<") = "<<m_vessels[2*k+1]->UpdateBndCondExpansion(l)->GetCoeff(0)<<endl;
			}
		}
		
		/*for (int k=0; k<m_domainsize; k++)
		{
			cout << "m_bndCondExpansions after linking in domain "<<k<<endl;
			for (int l=0; l<2; l++)
			{
				cout <<" A bndCondExp("<<l<<") = "<<m_vessels[2*k]->UpdateBndCondExpansion(k*2+l)->GetCoeff(0)<<"\t";
				cout <<" u bndCondExp("<<l<<") = "<<m_vessels[2*k+1]->UpdateBndCondExpansion(k*2+l)->GetCoeff(0)<<endl;
			}
		}*/
		
	}
	
	
		
	
	
	/**
	 * Bifurcation Riemann problem called to link three subdomains by a Bifurcation
	 */
	void PulseWaveSystem::BifurcationRiemann1_to_2(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu, 
												   Array<OneD, NekDouble> &beta, Array<OneD, NekDouble> &A_0)
	{
		NekDouble rho = m_rho;
		Array<OneD, NekDouble> W(3);
		Array<OneD, NekDouble> W_Au(3);
		Array<OneD, NekDouble> P_Au(3);
		Array<OneD, NekDouble> f(6);
		Array<OneD, NekDouble> g(6);
		Array<OneD, NekDouble> tmp(6);
		Array<OneD, Array<OneD, NekDouble> > inv_J(6);
		for (int i=0; i<6; i++)
		{
			inv_J[i] = Array<OneD, NekDouble> (6);
		}
		NekDouble k = 0.0;
		NekDouble k1 = 0.0;
		NekDouble k2 = 0.0;
		NekDouble k3 = 0.0;
		
		int proceed = 1;
		int iter = 0;
		int MAX_ITER = 7;
		
		//calculated from input
		W[0] = uu[0] + 4*sqrt(beta[0]/(2*rho))*(sqrt(sqrt(Au[0])) - sqrt(sqrt(A_0[0])));
		W[1] = uu[1] - 4*sqrt(beta[1]/(2*rho))*(sqrt(sqrt(Au[1])) - sqrt(sqrt(A_0[1])));
		W[2] = uu[2] - 4*sqrt(beta[2]/(2*rho))*(sqrt(sqrt(Au[2])) - sqrt(sqrt(A_0[2])));
		
		//Tolerances for the algorithm
		NekDouble Tol = 1.0e-10;
		NekDouble Small = 1.0e-200;
		
		
		//Newton Iteration
		while ((proceed) && (iter < MAX_ITER))
		{
			iter = iter+1;
			
			//Calculate the constraint vector, six equations:
			//3 characteristic variables, mass conservation, 
			//total pressure
			W_Au[0] = 4*sqrt(beta[0]/(2*rho))*(sqrt(sqrt(Au[0])) - sqrt(sqrt(A_0[0])));
			W_Au[1] = 4*sqrt(beta[1]/(2*rho))*(sqrt(sqrt(Au[1])) - sqrt(sqrt(A_0[1])));
			W_Au[2] = 4*sqrt(beta[2]/(2*rho))*(sqrt(sqrt(Au[2])) - sqrt(sqrt(A_0[2])));
			
			P_Au[0] = beta[0]*(sqrt(Au[0]) - sqrt(A_0[0]));
			P_Au[1] = beta[1]*(sqrt(Au[1]) - sqrt(A_0[1]));
			P_Au[2] = beta[2]*(sqrt(Au[2]) - sqrt(A_0[2]));
			
			f[0] = uu[0] + W_Au[0] - W[0];
			f[1] = uu[1] - W_Au[1] - W[1];
			f[2] = uu[2] - W_Au[2] - W[2];
			f[3] = Au[0]*uu[0] - Au[1]*uu[1] - Au[2]*uu[2];
			f[4] = uu[0]*uu[0] + 2.0/rho*P_Au[0] - uu[1]*uu[1] - 2.0/rho*P_Au[1];
			f[5] = uu[0]*uu[0] + 2.0/rho*P_Au[0] - uu[2]*uu[2] - 2.0/rho*P_Au[2];
			
			
			// Calculate the wave speed at each vessel
			NekDouble c1 = sqrt(beta[0]/(2*rho))*sqrt(sqrt(Au[0]));
			NekDouble c2 = sqrt(beta[1]/(2*rho))*sqrt(sqrt(Au[1]));
			NekDouble c3 = sqrt(beta[2]/(2*rho))*sqrt(sqrt(Au[2]));
			
			
			// Inverse Jacobian matrix J(x[n])^(-1), is already inverted here analytically
			k = c1*Au[1]*c3+Au[0]*c3*c2+Au[2]*c1*c2;
			k1 = (c1-uu[0])*k;
			inv_J[0][0] = (-c2*uu[0]*c3*Au[0]+Au[2]*c2*c1*c1+Au[1]*c1*c1*c3)/k1;
			inv_J[0][1] = Au[1]*(c2-uu[1])*c1*c3/k1;
			inv_J[0][2] = Au[2]*(c3-uu[2])*c1*c2/k1;
			inv_J[0][3] = c1*c2*c3/k1;
			inv_J[0][4] = -0.5*c1*Au[1]*c3/k1;
			inv_J[0][5] = -0.5*Au[2]*c1*c2/k1;
			
			k2 = (c2+uu[1])*k;
			inv_J[1][0] = Au[0]*(c1+uu[0])*c2*c3/k2;
			inv_J[1][1] = (c1*uu[1]*c3*Au[1]+Au[2]*c1*c2*c2+c3*c2*c2*Au[0])/k2;
			inv_J[1][2] = -Au[2]*(c3-uu[2])*c1*c2/k2;
			inv_J[1][3] = -c1*c2*c3/k2;
			inv_J[1][4] = -0.5*(c1*Au[2]+Au[0]*c3)*c2/k2;
			inv_J[1][5] = 0.5*Au[2]*c1*c2/k2;
			
			k3 = (c3+uu[2])*k;
			inv_J[2][0] = Au[0]*(c1+uu[0])*c2*c3/k3;
			inv_J[2][1] = -Au[1]*(c2-uu[1])*c1*c3/k3;
			inv_J[2][2] = (c1*c2*uu[2]*Au[2]+c1*Au[1]*c3*c3+c2*c3*c3*Au[0])/k3;
			inv_J[2][3] = -c1*c2*c3/k3;
			inv_J[2][4] = 0.5*c1*Au[1]*c3/k3;
			inv_J[2][5] = -0.5*(Au[1]*c1+c2*Au[0])*c3/k3;
			
			inv_J[3][0] = Au[0]*(Au[0]*c3*c2-uu[0]*c3*Au[1]-uu[0]*c2*Au[2])/k1;
			inv_J[3][1] = -Au[0]*Au[1]*(c2-uu[1])*c3/k1;
			inv_J[3][2] = -Au[0]*Au[2]*(c3-uu[2])*c2/k1;
			inv_J[3][3] = -Au[0]*c3*c2/k1;
			inv_J[3][4] = 0.5*Au[0]*Au[1]*c3/k1;
			inv_J[3][5] = 0.5*Au[0]*c2*Au[2]/k1;
			
			inv_J[4][0] = Au[0]*Au[1]*(c1+uu[0])*c3/k2;
			inv_J[4][1] = -Au[1]*(c1*Au[1]*c3+c1*uu[1]*Au[2]+c3*uu[1]*Au[0])/k2;
			inv_J[4][2] = -Au[2]*Au[1]*(c3-uu[2])*c1/k2;
			inv_J[4][3] = -c1*Au[1]*c3/k2;
			inv_J[4][4] = -0.5*Au[1]*(c1*Au[2]+Au[0]*c3)/k2;
			inv_J[4][5] = 0.5*Au[2]*Au[1]*c1/k2;
			
			inv_J[5][0] = Au[0]*Au[2]*(c1+uu[0])*c2/k3;
			inv_J[5][1] = -Au[2]*Au[1]*(c2-uu[1])*c1/k3;
			inv_J[5][2] = -Au[2]*(Au[2]*c1*c2+c1*uu[2]*Au[1]+c2*uu[2]*Au[0])/k3;
			inv_J[5][3] = -Au[2]*c1*c2/k3;
			inv_J[5][4] = 0.5*Au[2]*Au[1]*c1/k3;
			inv_J[5][5] = -0.5*Au[2]*(Au[1]*c1+c2*Au[0])/k3;
			
			
			// Solve the system by multiplying the Jacobian with the vector f:
			// g = (inv_J)*f
			for (int j=0; j<6; j++)
			{
				tmp[j] =0.0;
				g[j] = 0.0;
			}
			
			for (int j=0; j<6; j++)
			{
				for (int i=0; i<6; i++)
				{
					tmp[j] = inv_J[j][i]*f[i];
					g[j] += tmp[j];
				}
			}
			
			
			// Update the solution: x_new = x_old - dx
			uu[0] = uu[0] - g[0];
			uu[1] = uu[1] - g[1];
			uu[2] = uu[2] - g[2];
			Au[0] = Au[0] - g[3];
			Au[1] = Au[1] - g[4];
			Au[2] = Au[2] - g[5];
			
			// Check if the error of the solution is smaller than Tol
			if ((g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3]+ g[4]*g[4] + g[5]*g[5]) < Tol)
			{
				proceed = 0;
			}
			
			// Check if solver converges
			if (iter >= MAX_ITER)
			{
				ASSERTL0(false,"Riemann solver did not converge");
			}
			
		}
		
		cout << "-------------------------------------------------------------------------"<<endl;
		cout << "| Bifurcation Riemann solver converged, iterations: "<<iter<<"\t\t\t\t|"<< scientific<<endl;
		cout << "| uu[parent] = "<<uu[0]<<"   uu[daughter1] = "<<uu[1]<<"   uu[daughter2] = "<<uu[2]<<"\t|"<<endl;
		cout << "| Au[parent] = "<<Au[0]<<"   Au[daughter1] = "<<Au[1]<<"   Au[daughter2] = "<<Au[2]<<"\t|"<<endl;
		
		
		// Checks
		// whether the flow rate over the bifurcation matches
		Array<OneD, NekDouble> Q(4);
		Q[0] = Au[0]*uu[0];
		Q[1] = Au[1]*uu[1];
		Q[2] = Au[2]*uu[2];
		Q[3] = Q[0]-Q[1]-Q[2];
		cout << "| mass difference = "<<Q[3]<<"\t\t\t\t\t\t\t|"<<endl;
		
		// whether the total pressure is conserved
		Array<OneD, NekDouble> P(5);
		P[0] = uu[0]*uu[0] + 2.0/rho*(beta[0]*(sqrt(Au[0]) - sqrt(A_0[0])));
		P[1] = uu[1]*uu[1] + 2.0/rho*(beta[1]*(sqrt(Au[1]) - sqrt(A_0[1])));
		P[2] = uu[2]*uu[2] + 2.0/rho*(beta[2]*(sqrt(Au[2]) - sqrt(A_0[2])));
		P[3] = P[0]-P[1];
		P[4] = P[0]-P[2]; 
		cout << "| pressure difference 1-2 = "<<P[3]<<"\tpressure difference 1-3 = "<<P[4]<<"\t|"<<endl;
		cout << "-------------------------------------------------------------------------"<<endl;
		cout.setf(ios::floatfield);
	}
	
	
	/**
	 * \brief
	 * Solves the Riemann problem at an interdomain junction by assuming subsonic flow at
	 * both sides of the boundary and by applying conservation of mass and continuity of the
	 * total pressure \f$ \frac{p}{rho} + \frac{u^{2}}{2}. \f$
	 */
	void PulseWaveSystem::JunctionRiemann(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu,
										  Array<OneD, NekDouble> &beta, Array<OneD, NekDouble> &A_0)
	{		
		NekDouble rho = m_rho;
		Array<OneD, NekDouble> W(2);
		Array<OneD, NekDouble> W_Au(2);
		Array<OneD, NekDouble> P_Au(2);
		Array<OneD, NekDouble> f(4);
		Array<OneD, NekDouble> g(4);
		Array<OneD, NekDouble> tmp(4);
		Array<OneD, Array<OneD, NekDouble> > inv_J(4);
		for (int i=0; i<4; i++)
		{
			inv_J[i] = Array<OneD, NekDouble> (4);
		}
		NekDouble k = 0.0;
		NekDouble k1 = 0.0;
		NekDouble k2 = 0.0;
		
		int proceed = 1;
		int iter = 0;
		int MAX_ITER = 7;
		NekDouble Tol = 1.0e-10;
		
		// Calculated from input
		W[0] = uu[0] + 4*sqrt(beta[0]/(2*rho))*(sqrt(sqrt(Au[0])) - sqrt(sqrt(A_0[0])));
		W[1] = uu[1] - 4*sqrt(beta[1]/(2*rho))*(sqrt(sqrt(Au[1])) - sqrt(sqrt(A_0[1])));
		
		while((proceed) && (iter < MAX_ITER))
		{
			iter = iter+1;
			
			// Calculate the constraint vector, 4 equations:
			// 2 characteristic variables, mass conservation, 
			// total pressure
			W_Au[0] = 4*sqrt(beta[0]/(2*rho))*(sqrt(sqrt(Au[0])) - sqrt(sqrt(A_0[0])));
			W_Au[1] = 4*sqrt(beta[1]/(2*rho))*(sqrt(sqrt(Au[1])) - sqrt(sqrt(A_0[1])));
			
			P_Au[0] = beta[0]*(sqrt(Au[0]) - sqrt(A_0[0]));
			P_Au[1] = beta[1]*(sqrt(Au[1]) - sqrt(A_0[1]));
			
			f[0] = uu[0] + W_Au[0] - W[0];
			f[1] = uu[1] - W_Au[1] - W[1];
			f[2] = Au[0]*uu[0] - Au[1]*uu[1];
			f[3] = uu[0]*uu[0] + 2.0/rho*P_Au[0] - uu[1]*uu[1] - 2.0/rho*P_Au[1];
			
			// Calculate the wave speed at each vessel
			NekDouble cl = sqrt(beta[0]/(2*rho))*sqrt(sqrt(Au[0]));
			NekDouble cr = sqrt(beta[1]/(2*rho))*sqrt(sqrt(Au[1]));
			
			// Inverse Jacobian matrix J(x[n])^(-1), is already inverted here analytically
			k = (cl*Au[1]+Au[0]*cr);
			k1 = (cl-uu[0])*k;
			inv_J[0][0] = (Au[1]*cl*cl-cr*uu[0]*Au[0])/k1;
			inv_J[0][1] = Au[1]*(cr-uu[1])*cl/k1;
			inv_J[0][2] = cl*cr/k1;
			inv_J[0][3] = -0.5*cl*Au[1]/k1;
			
			k2 = (cr+uu[1])*k;
			inv_J[1][0] = Au[0]*(cl+uu[0])*cr/k2;
			inv_J[1][1] = (cl*uu[1]*Au[1]+cr*cr*Au[0])/k2;
			inv_J[1][2] = -cl*cr/k2;
			inv_J[1][3] = -0.5*Au[0]*cr/k2;
			
			inv_J[2][0] = Au[0]*(Au[0]*cr-uu[0]*Au[1])/k1;
			inv_J[2][1] = -Au[0]*Au[1]*(cr-uu[1])/k1;
			inv_J[2][2] = -Au[0]*cr/k1;
			inv_J[2][3] = 0.5*Au[1]*Au[0]/k1;
			
			inv_J[3][0] = Au[0]*Au[1]*(cl+uu[0])/k2;
			inv_J[3][1] = -Au[1]*(cl*Au[1]+uu[1]*Au[0])/k2;
			inv_J[3][2] = -cl*Au[1]/k2;
			inv_J[3][3] = -0.5*Au[1]*Au[0]/k2;
			
			// Solve the system by multiplying the Jacobian with the vector f:
			// g = (inv_J)*f
			for (int j=0; j<4; j++)
			{
				tmp[j] =0.0;
				g[j] = 0.0;
			}
			
			for (int j=0; j<4; j++)
			{
				for (int i=0; i<4; i++)
				{
					tmp[j] = inv_J[j][i]*f[i];
					g[j] += tmp[j];
				}
			}
			
			// Update solution: x_new = x_old - dx
			uu[0] -= g[0];
			uu[1] -= g[1];
			Au[0] -= g[2];
			Au[1] -= g[3];
			
			// Check if the error of the solution is smaller than Tol.
			if((g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3]) < Tol)      
				proceed = 0;
		}
		
		if(iter >= MAX_ITER)
		{
			ASSERTL0(false,"Riemann solver did not converge");
		}
		
		cout << "----------------------------------------------------------------"<<endl;
		cout << "| Junction Riemann solver converged, iterations: "<<iter<<"\t\t|"<< scientific<<endl;
		cout << "| uu[parent] = "<<uu[0]<<"   uu[daughter1] = "<<uu[1]<<"\t|"<<endl;
		cout << "| Au[parent] = "<<Au[0]<<"   Au[daughter1] = "<<Au[1]<<"\t|"<<endl;
		
		
		// Checks
		// whether the flow rate over the bifurcation matches
		Array<OneD, NekDouble> Q(4);
		Q[0] = Au[0]*uu[0];
		Q[1] = Au[1]*uu[1];
		Q[3] = Q[0]-Q[1];
		cout << "| mass difference = "<<Q[3]<<"\t\t\t\t|"<<endl;
		
		// whether the total pressure is conserved
		Array<OneD, NekDouble> P(5);
		P[0] = uu[0]*uu[0] + 2.0/rho*(beta[0]*(sqrt(Au[0]) - sqrt(A_0[0])));
		P[1] = uu[1]*uu[1] + 2.0/rho*(beta[1]*(sqrt(Au[1]) - sqrt(A_0[1])));
		P[3] = P[0]-P[1];
		cout << "| pressure difference 1-2 = "<<P[3]<<"\t\t\t|"<<endl;
		cout << "----------------------------------------------------------------"<<endl;
		cout.setf(ios::floatfield);
		
	}
	
	
	
	/**
	 * Gets the Area at static equilibrium specified in the inputfile
	 * in case of multiple domains different areas can be specified by
	 * using A[domain] in the inputfile
	 */
	void PulseWaveSystem::StaticArea(void)
    {
		m_A_0global = Array<OneD, Array<OneD, NekDouble> >(m_domainsize);
		int nq = 0;
		
		for (int omega = 0; omega < m_domainsize; omega++)
		{
			if (m_domainsize == 1)
				nq = m_fields[2*omega]->GetNpoints();
			else 
				nq = m_vessels[2*omega]->GetNpoints();
			
			m_A_0global[omega] = Array<OneD, NekDouble>(nq); 
			
			std::string velStr[1];
			std::string root("A_0[");
			std::string close("]");
			std:stringstream os;
			std::string omega_str;
			
			os << omega;
			omega_str = os.str();			
			velStr[0] = root + omega_str + close;
			
			LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction("A_0",velStr[0]);
			EvaluateFunction(m_A_0global[omega],ifunc,0.0);
		}
    }
	
	
	/**
	 * Gets the Material Properties of the artery specified in the inputfile.
	 * in case of multiple domains different areas can be specified by using
	 * E0[domain] in the inputfile 
	 */
	void PulseWaveSystem::MaterialProperties()
    {
		m_betaglobal = Array<OneD, Array<OneD, NekDouble> >(m_domainsize);
		int nq = 0;
		
		for (int omega = 0; omega < m_domainsize; omega++)
		{
			if (m_domainsize == 1)
				nq = m_fields[2*omega]->GetNpoints();
			else 
				nq = m_vessels[2*omega]->GetNpoints();
			
			m_betaglobal[omega] = Array<OneD, NekDouble>(nq); 			
			std::string velStr[1];
			
			// Two cases: one can either specify E0 or beta in the inputfile
			// Case 1: Read E0
			/*for (int j = 0; j<m_betaglobal[omega].num_elements(); j++)
			{
				m_betaglobal[omega][j] = sqrt(3.1415)*m_h0*m_betaglobal[omega][j]/((1-m_nue*m_nue)*m_A_0global[omega][j]);
			}*/
			
			// Case 2: read beta directly
			if(m_session->GetFunction("MaterialProperties","beta[0]"))
			{
				std::string root("beta[");
				std::string close("]");
				std:stringstream os;
				std::string omega_str;
				
				os << omega;
				omega_str = os.str();			
				velStr[0] = root + omega_str + close;
				
				LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction("MaterialProperties",velStr[0]);
				EvaluateFunction(m_betaglobal[omega],ifunc,0.0);
			}
		}
	}
	
	
	/**
     * Prepare the multifield output; put the multiple 
	 * subdomains together to one domain to use the 
	 * standard routines for writing the output
     */
    void PulseWaveSystem::PrepareMultidomainOutput(void)
    {		
		Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(2);
		int totCoeffs = 0;
		for(int omega = 0; omega < m_domainsize; omega++)
		{
			totCoeffs += m_vessels[2*omega]->UpdateCoeffs().num_elements();
		}
		for (int i = 0; i < 2; i++)
		{
			fieldcoeffs[i] = Array<OneD, NekDouble> (totCoeffs);
		}
		
		int offset =0;
		for(int omega = 0; omega < m_domainsize; omega++)
		{
			for (int i = 0; i < m_vessels[2*omega]->UpdateCoeffs().num_elements(); i++)
			{
				fieldcoeffs[0][i+offset] = m_vessels[2*omega]->UpdateCoeffs()[i];
				fieldcoeffs[1][i+offset] = m_vessels[2*omega+1]->UpdateCoeffs()[i];
			}
			offset += m_vessels[2*omega]->UpdateCoeffs().num_elements();
		}
		
		m_outfields[0]->UpdateCoeffs() = fieldcoeffs[0];					
		m_outfields[1]->UpdateCoeffs() = fieldcoeffs[1];
		
		m_fields[0] = m_outfields[0];
		m_fields[1] = m_outfields[1];
    }
	
	
	/**
	 * Write the .fld file at the end of the simulation
	 */
	void PulseWaveSystem::v_Output(void)
    {
        std::string outname = m_sessionName + ".fld";
		if (m_domainsize > 1)
		{
			PrepareMultidomainOutput();
		}
        WriteFld(outname); 
	}
	
	
	/**
	 * Compute the error in the L2-norm
	 * @param   field           The field to compare.
	 * @param   exactsoln       The exact solution to compare with.
	 * @param   Normalised      Normalise L2-error.
	 * @returns                 Error in the L2-norm.
	 */
	NekDouble PulseWaveSystem::v_L2Error(unsigned int field,
										 const Array<OneD, NekDouble> &exactsoln,
										 bool Normalised)
	{    		
		NekDouble L2error = -1.0;
		if (m_domainsize > 1)
		{
			for (int omega = 0; omega < m_domainsize; omega++)
			{
				m_fields[field] = m_vessels[field];
				
				if(m_NumQuadPointsError == 0)
				{
					if(m_fields[field]->GetPhysState() == false)
					{
						m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
												  m_fields[field]->UpdatePhys());
					}
					
					if(exactsoln.num_elements())
					{
						L2error = m_fields[field]->L2(exactsoln);
					}
					else if (m_session->DefinesFunction("ExactSolution"))
					{
						Array<OneD, NekDouble> exactsoln(m_fields[field]->GetNpoints());
						
						LibUtilities::EquationSharedPtr vEqu
						= m_session->GetFunction("ExactSolution",field);
						EvaluateFunction(exactsoln,vEqu,m_time);
						
						L2error = m_fields[field]->L2(exactsoln);
					}
					else
					{
						L2error = m_fields[field]->L2();
					}
				}
				else
				{
					ASSERTL0(false,"ErrorExtraPoints not allowed for this solver");
				}
			}
			return L2error;
		}
		else
		{
			return EquationSystem::v_L2Error(field, exactsoln);
		}
	}
	
	
	/**
	 * Compute the error in the L_inf-norm
	 * @param   field           The field to compare.
	 * @param   exactsoln       The exact solution to compare with.
	 * @returns                 Error in the L_inft-norm.
	 */
	NekDouble PulseWaveSystem::v_LinfError(unsigned int field,
										   const Array<OneD, NekDouble> &exactsoln)
	{
		NekDouble Linferror = -1.0;		
		if (m_domainsize > 1)
		{
			for (int omega = 0; omega < m_domainsize; omega++)
			{
				m_fields[field] = m_vessels[field];
			
				if(m_NumQuadPointsError == 0)
				{
					if(m_fields[field]->GetPhysState() == false)
					{
						m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
												  m_fields[field]->UpdatePhys());
					}
				
					if(exactsoln.num_elements())
					{
						Linferror = m_fields[field]->Linf(exactsoln);
					}
					else if (m_session->DefinesFunction("ExactSolution"))
					{
						Array<OneD, NekDouble> exactsoln(m_fields[field]->GetNpoints());
					
						LibUtilities::EquationSharedPtr vEqu
						= m_session->GetFunction("ExactSolution", field);
						EvaluateFunction(exactsoln,vEqu,m_time);
					
						Linferror = m_fields[field]->Linf(exactsoln);
					}
					else
					{
						Linferror = 0.0;
					}
				
				}
				else
				{
					ASSERTL0(false,"ErrorExtraPoints not allowed for this solver");
				}
			}
			return Linferror;
		}
		else 
		{
			return EquationSystem::v_LinfError(field, exactsoln);
		}
	}
	
    
}
