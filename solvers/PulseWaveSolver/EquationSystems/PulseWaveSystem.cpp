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
// Description: Generic timestepping for Pulse Wave Solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <MultiRegions/ContField1D.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>
#include <LibUtilities/BasicUtils/Timer.h>


namespace Nektar
{
    
    /**
     *  @class PulseWaveSystem
     *
     *  Initialises the arterial subdomains in m_vessels and sets up
     *  all domain-linking conditions (bifurcations, junctions,
     *  merging flows). Detects the network structure and assigns
     *  boundary conditons. Also provides the underlying timestepping
     *  framework for pulse wave solvers including the general
     *  timestepping routines.
     */
    
    /**
     *  Processes SolverInfo parameters from the session file and sets
     *  up timestepping-specific code.
     *
     *  @param   m_Session        Session object to read parameters from.
     */
    PulseWaveSystem::PulseWaveSystem(const LibUtilities::SessionReaderSharedPtr& m_session)
        : UnsteadySystem(m_session)
    {
    }
    
    /**
     *  Destructor
     */
    PulseWaveSystem::~PulseWaveSystem()
    {
    }
    
    /** 
     *  Initialisation routine for multidomain solver. Sets up the
     *  expansions for every arterial segment (m_vessels) and for one
     *  complete field m_outfield which is needed to write the
     *  postprocessing output. Also determines which upwind strategy
     *  is used (currently only upwinding scheme available) and reads
     *  blodd flow specific parameters from the inputfile
     * 
     */
    void PulseWaveSystem::v_InitObject()
    {       
        m_filename = m_session->GetFilename();
	
        // Save the basename of input file name for output details.
        m_sessionName = m_filename;
        m_sessionName = m_sessionName.substr(0, m_sessionName.find_last_of("."));
        
        // Read the geometry and the expansion information
        m_graph = SpatialDomains::MeshGraph::Read(m_session);
        m_nDomains = m_graph->GetDomain().size();
		
        // Determine projectiontype
        ASSERTL0(m_session->MatchSolverInfo("Projection","DisContinuous"),"Pulse solver only set up for Discontinuous projections");
        m_projectionType = MultiRegions::eDiscontinuous;
        ASSERTL0(m_graph->GetMeshDimension() == 1,"Pulse wave solver only set up for expansion dimension equal to 1");
        
        int i;
        m_nVariables =  m_session->GetVariables().size();
	
        m_fields  = Array<OneD, MultiRegions::ExpListSharedPtr> (m_nVariables);
        m_vessels = Array<OneD, MultiRegions::ExpListSharedPtr> (m_nVariables*m_nDomains);

        m_spacedim = m_graph->GetSpaceDimension();
        m_expdim   = m_graph->GetMeshDimension();
        m_HomoDirec			= 0;
        m_useFFT			= false;
        m_dealiasing	  	        = false;
        m_specHP_dealiasing             = false;
        m_HomogeneousType   = eNotHomogeneous;
                        
        const std::vector<SpatialDomains::CompositeMap> domain = m_graph->GetDomain();
			
        SpatialDomains::BoundaryConditions Allbcs(m_session, m_graph);

        // Set up domains and put geometry to be only one space dimension. 
        int cnt = 0;
        bool SetToOneSpaceDimension = true;

        if(m_session->DefinesCmdLineArgument("SetToOneSpaceDimension"))
        {
            std::string cmdline = m_session->GetCmdLineArgument("SetToOneSpaceDimension");
            if(boost::to_upper_copy(cmdline) == "FALSE")
            {
                SetToOneSpaceDimension  = false;
            }
        }

        for(i = 0 ; i < m_nDomains; ++i)
        {
            for(int j = 0; j < m_nVariables; ++j)
            {
                m_vessels[cnt++] = MemoryManager<MultiRegions::DisContField1D>
                    ::AllocateSharedPtr(m_session, m_graph, domain[i],
                                        Allbcs,
                                        m_session->GetVariable(j),
                                        SetToOneSpaceDimension);
            }
        }

        // Reset coeff and phys space to be continuous over all domains
        int totcoeffs = 0;
        int totphys   = 0;
        for(i = 0; i < m_nDomains; ++i)
        {
            totcoeffs += m_vessels[i*m_nVariables]->GetNcoeffs();
            totphys   += m_vessels[i*m_nVariables]->GetTotPoints();
        }
        

        for(int n = 0; n < m_nVariables; ++n)
        {
            Array<OneD, NekDouble> coeffs(totcoeffs,0.0);
            Array<OneD, NekDouble> phys(totphys,0.0);
            Array<OneD, NekDouble> tmpcoeffs,tmpphys;

            m_vessels[n]->SetCoeffsArray(coeffs);
            m_vessels[n]->SetPhysArray(phys);

            int cnt  = m_vessels[n]->GetNcoeffs();
            int cnt1 = m_vessels[n]->GetTotPoints();

            for(i = 1; i < m_nDomains; ++i)
            {
                m_vessels[i*m_nVariables+n]->SetCoeffsArray(tmpcoeffs = coeffs + cnt);
                m_vessels[i*m_nVariables+n]->SetPhysArray   (tmpphys = phys + cnt1);
                cnt  += m_vessels[i*m_nVariables+n]->GetNcoeffs();
                cnt1 += m_vessels[i*m_nVariables+n]->GetTotPoints();                
            }
        }
        
        // Set Default Parameter
        m_session->LoadParameter("Time",     m_time, 0.0);
        m_session->LoadParameter("TimeStep", m_timestep, 0.01);
        m_session->LoadParameter("NumSteps", m_steps, 0);
        m_session->LoadParameter("IO_CheckSteps", m_checksteps, m_steps);
        m_session->LoadParameter("FinTime", m_fintime, 0);
        m_session->LoadParameter("NumQuadPointsError", m_NumQuadPointsError, 0);
		
        m_fields[0] = m_vessels[0];
        m_fields[1] = m_vessels[1];
	
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
        // Load solver- specific parameters
        m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
        // Load blood density
        m_session->LoadParameter("rho", m_rho, 0.5);
        // Load external pressure
        m_session->LoadParameter("pext", m_pext, 0.0);
        
        int nq = 0;
        /**
         *  Gets the Material Properties of each arterial segment
         *  specified in the inputfile from section MaterialProperties
	 *  Also gets the Area at static equilibrium A_0 specified in the
	 *  inputfile. 
         *
         * Having found these points also extract the values at the
         * trace points and the normal direction consistent with the
         * left adjacent definition of Fwd and Bwd
	 */
        m_beta       = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
        m_beta_trace = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
        m_A_0        = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
        m_A_0_trace  = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
        m_trace_fwd_normal = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
        
        for (int omega = 0; omega < m_nDomains; omega++)
        {
            nq = m_vessels[2*omega]->GetNpoints();
            m_fields[0] = m_vessels[2*omega];
            
            m_beta[omega] = Array<OneD, NekDouble>(nq);
            EvaluateFunction("beta", m_beta[omega],"MaterialProperties",0.0,omega);

            m_A_0[omega] = Array<OneD, NekDouble>(nq); 
            EvaluateFunction("A_0", m_A_0[omega],"A_0",0.0,omega);

            int nqTrace = GetTraceTotPoints();

            m_beta_trace[omega] = Array<OneD, NekDouble>(nqTrace);
            m_fields[0]->ExtractTracePhys(m_beta[omega],m_beta_trace[omega]);
            
            m_A_0_trace[omega] = Array<OneD, NekDouble>(nqTrace);
            m_fields[0]->ExtractTracePhys(m_A_0[omega],m_A_0_trace[omega]);


            if(SetToOneSpaceDimension)
            {
                m_trace_fwd_normal[omega] = Array<OneD, NekDouble>(nqTrace,0.0);
                
                MultiRegions::ExpListSharedPtr trace = m_fields[0]->GetTrace(); 
                int nelmt_trace = trace->GetExpSize();
                
                Array<OneD, Array<OneD, NekDouble> > normals(nelmt_trace);
            
                for(int i = 0 ; i < nelmt_trace; ++i)
                {
                    normals[i] = m_trace_fwd_normal[omega]+i;
                }
 
                // need to set to 1 for consistency since boundary
                // conditions may not have coordim=1
               trace->GetExp(0)->GetGeom()->SetCoordim(1); 
                
                trace->GetNormals(normals);
            }
        }
    }
	
    
    /**
     *  Initialisation routine for multiple subdomain case. Sets the
     *  initial conditions for all arterial subdomains read from the
     *  inputfile. Sets the material properties and the A_0 area for
     *  all subdomains and fills the domain-linking boundary
     *  conditions with the initial values of their domain. 
     */
    void PulseWaveSystem::v_DoInitialise()
    {
        
        if (m_session->GetComm()->GetRank() == 0)
        {
            cout << "Initial Conditions:" << endl;
        }
	
        /* Loop over all subdomains to initialize all with the Initial
         * Conditions read from the inputfile*/
        for (int omega = 0; omega < m_nDomains; omega++)
        {
            m_fields[0] = m_vessels[2*omega];
            m_fields[1] = m_vessels[2*omega+1];
            
            if (m_session->GetComm()->GetRank() == 0)
            {
                cout << "Subdomain = " <<omega<<endl;
            }
            
            SetInitialConditions(0.0,0,omega);
        }	
        // Reset to first definition
        m_fields[0] = m_vessels[0];
        m_fields[1] = m_vessels[1];
            
    }

    /**
     * NEEDS Updating:
     *
     *  DoSolve routine for PulseWavePropagation with multiple
     *  subdomains taken from UnsteadySystem and modified for
     *  multidomain case. Initialises the time integration scheme (as
     *  specified in the session file), and perform the time
     *  integration. Within the timestepping loop the following is
     *  done: 1. Link all arterial segments according to the network
     *  structure, solve the Riemann problem between different
     *  arterial segments and assign the values to the boundary
     *  conditions (LinkSubdomains) 2. Every arterial segment is
     *  solved independentl for this timestep. This is done by handing
     *  the solution vector \f$ \mathbf{u} \f$ and the right hand side
     *  m_ode, which is the PulseWavePropagation class in this example
     *  over to the time integration scheme
     */
    void PulseWaveSystem::v_DoSolve()
    {
        NekDouble IntegrationTime = 0.0;
        int i,n,nchk = 1;
	
        Array<OneD, Array<OneD,NekDouble> >  fields(m_nVariables);			
        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
	
        for(int i = 0; i < m_nVariables; ++i)
        {
            fields[i]  = m_vessels[i]->UpdatePhys();
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
        case LibUtilities::eForwardEuler:
        case LibUtilities::eClassicalRungeKutta4:
        case LibUtilities::eRungeKutta2_ModifiedEuler:
        case LibUtilities::eRungeKutta2_ImprovedEuler:
        case LibUtilities::eAdamsBashforthOrder1:
            {
                int numMultiSteps = 1;

                IntScheme = Array<OneD, LibUtilities::
                    TimeIntegrationSchemeSharedPtr>(numMultiSteps);
                    
                LibUtilities::TimeIntegrationSchemeKey IntKey(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
                
                u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
        default:
            {
                ASSERTL0(false,"populate switch statement for integration scheme");
            }
        }
        
        // Time loop
        for(n = 0; n < m_steps; ++n)
        {				
            Timer timer;
            timer.Start();
            
#if 0 
            /* Link the domains:
             * 1. Set up boundary conditions for all subdomains
             * then modify the BCs in case of Junction or Bifurcation*/
            LinkSubdomains(fields);				
#endif            

            fields = IntScheme[0]->TimeIntegrate(m_timestep,u,m_ode);
                
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
            if(!((n+1)%m_checksteps))
            {
                for (i = 0; i < m_nVariables; ++i)
                {
                    int cnt = 0;
                    for (int omega = 0; omega < m_nDomains; omega++)
                    {
                        m_vessels[omega*m_nVariables+i]->FwdTrans(fields[i]+cnt,
                                     m_vessels[omega*m_nVariables+i]->UpdateCoeffs());
                        cnt += m_vessels[omega*m_nVariables+i]->GetTotPoints();
                    }
                }
                CheckPoint_Output(nchk++);
            }
            
        }//end of timeintegration
	
        //Copy Array To Vessel Phys Fields
        for(int i = 0; i < m_nVariables; ++i)
        {
            Vmath::Vcopy(fields[i].num_elements(), fields[i],1,m_vessels[i]->UpdatePhys(),1);
        }
        
        cout <<"\nCFL safety factor     : " 
             << m_cflSafetyFactor << endl;
        cout <<"Time-integration timing : " 
             << IntegrationTime << " s" << endl << endl;
    }	
	
#if 0 

    /**
     *  Links the subdomains for special network boundary
     *  conditions such as "Bifurcation", "Junction" and "Merging
     *  Flow" conditions. Calculates the upwinded Au & uu at
     *  boundaries connected to other subdomains by solving the
     *  appropriate Riemann problem. The solution satisfies
     *  conservation of mass and continuity of the total pressure.
     */
	void PulseWaveSystem::LinkSubdomains(Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &fields)
        {		
		Array<OneD, NekDouble> Au(3);
		Array<OneD, NekDouble> uu(3);
		Array<OneD, NekDouble> beta(3);
		Array<OneD, NekDouble> A_0(3);
		
		int nel_p = 0;
		int nel_d1 = 0;
		int nel_d2 = 0;
		int d1 = 0;
		int d2 = 0;
		int p_BCExp = 0;
		int d1_BCExp = 0;
		int d2_BCExp = 0;
		
		// Set the values of all boundary conditions
		for (int k = 0; k < m_vessels.num_elements(); k++)
		{
                    m_vessels[k]->EvaluateBoundaryConditions(m_time);
		}
		
		// Detect special network boundary conditions
		for (int omega = 0; omega < m_nDomains; omega++)
		{
                    // "Bifurcation": Check if the endpoint of the domain is a "Bifurcation" condition
                    if (m_vessels[2*omega]->GetBndConditions()[1]->GetBoundaryConditionType() == SpatialDomains::eBifurcation)
                    {
                        // Parent vessel
                        nel_p   = fields[omega][0].num_elements()-1;
                        Au[0]   = fields[omega][0][nel_p];
                        uu[0]   = fields[omega][1][nel_p];
                        beta[0] = m_betaglobal[omega][nel_p];
                        A_0[0]  = m_A_0global[omega][nel_p];
                        p_BCExp = 1;
			
                        // Daughter vessel 1
                        d1      = m_vessels[2*omega]->GetBndConditions()[1]->GetDaughter1();
                        Au[1]   = fields[d1][0][0];
                        uu[1]   = fields[d1][1][0];
                        beta[1] = m_betaglobal[d1][0];
                        A_0[1]  = m_A_0global[d1][0];
                        d1_BCExp = 0;
                        
                        // Daughter vessel 2
                        d2 = m_vessels[2*omega]->GetBndConditions()[1]->GetDaughter2();
                        Au[2] = fields[d2][0][0];
                        uu[2] = fields[d2][1][0];
                        beta[2] = m_betaglobal[d2][0];
                        A_0[2] = m_A_0global[d2][0];
                        d2_BCExp = 0;
			
                        // Solve the Riemann problem for a bifurcation
                        BifurcationRiemann1_to_2(Au, uu, beta, A_0);
			
                        // Store the values into the right positions:
                        // Parent vessel
                        m_vessels[2*omega]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = Au[0];
                        m_vessels[2*omega+1]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = uu[0];
                        
                        // Daughter vessel 1
                        m_vessels[2*d1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = Au[1];
                        m_vessels[2*d1+1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = uu[1];
                        
                        // Daughter vessel 2
                        m_vessels[2*d2]->UpdateBndCondExpansion(d2_BCExp)->UpdateCoeffs()[0] = Au[2];
                        m_vessels[2*d2+1]->UpdateBndCondExpansion(d2_BCExp)->UpdateCoeffs()[0] = uu[2];
                        
                    }
                    
                    // "Junction": Check if the endpoint of the domain is a "Junction" condition
                    if (m_vessels[2*omega]->GetBndConditions()[1]->GetBoundaryConditionType() == SpatialDomains::eJunction)
                    {
                        //Parent vessel
                        nel_p = fields[omega][0].num_elements()-1;
                        Au[0] = fields[omega][0][nel_p];
                        uu[0] = fields[omega][1][nel_p];
                        beta[0] = m_betaglobal[omega][nel_p];
                        A_0[0] = m_A_0global[omega][nel_p];
			
                        //Daughter vessel 1
                        d1 = m_vessels[2*omega]->GetBndConditions()[1]->GetDaughter1();
                        Au[1] = fields[d1][0][0];
                        uu[1] = fields[d1][1][0];
                        beta[1] = m_betaglobal[d1][0];
                        A_0[1] = m_A_0global[d1][0];
			
                        // Solve the Riemann problem for a junction
                        JunctionRiemann(Au, uu, beta, A_0);
			
                        // Store the values into the right positions:
                        // Parent vessel
                        p_BCExp = 1;
                        m_vessels[2*omega]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = Au[0];
                        m_vessels[2*omega+1]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = uu[0];
                        
                        // Daughter vessel 1
                        d1_BCExp = 0;
                        m_vessels[2*d1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = Au[1];
                        m_vessels[2*d1+1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = uu[1];
                        
                    }
                    
                    // "Merging Flow": Check if the startpoint of the domain is a "Merging Flow" conditon
                    if (m_vessels[2*omega]->GetBndConditions()[0]->GetBoundaryConditionType() == SpatialDomains::eMerging)
                    {
                        //cout << "\nFound the Merging Flow at m_bndcondexp["<<2*omega<<"]"<<" at the beginning of domain "<<omega<<endl;
			
                        // Parent vessel (merging vessel)
                        Au[0] = fields[omega][0][0];
                        uu[0] = fields[omega][1][0];
                        beta[0] = m_betaglobal[omega][0];
                        A_0[0] = m_A_0global[omega][0];
                        p_BCExp = 0;
			
                        // Daughter vessel 1
                        d1 = m_vessels[2*omega]->GetBndConditions()[0]->GetDaughter1();
                        nel_d1 = fields[d1][0].num_elements()-1; 
                        Au[1] = fields[d1][0][nel_d1];
                        uu[1] = fields[d1][1][nel_d1];
                        beta[1] = m_betaglobal[d1][0];
                        A_0[1] = m_A_0global[d1][0];
                        d1_BCExp = 1;
			
                        // Daughter vessel 2
                        d2 = m_vessels[2*omega]->GetBndConditions()[0]->GetDaughter2();
                        nel_d2 = fields[d2][0].num_elements()-1; 
                        Au[2] = fields[d2][0][nel_d2];
                        uu[2] = fields[d2][1][nel_d2];
                        beta[2] = m_betaglobal[d2][0];
                        A_0[2] = m_A_0global[d2][0];
                        d2_BCExp = 1;
			
			
                        // Solve the Riemann problem for a bifurcation
                        MergingRiemann2_to_1(Au, uu, beta, A_0);
			
			
                        // Store the values into the right positions:
                        // Parent vessel
                        m_vessels[2*omega]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = Au[0];
                        m_vessels[2*omega+1]->UpdateBndCondExpansion(p_BCExp)->UpdateCoeffs()[0] = uu[0];
                        
                        // Daughter vessel 1
                        m_vessels[2*d1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = Au[1];
                        m_vessels[2*d1+1]->UpdateBndCondExpansion(d1_BCExp)->UpdateCoeffs()[0] = uu[1];
                        
                        // Daughter vessel 2
                        m_vessels[2*d2]->UpdateBndCondExpansion(d2_BCExp)->UpdateCoeffs()[0] = Au[2];
                        m_vessels[2*d2+1]->UpdateBndCondExpansion(d2_BCExp)->UpdateCoeffs()[0] = uu[2];
                    }
                    
		}
                
    }
    
    
    /**
     *  Solves the Riemann problem at a bifurcation by assuming
     *  subsonic flow at both sides of the boundary and by
     *  applying conservation of mass and continuity of the total
     *  pressure \f$ \frac{p}{rho} + \frac{u^{2}}{2}. \f$ The
     *  other 3 missing equations come from the characteristic
     *  variables. For further information see "Pulse
     *  WavePropagation in the human vascular system" Section
     *  3.4.4
     */
    void PulseWaveSystem::BifurcationRiemann(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu, 
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
            
            // Calculated from input
            W[0] = uu[0] + 4*sqrt(beta[0]/(2*rho))*(sqrt(sqrt(Au[0])) - sqrt(sqrt(A_0[0])));
            W[1] = uu[1] - 4*sqrt(beta[1]/(2*rho))*(sqrt(sqrt(Au[1])) - sqrt(sqrt(A_0[1])));
            W[2] = uu[2] - 4*sqrt(beta[2]/(2*rho))*(sqrt(sqrt(Au[2])) - sqrt(sqrt(A_0[2])));
            
            // Tolerances for the algorithm
            NekDouble Tol = 1.0e-10;
            
            // Newton Iteration
            while ((proceed) && (iter < MAX_ITER))
            {
                iter = iter+1;
                
                // Calculate the constraint vector, six equations:
                // 3 characteristic variables, mass conservation, 
                // total pressure
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
                        ASSERTL0(false,"Riemann solver for Bifurcation did not converge");
                    }
		}
	}
	
	
	/**
	 *  Solves the Riemann problem at an interdomain junction by assuming subsonic flow at
	 *  both sides of the boundary and by applying conservation of mass and continuity of the
	 *  total pressure \f$ \frac{p}{rho} + \frac{u^{2}}{2}. \f$ The other 2 missing equations
	 *  come from the characteristic variables. For further information see "Pulse WavePropagation
	 *  in the human vascular system" Section 3.4.
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
			ASSERTL0(false,"Riemann solver for Junction did not converge");
        }
    }
    
#endif            

	
    /**
     *  Writes the .fld file at the end of the simulation. Similar to the normal
     *  v_Output however the Multidomain output has to be prepared.
     */
    void PulseWaveSystem::v_Output(void)
    {
        /**
         * Write the field data to file. The file is named according to the session
         * name with the extension .fld appended.
         */
        std::string outname = m_sessionName;
        if (m_comm->GetSize() > 1)
        {
            outname += "_P"+boost::lexical_cast<std::string>(m_comm->GetRank());
        }
        outname += ".fld";


        WriteVessels(outname);	
    }


    /**
     *  Writes the .fld file at the end of the simulation. Similar to the normal
     *  v_Output however the Multidomain output has to be prepared.
     */
    void PulseWaveSystem::CheckPoint_Output(const int n)
    {
        std::stringstream outname;
        outname << m_sessionName << "_" << n;
        
        if (m_comm->GetSize() > 1)
        {
            outname << "_P" << m_comm->GetRank();
        }
        outname << ".chk";

        WriteVessels(outname.str());	
    }


    /**
     * Writes the field data to a file with the given filename.
     * @param   outname         Filename to write to.
     */
    void PulseWaveSystem::WriteVessels(const std::string &outname)
    {

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
        std::vector<std::string> variables = m_session->GetVariables();

        for(int n = 0; n < m_nDomains; ++n)
        {
            m_vessels[n*m_nVariables]->GetFieldDefinitions(FieldDef);
        }
        
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
        
        int nFieldDefPerDomain = FieldDef.size()/m_nDomains;
        int cnt;
        // Copy Data into FieldData and set variable
        for(int n = 0; n < m_nDomains; ++n)
        {
            for(int j = 0; j < m_nVariables; ++j)
            {
                for(int i = 0; i < nFieldDefPerDomain; ++i)
                {
                    cnt = n*nFieldDefPerDomain+i;
                    FieldDef[cnt]->m_fields.push_back(variables[j]);
                    m_vessels[n*m_nVariables]->AppendFieldData(FieldDef[cnt], FieldData[cnt], m_vessels[n*m_nVariables+j]->UpdateCoeffs());
                }
            }            
        }
        
        // Update time in field info if required
        if(m_fieldMetaDataMap.find("Time") != m_fieldMetaDataMap.end())
        {
            m_fieldMetaDataMap["Time"] =  m_time; 
        }
        
        //LibUtilities::CombineFields(FieldDef, FieldData);
        
        LibUtilities::Write(outname, FieldDef, FieldData, m_fieldMetaDataMap);
    }    

    /* Compute the error in the L2-norm
     * @param   field           The field to compare.
     * @param   exactsoln       The exact solution to compare with.
     * @param   Normalised      Normalise L2-error.
     * @returns                 Error in the L2-norm.
     */
    NekDouble PulseWaveSystem::v_L2Error(unsigned int field,
                                         const Array<OneD, NekDouble> &exactsoln,
                                         bool Normalised)
    {    		
        NekDouble L2error = 0.0; 
        NekDouble L2error_dom;
        NekDouble Vol     = 0.0;

        if(m_NumQuadPointsError == 0)
        {
            for (int omega = 0; omega < m_nDomains; omega++)
            {
                int vesselid = field + omega*m_nVariables;
                
                if(m_vessels[vesselid]->GetPhysState() == false)
                {
                    m_vessels[vesselid]->BwdTrans(m_vessels[vesselid]->GetCoeffs(),
                                                  m_vessels[vesselid]->UpdatePhys());
                }
                
                if(exactsoln.num_elements())
                {
                    L2error_dom = m_vessels[vesselid]->L2(exactsoln);
                }
                else if (m_session->DefinesFunction("ExactSolution"))
                {
                    Array<OneD, NekDouble> exactsoln(m_vessels[vesselid]->GetNpoints());
                    
                    LibUtilities::EquationSharedPtr vEqu
                        = m_session->GetFunction("ExactSolution",field,omega);
                    EvaluateFunction(m_session->GetVariable(field),exactsoln,"ExactSolution",
                                     m_time);
                    
                    L2error_dom = m_vessels[vesselid]->L2(exactsoln);
                    
                }
                else
                {
                    L2error_dom = m_vessels[vesselid]->L2();
                }

                L2error += L2error_dom*L2error_dom;
                
                if(Normalised == true)
                {
                    Array<OneD, NekDouble> one(m_vessels[vesselid]->GetNpoints(), 1.0);
                    
                    Vol += m_vessels[vesselid]->PhysIntegral(one);
                }
                
            }
        }
        else
        {
            ASSERTL0(false,"Not set up");
        }
        
        
        if(Normalised == true)
        {
            m_comm->AllReduce(Vol, LibUtilities::ReduceSum);
        
            L2error = sqrt(L2error/Vol);
        }
        else
        {
            L2error = sqrt(L2error);
        }            

        return L2error;
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
        NekDouble LinferrorDom, Linferror = -1.0;		
        
        for (int omega = 0; omega < m_nDomains; omega++)
        {
            int vesselid = field + omega*m_nVariables;

            if(m_NumQuadPointsError == 0)
            {
                if(m_vessels[vesselid]->GetPhysState() == false)
                {
                    m_vessels[vesselid]->BwdTrans(m_vessels[vesselid]->GetCoeffs(),
                                                  m_vessels[vesselid]->UpdatePhys());
                }
		
                if(exactsoln.num_elements())
                {
                    LinferrorDom = m_vessels[vesselid]->Linf(exactsoln);
                }
                else if (m_session->DefinesFunction("ExactSolution"))
                {
                    Array<OneD, NekDouble> exactsoln(m_vessels[vesselid]->GetNpoints());
                    
                    EvaluateFunction(m_session->GetVariable(field),exactsoln,"ExactSolution",
                                     m_time);
                    
                    LinferrorDom = m_vessels[vesselid]->Linf(exactsoln);
                }
                else
                {
                    LinferrorDom = 0.0;
                }

                Linferror = (Linferror > LinferrorDom)? Linferror:LinferrorDom;
                
            }
            else
            {
                ASSERTL0(false,"ErrorExtraPoints not allowed for this solver");
            }
        }
        return Linferror;
    }
    
    void PulseWaveSystem::CalcCharacteristicVariables(int omega)
    {
        int nq = m_vessels[omega]->GetTotPoints();
        Array<OneD, NekDouble> A = m_vessels[omega]->UpdatePhys();
        Array<OneD, NekDouble> u = m_vessels[omega+1]->UpdatePhys();
        Array<OneD, NekDouble> c(nq);

        //Calc 4*c
        Vmath::Vsqrt(nq,A,1,c,1);
        Vmath::Vmul(nq,m_beta[omega],1,c,1,c,1);
        Vmath::Smul(nq,16.0/(2*m_rho),c,1,c,1);
        Vmath::Vsqrt(nq,c,1,c,1);

        // Characteristics
        Vmath::Vadd(nq,u,1,c,1,A,1);
        Vmath::Vsub(nq,u,1,c,1,u,1);
    }



}
