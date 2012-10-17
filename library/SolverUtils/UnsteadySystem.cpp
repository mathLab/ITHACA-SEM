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
// Description: Generic timestepping for Unsteady solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/UnsteadySystem.h>

namespace Nektar
{	
    namespace SolverUtils
    {
        NekDouble IntegrationTime = 0.0;
        /**
         * @class UnsteadySystem
         *
         * Provides the underlying timestepping framework for unsteady solvers
         * including the general timestepping routines. This class is not intended
         * to be directly instantiated, but rather is a base class on which to
         * define unsteady solvers.
         *
         * For details on implementing unsteady solvers see
         * \ref sectionADRSolverModuleImplementation here
         */

        /**
         * Processes SolverInfo parameters from the session file and sets up
         * timestepping-specific code.
         * @param   pSession        Session object to read parameters from.
         */
        UnsteadySystem::UnsteadySystem(
            const LibUtilities::SessionReaderSharedPtr& pSession)
            : EquationSystem(pSession),
              m_infosteps(10)

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
            m_session->MatchSolverInfo("REACTIONADVANCEMENT", "Explicit",
                                       m_explicitReaction, true);

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

            // Load generic input parameters
            m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
            m_session->LoadParameter("CFL", m_cfl, 0.0);

            // Set up filters
            LibUtilities::FilterMap::const_iterator x;
            LibUtilities::FilterMap f = m_session->GetFilters();
            for (x = f.begin(); x != f.end(); ++x)
            {
                m_filters.push_back(GetFilterFactory().CreateInstance(x->first, m_session, x->second));
            }


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
            int i,n,nchk = 1;
            int ncoeffs = m_fields[0]->GetNcoeffs();
            int npoints = m_fields[0]->GetNpoints();
            int nvariables = 0;

            if (m_intVariables.empty())
            {
                for (i = 0; i < m_fields.num_elements(); ++i)
                {
                    m_intVariables.push_back(i);
                }
                nvariables = m_fields.num_elements();
            }
            else
            {
                nvariables = m_intVariables.size();
            }

            // Set up wrapper to fields data storage.
            Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
            Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);
		
            for(i = 0; i < nvariables; ++i)
            {
                fields[i]  = m_fields[m_intVariables[i]]->UpdatePhys();
                m_fields[m_intVariables[i]]->SetPhysState(false);
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
                    LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eForwardEuler);

                    // Used for all other time steps
                    LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
                    IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                    IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                    
                    // Initialise the scheme for the actual time integration scheme
                    u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                    break;
                }
            case LibUtilities::eIMEXOrder2:
                {
                    numMultiSteps = 2;
                    
                    IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                    // Used in the first time step to initalize the scheme
                    LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eIMEXOrder1);
                    
                    // Used for all other time steps
                    LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
                    IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                    IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                    
                    // Initialise the scheme for the actual time integration scheme
                    u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                    break;
                }    
            case LibUtilities::eIMEXOrder3:
                {
                    numMultiSteps = 3;
                    
                    IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                    // Used in the first time step to initalize the scheme
                    LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eIMEXdirk_3_4_3);
                    LibUtilities::TimeIntegrationSchemeKey IntKey1(LibUtilities::eIMEXdirk_3_4_3);

                    // Used for all other time steps
                    LibUtilities::TimeIntegrationSchemeKey IntKey2(m_timeIntMethod);
                    IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                    IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                    IntScheme[2] = LibUtilities::TimeIntegrationSchemeManager()[IntKey2];

                    // Initialise the scheme for the actual time integration scheme
                    u = IntScheme[2]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                    break;
                }
                default:
                {
                    ASSERTL0(false,"populate switch statement for integration scheme");
                }
            }

            std::vector<FilterSharedPtr>::iterator x;
            for (x = m_filters.begin(); x != m_filters.end(); ++x)
            {
                (*x)->Initialise(m_fields, m_time);
            }
            
            //========================================================================
            // TIME LOOP =============================================================
            //========================================================================
            
            //========================================================================
            // WITH CFL CONTROLL CFL control has been implemented so
            // far just for profiling reasons (just for CFLTester) It
            // has to be generalised
            // =======================================================================
            
            if(m_cfl>0.0)
            {
                const Array<OneD,int> ExpOrder = GetNumExpModesPerExp();
		
                NekDouble TimeStability;
				
                switch(m_timeIntMethod)
                {
                case LibUtilities::eForwardEuler:
                case LibUtilities::eClassicalRungeKutta4:
                    {
                        TimeStability = 2.784;
                        break;
                    }
                case LibUtilities::eAdamsBashforthOrder1:
                case LibUtilities::eRungeKutta2_ModifiedEuler:
                case LibUtilities::eRungeKutta2_ImprovedEuler:
                    {
                        TimeStability = 2.0;
                        break;
                    }
                case LibUtilities::eAdamsBashforthOrder2:
                    {
                        TimeStability = 1.0;
                        break;
                    }
                default:
                    {
                        ASSERTL0(false,"No CFL control implementation for this time integration scheme");
                    }
                } 
		
                NekDouble CheckpointTime = 0.0;
                NekDouble QuarterOfLoop  = 0.25;
                NekDouble CFLtimestep;
		
                int number_of_checkpoints = ceil(m_fintime/QuarterOfLoop)+1;
			
                Array<OneD, NekDouble>   L2errors(number_of_checkpoints);
                Array<OneD, NekDouble>   LIerrors(number_of_checkpoints);
                Array<OneD, NekDouble>   TimeLevels(number_of_checkpoints);
		
                int checkpoints_cnt = 0;
		
                L2errors[checkpoints_cnt]         = L2Error(0);
                LIerrors[checkpoints_cnt]         = LinfError(0);
                TimeLevels[checkpoints_cnt]       = m_time;
                
                CheckpointTime += QuarterOfLoop;
                checkpoints_cnt++;
		
                // This function is implemented for the CFLTester only
                m_timestep = GetTimeStep(ExpOrder[0],m_cfl,TimeStability);
		
                CFLtimestep = m_timestep;
		
                Timer timer;
			
                int step = 0;
			
                while(m_time < m_fintime)
                {
                    timer.Start();
                    
                    // Increment time.
                    if(m_time + m_timestep > m_fintime && m_fintime > 0.0)
                    {
                        m_timestep = m_fintime - m_time;
                    }

                    // Integrate over timestep.
                    if( step < numMultiSteps-1)
                    {
                        // Use initialisation schemes if time step is less than the
                        // number of steps in the scheme.
                        fields = IntScheme[step]->TimeIntegrate(m_timestep,u,m_ode);
                    }
                    else
                    {
                        fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,m_ode);
                    }
                    
                    
                    m_time += m_timestep;
				
                    timer.Stop();
				
                    IntegrationTime += timer.TimePerTest(1);
				
                    if(m_time <= CheckpointTime && (m_time + m_timestep) > CheckpointTime)
                    {
                        m_fields[0]->FwdTrans(fields[0],m_fields[0]->UpdateCoeffs());
                        m_fields[0]->UpdatePhys() = fields[0];
                        L2errors[checkpoints_cnt] = L2Error(0);
                        LIerrors[checkpoints_cnt] = LinfError(0);
                        TimeLevels[checkpoints_cnt] = CheckpointTime;
                        CheckpointTime += QuarterOfLoop;
                        checkpoints_cnt++;
                    }				
                    // step advance
                    step++;
                }
		
                // At the end of the time integration, store final solution.
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->UpdatePhys() = fields[i];
                }
		
                cout <<"\nCFL number       : " << m_cfl                     << endl;
                cout <<"\nCFL time-step    : " << CFLtimestep               << endl;
                cout <<"\nTime-integration : " << IntegrationTime << " s"   << endl;
                for(int i = 0; i < number_of_checkpoints; i++)
                {
                    cout <<"Time : "<< TimeLevels[i] << "\tL2 error : " << L2errors[i] << "\tLI error : "  << LIerrors[i] << endl;
                }
            }
            
            //======================================================================
            //WITHOUT CFL CONTROL
            //======================================================================
            else
            {
                for(n = 0; n < m_steps; ++n)
                {
                    Timer timer;
                    timer.Start();
                    
                    /// Integrate over timestep
                    if( n < numMultiSteps-1)
                    {
                        /// Use initialisation schemes if time step is 
                        /// less than the number of steps in the scheme
                        fields = IntScheme[n]->TimeIntegrate(m_timestep,u,m_ode);
                    }
                    else
                    {
                        fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,m_ode);
                    }
				
                    m_time += m_timestep;
				
                    timer.Stop();
                    IntegrationTime += timer.TimePerTest(1);
				
                    /// Write out status information
                    if(m_session->GetComm()->GetRank() == 0
                       && !((n+1)%m_infosteps))
                    {
                        cout << "Steps: "           << n+1
                             << "\t Time: "         << m_time
                             << "\t Time-step: "    << m_timestep << "\t" << endl;
                    }

                    /// Transform data into coefficient space
                    for (i = 0; i < m_intVariables.size(); ++i)
                    {
                        m_fields[m_intVariables[i]]->FwdTrans_IterPerExp(fields[i],
                                                                         m_fields[m_intVariables[i]]->UpdateCoeffs());
                        //###### Vmath::Vcopy(npoints, fields[i], 1, m_fields[m_intVariables[i]]->UpdatePhys(), 1);
                        m_fields[m_intVariables[i]]->SetPhysState(false);
                    }

                    std::vector<FilterSharedPtr>::iterator x;
                    for (x = m_filters.begin(); x != m_filters.end(); ++x)
                    {
                        (*x)->Update(m_fields, m_time);
                    }

                    /// Write out checkpoint files.
                    if(m_checksteps&&n&&(!((n+1)%m_checksteps)))
                    {
                        Checkpoint_Output(nchk++);
                    }
                    /// Step advance
                }
			
                cout <<"\nCFL number              : " << m_cfl << endl;
                cout <<"Time-integration timing : " << IntegrationTime << " s" << endl << endl;

                // At the end of the time integration, store final solution.
                for(i = 0; i < m_intVariables.size(); ++i)
                {
	            if(m_fields[m_intVariables[i]]->GetPhysState() == false)
	            {
	                m_fields[m_intVariables[i]]->BwdTrans(m_fields[m_intVariables[i]]->GetCoeffs(),m_fields[m_intVariables[i]]->UpdatePhys());
	            }
                    m_fields[m_intVariables[i]]->UpdatePhys() = fields[i];
                }

                std::vector<FilterSharedPtr>::iterator x;
                for (x = m_filters.begin(); x != m_filters.end(); ++x)
                {
                    (*x)->Finalise(m_fields, m_time);
                }
            }
            //===========================================================================================
            // END OF TIME LOOP =========================================================================
            //===========================================================================================
            
            /// Print for 1D problems
            if(m_spacedim == 1)
            {
                v_AppendOutput1D(fields);   
            }
        }
        
        

        void UnsteadySystem::v_DoInitialise()
        {
            SetInitialConditions();
        }
        

        
        void UnsteadySystem::v_PrintSummary(std::ostream &out)
        {
            EquationSystem::v_PrintSummary(out);
            out << "\tAdvection       : " << (m_explicitAdvection ? "explicit" : "implicit") << endl;
            out << "\tDiffusion       : " << (m_explicitDiffusion ? "explicit" : "implicit") << endl;
            if (m_session->GetSolverInfo("EQTYPE")== "SteadyAdvectionDiffusionReaction")
            {
                out << "\tReaction        : " << (m_explicitReaction  ? "explicit" : "implicit") << endl;
            }
            out << "\tIntegration Type: " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod]<< endl;
            out << "\tTime Step       : " << m_timestep                                             << endl;
            out << "\tNo. of Steps    : " << m_steps                                                << endl;
            out << "\tCheckpoints     : " << m_checksteps << " steps"                               << endl;
        }
        
        
        
        void UnsteadySystem::v_AppendOutput1D(Array<OneD, Array<OneD, NekDouble> > &solution1D)
        {
            /// Coordinates of the quadrature points in the real physical space
            Array<OneD,NekDouble> x(GetNpoints());
            Array<OneD,NekDouble> y(GetNpoints());
            Array<OneD,NekDouble> z(GetNpoints());
            m_fields[0]->GetCoords(x, y, z);
            
            /// Print out the solution in a txt file
            ofstream outfile;
            outfile.open("solution1D.txt");
            for(int i = 0; i < GetNpoints(); i++)
            {
                outfile << scientific 
                << setw (17) 
                << setprecision(16) 
                << x[i]
                << "  " 
                << solution1D[0][i] 
                << endl;
            }
            outfile << endl << endl;
            outfile.close();
        }



        void UnsteadySystem::v_NumericalFlux(
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            ASSERTL0(false, "This function is not implemented for this equation.");
        }



        void UnsteadySystem::v_NumericalFlux(
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &numfluxX,
            Array<OneD, Array<OneD, NekDouble> > &numfluxY )
        {
            ASSERTL0(false, "This function is not implemented for this equation.");
        }



        void UnsteadySystem::v_NumFluxforScalar(
            Array<OneD, Array<OneD, NekDouble>  > &ufield,
            Array<OneD, Array<OneD, Array<OneD,NekDouble> > > &uflux)
        {
            int i,j;
            int nTraceNumPoints = GetTraceNpoints();
            int nvariables = m_fields.num_elements();
            int nqvar = uflux.num_elements();

            Array<OneD, NekDouble > Fwd(nTraceNumPoints);
            Array<OneD, NekDouble > Bwd(nTraceNumPoints);
            Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
            Array<OneD, NekDouble > fluxtemp (nTraceNumPoints,0.0);

            // Get the sign of (v \cdot n), v = an arbitrary vector

            //  Evaulate upwind flux of uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
            for (j = 0; j < nqvar; ++j)
            {
                for(i = 0; i < nvariables ; ++i)
                {
                    //  Compute Forward and Backward value of ufield of i direction

                    m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);

                    // if Vn >= 0, flux = uFwd, i.e.,
                    //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                    //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd

                    // else if Vn < 0, flux = uBwd, i.e.,
                    //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                    //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

                    m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],Fwd,Bwd,fluxtemp);

                    // Imposing weak boundary condition with flux
                    // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                    //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                    //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd

                    // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                    //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                    //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd

                    if(m_fields[0]->GetBndCondExpansions().num_elements())
                    {
                        WeakPenaltyforScalar(i,ufield[i],fluxtemp);
                    }

                    // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                    // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                    // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

                    // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n} ), i. e,
                    // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                    // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

                    Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,fluxtemp,1,uflux[j][i],1);

                }
            }
        }



        void UnsteadySystem::v_NumFluxforVector(
            Array<OneD, Array<OneD, NekDouble> > &ufield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            Array<OneD, Array<OneD, NekDouble> > &qflux)
        {
            int nTraceNumPoints = GetTraceNpoints();
            int nvariables = m_fields.num_elements();
            int nqvar = qfield.num_elements();

            NekDouble C11 = 1.0;
            Array<OneD, NekDouble > Fwd(nTraceNumPoints);
            Array<OneD, NekDouble > Bwd(nTraceNumPoints);
            Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

            Array<OneD, NekDouble > qFwd(nTraceNumPoints);
            Array<OneD, NekDouble > qBwd(nTraceNumPoints);
            Array<OneD, NekDouble > qfluxtemp(nTraceNumPoints,0.0);

            Array<OneD, NekDouble > uterm(nTraceNumPoints);

            // Evaulate upwind flux of qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
            for(int i = 0; i < nvariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
                for(int j = 0; j < nqvar; ++j)
                {
                    //  Compute Forward and Backward value of ufield of jth         direction
                    m_fields[i]->GetFwdBwdTracePhys(qfield[j][i],qFwd,qBwd);

                    // if Vn >= 0, flux = uFwd, i.e.,
                    //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick qflux = qBwd = q+
                    //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick qflux = qBwd = q-

                    // else if Vn < 0, flux = uBwd, i.e.,
                    //  edge::eForward, if V*n<0 <=> V*n_F<0, pick qflux = qFwd = q-
                    //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick qflux = qFwd   =q+

                    m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],qBwd,qFwd,qfluxtemp);
                    Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,qfluxtemp,1,qfluxtemp,1);

                    // Generate Stability term = - C11 ( u- - u+ )
                    m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);
                    Vmath::Vsub(nTraceNumPoints,Fwd,1,Bwd,1,uterm,1);
                    Vmath::Smul(nTraceNumPoints,-1.0*C11,uterm,1,uterm,1);

                    //  Flux = {Fwd,Bwd}*(nx,ny,nz) + uterm*(nx,ny)
                    Vmath::Vadd(nTraceNumPoints,uterm,1,qfluxtemp,1,qfluxtemp,1);

                    // Imposing weak boundary condition with flux
                    if(m_fields[0]->GetBndCondExpansions().num_elements())
                    {
                        WeakPenaltyforVector(i,j,qfield[j][i],qfluxtemp,C11);
                    }

                    // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                    // n_xi = n_x*tan_xi_x + n_y*tan_xi_y + n_z*tan_xi_z
                    // n_xi = n_x*tan_eta_x + n_y*tan_eta_y + n_z*tan_eta_z
                    Vmath::Vadd(nTraceNumPoints,qfluxtemp,1,qflux[i],1,qflux[i],1);
                }
            }
        }



        void UnsteadySystem::v_GetFluxVector(const int i, const int j,
                                             Array<OneD, Array<OneD, NekDouble> > &physfield,
                                             Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            for(int k = 0; k < flux.num_elements(); ++k)
            {
                Vmath::Zero(GetNpoints(),flux[k],1);
            }
            Vmath::Vcopy(GetNpoints(),physfield[i],1,flux[j],1);
        }



        void UnsteadySystem::WeakPenaltyforScalar(const int var,
                                                  const Array<OneD, const NekDouble> &physfield,
                                                  Array<OneD, NekDouble> &penaltyflux,
                                                  NekDouble time)
        {
            int i, j, e, npoints, id1, id2;
            // Number of boundary regions
            int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
            int Nfps, numBDEdge;
            int nTraceNumPoints = GetTraceNpoints();
            int cnt = 0;

            Array<OneD, NekDouble > uplus(nTraceNumPoints);

            m_fields[var]->ExtractTracePhys(physfield,uplus);
            for(i = 0; i < nbnd; ++i)
            {
                // Number of boundary expansion related to that region
                numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();
                // Evaluate boundary values g_D or g_N from input files
                LibUtilities::EquationSharedPtr ifunc =  m_session->GetFunction("InitialConditions", 0);
                npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();
                Array<OneD,NekDouble> BDphysics(npoints);
                Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);

                m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
                ifunc->Evaluate(x0,x1,x2,time,BDphysics);

                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < numBDEdge ; ++e)
                {
                    // Number of points on the expansion
                    Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0) ;
                    id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

                    // For Dirichlet boundary condition: uflux = g_D
                    if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vcopy(Nfps,&BDphysics[id1],1,&penaltyflux[id2],1);
                    }
                    // For Neumann boundary condition: uflux = u+
                    else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(Nfps,&uplus[id2],1,&penaltyflux[id2],1);
                    }
                }
            }
        }


        /**
         * Diffusion: Imposing weak boundary condition for q with flux
         *  uflux = g_D  on Dirichlet boundary condition
         *  uflux = u_Fwd  on Neumann boundary condition
         */
        void UnsteadySystem::WeakPenaltyforVector(
            const int var,
            const int dir,
            const Array<OneD, const NekDouble> &physfield,
            Array<OneD, NekDouble> &penaltyflux,
            NekDouble C11,
            NekDouble time)
        {
            int i, j, e, npoints, id1, id2;
            int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
            int numBDEdge, Nfps;
            int nTraceNumPoints = GetTraceNpoints();
            Array<OneD, NekDouble > uterm(nTraceNumPoints);
            Array<OneD, NekDouble > qtemp(nTraceNumPoints);
            int cnt = 0;

            m_fields[var]->ExtractTracePhys(physfield,qtemp);

            for(i = 0; i < nbnd; ++i)
            {
                numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();
                // Evaluate boundary values g_D or g_N from input files
                LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction("InitialConditions", 0);
                npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();

                Array<OneD,NekDouble> BDphysics(npoints);
                Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);

                m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
                ifunc->Evaluate(x0,x1,x2,time,BDphysics);

                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < numBDEdge ; ++e)
                {
                    Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);

                    id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

                    // For Dirichlet boundary condition: qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
                    if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&qtemp[id2],1, &penaltyflux[id2],1);
                    }
                    // For Neumann boundary condition: qflux = g_N
                    else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&BDphysics[id1],1,&penaltyflux[id2],1);
                    }
                }
            }
        }
	
	/**
	 *  This function calculate the proper time-step to keep the
	 *  problem stable.  It has been implemented to deal with an
	 *  explict treatment of the advection term.  In case of an
	 *  explicit treatment of the diffusion term a
	 *  re-implementation is required.  The actual implementation
	 *  can be found inside each equation class.
	 *
	 * @param ExpOrder          the expansion order we are using (P)
	 * @param CFL               the CFL number we want to impose (<1)
	 @ @param timeCFL           the stability coefficient, a combination of the spatial/temporal discretisation stability region
        */
	NekDouble UnsteadySystem::GetTimeStep(const Array<OneD,int> ExpOrder, 
                                              const Array<OneD,NekDouble> CFL, 
                                              NekDouble timeCFL)
	{
            NekDouble TimeStep = v_GetTimeStep(ExpOrder, CFL, timeCFL);
		
            return TimeStep;
	}
	
	/**
	*  This function calculate the proper time-step to keep the problem stable.
	*  It has been implemented to deal with an explict treatment of the advection term.
	*  In case of an explicit treatment of the diffusion term a re-implementation is required.
	*  The actual implementation can be found inside each equation class.
	*
	* @param ExpOrder          the expansion order we are using (P)
	* @param CFL               the CFL number we want to impose (<1)
	@ @param timeCFL           the stability coefficient of the time-integration scheme
	*/
	NekDouble UnsteadySystem::GetTimeStep(int ExpOrder, NekDouble CFL, NekDouble TimeStability)
	{
		NekDouble TimeStep = v_GetTimeStep(ExpOrder,CFL,TimeStability);
		
		return TimeStep;
	}
	
	/**
	 * See GetTimeStep. 
	 * This is the virtual fuction to redirect the implementation to the proper class.
	 */
	NekDouble UnsteadySystem::v_GetTimeStep(const Array<OneD,int> ExpOrder, 
                                                const Array<OneD,NekDouble> CFL, NekDouble timeCFL)
        {
            ASSERTL0(false, "v_GetTimeStep is not implemented in the base class (UnsteadySystem). Check if your equation class has its own implementation");
            
            return 0.0;
        }
	
	/**
	 * See GetTimeStep. 
	 * This is the virtual fuction to redirect the implementation to the proper class.
	 */
	NekDouble UnsteadySystem::v_GetTimeStep(int ExpOrder, NekDouble CFL, NekDouble TimeStability)
	{
		ASSERTL0(false, "v_GetTimeStep is not implemented in the base class (UnsteadySystem). Check if your equation class has its own implementation");
                
		return 0.0;
	}
	
        
	Array<OneD,NekDouble> UnsteadySystem::GetStdVelocity(const Array<OneD, Array<OneD,NekDouble> > inarray)
	{
            // Checking if the problem is 2D
            ASSERTL0(m_expdim>=2,"Method not implemented for 1D");
	
            int nTotQuadPoints  = GetTotPoints();
            int n_element  = m_fields[0]->GetExpSize();       // number of element in the mesh
            int nvel = inarray.num_elements();
            int npts = 0;
            
            NekDouble pntVelocity;
            
            // Getting the standard velocity vector on the 2D normal space
            Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
		
            Array<OneD, NekDouble> stdV(n_element,0.0);
            for (int i = 0; i < nvel; ++i)
            {
                stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints);
            }
		
            if(nvel == 2)
            {
                for(int el = 0; el < n_element; ++el)
                { 
                    
                    int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                    
                    Array<OneD, const NekDouble> jac  = m_fields[0]->GetExp(el)->GetGeom2D()->GetJac();
                    Array<TwoD, const NekDouble> gmat = m_fields[0]->GetExp(el)->GetGeom2D()->GetGmat();

                    if(m_fields[0]->GetExp(el)->GetGeom2D()->GetGtype() == SpatialDomains::eDeformed)
                    {
                        for(int i=0; i<n_points; i++)
                        {
                            stdVelocity[0][i] = gmat[0][i]*inarray[0][i] + gmat[2][i]*inarray[1][i];
                            stdVelocity[1][i] = gmat[1][i]*inarray[0][i] + gmat[3][i]*inarray[1][i];
                        }
                    }
                    else
                    {
                        for(int i=0; i<n_points; i++)
                        {
                            stdVelocity[0][i] = gmat[0][0]*inarray[0][i] + gmat[2][0]*inarray[1][i];
                            stdVelocity[1][i] = gmat[1][0]*inarray[0][i] + gmat[3][0]*inarray[1][i];
                        }
                    }
                    

                    for(int i=0; i<n_points; i++)
                    {
                        pntVelocity = sqrt(stdVelocity[0][i]*stdVelocity[0][i] + stdVelocity[1][i]*stdVelocity[1][i]);
                        if(pntVelocity>stdV[el])
                        {
                            stdV[el] = pntVelocity;
                        }
                        
                    }
                }
            }
            else
            {
                for(int el = 0; el < n_element; ++el)
                { 
                    
                    int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                    
                    Array<OneD, const NekDouble> jac  = m_fields[0]->GetExp(el)->GetGeom3D()->GetJac();
                    Array<TwoD, const NekDouble> gmat = m_fields[0]->GetExp(el)->GetGeom3D()->GetGmat();

                    if(m_fields[0]->GetExp(el)->GetGeom3D()->GetGtype() == SpatialDomains::eDeformed)
                    {
                        for(int i=0; i<n_points; i++)
                        {
                            stdVelocity[0][i] = gmat[0][i]*inarray[0][i] + gmat[3][i]*inarray[1][i] + gmat[6][i]*inarray[2][i];
                            stdVelocity[1][i] = gmat[1][i]*inarray[0][i] + gmat[4][i]*inarray[1][i] + gmat[7][i]*inarray[2][i];
                            stdVelocity[2][i] = gmat[2][i]*inarray[0][i] + gmat[5][i]*inarray[1][i] + gmat[8][i]*inarray[2][i];
                        }
                    }
                    else
                    {
                        Array<OneD, const NekDouble> jac  = m_fields[0]->GetExp(el)->GetGeom3D()->GetJac();
                        Array<TwoD, const NekDouble> gmat = m_fields[0]->GetExp(el)->GetGeom3D()->GetGmat();

                        for(int i=0; i<n_points; i++)
                        {
                            stdVelocity[0][i] = gmat[0][0]*inarray[0][i] + gmat[3][0]*inarray[1][i] + gmat[6][0]*inarray[2][i];
                            stdVelocity[1][i] = gmat[1][0]*inarray[0][i] + gmat[4][0]*inarray[1][i] + gmat[7][0]*inarray[2][i];
                            stdVelocity[2][i] = gmat[2][0]*inarray[0][i] + gmat[5][0]*inarray[1][i] + gmat[8][0]*inarray[2][i];
                        }
                    }
                    
                    for(int i=0; i<n_points; i++)
                    {
                        pntVelocity = sqrt(stdVelocity[0][i]*stdVelocity[0][i] + stdVelocity[1][i]*stdVelocity[1][i] + stdVelocity[2][i]*stdVelocity[2][i]);
                        if(pntVelocity>stdV[el])
                        {
                            stdV[el] = pntVelocity;
                        }
                    }
                }
            }
		
            return stdV;
	}
	
	NekDouble UnsteadySystem::GetStabilityLimit(int n)
	{
            if(n>20)
	    {
                ASSERTL0(false,"Illegal modes dimension for CFL calculation (P has to be less then 20)");
            }
		
            NekDouble CFLDG[21] = {2,6,11.8424,19.1569,27.8419,37.8247,49.0518,61.4815,75.0797,89.8181,105.67,122.62,140.64,159.73,179.85,201.01,223.18,246.36,270.53,295.69,321.83}; //CFLDG 1D [0-20]
            NekDouble CFLCG[2]  = {1.0,1.0};
            NekDouble CFL;
		
            if (m_projectionType == MultiRegions::eDiscontinuous)
            {
                CFL = CFLDG[n];
            }
            else 
            {
                ASSERTL0(false,"Continuos Galerkin stability coefficients not introduced yet.");
            }
		
            return CFL;
	}
	
	Array<OneD,NekDouble> UnsteadySystem::GetStabilityLimitVector(const Array<OneD,int> &ExpOrder)
	{
            int i;
            Array<OneD,NekDouble> returnval(m_fields[0]->GetExpSize(),0.0);
            for(i =0; i<m_fields[0]->GetExpSize(); i++)
            {
                returnval[i] = GetStabilityLimit(ExpOrder[i]);
            }
            return returnval;
	}
    }
}
