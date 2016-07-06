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

#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/UnsteadySystem.h>

using namespace std;

namespace Nektar
{	
    namespace SolverUtils
    {
        /**
         * @class UnsteadySystem
         *
         * Provides the underlying timestepping framework for unsteady solvers
         * including the general timestepping routines. This class is not 
         * intended to be directly instantiated, but rather is a base class 
         * on which to define unsteady solvers.
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

        /**
         * Initialization object for UnsteadySystem class.
         */
        void UnsteadySystem::v_InitObject()
        {
            EquationSystem::v_InitObject();
            
            m_initialStep = 0;

            // Load SolverInfo parameters
            m_session->MatchSolverInfo("DIFFUSIONADVANCEMENT","Explicit",
                                       m_explicitDiffusion,true);
            m_session->MatchSolverInfo("ADVECTIONADVANCEMENT","Explicit",
                                       m_explicitAdvection,true);
            m_session->MatchSolverInfo("REACTIONADVANCEMENT", "Explicit",
                                       m_explicitReaction, true);

            // For steady problems, we do not initialise the time integration
            if (m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"))
            {
                m_intScheme = LibUtilities::GetTimeIntegrationWrapperFactory().
                    CreateInstance(m_session->GetSolverInfo(
                                       "TIMEINTEGRATIONMETHOD"));

                // Load generic input parameters
                m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
                m_session->LoadParameter("CFL", m_cflSafetyFactor, 0.0);

                // Set up time to be dumped in field information
                m_fieldMetaDataMap["Time"] =
                        boost::lexical_cast<std::string>(m_time);
            }

            // By default attempt to forward transform initial condition.
            m_homoInitialFwd = true;

            // Set up filters
            LibUtilities::FilterMap::const_iterator x;
            LibUtilities::FilterMap f = m_session->GetFilters();
            for (x = f.begin(); x != f.end(); ++x)
            {
                m_filters.push_back(GetFilterFactory().CreateInstance(
                                                                x->first, 
                                                                m_session, 
                                                                x->second));
            }
        }
        
        /**
         * Destructor for the class UnsteadyAdvection.
         */
        UnsteadySystem::~UnsteadySystem()
        {
        }

        /**
         * @brief Returns the maximum time estimator for CFL control.
         */
        NekDouble UnsteadySystem::MaxTimeStepEstimator()
        {
            NekDouble TimeStability = 0.0;
            switch(m_intScheme->GetIntegrationMethod())
            {
                case LibUtilities::eForwardEuler:
                case LibUtilities::eClassicalRungeKutta4:
                case LibUtilities::eRungeKutta4:
                {
                    TimeStability = 2.784;
                    break;
                }
                case LibUtilities::eAdamsBashforthOrder1:
                case LibUtilities::eMidpoint:
                case LibUtilities::eRungeKutta2:
                case LibUtilities::eRungeKutta2_ImprovedEuler:
                case LibUtilities::eRungeKutta2_SSP:
                case LibUtilities::eRungeKutta3_SSP:
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
                    ASSERTL0(
                        false,
                        "No CFL control implementation for this time"
                        "integration scheme");
                }
            }
            return TimeStability;
        }
        
        /**
         * @brief Initialises the time integration scheme (as specified in the 
         * session file), and perform the time integration.
         */
        void UnsteadySystem::v_DoSolve()
        {
            ASSERTL0(m_intScheme != 0, "No time integration scheme.");

            int i = 1;
            int nvariables = 0;
            int nfields = m_fields.num_elements();

            if (m_intVariables.empty())
            {
                for (i = 0; i < nfields; ++i)
                {
                    m_intVariables.push_back(i);
                }
                nvariables = nfields;
            }
            else
            {
                nvariables = m_intVariables.size();
            }

            // Integrate in wave-space if using homogeneous1D
            if(m_HomogeneousType == eHomogeneous1D && m_homoInitialFwd)
            {
                for(i = 0; i < nfields; ++i)
                {
                    m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),
                                                     m_fields[i]->UpdatePhys());
                    m_fields[i]->SetWaveSpace(true);
                    m_fields[i]->SetPhysState(false);
                }
            }

            // Set up wrapper to fields data storage.
            Array<OneD, Array<OneD, NekDouble> > fields(nvariables);
            Array<OneD, Array<OneD, NekDouble> > tmp   (nvariables);
            
            // Order storage to list time-integrated fields first.
            for(i = 0; i < nvariables; ++i)
            {
                fields[i] = m_fields[m_intVariables[i]]->GetPhys();
                m_fields[m_intVariables[i]]->SetPhysState(false);
            }
            
            // Initialise time integration scheme
            m_intSoln = m_intScheme->InitializeScheme(
                m_timestep, fields, m_time, m_ode);

            // Initialise filters
            std::vector<FilterSharedPtr>::iterator x;
            for (x = m_filters.begin(); x != m_filters.end(); ++x)
            {
                (*x)->Initialise(m_fields, m_time);
            }

            // Ensure that there is no conflict of parameters
            if(m_cflSafetyFactor > 0.0)
            {
                // Check final condition
                ASSERTL0(m_fintime == 0.0 || m_steps == 0,
                         "Final condition not unique: "
                         "fintime > 0.0 and Nsteps > 0");
                
                // Check timestep condition
                ASSERTL0(m_timestep == 0.0, 
                         "Timestep not unique: timestep > 0.0 & CFL > 0.0");
            }

            // Check uniqueness of checkpoint output
            ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                     (m_checktime >  0.0 && m_checksteps == 0) || 
                     (m_checktime == 0.0 && m_checksteps >  0),
                     "Only one of IO_CheckTime and IO_CheckSteps "
                     "should be set!");

            Timer     timer;
            bool      doCheckTime   = false;
            int       step          = m_initialStep;
            NekDouble intTime       = 0.0;
            NekDouble lastCheckTime = 0.0;
            NekDouble cpuTime       = 0.0;
            NekDouble elapsed       = 0.0;

            while (step   < m_steps ||
                   m_time < m_fintime - NekConstants::kNekZeroTol)
            {
                if (m_cflSafetyFactor)
                {
                    m_timestep = GetTimeStep(fields);
                    
                    // Ensure that the final timestep finishes at the final
                    // time, or at a prescribed IO_CheckTime.
                    if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
                    {
                        m_timestep = m_fintime - m_time;
                    }
                    else if (m_checktime && 
                             m_time + m_timestep - lastCheckTime >= m_checktime)
                    {
                        lastCheckTime += m_checktime;
                        m_timestep     = lastCheckTime - m_time;
                        doCheckTime    = true;
                    }
                }
                
                // Perform any solver-specific pre-integration steps
                timer.Start();
                if (v_PreIntegrate(step))
                {
                    break;
                }

                fields = m_intScheme->TimeIntegrate(
                    step, m_timestep, m_intSoln, m_ode);
                timer.Stop();

                m_time  += m_timestep;
                elapsed  = timer.TimePerTest(1);
                intTime += elapsed;
                cpuTime += elapsed;
		
                // Write out status information
                if (m_session->GetComm()->GetRank() == 0 && 
                    !((step+1) % m_infosteps))
                {
                    cout << "Steps: " << setw(8)  << left << step+1 << " "
                         << "Time: "  << setw(12) << left << m_time;

                    if (m_cflSafetyFactor)
                    {
                        cout << " Time-step: " << setw(12)
                             << left << m_timestep;
                    }

                    stringstream ss;
                    ss << cpuTime << "s";
                    cout << " CPU Time: " << setw(8) << left
                         << ss.str() << endl;
                    cpuTime = 0.0;
                }

                // Transform data into coefficient space
                for (i = 0; i < nvariables; ++i)
                {
                    m_fields[m_intVariables[i]]->SetPhys(fields[i]);
                    m_fields[m_intVariables[i]]->FwdTrans_IterPerExp(
                        fields[i],
                        m_fields[m_intVariables[i]]->UpdateCoeffs());
                    m_fields[m_intVariables[i]]->SetPhysState(false);
                }

                // Perform any solver-specific post-integration steps
                if (v_PostIntegrate(step))
                {
                    break;
                }

                // search for NaN and quit if found
                bool nanFound = false;
                for (i = 0; i < nvariables; ++i)
                {
                    if (Vmath::Nnan(fields[i].num_elements(), fields[i], 1) > 0)
                    {
                        cout << "NaN found in variable \""
                             << m_session->GetVariable(i)
                             << "\", terminating" << endl;
                        nanFound = true;
                    }
                }

                if (nanFound)
                {
                    break;
                }

                // Update filters
                std::vector<FilterSharedPtr>::iterator x;
                for (x = m_filters.begin(); x != m_filters.end(); ++x)
                {
                    (*x)->Update(m_fields, m_time);
                }

                // Write out checkpoint files
                if ((m_checksteps && step && !((step + 1) % m_checksteps)) ||
                     doCheckTime)
                {
                    if(m_HomogeneousType == eHomogeneous1D)
                    {
                        vector<bool> transformed(nfields, false);
                        for(i = 0; i < nfields; i++)
                        {
                            if (m_fields[i]->GetWaveSpace())
                            {
                                m_fields[i]->SetWaveSpace(false);
                                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                                      m_fields[i]->UpdatePhys());
                                m_fields[i]->SetPhysState(true);
                                transformed[i] = true;
                            }
                        }
                        Checkpoint_Output(m_nchk++);
                        for(i = 0; i < nfields; i++)
                        {
                            if (transformed[i])
                            {
                                m_fields[i]->SetWaveSpace(true);
                                m_fields[i]->HomogeneousFwdTrans(
                                    m_fields[i]->GetPhys(),
                                    m_fields[i]->UpdatePhys());
                                m_fields[i]->SetPhysState(false);
                            }
                        }
                    }
                    else
                    {
                        Checkpoint_Output(m_nchk++);
                    }
                    doCheckTime = false;
                }

                // Step advance
                ++step;
            }
            
            // Print out summary statistics
            if (m_session->GetComm()->GetRank() == 0)
            {
                if (m_cflSafetyFactor > 0.0)
                {
                    cout << "CFL safety factor : " << m_cflSafetyFactor << endl
                         << "CFL time-step     : " << m_timestep        << endl;
                }

                if (m_session->GetSolverInfo("Driver") != "SteadyState")
                {
                    cout << "Time-integration  : " << intTime  << "s"   << endl;
                }
            }
            
            // If homogeneous, transform back into physical space if necessary.
            if(m_HomogeneousType == eHomogeneous1D)
            {
                for(i = 0; i < nfields; i++)
                {
                    if (m_fields[i]->GetWaveSpace())
                    {
                        m_fields[i]->SetWaveSpace(false);
                        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                              m_fields[i]->UpdatePhys());
                        m_fields[i]->SetPhysState(true);
                    }
                }
            }
            else
            {
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[m_intVariables[i]]->SetPhys(fields[i]);
                    m_fields[m_intVariables[i]]->SetPhysState(true);
                }
            }

            // Finalise filters
            for (x = m_filters.begin(); x != m_filters.end(); ++x)
            {
                (*x)->Finalise(m_fields, m_time);
            }
            
            // Print for 1D problems
            if(m_spacedim == 1)
            {
                v_AppendOutput1D(fields);   
            }
        }
        
        /**
         * @brief Sets the initial conditions.
         */
        void UnsteadySystem::v_DoInitialise()
        {
            CheckForRestartTime(m_time);
            SetBoundaryConditions(m_time);
            SetInitialConditions(m_time);
        }
        
        /**
         * @brief Prints a summary with some information regards the 
         * time-stepping.
         */
        void UnsteadySystem::v_GenerateSummary(SummaryList& s)
        {
            EquationSystem::v_GenerateSummary(s);
            AddSummaryItem(s, "Advection",
                           m_explicitAdvection ? "explicit" : "implicit");

            if(m_session->DefinesSolverInfo("AdvectionType"))
            {
                AddSummaryItem(s, "AdvectionType",
                               m_session->GetSolverInfo("AdvectionType"));
            }

            AddSummaryItem(s, "Diffusion",
                           m_explicitDiffusion ? "explicit" : "implicit");

            if (m_session->GetSolverInfo("EQTYPE")
                    == "SteadyAdvectionDiffusionReaction")
            {
                AddSummaryItem(s, "Reaction",
                               m_explicitReaction  ? "explicit" : "implicit");
            }

            AddSummaryItem(s, "Time Step", m_timestep);
            AddSummaryItem(s, "No. of Steps", m_steps);
            AddSummaryItem(s, "Checkpoints (steps)", m_checksteps);
            AddSummaryItem(s, "Integration Type",
                           LibUtilities::TimeIntegrationMethodMap[
                               m_intScheme->GetIntegrationMethod()]);
        }
        
        /**
         * Stores the solution in a file for 1D problems only. This method has 
         * been implemented to facilitate the post-processing for 1D problems.
         */
        void UnsteadySystem::v_AppendOutput1D(
            Array<OneD, Array<OneD, NekDouble> > &solution1D)
        {
            // Coordinates of the quadrature points in the real physical space
            Array<OneD,NekDouble> x(GetNpoints());
            Array<OneD,NekDouble> y(GetNpoints());
            Array<OneD,NekDouble> z(GetNpoints());
            m_fields[0]->GetCoords(x, y, z);
            
            // Print out the solution in a txt file
            ofstream outfile;
            outfile.open("solution1D.txt");
            for(int i = 0; i < GetNpoints(); i++)
            {
                outfile << scientific << setw (17) << setprecision(16) << x[i]
                        << "  " << solution1D[0][i] << endl;
            }
            outfile << endl << endl;
            outfile.close();
        }

        void UnsteadySystem::v_NumericalFlux(
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            ASSERTL0(false, 
                     "This function is not implemented for this equation.");
        }

        void UnsteadySystem::v_NumericalFlux(
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &numfluxX,
            Array<OneD, Array<OneD, NekDouble> > &numfluxY )
        {
            ASSERTL0(false, 
                     "This function is not implemented for this equation.");
        }

        void UnsteadySystem::v_NumFluxforScalar(
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            int i, j;
            int nTraceNumPoints = GetTraceNpoints();
            int nvariables      = m_fields.num_elements();
            int nqvar           = uflux.num_elements();

            Array<OneD, NekDouble > Fwd     (nTraceNumPoints);
            Array<OneD, NekDouble > Bwd     (nTraceNumPoints);
            Array<OneD, NekDouble > Vn      (nTraceNumPoints, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTraceNumPoints, 0.0);

            // Get the sign of (v \cdot n), v = an arbitrary vector

            // Evaulate upwind flux:
            // uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
            for (j = 0; j < nqvar; ++j)
            {
                for (i = 0; i < nvariables ; ++i)
                {
                    // Compute Fwd and Bwd value of ufield of i direction
                    m_fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);

                    // if Vn >= 0, flux = uFwd, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd

                    // else if Vn < 0, flux = uBwd, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

                    m_fields[i]->GetTrace()->Upwind(m_traceNormals[j], 
                                                    Fwd, Bwd, fluxtemp);

                    // Imposing weak boundary condition with flux
                    // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd

                    // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd

                    if(m_fields[0]->GetBndCondExpansions().num_elements())
                    {
                        WeakPenaltyforScalar(i, ufield[i], fluxtemp);
                    }

                    // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n}), 
                    // i.e,
                    // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n})
                    // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n})

                    // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n}), 
                    // i.e,
                    // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n})
                    // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n})

                    Vmath::Vmul(nTraceNumPoints, 
                                m_traceNormals[j], 1, 
                                fluxtemp, 1, 
                                uflux[j][i], 1);
                }
            }
        }

        
        
        void UnsteadySystem::v_NumFluxforVector(
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            int nTraceNumPoints = GetTraceNpoints();
            int nvariables = m_fields.num_elements();
            int nqvar = qfield.num_elements();

            NekDouble C11 = 1.0;
            Array<OneD, NekDouble > Fwd(nTraceNumPoints);
            Array<OneD, NekDouble > Bwd(nTraceNumPoints);
            Array<OneD, NekDouble > Vn (nTraceNumPoints, 0.0);

            Array<OneD, NekDouble > qFwd     (nTraceNumPoints);
            Array<OneD, NekDouble > qBwd     (nTraceNumPoints);
            Array<OneD, NekDouble > qfluxtemp(nTraceNumPoints, 0.0);

            Array<OneD, NekDouble > uterm(nTraceNumPoints);

            // Evaulate upwind flux:
            // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
            for (int i = 0; i < nvariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTraceNumPoints, 0.0);
                for (int j = 0; j < nqvar; ++j)
                {
                    //  Compute Fwd and Bwd value of ufield of jth direction
                    m_fields[i]->GetFwdBwdTracePhys(qfield[j][i],qFwd,qBwd);

                    // if Vn >= 0, flux = uFwd, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick 
                    // qflux = qBwd = q+
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick 
                    // qflux = qBwd = q-

                    // else if Vn < 0, flux = uBwd, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick 
                    // qflux = qFwd = q-
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick 
                    // qflux = qFwd = q+

                    m_fields[i]->GetTrace()->Upwind(m_traceNormals[j], 
                                                    qBwd, qFwd, 
                                                    qfluxtemp);
                    
                    Vmath::Vmul(nTraceNumPoints, 
                                m_traceNormals[j], 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);

                    // Generate Stability term = - C11 ( u- - u+ )
                    m_fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);
                    
                    Vmath::Vsub(nTraceNumPoints, 
                                Fwd, 1, Bwd, 1, 
                                uterm, 1);
                    
                    Vmath::Smul(nTraceNumPoints, 
                                -1.0 * C11, uterm, 1, 
                                uterm, 1);

                    // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
                    Vmath::Vadd(nTraceNumPoints, 
                                uterm, 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);

                    // Imposing weak boundary condition with flux
                    if (m_fields[0]->GetBndCondExpansions().num_elements())
                    {
                        WeakPenaltyforVector(i, j, 
                                             qfield[j][i], 
                                             qfluxtemp, 
                                             C11);
                    }

                    // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                    // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
                    // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
                    Vmath::Vadd(nTraceNumPoints, 
                                qfluxtemp, 1, 
                                qflux[i], 1, 
                                qflux[i], 1);
                }
            }
        }

        void UnsteadySystem::CheckForRestartTime(NekDouble &time)
        {
            if (m_session->DefinesFunction("InitialConditions"))
            {
                for (int i = 0; i < m_fields.num_elements(); ++i)
                {
                    LibUtilities::FunctionType vType;

                    vType = m_session->GetFunctionType(
                        "InitialConditions", m_session->GetVariable(i));

                    if (vType == LibUtilities::eFunctionTypeFile)
                    {
                        std::string filename
                            = m_session->GetFunctionFilename(
                                "InitialConditions", m_session->GetVariable(i));

                        fs::path pfilename(filename);

                        // redefine path for parallel file which is in directory
                        if(fs::is_directory(pfilename))
                        {
                            fs::path metafile("Info.xml");
                            fs::path fullpath = pfilename / metafile;
                            filename = LibUtilities::PortablePath(fullpath);
                        }
                        m_fld->ImportFieldMetaData(filename, m_fieldMetaDataMap);

                        // check to see if time defined
                        if (m_fieldMetaDataMap !=
                                LibUtilities::NullFieldMetaDataMap)
                        {
                            LibUtilities::FieldMetaDataMap::iterator iter; 
                            
                            iter = m_fieldMetaDataMap.find("Time");
                            if (iter != m_fieldMetaDataMap.end())
                            {
                                time = boost::lexical_cast<NekDouble>(
                                    iter->second);
                            }
                        }
                        
                        break;
                    }
                }
            }
        }
        
        void UnsteadySystem::WeakPenaltyforScalar(
            const int var,
            const Array<OneD, const NekDouble> &physfield,
                  Array<OneD,       NekDouble> &penaltyflux,
            NekDouble time)
        {
            int i, e, npoints, id1, id2;
            
            // Number of boundary regions
            int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
            int Nfps, numBDEdge;
            int nTraceNumPoints = GetTraceNpoints();
            int cnt = 0;

            Array<OneD, NekDouble > uplus(nTraceNumPoints);

            m_fields[var]->ExtractTracePhys(physfield, uplus);
            for (i = 0; i < nbnd; ++i)
            {
                // Number of boundary expansion related to that region
                numBDEdge = m_fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                
                // Evaluate boundary values g_D or g_N from input files
                LibUtilities::EquationSharedPtr ifunc = 
                    m_session->GetFunction("InitialConditions", 0);
                
                npoints = m_fields[var]->
                    GetBndCondExpansions()[i]->GetNpoints();
                
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
                    Nfps = m_fields[var]->
                        GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    id1 = m_fields[var]->
                        GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = m_fields[0]->GetTrace()->
                        GetPhys_Offset(m_fields[0]->GetTraceMap()->
                                        GetBndCondTraceToGlobalTraceMap(cnt++));

                    // For Dirichlet boundary condition: uflux = g_D
                    if (m_fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vcopy(Nfps, 
                                     &BDphysics[id1], 1, 
                                     &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: uflux = u+
                    else if ((m_fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(Nfps, 
                                     &uplus[id2], 1, 
                                     &penaltyflux[id2], 1);
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
            int i, e, npoints, id1, id2;
            int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
            int numBDEdge, Nfps;
            int nTraceNumPoints = GetTraceNpoints();
            Array<OneD, NekDouble > uterm(nTraceNumPoints);
            Array<OneD, NekDouble > qtemp(nTraceNumPoints);
            int cnt = 0;

            m_fields[var]->ExtractTracePhys(physfield,qtemp);

            for (i = 0; i < nbnd; ++i)
            {
                numBDEdge = m_fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                
                // Evaluate boundary values g_D or g_N from input files
                LibUtilities::EquationSharedPtr ifunc = 
                    m_session->GetFunction("InitialConditions", 0);
                
                npoints = m_fields[var]->
                    GetBndCondExpansions()[i]->GetNpoints();

                Array<OneD,NekDouble> BDphysics(npoints);
                Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);

                m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
                ifunc->Evaluate(x0,x1,x2,time,BDphysics);

                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < numBDEdge ; ++e)
                {
                    Nfps = m_fields[var]->
                        GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);

                    id1 = m_fields[var]->
                        GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = m_fields[0]->GetTrace()->
                        GetPhys_Offset(m_fields[0]->GetTraceMap()->
                                       GetBndCondTraceToGlobalTraceMap(cnt++));

                    // For Dirichlet boundary condition: 
                    //qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
                    if(m_fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vmul(Nfps, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &qtemp[id2], 1, 
                                    &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: qflux = g_N
                    else if((m_fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vmul(Nfps,
                                    &m_traceNormals[dir][id2], 1, 
                                    &BDphysics[id1], 1, 
                                    &penaltyflux[id2], 1);
                    }
                }
            }
        }
	
        /**
         * @brief Return the timestep to be used for the next step in the
         * time-marching loop.
         *
         * This function can be overloaded to facilitate solver which utilise a
         * CFL (or other) parameter to determine a maximum timestep under which
         * the problem remains stable.
         */
        NekDouble UnsteadySystem::GetTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray)
        {
            return v_GetTimeStep(inarray);
        }
        
        /**
         * @brief Return the timestep to be used for the next step in the
         * time-marching loop.
         *
         * @see UnsteadySystem::GetTimeStep
         */
        NekDouble UnsteadySystem::v_GetTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray)
        {
            ASSERTL0(false, "Not defined for this class");
            return 0.0;
        }

        bool UnsteadySystem::v_PreIntegrate(int step)
        {
            return false;
        }

        bool UnsteadySystem::v_PostIntegrate(int step)
        {
            return false;
        }

        bool UnsteadySystem::v_SteadyStateCheck(int step)
        {
            return false;
        }

        void UnsteadySystem::SVVVarDiffCoeff(
            const Array<OneD, Array<OneD, NekDouble> >  vel,
                  StdRegions::VarCoeffMap              &varCoeffMap)
        {
            int phystot = m_fields[0]->GetTotPoints();
            int nvel = vel.num_elements();

            Array<OneD, NekDouble> varcoeff(phystot),tmp;

            // calculate magnitude of v
            Vmath::Vmul(phystot,vel[0],1,vel[0],1,varcoeff,1);
            for(int n = 1; n < nvel; ++n)
            {
                Vmath::Vvtvp(phystot,vel[n],1,vel[n],1,varcoeff,1,varcoeff,1);
            }
            Vmath::Vsqrt(phystot,varcoeff,1,varcoeff,1);

            for(int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
            {
                int offset = m_fields[0]->GetPhys_Offset(i);
                int nq = m_fields[0]->GetExp(i)->GetTotPoints();
                Array<OneD, NekDouble> unit(nq,1.0);

                int nmodes = 0;

                for(int n = 0; n < m_fields[0]->GetExp(i)->GetNumBases(); ++n)
                {
                    nmodes = max(nmodes,
                                 m_fields[0]->GetExp(i)->GetBasisNumModes(n));
                }

                NekDouble h = m_fields[0]->GetExp(i)->Integral(unit);
                h = pow(h,(NekDouble) (1.0/nvel))/((NekDouble) nmodes);

                Vmath::Smul(nq,h,varcoeff+offset,1,tmp = varcoeff+offset,1);
            }

            // set up map with eVarCoffLaplacian key
            varCoeffMap[StdRegions::eVarCoeffLaplacian] = varcoeff;
        }
    }
}
