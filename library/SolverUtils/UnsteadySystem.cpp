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
using namespace std;

#include <boost/core/ignore_unused.hpp>
#include <boost/format.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/UnsteadySystem.h>

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
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph)
            : EquationSystem(pSession, pGraph),
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

            m_session->LoadParameter("CheckAbortSteps", m_abortSteps, 1);
            // Steady state tolerance
            m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 0.0);
            // Frequency for checking steady state
            m_session->LoadParameter("SteadyStateSteps",
                                          m_steadyStateSteps, 1);

            // For steady problems, we do not initialise the time integration
            if (m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"))
            {
                std::string methodName = m_session->GetSolverInfo(
                                                "TIMEINTEGRATIONMETHOD" );
                m_intScheme = LibUtilities::GetTimeIntegrationSchemeFactory()
		    .CreateInstance( methodName, "", 0,
				     std::vector<NekDouble>());

                // Load generic input parameters
                m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
                m_session->LoadParameter("IO_FiltersInfoSteps",
                    m_filtersInfosteps, 10.0 * m_infosteps);
                m_session->LoadParameter("CFL", m_cflSafetyFactor, 0.0);
                m_session->LoadParameter("CFLEnd", m_CFLEnd, 0.0);
                m_session->LoadParameter("CFLGrowth", m_CFLGrowth, 1.0);

                // Time tolerance between filter update time and time integration
                m_session->LoadParameter("FilterTimeWarning",
                                         m_filterTimeWarning, 1);

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
                else
                {
                    ASSERTL0(m_timestep != 0.0,
                             "Need to set either TimeStep or CFL");
                }

                // Ensure that there is no conflict of parameters
                if (m_CFLGrowth > 1.0)
                {
                    // Check final condition
                    ASSERTL0(m_CFLEnd >= m_cflSafetyFactor,
                             "m_CFLEnd >= m_cflSafetyFactor required");
                }

                // Set up time to be dumped in field information
                m_fieldMetaDataMap["Time"] =
                        boost::lexical_cast<std::string>(m_time);
            }

            // By default attempt to forward transform initial condition.
            m_homoInitialFwd = true;

            // Set up filters
            for (auto &x : m_session->GetFilters())
            {
                m_filters.push_back(make_pair(x.first, GetFilterFactory().CreateInstance(
                        x.first, m_session, shared_from_this(), x.second)));
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
            return m_intScheme->GetTimeStability();
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
            int nfields = m_fields.size();

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
            if(m_HomogeneousType != eNotHomogeneous && m_homoInitialFwd)
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

            // Order storage to list time-integrated fields first.
            for(i = 0; i < nvariables; ++i)
            {
                fields[i] = m_fields[m_intVariables[i]]->GetPhys();
                m_fields[m_intVariables[i]]->SetPhysState(false);
            }

            // Initialise time integration scheme
            m_intSoln = m_intScheme->InitializeScheme( m_timestep, fields,
                                                       m_time, m_ode );

            // Initialise filters
            for( auto &x : m_filters )
            {
                x.second->Initialise(m_fields, m_time);
            }

            LibUtilities::Timer         timer;
            bool      doCheckTime       = false;
            int       step              = m_initialStep;
            int       stepCounter       = 0;
            NekDouble intTime           = 0.0;
            NekDouble lastCheckTime     = 0.0;
            NekDouble cpuTime           = 0.0;
            NekDouble cpuPrevious       = 0.0;
            NekDouble elapsed           = 0.0;
            NekDouble totFilterTime     = 0.0;

            Array<OneD, int> abortFlags(2, 0);
            string    abortFile     = "abort";
            if (m_session->DefinesSolverInfo("CheckAbortFile"))
            {
                abortFile = m_session->GetSolverInfo("CheckAbortFile");
            }

            m_timestepMax = m_timestep;
            while ((step   < m_steps ||
                   m_time < m_fintime - NekConstants::kNekZeroTol) &&
                   abortFlags[1] == 0)
            {
                if (m_CFLGrowth > 1.0 && m_cflSafetyFactor < m_CFLEnd)
                {
                    m_cflSafetyFactor = min(
                        m_CFLEnd, m_CFLGrowth * m_cflSafetyFactor);
                }
                
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
                    stepCounter, m_timestep, m_intSoln, m_ode);
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
                    cpuPrevious = cpuTime;
                    cpuTime = 0.0;
                }

                // Transform data into coefficient space
                for (i = 0; i < nvariables; ++i)
                {
                    // copy fields into ExpList::m_phys and assign the new
                    // array to fields
                    m_fields[m_intVariables[i]]->SetPhys(fields[i]);
                    fields[i] = m_fields[m_intVariables[i]]->UpdatePhys();
                    if( v_RequireFwdTrans() )
                    {
                        m_fields[m_intVariables[i]]->FwdTrans_IterPerExp(
                            fields[i],
                            m_fields[m_intVariables[i]]->UpdateCoeffs());
                    }
                    m_fields[m_intVariables[i]]->SetPhysState(false);
                }
                
                // Perform any solver-specific post-integration steps
                if (v_PostIntegrate(step))
                {
                    break;
                }

                // Check for steady-state
                if (m_steadyStateTol > 0.0 && (!((step+1)%m_steadyStateSteps)) )
                {
                    if (CheckSteadyState(step))
                    {
                        if (m_comm->GetRank() == 0)
                        {
                            cout << "Reached Steady State to tolerance "
                                 << m_steadyStateTol << endl;
                        }
                        break;
                    }
                }

                // test for abort conditions (nan, or abort file)
                if (m_abortSteps && !((step+1) % m_abortSteps) )
                {
                    abortFlags[0] = 0;
                    for (i = 0; i < nvariables; ++i)
                    {
                        if (Vmath::Nnan(fields[i].size(),
                                fields[i], 1) > 0)
                        {
                            abortFlags[0] = 1;
                        }
                    }

                    //rank zero looks for abort file and deltes it
                    //if it exists. The communicates the abort
                    if(m_session->GetComm()->GetRank() == 0)
                    {
                        if(boost::filesystem::exists(abortFile))
                        {
                            boost::filesystem::remove(abortFile);
                            abortFlags[1] = 1;
                        }
                    }

                    m_session->GetComm()->AllReduce(abortFlags,
                                LibUtilities::ReduceMax);

                    ASSERTL0 (!abortFlags[0],
                                "NaN found during time integration.");
                }

                // Update filters
                for (auto &x : m_filters)
                {
                    timer.Start();
                    x.second->Update(m_fields, m_time);
                    timer.Stop();
                    elapsed = timer.TimePerTest(1);
                    totFilterTime += elapsed;

                    // Write out individual filter status information
                    if(m_session->GetComm()->GetRank() == 0 &&
                    !((step+1) % m_filtersInfosteps) && !m_filters.empty() &&
                    m_session->DefinesCmdLineArgument("verbose"))
                    {
                        stringstream s0;
                        s0 << x.first << ":";
                        stringstream s1;
                        s1 << elapsed << "s";
                        stringstream s2;
                        s2 << elapsed / cpuPrevious * 100 << "%";
                        cout << "CPU time for filter " << setw(25) << left
                            << s0.str() << setw(12) << left << s1.str() <<
                            endl << "\t Percentage of time integration:     "
                             << setw(10) << left << s2.str() << endl;
                    }
                }

                // Write out overall filter status information
                if (m_session->GetComm()->GetRank() == 0 &&
                    !((step+1) % m_filtersInfosteps) && !m_filters.empty())
                 {
                    stringstream ss;
                    ss << totFilterTime << "s";
                    cout << "Total filters CPU Time:\t\t\t     " << setw(10)
                        << left << ss.str() << endl;
                 }
                totFilterTime = 0.0;

                // Write out checkpoint files
                if ((m_checksteps && !((step + 1) % m_checksteps)) ||
                     doCheckTime)
                {
                    if(m_HomogeneousType != eNotHomogeneous)
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
                        Checkpoint_Output(m_nchk);
                        m_nchk++;
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
                        Checkpoint_Output(m_nchk);
                        m_nchk++;
                    }
                    doCheckTime = false;
                }

                // Step advance
                ++step;
                ++stepCounter;
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
            if(m_HomogeneousType != eNotHomogeneous)
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
            for (auto &x : m_filters)
            {
                x.second->Finalise(m_fields, m_time);
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
            CheckForRestartTime(m_time, m_nchk);
            SetBoundaryConditions(m_time);
            SetInitialConditions(m_time);
            InitializeSteadyState();
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

            AddSummaryItem( s, "Time Step", m_timestep );
            AddSummaryItem( s, "No. of Steps", m_steps );
            AddSummaryItem( s, "Checkpoints (steps)", m_checksteps );
            AddSummaryItem( s, "Integration Type", m_intScheme->GetName() );
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

        void UnsteadySystem::CheckForRestartTime(NekDouble &time, int &nchk)
        {
            if (m_session->DefinesFunction("InitialConditions"))
            {
                for (int i = 0; i < m_fields.size(); ++i)
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
                        LibUtilities::FieldIOSharedPtr fld =
                            LibUtilities::FieldIO::CreateForFile(
                                m_session, filename);
                        fld->ImportFieldMetaData(filename, m_fieldMetaDataMap);

                        // check to see if time defined
                        if (m_fieldMetaDataMap !=
                                LibUtilities::NullFieldMetaDataMap)
                        {
                            auto iter = m_fieldMetaDataMap.find("Time");
                            if (iter != m_fieldMetaDataMap.end())
                            {
                                time = boost::lexical_cast<NekDouble>(
                                    iter->second);
                            }

                            iter = m_fieldMetaDataMap.find("ChkFileNum");
                            if (iter != m_fieldMetaDataMap.end())
                            {
                                nchk = boost::lexical_cast<NekDouble>(
                                    iter->second);
                            }
                        }

                        break;
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
            boost::ignore_unused(inarray);
            NEKERROR(ErrorUtil::efatal, "Not defined for this class");
            return 0.0;
        }

        bool UnsteadySystem::v_PreIntegrate(int step)
        {
            boost::ignore_unused(step);
            return false;
        }

        bool UnsteadySystem::v_PostIntegrate(int step)
        {
            boost::ignore_unused(step);
            return false;
        }

        void UnsteadySystem::SVVVarDiffCoeff(
            const Array<OneD, Array<OneD, NekDouble> >  vel,
                  StdRegions::VarCoeffMap              &varCoeffMap)
        {
            int phystot = m_fields[0]->GetTotPoints();
            int nvel = vel.size();

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

        void UnsteadySystem::InitializeSteadyState()
        {
            if (m_steadyStateTol > 0.0)
            {
                const int nPoints = m_fields[0]->GetTotPoints();
                m_previousSolution = Array<OneD, Array<OneD, NekDouble> > (
                            m_fields.size());

                for (int i = 0; i < m_fields.size(); ++i)
                {
                    m_previousSolution[i] = Array<OneD, NekDouble>(nPoints);
                    Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1,
                                          m_previousSolution[i], 1);
                }

                if (m_comm->GetRank() == 0)
                {
                    std::string fName = m_session->GetSessionName() +
                        std::string(".res");
                    m_errFile.open(fName.c_str());
                    m_errFile << setw(26) << left << "# Time";

                    for (int i = 0; i < m_fields.size(); ++i)
                    {
                        m_errFile << setw(26) << m_session->GetVariables()[i];
                    }

                    m_errFile << endl;
                }
            }
        }

        /**
        * @brief Calculate whether the system has reached a steady state by
        * observing residuals to a user-defined tolerance.
        */
        bool UnsteadySystem::CheckSteadyState(int step)
        {
            const int nPoints = GetTotPoints();
            const int nFields = m_fields.size();

            // Holds L2 errors.
            Array<OneD, NekDouble> L2       (nFields);
            Array<OneD, NekDouble> residual (nFields);
            Array<OneD, NekDouble> reference(nFields);

            for (int i = 0; i < nFields; ++i)
            {
                Array<OneD, NekDouble> tmp(nPoints);

                Vmath::Vsub(nPoints, m_fields[i]->GetPhys(), 1,
                                     m_previousSolution[i], 1, tmp, 1);
                Vmath::Vmul(nPoints, tmp, 1, tmp, 1, tmp, 1);
                residual[i] = Vmath::Vsum(nPoints, tmp, 1);

                Vmath::Vmul(nPoints, m_previousSolution[i], 1,
                                     m_previousSolution[i], 1, tmp, 1);
                reference[i] = Vmath::Vsum(nPoints, tmp, 1);
            }

            m_comm->AllReduce(residual , LibUtilities::ReduceSum);
            m_comm->AllReduce(reference, LibUtilities::ReduceSum);

            // L2 error
            for (int i = 0; i < nFields; ++i)
            {
                reference[i] = (reference[i] == 0) ? 1 : reference[i];
                L2[i] = sqrt(residual[i] / reference[i]);
            }

            if (m_comm->GetRank() == 0 && ((step+1) % m_infosteps == 0))
            {
                // Output time
                m_errFile << boost::format("%25.19e") % m_time;

                // Output residuals
                for (int i = 0; i < nFields; ++i)
                {
                    m_errFile << " " << boost::format("%25.19e") % L2[i];
                }

                m_errFile << endl;
            }

            // Calculate maximum L2 error
            NekDouble maxL2 = Vmath::Vmax(nFields, L2, 1);

            if (m_session->DefinesCmdLineArgument("verbose") &&
                m_comm->GetRank() == 0 && ((step+1) % m_infosteps == 0))
            {
                cout << "-- Maximum L^2 residual: " << maxL2 << endl;
            }

            if (maxL2 <= m_steadyStateTol)
            {
                return true;
            }

            for (int i = 0; i < m_fields.size(); ++i)
            {
                Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1,
                                      m_previousSolution[i], 1);
            }

            return false;
        }
    }
}
