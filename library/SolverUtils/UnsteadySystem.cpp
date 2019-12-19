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

#include <boost/format.hpp>

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

            m_session->LoadParameter("CheckNanSteps", m_nanSteps, 1);
            // Steady state tolerance
            m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 0.0);
            // Frequency for checking steady state
            m_session->LoadParameter("SteadyStateSteps",
                                          m_steadyStateSteps, 1);

            // For steady problems, we do not initialise the time integration
            if (m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"))
            {
                m_intScheme = LibUtilities::GetTimeIntegrationWrapperFactory().
                    CreateInstance(m_session->GetSolverInfo(
                                       "TIMEINTEGRATIONMETHOD"));

                // Load generic input parameters
                m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
                m_session->LoadParameter("CFL", m_cflSafetyFactor, 0.0);
                m_session->LoadParameter("CFL_End", m_CFLEnd, 0.0);
                m_session->LoadParameter("CFLGrowth", m_CFLGrowth, 1.0);

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
                if(m_CFLGrowth > 1.0)
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
                m_filters.push_back(
                    GetFilterFactory().CreateInstance(
                        x.first, m_session, x.second));
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
                case LibUtilities::eBackwardEuler:
                case LibUtilities::eDIRKOrder2:
                case LibUtilities::eDIRKOrder3:
                case LibUtilities::eDIRKOrder3Stage5:
                case LibUtilities::eDIRKOrder4Stage6:
                {
                    TimeStability = 2.0;
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

            int nwidthcolm = 10;

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
            m_intSoln = m_intScheme->InitializeScheme(
                m_timestep, fields, m_time, m_ode);

            // Initialise filters
            for (auto &x : m_filters)
            {
                x->setFilterOperators(m_FilterOperators);
                x->Initialise(m_fields, m_time);
            }

            LibUtilities::Timer     timer;
            bool      doCheckTime   = false;
            int       step          = m_initialStep;
            int       stepCounter   = 0;
            int       restartStep   = -1;
            NekDouble intTime       = 0.0;
            NekDouble lastCheckTime = 0.0;
            NekDouble cpuTime       = 0.0;
            NekDouble elapsed       = 0.0;

            m_TotNewtonIts  = 0;
            m_TotGMRESIts   = 0;
            m_TotOdeRHS     = 0;
            m_TotImpStages  = 0;

            NekDouble tmp_cflSafetyFactor = m_cflSafetyFactor;
            m_CalcuPrecMatCounter = m_PrcdMatFreezNumb;
            bool flagFreezeCFL = false;

            m_timestepMax = m_timestep;
            while (step   < m_steps ||
                   m_time < m_fintime - NekConstants::kNekZeroTol)
            {
                restartStep++;

                // cout    <<" m_TotLinItePerStep= "<<m_TotLinItePerStep
                //         <<" m_StagesPerStep= "<<m_StagesPerStep
                //         <<" m_maxLinItePerNewton= "<<m_maxLinItePerNewton<<endl;
                if(m_CFLEnd>tmp_cflSafetyFactor)
                {
                    if( m_steadyStateTol > 0.0 &&
                        (NekDouble(m_TotLinItePerStep)/NekDouble(m_StagesPerStep)>0.5*NekDouble(m_maxLinItePerNewton)))
                    {
                        // cout <<"WARNINGL1(false,tmp_cflSafetyFactor *= 0.9; );"<<endl;
                        
                        tmp_cflSafetyFactor = 0.9*tmp_cflSafetyFactor;
                        flagFreezeCFL = true;
                        WARNINGL1(false," tmp_cflSafetyFactor *= 0.9; ");
                    }
                    else if(flagFreezeCFL)
                    {
                        flagFreezeCFL = false;
                    }
                    else
                    {
                        if (m_steadyStateTol > 0.0 && (!((step+1)%m_steadyStateSteps)) )
                        {
                            if(restartStep>1)
                            {
                                tmp_cflSafetyFactor = min(m_CFLEnd,std::pow(m_steadyStateRes0/m_steadyStateRes,0.7)*tmp_cflSafetyFactor);
                            }
                        }
                        else
                        {
                            if(m_CFLGrowth > 1.0&&m_cflSafetyFactor<m_CFLEnd)
                            {
                                tmp_cflSafetyFactor = min(m_CFLEnd,m_CFLGrowth*tmp_cflSafetyFactor);
                            }
                        }
                    }
                }
                // Flag to update AV
                m_calcuPhysicalAV = true;

                // Frozen preconditioner checks
                if(    (m_CalcuPrecMatCounter>=m_PrcdMatFreezNumb)
                    ||(m_time + m_timestep > m_fintime && m_fintime > 0.0)
                    ||(m_checktime && m_time + m_timestep - lastCheckTime >= m_checktime))
                {

                    m_CalcuPrecMatFlag      =   true;
                    m_CalcuPrecMatCounter   =   0;
                    m_cflSafetyFactor       =   tmp_cflSafetyFactor;

                    if (m_cflSafetyFactor)
                    {
                        m_timestep = GetTimeStep(fields);
                    }

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

                m_CalcuPrecMatCounter++;
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
                if (m_TimeIncrementFactor>1.0)
                {
                    NekDouble timeincrementFactor = m_TimeIncrementFactor;
                    m_timestep  *=  timeincrementFactor;

                    if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
                    {
                        m_timestep = m_fintime - m_time;
                    }
                }
                // if(m_CalcuPrecMatCounter>=m_PrcdMatFreezNumb)
                // {

                //     m_CalcuPrecMatFlag      =   true;
                //     m_CalcuPrecMatCounter   =   0;
                // }
                // m_CalcuPrecMatCounter++;

#endif

                // Perform any solver-specific pre-integration steps
                timer.Start();
                if (v_PreIntegrate(step))
                {
                    break;
                }

                m_StagesPerStep = 0;
                m_TotLinItePerStep = 0;

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
                    cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
                         << "Steps: " << setw(8)  << left << step+1 << " "
                         << "Time: "  << setw(8) << left << m_time;

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
                    if (m_cflSafetyFactor||m_TimeIncrementFactor>1.0)
#else
                    if (m_cflSafetyFactor)
#endif
                    {
                        cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
                             << " CFL: " << m_cflSafetyFactor
                             << " NonAcous-CFL : " << m_cflNonAcoustic
                             << " Time-step: " << m_timestep;
                    }

                    stringstream ss;
                    ss << cpuTime << "s";
                    cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
                         << " CPU Time: " << left
                         << ss.str();

                    cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
                         <<" INT Time: "<< intTime<<"s"<<endl;
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
                    if(m_flagImplItsStatistcs)
                    {
                        cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
                             << "       &&" 
                             << " TotImpStages= " << m_TotImpStages 
                             << " TotNewtonIts= " << m_TotNewtonIts
                             << " TotGMRESIts = " << m_TotGMRESIts  
                             << " TotOdeRHS   = " << m_TotOdeRHS    
                             <<endl;
                    }
#endif

                    cpuTime = 0.0;
                }

                // Transform data into coefficient space
                for (i = 0; i < nvariables; ++i)
                {
                    m_fields[m_intVariables[i]]->SetPhys(fields[i]);
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
                    if (CheckSteadyState(step,intTime))
                    {
                        if (m_comm->GetRank() == 0)
                        {
                            cout << "Reached Steady State to tolerance "
                                 << m_steadyStateTol << endl;
                        }
                        break;
                    }
                }

                // search for NaN and quit if found
                if (m_nanSteps && !((step+1) % m_nanSteps) )
                {
                    int nanFound = 0;
                    for (i = 0; i < nvariables; ++i)
                    {
                        if (Vmath::Nnan(fields[i].num_elements(),
                                fields[i], 1) > 0)
                        {
                            nanFound = 1;
                        }
                    }
                    m_session->GetComm()->AllReduce(nanFound,
                                LibUtilities::ReduceMax);
                    ASSERTL0 (!nanFound,
                                "NaN found during time integration.");
                }
                // Update filters
                for (auto &x : m_filters)
                {
                    x->Update(m_fields, m_time);
                }

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
                    cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
                         << "CFL safety factor : " << m_cflSafetyFactor << endl
                         << "CFL time-step     : " << m_timestep        << endl;
                }

                if (m_session->GetSolverInfo("Driver") != "SteadyState")
                {
                    cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
                         << "Time-integration  : " << intTime  << "s"   << endl;
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
                x->Finalise(m_fields, m_time);
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

        void UnsteadySystem::CheckForRestartTime(NekDouble &time, int &nchk)
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

        void UnsteadySystem::InitializeSteadyState()
        {
            if (m_steadyStateTol > 0.0)
            {
                const int nPoints = m_fields[0]->GetTotPoints();
                m_previousSolution = Array<OneD, Array<OneD, NekDouble> > (
                            m_fields.num_elements());

                for (int i = 0; i < m_fields.num_elements(); ++i)
                {
                    m_previousSolution[i] = Array<OneD, NekDouble>(nPoints);
                    Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1,
                                          m_previousSolution[i], 1);
                }

                if (m_comm->GetRank() == 0)
                {
                    std::string fName = m_session->GetSessionName() +
                        std::string(".resd");
                    m_errFile.open(fName.c_str());
                    m_errFile << setw(26) << left << "# Time";

                    m_errFile << setw(26) << left << "CPU_Time";

                    m_errFile << setw(26) << left << "Step";

                    for (int i = 0; i < m_fields.num_elements(); ++i)
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
            return CheckSteadyState(step,0.0);
        }

        bool UnsteadySystem::CheckSteadyState(int step, NekDouble totCPUTime)
        {
            const int nPoints = GetTotPoints();
            const int nFields = m_fields.num_elements();

            // Holds L2 errors.
            Array<OneD, NekDouble> L2       (nFields);

            SteadyStateResidual(step, L2);

            if (m_comm->GetRank() == 0 && ( ((step+1) % m_infosteps == 0)||((step== m_initialStep)) ))
            {
                // Output time
                m_errFile << boost::format("%25.19e") % m_time;

                m_errFile << " "<< boost::format("%25.19e") % totCPUTime;

                int stepp = step +1;

                m_errFile << " "<< boost::format("%25.19e") % stepp;

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

            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1,
                                      m_previousSolution[i], 1);
            }

            m_steadyStateRes0 = m_steadyStateRes;
            m_steadyStateRes = maxL2;

            return false;
        }

        void UnsteadySystem::v_SteadyStateResidual(
            int                         step, 
            Array<OneD, NekDouble>      &L2)
        {
            const int nPoints = GetTotPoints();
            const int nFields = m_fields.num_elements();

            // Holds L2 errors.
            Array<OneD, NekDouble> RHSL2    (nFields);
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
        }
    }
}
