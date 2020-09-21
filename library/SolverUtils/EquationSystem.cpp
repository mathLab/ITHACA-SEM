////////////////////////////////////////////////////////////////////////////////
//
// File EquationSystem.cpp
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
// Description: Main wrapper class for Advection Diffusion Reaction Solver
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <FieldUtils/Interpolator.h>
#include <SolverUtils/EquationSystem.h>

#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>

#include <MultiRegions/ExpList.h>       
#include <MultiRegions/ExpList2D.h>     // for ExpList2D, etc
#include <MultiRegions/ExpList3D.h>     // for ExpList3D
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>

#include <GlobalMapping/Mapping.h>

#include <boost/format.hpp>

#include <iostream>
#include <string>

using namespace std;

namespace Nektar
{
    namespace SolverUtils
    {

        std::string EquationSystem::equationSystemTypeLookupIds[2] = {
            LibUtilities::SessionReader::RegisterEnumValue("DEALIASING",
                "True", 0),
            LibUtilities::SessionReader::RegisterEnumValue("DEALIASING",
                "False", 1)};

        /**
         * @class EquationSystem
         *
         * This class is a base class for all solver implementations. It
         * provides the underlying generic functionality and interface for
         * solving equations.
         *
         * To solve a steady-state equation, create a derived class from this
         * class and reimplement the virtual functions to provide custom
         * implementation for the problem.
         *
         * To solve unsteady problems, derive from the UnsteadySystem class
         * instead which provides general time integration.
         */
        EquationSystemFactory& GetEquationSystemFactory()
        {
            static EquationSystemFactory instance;
            return instance;
        }

        /**
         * This constructor is protected as the objects of this class are never
         * instantiated directly.
         * @param   pSession The session reader holding problem parameters.
         */
        EquationSystem::EquationSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph)
            : m_comm (pSession->GetComm()),
              m_session (pSession),
              m_graph (pGraph),
              m_lambda (0),
              m_fieldMetaDataMap(LibUtilities::NullFieldMetaDataMap)
        {
            // set up session names in fieldMetaDataMap
            const vector<std::string> filenames = m_session->GetFilenames();

            for(int i = 0; i < filenames.size(); ++i)
            {
                string sessionname = "SessionName";
                sessionname += boost::lexical_cast<std::string>(i);
                m_fieldMetaDataMap[sessionname] = filenames[i];
                m_fieldMetaDataMap["ChkFileNum"] =
                        boost::lexical_cast<std::string>(0);
            }

        }

        /**
         * @brief Initialisation object for EquationSystem.
         */
        void EquationSystem::v_InitObject()
        {
            // Save the basename of input file name for output details
            m_sessionName = m_session->GetSessionName();

            // Instantiate a field reader/writer
            m_fld = LibUtilities::FieldIO::CreateDefault(m_session);

            // Also read and store the boundary conditions
            m_boundaryConditions =
                MemoryManager<SpatialDomains::BoundaryConditions>::
                    AllocateSharedPtr(m_session, m_graph);

            // Set space dimension for use in class
            m_spacedim = m_graph->GetSpaceDimension();

            // Setting parameteres for homogenous problems
            m_HomoDirec             = 0;
            m_useFFT                = false;
            m_homogen_dealiasing    = false;
            m_singleMode            = false;
            m_halfMode              = false;
            m_multipleModes         = false;
            m_HomogeneousType       = eNotHomogeneous;

            if (m_session->DefinesSolverInfo("HOMOGENEOUS"))
            {
                std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
                m_spacedim          = 3;

                if ((HomoStr == "HOMOGENEOUS1D") || (HomoStr == "Homogeneous1D")
                    || (HomoStr == "1D") || (HomoStr == "Homo1D"))
                {
                    m_HomogeneousType = eHomogeneous1D;
                    m_session->LoadParameter("LZ", m_LhomZ);
                    m_HomoDirec       = 1;

                    if(m_session->DefinesSolverInfo("ModeType"))
                    {
                        m_session->MatchSolverInfo("ModeType", "SingleMode",
                                                   m_singleMode, false);
                        m_session->MatchSolverInfo("ModeType", "HalfMode",
                                                   m_halfMode, false);
                        m_session->MatchSolverInfo("ModeType", "MultipleModes",
                                                   m_multipleModes, false);
                    }

                    // Stability Analysis flags
                    if (m_session->DefinesSolverInfo("ModeType"))
                    {
                        if(m_singleMode)
                        {
                            m_npointsZ = 2;
                        }
                        else if(m_halfMode)
                        {
                            m_npointsZ = 1;
                        }
                        else if(m_multipleModes)
                        {
                            m_npointsZ = m_session->GetParameter("HomModesZ");
                        }
                        else
                        {
                            ASSERTL0(false, "SolverInfo ModeType not valid");
                        }
                    }
                    else
                    {
                        m_npointsZ = m_session->GetParameter("HomModesZ");
                    }
                }

                if ((HomoStr == "HOMOGENEOUS2D") || (HomoStr == "Homogeneous2D")
                    || (HomoStr == "2D") || (HomoStr == "Homo2D"))
                {
                    m_HomogeneousType = eHomogeneous2D;
                    m_session->LoadParameter("HomModesY",   m_npointsY);
                    m_session->LoadParameter("LY",          m_LhomY);
                    m_session->LoadParameter("HomModesZ",   m_npointsZ);
                    m_session->LoadParameter("LZ",          m_LhomZ);
                    m_HomoDirec       = 2;
                }

                if ((HomoStr == "HOMOGENEOUS3D") || (HomoStr == "Homogeneous3D")
                    || (HomoStr == "3D") || (HomoStr == "Homo3D"))
                {
                    m_HomogeneousType = eHomogeneous3D;
                    m_session->LoadParameter("HomModesY",   m_npointsY);
                    m_session->LoadParameter("LY",          m_LhomY);
                    m_session->LoadParameter("HomModesZ",   m_npointsZ);
                    m_session->LoadParameter("LZ",          m_LhomZ);
                    m_HomoDirec       = 2;
                }

                m_session->MatchSolverInfo("USEFFT", "FFTW", m_useFFT, false);

                m_session->MatchSolverInfo("DEALIASING", "True",
                                           m_homogen_dealiasing, false);
            }
            else
            {
                // set to default value so can use to identify 2d or 3D
                // (homogeneous) expansions
                m_npointsZ = 1;
            }

            m_session->MatchSolverInfo("SPECTRALHPDEALIASING", "True",
                                       m_specHP_dealiasing, false);
            if (m_specHP_dealiasing == false)
            {
                m_session->MatchSolverInfo("SPECTRALHPDEALIASING", "On",
                                           m_specHP_dealiasing, false);
            }

            // Options to determine type of projection from file or directly
            // from constructor
            if (m_session->DefinesSolverInfo("PROJECTION"))
            {
                std::string ProjectStr = m_session->GetSolverInfo("PROJECTION");

                if ((ProjectStr == "Continuous") || (ProjectStr == "Galerkin") ||
                   (ProjectStr == "CONTINUOUS") || (ProjectStr == "GALERKIN"))
                {
                    m_projectionType = MultiRegions::eGalerkin;
                }
                else if ((ProjectStr == "MixedCGDG") ||
                        (ProjectStr == "Mixed_CG_Discontinuous"))
                {
                    m_projectionType = MultiRegions::eMixed_CG_Discontinuous;
                }
                else if(ProjectStr == "DisContinuous")
                {
                    m_projectionType = MultiRegions::eDiscontinuous;
                }
                else
                {
                    ASSERTL0(false,"PROJECTION value not recognised");
                }
            }
            else
            {
                cerr << "Projection type not specified in SOLVERINFO,"
                    "defaulting to continuous Galerkin" << endl;
                m_projectionType = MultiRegions::eGalerkin;
            }

            // Enforce singularity check for some problems
            m_checkIfSystemSingular = v_GetSystemSingularChecks();

            int i;
            int nvariables = m_session->GetVariables().size();
            bool DeclareCoeffPhysArrays = true;


            m_fields   = Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);
            m_spacedim = m_graph->GetSpaceDimension()+m_HomoDirec;
            m_expdim   = m_graph->GetMeshDimension();

            /// Continuous field
            if (m_projectionType == MultiRegions::eGalerkin ||
               m_projectionType == MultiRegions::eMixed_CG_Discontinuous)
            {
                switch(m_expdim)
                {
                case 1:
                    {
                        if (m_HomogeneousType == eHomogeneous2D
                            || m_HomogeneousType == eHomogeneous3D)
                        {
                            const LibUtilities::PointsKey PkeyY
                                (m_npointsY, LibUtilities::eFourierEvenlySpaced);
                            const LibUtilities::BasisKey  BkeyY
                                (LibUtilities::eFourier, m_npointsY, PkeyY);
                            const LibUtilities::PointsKey PkeyZ
                                (m_npointsZ, LibUtilities::eFourierEvenlySpaced);
                            const LibUtilities::BasisKey
                                BkeyZ(LibUtilities::eFourier, m_npointsZ, PkeyZ);

                            for (i = 0; i < m_fields.size(); i++)
                            {
                                m_fields[i] = MemoryManager<MultiRegions
                                    ::ContField3DHomogeneous2D>
                                        ::AllocateSharedPtr(
                                            m_session, BkeyY, BkeyZ, m_LhomY,
                                            m_LhomZ, m_useFFT,
                                            m_homogen_dealiasing, m_graph,
                                            m_session->GetVariable(i));
                            }
                        }
                        else
                        {
                            for (i = 0; i < m_fields.size(); i++)
                            {
                                m_fields[i] = MemoryManager
                                <MultiRegions::ContField1D>::
                                    AllocateSharedPtr(
                                        m_session, m_graph,
                                        m_session->GetVariable(i));
                            }
                        }
                        break;
                    }
                case 2:
                    {
                        if (m_HomogeneousType == eHomogeneous1D)
                        {
                            // Fourier single mode stability analysis
                            if (m_singleMode)
                            {
                                const LibUtilities::PointsKey PkeyZ(
                                    m_npointsZ,
                                    LibUtilities::eFourierSingleModeSpaced);

                                const LibUtilities::BasisKey  BkeyZ(
                                    LibUtilities::eFourierSingleMode,
                                    m_npointsZ,
                                    PkeyZ);

                                for(i = 0; i < m_fields.size(); i++)
                                {
                                    m_fields[i] = MemoryManager<MultiRegions
                                        ::ContField3DHomogeneous1D>
                                            ::AllocateSharedPtr(
                                                m_session, BkeyZ, m_LhomZ,
                                                m_useFFT, m_homogen_dealiasing,
                                                m_graph,
                                                m_session->GetVariable(i),
                                                m_checkIfSystemSingular[i]);
                                }
                            }
                            // Half mode stability analysis
                            else if(m_halfMode)
                            {
                                const LibUtilities::PointsKey PkeyZ(
                                    m_npointsZ,
                                    LibUtilities::eFourierSingleModeSpaced);

                                const LibUtilities::BasisKey  BkeyZR(
                                    LibUtilities::eFourierHalfModeRe,
                                    m_npointsZ, PkeyZ);

                                const LibUtilities::BasisKey  BkeyZI(
                                    LibUtilities::eFourierHalfModeIm,
                                    m_npointsZ, PkeyZ);


                                for (i = 0; i < m_fields.size(); i++)
                                {
                                    if(m_session->GetVariable(i).compare("w")
                                            == 0)
                                    {
                                        m_fields[i] = MemoryManager<MultiRegions
                                            ::ContField3DHomogeneous1D>
                                                ::AllocateSharedPtr(
                                                    m_session, BkeyZI, m_LhomZ,
                                                    m_useFFT,
                                                    m_homogen_dealiasing,
                                                    m_graph,
                                                    m_session->GetVariable(i),
                                                    m_checkIfSystemSingular[i]);
                                    }
                                    else
                                    {
                                        m_fields[i] = MemoryManager<MultiRegions
                                            ::ContField3DHomogeneous1D>
                                                ::AllocateSharedPtr(
                                                    m_session, BkeyZR, m_LhomZ,
                                                    m_useFFT, m_homogen_dealiasing,
                                                    m_graph,
                                                    m_session->GetVariable(i),
                                                    m_checkIfSystemSingular[i]);
                                    }


                                }
                            }
                            // Normal homogeneous 1D
                            else
                            {
                                const LibUtilities::PointsKey PkeyZ(
                                    m_npointsZ,
                                    LibUtilities::eFourierEvenlySpaced);
                                const LibUtilities::BasisKey  BkeyZ(
                                    LibUtilities::eFourier, m_npointsZ, PkeyZ);

                                for (i = 0; i < m_fields.size(); i++)
                                {
                                    m_fields[i] = MemoryManager<MultiRegions
                                        ::ContField3DHomogeneous1D>
                                            ::AllocateSharedPtr(
                                                m_session, BkeyZ, m_LhomZ,
                                                m_useFFT, m_homogen_dealiasing,
                                                m_graph,
                                                m_session->GetVariable(i),
                                                m_checkIfSystemSingular[i]);
                                }
                            }
                        }
                        else
                        {
                            i = 0;
                            MultiRegions::ContField2DSharedPtr firstfield;
                            firstfield = MemoryManager<MultiRegions::
                                ContField2D>::AllocateSharedPtr(
                                    m_session, m_graph,
                                    m_session->GetVariable(i),
                                    DeclareCoeffPhysArrays,
                                    m_checkIfSystemSingular[0]);
                            m_fields[0] = firstfield;
                            for (i = 1; i < m_fields.size(); i++)
                            {
                                if (m_graph->
                                      SameExpansions(m_session->GetVariable(0),
                                                     m_session->GetVariable(i)))
                                {
                                    m_fields[i] = MemoryManager<MultiRegions::
                                        ContField2D>::AllocateSharedPtr(
                                            *firstfield, m_graph,
                                            m_session->GetVariable(i),
                                            DeclareCoeffPhysArrays,
                                            m_checkIfSystemSingular[i]);
                                }
                                else
                                {
                                    m_fields[i] = MemoryManager<MultiRegions
                                        ::ContField2D>::AllocateSharedPtr(
                                            m_session, m_graph,
                                            m_session->GetVariable(i),
                                            DeclareCoeffPhysArrays,
                                            m_checkIfSystemSingular[i]);
                                }
                            }

                            if (m_projectionType ==
                               MultiRegions::eMixed_CG_Discontinuous)
                            {
                                /// Setting up the normals
                                m_traceNormals =
                                    Array<OneD, Array<OneD, NekDouble> >
                                                                   (m_spacedim);

                                for (i = 0; i < m_spacedim; ++i)
                                {
                                    m_traceNormals[i] = Array<OneD, NekDouble>
                                                            (GetTraceNpoints());
                                }

                                m_fields[0]->GetTrace()->
                                    GetNormals(m_traceNormals);
                            }

                        }

                        break;
                    }
                case 3:
                    {
                        i = 0;
                        MultiRegions::ContField3DSharedPtr firstfield =
                            MemoryManager<MultiRegions::ContField3D>
                            ::AllocateSharedPtr(m_session, m_graph,
                                                m_session->GetVariable(i),
                                                m_checkIfSystemSingular[i]);

                        m_fields[0] = firstfield;
                        for (i = 1; i < m_fields.size(); i++)
                        {
                            if(m_graph->SameExpansions(
                                        m_session->GetVariable(0),
                                        m_session->GetVariable(i)))
                            {
                                m_fields[i] = MemoryManager<MultiRegions
                                    ::ContField3D>::AllocateSharedPtr(
                                        *firstfield, m_graph,
                                        m_session->GetVariable(i),
                                        m_checkIfSystemSingular[i]);
                            }
                            else
                            {
                                m_fields[i] = MemoryManager<MultiRegions
                                    ::ContField3D>::AllocateSharedPtr(
                                        m_session, m_graph,
                                        m_session->GetVariable(i),
                                        m_checkIfSystemSingular[i]);
                            }
                        }

                        if (m_projectionType ==
                           MultiRegions::eMixed_CG_Discontinuous)
                        {
                            /// Setting up the normals
                            m_traceNormals =
                                Array<OneD, Array<OneD, NekDouble> >
                                                                (m_spacedim);
                            for(i = 0; i < m_spacedim; ++i)
                            {
                                m_traceNormals[i] =
                                    Array<OneD, NekDouble> (GetTraceNpoints());
                            }

                            m_fields[0]->GetTrace()->GetNormals(m_traceNormals);
                            // Call the trace on all fields to ensure DG setup.
                            for(i = 1; i < m_fields.size(); ++i)
                            {
                                m_fields[i]->GetTrace();
                            }
                        }
                        break;
                    }
                default:
                    ASSERTL0(false,"Expansion dimension not recognised");
                    break;
                }
            }
            // Discontinuous field
            else
            {
                switch(m_expdim)
                {
                case 1:
                    {
                        if (m_HomogeneousType == eHomogeneous2D
                            || m_HomogeneousType == eHomogeneous3D)
                        {
                            const LibUtilities::PointsKey PkeyY(
                                m_npointsY, LibUtilities::eFourierEvenlySpaced);
                            const LibUtilities::BasisKey  BkeyY(
                                LibUtilities::eFourier, m_npointsY, PkeyY);
                            const LibUtilities::PointsKey PkeyZ(
                                m_npointsZ, LibUtilities::eFourierEvenlySpaced);
                            const LibUtilities::BasisKey  BkeyZ(
                                LibUtilities::eFourier, m_npointsZ, PkeyZ);

                            for (i = 0; i < m_fields.size(); i++)
                            {
                                m_fields[i] = MemoryManager<MultiRegions
                                    ::DisContField3DHomogeneous2D>
                                        ::AllocateSharedPtr(
                                            m_session, BkeyY, BkeyZ, m_LhomY,
                                            m_LhomZ, m_useFFT,
                                            m_homogen_dealiasing, m_graph,
                                            m_session->GetVariable(i));
                            }
                        }
                        else
                        {
                            for (i = 0; i < m_fields.size(); i++)
                            {
                                m_fields[i] = MemoryManager<MultiRegions::
                                    DisContField1D>::AllocateSharedPtr(
                                        m_session, m_graph,
                                        m_session->GetVariable(i));
                            }
                        }

                        break;
                    }
                case 2:
                    {
                        if(m_HomogeneousType == eHomogeneous1D)
                        {
                            const LibUtilities::PointsKey PkeyZ(
                                m_npointsZ,LibUtilities::eFourierEvenlySpaced);
                            const LibUtilities::BasisKey BkeyZ(
                                LibUtilities::eFourier, m_npointsZ,PkeyZ);

                            for (i = 0; i < m_fields.size(); i++)
                            {
                                m_fields[i] = MemoryManager<MultiRegions
                                    ::DisContField3DHomogeneous1D>
                                        ::AllocateSharedPtr(
                                            m_session, BkeyZ, m_LhomZ, m_useFFT,
                                            m_homogen_dealiasing, m_graph,
                                            m_session->GetVariable(i));
                            }
                        }
                        else
                        {
                            for (i = 0; i < m_fields.size(); i++)
                            {
                                m_fields[i] = MemoryManager<MultiRegions::
                                    DisContField2D>::AllocateSharedPtr(
                                        m_session, m_graph,
                                        m_session->GetVariable(i));
                            }
                        }

                        break;
                    }
                case 3:
                    {
                        if (m_HomogeneousType == eHomogeneous3D)
                        {
                            ASSERTL0(false,
                              "3D fully periodic problems not implemented yet");
                        }
                        else
                        {
                            for (i = 0; i < m_fields.size(); i++)
                            {
                                m_fields[i] = MemoryManager<MultiRegions::
                                    DisContField3D>::AllocateSharedPtr(
                                        m_session, m_graph,
                                        m_session->GetVariable(i));
                            }
                        }
                        break;
                    }
                    default:
                        ASSERTL0(false, "Expansion dimension not recognised");
                        break;
                }

                // Setting up the normals
                m_traceNormals =
                    Array<OneD, Array<OneD, NekDouble> >(m_spacedim);

                for (i = 0; i < m_spacedim; ++i)
                {
                    m_traceNormals[i] =
                        Array<OneD, NekDouble> (GetTraceNpoints(), 0.0);
                }

                m_fields[0]->GetTrace()->GetNormals(m_traceNormals);
            }

            // Set Default Parameter
            m_session->LoadParameter("Time",          m_time,       0.0);
            m_session->LoadParameter("TimeStep",      m_timestep,   0.0);
            m_session->LoadParameter("NumSteps",      m_steps,      0);
            m_session->LoadParameter("IO_CheckSteps", m_checksteps, 0);
            m_session->LoadParameter("IO_CheckTime",  m_checktime,  0.0);
            m_session->LoadParameter("FinTime",       m_fintime,    0);
            m_session->LoadParameter("NumQuadPointsError",
                                     m_NumQuadPointsError, 0);

            // Check uniqueness of checkpoint output
            ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
                     (m_checktime >  0.0 && m_checksteps == 0) ||
                     (m_checktime == 0.0 && m_checksteps >  0),
                     "Only one of IO_CheckTime and IO_CheckSteps "
                     "should be set!");
                     
            m_nchk = 0;

            // Zero all physical fields initially
            ZeroPhysFields();
        }

        /**
         * @brief Destructor for class EquationSystem.
         */
        EquationSystem::~EquationSystem()
        {
            LibUtilities::NekManager<LocalRegions::MatrixKey,
                DNekScalMat, LocalRegions::MatrixKey::opLess>::ClearManager();
            LibUtilities::NekManager<LocalRegions::MatrixKey,
                DNekScalBlkMat, LocalRegions::MatrixKey::opLess>::ClearManager();
        }

        SessionFunctionSharedPtr EquationSystem::GetFunction(
            std::string name,
            const MultiRegions::ExpListSharedPtr &field,
            bool cache)
        {
            MultiRegions::ExpListSharedPtr vField = field;
            if (!field)
            {
                vField = m_fields[0];
            }

            if (cache)
            {
                if ((m_sessionFunctions.find(name) == m_sessionFunctions.end())
                    || (m_sessionFunctions[name]->GetSession() != m_session)
                    || (m_sessionFunctions[name]->GetExpansion() != vField)
                )
                {
                    m_sessionFunctions[name] =
                        MemoryManager<SessionFunction>::AllocateSharedPtr(
                            m_session, vField, name, cache);
                }

                return m_sessionFunctions[name];
            }
            else
            {
                return SessionFunctionSharedPtr(
                        new SessionFunction(m_session,vField, name, cache));
            }
        }

        /**
         * If boundary conditions are time-dependent, they will be evaluated at
         * the time specified.
         * @param   time            The time at which to evaluate the BCs
         */
        void EquationSystem::SetBoundaryConditions(NekDouble time)
        {
            std::string varName;
            int nvariables = m_fields.size();
            for (int i = 0; i < nvariables; ++i)
            {
                varName = m_session->GetVariable(i);
                m_fields[i]->EvaluateBoundaryConditions(time, varName);
            }
        }

        /**
         * Compute the error in the L2-norm.
         * @param   field           The field to compare.
         * @param   exactsoln       The exact solution to compare with.
         * @param   Normalised      Normalise L2-error.
         * @returns                 Error in the L2-norm.
         */
        NekDouble EquationSystem::v_L2Error(
            unsigned int field,
            const Array<OneD, NekDouble> &exactsoln,
            bool Normalised)
        {
            NekDouble L2error = -1.0;

            if (m_NumQuadPointsError == 0)
            {
                if (m_fields[field]->GetPhysState() == false)
                {
                    m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
                                              m_fields[field]->UpdatePhys());
                }

                if (exactsoln.size())
                {
                    L2error = m_fields[field]->L2(
                            m_fields[field]->GetPhys(), exactsoln);
                }
                else if (m_session->DefinesFunction("ExactSolution"))
                {
                    Array<OneD, NekDouble>
                        exactsoln(m_fields[field]->GetNpoints());

                    GetFunction("ExactSolution")->Evaluate(
                            m_session->GetVariable(field), exactsoln, m_time);

                    L2error = m_fields[field]->L2(
                            m_fields[field]->GetPhys(), exactsoln);
                }
                else
                {
                    L2error = m_fields[field]->L2(m_fields[field]->GetPhys());
                }

                if (Normalised == true)
                {
                    Array<OneD, NekDouble> one(m_fields[field]->GetNpoints(),
                                               1.0);

                    NekDouble Vol = m_fields[field]->PhysIntegral(one);
                    m_comm->AllReduce(Vol, LibUtilities::ReduceSum);

                    L2error = sqrt(L2error*L2error/Vol);
                }
            }
            else
            {
                Array<OneD,NekDouble> L2INF(2);
                L2INF = ErrorExtraPoints(field);
                L2error = L2INF[0];
            }
            return L2error;
        }

        /**
         * Compute the error in the L_inf-norm
         * @param   field           The field to compare.
         * @param   exactsoln       The exact solution to compare with.
         * @returns                 Error in the L_inft-norm.
         */
        NekDouble EquationSystem::v_LinfError (
                unsigned int field,
                const Array<OneD, NekDouble> &exactsoln)
        {
            NekDouble Linferror = -1.0;

            if (m_NumQuadPointsError == 0)
            {
                if (m_fields[field]->GetPhysState() == false)
                {
                    m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
                                              m_fields[field]->UpdatePhys());
                }

                if (exactsoln.size())
                {
                    Linferror = m_fields[field]->Linf(
                            m_fields[field]->GetPhys(), exactsoln);
                }
                else if (m_session->DefinesFunction("ExactSolution"))
                {
                    Array<OneD, NekDouble>
                        exactsoln(m_fields[field]->GetNpoints());

                    GetFunction("ExactSolution")->Evaluate(
                            m_session->GetVariable(field), exactsoln, m_time);

                    Linferror = m_fields[field]->Linf(
                            m_fields[field]->GetPhys(), exactsoln);
                }
                else
                {
                    Linferror = m_fields[field]->Linf(
                            m_fields[field]->GetPhys());
                }
            }
            else
            {
                Array<OneD,NekDouble> L2INF(2);
                L2INF = ErrorExtraPoints(field);
                Linferror = L2INF[1];
            }

            return Linferror;
        }

        /**
         * Compute the error in the L2-norm, L-inf for a larger number of
         * quadrature points.
         * @param   field              The field to compare.
         * @returns                    Error in the L2-norm and L-inf norm.
         */
        Array<OneD,NekDouble> EquationSystem::ErrorExtraPoints(
            unsigned int field)
        {
            int NumModes = GetNumExpModes();
            Array<OneD,NekDouble> L2INF(2);

            const LibUtilities::PointsKey PkeyT1(
                m_NumQuadPointsError,LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey PkeyT2(
                m_NumQuadPointsError, LibUtilities::eGaussRadauMAlpha1Beta0);
            const LibUtilities::PointsKey PkeyQ1(
                m_NumQuadPointsError, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey PkeyQ2(
                m_NumQuadPointsError, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::BasisKey  BkeyT1(
                LibUtilities::eModified_A,NumModes, PkeyT1);
            const LibUtilities::BasisKey  BkeyT2(
                LibUtilities::eModified_B, NumModes, PkeyT2);
            const LibUtilities::BasisKey  BkeyQ1(
                LibUtilities::eModified_A, NumModes, PkeyQ1);
            const LibUtilities::BasisKey  BkeyQ2(
                LibUtilities::eModified_A, NumModes, PkeyQ2);

            MultiRegions::ExpList2DSharedPtr ErrorExp =
                MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(
                        m_session, BkeyT1, BkeyT2, BkeyQ1, BkeyQ2, m_graph);

            int ErrorCoordim = ErrorExp->GetCoordim(0);
            int ErrorNq      = ErrorExp->GetTotPoints();

            Array<OneD,NekDouble> ErrorXc0(ErrorNq, 0.0);
            Array<OneD,NekDouble> ErrorXc1(ErrorNq, 0.0);
            Array<OneD,NekDouble> ErrorXc2(ErrorNq, 0.0);

            switch(ErrorCoordim)
            {
                case 1:
                    ErrorExp->GetCoords(ErrorXc0);
                    break;
                case 2:
                    ErrorExp->GetCoords(ErrorXc0, ErrorXc1);
                    break;
                case 3:
                    ErrorExp->GetCoords(ErrorXc0, ErrorXc1, ErrorXc2);
                    break;
            }
            LibUtilities::EquationSharedPtr exSol =
                m_session->GetFunction("ExactSolution", field);

            // Evaluate the exact solution
            Array<OneD,NekDouble> ErrorSol(ErrorNq);

            exSol->Evaluate(ErrorXc0,ErrorXc1,ErrorXc2,m_time,ErrorSol);

            // Calcualte spectral/hp approximation on the quadrature points
            // of this new expansion basis
            ErrorExp->BwdTrans_IterPerExp(m_fields[field]->GetCoeffs(),
                                          ErrorExp->UpdatePhys());

            L2INF[0] = ErrorExp->L2  (ErrorExp->GetPhys(), ErrorSol);
            L2INF[1] = ErrorExp->Linf(ErrorExp->GetPhys(), ErrorSol);

            return L2INF;
        }


        /**
         * Set the physical fields based on a restart file, or a function
         * describing the initial condition given in the session.
         * @param  initialtime           Time at which to evaluate the function.
         * @param  dumpInitialConditions Write the initial condition to file?
         */
        void EquationSystem::v_SetInitialConditions(NekDouble initialtime,
                                                    bool dumpInitialConditions,
                                                    const int domain)
        {
            boost::ignore_unused(initialtime);

            if (m_session->GetComm()->GetRank() == 0)
            {
                cout << "Initial Conditions:" << endl;
            }

            if (m_session->DefinesFunction("InitialConditions"))
            {
                GetFunction("InitialConditions")->Evaluate(
                        m_session->GetVariables(), m_fields, m_time, domain);
                // Enforce C0 Continutiy of initial condiiton
		if((m_projectionType == MultiRegions::eGalerkin)||
		    (m_projectionType == MultiRegions::eMixed_CG_Discontinuous))
		{
		     for (int i = 0; i < m_fields.size(); ++i)
		     {
                         m_fields[i]->LocalToGlobal();
                         m_fields[i]->GlobalToLocal();
		         m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
		                               m_fields[i]->UpdatePhys());
		      }
                } 

                if (m_session->GetComm()->GetRank() == 0)
                {

                    for (int i = 0; i < m_fields.size(); ++i)
                    {
                        std::string varName = m_session->GetVariable(i);
                        cout << "  - Field " << varName << ": "
                             << GetFunction("InitialConditions")->Describe(varName, domain)
                             << endl;
                    }
                }
            }
            else
            {
                int nq = m_fields[0]->GetNpoints();
                for (int i = 0; i < m_fields.size(); i++)
                {
                    Vmath::Zero(nq, m_fields[i]->UpdatePhys(), 1);
                    m_fields[i]->SetPhysState(true);
                    Vmath::Zero(m_fields[i]->GetNcoeffs(),
                                m_fields[i]->UpdateCoeffs(), 1);
                    if (m_session->GetComm()->GetRank() == 0)
                    {
                        cout << "  - Field "    << m_session->GetVariable(i)
                             << ": 0 (default)" << endl;
                    }
                }

            }

            if (dumpInitialConditions && m_checksteps)
            {
                Checkpoint_Output(m_nchk);
                m_nchk++;
            }
        }

        void EquationSystem::v_EvaluateExactSolution(
            unsigned int field,
            Array<OneD, NekDouble> &outfield,
            const NekDouble time)
        {
            ASSERTL0 (outfield.size() == m_fields[field]->GetNpoints(),
                      "ExactSolution array size mismatch.");
            Vmath::Zero(outfield.size(), outfield, 1);
            if (m_session->DefinesFunction("ExactSolution"))
            {
                GetFunction("ExactSolution")->Evaluate(
                        m_session->GetVariable(field), outfield, time);
            }
        }


        /**
         * By default, nothing needs initialising at the EquationSystem level.
         */
        void EquationSystem::v_DoInitialise()
        {

        }

        /**
         *
         */
        void EquationSystem::v_DoSolve()
        {

        }

        /**
         * Virtual function to define if operator in DoSolve is
         * negated with regard to the strong form. This is currently
         * only used in Arnoldi solves. Default is false.
         */
        bool EquationSystem::v_NegatedOp(void)
        {
            return false;
        }

        /**
         *
         */
        void EquationSystem::v_TransCoeffToPhys()
        {

        }

        /**
         *
         */
        void EquationSystem::v_TransPhysToCoeff()
        {

        }


        /// Virtual function for generating summary information.
        void EquationSystem::v_GenerateSummary(SummaryList& l)
        {
            SessionSummary(l);
        }


        /**
         * Write the field data to file. The file is named according to the
         * session name with the extension .fld appended.
         */
        void EquationSystem::v_Output(void)
        {
            WriteFld(m_sessionName + ".fld");
        }

        /**
         * Zero the physical fields.
         */
        void EquationSystem::ZeroPhysFields(void)
        {
            for (int i = 0; i < m_fields.size(); i++)
            {
                Vmath::Zero(m_fields[i]->GetNpoints(),
                            m_fields[i]->UpdatePhys(),1);
            }
        }

        /**
         * FwdTrans the m_fields members
         */
        void EquationSystem::FwdTransFields(void)
        {
            for (int i = 0; i < m_fields.size(); i++)
            {
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                      m_fields[i]->UpdateCoeffs());
                m_fields[i]->SetPhysState(false);
            }
        }

        /**
         * Write the n-th checkpoint file.
         * @param   n   The index of the checkpoint file.
         */
        void EquationSystem::Checkpoint_Output(const int n)
        {
            std::string outname =  m_sessionName +  "_" +
                boost::lexical_cast<std::string>(n);
            WriteFld(outname + ".chk");
        }

        /**
         * Write the n-th checkpoint file.
         * @param   n   The index of the checkpoint file.
         */
        void EquationSystem::Checkpoint_Output(
            const int n,
            MultiRegions::ExpListSharedPtr &field,
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string> &variables)
        {
            std::string outname =  m_sessionName +  "_" +
                boost::lexical_cast<std::string>(n);
            WriteFld(outname, field, fieldcoeffs, variables);
        }

        /**
         * Write the n-th base flow into a .chk file
         * @param   n   The index of the base flow file.
         */
        void EquationSystem::Checkpoint_BaseFlow(const int n)
        {
            std::string outname =  m_sessionName +  "_BaseFlow_" +
                boost::lexical_cast<std::string>(n);

            WriteFld(outname + ".chk");
        }

        /**
         * Writes the field data to a file with the given filename.
         * @param   outname     Filename to write to.
         */
        void EquationSystem::WriteFld(const std::string &outname)
        {
            std::vector<Array<OneD, NekDouble> > fieldcoeffs(
                m_fields.size());
            std::vector<std::string> variables(m_fields.size());

            for (int i = 0; i < m_fields.size(); ++i)
            {
                if (m_fields[i]->GetNcoeffs() == m_fields[0]->GetNcoeffs())
                {
                    fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
                }
                else
                {
                    fieldcoeffs[i] = Array<OneD,NekDouble>(m_fields[0]->
                                                           GetNcoeffs());
                    m_fields[0]->ExtractCoeffsToCoeffs(m_fields[i],
                                                       m_fields[i]->GetCoeffs(),
                                                       fieldcoeffs[i]);
                }
                variables[i] = m_boundaryConditions->GetVariable(i);
            }

            ExtraFldOutput(fieldcoeffs, variables);

            WriteFld(outname, m_fields[0], fieldcoeffs, variables);
        }



        /**
         * Writes the field data to a file with the given filename.
         * @param   outname         Filename to write to.
         * @param   field           ExpList on which data is based.
         * @param   fieldcoeffs     An array of array of expansion coefficients.
         * @param   variables       An array of variable names.
         */
        void EquationSystem::WriteFld(
            const std::string                    &outname,
            MultiRegions::ExpListSharedPtr       &field,
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string>             &variables)
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = field->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

            // Copy Data into FieldData and set variable
            for(int j = 0; j < fieldcoeffs.size(); ++j)
            {
                for(int i = 0; i < FieldDef.size(); ++i)
                {
                    // Could do a search here to find correct variable
                    FieldDef[i]->m_fields.push_back(variables[j]);
                    field->AppendFieldData(FieldDef[i], FieldData[i],
                                           fieldcoeffs[j]);
                }
            }

            // Update time in field info if required
            if(m_fieldMetaDataMap.find("Time") != m_fieldMetaDataMap.end())
            {
                m_fieldMetaDataMap["Time"] =
                        boost::lexical_cast<std::string>(m_time);
            }

            // Update step in field info if required
            if(m_fieldMetaDataMap.find("ChkFileNum") != m_fieldMetaDataMap.end())
            {
                m_fieldMetaDataMap["ChkFileNum"] =
                        boost::lexical_cast<std::string>(m_nchk);
            }

            // If necessary, add mapping information to metadata
            //      and output mapping coordinates
            Array<OneD, MultiRegions::ExpListSharedPtr> fields(1);
            fields[0] = field;
            GlobalMapping::MappingSharedPtr mapping =
                    GlobalMapping::Mapping::Load(m_session, fields);
            LibUtilities::FieldMetaDataMap fieldMetaDataMap(m_fieldMetaDataMap);
            mapping->Output( fieldMetaDataMap, outname);

#ifdef NEKTAR_DISABLE_BACKUPS
            bool backup = false;
#else
            bool backup = true;
#endif

            m_fld->Write(outname, FieldDef, FieldData, fieldMetaDataMap, backup);
        }


        /**
         * Import field from infile and load into \a m_fields. This routine will
         * also perform a \a BwdTrans to ensure data is in both the physical and
         * coefficient storage.
         * @param   infile  Filename to read.
         */
        void EquationSystem::ImportFld(
            const std::string &infile,
            Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;
            LibUtilities::FieldIOSharedPtr field_fld =
                LibUtilities::FieldIO::CreateForFile(m_session, infile);
            field_fld->Import(infile,FieldDef,FieldData);

            // Copy FieldData into m_fields
            for(int j = 0; j < pFields.size(); ++j)
            {
                Vmath::Zero(pFields[j]->GetNcoeffs(),
                            pFields[j]->UpdateCoeffs(),1);

                for(int i = 0; i < FieldDef.size(); ++i)
                {
                    ASSERTL1(FieldDef[i]->m_fields[j] ==
                             m_session->GetVariable(j),
                             std::string("Order of ") + infile
                             + std::string(" data and that defined in "
                                           "m_boundaryconditions differs"));

                    pFields[j]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                                    FieldDef[i]->m_fields[j],
                                                    pFields[j]->UpdateCoeffs());
                }
                pFields[j]->BwdTrans(pFields[j]->GetCoeffs(),
                                     pFields[j]->UpdatePhys());
            }
        }



        /**
         * Import field from infile and load into \a m_fields. This routine will
         * also perform a \a BwdTrans to ensure data is in both the physical and
         * coefficient storage.
         * @param   infile  Filename to read.
         * If optionan \a ndomains is specified it assumes we loop over nodmains
         * for each nvariables.
         */
        void EquationSystem::ImportFldToMultiDomains(
                    const std::string &infile,
                    Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                    const int ndomains)
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;

            LibUtilities::Import(infile,FieldDef,FieldData);

            int nvariables = GetNvariables();

            ASSERTL0(ndomains*nvariables == pFields.size(),
                "Number of fields does not match the number of variables and domains");

            // Copy FieldData into m_fields
            for(int j = 0; j < ndomains; ++j)
            {
                for(int i = 0; i < nvariables; ++i)
                {
                    Vmath::Zero(pFields[j*nvariables+i]->GetNcoeffs(),
                            pFields[j*nvariables+i]->UpdateCoeffs(),1);

                    for(int n = 0; n < FieldDef.size(); ++n)
                    {
                        ASSERTL1(FieldDef[n]->m_fields[i] == m_session->GetVariable(i),
                                 std::string("Order of ") + infile
                                 + std::string(" data and that defined in "
                                               "m_boundaryconditions differs"));

                        pFields[j*nvariables+i]->ExtractDataToCoeffs(
                                FieldDef[n], FieldData[n],
                                FieldDef[n]->m_fields[i],
                                pFields[j*nvariables+i]->UpdateCoeffs());
                    }
                    pFields[j*nvariables+i]->BwdTrans(
                                pFields[j*nvariables+i]->GetCoeffs(),
                                pFields[j*nvariables+i]->UpdatePhys());
                }
            }
        }

        /**
         * Import field from infile and load into \a pField. This routine will
         * also perform a \a BwdTrans to ensure data is in both the physical and
         * coefficient storage.
         */
        void EquationSystem::ImportFld(
            const std::string &infile,
            MultiRegions::ExpListSharedPtr &pField,
            std::string &pFieldName)
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;

            LibUtilities::FieldIOSharedPtr field_fld =
                LibUtilities::FieldIO::CreateForFile(m_session, infile);
            field_fld->Import(infile,FieldDef,FieldData);
            int idx = -1;

            Vmath::Zero(pField->GetNcoeffs(),pField->UpdateCoeffs(),1);

            for(int i = 0; i < FieldDef.size(); ++i)
            {
                // find the index of the required field in the file.
                for(int j = 0; j < FieldData.size(); ++j)
                {
                    if (FieldDef[i]->m_fields[j] == pFieldName)
                    {
                        idx = j;
                    }
                }
                ASSERTL1(idx >= 0, "Field " + pFieldName + " not found.");

                pField->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                            FieldDef[i]->m_fields[idx],
                                            pField->UpdateCoeffs());
            }
            pField->BwdTrans(pField->GetCoeffs(), pField->UpdatePhys());
        }

        /**
         * Import field from infile and load into the array \a coeffs.
         *
         * @param infile   Filename to read.
         * @param fieldStr an array of string identifying fields to be imported
         * @param coeffs   array of array of coefficients to store imported data
         */
        void EquationSystem::ImportFld(
            const std::string &infile,
            std::vector< std::string> &fieldStr,
            Array<OneD, Array<OneD, NekDouble> > &coeffs)
        {

            ASSERTL0(fieldStr.size() <= coeffs.size(),
                     "length of fieldstr should be the same as pFields");

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;

            LibUtilities::FieldIOSharedPtr field_fld =
                LibUtilities::FieldIO::CreateForFile(m_session, infile);
            field_fld->Import(infile,FieldDef,FieldData);

            // Copy FieldData into m_fields
            for(int j = 0; j < fieldStr.size(); ++j)
            {
                Vmath::Zero(coeffs[j].size(),coeffs[j],1);
                for(int i = 0; i < FieldDef.size(); ++i)
                {
                    m_fields[0]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                                     fieldStr[j], coeffs[j]);
                }
            }
        }

        /**
         * Write out a summary of the session data.
         * @param   out         Output stream to write data to.
         */
        void EquationSystem::SessionSummary(SummaryList& s)
        {
            AddSummaryItem(s, "EquationType",
                            m_session->GetSolverInfo("EQTYPE"));
            AddSummaryItem(s, "Session Name", m_sessionName);
            AddSummaryItem(s, "Spatial Dim.", m_spacedim);
            AddSummaryItem(s, "Max SEM Exp. Order",
                            m_fields[0]->EvalBasisNumModesMax());

            if (m_session->GetComm()->GetSize() > 1)
            {
                AddSummaryItem(s, "Num. Processes",
                        m_session->GetComm()->GetSize());
            }

            if(m_HomogeneousType == eHomogeneous1D)
            {
                AddSummaryItem(s, "Quasi-3D", "Homogeneous in z-direction");
                AddSummaryItem(s, "Expansion Dim.", m_expdim + 1);
                AddSummaryItem(s, "Num. Hom. Modes (z)", m_npointsZ);
                AddSummaryItem(s, "Hom. length (LZ)", m_LhomZ);
                AddSummaryItem(s, "FFT Type", m_useFFT ? "FFTW" : "MVM");
                if (m_halfMode)
                {
                    AddSummaryItem(s, "ModeType", "Half Mode");
                }
                else if (m_singleMode)
                {
                    AddSummaryItem(s, "ModeType", "Single Mode");
                }
                else if (m_multipleModes)
                {
                    AddSummaryItem(s, "ModeType", "Multiple Modes");
                }
            }
            else if(m_HomogeneousType == eHomogeneous2D)
            {
                AddSummaryItem(s, "Quasi-3D", "Homogeneous in yz-plane");
                AddSummaryItem(s, "Expansion Dim.", m_expdim + 2);
                AddSummaryItem(s, "Num. Hom. Modes (y)", m_npointsY);
                AddSummaryItem(s, "Num. Hom. Modes (z)", m_npointsZ);
                AddSummaryItem(s, "Hom. length (LY)", m_LhomY);
                AddSummaryItem(s, "Hom. length (LZ)", m_LhomZ);
                AddSummaryItem(s, "FFT Type", m_useFFT ? "FFTW" : "MVM");
            }
            else
            {
                AddSummaryItem(s, "Expansion Dim.", m_expdim);
            }

            if (m_session->DefinesSolverInfo("UpwindType"))
            {
                AddSummaryItem(s, "Riemann Solver",
                                  m_session->GetSolverInfo("UpwindType"));
            }

            if (m_session->DefinesSolverInfo("AdvectionType"))
            {
                std::string AdvectionType;
                AdvectionType = m_session->GetSolverInfo("AdvectionType");
                AddSummaryItem(s, "Advection Type", GetAdvectionFactory().
                               GetClassDescription(AdvectionType));
            }

            if (m_projectionType == MultiRegions::eGalerkin)
            {
                AddSummaryItem(s, "Projection Type", "Continuous Galerkin");
            }
            else if (m_projectionType == MultiRegions::eDiscontinuous)
            {
                AddSummaryItem(s, "Projection Type", "Discontinuous Galerkin");
            }
            else if (m_projectionType == MultiRegions::eMixed_CG_Discontinuous)
            {
                AddSummaryItem(s, "Projection Type",
                                "Mixed Continuous Galerkin and Discontinuous");
            }

            if (m_session->DefinesSolverInfo("DiffusionType"))
            {
                std::string DiffusionType;
                DiffusionType = m_session->GetSolverInfo("DiffusionType");
                AddSummaryItem(s, "Diffusion Type", GetDiffusionFactory().
                               GetClassDescription(DiffusionType));
            }
        }

        Array<OneD, bool> EquationSystem::v_GetSystemSingularChecks()
        {
            return Array<OneD, bool>(m_session->GetVariables().size(), false);
        }

        MultiRegions::ExpListSharedPtr EquationSystem::v_GetPressure()
        {
            ASSERTL0(false, "This function is not valid for the Base class");
            MultiRegions::ExpListSharedPtr null;
            return null;
        }

        void EquationSystem::v_ExtraFldOutput(
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string>             &variables)
        {
            boost::ignore_unused(fieldcoeffs, variables);
        }

    }
}
