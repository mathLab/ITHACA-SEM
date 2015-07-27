///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.cpp
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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/HexExp.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

namespace Nektar
{
    string CompressibleFlowSystem::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "CompressibleFlowSystem",
            CompressibleFlowSystem::create,
            "Auxiliary functions for the compressible flow system.");

    CompressibleFlowSystem::CompressibleFlowSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void CompressibleFlowSystem::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");

        // Do not forwards transform initial condition
        m_homoInitialFwd = false;

        // Set up locations of velocity vector.
        m_vecLocs = Array<OneD, Array<OneD, NekDouble> >(1);
        m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_vecLocs[0][i] = 1 + i;
        }

        // Get gamma parameter from session file.
        ASSERTL0(m_session->DefinesParameter("Gamma"),
                 "Compressible flow sessions must define a Gamma parameter.");
        m_session->LoadParameter("Gamma", m_gamma, 1.4);

        // Get E0 parameter from session file.
        ASSERTL0(m_session->DefinesParameter("pInf"),
                 "Compressible flow sessions must define a pInf parameter.");
        m_session->LoadParameter("pInf", m_pInf, 101325);

        // Get rhoInf parameter from session file.
        ASSERTL0(m_session->DefinesParameter("rhoInf"),
                 "Compressible flow sessions must define a rhoInf parameter.");
        m_session->LoadParameter("rhoInf", m_rhoInf, 1.225);

        // Get uInf parameter from session file.
        ASSERTL0(m_session->DefinesParameter("uInf"),
                 "Compressible flow sessions must define a uInf parameter.");
        m_session->LoadParameter("uInf", m_uInf, 0.1);

        m_UInf = m_uInf;

        // Get vInf parameter from session file.
        if (m_spacedim == 2 || m_spacedim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("vInf"),
                     "Compressible flow sessions must define a vInf parameter"
                     "for 2D/3D problems.");
            m_session->LoadParameter("vInf", m_vInf, 0.0);
            m_UInf = sqrt(m_uInf*m_uInf + m_vInf*m_vInf);
        }

        // Get wInf parameter from session file.
        if (m_spacedim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("wInf"),
                     "Compressible flow sessions must define a wInf parameter"
                     "for 3D problems.");
            m_session->LoadParameter("wInf", m_wInf, 0.0);
            m_UInf = sqrt(m_uInf*m_uInf + m_vInf*m_vInf + m_wInf*m_wInf);
        }

        m_session->LoadParameter ("GasConstant",   m_gasConstant,   287.058);
        m_session->LoadParameter ("Twall",         m_Twall,         300.15);
        m_session->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
        m_session->LoadParameter ("mu",            m_mu,            1.78e-05);
        m_session->LoadParameter ("Skappa",        m_Skappa,        -2.048);
        m_session->LoadParameter ("Kappa",         m_Kappa,         0.0);
        m_session->LoadParameter ("mu0",           m_mu0,           1.0);
        m_session->LoadParameter ("FL",            m_FacL,          0.0);
        m_session->LoadParameter ("FH",            m_FacH,          0.0);
        m_session->LoadParameter ("hFactor",       m_hFactor,       1.0);
        m_session->LoadParameter ("epsMax",        m_eps_max,       1.0);
        m_session->LoadParameter ("C1",            m_C1,            3.0);
        m_session->LoadParameter ("C2",            m_C2,            5.0);
        m_session->LoadSolverInfo("ShockCaptureType",
                                  m_shockCaptureType,    "Off");
        m_session->LoadParameter ("thermalConductivity",
                                  m_thermalConductivity, 0.0257);

        m_EqTypeStr = m_session->GetSolverInfo("EQTYPE");

        m_Cp      = m_gamma / (m_gamma - 1.0) * m_gasConstant;
        m_Prandtl = m_Cp * m_mu / m_thermalConductivity;

        m_session->LoadParameter("amplitude",      m_amplitude,      0.001);
        m_session->LoadParameter("omega",          m_omega,          1.0);
        m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 0.0);

        // Forcing terms for the sponge region
        m_forcing = SolverUtils::Forcing::Load(m_session, m_fields,
                                               m_fields.num_elements());

        // Loop over Boundary Regions for PressureOutflowFileBC
        int nvariables = m_fields.num_elements();
        Array<OneD, Array<OneD, NekDouble> > tmpStorage(nvariables);
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // PressureOutflowFile Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                "PressureOutflowFile")
            {
                int numBCPts = m_fields[0]->
                    GetBndCondExpansions()[n]->GetNpoints();
                m_pressureStorage = Array<OneD, NekDouble>(numBCPts, 0.0);
                for (int i = 0; i < nvariables; ++i)
                {
                    tmpStorage[i] = Array<OneD, NekDouble>(numBCPts, 0.0);

                    Vmath::Vcopy(
                        numBCPts,
                        m_fields[i]->GetBndCondExpansions()[n]->GetPhys(), 1,
                        tmpStorage[i], 1);
                }
                GetPressure(tmpStorage, m_pressureStorage);
            }
        }

        // Loop over Boundary Regions for PressureInflowFileBC
        m_fieldStorage = Array<OneD, Array<OneD, NekDouble> > (nvariables);
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // PressureInflowFile Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                "PressureInflowFile")
            {
                int numBCPts = m_fields[0]->
                    GetBndCondExpansions()[n]->GetNpoints();
                for (int i = 0; i < nvariables; ++i)
                {
                    m_fieldStorage[i] = Array<OneD, NekDouble>(numBCPts, 0.0);
                    Vmath::Vcopy(
                        numBCPts,
                        m_fields[i]->GetBndCondExpansions()[n]->GetPhys(), 1,
                        m_fieldStorage[i], 1);
                }
            }
        }

        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Continuous field
            case MultiRegions::eGalerkin:
            {
                ASSERTL0(false, "Continuous field not supported.");
                break;
            }
            // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                string advName, diffName, riemName;

                // Setting up advection and diffusion operators
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

                m_advection = SolverUtils::GetAdvectionFactory()
                                            .CreateInstance(advName, advName);
                m_diffusion = SolverUtils::GetDiffusionFactory()
                                            .CreateInstance(diffName, diffName);

                if (m_specHP_dealiasing)
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVectorDeAlias, this);
                    m_diffusion->SetFluxVectorNS(
                        &CompressibleFlowSystem::GetViscousFluxVectorDeAlias,
                        this);
                }
                else
                {
                    m_advection->SetFluxVector  (&CompressibleFlowSystem::
                                                  GetFluxVector, this);
                    m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                                                  GetViscousFluxVector, this);
                }

                if (m_shockCaptureType=="Smooth" && m_EqTypeStr=="EulerADCFE")
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVector, this);

                    m_diffusion->SetArtificialDiffusionVector(
                        &CompressibleFlowSystem::GetSmoothArtificialViscosity, this);
                }
                if (m_shockCaptureType=="NonSmooth" && m_EqTypeStr=="EulerADCFE")
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVector, this);

                    m_diffusion->SetArtificialDiffusionVector(
                        &CompressibleFlowSystem::GetArtificialDynamicViscosity, this);
                }

                // Setting up Riemann solver for advection operator
                m_session->LoadSolverInfo("UpwindType", riemName, "Average");

                m_riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                            .CreateInstance(riemName);

                // Setting up upwind solver for diffusion operator
                m_riemannSolverLDG = SolverUtils::GetRiemannSolverFactory()
                                                .CreateInstance("UpwindLDG");

                // Setting up parameters for advection operator Riemann solver
                m_riemannSolver->SetParam (
                    "gamma",   &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolver->SetAuxVec(
                    "vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
                m_riemannSolver->SetVector(
                    "N",       &CompressibleFlowSystem::GetNormals, this);

                // Setting up parameters for diffusion operator Riemann solver
                m_riemannSolverLDG->SetParam (
                    "gamma",   &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolverLDG->SetVector(
                    "vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
                m_riemannSolverLDG->SetVector(
                    "N",       &CompressibleFlowSystem::GetNormals, this);

                // Concluding initialisation of advection / diffusion operators
                m_advection->SetRiemannSolver   (m_riemannSolver);
                m_diffusion->SetRiemannSolver   (m_riemannSolverLDG);
                m_advection->InitObject         (m_session, m_fields);
                m_diffusion->InitObject         (m_session, m_fields);
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }
    }

    /**
     * @brief Destructor for CompressibleFlowSystem class.
     */
    CompressibleFlowSystem::~CompressibleFlowSystem()
    {

    }

    /**
     * @brief Print out a summary with some relevant information.
     */
    void CompressibleFlowSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        UnsteadySystem::v_GenerateSummary(s);
    }

    /**
     * @brief Set boundary conditions which can be: 
     * a) Wall and Symmerty BCs implemented at CompressibleFlowSystem level
     *          since they are compressible solver specific;
     * b) Time dependent BCs.
     * 
     * @param inarray: fields array.
     * @param time:    time.
     */
    void CompressibleFlowSystem::SetCommonBC(
                      const std::string &userDefStr,
                      const int n,
                      const NekDouble time,
                      int &cnt,
                      Array<OneD, Array<OneD, NekDouble> > &inarray)
    {
        std::string varName;
        int nvariables = m_fields.num_elements();

        if(!userDefStr.empty())
        {
            if(boost::iequals(userDefStr,"Wall"))
            {
                WallBC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"WallViscous") ||
                    boost::iequals(userDefStr,"WallAdiabatic"))
            {
                // Wall Boundary Condition
                WallViscousBC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"Symmetry"))
            {
                // Symmetric Boundary Condition
                SymmetryBC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"RiemannInvariant"))
            {
                // Riemann invariant characteristic Boundary Condition
                RiemannInvariantBC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"PressureOutflowNonReflective"))
            {
                // Pressure outflow non-reflective Boundary Condition
                PressureOutflowNonReflectiveBC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"PressureOutflow"))
            {
                // Pressure outflow Boundary Condition
                PressureOutflowBC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"PressureOutflowFile"))
            {
                // Pressure outflow Boundary Condition from file 
                PressureOutflowFileBC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"PressureInflowFile"))
            {
                // Pressure inflow Boundary Condition from file
                PressureInflowFileBC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"ExtrapOrder0"))
            {
                // Extrapolation of the data at the boundaries
                ExtrapOrder0BC(n, cnt, inarray);
            }
            else if(boost::iequals(userDefStr,"TimeDependent"))
            {
                for (int i = 0; i < nvariables; ++i)
                {
                    varName = m_session->GetVariable(i);
                    m_fields[i]->EvaluateBoundaryConditions(time, varName);
                }
            }
            else
            {
                string errmsg = "Unrecognised boundary condition: ";
                errmsg += userDefStr;
                ASSERTL0(false,errmsg.c_str());
            }
        }
    }

    /**
     * @brief Wall boundary conditions for compressible flow problems.
     */
    void CompressibleFlowSystem::WallBC(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();

        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Adjust the physical values of the trace to take
        // user defined boundaries into account
        int e, id1, id2, nBCEdgePts, eMax;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            // Boundary condition for epsilon term.
            if (nVariables == m_spacedim+3)
            {
                NekDouble factor  = 1.0;
                NekDouble factor2 = 1.0;

                Array<OneD, NekDouble > tmp2(nBCEdgePts, 0.0);
                Vmath::Smul(nBCEdgePts,
                            factor,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1);

                Vmath::Vsub(nBCEdgePts,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1,
                            &Fwd[nVariables-1][id2], 1);

                Vmath::Smul(nBCEdgePts,
                            factor2,
                            &Fwd[nVariables-1][id2], 1,
                            &Fwd[nVariables-1][id2], 1);

            }

            // For 2D/3D, define: v* = v - 2(v.n)n
            Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);

            // Calculate (v.n)
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &Fwd[1+i][id2], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp[0], 1,
                             &tmp[0], 1);
            }

            // Calculate 2.0(v.n)
            Vmath::Smul(nBCEdgePts, -2.0, &tmp[0], 1, &tmp[0], 1);

            // Calculate v* = v - 2.0(v.n)n
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &tmp[0], 1,
                             &m_traceNormals[i][id2], 1,
                             &Fwd[1+i][id2], 1,
                             &Fwd[1+i][id2], 1);
            }

            // Copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1], 1);
            }
        }
    }

    /**
     * @brief Wall boundary conditions for viscous compressible flow problems.
     */
    void CompressibleFlowSystem::WallViscousBC(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();

        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Take into account that for PDE based shock capturing, eps = 0 at the
        // wall. Adjust the physical values of the trace to take user defined
        // boundaries into account
        int e, id1, id2, nBCEdgePts, eMax;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            if (nVariables == m_spacedim+3)
            {
                NekDouble factor  = 0.0;
                NekDouble factor2 = 1.0;

                Array<OneD, NekDouble > tmp2(nBCEdgePts, 0.0);
                Vmath::Smul(nBCEdgePts,
                            factor,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1);

                Vmath::Vsub(nBCEdgePts,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1,
                            &Fwd[nVariables-1][id2], 1);

                Vmath::Smul(nBCEdgePts,
                            factor2,
                            &Fwd[nVariables-1][id2], 1,
                            &Fwd[nVariables-1][id2], 1);
            }

            for (i = 0; i < m_spacedim; i++)
            {
                Vmath::Neg(nBCEdgePts, &Fwd[i+1][id2], 1);
            }

            // Copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                               UpdatePhys())[id1], 1);
            }
        }
    }

    /**
     * @brief Symmetry boundary conditions for compressible flow problems.
     */
    void CompressibleFlowSystem::SymmetryBC(
        int                                      bcRegion,
        int                                      cnt,
        Array<OneD, Array<OneD, NekDouble> >    &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();

        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();

        // Get physical values of the forward trace (from exp to phys)
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Take into account that for PDE based shock capturing, eps = 0 at the
        // wall.
        int e, id1, id2, nBCEdgePts, eMax;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            if (nVariables == m_spacedim+3)
            {
                NekDouble factor  = 0.0;
                NekDouble factor2 = 1.0;

                Array<OneD, NekDouble > tmp2(nBCEdgePts, 0.0);
                Vmath::Smul(nBCEdgePts,
                            factor,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1);

                Vmath::Vsub(nBCEdgePts,
                            &Fwd[nVariables-1][id2], 1,
                            &tmp2[0], 1,
                            &Fwd[nVariables-1][id2], 1);

                Vmath::Smul(nBCEdgePts,
                            factor2,
                            &Fwd[nVariables-1][id2], 1,
                            &Fwd[nVariables-1][id2], 1);
            }

            // For 2D/3D, define: v* = v - 2(v.n)n
            Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);

            // Calculate (v.n)
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &Fwd[1+i][id2], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp[0], 1,
                             &tmp[0], 1);
            }

            // Calculate 2.0(v.n)
            Vmath::Smul(nBCEdgePts, -2.0, &tmp[0], 1, &tmp[0], 1);

            // Calculate v* = v - 2.0(v.n)n
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nBCEdgePts,
                             &tmp[0], 1,
                             &m_traceNormals[i][id2], 1,
                             &Fwd[1+i][id2], 1,
                             &Fwd[1+i][id2], 1);
            }

            // Copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nVariables; ++i)
            {
                Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                               UpdatePhys())[id1], 1);
            }
        }
    }

    /**
     * @brief Outflow characteristic boundary conditions for compressible
     * flow problems.
     */
    void CompressibleFlowSystem::RiemannInvariantBC(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nDimensions = m_spacedim;

        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();

        NekDouble gamma            = m_gamma;
        NekDouble gammaInv         = 1.0 / gamma;
        NekDouble gammaMinusOne    = gamma - 1.0;
        NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;

        Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2 (nTracePts, 0.0);
        Array<OneD, NekDouble> VnInf(nTracePts, 0.0);
        Array<OneD, NekDouble> velInf(nDimensions, 0.0);

        // Computing the normal velocity for characteristics coming
        // from outside the computational domain
        velInf[0] = m_uInf;
        Vmath::Smul(nTracePts, m_uInf, m_traceNormals[0], 1, VnInf, 1);
        if (nDimensions == 2 || nDimensions == 3)
        {
            velInf[1] = m_vInf;
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[1], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nDimensions == 3)
        {
            velInf[2] = m_wInf;
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[2], 1, tmp2, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp2, 1, VnInf, 1);
        }

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Computing the normal velocity for characteristics coming
        // from inside the computational domain
        Array<OneD, NekDouble > Vn (nTracePts, 0.0);
        Array<OneD, NekDouble > Vel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
        }

        // Computing the absolute value of the velocity in order to compute the
        // Mach number to decide whether supersonic or subsonic
        Array<OneD, NekDouble > absVel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp1, 1);
            Vmath::Vmul(nTracePts, tmp1, 1, tmp1, 1, tmp1, 1);
            Vmath::Vadd(nTracePts, tmp1, 1, absVel, 1, absVel, 1);
        }
        Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);

        // Get speed of sound
        Array<OneD, NekDouble > soundSpeed(nTracePts);
        Array<OneD, NekDouble > pressure  (nTracePts);

        for (i = 0; i < nTracePts; i++)
        {
            if (m_spacedim == 1)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[2][i] -
                               0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i]));
            }
            else if (m_spacedim == 2)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[3][i] -
                               0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                               Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
            }
            else
            {
                pressure[i] = (gammaMinusOne) * (Fwd[4][i] -
                               0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                               Fwd[2][i] * Fwd[2][i] / Fwd[0][i] +
                               Fwd[3][i] * Fwd[3][i] / Fwd[0][i]));
            }

            soundSpeed[i] = sqrt(gamma * pressure[i] / Fwd[0][i]);
        }

        // Get Mach
        Array<OneD, NekDouble > Mach(nTracePts, 0.0);
        Vmath::Vdiv(nTracePts, Vn, 1, soundSpeed, 1, Mach, 1);
        Vmath::Vabs(nTracePts, Mach, 1, Mach, 1);

        // Auxiliary variables
        int eMax;
        int e, id1, id2, nBCEdgePts, pnt;
        NekDouble cPlus, rPlus, cMinus, rMinus, VDBC, VNBC;
        Array<OneD, NekDouble> velBC(nDimensions, 0.0);
        Array<OneD, NekDouble> rhoVelBC(nDimensions, 0.0);
        NekDouble rhoBC, EBC, cBC, sBC, pBC;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        // Loop on bcRegions
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();

            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            // Loop on the points of the bcRegion
            for (i = 0; i < nBCEdgePts; i++)
            {
                pnt = id2+i;

                // Impose inflow Riemann invariant
                if (Vn[pnt] <= 0.0)
                {
                    // Subsonic flows
                    if (Mach[pnt] < 1.00)
                    {
                        // + Characteristic from inside
                        cPlus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                        rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;

                        // - Characteristic from boundary
                        cMinus = sqrt(gamma * m_pInf / m_rhoInf);
                        rMinus = VnInf[pnt] - 2.0 * cMinus * gammaMinusOneInv;
                    }
                    else
                    {
                        // + Characteristic from inside
                        cPlus = sqrt(gamma * m_pInf / m_rhoInf);
                        rPlus = VnInf[pnt] + 2.0 * cPlus * gammaMinusOneInv;

                        // + Characteristic from inside
                        cMinus = sqrt(gamma * m_pInf / m_rhoInf);
                        rMinus = VnInf[pnt] - 2.0 * cPlus * gammaMinusOneInv;
                    }

                    // Riemann boundary variables
                    VNBC = 0.5 * (rPlus + rMinus);
                    cBC = 0.25 * gammaMinusOne * (rPlus - rMinus);
                    VDBC = VNBC - VnInf[pnt];

                    // Thermodynamic boundary variables
                    sBC = m_pInf / (pow(m_rhoInf, gamma));
                    rhoBC = pow((cBC * cBC) / (gamma * sBC), gammaMinusOneInv);
                    pBC = rhoBC * cBC * cBC * gammaInv;

                    // Kinetic energy initialiasation
                    NekDouble EkBC = 0.0;

                    // Boundary velocities
                    for ( j = 0; j < nDimensions; ++j)
                    {
                        velBC[j] = velInf[j] + VDBC * m_traceNormals[j][pnt];
                        rhoVelBC[j] = rhoBC * velBC[j];
                        EkBC += 0.5 * rhoBC * velBC[j]*velBC[j];
                    }

                    // Boundary energy
                    EBC = pBC * gammaMinusOneInv + EkBC;

                    // Imposing Riemann Invariant boundary conditions
                    (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = rhoBC;
                    for (j = 0; j < nDimensions; ++j)
                    {
                        (m_fields[j+1]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = rhoVelBC[j];
                    }
                    (m_fields[nDimensions+1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = EBC;

                }
                else // Impose outflow Riemann invariant
                {
                    // Subsonic flows
                    if (Mach[pnt] < 1.00)
                    {
                        // + Characteristic from inside
                        cPlus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                        rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;

                        // - Characteristic from boundary
                        cMinus = sqrt(gamma * m_pInf / m_rhoInf);
                        rMinus = VnInf[pnt] - 2.0 * cMinus * gammaMinusOneInv;
                    }
                    else
                    {
                        // + Characteristic from inside
                        cPlus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                        rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;

                        // + Characteristic from inside
                        cMinus = sqrt(gamma * pressure[pnt] / Fwd[0][pnt]);
                        rMinus = Vn[pnt] - 2.0 * cPlus * gammaMinusOneInv;
                    }

                    // Riemann boundary variables
                    VNBC = 0.5 * (rPlus + rMinus);
                    cBC = 0.25 * gammaMinusOne * (rPlus - rMinus);
                    VDBC = VNBC - Vn[pnt];

                    // Thermodynamic boundary variables
                    sBC = pressure[pnt] / (pow(Fwd[0][pnt], gamma));
                    rhoBC = pow((cBC * cBC) / (gamma * sBC), gammaMinusOneInv);
                    pBC = rhoBC * cBC * cBC * gammaInv;

                    // Kinetic energy initialiasation
                    NekDouble EkBC = 0.0;

                    // Boundary velocities
                    for ( j = 0; j < nDimensions; ++j)
                    {
                        velBC[j] = Fwd[j+1][pnt] / Fwd[0][pnt] +
                                    VDBC * m_traceNormals[j][pnt];
                        rhoVelBC[j] = rhoBC * velBC[j];
                        EkBC += 0.5 * rhoBC * velBC[j]*velBC[j];
                    }

                    // Boundary energy
                    EBC = pBC * gammaMinusOneInv + EkBC;

                    // Imposing Riemann Invariant boundary conditions
                    (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = rhoBC;
                    for (j = 0; j < nDimensions; ++j)
                    {
                        (m_fields[j+1]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = rhoVelBC[j];
                    }
                    (m_fields[nDimensions+1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = EBC;
                }
            }
        }
    }

    /**
     * @brief Pressure outflow non-reflective
     * boundary conditions for compressible flow
     * problems.
     */
    void CompressibleFlowSystem::PressureOutflowNonReflectiveBC(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nDimensions = m_spacedim;

        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();

        NekDouble gamma            = m_gamma;
        NekDouble gammaMinusOne    = gamma - 1.0;
        NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;

        Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2 (nTracePts, 0.0);
        Array<OneD, NekDouble> VnInf(nTracePts, 0.0);
        Array<OneD, NekDouble> velInf(nDimensions, 0.0);

        // Computing the normal velocity for characteristics coming
        // from outside the computational domain
        velInf[0] = m_uInf;
        Vmath::Smul(nTracePts, m_uInf, m_traceNormals[0], 1, VnInf, 1);
        if (nDimensions == 2 || nDimensions == 3)
        {
            velInf[1] = m_vInf;
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[1], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nDimensions == 3)
        {
            velInf[2] = m_wInf;
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[2], 1, tmp2, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp2, 1, VnInf, 1);
        }

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Computing the normal velocity for characteristics coming
        // from inside the computational domain
        Array<OneD, NekDouble > Vn (nTracePts, 0.0);
        Array<OneD, NekDouble > Vel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
        }

        // Computing the absolute value of the velocity in order to compute the
        // Mach number to decide whether supersonic or subsonic
        Array<OneD, NekDouble > absVel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp1, 1);
            Vmath::Vmul(nTracePts, tmp1, 1, tmp1, 1, tmp1, 1);
            Vmath::Vadd(nTracePts, tmp1, 1, absVel, 1, absVel, 1);
        }
        Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);

        // Get speed of sound
        Array<OneD, NekDouble > soundSpeed(nTracePts);
        Array<OneD, NekDouble > pressure  (nTracePts);

        for (i = 0; i < nTracePts; i++)
        {
            if (m_spacedim == 1)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[2][i] -
                                    0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i]));
            }
            else if (m_spacedim == 2)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[3][i] -
                                    0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                           Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
            }
            else
            {
                pressure[i] = (gammaMinusOne) * (Fwd[4][i] -
                                    0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                           Fwd[2][i] * Fwd[2][i] / Fwd[0][i] +
                                           Fwd[3][i] * Fwd[3][i] / Fwd[0][i]));
            }

            soundSpeed[i] = sqrt(gamma * pressure[i] / Fwd[0][i]);
        }

        // Get Mach
        Array<OneD, NekDouble > Mach(nTracePts, 0.0);
        Vmath::Vdiv(nTracePts, Vn, 1, soundSpeed, 1, Mach, 1);
        Vmath::Vabs(nTracePts, Mach, 1, Mach, 1);

        // Auxiliary variables
        int e, id1, id2, npts, pnt;
        NekDouble rhoeb;

        // Loop on the bcRegions
        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->
             GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
               GetExp(e)->GetTotPoints();
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e) ;
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            // Loop on points of bcRegion 'e'
            for (i = 0; i < npts; i++)
            {
                pnt = id2+i;

                // Subsonic flows
                if (Mach[pnt] < 0.99)
                {
                    // Kinetic energy calculation
                    NekDouble Ek = 0.0;
                    for (j = 1; j < nVariables-1; ++j)
                    {
                        Ek += 0.5 * (Fwd[j][pnt] * Fwd[j][pnt]) / Fwd[0][pnt];
                    }

                    rhoeb = m_pInf * gammaMinusOneInv + Ek;

                    // Partial extrapolation for subsonic cases
                    for (j = 0; j < nVariables-1; ++j)
                    {
                        (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = Fwd[j][pnt];
                    }

                    (m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = 2.0 * rhoeb - Fwd[nVariables-1][pnt];
                }
                // Supersonic flows
                else
                {
                    for (j = 0; j < nVariables; ++j)
                    {
                        // Extrapolation for supersonic cases
                        (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = Fwd[j][pnt];
                    }
                }
            }
        }
    }

    /**
     * @brief Pressure outflow boundary conditions for compressible flow
     * problems.
     */
    void CompressibleFlowSystem::PressureOutflowBC(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nDimensions = m_spacedim;

        const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

        NekDouble gamma            = m_gamma;
        NekDouble gammaMinusOne    = gamma - 1.0;
        NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;

        Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2 (nTracePts, 0.0);
        Array<OneD, NekDouble> VnInf(nTracePts, 0.0);
        Array<OneD, NekDouble> velInf(nDimensions, 0.0);

        // Computing the normal velocity for characteristics coming
        // from outside the computational domain
        velInf[0] = m_uInf;
        Vmath::Smul(nTracePts, m_uInf, m_traceNormals[0], 1, VnInf, 1);
        if (nDimensions == 2 || nDimensions == 3)
        {
            velInf[1] = m_vInf;
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[1], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nDimensions == 3)
        {
            velInf[2] = m_wInf;
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[2], 1, tmp2, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp2, 1, VnInf, 1);
        }

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Computing the normal velocity for characteristics coming
        // from inside the computational domain
        Array<OneD, NekDouble > Vn (nTracePts, 0.0);
        Array<OneD, NekDouble > Vel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
        }

        // Computing the absolute value of the velocity in order to compute the
        // Mach number to decide whether supersonic or subsonic
        Array<OneD, NekDouble > absVel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp1, 1);
            Vmath::Vmul(nTracePts, tmp1, 1, tmp1, 1, tmp1, 1);
            Vmath::Vadd(nTracePts, tmp1, 1, absVel, 1, absVel, 1);
        }
        Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);

        // Get speed of sound
        Array<OneD, NekDouble > soundSpeed(nTracePts);
        Array<OneD, NekDouble > pressure  (nTracePts);

        for (i = 0; i < nTracePts; i++)
        {
            if (m_spacedim == 1)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[2][i] -
                                0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i]));
            }
            else if (m_spacedim == 2)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[3][i] -
                                0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                       Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
            }
            else
            {
                pressure[i] = (gammaMinusOne) * (Fwd[4][i] -
                                0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                       Fwd[2][i] * Fwd[2][i] / Fwd[0][i] +
                                       Fwd[3][i] * Fwd[3][i] / Fwd[0][i]));
            }

            soundSpeed[i] = sqrt(gamma * pressure[i] / Fwd[0][i]);
        }

        // Get Mach
        Array<OneD, NekDouble > Mach(nTracePts, 0.0);
        Vmath::Vdiv(nTracePts, Vn, 1, soundSpeed, 1, Mach, 1);
        Vmath::Vabs(nTracePts, Mach, 1, Mach, 1);

        // Auxiliary variables
        int e, id1, id2, npts, pnt;
        NekDouble rhoeb;

        // Loop on the bcRegions
        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->
             GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetTotPoints();
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e) ;
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            // Loop on points of bcRegion 'e'
            for (i = 0; i < npts; i++)
            {
                pnt = id2+i;

                // Subsonic flows
                if (Mach[pnt] < 0.99)
                {
                    // Kinetic energy calculation
                    NekDouble Ek = 0.0;
                    for (j = 1; j < nVariables-1; ++j)
                    {
                        Ek += 0.5 * (Fwd[j][pnt] * Fwd[j][pnt]) / Fwd[0][pnt];
                    }

                    rhoeb = m_pInf * gammaMinusOneInv + Ek;

                    // Partial extrapolation for subsonic cases
                    for (j = 0; j < nVariables-1; ++j)
                    {
                        (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = Fwd[j][pnt];
                    }

                    (m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = rhoeb;
                }
                // Supersonic flows
                else
                {
                    for (j = 0; j < nVariables; ++j)
                    {
                        // Extrapolation for supersonic cases
                        (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = Fwd[j][pnt];
                    }
                }
            }
        }
    }


    /**
     * @brief Pressure outflow boundary conditions for compressible flow
     * problems.
     */
    void CompressibleFlowSystem::PressureOutflowFileBC(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nDimensions = m_spacedim;

        const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

        NekDouble gamma            = m_gamma;
        NekDouble gammaMinusOne    = gamma - 1.0;
        NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;

        Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2 (nTracePts, 0.0);
        Array<OneD, NekDouble> VnInf(nTracePts, 0.0);
        Array<OneD, NekDouble> velInf(nDimensions, 0.0);

        // Computing the normal velocity for characteristics coming
        // from outside the computational domain
        velInf[0] = m_uInf;
        Vmath::Smul(nTracePts, m_uInf, m_traceNormals[0], 1, VnInf, 1);
        if (nDimensions == 2 || nDimensions == 3)
        {
            velInf[1] = m_vInf;
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[1], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nDimensions == 3)
        {
            velInf[2] = m_wInf;
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[2], 1, tmp2, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp2, 1, VnInf, 1);
        }

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);

        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Computing the normal velocity for characteristics coming
        // from inside the computational domain
        Array<OneD, NekDouble > Vn (nTracePts, 0.0);
        Array<OneD, NekDouble > Vel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
        }

        // Computing the absolute value of the velocity in order to compute the
        // Mach number to decide whether supersonic or subsonic
        Array<OneD, NekDouble > absVel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp1, 1);
            Vmath::Vmul(nTracePts, tmp1, 1, tmp1, 1, tmp1, 1);
            Vmath::Vadd(nTracePts, tmp1, 1, absVel, 1, absVel, 1);
        }
        Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);

        // Get speed of sound
        Array<OneD, NekDouble > soundSpeed (nTracePts, 0.0);
        Array<OneD, NekDouble > pressure   (nTracePts, 0.0);

        for (i = 0; i < nTracePts; i++)
        {
            if (m_spacedim == 1)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[2][i] -
                    0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i]));
            }
            else if (m_spacedim == 2)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[3][i] -
                    0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                           Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
            }
            else
            {
                pressure[i] = (gammaMinusOne) * (Fwd[4][i] -
                    0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                           Fwd[2][i] * Fwd[2][i] / Fwd[0][i] +
                           Fwd[3][i] * Fwd[3][i] / Fwd[0][i]));
            }

            soundSpeed[i] = sqrt(gamma * pressure[i] / Fwd[0][i]);
        }

        // Get Mach
        Array<OneD, NekDouble > Mach(nTracePts, 0.0);
        Vmath::Vdiv(nTracePts, Vn, 1, soundSpeed, 1, Mach, 1);
        Vmath::Vabs(nTracePts, Mach, 1, Mach, 1);

        // Auxiliary variables
        int e, id1, id2, npts, pnt;
        NekDouble rhoeb;

        // Loop on the bcRegions
        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->
             GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                                                GetExp(e)->GetTotPoints();
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                                                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            // Loop on points of bcRegion 'e'
            for (i = 0; i < npts; i++)
            {
                pnt = id2+i;

                // Subsonic flows
                if (Mach[pnt] < 0.99)
                {
                    // Kinetic energy calculation
                    NekDouble Ek = 0.0;
                    for (j = 1; j < nVariables-1; ++j)
                    {
                        Ek += 0.5 * (Fwd[j][pnt] * Fwd[j][pnt]) / Fwd[0][pnt];
                    }

                    rhoeb = m_pressureStorage[id1+i] * gammaMinusOneInv + Ek;

                    // Partial extrapolation for subsonic cases
                    for (j = 0; j < nVariables-1; ++j)
                    {
                        (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = Fwd[j][pnt];
                    }

                    (m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = rhoeb;
                }
                // Supersonic flows
                else
                {
                    for (j = 0; j < nVariables; ++j)
                    {
                        // Extrapolation for supersonic cases
                        (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = Fwd[j][pnt];
                    }
                }
            }
        }
    }


    /**
     * @brief Pressure inflow boundary conditions for compressible flow
     * problems where either the density and the velocities are assigned from a
     * file or the full state is assigned from a file (depending on the problem
     * type, either subsonic or supersonic).
     */
    void CompressibleFlowSystem::PressureInflowFileBC(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nDimensions = m_spacedim;

        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();

        NekDouble gamma            = m_gamma;
        NekDouble gammaMinusOne    = gamma - 1.0;
        NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;

        Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2 (nTracePts, 0.0);
        Array<OneD, NekDouble> VnInf(nTracePts, 0.0);
        Array<OneD, NekDouble> velInf(nDimensions, 0.0);

        // Computing the normal velocity for characteristics coming
        // from outside the computational domain
        velInf[0] = m_uInf;
        Vmath::Smul(nTracePts, m_uInf, m_traceNormals[0], 1, VnInf, 1);
        if (nDimensions == 2 || nDimensions == 3)
        {
            velInf[1] = m_vInf;
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[1], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nDimensions == 3)
        {
            velInf[2] = m_wInf;
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[2], 1, tmp2, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp2, 1, VnInf, 1);
        }

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);

        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        // Computing the normal velocity for characteristics coming
        // from inside the computational domain
        Array<OneD, NekDouble > Vn (nTracePts, 0.0);
        Array<OneD, NekDouble > Vel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
        }

        // Computing the absolute value of the velocity in order to compute the
        // Mach number to decide whether supersonic or subsonic
        Array<OneD, NekDouble > absVel(nTracePts, 0.0);
        for (i = 0; i < nDimensions; ++i)
        {
            Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp1, 1);
            Vmath::Vmul(nTracePts, tmp1, 1, tmp1, 1, tmp1, 1);
            Vmath::Vadd(nTracePts, tmp1, 1, absVel, 1, absVel, 1);
        }
        Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);

        // Get speed of sound
        Array<OneD, NekDouble > soundSpeed (nTracePts, 0.0);
        Array<OneD, NekDouble > pressure   (nTracePts, 0.0);

        for (i = 0; i < nTracePts; i++)
        {
            if (m_spacedim == 1)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[2][i] -
                                0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i]));
            }
            else if (m_spacedim == 2)
            {
                pressure[i] = (gammaMinusOne) * (Fwd[3][i] -
                                0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                       Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
            }
            else
            {
                pressure[i] = (gammaMinusOne) * (Fwd[4][i] -
                                0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                       Fwd[2][i] * Fwd[2][i] / Fwd[0][i] +
                                       Fwd[3][i] * Fwd[3][i] / Fwd[0][i]));
            }

            soundSpeed[i] = sqrt(gamma * pressure[i] / Fwd[0][i]);
        }

        // Get Mach
        Array<OneD, NekDouble > Mach(nTracePts, 0.0);
        Vmath::Vdiv(nTracePts, Vn, 1, soundSpeed, 1, Mach, 1);
        Vmath::Vabs(nTracePts, Mach, 1, Mach, 1);

        // Auxiliary variables
        int e, id1, id2, npts, pnt;
        NekDouble rhoeb;

        // Loop on the bcRegions
        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->
             GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetTotPoints();
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            // Loop on points of bcRegion 'e'
            for (i = 0; i < npts; i++)
            {
                pnt = id2+i;

                // Subsonic flows
                if (Mach[pnt] < 0.99)
                {
                    // Partial extrapolation for subsonic cases
                    for (j = 0; j < nVariables-1; ++j)
                    {
                        (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = m_fieldStorage[j][id1+i];
                    }

                    // Kinetic energy calculation
                    NekDouble Ek = 0.0;
                    for (j = 1; j < nVariables-1; ++j)
                    {
                        Ek += 0.5 * (m_fieldStorage[j][id1+i] *
                                     m_fieldStorage[j][id1+i]) /
                                        m_fieldStorage[0][id1+i];
                    }
                    rhoeb = gammaMinusOneInv * pressure[pnt] + Ek;

                    (m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = rhoeb;
                }
                // Supersonic flows
                else
                {
                    for (j = 0; j < nVariables; ++j)
                    {
                        // Extrapolation for supersonic cases
                        (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                         UpdatePhys())[id1+i] = Fwd[j][pnt];
                    }
                }
            }
        }
    }


    /**
     * @brief Extrapolation of order 0 for all the variables such that,
     * at the boundaries, a trivial Riemann problem is solved.
     */
    void CompressibleFlowSystem::ExtrapOrder0BC(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i, j;
        int e, pnt;
        int id1, id2, nBCEdgePts;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        int nDimensions = m_spacedim;

        const Array<OneD, const int> &traceBndMap
            = m_fields[0]->GetTraceBndMap();

        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        int eMax;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        // Loop on bcRegions
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e) ;
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

            // Loop on points of bcRegion 'e'
            for (i = 0; i < nBCEdgePts; i++)
            {
                pnt = id2+i;

                // Setting up bcs for density
                (m_fields[0]->GetBndCondExpansions()[bcRegion]->
                    UpdatePhys())[id1+i] = Fwd[0][pnt];

                // Setting up bcs for velocities
                for (j = 1; j <=nDimensions; ++j)
                {
                    (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[j][pnt];
                }

                // Setting up bcs for energy
                (m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->
                    UpdatePhys())[id1+i] = Fwd[nVariables-1][pnt];
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physfield[i+1], 1, flux[0][i], 1);
        }

        GetVelocityVector(physfield, velocity);
        GetPressure      (physfield, velocity, pressure);

        // Flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield[i+1], 1,
                            flux[i+1][j], 1);
            }

            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux[i+1][i], 1, pressure, 1, flux[i+1][i], 1);
        }

        // Flux vector for energy.
        Vmath::Vadd(nq, physfield[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);

        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux[m_spacedim+1][j], 1);
        }

        // For the smooth viscosity model
        if (nVariables == m_spacedim+3)
        {
            // Add a zero row for the advective fluxes
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Zero(nq, flux[m_spacedim+2][j], 1);
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations
     * by using the de-aliasing technique.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        Array<OneD, Array<OneD, NekDouble> > physfield_interp(nVariables);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > flux_interp(
                                                            nVariables);

        for (i = 0; i < nVariables; ++ i)
        {
            physfield_interp[i] = Array<OneD, NekDouble>(nq);
            flux_interp[i] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], physfield_interp[i]);

            for (j = 0; j < m_spacedim; ++j)
            {
                flux_interp[i][j] = Array<OneD, NekDouble>(nq);
            }
        }

        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);

            // Galerkin project solution back to original space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale, physfield_interp[i+1], flux[0][i]);
        }

        GetVelocityVector(physfield_interp, velocity);
        GetPressure      (physfield_interp, velocity, pressure);

        // Evaluation of flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield_interp[i+1], 1,
                            flux_interp[i+1][j], 1);
            }

            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux_interp[i+1][i], 1, pressure,1,
                        flux_interp[i+1][i], 1);
        }

        // Galerkin project solution back to origianl space
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale, flux_interp[i+1][j], flux[i+1][j]);
            }
        }

        // Evaluation of flux vector for energy
        Vmath::Vadd(nq, physfield_interp[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);

        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux_interp[m_spacedim+1][j], 1);

            // Galerkin project solution back to origianl space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale,
                flux_interp[m_spacedim+1][j],
                flux[m_spacedim+1][j]);
        }
    }


    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void CompressibleFlowSystem::GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int j, k;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();

        // Stokes hypotesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > mu2                (nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, &mu[0], 1, &thermalConductivity[0], 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, &mu[0], 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        &thermalConductivity[0], 1);
        }

        // Computing diagonal terms of viscous stress tensor
        Array<OneD, Array<OneD, NekDouble> > tmp(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > Sgg(m_spacedim);

        // mu2 = 2 * mu
        Vmath::Smul(nPts, 2.0, &mu[0], 1, &mu2[0], 1);

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, &divVel[0], 1, &derivativesO1[j][j][0], 1,
                        &divVel[0], 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, &divVel[0], 1, &divVel[0], 1);
        Vmath::Vmul(nPts, &mu[0], 1, &divVel[0], 1, &divVel[0], 1);

        // Diagonal terms of viscous stress tensor (Sxx, Syy, Szz)
        // Sjj = 2 * mu * du_j/dx_j - (2 / 3) * mu * sum_j(du_j/dx_j)
        for (j = 0; j < m_spacedim; ++j)
        {
            tmp[j] = Array<OneD, NekDouble>(nPts, 0.0);
            Sgg[j] = Array<OneD, NekDouble>(nPts, 0.0);

            Vmath::Vmul(nPts, &mu2[0], 1, &derivativesO1[j][j][0], 1,
                        &tmp[j][0], 1);

            Vmath::Vadd(nPts, &tmp[j][0], 1, &divVel[0], 1, &Sgg[j][0], 1);
        }

        // Extra diagonal terms of viscous stress tensor (Sxy, Sxz, Syz)
        // Note: they exist for 2D and 3D problems only
        Array<OneD, NekDouble > Sxy(nPts, 0.0);
        Array<OneD, NekDouble > Sxz(nPts, 0.0);
        Array<OneD, NekDouble > Syz(nPts, 0.0);

        if (m_spacedim == 2)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][1][0], 1,
                        &derivativesO1[1][0][0], 1, &Sxy[0], 1);

            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][1][0], 1,
                        &derivativesO1[1][0][0], 1, &Sxy[0], 1);

            // Sxz = (du/dz + dw/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][2][0], 1,
                        &derivativesO1[2][0][0], 1, &Sxz[0], 1);

            // Syz = (dv/dz + dw/dy)
            Vmath::Vadd(nPts, &derivativesO1[1][2][0], 1,
                        &derivativesO1[2][1][0], 1, &Syz[0], 1);

            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);

            // Sxz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxz[0], 1, &Sxz[0], 1);

            // Syz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Syz[0], 1, &Syz[0], 1);
        }

        // Energy-related terms
        Array<OneD, NekDouble > STx(nPts, 0.0);
        Array<OneD, NekDouble > STy(nPts, 0.0);
        Array<OneD, NekDouble > STz(nPts, 0.0);

        // Building the viscous flux vector

        // Viscous flux vector for the rho equation
        for (k = 0; k < m_spacedim; ++k)
        {
            Vmath::Zero(nPts, viscousTensor[k][0], 1);
        }

        if (m_spacedim == 1)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);

            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);

            // k * dT/dx
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[0][1][0], 1, &tmp1[0], 1);

            // STx = u * Sxx + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
        }
        else if (m_spacedim == 2)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);

            // Computation of STx

            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);

            // v * Sxy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sxy[0], 1, &tmp1[0], 1);

            // k * dT/dx
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[0][2][0], 1, &tmp2[0], 1);

            // STx = u * Sxx + v * Sxy + K * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);

            // Computation of STy

            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);

            // v * Syy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sgg[1][0], 1, &STy[0], 1);

            // u * Sxy
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sxy[0], 1, &tmp1[0], 1);

            // k * dT/dy
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[1][2][0], 1, &tmp2[0], 1);

            // STy = v * Syy + u * Sxy + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            Array<OneD, NekDouble > tmp3(nPts, 0.0);

            // Computation of STx

            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);

            // v * Sxy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sxy[0], 1, &tmp1[0], 1);

            // v * Sxz
            Vmath::Vmul(nPts, &physfield[2][0], 1, &Sxz[0], 1, &tmp2[0], 1);

            // k * dT/dx
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[0][3][0], 1, &tmp3[0], 1);

            // STx = u * Sxx + v * Sxy + w * Sxz + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp3[0], 1, &STx[0], 1);

            // Computation of STy

            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);

            // v * Syy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sgg[1][0], 1, &STy[0], 1);

            // u * Sxy
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sxy[0], 1, &tmp1[0], 1);

            // w * Syz
            Vmath::Vmul(nPts, &physfield[2][0], 1, &Syz[0], 1, &tmp2[0], 1);

            // k * dT/dy
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[1][3][0], 1, &tmp3[0], 1);

            // STy = v * Syy + u * Sxy + w * Syz + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp3[0], 1, &STy[0], 1);

            // Computation of STz

            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);

            // w * Szz
            Vmath::Vmul(nPts, &physfield[2][0], 1, &Sgg[2][0], 1, &STz[0], 1);

            // u * Sxz
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sxz[0], 1, &tmp1[0], 1);

            // v * Syz
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Syz[0], 1, &tmp2[0], 1);

            // k * dT/dz
            Vmath::Vmul(nPts, &thermalConductivity[0], 1,
                        &derivativesO1[2][3][0], 1, &tmp3[0], 1);

            // STz = w * Szz + u * Sxz + v * Syz + K * dT/dz
            Vmath::Vadd(nPts, &STz[0], 1, &tmp1[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp2[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp3[0], 1, &STz[0], 1);
        }

        switch (m_spacedim)
        {
            case 1:
            {
                // f_11v = f_rho = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);

                // f_21v = f_rhou
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);

                // f_31v = f_E
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][2][0], 1);
                break;
            }
            case 2:
            {
                // f_11v = f_rho1 = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nPts, &viscousTensor[1][0][0], 1);

                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[1][1][0], 1);

                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor[1][2][0], 1);

                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][3][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor[1][3][0], 1);
                break;
            }
            case 3:
            {
                // f_11v = f_rho1 = 0
                Vmath::Zero(nPts, &viscousTensor[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nPts, &viscousTensor[1][0][0], 1);
                // f_13v = f_rho3 = 0
                Vmath::Zero(nPts, &viscousTensor[2][0][0], 1);

                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[1][1][0], 1);
                // f_23v = f_rhou3
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor[2][1][0], 1);

                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor[1][2][0], 1);
                // f_33v = f_rhov3
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor[2][2][0], 1);

                // f_31v = f_rhow1
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor[0][3][0], 1);
                // f_32v = f_rhow2
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor[1][3][0], 1);
                // f_33v = f_rhow3
                Vmath::Vcopy(nPts, &Sgg[2][0], 1, &viscousTensor[2][3][0], 1);

                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][4][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor[1][4][0], 1);
                // f_43v = f_E3
                Vmath::Vcopy(nPts, &STz[0], 1, &viscousTensor[2][4][0], 1);
                break;
            }
            default:
            {
                ASSERTL0(false, "Illegal expansion dimension");
            }
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void CompressibleFlowSystem::GetViscousFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
#if 0
        int i, j, k;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();

        int variables_phys = physfield.num_elements();

        // Factor to rescale 1d points in dealiasing.
        NekDouble OneDptscale = 2;

        // Get number of points to dealias a cubic non-linearity
        nPts = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        int nVariables_aux = derivativesO1[0].num_elements();

        Array<OneD, Array<OneD, NekDouble> > physfield_interp(variables_phys);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > derivativesO1_interp(
                                                                    m_spacedim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >viscousTensor_interp(
                                                                    m_spacedim);

        for (i = 0; i < m_spacedim; ++ i)
        {
            viscousTensor_interp[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                    nVariables);
            for (j = 0; j < nVariables; ++j)
            {
                viscousTensor_interp[i][j] = Array<OneD, NekDouble>(nPts);
            }
        }

        // Stokes hypotesis
        NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu         (nPts, 0.0);
        Array<OneD, NekDouble > mu2        (nPts, 0.0);
        Array<OneD, NekDouble > divVel     (nPts, 0.0);
        Array<OneD, NekDouble > pressure   (nPts, 0.0);
        Array<OneD, NekDouble > temperature(nPts, 0.0);

        for (i = 0; i < nVariables; ++i)
        {
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], fields_interp[i]);
        }

        for (i = 0; i < variables_phys; ++i)
        {
            physfield_interp[i] = Array<OneD, NekDouble> (nPts);

            // Interpolation to higher space
            m_fields[0]->PhysInterp1DScaled(OneDptscale, physfield[i],
                                            physfield_interp[i]);
        }

        for (i = 0; i < m_spacedim; ++i)
        {
            derivativesO1_interp[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                            nVariables_aux);
            for (j = 0; j < nVariables_aux; ++j)
            {
                derivativesO1_interp[i][j] = Array<OneD, NekDouble>(nPts);
                m_fields[0]->PhysInterp1DScaled(
                    OneDptscale, derivativesO1[i][j],
                    derivativesO1_interp[i][j]);
            }
        }

        // Thermodynamic related quantities
        GetPressure(fields_interp, pressure);
        GetTemperature(fields_interp, pressure, temperature);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            GetDynamicViscosity(fields_interp[variables_phys-1], mu);
        }
        else
        {
            Vmath::Sadd(nPts, m_mu, &mu[0], 1, &mu[0], 1);
        }

        // Computing diagonal terms of viscous stress tensor
        Array<OneD, Array<OneD, NekDouble> > tmp(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > Sgg(m_spacedim);

        // mu2 = 2 * mu
        Vmath::Smul(nPts, 2.0, &mu[0], 1, &mu2[0], 1);

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, &divVel[0], 1, &derivativesO1_interp[j][j][0], 1,
                        &divVel[0], 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, &divVel[0], 1, &divVel[0], 1);
        Vmath::Vmul(nPts, &mu[0], 1, &divVel[0], 1, &divVel[0], 1);

        // Digonal terms of viscous stress tensor (Sxx, Syy, Szz)
        // Sjj = 2 * mu * du_j/dx_j - (2 / 3) * mu * sum_j(du_j/dx_j)
        for (j = 0; j < m_spacedim; ++j)
        {
            tmp[j] = Array<OneD, NekDouble>(nPts, 0.0);
            Sgg[j] = Array<OneD, NekDouble>(nPts, 0.0);

            Vmath::Vmul(nPts, &mu2[0], 1, &derivativesO1_interp[j][j][0], 1,
                        &tmp[j][0], 1);

            Vmath::Vadd(nPts, &tmp[j][0], 1, &divVel[0], 1, &Sgg[j][0], 1);
        }

        // Extra diagonal terms of viscous stress tensor (Sxy, Sxz, Syz)
        // Note: they exist for 2D and 3D problems only
        Array<OneD, NekDouble > Sxy(nPts, 0.0);
        Array<OneD, NekDouble > Sxz(nPts, 0.0);
        Array<OneD, NekDouble > Syz(nPts, 0.0);

        if (m_spacedim == 2)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1_interp[0][1][0], 1,
                        &derivativesO1_interp[1][0][0], 1, &Sxy[0], 1);

            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1_interp[0][1][0], 1,
                        &derivativesO1_interp[1][0][0], 1, &Sxy[0], 1);

            // Sxz = (du/dz + dw/dx)
            Vmath::Vadd(nPts, &derivativesO1_interp[0][2][0], 1,
                        &derivativesO1_interp[2][0][0], 1, &Sxz[0], 1);

            // Syz = (dv/dz + dw/dy)
            Vmath::Vadd(nPts, &derivativesO1_interp[1][2][0], 1,
                        &derivativesO1_interp[2][1][0], 1, &Syz[0], 1);

            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);

            // Sxz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxz[0], 1, &Sxz[0], 1);

            // Syz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Syz[0], 1, &Syz[0], 1);
        }

        // Energy-related terms
        Array<OneD, NekDouble > STx(nPts, 0.0);
        Array<OneD, NekDouble > STy(nPts, 0.0);
        Array<OneD, NekDouble > STz(nPts, 0.0);
        // Building the viscous flux vector
        if (i == 0)
        {
            // Viscous flux vector for the rho equation
            for (k = 0; k < m_spacedim; ++k)
            {
                Vmath::Zero(nPts, viscousTensor_interp[k][i], 1);
            }
        }

        if (m_spacedim == 1)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);

            // u * Sxx
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);

            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[0][1][0], 1, &tmp1[0], 1);

            // STx = u * Sxx + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
        }
        else if (m_spacedim == 2)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);

            // Computation of STx

            // u * Sxx
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);

            // v * Sxy
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);

            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[0][2][0], 1, &tmp2[0], 1);

            // STx = u * Sxx + v * Sxy + K * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);

            // Computation of STy

            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);

            // v * Syy
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Sgg[1][0], 1, &STy[0], 1);

            // u * Sxy
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);

            // k * dT/dy
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[1][2][0], 1, &tmp2[0], 1);

            // STy = v * Syy + u * Sxy + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
        }
        else if (m_spacedim == 3)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            Array<OneD, NekDouble > tmp3(nPts, 0.0);

            // Computation of STx

            // u * Sxx
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sgg[0][0], 1, &STx[0], 1);

            // v * Sxy
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);

            // v * Sxy
            Vmath::Vmul(nPts, &physfield_interp[2][0], 1,
                        &Sxz[0], 1, &tmp2[0], 1);

            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[0][3][0], 1, &tmp3[0], 1);

            // STx = u * Sxx + v * Sxy + w * Sxz + (K / mu) * dT/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp3[0], 1, &STx[0], 1);

            // Computation of STy

            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);

            // v * Syy
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Sgg[1][0], 1, &STy[0], 1);

            // u * Sxy
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sxy[0], 1, &tmp1[0], 1);

            // w * Syz
            Vmath::Vmul(nPts, &physfield_interp[2][0], 1,
                        &Syz[0], 1, &tmp2[0], 1);

            // k * dT/dy
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[1][3][0], 1, &tmp3[0], 1);

            // STy = v * Syy + u * Sxy + w * Syz + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp3[0], 1, &STy[0], 1);

            // Computation of STz

            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);

            // w * Szz
            Vmath::Vmul(nPts, &physfield_interp[2][0], 1,
                        &Sgg[2][0], 1, &STz[0], 1);

            // u * Sxz
            Vmath::Vmul(nPts, &physfield_interp[0][0], 1,
                        &Sxz[0], 1, &tmp1[0], 1);

            // v * Syz
            Vmath::Vmul(nPts, &physfield_interp[1][0], 1,
                        &Syz[0], 1, &tmp2[0], 1);

            // k * dT/dz
            Vmath::Smul(nPts, m_thermalConductivity,
                        &derivativesO1_interp[2][3][0], 1, &tmp3[0], 1);

            // STz = w * Szz + u * Sxz + v * Syz + K * dT/dz
            Vmath::Vadd(nPts, &STz[0], 1, &tmp1[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp2[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp3[0], 1, &STz[0], 1);
        }

        switch (m_spacedim)
        {
            case 1:
            {

                int nq = physfield[0].num_elements();
                // f_11v = f_rho = 0
                Vmath::Zero(nq, &viscousTensor_interp[0][0][0], 1);

                // f_21v = f_rhou
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor_interp[0][1][0], 1);

                // f_31v = f_E
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor_interp[0][2][0], 1);
                break;
            }
            case 2:
            {
                int nq = physfield[0].num_elements();
                // f_11v = f_rho1 = 0
                Vmath::Zero(nq, &viscousTensor_interp[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nq, &viscousTensor_interp[1][0][0], 1);

                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor_interp[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor_interp[1][1][0], 1);

                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor_interp[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor_interp[1][2][0], 1);

                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor_interp[0][3][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor_interp[1][3][0], 1);
                break;
            }
            case 3:
            {
                int nq = physfield[0].num_elements();
                // f_11v = f_rho1 = 0
                Vmath::Zero(nq, &viscousTensor_interp[0][0][0], 1);
                // f_12v = f_rho2 = 0
                Vmath::Zero(nq, &viscousTensor_interp[1][0][0], 1);
                // f_13v = f_rho3 = 0
                Vmath::Zero(nq, &viscousTensor_interp[2][0][0], 1);

                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor_interp[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor_interp[1][1][0], 1);
                // f_23v = f_rhou3
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor_interp[2][1][0], 1);

                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor_interp[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor_interp[1][2][0], 1);
                // f_33v = f_rhov3
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor_interp[2][2][0], 1);

                // f_31v = f_rhow1
                Vmath::Vcopy(nPts, &Sxz[0],    1, &viscousTensor_interp[0][3][0], 1);
                // f_32v = f_rhow2
                Vmath::Vcopy(nPts, &Syz[0],    1, &viscousTensor_interp[1][3][0], 1);
                // f_33v = f_rhow3
                Vmath::Vcopy(nPts, &Sgg[2][0], 1, &viscousTensor_interp[2][3][0], 1);

                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor_interp[0][4][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor_interp[1][4][0], 1);
                // f_43v = f_E3
                Vmath::Vcopy(nPts, &STz[0], 1, &viscousTensor_interp[2][4][0], 1);
                break;
            }
            default:
            {
                ASSERTL0(false, "Illegal expansion dimension");
            }
        }

        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 1; j < nVariables; ++j)
            {
                // Galerkin project solution back to origianl space
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale,
                    viscousTensor_interp[i][j],
                    viscousTensor[i][j]);
            }
        }
#endif
}

    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal
     * gas law.
     *
     * @param physfield  Input momentum.
     * @param pressure   Computed pressure field.
     */
    void CompressibleFlowSystem::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD,                   NekDouble>   &pressure)
    {
        int       nBCEdgePts  = physfield[0].num_elements();
        NekDouble alpha = -0.5;

        // Calculate ||rho v||^2
        Vmath::Vmul(nBCEdgePts, physfield[1], 1, physfield[1], 1, pressure, 1);
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, physfield[1+i], 1, physfield[1+i], 1,
                               pressure,       1, pressure,       1);
        }
        // Divide by rho to get rho*||v||^2
        Vmath::Vdiv (nBCEdgePts, pressure, 1, physfield[0], 1, pressure, 1);
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(nBCEdgePts, alpha,
                     pressure, 1, physfield[m_spacedim+1], 1, pressure, 1);
        // Multiply by (gamma-1)
        Vmath::Smul (nBCEdgePts, m_gamma-1, pressure, 1, pressure, 1);
    }

    /**
     * @brief Compute the enthalpy term \f$ H = E + p/rho \$.
     */
    void CompressibleFlowSystem::GetEnthalpy(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &pressure,
                  Array<OneD,                   NekDouble>   &enthalpy)
    {
        int npts  = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> tmp(npts, 0.0);

        // Calculate E = rhoE/rho
        Vmath::Vdiv(npts, physfield[m_spacedim+1], 1, physfield[0], 1, tmp, 1);
        // Calculate p/rho
        Vmath::Vdiv(npts, pressure, 1, physfield[0], 1, enthalpy, 1);
        // Calculate H = E + p/rho
        Vmath::Vadd(npts, tmp, 1, enthalpy, 1, enthalpy, 1);
    }

    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal
     * gas law.
     *
     * This is a slightly optimised way to calculate the pressure field which
     * avoids division by the density field if the velocity field has already
     * been calculated.
     *
     * @param physfield  Input momentum.
     * @param velocity   Velocity vector.
     * @param pressure   Computed pressure field.
     */
    void CompressibleFlowSystem::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        const Array<OneD, const Array<OneD, NekDouble> > &velocity,
              Array<OneD,                   NekDouble>   &pressure)
    {
        int nBCEdgePts = physfield[0].num_elements();
        NekDouble alpha = -0.5;

        // Calculate ||\rho v||^2.
        Vmath::Vmul (nBCEdgePts, velocity[0], 1, physfield[1], 1, pressure, 1);
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, velocity[i], 1, physfield[1+i], 1,
                               pressure,    1, pressure,       1);
        }
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(nBCEdgePts,     alpha,
                     pressure, 1, physfield[m_spacedim+1], 1, pressure, 1);
        // Multiply by (gamma-1).
        Vmath::Smul (nBCEdgePts, m_gamma-1, pressure, 1, pressure, 1);
    }

    /**
     * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
     * \f$ \rho\mathbf{v} \f$.
     *
     * @param physfield  Momentum field.
     * @param velocity   Velocity field.
     */
    void CompressibleFlowSystem::GetVelocityVector(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &velocity)
    {
        const int nBCEdgePts = physfield[0].num_elements();

        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vdiv(nBCEdgePts, physfield[1+i], 1, physfield[0], 1,
                              velocity[i],    1);
        }
    }

    /**
     * @brief Compute the temperature \f$ T = p/\rho R \f$.
     *
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param temperature  The resulting temperature \f$ T \f$.
     */
    void CompressibleFlowSystem::GetTemperature(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        Array<OneD,                         NekDouble  > &pressure,
        Array<OneD,                         NekDouble  > &temperature)
    {
        const int nq = physfield[0].num_elements();

        Vmath::Vdiv(nq, pressure, 1, physfield[0], 1, temperature, 1);
        Vmath::Smul(nq, 1.0/m_gasConstant, temperature, 1, temperature, 1);
    }

    /**
     * @brief Compute the sound speed \f$ c = sqrt(\gamma p/\rho) \f$.
     *
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param soundspeed   The resulting sound speed \f$ c \f$.
     */
    void CompressibleFlowSystem::GetSoundSpeed(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD,             NekDouble  > &pressure,
              Array<OneD,             NekDouble  > &soundspeed)
    {
        const int nq = m_fields[0]->GetTotPoints();
        Vmath::Vdiv (nq, pressure, 1, physfield[0], 1, soundspeed, 1);
        Vmath::Smul (nq, m_gamma, soundspeed, 1, soundspeed, 1);
        Vmath::Vsqrt(nq, soundspeed, 1, soundspeed, 1);
    }

    /**
     * @brief Compute the mach number \f$ M = \| \mathbf{v} \|^2 / c \f$.
     *
     * @param physfield    Input physical field.
     * @param soundfield   The speed of sound corresponding to physfield.
     * @param mach         The resulting mach number \f$ M \f$.
     */
    void CompressibleFlowSystem::GetMach(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD,             NekDouble  > &soundspeed,
        Array<OneD,             NekDouble  > &mach)
    {
        const int nq = m_fields[0]->GetTotPoints();

        Vmath::Vmul(nq, physfield[1], 1, physfield[1], 1, mach, 1);

        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nq, physfield[1+i], 1, physfield[1+i], 1,
                             mach,           1, mach,           1);
        }

        Vmath::Vdiv(nq, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(nq, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vsqrt(nq, mach, 1, mach, 1);

        Vmath::Vdiv(nq, mach, 1, soundspeed,   1, mach, 1);
    }

    /**
     * @brief Compute the dynamic viscosity using the Sutherland's law
     * \f$ \mu = \mu_star * (T / T_star)^3/2 * (T_star + 110) / (T + 110) \f$,
     * where:   \mu_star = 1.7894 * 10^-5 Kg / (m * s)
     *          T_star   = 288.15 K
     *
     * @param physfield    Input physical field.
     * @param mu           The resulting dynamic viscosity.
     */
    void CompressibleFlowSystem::GetDynamicViscosity(
        const Array<OneD, const NekDouble> &temperature,
              Array<OneD,       NekDouble> &mu)
    {
        const int nPts    = temperature.num_elements();
        NekDouble mu_star = m_mu;
        NekDouble T_star  = m_pInf / (m_rhoInf * m_gasConstant);
        NekDouble ratio;

        for (int i = 0; i < nPts; ++i)
        {
            ratio = temperature[i] / T_star;
            mu[i] = mu_star * ratio * sqrt(ratio) *
                    (T_star + 110.0) / (temperature[i] + 110.0);
        }
    }

    /**
     * @brief Perform post-integration checks, presently just to check steady
     * state behaviour.
     */
    bool CompressibleFlowSystem::v_PostIntegrate(int step)
    {
        if (m_steadyStateTol > 0.0)
        {
            bool doOutput = step % m_infosteps == 0;
            if (CalcSteadyState(doOutput))
            {
                if (m_comm->GetRank() == 0)
                {
                    cout << "Reached Steady State to tolerance "
                         << m_steadyStateTol << endl;
                }
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Calculate whether the system has reached a steady state by
     * observing residuals to a user-defined tolerance.
     */
    bool CompressibleFlowSystem::CalcSteadyState(bool output)
    {
        const int nPoints = GetTotPoints();
        const int nFields = m_fields.num_elements();

        // Holds L2 errors.
        Array<OneD, NekDouble> L2      (nFields);
        Array<OneD, NekDouble> residual(nFields);

        for (int i = 0; i < nFields; ++i)
        {
            Array<OneD, NekDouble> diff(nPoints);

            Vmath::Vsub(nPoints, m_fields[i]->GetPhys(), 1, m_un[i], 1, diff, 1);
            Vmath::Vmul(nPoints, diff, 1, diff, 1, diff, 1);
            residual[i] = Vmath::Vsum(nPoints, diff, 1);
        }

        m_comm->AllReduce(residual, LibUtilities::ReduceSum);

        // L2 error
        L2[0] = sqrt(residual[0]) / m_rhoInf;

        for (int i = 1; i < nFields-1; ++i)
        {
            L2[i] = sqrt(residual[i]) / m_UInf / m_rhoInf;
        }

        NekDouble Einf = m_pInf / (m_gamma-1.0) + 0.5 * m_rhoInf * m_UInf;
        L2[nFields-1] = sqrt(residual[nFields-1]) / Einf;

        if (m_comm->GetRank() == 0 && output)
        {
            // Output time
            m_errFile << setprecision(8) << setw(17) << scientific << m_time;

            // Output residuals
            for (int i = 0; i < nFields; ++i)
            {
                m_errFile << setprecision(11) << setw(22) << scientific
                          << L2[i];
            }

            m_errFile << endl;
        }

        // Calculate maximum L2 error
        NekDouble maxL2 = Vmath::Vmax(nFields, L2, 1);

        if (m_session->DefinesCmdLineArgument("verbose") &&
            m_comm->GetRank() == 0 && output)
        {
            cout << "-- Maximum L^2 residual: " << maxL2 << endl;
        }

        if (maxL2 <= m_steadyStateTol)
        {
            return true;
        }

        return false;
    }

    /**
     * @brief Calculate entropy.
     */
    void CompressibleFlowSystem::GetEntropy(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        const Array<OneD, const NekDouble>               &pressure,
        const Array<OneD, const NekDouble>               &temperature,
              Array<OneD,       NekDouble>               &entropy)
    {
        NekDouble entropy_sum = 0.0;
        const int npts = m_fields[0]->GetTotPoints();
        const NekDouble temp_inf = m_pInf/(m_rhoInf*m_gasConstant);;
        Array<OneD, NekDouble> L2entropy(npts, 0.0);

        for (int i = 0; i < npts; ++i)
        {
            entropy[i] = m_gamma / (m_gamma - 1.0) * m_gasConstant *
                            log(temperature[i]/temp_inf) - m_gasConstant *
                            log(pressure[i] / m_pInf);
        }

        Vmath::Vmul(npts, entropy, 1, entropy, 1, L2entropy, 1);

        entropy_sum = Vmath::Vsum(npts, L2entropy, 1);

        entropy_sum = sqrt(entropy_sum);

        std::ofstream m_file( "L2entropy.txt", std::ios_base::app);

        m_file << setprecision(16) << scientific << entropy_sum << endl;
        //m_file << Vmath::Vmax(entropy.num_elements(),entropy,1) << endl;

        m_file.close();
    }

    /**
     * @brief Calculate the maximum timestep subject to CFL restrictions.
     */
    NekDouble CompressibleFlowSystem::v_GetTimeStep(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray)
    {
        int n;
        int nElements = m_fields[0]->GetExpSize();
        const Array<OneD, int> ExpOrder = GetNumExpModesPerExp();

        Array<OneD, NekDouble> tstep      (nElements, 0.0);
        Array<OneD, NekDouble> stdVelocity(nElements);

        // Get standard velocity to compute the time-step limit
        GetStdVelocity(inarray, stdVelocity);

        // Factors to compute the time-step limit
        NekDouble minLength = 0.0;
        NekDouble alpha     = MaxTimeStepEstimator();
        NekDouble cLambda   = 0.2; // Spencer book-317

        // Loop over elements to compute the time-step limit for each element
        for(n = 0; n < nElements; ++n)
        {
            int npoints = m_fields[0]->GetExp(n)->GetTotPoints();
            Array<OneD, NekDouble> one2D(npoints, 1.0);
            NekDouble Area = m_fields[0]->GetExp(n)->Integral(one2D);

            minLength = sqrt(Area);
            if (m_fields[0]->GetExp(n)->as<LocalRegions::TriExp>())
            {
                minLength *= 2.0;
            }

            tstep[n] = m_cflSafetyFactor * alpha * minLength
                     / (stdVelocity[n] * cLambda
                        * (ExpOrder[n] - 1) * (ExpOrder[n] - 1));
        }

        // Get the minimum time-step limit and return the time-step
        NekDouble TimeStep = Vmath::Vmin(nElements, tstep, 1);
        m_comm->AllReduce(TimeStep, LibUtilities::ReduceMin);
        return TimeStep;
    }

    /**
     * @brief Set up logic for residual calculation.
     */
    void CompressibleFlowSystem::v_SetInitialConditions(
        NekDouble initialtime,
        bool      dumpInitialConditions,
        const int domain)
    {
        if (m_session->DefinesParameter("SteadyStateTol"))
        {
            const int nPoints = m_fields[0]->GetTotPoints();
            m_un = Array<OneD, Array<OneD, NekDouble> > (
                m_fields.num_elements());

            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                cout << Vmath::Vmax(nPoints, m_fields[i]->GetPhys(), 1) << endl;
                m_un[i] = Array<OneD, NekDouble>(nPoints);
                Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1, m_un[i], 1);
            }

            if (m_comm->GetRank() == 0)
            {
                std::string fName = m_session->GetSessionName() +
                    std::string(".res");
                m_errFile.open(fName.c_str());
                m_errFile << "# "
                          << setw(15) << left << "Time"
                          << setw(22) << left << "rho";

                std::string velFields[3] = {"u", "v", "w"};

                for (int i = 0; i < m_fields.num_elements()-2; ++i)
                {
                    m_errFile << setw(22) << "rho"+velFields[i];
                }

                m_errFile << setw(22) << left << "E" << endl;
            }
        }
    }


    /**
     * @brief Compute the advection velocity in the standard space
     * for each element of the expansion.
     *
     * @param inarray    Momentum field.
     * @param stdV       Standard velocity field.
     */
    void CompressibleFlowSystem::GetStdVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,                   NekDouble>   &stdV)
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int nBCEdgePts           = 0;

        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_spacedim);
        Array<OneD, NekDouble>               pressure   (nTotQuadPoints);
        Array<OneD, NekDouble>               soundspeed (nTotQuadPoints);
        LibUtilities::PointsKeyVector        ptsKeys;

        // Zero output array
        Vmath::Zero(stdV.num_elements(), stdV, 1);

        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity   [i] = Array<OneD, NekDouble>(nTotQuadPoints);
            stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }

        GetVelocityVector(inarray, velocity);
        GetPressure      (inarray, velocity, pressure);
        GetSoundSpeed    (inarray, pressure, soundspeed);

        for(int el = 0; el < n_element; ++el)
        {
            ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();

            // Possible bug: not multiply by jacobian??
            const SpatialDomains::GeomFactorsSharedPtr metricInfo =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo();
            const Array<TwoD, const NekDouble> &gmat =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()
                                                  ->GetDerivFactors(ptsKeys);

            int nq = m_fields[0]->GetExp(el)->GetTotPoints();

            if(metricInfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // d xi/ dx = gmat = 1/J * d x/d xi
                for (int i = 0; i < m_spacedim; ++i)
                {
                    Vmath::Vmul(nq, gmat[i], 1, velocity[0], 1,
                                stdVelocity[i], 1);
                    for (int j = 1; j < m_spacedim; ++j)
                    {
                        Vmath::Vvtvp(nq, gmat[m_spacedim*j+i], 1, velocity[j],
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            else
            {
                for (int i = 0; i < m_spacedim; ++i)
                {
                    Vmath::Smul(nq, gmat[i][0], velocity[0], 1,
                                stdVelocity[i], 1);
                    for (int j = 1; j < m_spacedim; ++j)
                    {
                        Vmath::Svtvp(nq, gmat[m_spacedim*j+i][0], velocity[j],
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }

            for (int i = 0; i < nq; ++i)
            {
                NekDouble pntVelocity = 0.0;
                for (int j = 0; j < m_spacedim; ++j)
                {
                    pntVelocity += stdVelocity[j][i]*stdVelocity[j][i];
                }
                pntVelocity = sqrt(pntVelocity) + soundspeed[nBCEdgePts];
                if (pntVelocity > stdV[el])
                {
                    stdV[el] = pntVelocity;
                }
                nBCEdgePts++;
            }
        }
    }

    /**
     * @brief Set the denominator to compute the time step when a cfl
     * control is employed. This function is no longer used but is still
     * here for being utilised in the future.
     *
     * @param n   Order of expansion element by element.
     */
    NekDouble CompressibleFlowSystem::GetStabilityLimit(int n)
    {
        ASSERTL0(n <= 20, "Illegal modes dimension for CFL calculation "
                          "(P has to be less then 20)");

        NekDouble CFLDG[21] = {  2.0000,   6.0000,  11.8424,  19.1569,
                                27.8419,  37.8247,  49.0518,  61.4815,
                                75.0797,  89.8181, 105.6700, 122.6200,
                               140.6400, 159.7300, 179.8500, 201.0100,
                               223.1800, 246.3600, 270.5300, 295.6900,
                               321.8300}; //CFLDG 1D [0-20]
        NekDouble CFL = 0.0;

        if (m_projectionType == MultiRegions::eDiscontinuous)
        {
            CFL = CFLDG[n];
        }
        else
        {
            ASSERTL0(false, "Continuous Galerkin stability coefficients "
                            "not introduced yet.");
        }

        return CFL;
    }

    /**
     * @brief Compute the vector of denominators to compute the time step
     * when a cfl control is employed. This function is no longer used but
     * is still here for being utilised in the future.
     *
     * @param ExpOrder   Order of expansion element by element.
     */
    Array<OneD, NekDouble> CompressibleFlowSystem::GetStabilityLimitVector(
        const Array<OneD,int> &ExpOrder)
    {
        int i;
        Array<OneD,NekDouble> returnval(m_fields[0]->GetExpSize(), 0.0);
        for (i =0; i<m_fields[0]->GetExpSize(); i++)
        {
            returnval[i] = GetStabilityLimit(ExpOrder[i]);
        }
        return returnval;
    }

    void CompressibleFlowSystem::GetSensor(
        const Array<OneD, const Array<OneD, NekDouble> > &physarray,
              Array<OneD,                   NekDouble>   &Sensor,
              Array<OneD,                   NekDouble>   &SensorKappa)
    {
        int e, NumModesElement, nQuadPointsElement;
        int nTotQuadPoints  = GetTotPoints();
        int nElements       = m_fields[0]->GetExpSize();

        // Find solution (SolP) at p = P;
        // The input array (physarray) is the solution at p = P;

        Array<OneD,int> ExpOrderElement = GetNumExpModesPerExp();

        Array<OneD, NekDouble> SolP    (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> SolPmOne(nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> SolNorm (nTotQuadPoints, 0.0);

        Vmath::Vcopy(nTotQuadPoints, physarray[0], 1, SolP, 1);

        int CoeffsCount = 0;

        for (e = 0; e < nElements; e++)
        {
            NumModesElement        = ExpOrderElement[e];
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            int nCoeffsElement     = m_fields[0]->GetExp(e)->GetNcoeffs();
            int numCutOff          = NumModesElement - 1;

            // Set-up of the Orthogonal basis for a Quadrilateral element which
            // is needed to obtain thesolution at P =  p - 1;

            Array<OneD, NekDouble> SolPElementPhys  (nQuadPointsElement, 0.0);
            Array<OneD, NekDouble> SolPElementCoeffs(nCoeffsElement,     0.0);

            Array<OneD, NekDouble> SolPmOneElementPhys(nQuadPointsElement, 0.0);
            Array<OneD, NekDouble> SolPmOneElementCoeffs(nCoeffsElement, 0.0);

            // create vector the save the solution points per element at P = p;

            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPElementPhys[i] = SolP[CoeffsCount+i];
            }

            m_fields[0]->GetExp(e)->FwdTrans(SolPElementPhys,
                                             SolPElementCoeffs);

            // ReduceOrderCoeffs reduces the polynomial order of the solution
            // that is represented by the coeffs given as an inarray. This is
            // done by projecting the higher order solution onto the orthogonal
            // basis and padding the higher order coefficients with zeros.

            m_fields[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff,
                                                      SolPElementCoeffs,
                                                      SolPmOneElementCoeffs);

            m_fields[0]->GetExp(e)->BwdTrans(SolPmOneElementCoeffs,
                                             SolPmOneElementPhys);

            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPmOne[CoeffsCount+i] = SolPmOneElementPhys[i];
            }

            NekDouble SolPmeanNumerator   = 0.0;
            NekDouble SolPmeanDenumerator = 0.0;

            // Determining the norm of the numerator of the Sensor

            Vmath::Vsub(nQuadPointsElement,
                        SolPElementPhys, 1,
                        SolPmOneElementPhys, 1,
                        SolNorm, 1);

            Vmath::Vmul(nQuadPointsElement,
                        SolNorm, 1,
                        SolNorm, 1,
                        SolNorm, 1);

            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPmeanNumerator   += SolNorm[i];
                SolPmeanDenumerator += SolPElementPhys[i];
            }

            for (int i = 0; i < nQuadPointsElement; ++i)
            {
                Sensor[CoeffsCount+i] =
                    sqrt(SolPmeanNumerator / nQuadPointsElement) /
                    sqrt(SolPmeanDenumerator / nQuadPointsElement);

                Sensor[CoeffsCount+i] = log10(Sensor[CoeffsCount+i]);
            }
            CoeffsCount += nQuadPointsElement;
        }

        CoeffsCount = 0.0;

        for (e = 0; e < nElements; e++)
        {
            NumModesElement    = ExpOrderElement[e];
            NekDouble ThetaS   = m_mu0;
            NekDouble Phi0     = m_Skappa;
            NekDouble DeltaPhi = m_Kappa;
            nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();

            for (int i = 0; i < nQuadPointsElement; i++)
            {
                if (Sensor[CoeffsCount+i] <= (Phi0 - DeltaPhi))
                {
                    SensorKappa[CoeffsCount+i] = 0;
                }
                else if(Sensor[CoeffsCount+i] >= (Phi0 + DeltaPhi))
                {
                    SensorKappa[CoeffsCount+i] = ThetaS;
                }
                else if(abs(Sensor[CoeffsCount+i]-Phi0) < DeltaPhi)
                {
                    SensorKappa[CoeffsCount+i] =
                        ThetaS / 2 * (1 + sin(M_PI * (Sensor[CoeffsCount+i] -
                                                      Phi0) / (2 * DeltaPhi)));
                }
            }

            CoeffsCount += nQuadPointsElement;
        }

    }

    void CompressibleFlowSystem::GetForcingTerm(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > outarrayForcing)
    {
        const int nPts = m_fields[0]->GetTotPoints();
        const int nvariables = m_fields.num_elements();
        const int nElements = m_fields[0]->GetExpSize();

        Array<OneD,  NekDouble>  Sensor(nPts, 0.0);
        Array<OneD,  NekDouble>  SensorKappa(nPts, 0.0);
        Array <OneD, NekDouble > Lambda(nPts, 0.0);
        Array <OneD, NekDouble > Tau(nPts, 1.0);
        Array <OneD, NekDouble > soundspeed(nPts, 0.0);
        Array <OneD, NekDouble > pressure(nPts, 0.0);
        Array <OneD, NekDouble > temperature(nPts, 0.0);
        Array <OneD, NekDouble > absVelocity(nPts, 0.0);
        Array <OneD, NekDouble > hel(nPts, 0.0);
        Array <OneD, NekDouble > h_minmin(m_spacedim, 0.0);

        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();
        Array<OneD, NekDouble> pOrder (nPts, 0.0);

        // Thermodynamic related quantities
        GetPressure(inarray, pressure);
        GetTemperature(inarray, pressure, temperature);
        GetSoundSpeed(inarray, pressure, soundspeed);
        GetAbsoluteVelocity(inarray, absVelocity);
        GetSensor(inarray, Sensor, SensorKappa);

        // Determine the maximum wavespeed
        Vmath::Vadd(nPts, absVelocity, 1, soundspeed, 1, Lambda, 1);

        NekDouble LambdaMax = Vmath::Vmax(nPts, Lambda, 1);

        int PointCount = 0;

        for (int e = 0; e < nElements; e++)
        {
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();

            for (int n = 0; n < nQuadPointsElement; n++)
            {
                pOrder[n + PointCount] = pOrderElmt[e];

                // order 1.0e-06
                Tau[n + PointCount] =
                    1.0 / (m_C1*pOrder[n + PointCount]*LambdaMax);

                outarrayForcing[nvariables-1][n + PointCount] =
                    1 / Tau[n + PointCount] * (m_hFactor * LambdaMax /
                                        pOrder[n + PointCount] *
                                        SensorKappa[n + PointCount] -
                                        inarray[nvariables-1][n + PointCount]);
            }
            PointCount += nQuadPointsElement;
        }
    }

    void CompressibleFlowSystem::GetElementDimensions(
        Array<OneD,       Array<OneD, NekDouble> > &outarray,
        Array<OneD,                   NekDouble >  &hmin)
    {
        // So far, this function is only implemented for quads
        const int nElements = m_fields[0]->GetExpSize();

        SpatialDomains::HexGeomSharedPtr    ElHexGeom;

        NekDouble hx = 0.0;
        NekDouble hy = 0.0;

        for (int e = 0; e < nElements; e++)
        {
            NekDouble nedges = m_fields[0]->GetExp(e)->GetNedges();
            Array <OneD, NekDouble> L1(nedges, 0.0);

            for (int j = 0; j < nedges; ++j)
            {

                NekDouble x0 = 0.0;
                NekDouble y0 = 0.0;
                NekDouble z0 = 0.0;

                NekDouble x1 = 0.0;
                NekDouble y1 = 0.0;
                NekDouble z1 = 0.0;

                if (boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                                            m_fields[0]->GetExp(e)->GetGeom()))
                {
                    SpatialDomains::QuadGeomSharedPtr ElQuadGeom =
                        boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                                            m_fields[0]->GetExp(e)->GetGeom());

                    ElQuadGeom->GetEdge(j)->GetVertex(0)->GetCoords(x0, y0, z0);
                    ElQuadGeom->GetEdge(j)->GetVertex(1)->GetCoords(x1, y1, z1);

                    L1[j] = sqrt(pow((x0-x1),2)+pow((y0-y1),2)+pow((z0-z1),2));
                }
                else
                {
                    ASSERTL0(false, "GetElementDimensions() is only "
                                    "implemented for quadrilateral elements");
                }
            }
            // determine the minimum length in x and y direction
            // still have to find a better estimate when dealing
            // with unstructured meshes
            if(boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                                            m_fields[0]->GetExp(e)->GetGeom()))
            {
                hx = min(L1[0], L1[2]);
                hy = min(L1[1], L1[3]);

                outarray[0][e] = hx;
                outarray[1][e] = hy;
            }
        }

        hmin[0] = Vmath::Vmin(outarray[0].num_elements(), outarray[0], 1);
        hmin[1] = Vmath::Vmin(outarray[1].num_elements(), outarray[1], 1);
    }

    void CompressibleFlowSystem::GetAbsoluteVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,                   NekDouble>   &Vtot)
    {
        int nTotQuadPoints = GetTotPoints();

        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);

        Vmath::Zero(Vtot.num_elements(), Vtot, 1);

        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nTotQuadPoints);
        }

        GetVelocityVector(inarray, velocity);

        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nTotQuadPoints,
                         velocity[i], 1,
                         velocity[i], 1,
                         Vtot, 1,
                         Vtot, 1);
        }

        Vmath::Vsqrt(Vtot.num_elements(),Vtot,1,Vtot,1);
    }

    void CompressibleFlowSystem::GetSmoothArtificialViscosity(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD,             NekDouble  > &eps_bar)
    {
        int nvariables = physfield.num_elements();
        int nPts       = m_fields[0]->GetTotPoints();

        Array<OneD, NekDouble > pressure   (nPts, 0.0);
        Array<OneD, NekDouble > temperature(nPts, 0.0);
        Array<OneD, NekDouble > sensor     (nPts, 0.0);
        Array<OneD, NekDouble > SensorKappa(nPts, 0.0);
        Array<OneD, NekDouble > absVelocity(nPts, 0.0);
        Array<OneD, NekDouble > soundspeed (nPts, 0.0);
        Array<OneD, NekDouble > Lambda     (nPts, 0.0);
        Array<OneD, NekDouble > mu_var     (nPts, 0.0);
        Array<OneD, NekDouble > h_minmin   (m_spacedim, 0.0);
        Vmath::Zero(nPts, eps_bar, 1);

        // Thermodynamic related quantities
        GetPressure(physfield, pressure);
        GetTemperature(physfield, pressure, temperature);
        GetSoundSpeed(physfield, pressure, soundspeed);
        GetAbsoluteVelocity(physfield, absVelocity);
        GetSensor(physfield, sensor, SensorKappa);

        // Determine the maximum wavespeed
        Vmath::Vadd(nPts, absVelocity, 1, soundspeed, 1, Lambda, 1);

        // Determine hbar = hx_i/h
        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();

        NekDouble ThetaH = m_FacH;
        NekDouble ThetaL = m_FacL;

        NekDouble Phi0     = (ThetaH+ThetaL)/2;
        NekDouble DeltaPhi = ThetaH-Phi0;

        Vmath::Zero(eps_bar.num_elements(), eps_bar, 1);

            /*Vmath::Smul(eps_bar.num_elements(),
                    m_eps_max,
                    &physfield[nvariables-1][0], 1,
                    &eps_bar[0], 1);*/

        for (int e = 0; e < eps_bar.num_elements(); e++)
        {
            if (physfield[nvariables-1][e] <= (Phi0 - DeltaPhi))
            {
                eps_bar[e] = 0;
            }
            else if(physfield[nvariables-1][e] >= (Phi0 + DeltaPhi))
            {
                eps_bar[e] = m_mu0;
            }
            else if(abs(physfield[nvariables-1][e]-Phi0) < DeltaPhi)
            {
                eps_bar[e] = m_mu0/2*(1+sin(M_PI*
                (physfield[nvariables-1][e]-Phi0)/(2*DeltaPhi)));
            }
        }

    }

    void CompressibleFlowSystem::GetArtificialDynamicViscosity(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD,             NekDouble  > &mu_var)
    {
        const int nElements = m_fields[0]->GetExpSize();
        int PointCount      = 0;
        int nTotQuadPoints  = GetTotPoints();

        Array<OneD, NekDouble> S_e        (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> se         (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> Sensor     (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> SensorKappa(nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> absVelocity(nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> soundspeed (nTotQuadPoints, 0.0);
        Array<OneD, NekDouble> pressure   (nTotQuadPoints, 0.0);

        GetAbsoluteVelocity(physfield, absVelocity);
        GetPressure        (physfield, pressure);
        GetSoundSpeed      (physfield, pressure, soundspeed);
        GetSensor          (physfield, Sensor, SensorKappa);

        Array<OneD, int> pOrderElmt = GetNumExpModesPerExp();
        Array<OneD, NekDouble> Lambda(nTotQuadPoints, 1.0);
        Vmath::Vadd(nTotQuadPoints, absVelocity, 1, soundspeed, 1, Lambda, 1);

        for (int e = 0; e < nElements; e++)
        {
            // Threshold value specified in C. Biottos thesis.  Based on a 1D
            // shock tube problem S_k = log10(1/p^4). See G.E. Barter and
            // D.L. Darmofal. Shock Capturing with PDE-based artificial
            // diffusion for DGFEM: Part 1 Formulation, Journal of Computational
            // Physics 229 (2010) 1810-1827 for further reference

            // Adjustable depending on the coarsness of the mesh. Might want to
            // move this variable into the session file

            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            Array <OneD, NekDouble> one2D(nQuadPointsElement, 1.0);

            for (int n = 0; n < nQuadPointsElement; n++)
            {
                NekDouble mu_0 = m_mu0;

                if (Sensor[n+PointCount] < (m_Skappa-m_Kappa))
                {
                    mu_var[n+PointCount] = 0;
                }
                else if (Sensor[n+PointCount] >= (m_Skappa-m_Kappa)
                        && Sensor[n+PointCount] <= (m_Skappa+m_Kappa))
                {
                    mu_var[n+PointCount] = mu_0 * (0.5 * (1 + sin(
                                           M_PI * (Sensor[n+PointCount] -
                                                   m_Skappa - m_Kappa) /
                                                            (2*m_Kappa))));
                }
                else if (Sensor[n+PointCount] > (m_Skappa+m_Kappa))
                {
                    mu_var[n+PointCount] = mu_0;
                }
            }

            PointCount += nQuadPointsElement;
        }
    }

    void CompressibleFlowSystem::SetVarPOrderElmt(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD,                   NekDouble  > &PolyOrder)
    {
        int e;
        NekDouble s_ds, s_sm, s_fl;

        int nElements = m_fields[0]->GetExpSize();
        int npts      = m_fields[0]->GetTotPoints();

        Array<OneD, NekDouble > Sensor     (npts, 0.0);
        Array<OneD, NekDouble > SensorKappa(npts, 0.0);
        Array<OneD, NekDouble > se         (npts, 0.0);

        GetSensor(physfield, Sensor, SensorKappa);

        int nQuadPointsElement  = 0;
        int npCount             = 0;
        int MinOrder            = 2;
        int MaxOrder            = 12;
        int MinOrderShock       = 4;

        std::ofstream m_file( "VariablePComposites.txt", std::ios_base::app);
        for (int e = 0; e < nElements; e++)
        {
            m_file << "<C ID=\"" << e+1 << "\"> Q[" << e << "] </C>"<< endl;
        }
        m_file.close();

        std::ofstream m_file2( "VariablePExpansions.txt", std::ios_base::app);

        for (e = 0; e < nElements; e++)
        {
            nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();

            // Define thresholds
            // Ideally, these threshold values could be given
            // in the Session File
            s_ds =  -5.0;
            //s_ds = s_0*log10(PolyOrder[e]);
            s_sm = -6;
            s_fl = -7;


            for (int i = 0; i < nQuadPointsElement; i++)
            {
                se[npCount + i] = (Sensor[npCount + i]);

                if (se[npCount + i] > s_ds)
                {
                    if (PolyOrder[npCount + i] > MinOrderShock)
                    {
                        PolyOrder[npCount + i] = PolyOrder[npCount + i] - 1;
                    }
                    else if(PolyOrder[e] < MinOrderShock)
                    {
                        PolyOrder[npCount + i] = PolyOrder[npCount + i] + 1;
                    }

                }
                else if (se[npCount + i] > s_sm && se[npCount + i] < s_ds)
                {
                    if (PolyOrder[npCount + i] < MaxOrder)
                    {
                        PolyOrder[npCount + i] = PolyOrder[npCount + i] + 2;
                    }
                }
                else if (se[npCount + i] > s_fl && se[npCount + i] < s_sm)
                {
                    PolyOrder[npCount + i] = PolyOrder[npCount + i] + 1;
                }
                else if (se[npCount + i] < s_fl)
                {
                    if (PolyOrder[npCount + i] > MinOrder)
                    {
                            PolyOrder[npCount + i] = PolyOrder[npCount + i];
                    }
                }
            }
            m_file2 << "<E COMPOSITE= \"C[" << e+1
                    << "]\" NUMMODES=\"" << PolyOrder[npCount + 1]
                    << "\" TYPE=\"MODIFIED\" FIELDS=\"rho,rhou,rhov,rhow,E\" />"
                    << endl;
            npCount += nQuadPointsElement;
        }

        m_file2.close();
    }

    void CompressibleFlowSystem::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();
        Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.num_elements());

        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            tmp[i] = m_fields[i]->GetPhys();
        }

        Array<OneD, NekDouble> pressure(nPhys), soundspeed(nPhys), mach(nPhys);
        Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys), smooth(nPhys);

        GetPressure  (tmp, pressure);
        GetSoundSpeed(tmp, pressure, soundspeed);
        GetMach      (tmp, soundspeed, mach);
        GetSensor    (tmp, sensor, SensorKappa);
        GetSmoothArtificialViscosity    (tmp, smooth);

        Array<OneD, NekDouble> pFwd(nCoeffs), sFwd(nCoeffs), mFwd(nCoeffs);
        Array<OneD, NekDouble> sensFwd(nCoeffs), smoothFwd(nCoeffs);

        m_fields[0]->FwdTrans(pressure,   pFwd);
        m_fields[0]->FwdTrans(soundspeed, sFwd);
        m_fields[0]->FwdTrans(mach,       mFwd);
        m_fields[0]->FwdTrans(sensor,     sensFwd);
        m_fields[0]->FwdTrans(smooth,     smoothFwd);

        variables.push_back  ("p");
        variables.push_back  ("a");
        variables.push_back  ("Mach");
        variables.push_back  ("Sensor");
        variables.push_back  ("SmoothVisc");
        fieldcoeffs.push_back(pFwd);
        fieldcoeffs.push_back(sFwd);
        fieldcoeffs.push_back(mFwd);
        fieldcoeffs.push_back(sensFwd);
        fieldcoeffs.push_back(smoothFwd);
    }
}
