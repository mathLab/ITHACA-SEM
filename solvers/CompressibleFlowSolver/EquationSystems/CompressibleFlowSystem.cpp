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
        m_velLoc = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_velLoc[i] = i+1;
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

        // Get vInf parameter from session file.
        if (m_spacedim == 2 || m_spacedim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("vInf"),
                     "Compressible flow sessions must define a vInf parameter"
                     "for 2D/3D problems.");
            m_session->LoadParameter("vInf", m_vInf, 0.0);
        }

        // Get wInf parameter from session file.
        if (m_spacedim == 3)
        {
            ASSERTL0(m_session->DefinesParameter("wInf"),
                     "Compressible flow sessions must define a wInf parameter"
                     "for 3D problems.");
            m_session->LoadParameter("wInf", m_wInf, 0.0);
        }

        m_session->LoadParameter ("GasConstant",   m_gasConstant,   287.058);
        m_session->LoadParameter ("Twall",         m_Twall,         300.15);
        m_session->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
        m_session->LoadSolverInfo("ShockCaptureType",   m_shockCaptureType, "Off");
        m_session->LoadParameter ("mu",            m_mu,            1.78e-05);
        m_session->LoadParameter ("Skappa",        m_Skappa,        -2.048);
        m_session->LoadParameter ("Kappa",         m_Kappa,         0.0);
        m_session->LoadParameter ("mu0",           m_mu0,           1.0);
        m_session->LoadParameter ("FL",            m_FacL,           0.0);
        m_session->LoadParameter ("FH",            m_FacH,           0.0);
        m_session->LoadParameter ("epsMax",            m_eps_max,           0.0);
        m_session->LoadParameter ("thermalConductivity",
                                  m_thermalConductivity, 0.0257);

        m_Cp      = m_gamma / (m_gamma - 1.0) * m_gasConstant;
        m_Prandtl = m_Cp * m_mu / m_thermalConductivity;

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

                if (m_shockCaptureType == "Off" && m_specHP_dealiasing)
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVectorDeAlias, this);
                    
                    m_diffusion->SetFluxVectorNS(
                                                 &CompressibleFlowSystem::
                                                 GetViscousFluxVectorDeAlias, this);
                }
                else if(m_shockCaptureType == "Off")
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVector, this);
                    
                    m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                                                 GetViscousFluxVector, this);
                }
                // Setting up flux vector for diffusion operator
                
                if (m_shockCaptureType == "Smooth")
                {
                    m_advection->SetFluxVector(&CompressibleFlowSystem::
                                               GetFluxVectorPDESC, this);
                    
                    m_diffusion->SetFluxVectorNS(&CompressibleFlowSystem::
                                            GetViscousFluxVectorPDESC, this);
                }
                if (m_shockCaptureType == "NonSmooth")
                {
                    m_diffusion->SetArtificialDiffusionVector(
                        &CompressibleFlowSystem::GetArtificialDynamicViscosity, this);
                }
                // Setting up Riemann solver for advection operator
                m_session->LoadSolverInfo("UpwindType", riemName, "Average");
                
                
                m_diffusion->SetArtificialDiffusionVector(&CompressibleFlowSystem::GetArtificialDynamicViscosity, this);

                m_riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                            .CreateInstance(riemName);

                // Setting up upwind solver for diffusion operator
                m_riemannSolverLDG = SolverUtils::GetRiemannSolverFactory()
                                                .CreateInstance("UpwindLDG");

                // Setting up parameters for advection operator Riemann solver
                m_riemannSolver->SetParam (
                    "gamma",  &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolver->SetAuxiliary(
                    "velLoc", &CompressibleFlowSystem::GetVelLoc,  this);
                m_riemannSolver->SetVector(
                    "N",      &CompressibleFlowSystem::GetNormals, this);

                // Setting up parameters for diffusion operator Riemann solver
                m_riemannSolverLDG->SetParam (
                    "gamma",  &CompressibleFlowSystem::GetGamma,   this);
                m_riemannSolverLDG->SetAuxiliary(
                    "velLoc", &CompressibleFlowSystem::GetVelLoc,  this);
                m_riemannSolverLDG->SetVector(
                    "N",      &CompressibleFlowSystem::GetNormals, this);

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

        // Adjust the physical values of the trace to
        // take user defined boundaries into account
        int e, id1, id2, nBCEdgePts, eMax;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

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
    
    void CompressibleFlowSystem::ArtificialViscosityBC(
                    int                                   bcRegion,
                    int                                   cnt,
                    Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nTracePts = GetTraceTotPoints();
        int nVariables = physarray.num_elements();
        const int nElements  = m_fields[0]->GetExpSize();
        const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();
    
        // Get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
        for (i = 0; i < nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts,0.0);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }
        
        int e, id1, id2, nBCEdgePts, eMax;
        
        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
    
        //Vmath::Zero(nTracePts, Fwd[nVariables-1], 1);
        
        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);
            
            // Copy boundary adjusted values into the boundary expansion
            
            Vmath::Vcopy(nBCEdgePts, &Fwd[nVariables-1][id2], 1,
                         &(m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1], 1);
        }
    }

    /**
     * @brief Simmetry boundary conditions for compressible flow problems.
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

        int e, id1, id2, nBCEdgePts, eMax;

        eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

        for (e = 0; e < eMax; ++e)
        {
            nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

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
            Vmath::Smul(nTracePts, m_vInf, m_traceNormals[0], 1, tmp1, 1);
            Vmath::Vadd(nTracePts, VnInf, 1, tmp1, 1, VnInf, 1);
        }
        if (nDimensions == 3)
        {
            velInf[2] = m_wInf;
            Vmath::Smul(nTracePts, m_wInf, m_traceNormals[0], 1, tmp2, 1);
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
                GetExp(e)->GetNumPoints(0);

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
                GetExp(e)->GetNumPoints(0);
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
    }
    
    void CompressibleFlowSystem::GetFluxVectorPDESC(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = m_fields[0]->GetTotPoints();
        
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
        
        // Add a zero row for the advective fluxes
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Zero(nq, flux[m_spacedim+2][j], 1);
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
        Array<OneD, NekDouble > mu    (nPts, 0.0);
        Array<OneD, NekDouble > mu2   (nPts, 0.0);
        Array<OneD, NekDouble > divVel(nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            GetDynamicViscosity(physfield[nVariables-2], mu);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, &mu[0], 1);
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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][1][0], 1,
                        &tmp1[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][2][0], 1,
                        &tmp2[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[1][2][0], 1,
                        &tmp2[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][3][0], 1,
                        &tmp3[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[1][3][0], 1,
                        &tmp3[0], 1);

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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[2][3][0], 1,
                        &tmp3[0], 1);

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
    
    void CompressibleFlowSystem::GetViscousFluxVectorPDESC(
      const Array<OneD, Array<OneD, NekDouble> >               &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j, k;
        int nvariables = m_fields.num_elements();
        int nPts       = m_fields[0]->GetTotPoints();
        
        // Stokes hypotesis
        NekDouble lambda = -0.66666;
        NekDouble C1     = 3.0;
        NekDouble C1C2   = 0.075;
        // Auxiliary variables
        Array<OneD, NekDouble > mu                  (nPts, 0.0);
        Array<OneD, NekDouble > mu2                 (nPts, 0.0);
        Array<OneD, NekDouble > divVel              (nPts, 0.0);
        Array<OneD, NekDouble > pressure            (nPts, 0.0);
        Array<OneD, NekDouble > temperature         (nPts, 0.0);
        Array <OneD, NekDouble > sensor             (nPts, 0.0);
        Array <OneD, NekDouble > SensorKappa        (nPts, 0.0);
        Array <OneD, NekDouble > absVelocity        (nPts, 0.0);
        Array <OneD, NekDouble > soundspeed         (nPts, 0.0);
        Array <OneD, NekDouble > Lambda             (nPts, 0.0);
        Array <OneD, NekDouble > mu_var             (nPts, 0.0);
        
        // Set up wrapper to fields data storage
        Array<OneD, Array<OneD, NekDouble> > fields(nvariables);
        
        // Reorder storage to list time-integrated fields first
        for (i = 0; i < nvariables; ++i)
        {
            fields[i] = m_fields[i]->UpdatePhys();
        }
        // Thermodynamic related quantities
        GetPressure(fields, pressure);
        GetTemperature(fields, pressure, temperature);
        GetSoundSpeed(fields, pressure, soundspeed);
        GetAbsoluteVelocity(fields, absVelocity);
        GetSensor(fields, sensor, SensorKappa);
        
        // Determine the maximum wavespeed
        Vmath::Vadd(nPts, absVelocity, 1, soundspeed, 1, Lambda, 1);
        NekDouble LambdaMax = Vmath::Vmax(nPts, Lambda, 1);
    
        // Determine hbar = hx_i/h
        const int nElements  = m_fields[0]->GetExpSize();
        
        // PUT CODE TO DETERMINE hx_i AND min(hx_i)
        //
        Array <OneD, Array <OneD, NekDouble > > h_av(m_spacedim);
        Array <OneD, NekDouble> h_minmin(m_spacedim, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            h_av[i] = Array <OneD, NekDouble > (nPts, 0.0);
        }
        Array <OneD, Array <OneD, NekDouble > > ElDim(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            ElDim[i] = Array <OneD, NekDouble > (nElements, 0.0);
        }
        
        GetElementDimensions(ElDim, h_minmin);
        
        NekDouble hxmin = 0.0;
        NekDouble hymin = 0.0;
        NekDouble hzmin = 0.0;
        
        NekDouble h_mean_sumx = 0.0;
        NekDouble h_mean_sumy = 0.0;
        NekDouble h_mean_sumz = 0.0;
        
        NekDouble h_mean_x = 0.0;
        NekDouble h_mean_y = 0.0;
        NekDouble h_mean_z = 0.0;
        NekDouble h_mean   = 0.0;
        
        NekDouble hminmin  = 0.0;
        
        hxmin = Vmath::Vmin(ElDim[0].num_elements(), &ElDim[0][0], 1);
        hymin = Vmath::Vmin(ElDim[1].num_elements(), &ElDim[1][0], 1);
        int PointCount = 0.0;
        
        if (m_spacedim == 2)
        {
            
        
            for (int i = 0; i < nElements; ++i)
            {
                h_mean_sumx += ElDim[0][i];
                h_mean_sumy += ElDim[1][i];
            }
            
            h_mean_x = h_mean_sumx/nElements;
            h_mean_y = h_mean_sumy/nElements;
            
            h_mean = (h_mean_x+h_mean_y)/2;
        }
        
        if (m_spacedim == 3)
        {
            for (int e = 0; e < nElements; e++)
            {
                int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
                
                for (int n = 0; n < nQuadPointsElement; n++)
                {
                    h_av[0][n + PointCount] = ElDim[0][e];
                    h_av[1][n + PointCount] = ElDim[1][e];
                    h_av[2][n + PointCount] = ElDim[2][e];
                }
                
                PointCount += nQuadPointsElement;
            }
            
            for (int i = 0; i < nElements; ++i)
            {
                h_mean_sumx += h_av[0][i];
                h_mean_sumy += h_av[1][i];
                h_mean_sumz += h_av[2][i];
            }
            
            h_mean_x = h_mean_sumx/nElements;
            h_mean_y = h_mean_sumy/nElements;
            h_mean_z = h_mean_sumy/nElements;
            
            h_mean = (h_mean_x+h_mean_y+h_mean_z)/3;
        }
        
        //
        
        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();
        Array<OneD, NekDouble> pOrder (nPts, 0.0);
        // Based on eps, determine eps_bar
        
        PointCount = 0;
        
        Array<OneD, NekDouble> eps_bar(nPts, 0.0);
        
        Vmath::Zero(nPts, eps_bar, 1);
        
        GetSmoothArtificialViscosity(fields, eps_bar);
        
        for (int e = 0; e < nElements; e++)
        {
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            
            for (int n = 0; n < nQuadPointsElement; n++)
            {
                pOrder[n + PointCount] = pOrderElmt[e];
            }
            PointCount += nQuadPointsElement;
        }
        
        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            GetDynamicViscosity(physfield[nvariables-3], mu);
        }
        else
        {
            Vmath::Sadd(nPts, m_mu, &mu[0], 1, &mu[0], 1);
        }
        
        // Computing diagonal terms of viscous stress tensor
        Array<OneD, Array<OneD, NekDouble> > tmp(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > tmpar(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > tmpar1(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > Sgg(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > Seps(m_spacedim);
        // mu2 = 2 * mu
        Vmath::Smul(nPts, 2.0, &mu[0], 1, &mu2[0], 1);
        
        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts,
                        &divVel[0], 1,
                        &derivativesO1[j][j][0], 1,
                        &divVel[0], 1);
        }
        
        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, &divVel[0], 1, &divVel[0], 1);
        Vmath::Vmul(nPts, &mu[0], 1, &divVel[0], 1, &divVel[0], 1);
        
        
        Array<OneD, Array<OneD, NekDouble> >vel     (m_spacedim);
        // Collect the velocity components
        for (i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble> (nPts, 0.0);
            Vmath::Vdiv(nPts, physfield[i+1], 1, physfield[0], 1, vel[i], 1);
        }
        
        // Apply chain rule to determine all mixed derivatives
        
        // Initiate the matrix containing the mixed derivatives
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > derivativesMix(m_spacedim);
        
        for (i = 0; i < m_spacedim; ++i)
        {
            derivativesMix[i] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim+1);
            
            for (j = 0; j < m_spacedim+1; ++j)
            {
                derivativesMix[i][j] = Array<OneD, NekDouble>(nPts, 0.0);
            }
        }
        
        for (i = 0; i < m_spacedim; ++i)
        {
            // mixed terms of velocity components and the energy
            for (j = 0; j < m_spacedim; ++j)
            {
                // d(rhou_i)/dx_i = u*drho/dx_j  + rho*du_i/dx_j
                // d(rhou)/dx = u*drho/dx
                Vmath::Vmul(nPts,
                            &vel[j][0], 1,
                            &derivativesO1[i][nvariables-2][0], 1,
                            &derivativesMix[i][j][0], 1);
                
                // d(rhou)/dx = u*drho/dx + rho*du/dx
                Vmath::Vvtvp(nPts,
                             &physfield[nvariables-2][0], 1,
                             &derivativesO1[i][j][0], 1,
                             &derivativesMix[i][j][0], 1,
                             &derivativesMix[i][j][0], 1);
            }
        }
        
        // Calculate the mixed term for the enthalpy
        
        for (i = 0; i < m_spacedim; ++i)
        {
            //d(rhoH)/dx = H*drho/dx
            Vmath::Vmul(nPts,
                        &physfield[nvariables][0], 1,
                        &derivativesO1[i][nvariables-2][0], 1,
                        &derivativesMix[i][m_spacedim][0], 1);
            
            // d(rhou)/dx = H*drho/dx + rho*dH/dx
            Vmath::Vvtvp(nPts,
                         &physfield[nvariables-2][0], 1,
                         &derivativesO1[i][m_spacedim][0], 1,
                         &derivativesMix[i][m_spacedim][0], 1,
                         &derivativesMix[i][m_spacedim][0], 1);
        }
        // Digonal terms of viscous stress tensor (Sxx, Syy, Szz)
        // Sjj = 2 * mu * du_j/dx_j - (2 / 3) * mu * sum_j(du_j/dx_j)
        for (j = 0; j < m_spacedim; ++j)
        {
            tmp[j]      = Array<OneD, NekDouble>(nPts, 0.0);
            tmpar[j]    = Array<OneD, NekDouble>(nPts, 0.0);
            Sgg[j]      = Array<OneD, NekDouble>(nPts, 0.0);
            
            Vmath::Vmul(nPts,
                        &mu2[0], 1,
                        &derivativesO1[j][j][0], 1,
                        &tmp[j][0], 1);
            
            Vmath::Vadd(nPts,
                        &tmp[j][0], 1,
                        &divVel[0], 1,
                        &Sgg[j][0], 1);
            
            // tmp = h_av*eps_bar
            Vmath::Vmul(nPts,
                        &h_av[j][0], 1,
                        &eps_bar[0], 1,
                        &tmpar[j][0], 1);
            
            // tmp = h_av*eps_bar*d(rhou_i)/dx_i
            Vmath::Vmul(nPts,
                        &eps_bar[0], 1,
                        &derivativesMix[j][j][0], 1,
                        &tmpar[j][0], 1);
        }
        
        // Determine terms for viscosity PDE C1C2 p lambda h_i^2/(mini(h_i) deps/dx_i
        
        for (j = 0; j < m_spacedim; ++j)
        {
            tmpar1[j] = Array<OneD, NekDouble>(nPts, 0.0);
            Seps[j] = Array<OneD, NekDouble>(nPts, 0.0);
            
            // tmp = h_i^2
            Vmath::Vmul(nPts,
                        &h_av[j][0], 1,
                        &h_av[j][0], 1,
                        &tmpar1[j][0], 1);
            
            // tmp = h_i^2/(min(h_i))
            // NekDouble hminimal      = min(h_minmin[0],h_minmin[1]);
            NekDouble hminimal_inv  = 1/h_minmin[j];
            
            Vmath::Smul(nPts, hminimal_inv,
                        &tmpar1[j][0], 1,
                        &tmpar1[j][0], 1);
            
            // tmp = LambdaMax * C1C2 * h_i^2/(min(h_i))
            
            Vmath::Smul(nPts, C1C2,
                        &tmpar1[j][0], 1,
                        &tmpar1[j][0], 1);
            
            // tmp = p * LambdaMax * C1C2 * h_i^2/(min(h_i))
            
            Vmath::Vmul(nPts,
                        &tmpar1[j][0], 1,
                        &pOrder[0], 1,
                        &tmpar1[j][0], 1);

            // tmp = p * C1C2 * h_i^2/(min(h_i))
            
            Vmath::Smul(nPts,
                        LambdaMax,
                        &tmpar1[j][0], 1,
                        &tmpar1[j][0], 1);
            
            // tmp = deps/dx_i * Lambda * p * C1C2 * h_i^2/(min(h_i))
            Vmath::Vmul(nPts,
                        &tmpar1[j][0], 1,
                        &derivativesO1[j][nvariables-1][0], 1,
                        &tmpar1[j][0], 1);
            
            
            Vmath::Vcopy(nPts, &tmpar1[j][0], 1, &Seps[j][0], 1);
            
            //Vmath::Vcopy(nPts,
            //             &derivativesO1[j][nvariables-1][0], 1,
            //             &Seps[j][0], 1);
        }
        
        // Extra diagonal terms of viscous stress tensor (Sxy, Sxz, Syz)
        // Note: they exist for 2D and 3D problems only
        Array<OneD, NekDouble > Sxy(nPts, 0.0);
        Array<OneD, NekDouble > Syx(nPts, 0.0);
        Array<OneD, NekDouble > Sxz(nPts, 0.0);
        Array<OneD, NekDouble > Syz(nPts, 0.0);
        
        Array<OneD, NekDouble > tmpar2(nPts, 0.0);
        Array<OneD, NekDouble > tmpar20(nPts, 0.0);
        if (m_spacedim == 2)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts,
                        &derivativesO1[0][1][0], 1,
                        &derivativesO1[1][0][0], 1,
                        &Sxy[0], 1);
            
            // Syx = (du/dx + du/dy)
            Vmath::Vadd(nPts,
                        &derivativesO1[1][0][0], 1,
                        &derivativesO1[0][1][0], 1,
                        &Syx[0], 1);
            
            // tmp = h_x*eps_bar
            Vmath::Vmul(nPts,
                        &h_av[0][0], 1,
                        &eps_bar[0], 1,
                        &tmpar2[0], 1);
            
            // tmp = h_x*eps_bar*d(rhov)/dx
            Vmath::Vmul(nPts,
                        &eps_bar[0], 1,
                        &derivativesMix[0][1][0], 1,
                        &tmpar2[0], 1);
            
            Vmath::Vmul(nPts,
                        &h_av[1][0], 1,
                        &eps_bar[0], 1,
                        &tmpar20[0], 1);
            
            // tmp = h_y*eps_bar*d(rhou)/dy
            Vmath::Vmul(nPts,
                        &eps_bar[0], 1,
                        &derivativesMix[1][0][0], 1,
                        &tmpar20[0], 1);
            
            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
            Vmath::Vmul(nPts, &mu[0], 1, &Syx[0], 1, &Syx[0], 1);
            // Sxy = mu * (du/dy + dv/dx) + h_x*eps_bar*d(rhov)/dx
            Vmath::Vadd(nPts, &tmpar2[0], 1, &Sxy[0], 1, &Sxy[0], 1);
            // Syx = mu * (du/dy + dv/dx) + h_y*eps_bar*d(rhou)/dy
            Vmath::Vadd(nPts, &tmpar20[0], 1, &Syx[0], 1, &Syx[0], 1);
        }
        else if (m_spacedim == 3)
        {
            // Sxy = (du/dy + dv/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][1][0], 1,
                        &derivativesO1[1][0][0], 1, &Sxy[0], 1);
            
            Vmath::Zero(nPts, &tmpar2[0], 1);
            
            // tmp = h_x*eps_bar
            Vmath::Vmul(nPts,
                        &h_av[0][0], 1,
                        &eps_bar[0], 1,
                        &tmpar2[0], 1);
            
            // tmp = h_x*eps_bar*d(rhov)/dx
            Vmath::Vmul(nPts,
                        &eps_bar[0], 1,
                        &derivativesMix[0][1][0], 1,
                        &tmpar2[0], 1);
            
            // Sxy = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);
            // Sxy = mu * (du/dy + dv/dx) + h_x*eps_bar*d(rhov)/dx
            Vmath::Vadd(nPts, &tmpar2[0], 1, &Sxy[0], 1, &Sxy[0], 1);
            
            // Sxz = (du/dz + dw/dx)
            Vmath::Vadd(nPts, &derivativesO1[0][2][0], 1,
                        &derivativesO1[2][0][0], 1, &Sxz[0], 1);
            
            Vmath::Zero(nPts, &tmpar2[0], 1);
            
            // tmp = h_x*eps_bar
            Vmath::Vmul(nPts,
                        &h_av[0][0], 1,
                        &eps_bar[0], 1,
                        &tmpar2[0], 1);
            
            // tmp = h_x*eps_bar*d(rhow)/dx
            Vmath::Vmul(nPts,
                        &eps_bar[0], 1,
                        &derivativesMix[0][2][0], 1,
                        &tmpar2[0], 1);
            
            // Sxz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Sxz[0], 1, &Sxz[0], 1);
            // Sxz = mu * (du/dy + dv/dx) + h_x*eps_bar*d(rhow)/dx
            Vmath::Vadd(nPts, &tmpar2[0], 1, &Sxy[0], 1, &Sxy[0], 1);
            
            // Syz = (dv/dz + dw/dy)
            Vmath::Vadd(nPts, &derivativesO1[1][2][0], 1,
                        &derivativesO1[2][1][0], 1, &Syz[0], 1);
            
            Vmath::Zero(nPts, &tmpar2[0], 1);
            
            // tmp = h_y*eps_bar
            Vmath::Vmul(nPts,
                        &h_av[1][0], 1,
                        &eps_bar[0], 1,
                        &tmpar2[0], 1);
            
            // tmp = h_y*eps_bar*d(rhow)/dy
            Vmath::Vmul(nPts,
                        &eps_bar[0], 1,
                        &derivativesMix[1][2][0], 1,
                        &tmpar2[0], 1);
            
            // Syz = mu * (du/dy + dv/dx)
            Vmath::Vmul(nPts, &mu[0], 1, &Syz[0], 1, &Syz[0], 1);
            // Syz = mu * (du/dy + dv/dx) + h_y*eps_bar*d(rhow)/dy
            Vmath::Vadd(nPts, &tmpar2[0], 1, &Sxy[0], 1, &Sxy[0], 1);
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
                Vmath::Zero(nPts, viscousTensor[k][i], 1);
            }
        }
        
        if (m_spacedim == 1)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);
            
            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][2][0], 1,
                        &tmp1[0], 1);
            
            // tmp = h_av*eps_bar
            Vmath::Vmul(nPts, &h_av[0][0], 1, &eps_bar[0], 1, &tmp2[0], 1);
            // tmp = h_av*eps_bar*d(rhoH)/dy
            Vmath::Vmul(nPts, &eps_bar[0], 1,
                        &derivativesMix[0][m_spacedim][0], 1,
                        &tmp2[0], 1);
            
            // STx = u * Sxx + (K / mu) * dT/dx +  + eps_bar*d(rhoH)/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
        }
        else if (m_spacedim == 2)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            Array<OneD, NekDouble > tmp3(nPts, 0.0);
            
            // Computation of STx
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sxy[0], 1, &tmp1[0], 1);
            
            // k * dT/dx
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][2][0], 1,
                        &tmp2[0], 1);
            
            // tmp = h_av*eps_bar
            Vmath::Vmul(nPts, &h_av[0][0], 1, &eps_bar[0], 1, &tmp3[0], 1);
            // tmp = h_av*eps_bar*d(rhoH)/dx
            Vmath::Vmul(nPts, &eps_bar[0], 1,
                        &derivativesMix[0][m_spacedim][0], 1,
                        &tmp3[0], 1);
            
            // STx = u * Sxx + v * Sxy + K * dT/dx + eps_bar*d(rhoH)/dx
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
            
            // k * dT/dy
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[1][2][0], 1,
                        &tmp2[0], 1);
            
            // tmp = h_av*eps_bar
            Vmath::Vmul(nPts, &h_av[1][0], 1, &eps_bar[0], 1, &tmp3[0], 1);
            // tmp = h_av*eps_bar*d(rhoH)/dy
            Vmath::Vmul(nPts,
                        &eps_bar[0], 1,
                        &derivativesMix[1][m_spacedim][0], 1,
                        &tmp3[0], 1);
            
            // STy = v * Syy + u * Sxy + K * dT/dy + eps_bar*d(rhoH)/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);

            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
            
            Vmath::Vadd(nPts, &STy[0], 1, &tmp3[0], 1, &STy[0], 1);
            
        }
        else if (m_spacedim == 3)
        {
            Array<OneD, NekDouble > tmp1(nPts, 0.0);
            Array<OneD, NekDouble > tmp2(nPts, 0.0);
            Array<OneD, NekDouble > tmp3(nPts, 0.0);
            Array<OneD, NekDouble > tmp4(nPts, 0.0);
            
            // Computation of STx
            
            // u * Sxx
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sgg[0][0], 1, &STx[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Sxy[0], 1, &tmp1[0], 1);
            
            // v * Sxy
            Vmath::Vmul(nPts, &physfield[2][0], 1, &Sxz[0], 1, &tmp2[0], 1);
            
            // k * dT/dx
            
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[0][2][0], 1,
                        &tmp3[0], 1);
            
            // tmp = h_av*eps_bar
            Vmath::Vmul(nPts, &h_av[0][0], 1, &eps_bar[0], 1, &tmp4[0], 1);
            // tmp = h_av*eps_bar*d(rhoH)/dx
            Vmath::Vmul(nPts, &eps_bar[0], 1, &derivativesMix[0][m_spacedim][0], 1,
                        &tmp4[0], 1);
            
            // STx = u * Sxx + v * Sxy + w * Sxz
            //     + (K / mu) * dT/dx + eps_bar*d(rhoH)/dx
            Vmath::Vadd(nPts, &STx[0], 1, &tmp1[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp2[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp3[0], 1, &STx[0], 1);
            Vmath::Vadd(nPts, &STx[0], 1, &tmp4[0], 1, &STx[0], 1);
            
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
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[1][2][0], 1,
                        &tmp3[0], 1);
            
            // tmp = h_av*eps_bar
            Vmath::Vmul(nPts, &h_av[1][0], 1, &eps_bar[0], 1, &tmp4[0], 1);
            // tmp = h_av*eps_bar*d(rhoH)/dy
            Vmath::Vmul(nPts, &eps_bar[0], 1, &derivativesMix[1][m_spacedim][0], 1,
                        &tmp4[0], 1);
            
            // STy = v * Syy + u * Sxy + w * Syz + K * dT/dy
            Vmath::Vadd(nPts, &STy[0], 1, &tmp1[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp2[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp3[0], 1, &STy[0], 1);
            Vmath::Vadd(nPts, &STy[0], 1, &tmp4[0], 1, &STy[0], 1);
            // Computation of STz
            
            // Re-initialise temporary arrays
            Vmath::Zero(nPts, &tmp1[0], 1);
            Vmath::Zero(nPts, &tmp2[0], 1);
            Vmath::Zero(nPts, &tmp3[0], 1);
            Vmath::Zero(nPts, &tmp4[0], 1);
            
            // w * Szz
            Vmath::Vmul(nPts, &physfield[2][0], 1, &Sgg[2][0], 1, &STz[0], 1);
            
            // u * Sxz
            Vmath::Vmul(nPts, &physfield[0][0], 1, &Sxz[0], 1, &tmp1[0], 1);
            
            // v * Syz
            Vmath::Vmul(nPts, &physfield[1][0], 1, &Syz[0], 1, &tmp2[0], 1);
            
            // k * dT/dz
            Vmath::Smul(nPts, m_thermalConductivity, &derivativesO1[2][2][0], 1,
                        &tmp3[0], 1);
            
            // tmp = h_av*eps_bar
            Vmath::Vmul(nPts, &h_av[2][0], 1, &eps_bar[0], 1, &tmp4[0], 1);
            // tmp = h_av*eps_bar*d(rhoH)/dy
            Vmath::Vmul(nPts, &eps_bar[0], 1, &derivativesMix[2][m_spacedim][0], 1,
                        &tmp4[0], 1);
            
            // STz = w * Szz + u * Sxz + v * Syz + K * dT/dz
            //     + (K / mu) * dT/dz + eps_bar*d(rhoH)/dz
            Vmath::Vadd(nPts, &STz[0], 1, &tmp1[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp2[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp3[0], 1, &STz[0], 1);
            Vmath::Vadd(nPts, &STz[0], 1, &tmp4[0], 1, &STz[0], 1);
        }
        
        switch (m_spacedim)
        {
            case 1:
            {
                Array<OneD, NekDouble> tmpar3(nPts, 0.0);
                // f_11v = f_rho = eps_bar*h_x drho/dx
                Vmath::Vmul(nPts, &h_av[0][0], 1, &eps_bar[0], 1, &tmpar3[0], 1);
                Vmath::Vmul(nPts,
                            &tmpar3[0], 1,
                            &derivativesO1[0][nvariables-2][0], 1,
                            &viscousTensor[0][0][0], 1);
                
                // f_21v = f_rhou
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                
                // f_31v = f_E
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][2][0], 1);
                
                // f_41v = f_Eps
                Vmath::Vcopy(nPts, &Seps[0][0], 1, &viscousTensor[0][3][0], 1);
                break;
            }
            case 2:
            {
                // f_11v = f_rho = eps_bar*h_x drho/dx
                Array<OneD, NekDouble> tmpar3(nPts, 0.0);
                
                Vmath::Vmul(nPts, &h_av[0][0], 1, &eps_bar[0], 1, &tmpar3[0], 1);
                
                Vmath::Vmul(nPts,
                            &tmpar3[0], 1,
                            &derivativesO1[0][nvariables-2][0], 1,
                            &viscousTensor[0][0][0], 1);
                
                // f_12v = f_rho2 = eps_bar*h_y drho/dy
                Vmath::Zero(nPts, &tmpar3[0], 1);
                
                Vmath::Vmul(nPts, &h_av[1][0], 1, &eps_bar[0], 1, &tmpar3[0], 1);
                Vmath::Vmul(nPts,
                            &tmpar3[0], 1,
                            &derivativesO1[1][nvariables-2][0], 1,
                            &viscousTensor[1][0][0], 1);
                
                // f_21v = f_rhou1
                Vmath::Vcopy(nPts, &Sgg[0][0], 1, &viscousTensor[0][1][0], 1);
                // f_22v = f_rhou2
                Vmath::Vcopy(nPts, &Sxy[0],    1, &viscousTensor[1][1][0], 1);
                
                // f_31v = f_rhov1
                Vmath::Vcopy(nPts, &Syx[0],    1, &viscousTensor[0][2][0], 1);
                // f_32v = f_rhov2
                Vmath::Vcopy(nPts, &Sgg[1][0], 1, &viscousTensor[1][2][0], 1);
                
                // f_41v = f_E1
                Vmath::Vcopy(nPts, &STx[0], 1, &viscousTensor[0][3][0], 1);
                // f_42v = f_E2
                Vmath::Vcopy(nPts, &STy[0], 1, &viscousTensor[1][3][0], 1);
                
                // f_51v = f_Eps1
                Vmath::Vcopy(nPts, &Seps[0][0], 1, &viscousTensor[0][4][0], 1);
                // f_52v = f_Eps2
                Vmath::Vcopy(nPts, &Seps[1][0], 1, &viscousTensor[1][4][0], 1);
                break;
            }
            case 3:
            {
                // f_11v = f_rho = eps_bar*h_x drho/dx
                Array<OneD, NekDouble> tmpar3(nPts, 0.0);
                
                Vmath::Vmul(nPts, &h_av[0][0], 1, &eps_bar[0], 1, &tmpar3[0], 1);
                
                Vmath::Vmul(nPts,
                            &tmpar3[0], 1,
                            &derivativesO1[0][nvariables-2][0], 1,
                            &viscousTensor[0][0][0], 1);
                
                // f_12v = f_rho2 = eps_bar*h_y drho/dy
                Vmath::Zero(nPts, &tmpar3[0], 1);
                
                Vmath::Vmul(nPts, &h_av[1][0], 1, &eps_bar[0], 1, &tmpar3[0], 1);
                Vmath::Vmul(nPts,
                            &tmpar3[0], 1,
                            &derivativesO1[1][nvariables-2][0], 1,
                            &viscousTensor[1][0][0], 1);
                
                // f_13v = f_rho3 = eps_bar*h_z drho/dz
                Vmath::Zero(nPts, &tmpar3[0], 1);
                
                Vmath::Vmul(nPts, &h_av[2][0], 1, &eps_bar[0], 1, &tmpar3[0], 1);
                Vmath::Vmul(nPts,
                            &tmpar3[0], 1,
                            &derivativesO1[2][nvariables-2][0], 1,
                            &viscousTensor[2][0][0], 1);
                
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
                
                // f_51v = f_Eps1
                Vmath::Vcopy(nPts, &Seps[0][0], 1, &viscousTensor[0][5][0], 1);
                // f_52v = f_Eps2
                Vmath::Vcopy(nPts, &Seps[1][0], 1, &viscousTensor[1][5][0], 1);
                // f_53v = f_Eps3
                Vmath::Vcopy(nPts, &Seps[2][0], 1, &viscousTensor[2][5][0], 1);
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

    void CompressibleFlowSystem::GetEnthalpy(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &pressure,
                  Array<OneD,                   NekDouble>   &enthalpy)
    {
        int       npts  = m_fields[0]->GetTotPoints();
        
        // Calculate p/rho;
        Vmath::Vdiv(npts, pressure, 1, physfield[0], 1, enthalpy, 1);
        
        // pressure <- E - 0.5*pressure
        Vmath::Vadd(npts, enthalpy, 1, physfield[m_spacedim+1], 1, enthalpy, 1);
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
        const int nBCEdgePts = m_fields[0]->GetTotPoints();

        Vmath::Vmul(nBCEdgePts, physfield[1], 1, physfield[1], 1, mach, 1);

        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, physfield[1+i], 1, physfield[1+i], 1,
                                     mach,           1, mach,           1);
        }

        Vmath::Vdiv(nBCEdgePts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(nBCEdgePts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(nBCEdgePts, mach, 1, soundspeed,   1, mach, 1);
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
            mu[i] = mu_star * pow(ratio, 1.50) *
                    (T_star + 110.0) / (temperature[i] + 110.0);
        }
    }

    /**
     * @brief Calcualte entropy.
     */
    void CompressibleFlowSystem::GetEntropy(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        const Array<OneD, const NekDouble>               &pressure,
        const Array<OneD, const NekDouble>               &temperature,
              Array<OneD,       NekDouble>               &entropy)
    {
        NekDouble entropyL2, entropy_sum = 0.0;
        const int npts = m_fields[0]->GetTotPoints();
        const NekDouble temp_inf = m_pInf/(m_rhoInf*m_gasConstant);;
        
        Array<OneD,       NekDouble> L2entropy(npts, 0.0);
        
        for (int i = 0; i < npts; ++i)
        {
            entropy[i] = m_gamma/(m_gamma-1.0)*m_gasConstant*log(temperature[i]/temp_inf) -
            m_gasConstant*log(pressure[i]/m_pInf);
            
        }
        
        Vmath::Vmul(npts,entropy,1,entropy,1,L2entropy,1);
        
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
        NekDouble minLength;
        NekDouble alpha   = MaxTimeStepEstimator();
        NekDouble cLambda = 0.2; // Spencer book-317

        // Loop over elements to compute the time-step limit for each element
        for(n = 0; n < nElements; ++n)
        {
            int npoints = m_fields[0]->GetExp(n)->GetTotPoints();
            Array<OneD, NekDouble> one2D(npoints, 1.0);
            NekDouble Area = m_fields[0]->GetExp(n)->Integral(one2D);

            if (boost::dynamic_pointer_cast<LocalRegions::TriExp>(
                    m_fields[0]->GetExp(n)))
            {
                minLength = 2.0 * sqrt(Area);
            }

            else if (boost::dynamic_pointer_cast<LocalRegions::QuadExp>(
                         m_fields[0]->GetExp(n)))
            {
                minLength = sqrt(Area);
            }
            else if (boost::dynamic_pointer_cast<LocalRegions::HexExp>(
                         m_fields[0]->GetExp(n)))
            {
                minLength = sqrt(Area);
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
        
        int i, e, nCoeffsElement, NumModesElement, NumModesCuttOff, NumModesCuttOff_Dir1, NumModesCuttOff_Dir2, nQuadPointsElement;
        NekDouble SensorNumerator, SensorDenominator;
        
        int nVariables      = m_fields.num_elements();
        int nTotQuadPoints  = GetTotPoints();
        int nElements       = m_fields[0]->GetExpSize();
        
        // Find solution (SolP) at p = P;
        // The input array (physarray) is the solution at p = P;
        
        Array<OneD,int> ExpOrderElement = GetNumExpModesPerExp();
        
        Array<OneD, NekDouble> SolP(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> SolPmOne(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> SolNorm(nTotQuadPoints,0.0);
        
        Vmath::Vcopy(nTotQuadPoints,physarray[0],1,SolP,1);
        
        // Apply filtering procedure in 2D to obtain the solution at p = P - 1;
        
        // Getting the sensor for 2D elements. For now, this is only implemented for Quadrilateral elements
        
        // filtering procedure for Quadrilateral elements described by Biotto page 94-95
        int CoeffsCount = 0;
        
        for (e = 0; e < nElements; e++)
        {
            NumModesElement         = ExpOrderElement[e];
            
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            int nCoeffsElement = m_fields[0]->GetExp(e)->GetNcoeffs();
            int numCutOff = NumModesElement - 1;
            
            // Set-up of the Orthogonal basis for a Quadrilateral element which is
            // needed to obtain thesolution at P =  p - 1;
            
            Array<OneD, NekDouble> SolPElementPhys(nQuadPointsElement,0.0);
            Array<OneD, NekDouble> SolPElementCoeffs(nCoeffsElement,0.0);
            
            Array<OneD, NekDouble> SolPmOneElementPhys(nQuadPointsElement,0.0);
            Array<OneD, NekDouble> SolPmOneElementCoeffs(nCoeffsElement,0.0);
            
            // create vector the save the solution points per element at P = p;
            
            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPElementPhys[i] = SolP[CoeffsCount+i];
            }
            
            m_fields[0]->GetExp(e)->FwdTrans(SolPElementPhys,
                                             SolPElementCoeffs);
            
            // ReduceOrderCoeffs reduces the polynomial order of the solution that
            // is represented by the coeffs given as an inarray. This is done by
            // projecting the higher order solution onto the orthogonal basis and
            // padding the higher order coefficients with zeros.
            
            m_fields[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff,
                                                      SolPElementCoeffs,
                                                      SolPmOneElementCoeffs);
            
            m_fields[0]->GetExp(e)->BwdTrans(SolPmOneElementCoeffs,
                                             SolPmOneElementPhys);
            
            for (int i = 0; i < nQuadPointsElement; i++)
            {
                SolPmOne[CoeffsCount+i] = SolPmOneElementPhys[i];
            }
            
            NekDouble SolPmeanNumerator = 0.0;
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
                Sensor[CoeffsCount+i] = sqrt(SolPmeanNumerator/nQuadPointsElement)
                                    /sqrt(SolPmeanDenumerator/nQuadPointsElement);
                
                Sensor[CoeffsCount+i] = log10(Sensor[CoeffsCount+i]);
            }
            CoeffsCount += nQuadPointsElement;
        }

        CoeffsCount = 0.0;
        
        for (e = 0; e < nElements; e++)
        {
            NumModesElement         = ExpOrderElement[e];
            NekDouble ThetaS        = m_mu0;
            NekDouble Phi0          = m_Skappa;
            NekDouble DeltaPhi      = m_Kappa;
            nQuadPointsElement      = m_fields[0]->GetExp(e)->GetTotPoints();
            
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
                    SensorKappa[CoeffsCount+i] = ThetaS/2*(1+sin(M_PI*
                                        (Sensor[CoeffsCount+i]-Phi0)
                                                    /(2*DeltaPhi)));
                }
            }
            
            CoeffsCount += nQuadPointsElement;
        }

    }
    
    void CompressibleFlowSystem::GetForcingTerm(
              const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > outarrayForcing)
    {
        const int nPts = m_fields[0]->GetTotPoints();
        const int nvariables = m_fields.num_elements();
        const int nElements = m_fields[0]->GetExpSize();
        
        NekDouble hxmin = 0.0;
        NekDouble hymin = 0.0;
        NekDouble hmin  = 0.0;
        
        int C1 = 3.0;
        
        Array<OneD,  NekDouble>  Sensor(nPts, 0.0);
        Array<OneD,  NekDouble>  SensorKappa(nPts, 0.0);
        Array <OneD, NekDouble > Lambda(nPts, 0.0);
        Array <OneD, NekDouble > Tau(nPts, 0.0);
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
        
        // Determine the spacial dimension approximation of the element
        Array <OneD, Array <OneD, NekDouble > > h_av(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            h_av[i] = Array <OneD, NekDouble > (nPts, 0.0);
        }
        Array <OneD, Array <OneD, NekDouble > > ElDim(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            ElDim[i] = Array <OneD, NekDouble > (nElements, 0.0);
        }
        //
        GetElementDimensions(ElDim, h_minmin);
        //
        int PointCount = 0.0;
        for (int e = 0; e < nElements; e++)
        {
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            
            for (int n = 0; n < nQuadPointsElement; n++)
            {
                h_av[0][n + PointCount] = ElDim[0][e];
                h_av[1][n + PointCount] = ElDim[1][e];
            }
            
            PointCount += nQuadPointsElement;
        }
        
        hxmin = Vmath::Vmin(ElDim[0].num_elements(), &ElDim[0][0], 1);
        hymin = Vmath::Vmin(ElDim[1].num_elements(), &ElDim[1][0], 1);
        hmin  = min(hxmin, hymin);
        PointCount = 0;
        for (int e = 0; e < nElements; e++)
        {
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            
            for (int n = 0; n < nQuadPointsElement; n++)
            {
                pOrder[n + PointCount] = pOrderElmt[e];
                
                Tau[n + PointCount] = hmin/(C1*pOrder[n + PointCount]*LambdaMax); // order 1.0e-06
                
                outarrayForcing[nvariables-1][n + PointCount] = 1/Tau[n + PointCount]*(hxmin/pOrder[n + PointCount]*LambdaMax*SensorKappa[n + PointCount]-inarray[nvariables-1][n + PointCount]);
            }
            PointCount += nQuadPointsElement;
        }
    }

    void CompressibleFlowSystem::GetElementDimensions(
                          Array<OneD,       Array<OneD, NekDouble> > &outarray,
                          Array<OneD,       NekDouble > &hmin)
    {
        
        const int nPts = m_fields[0]->GetTotPoints();
        const int nvariables = m_fields.num_elements();
        const int nElements = m_fields[0]->GetExpSize();
        
        SpatialDomains::QuadGeomSharedPtr   ElQuadGeom;
        SpatialDomains::TriGeomSharedPtr    ElTriGeom;
        SpatialDomains::HexGeomSharedPtr    ElHexGeom;
        SpatialDomains::TetGeomSharedPtr    ElTetGeom;
        SpatialDomains::PrismGeomSharedPtr  ElPrismGeom;
        
        NekDouble hx = 0.0;
        NekDouble hy = 0.0;
        NekDouble hz = 0.0;
        
        int PointCount = 0;
        
        for (int e = 0; e < nElements; e++)
        {
            NekDouble nedges = m_fields[0]->GetExp(e)->GetNedges();
            Array <OneD, NekDouble> L1(nedges, 0.0);
            
            for (int j = 0; j < nedges; ++j)
            {
                if (boost::dynamic_pointer_cast<LocalRegions::QuadExp>(
                                                m_fields[0]->GetExp(e)))
                {
                    ElQuadGeom = boost::dynamic_pointer_cast<
                    SpatialDomains::QuadGeom>(m_fields[0]->GetExp(e)->GetGeom());
                    
                    NekDouble x0 = 0.0;
                    NekDouble y0 = 0.0;
                    NekDouble z0 = 0.0;
                    
                    NekDouble x1 = 0.0;
                    NekDouble y1 = 0.0;
                    NekDouble z1 = 0.0;
                    
                    ElQuadGeom->GetEdge(j)->GetVertex(0)->GetCoords(x0,y0,z0);
                    ElQuadGeom->GetEdge(j)->GetVertex(1)->GetCoords(x1,y1,z1);
                    
                    L1[j] = sqrt(pow((x0-x1),2)+pow((y0-y1),2)+pow((z0-z1),2));
                }
                
                if (boost::dynamic_pointer_cast<LocalRegions::TriExp>(
                                                m_fields[0]->GetExp(e)))
                {
                    ElTriGeom = boost::dynamic_pointer_cast<
                    SpatialDomains::TriGeom>(m_fields[0]->GetExp(e)->GetGeom());
                    
                    NekDouble x0 = 0.0;
                    NekDouble y0 = 0.0;
                    NekDouble z0 = 0.0;
                    
                    NekDouble x1 = 0.0;
                    NekDouble y1 = 0.0;
                    NekDouble z1 = 0.0;
                    
                    ElTriGeom->GetEdge(j)->GetVertex(0)->GetCoords(x0,y0,z0);
                    ElTriGeom->GetEdge(j)->GetVertex(1)->GetCoords(x1,y1,z1);
                    
                    L1[j] = sqrt(pow((x0-x1),2)+pow((y0-y1),2)+pow((z0-z1),2));
                }
                
                if (boost::dynamic_pointer_cast<LocalRegions::HexExp>(
                                                m_fields[0]->GetExp(e)))
                {
                    ElHexGeom = boost::dynamic_pointer_cast<
                    SpatialDomains::HexGeom>(m_fields[0]->GetExp(e)->GetGeom());
                    
                    NekDouble x0 = 0.0;
                    NekDouble y0 = 0.0;
                    NekDouble z0 = 0.0;
                    
                    NekDouble x1 = 0.0;
                    NekDouble y1 = 0.0;
                    NekDouble z1 = 0.0;
                    
                    ElHexGeom->GetEdge(j)->GetVertex(0)->GetCoords(x0,y0,z0);
                    ElHexGeom->GetEdge(j)->GetVertex(1)->GetCoords(x1,y1,z1);
                }
            }
            if(boost::dynamic_pointer_cast<LocalRegions::QuadExp>(
                                                m_fields[0]->GetExp(e)));
            {
                hx = min(L1[0], L1[2]);
                hy = min(L1[1], L1[3]);
                    
                outarray[0][e] = hx;
                outarray[1][e] = hy;
            }
            if (boost::dynamic_pointer_cast<LocalRegions::TriExp>(
                                                m_fields[0]->GetExp(e)))
            {
                hx = Vmath::Vmin(nedges, &L1[0], 1);
                hy = Vmath::Vmin(nedges, &L1[0], 1);
                    
                outarray[0][e] = hx;
                outarray[1][e] = hy;
            }
            
            if(boost::dynamic_pointer_cast<LocalRegions::HexExp>(
                                                m_fields[0]->GetExp(e)));
            {
                hx = Vmath::Vmin(nedges, &L1[0], 1);
                hy = Vmath::Vmin(nedges, &L1[0], 1);
                hz = Vmath::Vmin(nedges, &L1[0], 1);
                
                outarray[0][e] = hx;
                outarray[1][e] = hy;
                outarray[2][e] = hy;
            }
            
        }
        
        if (m_spacedim == 2)
        {
            hmin[0] = Vmath::Vmin(outarray[0].num_elements(), outarray[0], 1);
            hmin[1] = Vmath::Vmin(outarray[1].num_elements(), outarray[1], 1);
        }
        
        if (m_spacedim == 3)
        {
            hmin[0] = Vmath::Vmin(outarray[0].num_elements(), outarray[0], 1);
            hmin[1] = Vmath::Vmin(outarray[1].num_elements(), outarray[1], 1);
            hmin[2] = Vmath::Vmin(outarray[2].num_elements(), outarray[2], 1);
        }
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
        
        int i, j, k;
        int nvariables = physfield.num_elements();
        int nPts       = m_fields[0]->GetTotPoints();
        
        
        Array<OneD, NekDouble > pressure            (nPts, 0.0);
        Array<OneD, NekDouble > temperature         (nPts, 0.0);
        Array <OneD, NekDouble > sensor             (nPts, 0.0);
        Array <OneD, NekDouble > SensorKappa        (nPts, 0.0);
        Array <OneD, NekDouble > absVelocity        (nPts, 0.0);
        Array <OneD, NekDouble > soundspeed         (nPts, 0.0);
        Array <OneD, NekDouble > Lambda             (nPts, 0.0);
        Array <OneD, NekDouble > mu_var             (nPts, 0.0);
        Array <OneD, NekDouble > h_minmin           (m_spacedim, 0.0);
        Vmath::Zero(nPts, eps_bar, 1);
        
        // Thermodynamic related quantities
        GetPressure(physfield, pressure);
        GetTemperature(physfield, pressure, temperature);
        GetSoundSpeed(physfield, pressure, soundspeed);
        GetAbsoluteVelocity(physfield, absVelocity);
        GetSensor(physfield, sensor, SensorKappa);
        
        // Determine the maximum wavespeed
        Vmath::Vadd(nPts, absVelocity, 1, soundspeed, 1, Lambda, 1);
        NekDouble LambdaMax = Vmath::Vmax(nPts, Lambda, 1);
        
        // Determine hbar = hx_i/h
        const int nElements  = m_fields[0]->GetExpSize();
        
        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();
        Array<OneD, NekDouble> pOrder (nPts, 0.0);
        
        
        Array <OneD, Array <OneD, NekDouble > > ElDim(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            ElDim[i] = Array <OneD, NekDouble > (nElements, 0.0);
        }
        
        GetElementDimensions(ElDim, h_minmin);
        
        NekDouble h_mean_sumx = 0.0;
        NekDouble h_mean_sumy = 0.0;
        
        NekDouble h_meanx     = 0.0;
        NekDouble h_meany     = 0.0;
        
        for (int i = 0; i < nElements; ++i)
        {
            h_mean_sumx += ElDim[0][i];
            h_mean_sumy += ElDim[1][i];
        }
        
        h_meanx =  h_mean_sumx/nElements;
        h_meany =  h_mean_sumy/nElements;
        
        NekDouble h_mean  = (h_meanx+h_meany)/2;
        
        NekDouble ThetaH = LambdaMax/2*h_mean;
        NekDouble ThetaL = 0.01*LambdaMax/2*h_mean;;
        
        NekDouble Phi0     = (ThetaH+ThetaL)/2;
        NekDouble DeltaPhi = ThetaH-Phi0;
        
        int PointCount = 0.0;
        
        Vmath::Zero(eps_bar.num_elements(), eps_bar, 1);
        
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
        
        Vmath::Smul(nPts, m_eps_max, eps_bar, 1, eps_bar, 1);
    }
    
    
    void CompressibleFlowSystem::GetArtificialDynamicViscosity(
                        const Array<OneD, Array<OneD, NekDouble> > &physfield,
                              Array<OneD,             NekDouble  > &mu_var)
    {
        const int npts       = m_fields[0]->GetTotPoints();
        const int nElements  = m_fields[0]->GetExpSize();
        const int ploc       = m_fields[0]->GetExpSize();
    
        int PointCount = 0;
        int nTotQuadPoints  = GetTotPoints();
        
        Array <OneD , NekDouble > S_e               (nTotQuadPoints, 0.0);
        Array <OneD , NekDouble > se                (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > Sensor             (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > SensorKappa        (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > absVelocity        (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > soundspeed         (nTotQuadPoints, 0.0);
        Array <OneD, NekDouble > pressure           (nTotQuadPoints, 0.0);
        
        GetAbsoluteVelocity   (physfield, absVelocity);
        GetPressure      (physfield, pressure);
        GetSoundSpeed    (physfield, pressure, soundspeed);
        GetSensor        (physfield, Sensor, SensorKappa);
        
        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();
        
        Array <OneD , NekDouble > Lambda(nTotQuadPoints, 1.0);
        Vmath::Vadd(nTotQuadPoints,absVelocity,1,soundspeed,1,Lambda,1);
        
        // Determining the maximum wave speed in the element
        NekDouble LambdaMax = Vmath::Vmax(nTotQuadPoints,Lambda,1);
        NekDouble StdVMax = Vmath::Vmax(nTotQuadPoints,absVelocity,1);
        
        for (int e = 0; e < nElements; e++)
        {
            // Threshold value specified in C. Biottos thesis.
            // Based on a 1D shock tube problem S_k = log10(1/p^4)
            // See  G.E. Barter and D.L. Darmofal. Shock Capturing with PDE-based
            // artificial diffusion for DGFEM: Part 1 Formulation, Journal od
            // Computational Physics 229 (2010) 1810-1827 for further reference
            
            //NekDouble S_Kappa = -4.0*log10(pOrderElmt[e]);
            //const double Kappa = 0.5;
            
            // Adjustable depending on the coarsness of the mesh. Might want to
            // move this variable into the session file
            
            int nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            Array <OneD, NekDouble> one2D(nQuadPointsElement, 1.0);
            NekDouble Area = m_fields[0]->GetExp(e)->Integral(one2D);
            
            for (int n = 0; n < nQuadPointsElement; n++)
            {
                NekDouble mu_0 = m_mu0;
                
                if (Sensor[n+PointCount] < (m_Skappa-m_Kappa))
                {
                    mu_var[n+PointCount] = 0;
                }
                else if(Sensor[n+PointCount] >= (m_Skappa-m_Kappa)
                        && Sensor[n+PointCount] <= (m_Skappa+m_Kappa))
                {
                    mu_var[n+PointCount] = mu_0*(0.5*(1+sin(
                                           M_PI*(Sensor[n+PointCount]-m_Skappa-m_Kappa)/(2*m_Kappa))));
                }
                else if(Sensor[n+PointCount] > (m_Skappa+m_Kappa))
                {
                    mu_var[n+PointCount] = mu_0;
                }
            }
            
            PointCount += nQuadPointsElement;
        }
    }

    void CompressibleFlowSystem::SetVarPOrderElmt(
                const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                      Array<OneD,                    NekDouble  > &PolyOrder)
    {
        int e, cnt;
        NekDouble s_0, s_ds, s_sm, s_fl;
        
        int nElements  = m_fields[0]->GetExpSize();
        int npts       = m_fields[0]->GetTotPoints();
        
        Array<OneD, NekDouble > Sensor           (npts, 0.0);
        Array<OneD, NekDouble > SensorKappa      (npts, 0.0);
        Array<OneD, NekDouble > se               (npts,0.0);
        
        GetSensor(physfield, Sensor, SensorKappa);
        
        
        int numfields = m_fields.num_elements();
        
        
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
            // Ideally, these threshold values could be given in the Session File
            
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
                else if(se[npCount + i] > s_sm && se[npCount + i] < s_ds)
                {
                    if (PolyOrder[npCount + i] < MaxOrder)
                    {
                        PolyOrder[npCount + i] = PolyOrder[npCount + i] + 2;
                    }
                }
                else if(se[npCount + i] > s_fl && se[npCount + i] < s_sm)
                {
                    PolyOrder[npCount + i] = PolyOrder[npCount + i] + 1;
                }
                else if(se[npCount + i] < s_fl)
                {
                    if (PolyOrder[npCount + i] > MinOrder)
                    {
                            PolyOrder[npCount + i] = PolyOrder[npCount + i];
                    }
                }
            }
            m_file2 << "<E COMPOSITE= \"C[" << e+1 << "]\" NUMMODES=\"" << PolyOrder[npCount + 1] << "\" TYPE=\"MODIFIED\" FIELDS=\"rho,rhou,rhov,rhow,E\" />" << endl;
            npCount += nQuadPointsElement;
        }
        
        m_file2.close();
    }
}

