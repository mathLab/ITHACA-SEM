///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFE.cpp
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
// Description: Navier Stokes equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/NavierStokesCFE.h>

using namespace std;

namespace Nektar
{
    string NavierStokesCFE::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "NavierStokesCFE", NavierStokesCFE::create,
            "NavierStokes equations in conservative variables.");

    NavierStokesCFE::NavierStokesCFE(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          CompressibleFlowSystem(pSession, pGraph)
    {
    }

    NavierStokesCFE::~NavierStokesCFE()
    {

    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void NavierStokesCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        // rest of initialisation is in this routine so it can also be called
        // in NavierStokesImplicitCFE initialisation
        InitObject_Explicit(); 
    }

    void NavierStokesCFE::InitObject_Explicit()
    {
        // Get gas constant from session file and compute Cp
        NekDouble gasConstant;
        m_session->LoadParameter ("GasConstant",   gasConstant,   287.058);
        m_Cp      = m_gamma / (m_gamma - 1.0) * gasConstant;
        m_Cv      = m_Cp / m_gamma;

        m_session->LoadParameter ("Twall", m_Twall, 300.15);

        // Viscosity
        int nPts = m_fields[0]->GetNpoints();
        m_session->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
        m_session->LoadParameter ("mu",            m_muRef,           1.78e-05);
        m_mu = Array<OneD, NekDouble>(nPts, m_muRef);
        if (m_ViscosityType == "Variable")
        {
            m_is_mu_variable = true;
        }

        // Thermal conductivity or Prandtl
        if ( m_session->DefinesParameter("thermalConductivity"))
        {
            ASSERTL0( !m_session->DefinesParameter("Pr"),
                 "Cannot define both Pr and thermalConductivity.");

            m_session->LoadParameter ("thermalConductivity",
                                        m_thermalConductivityRef);
            m_Prandtl = m_Cp * m_muRef / m_thermalConductivityRef;
        }
        else
        {
            m_session->LoadParameter ("Pr", m_Prandtl, 0.72);
            m_thermalConductivityRef = m_Cp * m_muRef / m_Prandtl;
        }
        m_thermalConductivity =
                Array<OneD, NekDouble>(nPts, m_thermalConductivityRef);

        // Artificial viscosity parameter
        m_session->LoadParameter("mu0", m_mu0, 1.0);

        // load smoothing tipe
        m_session->LoadSolverInfo("Smoothing", m_smoothing, "Off");
        if (m_smoothing == "C0")
        {
            m_C0ProjectExp = MemoryManager<MultiRegions::ContField>::
            AllocateSharedPtr(m_session,m_graph,m_session->GetVariable(0));
        }
        // load physical sensor type
        m_session->LoadSolverInfo("PhysicalSensorType", m_physicalSensorType,
            "Off");

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

        m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance(diffName, diffName);
        if ("InteriorPenalty" == diffName)
        {
            m_is_diffIP = true;
            SetBoundaryConditionsBwdWeight();
        }

        if ("LDGNS" == diffName||
            "LDGNS3DHomogeneous1D" == diffName)
        {
            m_diffusion->SetFluxPenaltyNS(&NavierStokesCFE::
                v_GetFluxPenalty, this);
        }

        if (m_specHP_dealiasing)
        {
            m_diffusion->SetFluxVectorNS(
                &NavierStokesCFE::v_GetViscousFluxVectorDeAlias,
                this);
        }
        else
        {
            m_diffusion->SetFluxVectorNS(&NavierStokesCFE::
                                          v_GetViscousFluxVector, this);
        }

        m_diffusion->SetDiffusionFluxCons(
            &NavierStokesCFE::GetViscousFluxVectorConservVar<false>, this);

        m_diffusion->SetDiffusionFluxConsTrace(
            &NavierStokesCFE::GetViscousFluxVectorConservVar<true>, this);

        m_diffusion->SetSpecialBndTreat(
            &NavierStokesCFE::SpecialBndTreat, this);

        m_diffusion->SetDiffusionSymmFluxCons(
            &NavierStokesCFE::GetViscousSymmtrFluxConservVar, this);

        if (m_shockCaptureType != "Off")
        {
            m_diffusion->SetArtificialDiffusionVector(
                &NavierStokesCFE::GetArtificialViscosity, this);
        }

        m_diffusion->SetCalcViscosity(
                &NavierStokesCFE::CalcViscosity, this);
        
        // Concluding initialisation of diffusion operator
        m_diffusion->InitObject         (m_session, m_fields);

    }

    void NavierStokesCFE::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        bool extraFields;
        m_session->MatchSolverInfo("OutputExtraFields","True",
                                   extraFields, true);
        if (extraFields)
        {
            const int nPhys   = m_fields[0]->GetNpoints();
            const int nCoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.size());

            for (int i = 0; i < m_fields.size(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }

            Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);
            Array<OneD, Array<OneD, NekDouble> > velFwd  (m_spacedim);
            for (int i = 0; i < m_spacedim; ++i)
            {
                velocity[i] = Array<OneD, NekDouble> (nPhys);
                velFwd[i]   = Array<OneD, NekDouble> (nCoeffs);
            }

            Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
            Array<OneD, NekDouble> entropy(nPhys);
            Array<OneD, NekDouble> soundspeed(nPhys), mach(nPhys);
            Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

            m_varConv->GetVelocityVector(tmp, velocity);
            m_varConv->GetPressure  (tmp, pressure);
            m_varConv->GetTemperature(tmp, temperature);
            m_varConv->GetEntropy   (tmp, entropy);
            m_varConv->GetSoundSpeed(tmp, soundspeed);
            m_varConv->GetMach      (tmp, soundspeed, mach);

            int sensorOffset;
            m_session->LoadParameter ("SensorOffset", sensorOffset, 1);
            m_varConv->GetSensor (m_fields[0], tmp, sensor, SensorKappa,
                                    sensorOffset);

            Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);
            Array<OneD, NekDouble> sFwd(nCoeffs);
            Array<OneD, NekDouble> aFwd(nCoeffs), mFwd(nCoeffs);
            Array<OneD, NekDouble> sensFwd(nCoeffs);

            string velNames[3] = {"u", "v", "w"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                m_fields[0]->FwdTrans_IterPerExp(velocity[i], velFwd[i]);
                variables.push_back(velNames[i]);
                fieldcoeffs.push_back(velFwd[i]);
            }

            m_fields[0]->FwdTrans_IterPerExp(pressure,   pFwd);
            m_fields[0]->FwdTrans_IterPerExp(temperature,TFwd);
            m_fields[0]->FwdTrans_IterPerExp(entropy,    sFwd);
            m_fields[0]->FwdTrans_IterPerExp(soundspeed, aFwd);
            m_fields[0]->FwdTrans_IterPerExp(mach,       mFwd);
            m_fields[0]->FwdTrans_IterPerExp(sensor,     sensFwd);

            variables.push_back  ("p");
            variables.push_back  ("T");
            variables.push_back  ("s");
            variables.push_back  ("a");
            variables.push_back  ("Mach");
            variables.push_back  ("Sensor");
            fieldcoeffs.push_back(pFwd);
            fieldcoeffs.push_back(TFwd);
            fieldcoeffs.push_back(sFwd);
            fieldcoeffs.push_back(aFwd);
            fieldcoeffs.push_back(mFwd);
            fieldcoeffs.push_back(sensFwd);

            if (m_artificialDiffusion)
            {
                // Get min h/p
                // m_artificialDiffusion->m_hOverP = GetElmtMinHP();
                // reuse pressure
                Array<OneD, NekDouble> sensorFwd(nCoeffs);
                m_artificialDiffusion->GetArtificialViscosity(tmp, pressure);
                m_fields[0]->FwdTrans_IterPerExp(pressure,   sensorFwd);

                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(sensorFwd);

            }

            if (m_shockCaptureType == "Physical")
            {
                // GetPhysicalAV(tmp);
                Array<OneD, NekDouble> muavFwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(m_muav,   muavFwd);
                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(muavFwd);

                // Debug Ducros
                // div square
                Array<OneD, NekDouble> dv2Fwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(m_diffusion->m_divVelSquare,
                    dv2Fwd);
                variables.push_back  ("divVelSquare");
                fieldcoeffs.push_back(dv2Fwd);
                // curl square
                Array<OneD, NekDouble> cv2Fwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(m_diffusion->m_curlVelSquare,
                    cv2Fwd);
                variables.push_back  ("curlVelSquare");
                fieldcoeffs.push_back(cv2Fwd);
                // Ducros
                Array<OneD, NekDouble> duc(nPhys,1.0);
                Ducros(duc);
                Array<OneD, NekDouble> ducFwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(duc, ducFwd);
                variables.push_back  ("Ducros");
                fieldcoeffs.push_back(ducFwd);

            }
        }
    }

    void NavierStokesCFE::v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              Array<OneD,       Array<OneD, NekDouble>> &outarray,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
    {
        size_t nvariables = inarray.size();
        size_t npoints    = GetNpoints();
        size_t nTracePts  = GetTraceTotPoints();

        // this should be preallocated
        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }

        // get artificial viscosity
        if (m_shockCaptureType == "Physical" && m_CalcPhysicalAV)
        {
            GetPhysicalAV(inarray);
        }

        if (m_is_diffIP)
        {
            if (m_bndEvaluateTime < 0.0)
            {
                NEKERROR(ErrorUtil::efatal, "m_bndEvaluateTime not setup");
            }
            m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff,
                                m_bndEvaluateTime,
                                pFwd, pBwd);
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
        }
        else
        {
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inFwd(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inBwd(nvariables-1);


            for (int i = 0; i < nvariables - 1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>{npoints};
                inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
                inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
            }

            // Extract pressure
            // (use inarrayDiff[0] as a temporary storage for the pressure)
            m_varConv->GetPressure(inarray, inarrayDiff[0]);

            // Extract temperature
            m_varConv->GetTemperature(inarray, inarrayDiff[nvariables - 2]);

            // Extract velocities
            m_varConv->GetVelocityVector(inarray, inarrayDiff);

            // Repeat calculation for trace space
            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                inFwd = NullNekDoubleArrayofArray;
                inBwd = NullNekDoubleArrayofArray;
            }
            else
            {
                m_varConv->GetPressure(pFwd, inFwd[0]);
                m_varConv->GetPressure(pBwd, inBwd[0]);

                m_varConv->GetTemperature(pFwd, inFwd[nvariables - 2]);
                m_varConv->GetTemperature(pBwd, inBwd[nvariables - 2]);

                m_varConv->GetVelocityVector(pFwd, inFwd);
                m_varConv->GetVelocityVector(pBwd, inBwd);
            }

            // Diffusion term in physical rhs form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff,
                                outarrayDiff, inFwd, inBwd);

            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
        }
    }



    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::v_GetViscousFluxVector(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
              TensorOfArray3D<NekDouble>                &derivativesO1,
              TensorOfArray3D<NekDouble>                &viscousTensor)
    {
        // Auxiliary variables
        size_t nScalar    = physfield.size();
        size_t nPts       = physfield[0].size();
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Update viscosity and thermal conductivity
        GetViscosityAndThermalCondFromTemp(physfield[nScalar-1], m_mu,
            m_thermalConductivity);

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            // Apply Ducros sensor
            if (m_physicalSensorType == "Ducros" && m_CalcPhysicalAV)
            {
                Ducros(m_muav);
            }
            // Apply approximate c0 smoothing
            if (m_smoothing == "C0" && m_CalcPhysicalAV)
            {
                C0Smooth(m_muav);
            }
            Vmath::Vadd(nPts, m_mu, 1, m_muav, 1, m_mu, 1);
            // Freeze AV for Implicit time stepping
            if (m_explicitDiffusion == false)
            {
                m_CalcPhysicalAV = false;
            }
        }

        // Velocity divergence
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, derivativesO1[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, m_mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, derivativesO1[i][j], 1,
                                  derivativesO1[j][i], 1,
                                  viscousTensor[i][j+1], 1);

                Vmath::Vmul(nPts, m_mu, 1,
                                  viscousTensor[i][j+1], 1,
                                  viscousTensor[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, viscousTensor[i][j+1], 1,
                                  divVel, 1,
                                  viscousTensor[i][j+1], 1);
                }
                else
                {
                    // Copy to make symmetric
                    Vmath::Vcopy(nPts, viscousTensor[i][j+1], 1,
                                       viscousTensor[j][i+1], 1);
                }
            }
        }

        // Terms for the energy equation
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (int j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, physfield[j], 1,
                               viscousTensor[i][j+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, m_thermalConductivity, 1,
                               derivativesO1[i][m_spacedim], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::v_GetViscousFluxVectorDeAlias(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
              TensorOfArray3D<NekDouble>                &derivativesO1,
              TensorOfArray3D<NekDouble>                &viscousTensor)
    {
        // Factor to rescale 1d points in dealiasing.
        NekDouble OneDptscale = 2;
        // Get number of points to dealias a cubic non-linearity
        size_t nScalar   = physfield.size();
        int nPts      = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        size_t nPts_orig = physfield[0].size();

        // Auxiliary variables
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Update viscosity and thermal conductivity
        GetViscosityAndThermalCondFromTemp(physfield[nScalar-1], m_mu,
            m_thermalConductivity);

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            // Apply Ducros sensor
            if (m_physicalSensorType == "Ducros" && m_CalcPhysicalAV)
            {
                Ducros(m_muav);
            }
            // Apply approximate c0 smoothing
            if (m_smoothing == "C0" && m_CalcPhysicalAV)
            {
                C0Smooth(m_muav);
            }
            Vmath::Vadd(nPts, m_mu, 1, m_muav, 1, m_mu, 1);
            // Freeze AV for Implicit time stepping
            if (m_explicitDiffusion == false)
            {
                m_CalcPhysicalAV = false;
            }
        }

        // Interpolate inputs and initialise interpolated output
        Array<OneD, Array<OneD, NekDouble> > vel_interp(m_spacedim);
        TensorOfArray3D<NekDouble>           deriv_interp(m_spacedim);
        TensorOfArray3D<NekDouble>           out_interp(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            // Interpolate velocity
            vel_interp[i]   = Array<OneD, NekDouble> (nPts);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], vel_interp[i]);

            // Interpolate derivatives
            deriv_interp[i] = Array<OneD, Array<OneD, NekDouble> >
                             (m_spacedim+1);
            for (int j = 0; j < m_spacedim+1; ++j)
            {
                deriv_interp[i][j] = Array<OneD, NekDouble> (nPts);
                m_fields[0]->PhysInterp1DScaled(
                    OneDptscale, derivativesO1[i][j], deriv_interp[i][j]);
            }

            // Output (start from j=1 since flux is zero for rho)
            out_interp[i] = Array<OneD, Array<OneD, NekDouble> > (m_spacedim+2);
            for (int j = 1; j < m_spacedim+2; ++j)
            {
                out_interp[i][j] = Array<OneD, NekDouble> (nPts);
            }
        }

        // Velocity divergence
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, deriv_interp[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, m_mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0 (no need to dealias)
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts_orig, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, deriv_interp[i][j], 1,
                                  deriv_interp[j][i], 1,
                                  out_interp[i][j+1], 1);

                Vmath::Vmul(nPts, m_mu, 1,
                                  out_interp[i][j+1], 1,
                                  out_interp[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, out_interp[i][j+1], 1,
                                  divVel, 1,
                                  out_interp[i][j+1], 1);
                }
                else
                {
                    // Make symmetric
                    out_interp[j][i+1] = out_interp[i][j+1];
                }
            }
        }

        // Terms for the energy equation
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, out_interp[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (int j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, vel_interp[j], 1,
                               out_interp[i][j+1], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, m_thermalConductivity, 1,
                               deriv_interp[i][m_spacedim], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
        }

        // Project to original space
        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = 1; j < m_spacedim+2; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale,
                    out_interp[i][j],
                    viscousTensor[i][j]);
            }
        }
    }


    /**
     * @brief For very special treatment. For general boundaries it does nothing
     * But for WallViscous and WallAdiabatic, special treatment is needed
     * because they get the same Bwd value, special treatment is needed for
     * boundary treatment of diffusion flux
     * Note: This special treatment could be removed by seperating
     * WallViscous and WallAdiabatic into two different classes.
     *
     */
    void NavierStokesCFE::SpecialBndTreat(
        Array<OneD, Array<OneD, NekDouble>> &consvar)
    {
        size_t nConvectiveFields = consvar.size();
        int ndens       = 0;
        int nengy       = nConvectiveFields - 1;

        Array<OneD, Array<OneD, NekDouble>> bndCons {nConvectiveFields};
        Array<OneD, NekDouble> bndTotEngy;
        Array<OneD, NekDouble> bndPressure;
        Array<OneD, NekDouble> bndRho;
        Array<OneD, NekDouble> bndIntEndy;
        int nLengthArray = 0;

        int cnt = 0;
        int nBndRegions = m_fields[nengy]->
            GetBndCondExpansions().size();
        for (int j = 0; j < nBndRegions; ++j)
        {
            if (m_fields[nengy]->GetBndConditions()[j]->
                GetBoundaryConditionType() ==
                SpatialDomains::ePeriodic)
            {
                continue;
            }

            size_t nBndEdges = m_fields[nengy]->
            GetBndCondExpansions()[j]->GetExpSize();
            for (int e = 0; e < nBndEdges; ++e)
            {
                size_t nBndEdgePts = m_fields[nengy]->
                GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                int id2 = m_fields[0]->GetTrace()->
                GetPhys_Offset(m_fields[0]->GetTraceMap()->
                            GetBndCondIDToGlobalTraceID(cnt++));

                // Imposing Temperature Twall at the wall
                if (boost::iequals(m_fields[nengy]->GetBndConditions()[j]->
                    GetUserDefined(), "WallViscous"))
                {
                    if (nBndEdgePts != nLengthArray)
                    {
                        for (int i = 0; i < nConvectiveFields; ++i)
                        {
                            bndCons[i] = Array<OneD, NekDouble>
                                        {nBndEdgePts, 0.0};
                        }
                        bndTotEngy  = Array<OneD, NekDouble> {nBndEdgePts, 0.0};
                        bndPressure = Array<OneD, NekDouble> {nBndEdgePts, 0.0};
                        bndRho      = Array<OneD, NekDouble> {nBndEdgePts, 0.0};
                        bndIntEndy  = Array<OneD, NekDouble> {nBndEdgePts, 0.0};
                        nLengthArray = nBndEdgePts;
                    }
                    else
                    {
                        Vmath::Zero(nLengthArray, bndPressure, 1);
                        Vmath::Zero(nLengthArray, bndRho     , 1);
                        Vmath::Zero(nLengthArray, bndIntEndy , 1);
                    }

                    Array<OneD, NekDouble> tmp;

                    for (int k = 0; k < nConvectiveFields; ++k)
                    {
                        Vmath::Vcopy(nBndEdgePts, tmp = consvar[k] + id2, 1,
                                    bndCons[k], 1);
                    }

                    m_varConv->GetPressure(bndCons, bndPressure);
                    Vmath::Fill(nLengthArray, m_Twall, bndTotEngy, 1);
                    m_varConv->GetRhoFromPT(bndPressure, bndTotEngy, bndRho);
                    m_varConv->GetEFromRhoP(bndRho, bndPressure, bndIntEndy);
                    m_varConv->GetDynamicEnergy(bndCons, bndTotEngy);

                    Vmath::Vvtvp(nBndEdgePts, bndIntEndy, 1, bndCons[ndens], 1,
                        bndTotEngy, 1, bndTotEngy, 1);

                    Vmath::Vcopy(nBndEdgePts,
                                bndTotEngy, 1,
                                tmp = consvar[nengy] + id2, 1);
                }
            }
        }
    }

    /**
     * @brief Calculate and return the ArtificialViscosity for shock-capturing.
     */
    void NavierStokesCFE::GetArtificialViscosity(
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, NekDouble>                    &muav)
    {
        m_artificialDiffusion->GetArtificialViscosity(inarray, muav);
    }

    /**
     * @brief Calculate and return the Symmetric flux in IP method.
     */
    void NavierStokesCFE::GetViscousSymmtrFluxConservVar(
        const int                                           nDim,
        const Array<OneD, Array<OneD, NekDouble> >          &inaverg,
        const Array<OneD, Array<OneD, NekDouble > >         &inarray,
        TensorOfArray3D<NekDouble>                          &outarray,
        Array< OneD, int >                                  &nonZeroIndex,
        const Array<OneD, Array<OneD, NekDouble> >          &normal)
    {
        size_t nConvectiveFields   = inarray.size();
        size_t nPts                = inaverg[nConvectiveFields - 1].size();
        nonZeroIndex = Array<OneD, int>{nConvectiveFields - 1, 0};
        for (int i = 0; i < nConvectiveFields - 1; ++i)
        {
            nonZeroIndex[i] =   i + 1;
        }

        std::vector<NekDouble> inAvgTmp(nConvectiveFields);
        std::vector<NekDouble> inTmp(nConvectiveFields);
        std::vector<NekDouble> outTmp(nConvectiveFields);
        for (int d = 0; d < nDim; ++d)
        {
            for (int nderiv = 0; nderiv < nDim; ++nderiv)
            {
                for (size_t p = 0; p < nPts; ++p)
                {
                    // rearrenge data
                    for (int f = 0; f < nConvectiveFields; ++f)
                    {
                        inAvgTmp[f] = inaverg[f][p];
                        inTmp[f] = inarray[f][p];
                    }
                    
                    // get temp
                    NekDouble temperature = m_varConv->GetTemperature(inTmp.data());
                    // get viscosity
                    NekDouble mu;
                    GetViscosityFromTempKernel(temperature, mu);

                    GetViscousFluxBilinearFormKernel(nDim, d, nderiv,
                        inAvgTmp.data(), inTmp.data(), mu, outTmp.data());

                    for (int f = 0; f < nConvectiveFields; ++f)
                    {
                        outarray[d][f][p] += normal[d][p] * outTmp[f];
                    }
                }
            }
        }
    }
    
     /**
    * @brief Calculate the physical artificial viscosity
    *
    * @param physfield  Input field.
    */
    void NavierStokesCFE::GetPhysicalAV(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield)
    {
        int nPts = physfield[0].size();
        int nElements = m_fields[0]->GetExpSize();
        Array <OneD, NekDouble > hOverP(nElements, 0.0);
        hOverP = GetElmtMinHP();

        // Determine the maximum wavespeed
        Array <OneD, NekDouble > Lambdas(nPts, 0.0);
        Array <OneD, NekDouble > soundspeed(nPts, 0.0);
        Array <OneD, NekDouble > absVelocity(nPts, 0.0);
        m_varConv->GetSoundSpeed(physfield, soundspeed);
        m_varConv->GetAbsoluteVelocity(physfield, absVelocity);

        Vmath::Vadd(nPts, absVelocity, 1, soundspeed, 1, Lambdas, 1);

        // Compute sensor based on rho
        Array<OneD, NekDouble> Sensor(nPts, 0.0);
        m_varConv->GetSensor(m_fields[0], physfield, Sensor, m_muav, 1);

        Array<OneD, NekDouble> tmp;
        for (int e = 0; e < nElements; e++)
        {
            int physOffset      = m_fields[0]->GetPhys_Offset(e);
            int nElmtPoints     = m_fields[0]->GetExp(e)->GetTotPoints();

            // Compute the maximum wave speed
            NekDouble LambdaElmt = Vmath::Vmax(nElmtPoints, tmp = Lambdas
                + physOffset, 1);

            // Compute average bounded density
            NekDouble rhoAve = Vmath::Vsum(nElmtPoints, tmp = physfield[0]
                + physOffset, 1);
            rhoAve = rhoAve / nElmtPoints;
            rhoAve = Smath::Smax(rhoAve , 1.0e-4, 1.0e+4);

            // Scale sensor by coeff, h/p, and density
            LambdaElmt *= m_mu0 * hOverP[e] * rhoAve;
            Vmath::Smul(nElmtPoints, LambdaElmt, tmp = m_muav + physOffset, 1,
                tmp = m_muav + physOffset, 1);
        }
    }
        
    void NavierStokesCFE::CalcViscosity(
        const Array<OneD, const Array<OneD, NekDouble>> &inaverg,
              Array<OneD, NekDouble>                    &mu)
    {
        int nConvectiveFields = inaverg.size();
        int nPts = inaverg[nConvectiveFields-1].size();

        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble> tmp(nPts,0.0);
            m_varConv->GetTemperature(inaverg,tmp);
            m_varConv->GetDynamicViscosity(tmp, mu);
        }
        else
        {
            //mu may be on volume or trace 
            Vmath::Fill(nPts, m_mu[0], mu, 1);
        }
    }
    /**
     * @brief Make field C0.
     *
     * @param field Input Field
     */
    void NavierStokesCFE::C0Smooth( Array<OneD, NekDouble> &field )
    {
        int nCoeffs = m_C0ProjectExp->GetNcoeffs();
        Array<OneD, NekDouble> muFwd(nCoeffs);
        Array<OneD, NekDouble> weights(nCoeffs, 1.0);
        // Assemble global expansion coefficients for viscosity
        m_C0ProjectExp->FwdTrans_IterPerExp(field,
            m_C0ProjectExp->UpdateCoeffs());
        m_C0ProjectExp->Assemble();
        Vmath::Vcopy(nCoeffs, m_C0ProjectExp->GetCoeffs(), 1, muFwd, 1);
        // Global coefficients
        Vmath::Vcopy(nCoeffs, weights, 1,
            m_C0ProjectExp->UpdateCoeffs(), 1);
        // This is the sign vector
        m_C0ProjectExp->GlobalToLocal();
        // Get weights
        m_C0ProjectExp->Assemble();
        // Divide
        Vmath::Vdiv(nCoeffs, muFwd, 1, m_C0ProjectExp->GetCoeffs(), 1,
            m_C0ProjectExp->UpdateCoeffs(), 1);
        // Get local coefficients
        m_C0ProjectExp->GlobalToLocal();
        // Get C0 field
        m_C0ProjectExp->BwdTrans_IterPerExp(
        m_C0ProjectExp->GetCoeffs(), field);
    }

    /**
    * @brief Applied Ducros (anti-vorticity) sensor.
    *
    * @param field Input Field
    */
    void NavierStokesCFE::Ducros( Array<OneD, NekDouble> &field )
    {
        int nPts = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> denDuc(nPts, NekConstants::kNekZeroTol);
        Array<OneD, NekDouble> ducros(nPts, 0.0);
        Vmath::Vadd(nPts, denDuc, 1, m_diffusion->m_divVelSquare, 1,
            denDuc, 1);
        Vmath::Vadd(nPts, denDuc, 1, m_diffusion->m_curlVelSquare, 1,
            denDuc, 1);
        Vmath::Vdiv(nPts, m_diffusion->m_divVelSquare, 1, denDuc, 1,
            ducros, 1);
        // Average in cell
        Array<OneD, NekDouble> tmp;
        for (int e = 0; e < m_fields[0]->GetExpSize(); e++)
        {
            int nElmtPoints     = m_fields[0]->GetExp(e)->GetTotPoints();
            int physOffset      = m_fields[0]->GetPhys_Offset(e);

            NekDouble eAve = Vmath::Vsum(nElmtPoints, 
                tmp = ducros + physOffset, 1);
            eAve = eAve / nElmtPoints;
            Vmath::Fill(nElmtPoints, eAve, tmp = ducros + physOffset, 1);
        }
        Vmath::Vmul(nPts, ducros, 1, field, 1, field, 1);
    }


    /**
     * @brief Return the penalty vector for the LDGNS diffusion problem.
     */
    void NavierStokesCFE::v_GetFluxPenalty(
        const Array<OneD, const Array<OneD, NekDouble>> &uFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &uBwd,
              Array<OneD,       Array<OneD, NekDouble>> &penaltyCoeff)
    {
        unsigned int nTracePts  = uFwd[0].size();

        // Compute average temperature
        unsigned int nVariables = uFwd.size();
        Array<OneD, NekDouble> tAve{nTracePts, 0.0};
        Vmath::Svtsvtp(nTracePts, 0.5, uFwd[nVariables-1], 1,
            0.5, uBwd[nVariables-1], 1, tAve, 1);

        // Get average viscosity and thermal conductivity
        Array<OneD, NekDouble> muAve{nTracePts, 0.0};
        Array<OneD, NekDouble> tcAve{nTracePts, 0.0};

        GetViscosityAndThermalCondFromTemp(tAve, muAve, tcAve);

        // Compute penalty term
        for (int i = 0; i < nVariables; ++i)
        {
            // Get jump of u variables
            Vmath::Vsub(nTracePts, uFwd[i], 1, uBwd[i], 1, penaltyCoeff[i], 1);
            // Multiply by variable coefficient = {coeff} ( u^+ - u^- )
            if ( i < nVariables-1 )
            {
                Vmath::Vmul(nTracePts, muAve, 1, penaltyCoeff[i], 1,
                    penaltyCoeff[i], 1);
            }
            else
            {
                Vmath::Vmul(nTracePts, tcAve, 1, penaltyCoeff[i], 1,
                    penaltyCoeff[i], 1);
            }
        }
    }

    
    /**
     * @brief Update viscosity
     * todo: add artificial viscosity here
     */
    void NavierStokesCFE::GetViscosityAndThermalCondFromTemp(
        const Array<OneD, NekDouble> &temperature,
              Array<OneD, NekDouble> &mu,
              Array<OneD, NekDouble> &thermalCond)
    {
        int nPts = temperature.size();

        for (size_t p = 0; p < nPts; ++p)
        {
            GetViscosityAndThermalCondFromTempKernel(temperature[p], mu[p],
                thermalCond[p]);
        }
    }

    /**
    * @brief Get trace of the physical artificial viscosity
    *
    */
    void NavierStokesCFE::GetTracePhysicalAV()
    {
        int nTracePts = m_fields[0]->GetTrace()->GetTotPoints();
        Array<OneD, NekDouble> Fwd(nTracePts,0.0);
        Array<OneD, NekDouble> Bwd(nTracePts,0.0);
        // BwdMuvar is left to be 0.0 according to DiffusionLDG.cpp
        m_fields[0]->GetFwdBwdTracePhys(m_muav,Fwd,Bwd, false);

        for(int k = 0; k < nTracePts; ++k)
        {
            m_muavTrace[k] = 0.5 * (Fwd[k] + Bwd[k]) ;
        }
    }
    
}
