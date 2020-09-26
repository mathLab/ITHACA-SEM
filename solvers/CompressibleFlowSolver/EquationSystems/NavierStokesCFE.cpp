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
            int nDim = m_fields[0]->GetCoordim(0);
            if (nDim == 2)
            {
                m_C0Project2DExp = MemoryManager<MultiRegions::ContField2D>::
                AllocateSharedPtr(m_session,m_graph,m_session->GetVariable(0));
            }
            else if(nDim == 3)
            {
                m_C0Project3DExp = MemoryManager<MultiRegions::ContField3D>::
                AllocateSharedPtr(m_session,m_graph,m_session->GetVariable(0));
            }
            else
            {
                ASSERTL0(false, "AV C0 smoothing not implemented in 1D.")
            }
        }
        // load physical sensor type
        m_session->LoadSolverInfo("PhysicalSensorType", m_physicalSensorType,
            "Off");


        string diffName, advName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

        m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance(diffName, diffName);
        if ("InteriorPenalty" == diffName)
        {
            SetBoundaryConditionsBwdWeight();
        }

        if (m_specHP_dealiasing)
        {
            m_diffusion->SetFluxVectorNS(
                &NavierStokesCFE::v_GetViscousFluxVectorDeAlias,
                this);

            m_diffusion->SetDiffusionFluxCons(
                &NavierStokesCFE::GetViscousFluxVectorConservVar, this);

            m_diffusion->SetFunctorDerivBndCond(
                &NavierStokesCFE::SetBoundaryConditionsDeriv, this);

            m_diffusion->SetSpecialBndTreat(
                &NavierStokesCFE::SpecialBndTreat, this);

            m_diffusion->SetDiffusionSymmFluxCons(
                &NavierStokesCFE::GetViscousSymmtrFluxConservVar, this);

            if (m_artificialDiffusion)
            {
                m_diffusion->SetArtificialDiffusionVector(
                    &NavierStokesCFE::GetArtificialViscosity, this);
            }
        }
        else
        {
            m_diffusion->SetFluxVectorNS(&NavierStokesCFE::
                                          v_GetViscousFluxVector, this);
            m_diffusion->SetDiffusionFluxCons(
                &NavierStokesCFE::GetViscousFluxVectorConservVar, this);

            m_diffusion->SetSpecialBndTreat(
                &NavierStokesCFE::SpecialBndTreat, this);

            m_diffusion->SetDiffusionSymmFluxCons(
                &NavierStokesCFE::GetViscousSymmtrFluxConservVar, this);

            if (m_artificialDiffusion)
            {
                m_diffusion->SetArtificialDiffusionVector(
                    &NavierStokesCFE::GetArtificialViscosity, this);
            }
        }

        m_diffusion->SetCalcViscosity(
                &NavierStokesCFE::CalcViscosity, this);

        // Concluding initialisation of diffusion operator
        m_diffusion->InitObject         (m_session, m_fields);

        //Only NavierStokes equation and using weakDG,LDGNS can temparary use the codes
        m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
        if(m_useUnifiedWeakIntegration)
        // Set up penalty term for LDGNS
        if ("LDGNS" == diffName||
            "LDGNS3DHomogeneous1D" == diffName)
        {
            m_diffusion->SetFluxPenaltyNS(&NavierStokesCFE::
                v_GetFluxPenalty, this);
        }

        m_GetdFlux_dDeriv_Array = Array<OneD, GetdFlux_dDeriv> (m_spacedim);
        switch (m_spacedim)
        {
        case 2:
            /* code */
            m_GetdFlux_dDeriv_Array[0] = std::bind(
            &NavierStokesCFE::GetdFlux_dQx_2D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);

            m_GetdFlux_dDeriv_Array[1] = std::bind(
            &NavierStokesCFE::GetdFlux_dQy_2D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);
            break;
        case 3:
            /* code */
            m_GetdFlux_dDeriv_Array[0] = std::bind(
            &NavierStokesCFE::GetdFlux_dQx_3D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);

            m_GetdFlux_dDeriv_Array[1] = std::bind(
            &NavierStokesCFE::GetdFlux_dQy_3D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);
            m_GetdFlux_dDeriv_Array[2] = std::bind(
            &NavierStokesCFE::GetdFlux_dQz_3D, this, std::placeholders::_1,
                                                     std::placeholders::_2,
                                                     std::placeholders::_3,
                                                     std::placeholders::_4);

            break;

        default:

            break;
        }
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
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.num_elements());

            for (int i = 0; i < m_fields.num_elements(); ++i)
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
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        size_t nvariables = inarray.size();
        size_t npoints    = GetNpoints();
        size_t nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }

        // get artificial viscosity
        if (m_shockCaptureType == "Physical" && m_calcuPhysicalAV)
        {
            GetPhysicalAV(inarray);
        }

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
        if ("InteriorPenalty" == diffName)
        {
            if (m_BndEvaluateTime < 0.0)
            {
                ASSERTL0(false, "m_BndEvaluateTime not setup");
            }
            m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff,
                                m_BndEvaluateTime,
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


    void NavierStokesCFE::v_DoDiffusion_coeff(
        const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
        Array<OneD, Array<OneD, NekDouble> >                &outarray,
        const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
    {
        size_t nvariables = inarray.size();
        size_t npoints    = GetNpoints();
        size_t ncoeffs    = GetNcoeffs();
        size_t nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff{nvariables};
        for (int i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>{ncoeffs, 0.0};
        }

        // get artificial viscosity
        if (m_shockCaptureType == "Physical" && m_calcuPhysicalAV)
        {
            GetPhysicalAV(inarray);
        }

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
        if ("InteriorPenalty" == diffName)
        {
            if (m_BndEvaluateTime < 0.0)
            {
                ASSERTL0(false, "m_BndEvaluateTime not setup");
            }
            m_diffusion->DiffuseCoeffs(nvariables, m_fields, inarray,
                                        outarrayDiff, m_BndEvaluateTime,
                                        pFwd, pBwd);
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(ncoeffs,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
        }
        else
        {
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff{nvariables - 1};
            Array<OneD, Array<OneD, NekDouble> > inFwd{nvariables - 1};
            Array<OneD, Array<OneD, NekDouble> > inBwd{nvariables-1};

            for (int i = 0; i < nvariables-1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>{npoints};
                inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
                inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
            }

            // Extract pressure
            //    (use inarrayDiff[0] as a temporary storage for the pressure)
            m_varConv->GetPressure(inarray, inarrayDiff[0]);

            // Extract temperature
            m_varConv->GetTemperature(inarray, inarrayDiff[nvariables-2]);

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
                m_varConv->GetPressure(pFwd,    inFwd[0]);
                m_varConv->GetPressure(pBwd,    inBwd[0]);

                m_varConv->GetTemperature(pFwd, inFwd[nvariables-2]);
                m_varConv->GetTemperature(pBwd, inBwd[nvariables-2]);

                m_varConv->GetVelocityVector(pFwd, inFwd);
                m_varConv->GetVelocityVector(pBwd, inBwd);
            }

            // Diffusion term in physical rhs form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff,
                                inFwd, inBwd);

            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }

            if (m_shockCaptureType != "Off")
            {
                m_artificialDiffusion->DoArtificialDiffusion_coeff(
                    inarray, outarray);
            }
        }
    }


    void NavierStokesCFE::v_DoDiffusionFlux(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>&VolumeFlux,
        Array<OneD, Array<OneD, NekDouble>>             &TraceFlux,
        const Array<OneD, Array<OneD, NekDouble>>       &pFwd,
        const Array<OneD, Array<OneD, NekDouble>>       &pBwd)
    {
        int nDim= m_fields[0]->GetCoordim(0);
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints,0.0);
        }

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
        if("InteriorPenalty"==diffName)
        {
            Array<OneD,Array<OneD, Array<OneD, NekDouble>>> inarrayDiffderivative(nDim);
            for (int i = 0; i < nDim; i++)
            {
                inarrayDiffderivative[i]=Array<OneD, Array<OneD, NekDouble>> (nvariables);
                for(int j=0;j<nvariables;j++)
                {
                    inarrayDiffderivative[i][j]=Array<OneD, NekDouble>(npoints,0.0);
                }
            }
            m_diffusion->DiffuseCalculateDerivative(nvariables,m_fields,inarray,inarrayDiffderivative,pFwd,pBwd);
            m_diffusion->DiffuseVolumeFlux(nvariables, m_fields, inarray,inarrayDiffderivative, VolumeFlux);
            m_diffusion->DiffuseTraceFlux(nvariables, m_fields, inarray,inarrayDiffderivative,VolumeFlux,TraceFlux,pFwd,pBwd);
        }
        else
        {
            Array<OneD, Array<OneD, NekDouble>> inarrayDiff;
            Array<OneD, Array<OneD, NekDouble>> inFwd;
            Array<OneD, Array<OneD, NekDouble>> inBwd;
            Array<OneD,Array<OneD, Array<OneD, NekDouble>>> inarrayDiffderivative(nDim);

            //Allocate memory
            inarrayDiff=Array<OneD, Array<OneD, NekDouble>> (nvariables - 1);
            inFwd=Array<OneD, Array<OneD, NekDouble>>(nvariables - 1);
            inBwd=Array<OneD, Array<OneD, NekDouble>> (nvariables - 1);
            inarrayDiffderivative=Array<OneD,Array<OneD, Array<OneD, NekDouble>>> (nDim);

            for (int i = 0; i < nvariables-1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>(npoints,0.0);
                inFwd[i]       = Array<OneD, NekDouble>(nTracePts,0.0);
                inBwd[i]       = Array<OneD, NekDouble>(nTracePts,0.0);
            }

            for (int i = 0; i < nDim; i++)
            {
                inarrayDiffderivative[i]=Array<OneD, Array<OneD, NekDouble>> (nvariables-1);
                for(int j=0;j<nvariables-1;j++)
                {
                    inarrayDiffderivative[i][j]=Array<OneD, NekDouble>(npoints,0.0);
                }
            }

            // Extract pressure
            //    (use inarrayDiff[0] as a temporary storage for the pressure)
            m_varConv->GetPressure(inarray, inarrayDiff[0]);

            // Extract temperature
            m_varConv->GetTemperature(inarray, inarrayDiff[nvariables - 2]);

            // Extract velocities
            m_varConv->GetVelocityVector(inarray, inarrayDiff);

            // Repeat calculation for trace space
            if (pFwd == NullNekDoubleArrayofArray || pBwd == NullNekDoubleArrayofArray)
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
            // To notice, needs to firstly calculate volumeflux, traceflux uses it.
            m_diffusion->DiffuseCalculateDerivative(nvariables,m_fields,inarrayDiff,inarrayDiffderivative,inFwd,inBwd);
            m_diffusion->DiffuseVolumeFlux(nvariables, m_fields, inarrayDiff,inarrayDiffderivative, VolumeFlux);
            m_diffusion->DiffuseTraceFlux(nvariables, m_fields, inarrayDiff,inarrayDiffderivative,VolumeFlux,TraceFlux, inFwd, inBwd);

            //Artificial Diffusion need to implement
            if (m_artificialDiffusion)
            {
                m_artificialDiffusion->DoArtificialDiffusion_coeff(
                    inarray, outarray);
            }
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::v_GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >              &physfield,
        TensorOfArray3D<NekDouble>                              &derivativesO1,
        TensorOfArray3D<NekDouble>                              &viscousTensor)
    {
        // Auxiliary variables
        size_t nScalar    = physfield.size();
        size_t nPts       = physfield[0].size();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, m_mu);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
        }

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            // Apply Ducros sensor
            if (m_physicalSensorType == "Ducros" && m_calcuPhysicalAV)
            {
                Ducros(m_muav);
            }
            // Apply approximate c0 smoothing
            if (m_smoothing == "C0" && m_calcuPhysicalAV)
            {
                C0Smooth(m_muav);
            }
            Vmath::Vadd(nPts, mu, 1, m_muav, 1, mu, 1);
            // Freeze AV for Implicit time stepping
            if (m_explicitDiffusion == false)
            {
                m_calcuPhysicalAV = false;
            }
        }

        // Set thermal conductivity
        NekDouble tRa = m_Cp / m_Prandtl;
        Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);

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
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              TensorOfArray3D<NekDouble>                         &derivativesO1,
              TensorOfArray3D<NekDouble>                         &viscousTensor)
    {
        // Factor to rescale 1d points in dealiasing.
        NekDouble OneDptscale = 2;
        // Get number of points to dealias a cubic non-linearity
        size_t nScalar   = physfield.size();
        int nPts      = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        size_t nPts_orig = physfield[0].size();

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, m_mu);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
        }

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            // Apply Ducros sensor
            if (m_physicalSensorType == "Ducros" && m_calcuPhysicalAV)
            {
                Ducros(m_muav);
            }
            // Apply approximate c0 smoothing
            if (m_smoothing == "C0" && m_calcuPhysicalAV)
            {
                C0Smooth(m_muav);
            }
            Vmath::Vadd(nPts, mu, 1, m_muav, 1, mu, 1);
            // Freeze AV for Implicit time stepping
            if (m_explicitDiffusion == false)
            {
                m_calcuPhysicalAV = false;
            }
        }

        // Set thermal conductivity
        NekDouble tRa = m_Cp / m_Prandtl;
        Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);

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
     * @brief Return the flux vector for the IP diffusion problem, based on
     * conservative variables
     */
    void NavierStokesCFE::GetViscousFluxVectorConservVar(
        const int                                              nDim,
        const Array<OneD, Array<OneD, NekDouble> >             &inarray,
        const TensorOfArray3D<NekDouble>                       &qfields,
        TensorOfArray3D<NekDouble>                             &outarray,
        Array< OneD, int >                                     &nonZeroIndex,
        const Array<OneD, Array<OneD, NekDouble> >             &normal,
        const Array<OneD, NekDouble>                           &ArtifDiffFactor)
    {
        size_t nConvectiveFields   = inarray.size();
        size_t nPts=inarray[0].size();
        int n_nonZero   =   nConvectiveFields - 1;
        TensorOfArray3D<NekDouble> fluxVec;
        Array<OneD, Array<OneD, NekDouble>> outtmp{nConvectiveFields};

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical" && m_calcuPhysicalAV)
        {
            // Apply Ducros sensor
            if (m_physicalSensorType == "Ducros")
            {
                Ducros(m_muav);
            }
            // Apply approximate c0 smoothing
            if (m_smoothing == "C0")
            {
                C0Smooth(m_muav);
            }
            // Update trace
            GetTracePhysicalAV();
            // Freeze AV for Implicit time stepping
            if (m_explicitDiffusion == false)
            {
                m_calcuPhysicalAV = false;
            }
        }

        for (int i = 0; i < nConvectiveFields; ++i)
        {
            outtmp[i]=Array<OneD, NekDouble>{nPts, 0.0};
        }

        for (int i = 0; i < outarray.size(); ++i)
        {
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                Vmath::Zero(nPts, outarray[i][j], 1);
            }
        }


        int nwspD1 = 2*m_spacedim+4;
        Array<OneD, Array<OneD, NekDouble > > auxVars (nwspD1);
        for(int i=0; i<nwspD1;i++)
        {
            auxVars[i]  =   Array<OneD, NekDouble > (nPts,0.0);
        }
        Array<OneD, NekDouble > mu(nPts,0.0);

        CalcAuxiVarForBilinearFom(nConvectiveFields,inarray,mu,auxVars);

        if (normal.size())
        {
            for (int nd = 0; nd < nDim; ++nd)
            {
                for (int nderiv = 0; nderiv < nDim; ++nderiv)
                {
                    GetViscousFluxBilinearForm(nConvectiveFields, nd, nderiv,
                                               inarray, qfields[nderiv],
                                               outtmp);

                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        Vmath::Vvtvp(nPts, &normal[nd][0], 1,
                                    &outtmp[j][0], 1,
                                    &outarray[0][j][0], 1,
                                    &outarray[0][j][0], 1);
                    }
                }
            }
        }
        else
        {
         fluxVec = outarray;
         
         for (int nd = 0; nd < nDim; ++nd)
         {
             for (int nderiv = 0; nderiv < nDim; ++nderiv)
             {
                    GetViscousFluxBilinearForm(nDim, nd, nderiv, inarray,
                                                qfields[nderiv], outtmp);

                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        Vmath::Vadd(nPts, &outtmp[j][0], 1,
                                    &fluxVec[nd][j][0], 1,
                                    &fluxVec[nd][j][0], 1);
                    }
                }
            }
        }

        if (ArtifDiffFactor.size())
        {
            n_nonZero   =   nConvectiveFields;

            if (normal.size())
            {
                Array<OneD, NekDouble> tmparray{nPts, 0.0};
                for (int i = 0; i < nDim; ++i)
                {
                    Vmath::Vmul(nPts, ArtifDiffFactor, 1, normal[i], 1,
                                tmparray, 1);
                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        Vmath::Vvtvp(nPts, tmparray, 1, qfields[i][j], 1,
                                    outarray[0][j], 1, outarray[0][j], 1);
                    }
                }
            }
            else
            {
                for (int i = 0; i < nDim; ++i)
                {
                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        Vmath::Vvtvp(nPts, ArtifDiffFactor, 1, qfields[i][j], 1,
                                    outarray[i][j], 1, outarray[i][j], 1);
                    }
                }
            }
        }

        nonZeroIndex = Array< OneD, int > {size_t(n_nonZero), 0};
        for (int i = 1; i < n_nonZero + 1; ++i)
        {
            nonZeroIndex[n_nonZero - i] =   nConvectiveFields - i;
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
              Array<OneD,       Array<OneD, NekDouble> >    &consvar)
    {
        size_t nConvectiveFields = consvar.size();
        int ndens       = 0;
        int nengy       = nConvectiveFields - 1;
        int nvelst      = ndens + 1;
        int nveled      = nengy;

        int cnt;
        int j, e;
        int id2;

        int nBndEdgePts, nBndEdges, nBndRegions;

        NekDouble InternalEnergy  =   m_eos->GetInternalEnergy(m_Twall);

        Array<OneD, NekDouble> wallTotEngy;



        Array<OneD, Array<OneD, NekDouble>> bndCons {nConvectiveFields};

        Array<OneD, NekDouble> bndTotEngy;
        Array<OneD, NekDouble> bndPressure;
        Array<OneD, NekDouble> bndRho;
        Array<OneD, NekDouble> bndIntEndy;

        int nLengthArray = 0;

        // Compute boundary conditions  for Energy
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

                    for (int k = 0; k < nConvectiveFields; ++k)
                    {
                        Vmath::Vcopy(nBndEdgePts, &consvar[k][id2], 1,
                                    &bndCons[k][0], 1);
                    }

                    m_varConv->GetPressure(bndCons, bndPressure);
                    Vmath::Fill(nLengthArray, m_Twall, bndTotEngy, 1);
                    m_varConv->GetRhoFromPT(bndPressure, bndTotEngy, bndRho);
                    m_varConv->GetEFromRhoP(bndRho, bndPressure, bndIntEndy);
                    m_varConv->GetDynamicEnergy(bndCons, bndTotEngy);

                    Vmath::Vvtvp(nBndEdgePts, &bndIntEndy[0], 1,
                                &bndCons[ndens][0], 1, &bndTotEngy[0], 1,
                                &bndTotEngy[0], 1);


                    Vmath::Vcopy(nBndEdgePts,
                                &bndTotEngy[0], 1,
                                &consvar[nengy][id2], 1);
                }
            }
        }
    }

    /**
     * @brief aplly Neuman boundary conditions on flux
     *        Currently only consider WallAdiabatic
     *
     */
    void NavierStokesCFE::ApplyFluxBndConds(
        const int                                           nConvectiveFields,
              Array<OneD,       Array<OneD, NekDouble> >    &flux)
    {
        int ndens       = 0;
        int nengy       = nConvectiveFields-1;
        int nvelst      = ndens + 1;
        int nveled      = nengy;

        int cnt;
        int j, e;
        int id2;

        int nBndEdgePts, nBndEdges, nBndRegions;

        int nLengthArray    =0;

        // Compute boundary conditions  for Energy
        cnt = 0;
        nBndRegions = m_fields[nengy]->
        GetBndCondExpansions().num_elements();
        for (j = 0; j < nBndRegions; ++j)
        {
            if (m_fields[nengy]->GetBndConditions()[j]->
                GetBoundaryConditionType() ==
                SpatialDomains::ePeriodic)
            {
                continue;
            }

            nBndEdges = m_fields[nengy]->
            GetBndCondExpansions()[j]->GetExpSize();
            for (e = 0; e < nBndEdges; ++e)
            {
                nBndEdgePts = m_fields[nengy]->
                GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                id2 = m_fields[0]->GetTrace()->
                GetPhys_Offset(m_fields[0]->GetTraceMap()->
                            GetBndCondTraceToGlobalTraceMap(cnt++));

                // Imposing Temperature Twall at the wall
                if (boost::iequals(m_fields[nengy]->GetBndConditions()[j]->
                    GetUserDefined(),"WallAdiabatic"))
                {
                    Vmath::Zero(nBndEdgePts, &flux[nengy][id2], 1);
                }
            }
        }
    }

    /**
     * @brief Calculate and return the ArtificialViscosity for shock-capturing.
     */
    void NavierStokesCFE::GetArtificialViscosity(
        const Array<OneD, Array<OneD, NekDouble> >  &inarray,
              Array<OneD,             NekDouble  >  &muav)
    {
        m_artificialDiffusion->GetArtificialViscosity(inarray, muav);
    }

    /**
     * @brief Calculate and return the Symmetric flux in IP method.
     */
    void NavierStokesCFE::GetViscousSymmtrFluxConservVar(
        const int                                           nSpaceDim,
        const Array<OneD, Array<OneD, NekDouble> >          &inaverg,
        const Array<OneD, Array<OneD, NekDouble > >         &inarray,
        TensorOfArray3D<NekDouble>                          &outarray,
        Array< OneD, int >                                  &nonZeroIndex,
        const Array<OneD, Array<OneD, NekDouble> >          &normals)
    {
        size_t nConvectiveFields   = inarray.size();
        size_t nPts                = inaverg[nConvectiveFields - 1].size();
        nonZeroIndex = Array<OneD, int>{nConvectiveFields - 1, 0};
        for (int i = 0; i < nConvectiveFields - 1; ++i)
        {
            nonZeroIndex[i] =   i + 1;
        }
        Array<OneD, Array<OneD, NekDouble> > outtmp{nConvectiveFields};
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            outtmp[i] =   Array<OneD, NekDouble> {nPts, 0.0};
        }

        int nwspD1 = 2*m_spacedim+4;
        Array<OneD, Array<OneD, NekDouble > > auxVars (nwspD1);
        for(int i=0; i<nwspD1;i++)
        {
            auxVars[i]  =   Array<OneD, NekDouble > (nPts,0.0);
        }
        Array<OneD, NekDouble > mu(nPts,0.0);

        CalcAuxiVarForBilinearFom(nConvectiveFields,inaverg,mu,auxVars);

        for (int nd = 0; nd < nSpaceDim; ++nd)
        {
            for (int nderiv = 0; nderiv < nSpaceDim; ++nderiv)
            {
                GetViscousFluxBilinearForm(nSpaceDim, nd, nderiv, inaverg,
                                            inarray, outtmp);

                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Vmath::Vvtvp(nPts, &outtmp[i][0], 1, &normals[nderiv][0], 1,
                                &outarray[nd][i][0], 1, &outarray[nd][i][0], 1);
                }
            }
        }
    }

    void NavierStokesCFE::CalcAuxiVarForBilinearFom(
        const int                                                       nConvectiveFields,
        const Array<OneD, const Array<OneD, NekDouble> >                &inaverg,
        Array<OneD, NekDouble>                                          &mu,
        Array<OneD, Array<OneD, NekDouble> >                            &auxVars)
    {
        int nPts                = inaverg[nConvectiveFields-1].num_elements();
        int nDim=m_spacedim;

        CalcViscosity(inaverg,mu);

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            Array<OneD, NekDouble> muav;
            if (m_fields[0]->GetTrace()->GetTotPoints()==nPts)
            {
                muav = m_muavTrace;
            }
            else
            {
                muav = m_muav;
            }
            Vmath::Vadd(nPts, mu, 1, muav, 1, mu, 1);
        }

        int nAuxVars_count = 0;

        //TODO: to get primary variable outside.(even in the beginning of DoODERhs)
        Array<OneD,Array<OneD,NekDouble>> u(nDim);
        Array<OneD,Array<OneD,NekDouble>> u2(nDim);
        for(int i=0;i<nDim;i++)
        {
            u[i]    =   auxVars[nAuxVars_count];
            nAuxVars_count++;
        }

        for(int i=0;i<nDim;i++)
        {
            u2[i]   =   auxVars[nAuxVars_count];
            nAuxVars_count++;
        }
        Array<OneD,NekDouble> q2    =   auxVars[nAuxVars_count];
        nAuxVars_count++;
        Array<OneD,NekDouble> E_minus_q2=   auxVars[nAuxVars_count];
        nAuxVars_count++;
        Array<OneD,NekDouble> orho=   auxVars[nAuxVars_count];
        nAuxVars_count++;
        Array<OneD,NekDouble>tmp=   auxVars[nAuxVars_count];
        nAuxVars_count++;

        int nDim_plus_one=nDim+1;
        Vmath::Zero(nPts,&q2[0],1);
        Vmath::Sdiv(nPts,1.0,&inaverg[0][0],1,&orho[0],1);
        for(int i=0;i<nDim;i++)

    /**
     * @brief Calculate diffusion flux using the Jacobian form.
     */
    void NavierStokesCFE::GetViscousFluxBilinearForm(
        const int                                           nSpaceDim,
        const int                                           FluxDirection,
        const int                                           DerivDirection,
        const Array<OneD, const Array<OneD, NekDouble> >    &inaverg,
        const Array<OneD, const Array<OneD, NekDouble> >    &injumpp,
        Array<OneD, Array<OneD, NekDouble> >                &outarray)
    {
        size_t nConvectiveFields   = inaverg.size();
        size_t nPts = inaverg[nConvectiveFields - 1].size();
        size_t nDim=nSpaceDim;

        Array<OneD, NekDouble > temperature        {nPts, 0.0};
        Array<OneD, NekDouble > mu                 {nPts, 0.0};
        Array<OneD, NekDouble > thermalConductivity        {nPts, 0.0};
        m_varConv->GetTemperature(inaverg, temperature);
        GetViscosityAndThermalCondFromTemp(temperature, mu,
                                            thermalConductivity);

        Array<OneD, Array<OneD, NekDouble>> outtmp = outarray;
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Zero(nPts, &outarray[i][0], 1);
        }

        Array<OneD, Array<OneD, NekDouble>> u{nDim};
        Array<OneD, Array<OneD, NekDouble>> u2{nDim};
        for (int i = 0; i < nDim; ++i)
        {
            u[i]=Array<OneD, NekDouble> {nPts, 0.0};
            u2[i]=Array<OneD, NekDouble> {nPts, 0.0};
        }
        Array<OneD, NekDouble> q2{nPts, 0.0};
        Array<OneD, NekDouble> E_minus_q2{nPts, 0.0};
        Array<OneD, NekDouble> orho{nPts, 0.0};
        Array<OneD, NekDouble>tmp {nPts, 0.0};
        Array<OneD, NekDouble>tmp1{nPts, 0.0};

        //Constants
        int nDim_plus_one = nDim + 1;
        int FluxDirection_plus_one = FluxDirection + 1;
        NekDouble gamma = m_gamma;
        NekDouble Pr = m_Prandtl;
        NekDouble gammaoPr = gamma / Pr;
        NekDouble one_minus_gammaoPr = 1.0 - gammaoPr;
        const NekDouble OneThird = 1. / 3.;
        const NekDouble TwoThird = 2. * OneThird;
        const NekDouble FourThird = 4. * OneThird;

        Vmath::Sdiv(nPts, 1.0, &inaverg[0][0], 1, &orho[0], 1);
        m_varConv->GetVelocityVector(inaverg, u);
        for (int i = 0; i < nDim; ++i)
        {
            Vmath::Vmul(nPts, &u[i][0], 1, &u[i][0], 1, &u2[i][0], 1);
            Vmath::Vadd(nPts, &q2[0], 1, &u2[i][0], 1, &q2[0], 1);
        }
        Vmath::Vmul(nPts, &inaverg[nDim_plus_one][0], 1,
                    &orho[0], 1, &E_minus_q2[0], 1);
        Vmath::Vsub(nPts, &E_minus_q2[0], 1, &q2[0], 1, &E_minus_q2[0], 1);
        Vmath::Vmul(nPts, &mu[0], 1, &orho[0], 1, &tmp[0], 1);
        int DerivDirection_plus_one = DerivDirection + 1;
        if (DerivDirection == FluxDirection)
        {
            Vmath::Svtvp(nPts, OneThird, &u2[FluxDirection][0], 1,
                        &q2[0], 1, &tmp1[0], 1);
            Vmath::Svtvp(nPts, gammaoPr, &E_minus_q2[0], 1,
                        &tmp1[0], 1,
                        &tmp1[0], 1);
            Vmath::Vmul(nPts, &tmp1[0], 1, &injumpp[0][0], 1, &tmp1[0], 1);
            //orho is tmperary array
            Vmath::Svtvm(nPts, gammaoPr, &injumpp[nDim_plus_one][0], 1,
                        &tmp1[0], 1, &orho[0], 1);

            for (int i = 0; i < nDim; ++i)
            {
                int i_plus_one = i + 1;
                //flux[rhou, rhov, rhow]
                Vmath::Vvtvm(nPts, &u[i][0], 1, &injumpp[0][0], 1,
                            &injumpp[i_plus_one][0], 1,
                            &outtmp[i_plus_one][0], 1);
                Vmath::Neg(nPts, &outtmp[i_plus_one][0], 1);
                Vmath::Vmul(nPts, &tmp[0], 1, &outtmp[i_plus_one][0], 1,
                            &outtmp[i_plus_one][0], 1);
                //flux rhoE
                Vmath::Smul(nPts, one_minus_gammaoPr, &u[i][0], 1, &tmp1[0], 1);
                Vmath::Vvtvp(nPts, &tmp1[0], 1, &injumpp[i_plus_one][0], 1,
                            &outtmp[nDim_plus_one][0], 1,
                            &outtmp[nDim_plus_one][0], 1);

                if (i == FluxDirection)
                {
                    Vmath::Smul(nPts, FourThird, &outtmp[i_plus_one][0], 1,
                                &outtmp[i_plus_one][0], 1);
                    Vmath::Smul(nPts, OneThird,
                                &u[FluxDirection][0], 1, &tmp1[0], 1);
                    Vmath::Vvtvp(nPts, &tmp1[0], 1,
                                &injumpp[FluxDirection_plus_one][0], 1,
                                &outtmp[nDim_plus_one][0], 1,
                                &outtmp[nDim_plus_one][0], 1);
                }
            }
            Vmath::Vadd(nPts, &orho[0], 1, &outtmp[nDim_plus_one][0], 1,
                        &outtmp[nDim_plus_one][0], 1);
            Vmath::Vmul(nPts, &tmp[0], 1, &outtmp[nDim_plus_one][0], 1,
                        &outtmp[nDim_plus_one][0], 1);

        }
        else
        {
            Vmath::Vvtvm(nPts, &u[DerivDirection][0], 1, &injumpp[0][0], 1,
                        &injumpp[DerivDirection_plus_one][0], 1, &tmp1[0], 1);
            Vmath::Smul(nPts, TwoThird, &tmp1[0], 1, &tmp1[0], 1);
            Vmath::Vmul(nPts, &tmp[0], 1, &tmp1[0], 1,
                        &outtmp[FluxDirection_plus_one][0], 1);

            Vmath::Vvtvm(nPts, &u[FluxDirection][0], 1, &injumpp[0][0], 1,
                        &injumpp[FluxDirection_plus_one][0], 1, &tmp1[0], 1);
            Vmath::Neg(nPts, &tmp1[0], 1);
            Vmath::Vmul(nPts, &tmp[0], 1, &tmp1[0], 1,
                        &outtmp[DerivDirection_plus_one][0], 1);

            Vmath::Smul(nPts, OneThird, &u[FluxDirection][0], 1, &tmp1[0], 1);
            Vmath::Vmul(nPts, &tmp1[0], 1, &u[DerivDirection][0], 1,
                        &tmp1[0], 1);
            Vmath::Vmul(nPts, &tmp1[0], 1, &injumpp[0][0], 1, &tmp1[0], 1);
            //previous orho as a tmperary memory because it is
            //non-used any more
            Vmath::Smul(nPts, TwoThird, &u[FluxDirection][0], 1, &orho[0], 1);
            Vmath::Vmul(nPts, &orho[0], 1,
                        &injumpp[DerivDirection_plus_one][0], 1,
                        &orho[0], 1);
            Vmath::Vadd(nPts, &tmp1[0], 1, &orho[0], 1, &tmp1[0], 1);
            Vmath::Neg(nPts, &tmp1[0], 1);
            Vmath::Vvtvp(nPts, &u[DerivDirection][0], 1,
                        &injumpp[FluxDirection_plus_one][0], 1, &tmp1[0], 1,
                        &tmp1[0], 1);
            Vmath::Vmul(nPts, &tmp[0], 1, &tmp1[0], 1,
                        &outtmp[nDim_plus_one][0], 1);
        }
    }
        
    void NavierStokesCFE::CalcViscosity(
        const Array<OneD, const Array<OneD, NekDouble> >                &inaverg,
        Array<OneD, NekDouble>                                          &mu)
    {
        int nConvectiveFields       = inaverg.num_elements();
        int nPts                = inaverg[nConvectiveFields-1].num_elements();
        int nDim=m_spacedim;

        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble> tmp(nPts,0.0);
            m_varConv->GetTemperature(inaverg,tmp);
            m_varConv->GetDynamicViscosity(tmp, mu);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
        }
    }

    /**
     * @brief return part of viscous Jacobian:
     * \todo flux derived with Qx=[drho_dx,drhou_dx,drhov_dx,drhoE_dx]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhoE]
     * Output: 2D 3*4 Matrix (flux with rho is zero)
     */
    void NavierStokesCFE::GetdFlux_dQx_2D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble E=U[3]*orho;
        NekDouble q2=u*u+v*v;
        NekDouble e=E-0.5*q2;
        NekDouble gamma=m_gamma;
        NekDouble Cp=m_Cp;
        NekDouble Cv=m_Cv;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr=  1.0/Pr;
        NekDouble tRa = m_Cp *oPr;
        NekDouble kappa=mu*tRa;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I Do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();
        int ncol = OutputMatrix->GetColumns();

        tmpArray[0+0*nrow]=tmp*(-FourThird*u*nx-v*ny);
        tmpArray[0+1*nrow]=tmp*(FourThird*nx);
        tmpArray[0+2*nrow]=tmp*ny;
        tmpArray[0+3*nrow]=0.0;
        tmpArray[1+0*nrow]=tmp*(-v*nx+TwoThird*u*ny);
        tmpArray[1+1*nrow]=tmp*(-TwoThird*ny);
        tmpArray[1+2*nrow]=tmp*nx;
        tmpArray[1+3*nrow]=0.0;
        tmpArray[2+0*nrow]=(FourThird*u*u+v*v+tmp2*(E-q2))*nx+OneThird*u*v*ny;
        tmpArray[2+0*nrow]=-tmp*(*OutputMatrix)(2,0);
        tmpArray[2+1*nrow]=(FourThird-tmp2)*u*nx-TwoThird*v*ny;
        tmpArray[2+1*nrow]=tmp*(*OutputMatrix)(2,1);
        tmpArray[2+2*nrow]=(1-tmp2)*v*nx+u*ny;
        tmpArray[2+2*nrow]=tmp*(*OutputMatrix)(2,2);
        tmpArray[2+3*nrow]=tmp*tmp2*nx;
    }

     /**
     * @brief return part of viscous Jacobian:
     * \todo flux derived with Qx=[drho_dy,drhou_dy,drhov_dy,drhoE_dy]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhoE]
     * Output: 2D 3*4 Matrix (flux with rho is zero)
     */
    void NavierStokesCFE::GetdFlux_dQy_2D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble E=U[3]*orho;
        NekDouble q2=u*u+v*v;
        NekDouble e=E-0.5*q2;
        NekDouble gamma=m_gamma;
        NekDouble Cp=m_Cp;
        NekDouble Cv=m_Cv;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr=  1.0/Pr;
        NekDouble tRa = m_Cp *oPr;
        NekDouble kappa=mu*tRa;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I Do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();
        int ncol = OutputMatrix->GetColumns();
           
        tmpArray[0+0*nrow]=tmp*(TwoThird*v*nx-u*ny);
        tmpArray[0+1*nrow]=tmp*ny;
        tmpArray[0+2*nrow]=tmp*(-TwoThird)*nx;
        tmpArray[0+3*nrow]=0.0;
        tmpArray[1+0*nrow]=tmp*(-u*nx-FourThird*v*ny);
        tmpArray[1+1*nrow]=tmp*nx;
        tmpArray[1+2*nrow]=tmp*(FourThird*ny);
        tmpArray[1+3*nrow]=0.0;
        tmpArray[2+0*nrow]=OneThird*u*v*nx+(FourThird*v*v+u*u+tmp2*(E-q2))*ny;
        tmpArray[2+0*nrow]=-tmp*(*OutputMatrix)(2,0);
        tmpArray[2+1*nrow]=(1-tmp2)*u*ny+v*nx;
        tmpArray[2+1*nrow]=tmp*(*OutputMatrix)(2,1);
        tmpArray[2+2*nrow]=(FourThird-tmp2)*v*ny-TwoThird*u*nx;
        tmpArray[2+2*nrow]=tmp*(*OutputMatrix)(2,2);
        tmpArray[2+3*nrow]=tmp*tmp2*ny;
    }

    /**
     * @brief return part of viscous Jacobian derived with Qx=[drho_dx,drhou_dx,drhov_dx,drhow_dx,drhoE_dx]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qx=[drho_dx,drhou_dx,drhov_dx,drhow_dx,drhoE_dx]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=0)= dF_dQx;
     */
    void NavierStokesCFE::GetdFlux_dQx_3D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble w=U[3]*orho;
        NekDouble E=U[4]*orho;
        NekDouble q2=u*u+v*v+w*w;
        NekDouble e=E-0.5*q2;
        NekDouble gamma=m_gamma;
        NekDouble Cp=m_Cp;
        NekDouble Cv=m_Cv;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;
        NekDouble tRa = m_Cp *oPr;
        NekDouble kappa=mu*tRa;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();
        int ncol = OutputMatrix->GetColumns();

        tmpArray[0+0*nrow]=tmpx*(-FourThird*u)+tmpy*(-v)+tmpz*(-w);
        tmpArray[0+1*nrow]=tmpx*FourThird;
        tmpArray[0+2*nrow]=tmpy;
        tmpArray[0+3*nrow]=tmpz;
        tmpArray[0+4*nrow]=0.0;
        tmpArray[1+0*nrow]=tmpx*(-v)+tmpy*(TwoThird*u);
        tmpArray[1+1*nrow]=tmpy*(-TwoThird);
        tmpArray[1+2*nrow]=tmpx;
        tmpArray[1+3*nrow]=0.0;
        tmpArray[1+4*nrow]=0.0;
        tmpArray[2+0*nrow]=tmpx*(-w)+tmpz*(TwoThird*u);
        tmpArray[2+1*nrow]=tmpz*(-TwoThird);
        tmpArray[2+2*nrow]=0.0;
        tmpArray[2+3*nrow]=tmpx;
        tmpArray[2+4*nrow]=0.0;
        tmpArray[3+0*nrow]=-tmpx*(FourThird*u*u+v*v+w*w+tmp2*(E-q2))+tmpy*(-OneThird*u*v)+tmpz*(-OneThird*u*w);
        tmpArray[3+1*nrow]=tmpx*(FourThird-tmp2)*u+tmpy*(-TwoThird*v)+tmpz*(-TwoThird*w);
        tmpArray[3+2*nrow]=tmpx*(1.0-tmp2)*v+tmpy*u;
        tmpArray[3+3*nrow]=tmpx*(1.0-tmp2)*w+tmpz*u;
        tmpArray[3+4*nrow]=tmpx*tmp2;
    }

    /**
     * @brief return part of viscous Jacobian derived with Qy=[drho_dy,drhou_dy,drhov_dy,drhow_dy,drhoE_dy]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qy=[drho_dy,drhou_dy,drhov_dy,drhow_dy,drhoE_dy]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=1)= dF_dQy;
     */
    void NavierStokesCFE::GetdFlux_dQy_3D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble w=U[3]*orho;
        NekDouble E=U[4]*orho;
        NekDouble q2=u*u+v*v+w*w;
        NekDouble e=E-0.5*q2;
        NekDouble gamma=m_gamma;
        NekDouble Cp=m_Cp;
        NekDouble Cv=m_Cv;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;
        NekDouble tRa = m_Cp *oPr;
        NekDouble kappa=mu*tRa;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();
        int ncol = OutputMatrix->GetColumns();

        tmpArray[0+0*nrow]=tmpx*(TwoThird*v)+tmpy*(-u);
        tmpArray[0+1*nrow]=tmpy;
        tmpArray[0+2*nrow]=tmpx*(-TwoThird);
        tmpArray[0+3*nrow]=0.0;
        tmpArray[0+4*nrow]=0.0;
        tmpArray[1+0*nrow]=tmpx*(-u)+tmpy*(-FourThird*v)+tmpz*(-w);
        tmpArray[1+1*nrow]=tmpx;
        tmpArray[1+2*nrow]=tmpy*FourThird;
        tmpArray[1+3*nrow]=tmpz;
        tmpArray[1+4*nrow]=0.0;
        tmpArray[2+0*nrow]=tmpy*(-w)+tmpz*(TwoThird*v);
        tmpArray[2+1*nrow]=0.0;
        tmpArray[2+2*nrow]=tmpz*(-TwoThird);
        tmpArray[2+3*nrow]=tmpy;
        tmpArray[2+4*nrow]=0.0;
        tmpArray[3+0*nrow]=tmpx*(-OneThird*u*v)-tmpy*(u*u+FourThird*v*v+w*w+tmp2*(E-q2))+tmpz*(-OneThird*v*w);
        tmpArray[3+1*nrow]=tmpx*v+tmpy*(1-tmp2)*u;
        tmpArray[3+2*nrow]=tmpx*(-TwoThird*u)+tmpy*(FourThird-tmp2)*v+tmpz*(-TwoThird*w);
        tmpArray[3+3*nrow]=tmpy*(1-tmp2)*w+tmpz*v;
        tmpArray[3+4*nrow]=tmpy*tmp2;
    }

    /**
     * @brief return part of viscous Jacobian derived with Qz=[drho_dz,drhou_dz,drhov_dz,drhow_dz,drhoE_dz]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qz=[drho_dz,drhou_dz,drhov_dz,drhow_dz,drhoE_dz]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=2)= dF_dQz;
     */
    void NavierStokesCFE::GetdFlux_dQz_3D(
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U,
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble orho=1.0/rho;
        NekDouble u=U[1]*orho;
        NekDouble v=U[2]*orho;
        NekDouble w=U[3]*orho;
        NekDouble E=U[4]*orho;
        NekDouble q2=u*u+v*v+w*w;
        NekDouble e=E-0.5*q2;
        NekDouble gamma=m_gamma;
        NekDouble Cp=m_Cp;
        NekDouble Cv=m_Cv;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;
        NekDouble tRa = m_Cp *oPr;
        NekDouble kappa=mu*tRa;
        //To notice, here is positive, which is consistent with
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu*orho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;
        NekDouble tmp2=gamma*oPr;
        NekDouble OneThird,TwoThird,FourThird;
        OneThird=1.0/3.0;
        TwoThird=2.0*OneThird;
        FourThird=4.0*OneThird;

        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();
        int ncol = OutputMatrix->GetColumns();

        tmpArray[0+0*nrow]=tmpx*(TwoThird*w)+tmpz*(-u);
        tmpArray[0+1*nrow]=tmpz;
        tmpArray[0+2*nrow]=0.0;
        tmpArray[0+3*nrow]=tmpx*(-TwoThird);
        tmpArray[0+4*nrow]=0.0;
        tmpArray[1+0*nrow]=tmpy*(TwoThird*w)+tmpz*(-v);
        tmpArray[1+1*nrow]=0.0;
        tmpArray[1+2*nrow]=tmpz;
        tmpArray[1+3*nrow]=tmpy*(-TwoThird);
        tmpArray[1+4*nrow]=0.0;
        tmpArray[2+0*nrow]=tmpx*(-u)+tmpy*(-v)+tmpz*(-FourThird*w);
        tmpArray[2+1*nrow]=tmpx;
        tmpArray[2+2*nrow]=tmpy;
        tmpArray[2+3*nrow]=tmpz*FourThird;
        tmpArray[2+4*nrow]=0.0;
        tmpArray[3+0*nrow]=tmpx*(-OneThird*u*w)+tmpy*(-OneThird*v*w)-tmpz*(u*u+v*v+FourThird*w*w+tmp2*(E-q2));
        tmpArray[3+1*nrow]=tmpx*w+tmpz*(1-tmp2)*u;
        tmpArray[3+2*nrow]=tmpy*w+tmpz*(1-tmp2)*v;
        tmpArray[3+3*nrow]=tmpx*(-TwoThird*u)+tmpy*(-TwoThird*v)+tmpz*(FourThird-tmp2)*w;
        tmpArray[3+4*nrow]=tmpz*tmp2;
    }

    /**
     * @brief return part of viscous Jacobian
     * Input:
     * normals:Point normals
     * mu: dynamicviscosity
     * dmu_dT: mu's derivative with T using Sutherland's law
     * U=[rho,rhou,rhov,rhoE]
     * Output: 3*4 Matrix (the flux about rho is zero)
     * OutputMatrix dFLux_dU,  the matrix sign is consistent with SIPG
     */
    void NavierStokesCFE::GetdFlux_dU_2D(
        const Array<OneD, NekDouble>                        &normals,
        const NekDouble                                     mu,
        const NekDouble                                     dmu_dT,
        const Array<OneD, NekDouble>                        &U,
        const Array<OneD, const Array<OneD, NekDouble> >    &qfield,
              DNekMatSharedPtr                              &OutputMatrix)
    {
        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();
        int ncol = OutputMatrix->GetColumns();

        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble U1=U[0];
        NekDouble U2=U[1];
        NekDouble U3=U[2];
        NekDouble U4=U[3];
        NekDouble dU1_dx=qfield[0][0];
        NekDouble dU2_dx=qfield[0][1];
        NekDouble dU3_dx=qfield[0][2];
        NekDouble dU4_dx=qfield[0][3];
        NekDouble dU1_dy=qfield[1][0];
        NekDouble dU2_dy=qfield[1][1];
        NekDouble dU3_dy=qfield[1][2];
        NekDouble dU4_dy=qfield[1][3];
        NekDouble gamma=m_gamma;
        NekDouble Cv=m_Cv;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;

        NekDouble orho1,orho2,orho3,orho4;
        NekDouble oCv=1./Cv;
        orho1=1.0/U1;
        orho2=orho1*orho1;
        orho3=orho1*orho2;
        orho4=orho2*orho2;

        //Assume Fn=mu*Sn
        //Sn=Sx*nx+Sy*ny

        NekDouble TwoThrid=2./3.;
        NekDouble FourThird=2.0*TwoThrid;
        NekDouble u=U2*orho1;
        NekDouble v=U3*orho1;
        NekDouble du_dx=orho1*(dU2_dx-u*dU1_dx);
        NekDouble dv_dx=orho1*(dU3_dx-v*dU1_dx);
        NekDouble du_dy=orho1*(dU2_dy-u*dU1_dy);
        NekDouble dv_dy=orho1*(dU3_dy-v*dU1_dy);
        NekDouble s12=FourThird*du_dx-TwoThrid*dv_dy;
        NekDouble s13=du_dy+dv_dx;
        NekDouble s22=s13;
        NekDouble s23=FourThird*dv_dy-TwoThrid*du_dx;
        NekDouble snx=s12*nx+s22*ny;
        NekDouble sny=s13*nx+s23*ny;
        NekDouble snv=snx*u+sny*v;
        NekDouble qx=-gamma*mu*oPr*(orho1*dU4_dx-U[3]*orho2*dU1_dx-u*(orho1*dU2_dx-U[1]*orho2*dU1_dx)-v*(orho1*dU3_dx-U[2]*orho2*dU1_dx));
        NekDouble qy=-gamma*mu*oPr*(orho1*dU4_dy-U[3]*orho2*dU1_dy-u*(orho1*dU2_dy-U[1]*orho2*dU1_dy)-v*(orho1*dU3_dy-U[2]*orho2*dU1_dy));
        NekDouble qn=qx*nx+qy*ny;

        //Term1 mu's derivative with U: dmu_dU*Sn
        Array<OneD,NekDouble> tmp(3,0.0);
        tmp[0]=snx;
        tmp[1]=sny;
        tmp[2]=snv-qn/mu;
        Array<OneD,NekDouble> dT_dU (4,0.0);
        dT_dU[0]=oCv*(-orho2*U4+orho3*U2*U2+orho3*U3*U3);
        dT_dU[1]=-oCv*orho2*U2;
        dT_dU[2]=-oCv*orho2*U3;
        dT_dU[3]=oCv*orho1;
        for(int i=0;i<3;i++)
            for(int j=0;j<4;j++)
            {
                tmpArray[i+j*nrow]=dmu_dT*dT_dU[j]*tmp[i];
            }
        }

        //Term 2 +mu*dSn_dU
        NekDouble du_dx_dU1,du_dx_dU2;
        NekDouble du_dy_dU1,du_dy_dU2;
        NekDouble dv_dx_dU1,dv_dx_dU3;
        NekDouble dv_dy_dU1,dv_dy_dU3;
        NekDouble ds12_dU1,ds12_dU2,ds12_dU3;
        NekDouble ds13_dU1,ds13_dU2,ds13_dU3;
        NekDouble ds22_dU1,ds22_dU2,ds22_dU3;
        NekDouble ds23_dU1,ds23_dU2,ds23_dU3;
        NekDouble dsnx_dU1,dsnx_dU2,dsnx_dU3;
        NekDouble dsny_dU1,dsny_dU2,dsny_dU3;
        NekDouble dsnv_dU1,dsnv_dU2,dsnv_dU3;

        du_dx_dU1=-orho2*dU2_dx+2*orho3*U2*dU1_dx;
        du_dx_dU2=-orho2*dU1_dx;
        du_dy_dU1=-orho2*dU2_dy+2*orho3*U2*dU1_dy;
        du_dy_dU2=-orho2*dU1_dy;
        dv_dx_dU1=-orho2*dU3_dx+2*orho3*U3*dU1_dx;
        dv_dx_dU3=du_dx_dU2;
        dv_dy_dU1=-orho2*dU3_dy+2*orho3*U3*dU1_dy;
        dv_dy_dU3=du_dy_dU2;
        ds12_dU1=FourThird*du_dx_dU1-TwoThrid*dv_dy_dU1;
        ds12_dU2=FourThird*du_dx_dU2;
        ds12_dU3=-TwoThrid*dv_dy_dU3;
        ds13_dU1=du_dy_dU1+dv_dx_dU1;
        ds13_dU2=du_dy_dU2;
        ds13_dU3=dv_dx_dU3;
        ds22_dU1=ds13_dU1;
        ds22_dU2=ds13_dU2;
        ds22_dU3=ds13_dU3;
        ds23_dU1=FourThird*dv_dy_dU1-TwoThrid*du_dx_dU1;
        ds23_dU2=-TwoThrid*du_dx_dU2;
        ds23_dU3=FourThird*dv_dy_dU3;
        dsnx_dU1=ds12_dU1*nx+ds22_dU1*ny;
        dsnx_dU2=ds12_dU2*nx+ds22_dU2*ny;
        dsnx_dU3=ds12_dU3*nx+ds22_dU3*ny;
        dsny_dU1=ds13_dU1*nx+ds23_dU1*ny;
        dsny_dU2=ds13_dU2*nx+ds23_dU2*ny;
        dsny_dU3=ds13_dU3*nx+ds23_dU3*ny;
        dsnv_dU1=u*dsnx_dU1+v*dsny_dU1-orho2*U2*snx-orho2*U3*sny;
        dsnv_dU2=u*dsnx_dU2+v*dsny_dU2+orho1*snx;
        dsnv_dU3=u*dsnx_dU3+v*dsny_dU3+orho1*sny;
        tmpArray[0+0*nrow]=tmpArray[0+0*nrow]+mu*dsnx_dU1;
        tmpArray[0+1*nrow]=tmpArray[0+1*nrow]+mu*dsnx_dU2;
        tmpArray[0+2*nrow]=tmpArray[0+2*nrow]+mu*dsnx_dU3;
        tmpArray[1+0*nrow]=tmpArray[1+0*nrow]+mu*dsny_dU1;
        tmpArray[1+1*nrow]=tmpArray[1+1*nrow]+mu*dsny_dU2;
        tmpArray[1+2*nrow]=tmpArray[1+2*nrow]+mu*dsny_dU3;
        tmpArray[2+0*nrow]=tmpArray[2+0*nrow]+mu*dsnv_dU1;
        tmpArray[2+1*nrow]=tmpArray[2+1*nrow]+mu*dsnv_dU2;
        tmpArray[2+2*nrow]=tmpArray[2+2*nrow]+mu*dsnv_dU3;

        //Consider +qn's effect (does not include mu's effect)
        NekDouble dqx_dU1,dqx_dU2,dqx_dU3,dqx_dU4;
        NekDouble dqy_dU1,dqy_dU2,dqy_dU3,dqy_dU4;
        NekDouble tmpx=-nx*mu*gamma*oPr;
        dqx_dU1=tmpx*(-orho2*dU4_dx+2*orho3*U4*dU1_dx+2*orho3*U2*dU2_dx-3*orho4*U2*U2*dU1_dx+2*orho3*U3*dU3_dx-3*orho4*U3*U3*dU1_dx);
        dqx_dU2=tmpx*(-orho2*dU2_dx+2*orho3*U2*dU1_dx);
        dqx_dU3=tmpx*(-orho2*dU3_dx+2*orho3*U3*dU1_dx);
        dqx_dU4=-tmpx*orho2*dU1_dx;
        NekDouble tmpy=-ny*mu*gamma*oPr;
        dqy_dU1=tmpy*(-orho2*dU4_dy+2*orho3*U4*dU1_dy+2*orho3*U2*dU2_dy-3*orho4*U2*U2*dU1_dy+2*orho3*U3*dU3_dy-3*orho4*U3*U3*dU1_dy);
        dqy_dU2=tmpy*(-orho2*dU2_dy+2*orho3*U2*dU1_dy);
        dqy_dU3=tmpy*(-orho2*dU3_dy+2*orho3*U3*dU1_dy);
        dqy_dU4=-tmpy*orho2*dU1_dy;
        tmpArray[2+0*nrow]=tmpArray[2+0*nrow]-dqx_dU1-dqy_dU1;
        tmpArray[2+1*nrow]=tmpArray[2+1*nrow]-dqx_dU2-dqy_dU2;
        tmpArray[2+2*nrow]=tmpArray[2+2*nrow]-dqx_dU3-dqy_dU3;
        tmpArray[2+3*nrow]=tmpArray[2+3*nrow]-dqx_dU4-dqy_dU4;
    }

     /**
     * @brief return part of viscous Jacobian
     * Input:
     * normals:Point normals
     * mu: dynamicviscosity
     * dmu_dT: mu's derivative with T using Sutherland's law
     * U=[rho,rhou,rhov,rhow,rhoE]
     * Output: 4*5 Matrix (the flux about rho is zero)
     * OutputMatrix dFLux_dU,  the matrix sign is consistent with SIPG
     */
    void NavierStokesCFE::GetdFlux_dU_3D(
        const Array<OneD, NekDouble>                        &normals,
        const NekDouble                                     mu,
        const NekDouble                                     dmu_dT,
        const Array<OneD, NekDouble>                        &U,
        const Array<OneD, const Array<OneD, NekDouble> >    &qfield,
              DNekMatSharedPtr                              &OutputMatrix)
    {
        Array<OneD, NekDouble> tmpArray;
        tmpArray = OutputMatrix->GetPtr();
        int nrow = OutputMatrix->GetRows();
        int ncol = OutputMatrix->GetColumns();

        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble U1=U[0];
        NekDouble U2=U[1];
        NekDouble U3=U[2];
        NekDouble U4=U[3];
        NekDouble U5=U[4];
        NekDouble dU1_dx=qfield[0][0];
        NekDouble dU2_dx=qfield[0][1];
        NekDouble dU3_dx=qfield[0][2];
        NekDouble dU4_dx=qfield[0][3];
        NekDouble dU5_dx=qfield[0][4];
        NekDouble dU1_dy=qfield[1][0];
        NekDouble dU2_dy=qfield[1][1];
        NekDouble dU3_dy=qfield[1][2];
        NekDouble dU4_dy=qfield[1][3];
        NekDouble dU5_dy=qfield[1][4];
        NekDouble dU1_dz=qfield[2][0];
        NekDouble dU2_dz=qfield[2][1];
        NekDouble dU3_dz=qfield[2][2];
        NekDouble dU4_dz=qfield[2][3];
        NekDouble dU5_dz=qfield[2][4];
        NekDouble gamma=m_gamma;
        NekDouble Cv=m_Cv;
        NekDouble Pr=  m_Prandtl;
        NekDouble oPr = 1.0/Pr;

        NekDouble orho1,orho2,orho3,orho4;
        NekDouble oCv=1./Cv;
        orho1=1.0/U1;
        orho2=orho1*orho1;
        orho3=orho1*orho2;
        orho4=orho2*orho2;

        //Assume Fn=mu*Sn
        //Sn=Sx*nx+Sy*ny+Sz*nz
        NekDouble TwoThrid=2./3.;
        NekDouble FourThird=2.0*TwoThrid;
        NekDouble tmp2=gamma*mu*oPr;
        NekDouble u=U2*orho1;
        NekDouble v=U3*orho1;
        NekDouble w=U4*orho1;
        NekDouble du_dx=orho1*(dU2_dx-u*dU1_dx);
        NekDouble dv_dx=orho1*(dU3_dx-v*dU1_dx);
        NekDouble dw_dx=orho1*(dU4_dx-w*dU1_dx);
        NekDouble du_dy=orho1*(dU2_dy-u*dU1_dy);
        NekDouble dv_dy=orho1*(dU3_dy-v*dU1_dy);
        NekDouble dw_dy=orho1*(dU4_dy-w*dU1_dy);
        NekDouble du_dz=orho1*(dU2_dz-u*dU1_dz);
        NekDouble dv_dz=orho1*(dU3_dz-v*dU1_dz);
        NekDouble dw_dz=orho1*(dU4_dz-w*dU1_dz);
        NekDouble s12=FourThird*du_dx-TwoThrid*dv_dy-TwoThrid*dw_dz;
        NekDouble s13=du_dy+dv_dx;
        NekDouble s14=dw_dx+du_dz;
        NekDouble s22=s13;
        NekDouble s23=FourThird*dv_dy-TwoThrid*du_dx-TwoThrid*dw_dz;
        NekDouble s24=dv_dz+dw_dy;
        NekDouble s32=s14;
        NekDouble s33=s24;
        NekDouble s34=FourThird*dw_dz-TwoThrid*du_dx-TwoThrid*dv_dy;
        NekDouble snx=s12*nx+s22*ny+s32*nz;
        NekDouble sny=s13*nx+s23*ny+s33*nz;
        NekDouble snz=s14*nz+s24*ny+s34*nz;
        NekDouble snv=snx*u+sny*v+snz*w;
        NekDouble qx=-tmp2*(orho1*dU5_dx-U5*orho2*dU1_dx-u*(orho1*dU2_dx-U2*orho2*dU1_dx)-v*(orho1*dU3_dx-U3*orho2*dU1_dx)-w*(orho1*dU4_dx-U4*orho2*dU1_dx));
        NekDouble qy=-tmp2*(orho1*dU5_dy-U5*orho2*dU1_dy-u*(orho1*dU2_dy-U2*orho2*dU1_dy)-v*(orho1*dU3_dy-U3*orho2*dU1_dy)-w*(orho1*dU4_dy-U4*orho2*dU1_dy));
        NekDouble qz=-tmp2*(orho1*dU5_dz-U5*orho2*dU1_dz-u*(orho1*dU2_dz-U2*orho2*dU1_dz)-v*(orho1*dU3_dz-U3*orho2*dU1_dz)-w*(orho1*dU4_dz-U4*orho2*dU1_dz));
        NekDouble qn=qx*nx+qy*ny+qz*nz;

        //Term1 mu's derivative with U: dmu_dU*Sn
        Array<OneD,NekDouble> tmp(4,0.0);
        tmp[0]=snx;
        tmp[1]=sny;
        tmp[2]=snz;
        tmp[3]=snv-qn/mu;
        Array<OneD,NekDouble> dT_dU (5,0.0);
        dT_dU[0]=oCv*(-orho2*U5+orho3*U2*U2+orho3*U3*U3+orho3*U4*U4);
        dT_dU[1]=-oCv*orho2*U2;
        dT_dU[2]=-oCv*orho2*U3;
        dT_dU[3]=-oCv*orho2*U4;
        dT_dU[4]=oCv*orho1;
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<5;j++)
            {
                tmpArray[i+j*nrow]=dmu_dT*dT_dU[j]*tmp[i];
            }
        }

        //Term 2 +mu*dSn_dU
        NekDouble du_dx_dU1,du_dx_dU2;
        NekDouble du_dy_dU1,du_dy_dU2;
        NekDouble du_dz_dU1,du_dz_dU2;
        NekDouble dv_dx_dU1,dv_dx_dU3;
        NekDouble dv_dy_dU1,dv_dy_dU3;
        NekDouble dv_dz_dU1,dv_dz_dU3;
        NekDouble dw_dx_dU1,dw_dx_dU4;
        NekDouble dw_dy_dU1,dw_dy_dU4;
        NekDouble dw_dz_dU1,dw_dz_dU4;
        NekDouble ds12_dU1,ds12_dU2,ds12_dU3,ds12_dU4;
        NekDouble ds13_dU1,ds13_dU2,ds13_dU3;
        NekDouble ds14_dU1,ds14_dU2,ds14_dU4;
        NekDouble ds22_dU1,ds22_dU2,ds22_dU3;
        NekDouble ds23_dU1,ds23_dU2,ds23_dU3,ds23_dU4;
        NekDouble ds24_dU1,ds24_dU3,ds24_dU4;
        NekDouble ds32_dU1,ds32_dU2,ds32_dU4;
        NekDouble ds33_dU1,ds33_dU3,ds33_dU4;
        NekDouble ds34_dU1,ds34_dU2,ds34_dU3,ds34_dU4;
        NekDouble dsnx_dU1,dsnx_dU2,dsnx_dU3,dsnx_dU4;
        NekDouble dsny_dU1,dsny_dU2,dsny_dU3,dsny_dU4;
        NekDouble dsnz_dU1,dsnz_dU2,dsnz_dU3,dsnz_dU4;
        NekDouble dsnv_dU1,dsnv_dU2,dsnv_dU3,dsnv_dU4;

        du_dx_dU1=-orho2*dU2_dx+2*orho3*U2*dU1_dx;
        du_dx_dU2=-orho2*dU1_dx;
        du_dy_dU1=-orho2*dU2_dy+2*orho3*U2*dU1_dy;
        du_dy_dU2=-orho2*dU1_dy;
        du_dz_dU1=-orho2*dU2_dz+2*orho3*U2*dU1_dz;
        du_dz_dU2=-orho2*dU1_dz;
        dv_dx_dU1=-orho2*dU3_dx+2*orho3*U3*dU1_dx;
        dv_dx_dU3=-orho2*dU1_dx;
        dv_dy_dU1=-orho2*dU3_dy+2*orho3*U3*dU1_dy;
        dv_dy_dU3=-orho2*dU1_dy;
        dv_dz_dU1=-orho2*dU3_dz+2*orho3*U3*dU1_dz;
        dv_dz_dU3=-orho2*dU1_dz;
        dw_dx_dU1=-orho2*dU4_dx+2*orho3*U4*dU1_dx;
        dw_dx_dU4=-orho2*dU1_dx;
        dw_dy_dU1=-orho2*dU4_dy+2*orho3*U4*dU1_dy;
        dw_dy_dU4=-orho2*dU1_dy;
        dw_dz_dU1=-orho2*dU4_dz+2*orho3*U4*dU1_dz;
        dw_dz_dU4=-orho2*dU1_dz;
        ds12_dU1=FourThird*du_dx_dU1-TwoThrid*dv_dy_dU1-TwoThrid*dw_dz_dU1;
        ds12_dU2=FourThird*du_dx_dU2;
        ds12_dU3=-TwoThrid*dv_dy_dU3;
        ds12_dU4=-TwoThrid*dw_dz_dU4;
        ds13_dU1=du_dy_dU1+dv_dx_dU1;
        ds13_dU2=du_dy_dU2;
        ds13_dU3=dv_dx_dU3;
        ds14_dU1=dw_dx_dU1+du_dz_dU1;
        ds14_dU2=du_dz_dU2;
        ds14_dU4=dw_dx_dU4;
        ds22_dU1=du_dy_dU1+dv_dx_dU1;
        ds22_dU2=du_dy_dU2;
        ds22_dU3=dv_dx_dU3;
        ds23_dU1=FourThird*dv_dy_dU1-TwoThrid*du_dx_dU1-TwoThrid*dw_dz_dU1;
        ds23_dU2=-TwoThrid*du_dx_dU2;
        ds23_dU3=FourThird*dv_dy_dU3;
        ds23_dU4=-TwoThrid*dw_dz_dU4;
        ds24_dU1=dv_dz_dU1+dw_dy_dU1;
        ds24_dU3=dv_dz_dU3;
        ds24_dU4=dw_dy_dU4;
        ds32_dU1=dw_dx_dU1+du_dz_dU1;
        ds32_dU2=du_dz_dU2;
        ds32_dU4=dw_dx_dU4;
        ds33_dU1=dv_dz_dU1+dw_dy_dU1;
        ds33_dU3=dv_dz_dU3;
        ds33_dU4=dw_dy_dU4;
        ds34_dU1=FourThird*dw_dz_dU1-TwoThrid*du_dx_dU1-TwoThrid*dv_dy_dU1;
        ds34_dU2=-TwoThrid*du_dx_dU2;
        ds34_dU3=-TwoThrid*dv_dy_dU3;
        ds34_dU4=FourThird*dw_dz_dU4;
        dsnx_dU1=ds12_dU1*nx+ds22_dU1*ny+ds32_dU1*nz;
        dsnx_dU2=ds12_dU2*nx+ds22_dU2*ny+ds32_dU2*nz;
        dsnx_dU3=ds12_dU3*nx+ds22_dU3*ny;
        dsnx_dU4=ds12_dU4*nx+ds32_dU4*nz;
        dsny_dU1=ds13_dU1*nx+ds23_dU1*ny+ds33_dU1*nz;
        dsny_dU2=ds13_dU2*nx+ds23_dU2*ny;
        dsny_dU3=ds13_dU3*nx+ds23_dU3*ny+ds33_dU3*nz;
        dsny_dU4=ds23_dU4*ny+ds33_dU4*nz;
        dsnz_dU1=ds14_dU1*nx+ds24_dU1*ny+ds34_dU1*nz;
        dsnz_dU2=ds14_dU2*nx+ds34_dU2*nz;
        dsnz_dU3=ds24_dU3*ny+ds34_dU3*nz;
        //? why there is value if 2D
        dsnz_dU4=ds14_dU4*nx+ds24_dU4*ny+ds34_dU4*nz;
        dsnv_dU1=u*dsnx_dU1+v*dsny_dU1+w*dsnz_dU1-orho2*U2*snx-orho2*U3*sny-orho2*U4*snz;
        dsnv_dU2=u*dsnx_dU2+v*dsny_dU2+w*dsnz_dU2+orho1*snx;
        dsnv_dU3=u*dsnx_dU3+v*dsny_dU3+w*dsnz_dU3+orho1*sny;
        dsnv_dU4=u*dsnx_dU4+v*dsny_dU4+w*dsnz_dU4+orho1*snz;
        tmpArray[0+0*nrow]=tmpArray[0+0*nrow]+mu*dsnx_dU1;
        tmpArray[0+1*nrow]=tmpArray[0+1*nrow]+mu*dsnx_dU2;
        tmpArray[0+2*nrow]=tmpArray[0+2*nrow]+mu*dsnx_dU3;
        tmpArray[0+3*nrow]=tmpArray[0+3*nrow]+mu*dsnx_dU4;
        tmpArray[1+0*nrow]=tmpArray[1+0*nrow]+mu*dsny_dU1;
        tmpArray[1+1*nrow]=tmpArray[1+1*nrow]+mu*dsny_dU2;
        tmpArray[1+2*nrow]=tmpArray[1+2*nrow]+mu*dsny_dU3;
        tmpArray[1+3*nrow]=tmpArray[1+3*nrow]+mu*dsny_dU4;
        tmpArray[2+0*nrow]=tmpArray[2+0*nrow]+mu*dsnz_dU1;
        tmpArray[2+1*nrow]=tmpArray[2+1*nrow]+mu*dsnz_dU2;
        tmpArray[2+2*nrow]=tmpArray[2+2*nrow]+mu*dsnz_dU3;
        tmpArray[2+3*nrow]=tmpArray[2+3*nrow]+mu*dsnz_dU4;
        tmpArray[3+0*nrow]=tmpArray[3+0*nrow]+mu*dsnv_dU1;
        tmpArray[3+1*nrow]=tmpArray[3+1*nrow]+mu*dsnv_dU2;
        tmpArray[3+2*nrow]=tmpArray[3+2*nrow]+mu*dsnv_dU3;
        tmpArray[3+3*nrow]=tmpArray[3+3*nrow]+mu*dsnv_dU4;

        //Consider heat flux qn's effect (does not include mu's effect)
        NekDouble dqx_dU1,dqx_dU2,dqx_dU3,dqx_dU4,dqx_dU5;
        NekDouble dqy_dU1,dqy_dU2,dqy_dU3,dqy_dU4,dqy_dU5;
        NekDouble dqz_dU1,dqz_dU2,dqz_dU3,dqz_dU4,dqz_dU5;
        NekDouble tmpx=-nx*tmp2;
        dqx_dU1=tmpx*(-orho2*dU5_dx+2*orho3*U5*dU1_dx+2*orho3*U2*dU2_dx-3*orho4*U2*U2*dU1_dx+2*orho3*U3*dU3_dx-3*orho4*U3*U3*dU1_dx+2*orho3*U4*dU4_dx-3*orho4*U4*U4*dU1_dx);
        dqx_dU2=tmpx*(-orho2*dU2_dx+2*orho3*U2*dU1_dx);
        dqx_dU3=tmpx*(-orho2*dU3_dx+2*orho3*U3*dU1_dx);
        dqx_dU4=tmpx*(-orho2*dU4_dx+2*orho3*U4*dU1_dx);
        dqx_dU5=-tmpx*orho2*dU1_dx;
        NekDouble tmpy=-ny*tmp2;
        dqy_dU1=tmpy*(-orho2*dU5_dy+2*orho3*U5*dU1_dy+2*orho3*U2*dU2_dy-3*orho4*U2*U2*dU1_dy+2*orho3*U3*dU3_dy-3*orho4*U3*U3*dU1_dy+2*orho3*U4*dU4_dy-3*orho4*U4*U4*dU1_dy);
        dqy_dU2=tmpy*(-orho2*dU2_dy+2*orho3*U2*dU1_dy);
        dqy_dU3=tmpy*(-orho2*dU3_dy+2*orho3*U3*dU1_dy);
        dqy_dU4=tmpy*(-orho2*dU4_dy+2*orho3*U4*dU1_dy);
        dqy_dU5=-tmpy*orho2*dU1_dy;
        NekDouble tmpz=-nz*tmp2;
        dqz_dU1=tmpz*(-orho2*dU5_dz+2*orho3*U5*dU1_dz+2*orho3*U2*dU2_dz-3*orho4*U2*U2*dU1_dz+2*orho3*U3*dU3_dz-3*orho4*U3*U3*dU1_dz+2*orho3*U4*dU4_dz-3*orho4*U4*U4*dU1_dz);
        dqz_dU2=tmpz*(-orho2*dU2_dz+2*orho3*U2*dU1_dz);
        dqz_dU3=tmpz*(-orho2*dU3_dz+2*orho3*U3*dU1_dz);
        dqz_dU4=tmpz*(-orho2*dU4_dz+2*orho3*U4*dU1_dz);
        dqz_dU5=-tmpz*orho2*dU1_dz;
        tmpArray[3+0*nrow]=tmpArray[3+0*nrow]-dqx_dU1-dqy_dU1-dqz_dU1;
        tmpArray[3+1*nrow]=tmpArray[3+1*nrow]-dqx_dU2-dqy_dU2-dqz_dU2;
        tmpArray[3+2*nrow]=tmpArray[3+2*nrow]-dqx_dU3-dqy_dU3-dqz_dU3;
        tmpArray[3+3*nrow]=tmpArray[3+3*nrow]-dqx_dU4-dqy_dU4-dqz_dU4;
        tmpArray[3+4*nrow]=tmpArray[3+4*nrow]-dqx_dU5-dqy_dU5-dqz_dU5;
    }

    void NavierStokesCFE::v_MinusDiffusionFluxJacDirctn(
        const int                                                       nDirctn,
        const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();
        int nPts            = m_fields[0]->GetTotPoints();
        int nSpaceDim           = m_graph->GetSpaceDimension();
        Array<OneD, NekDouble> normals;
        Array<OneD, Array<OneD, NekDouble> > normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        normals =   normal3D[nDirctn];

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, m_mu);
        Array<OneD, NekDouble > DmuDT              (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (nPts, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
            m_varConv->GetDmuDT(temperature,mu,DmuDT);
        }
        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            Vmath::Vadd(nPts, mu, 1, m_muav, 1, mu, 1);
            // Get numerical DmuDT
        }

        // What about thermal conductivity?


        NekDouble pointmu       = 0.0;
        NekDouble pointDmuDT    = 0.0;
        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> locDmuDT;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > pointDerv(nSpaceDim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > locDerv(nSpaceDim);
        for(int j = 0; j < nSpaceDim; j++)
        {
            pointDerv[j] = Array<OneD, NekDouble>(nConvectiveFields,0.0);
            locDerv[j]   = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
        }

        DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
                                ::AllocateSharedPtr(nConvectiveFields-1, nConvectiveFields,0.0);
        Array<OneD, NekDouble > PointFJac_data = PointFJac->GetPtr();

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
            int noffest             = GetPhys_Offset(nelmt);

            for(int j = 0; j < nConvectiveFields; j++)
            {
                locVars[j] = inarray[j]+noffest;
            }

            for(int j = 0; j < nSpaceDim; j++)
            {
                for(int k = 0; k < nConvectiveFields; k++)
                {
                    locDerv[j][k] = qfields[j][k]+noffest;
                }
            locmu       =   mu      + noffest;
            locDmuDT    =   DmuDT   + noffest;
            for(int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                for(int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }
                for(int j = 0; j < nSpaceDim; j++)
                {
                    for(int k = 0; k < nConvectiveFields; k++)
                    {
                        pointDerv[j][k] = locDerv[j][k][npnt];
                    }
                }

                pointmu     = locmu[npnt];
                pointDmuDT  = locDmuDT[npnt];

                GetDiffusionFluxJacPoint(pointVar,pointDerv,pointmu,pointDmuDT,normals,PointFJac);
                for (int j =0; j < nConvectiveFields; j++)
                {
                    int noffset = j*(nConvectiveFields-1);
                    for (int i =0; i < nConvectiveFields-1; i++)
                    {
                        ElmtJacArray[i+1][j][nDirctn][nelmt][npnt] -= PointFJac_data[noffset+i];
                    }
                }
            }
        }
    }


    void NavierStokesCFE::v_MinusDiffusionFluxJacDirctnElmt(
        const int                                                       nConvectiveFields,
        const int                                                       nElmtPnt,
        const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
        const Array<OneD, Array<OneD,  Array<OneD, NekDouble> > >       &locDerv,
        const Array<OneD, NekDouble>                                    &locmu,
        const Array<OneD, NekDouble>                                    &locDmuDT,
        const Array<OneD, NekDouble>                                    &normals,
        DNekMatSharedPtr                                                &wspMat,
        Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray)
    {
        int nSpaceDim           = m_graph->GetSpaceDimension();  

        NekDouble pointmu       = 0.0;
        NekDouble pointDmuDT    = 0.0;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > pointDerv(nSpaceDim);
        for(int j = 0; j < nSpaceDim; j++)
        {   
            pointDerv[j] = Array<OneD, NekDouble>(nConvectiveFields,0.0);
        }

        Array<OneD, NekDouble > wspMatData = wspMat->GetPtr();

        for(int npnt = 0; npnt < nElmtPnt; npnt++)
        {
            for(int j = 0; j < nConvectiveFields; j++)
            {
                pointVar[j] = locVars[j][npnt];
            }
            for(int j = 0; j < nSpaceDim; j++)
            {   
                for(int k = 0; k < nConvectiveFields; k++)
                {
                    pointDerv[j][k] = locDerv[j][k][npnt];
                }
            }

            pointmu     = locmu[npnt];
            pointDmuDT  = locDmuDT[npnt];

            GetDiffusionFluxJacPoint(pointVar,pointDerv,pointmu,pointDmuDT,normals,wspMat);
            for (int j =0; j < nConvectiveFields; j++)
            {
                int noffset = j*nConvectiveFields;

                Vmath::Vsub(nConvectiveFields-1,&PntJacArray[npnt][noffset+1],1,
                                                &wspMatData[noffset-j],1,
                                                &PntJacArray[npnt][noffset+1],1);
            }
        }
    }

    /**
     * @brief Return the penalty vector for the LDGNS diffusion problem.
     */
    void NavierStokesCFE::v_GetFluxPenalty(
        const Array<OneD, Array<OneD, NekDouble> >  &uFwd,
        const Array<OneD, Array<OneD, NekDouble> >  &uBwd,
              Array<OneD, Array<OneD, NekDouble> >  &penaltyCoeff)
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

    void NavierStokesCFE::v_GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals,
                 DNekMatSharedPtr                               &fluxJac)
    {
        switch (m_spacedim)
        {
        case 2:
            GetdFlux_dU_2D(normals,mu,DmuDT,conservVar,conseDeriv,fluxJac);
            break;

        case 3:
            GetdFlux_dU_3D(normals,mu,DmuDT,conservVar,conseDeriv,fluxJac);
            break;

        default:
            ASSERTL0(false, "v_GetDiffusionFluxJacPoint not coded");
            break;
        }
    }

    void NavierStokesCFE::v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr                            &explist,
        const Array<OneD, const Array<OneD, NekDouble> >                &normals,
        const int                                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray,
        const int                                                       nfluxDir)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    explist->GetExp();
        int ntotElmt            = (*expvect).size();
        int nPts                = explist->GetTotPoints();
        int nSpaceDim           = m_graph->GetSpaceDimension();

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (nPts, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
        }

        NekDouble pointmu       = 0.0;
        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        Array<OneD, NekDouble> pointnormals(nSpaceDim,0.0);
        Array<OneD, Array<OneD, NekDouble> > locnormal(nSpaceDim);

        DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
                                ::AllocateSharedPtr(nConvectiveFields-1, nConvectiveFields);
        Array<OneD, NekDouble > PointFJac_data = PointFJac->GetPtr();

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
            int noffest             = explist->GetPhys_Offset(nelmt);

            for(int j = 0; j < nConvectiveFields; j++)
            {
                locVars[j] = inarray[j]+noffest;
            }

            for(int j = 0; j < nSpaceDim; j++)
            {
                locnormal[j] = normals[j]+noffest;
            }

            locmu       =   mu      + noffest;
            for(int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                for(int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }
                for(int j = 0; j < nSpaceDim; j++)
                {
                    pointnormals[j] = locnormal[j][npnt];
                }

                pointmu     = locmu[npnt];

                // GetdFlux_dQx_2D(pointnormals,pointmu,pointVar,PointFJac);
                // functor(pointnormals,pointmu,pointVar,PointFJac);
                m_GetdFlux_dDeriv_Array[nDervDir](pointnormals,pointmu,pointVar,PointFJac);
                for (int j =0; j < nConvectiveFields; j++)
                {
                    // (*ElmtJac[nelmt][npnt])(0,j) =  0.0;
                    ElmtJacArray[0][j][nfluxDir][nelmt][npnt] = 0.0;
                }
                for (int j =0; j < nConvectiveFields; j++)
                {
                    int noffset = j*(nConvectiveFields-1);
                    for (int i =0; i < nConvectiveFields-1; i++)
                    {
                        ElmtJacArray[i+1][j][nfluxDir][nelmt][npnt] = PointFJac_data[noffset+i];
                    }
                }
            }
        }
    }

    void NavierStokesCFE::v_GetFluxDerivJacDirctnElmt(
        const int                                                       nConvectiveFields,
        const int                                                       nElmtPnt,
        const int                                                       nDervDir,
        const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
        const Array<OneD, NekDouble>                                    &locmu,
        const Array<OneD, Array<OneD, NekDouble> >                      &locnormal,
        DNekMatSharedPtr                                                &wspMat,
        Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray)
    {
        int nSpaceDim           = m_graph->GetSpaceDimension();  
        
        NekDouble pointmu       = 0.0;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, NekDouble> pointnormals(nSpaceDim,0.0);

        Array<OneD, NekDouble > wspMatData = wspMat->GetPtr();
                
        for(int npnt = 0; npnt < nElmtPnt; npnt++)
        {
            for(int j = 0; j < nConvectiveFields; j++)
            {
                pointVar[j] = locVars[j][npnt];
            }
            for(int j = 0; j < nSpaceDim; j++)
            {   
                pointnormals[j] = locnormal[j][npnt];
            }

            pointmu     = locmu[npnt];

            m_GetdFlux_dDeriv_Array[nDervDir](pointnormals,pointmu,pointVar,wspMat);
            Vmath::Zero(nConvectiveFields,&PntJacArray[npnt][0],nConvectiveFields);
            for (int j =0; j < nConvectiveFields; j++)
            {
                int noffset = j*(nConvectiveFields-1);
                Vmath::Vcopy((nConvectiveFields-1),&wspMatData[noffset],1,&PntJacArray[npnt][noffset+j+1],1);
            }
        }
    }
    
    void NavierStokesCFE::v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr                &explist,
        const Array<OneD, const Array<OneD, NekDouble> >    &normals,
        const int                                           nDervDir,
        const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
              Array<OneD, Array<OneD, DNekMatSharedPtr> >   &ElmtJac)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =
            explist->GetExp();
        int ntotElmt            = (*expvect).size();
        int nPts                = explist->GetTotPoints();
        int nSpaceDim           = m_graph->GetSpaceDimension();

        //Debug
        if(!ElmtJac.num_elements())
        {
            ElmtJac =   Array<OneD, Array<OneD, DNekMatSharedPtr> > (ntotElmt);
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
                ElmtJac[nelmt] =   Array<OneD, DNekMatSharedPtr>(nElmtPnt);
                for(int npnt = 0; npnt < nElmtPnt; npnt++)
                {
                    ElmtJac[nelmt][npnt] = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(nConvectiveFields, nConvectiveFields);
                }
            }
        }
        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, m_mu);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (nPts, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
        }

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            Array<OneD, NekDouble> muav;
            if (m_fields[0]->GetTrace()->GetTotPoints()==nPts)
            {
                muav = m_muavTrace;
            }
            else
            {
                muav = m_muav;
            }
            Vmath::Vadd(nPts, mu, 1, muav, 1, mu, 1);
        }

        // What about thermal conductivity?

        NekDouble pointmu       = 0.0;
        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        Array<OneD, NekDouble> pointnormals(nSpaceDim,0.0);
        Array<OneD, Array<OneD, NekDouble> > locnormal(nSpaceDim);

        DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
                                ::AllocateSharedPtr(nConvectiveFields-1, nConvectiveFields);
        Array<OneD, NekDouble > tmpMatinnData, tmpMatoutData;
        
        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
            int noffest             = explist->GetPhys_Offset(nelmt);
            
            for(int j = 0; j < nConvectiveFields; j++)
            {
                locVars[j] = inarray[j]+noffest;
            }
            
            for(int j = 0; j < nSpaceDim; j++)
            {
                locnormal[j] = normals[j]+noffest;
            }
            
            locmu       =   mu      + noffest;
            for(int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                for(int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }
                for(int j = 0; j < nSpaceDim; j++)
                {
                    pointnormals[j] = locnormal[j][npnt];
                        }
                
                pointmu     = locmu[npnt];
                
                // GetdFlux_dQx_2D(pointnormals,pointmu,pointVar,PointFJac);
                // functor(pointnormals,pointmu,pointVar,PointFJac);
                m_GetdFlux_dDeriv_Array[nDervDir](pointnormals,pointmu,pointVar,PointFJac);
                tmpMatinnData = PointFJac->GetPtr();
                tmpMatoutData = ElmtJac[nelmt][npnt]->GetPtr();
                
                Vmath::Fill(nConvectiveFields,0.0,&tmpMatoutData[0],
                            nConvectiveFields);
                for (int j =0; j < nConvectiveFields; j++)
                {
                    Vmath::Vcopy(nConvectiveFields-1,
                                 &tmpMatinnData[j*(nConvectiveFields-1)],1,
                                 &tmpMatoutData[1+j*nConvectiveFields],1);
                }
                // (*ElmtJac[nelmt][npnt]) =   (*PointFJac);
            }
        }
    }
    
    void NavierStokesCFE::v_CalphysDeriv(
        const Array<OneD, const Array<OneD, NekDouble> >          &inarray,
        Array<OneD,       Array<OneD, Array<OneD, NekDouble> > >  &qfield)
    {
        int nConvectiveFields = m_fields.num_elements();
        int npoints           = GetTotPoints();
        const Array<OneD, Array<OneD, NekDouble> >                  pFwd;
        const Array<OneD, Array<OneD, NekDouble> >                  pBwd;
        if(!qfield.num_elements())
        {
            qfield  =   Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(m_spacedim);
            for(int i = 0; i< m_spacedim; i++)
            {
                qfield[i]   =   Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for(int j = 0; j< nConvectiveFields; j++)
                {
                    qfield[i][j]   =   Array<OneD, NekDouble>(npoints,0.0);
                }
            }
        }
        m_diffusion->DiffuseCalculateDerivative(nConvectiveFields,m_fields,inarray,qfield,pFwd,pBwd);
    }

    void NavierStokesCFE::v_CalcMuDmuDT(
        const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
        Array<OneD, NekDouble>                                          &mu,
        Array<OneD, NekDouble>                                          &DmuDT)
    {
        int npoints = mu.num_elements();
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (npoints, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
            if(DmuDT.num_elements()>0)
            {
                m_varConv->GetDmuDT(temperature,mu,DmuDT);
            }
        }
        else
        {
            Vmath::Fill(npoints, m_mu, mu, 1);
            if(DmuDT.num_elements()>0)
            {
                Vmath::Zero(npoints, DmuDT, 1);
            }
        }
    }
}
