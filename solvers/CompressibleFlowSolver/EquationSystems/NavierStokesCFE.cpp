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

        string diffName;
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

            m_diffusion->SetSpecialBndTreat(
                &NavierStokesCFE::SpecialBndTreat, this);

            m_diffusion->SetDiffusionSymmFluxCons(
                &NavierStokesCFE::GetViscousSymmtrFluxConservVar, this);

            if (m_shockCaptureType != "Off")
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

            if (m_shockCaptureType != "Off")
            {
                m_diffusion->SetArtificialDiffusionVector(
                    &NavierStokesCFE::GetArtificialViscosity, this);
            }
        }

        // Set up penalty term for LDGNS
        if ("LDGNS" == diffName||
            "LDGNS3DHomogeneous1D" == diffName)
        {
            m_diffusion->SetFluxPenaltyNS(&NavierStokesCFE::
                v_GetFluxPenalty, this);
        }

        // Concluding initialisation of diffusion operator
        m_diffusion->InitObject         (m_session, m_fields);
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
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Update viscosity and thermal conductivity
        GetViscosityAndThermalCondFromTemp(physfield[nScalar-1], m_mu,
            m_thermalConductivity);

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
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Update viscosity and thermal conductivity
        GetViscosityAndThermalCondFromTemp(physfield[nScalar-1], m_mu,
            m_thermalConductivity);

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
        size_t nConvectiveFields = inarray.size();
        size_t nPts = inarray[0].size();
        int n_nonZero = nConvectiveFields - 1;

        TensorOfArray3D<NekDouble> fluxVec;
        Array<OneD, Array<OneD, NekDouble>> outtmp{nConvectiveFields};

        for (int i = 0; i < nConvectiveFields; ++i)
        {
            outtmp[i] = Array<OneD, NekDouble>{nPts, 0.0};
        }

        // is this zeroing necessary?
        for (int i = 0; i < outarray.size(); ++i)
        {
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                Vmath::Zero(nPts, outarray[i][j], 1);
            }
        }

        // unless moved to the npts loop these can be pre-allocated
        Array<OneD, NekDouble > temperature        {nPts, 0.0};
        Array<OneD, NekDouble > mu                 {nPts, 0.0};
        Array<OneD, NekDouble > thermalConductivity        {nPts, 0.0};
        m_varConv->GetTemperature(inarray, temperature);
        GetViscosityAndThermalCondFromTemp(temperature, mu,
                                            thermalConductivity);


        if (normal.size())
        {
            for (int nd = 0; nd < nDim; ++nd)
            {
                for (int nderiv = 0; nderiv < nDim; ++nderiv)
                {
                    GetViscousFluxBilinearForm(nDim, nd, nderiv, inarray,
                                                qfields[nderiv], mu, outtmp);

                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        Vmath::Vvtvp(nPts, normal[nd], 1, outtmp[j], 1,
                            outarray[0][j], 1, outarray[0][j], 1);
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
                                                qfields[nderiv], mu, outtmp);

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

                    Vmath::Vvtvp(nBndEdgePts, bndIntEndy, 1, bndCons[ndens], 1,
                        bndTotEngy, 1, bndTotEngy, 1);


                    Vmath::Vcopy(nBndEdgePts,
                                &bndTotEngy[0], 1,
                                &consvar[nengy][id2], 1);
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

        // unless moved to the npts loop these can be pre-allocated
        Array<OneD, NekDouble > temperature        {nPts, 0.0};
        Array<OneD, NekDouble > mu                 {nPts, 0.0};
        Array<OneD, NekDouble > thermalConductivity        {nPts, 0.0};
        m_varConv->GetTemperature(inaverg, temperature);
        GetViscosityAndThermalCondFromTemp(temperature, mu,
                                            thermalConductivity);


        Array<OneD, Array<OneD, NekDouble> > outtmp{nConvectiveFields};
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            outtmp[i] =   Array<OneD, NekDouble> {nPts, 0.0};
        }
        for (int nd = 0; nd < nSpaceDim; ++nd)
        {
            for (int nderiv = 0; nderiv < nSpaceDim; ++nderiv)
            {
                GetViscousFluxBilinearForm(nSpaceDim, nd, nderiv, inaverg,
                                            inarray, mu, outtmp);

                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Vmath::Vvtvp(nPts, outtmp[i], 1, normals[nderiv], 1,
                        outarray[nd][i], 1, outarray[nd][i], 1);
                }
            }
        }
    }

    /**
     * @brief Calculate diffusion flux using the Jacobian form.
     *
     * @param
     * outarray [nvars][npts]
     */
    void NavierStokesCFE::GetViscousFluxBilinearForm(
        const int                                           nSpaceDim,
        const int                                           FluxDirection,
        const int                                           DerivDirection,
        const Array<OneD, const Array<OneD, NekDouble> >    &inaverg,
        const Array<OneD, const Array<OneD, NekDouble> >    &injumpp,
        const Array<OneD, NekDouble>                        &mu,
        Array<OneD, Array<OneD, NekDouble> >                &outarray)
    {
        size_t nConvectiveFields   = inaverg.size();
        size_t nPts = inaverg[nConvectiveFields - 1].size();
        size_t nDim = nSpaceDim;

        // Constants
        size_t nDim_plus_one = nDim + 1;
        size_t FluxDirection_plus_one = FluxDirection + 1;

        NekDouble gammaoPr = m_gamma / m_Prandtl;
        NekDouble one_minus_gammaoPr = 1.0 - gammaoPr;

        constexpr NekDouble OneThird = 1. / 3.;
        constexpr NekDouble TwoThird = 2. * OneThird;
        constexpr NekDouble FourThird = 4. * OneThird;

        size_t DerivDirection_plus_one = DerivDirection + 1;
        if (DerivDirection == FluxDirection)
        {
            // Loop over the points
            for (size_t p = 0; p < nPts; ++p)
            {
                // necessary???
                outarray[0][p] = 0.0; // store 1x

                // load 1/rho
                NekDouble oneOrho = 1.0 / inaverg[0][p]; // load 1x
                // get vel, vel^2, sum of vel^2
                std::array<NekDouble, 3> u{}, u2{};
                NekDouble u2sum{};
                for (size_t d = 0; d < nDim; ++d)
                {
                    u[d] = inaverg[d+1][p] * oneOrho; // load 1x
                    u2[d] = u[d]*u[d];
                    u2sum += u2[d];
                }

                // get E - sum v^2
                NekDouble E_minus_u2sum = inaverg[nDim_plus_one][p]; // load 1x
                E_minus_u2sum *= oneOrho;
                E_minus_u2sum -= u2sum;

                // get nu = mu/rho
                NekDouble nu = mu[p] * oneOrho; // load 1x


                // ^^^^ above is almost the same for both loops

                NekDouble tmp1 = OneThird * u2[FluxDirection] + u2sum;
                tmp1 += gammaoPr * E_minus_u2sum;
                tmp1 *= injumpp[0][p]; //load 1x

                NekDouble tmp2 = gammaoPr * injumpp[nDim_plus_one][p] - tmp1; //load 1x

                // local var for energy output
                NekDouble outTmpE = 0.0;
                for (size_t d = 0; d < nDim; ++d)
                {
                    size_t d_plus_one = d + 1;
                    //flux[rhou, rhov, rhow]
                    NekDouble outTmpD = injumpp[d_plus_one][p] - u[d] * injumpp[0][p];
                    outTmpD *= nu;
                    //flux rhoE
                    NekDouble tmp3 = one_minus_gammaoPr * u[d];
                    outTmpE +=  tmp3 * injumpp[d_plus_one][p];

                    if (d == FluxDirection)
                    {
                        outTmpD *= FourThird;
                        NekDouble tmp4 = OneThird * u[FluxDirection];
                        outTmpE += tmp4 * injumpp[FluxDirection_plus_one][p];
                    }

                    outarray[d_plus_one][p] = outTmpD; //store 1x
                }

                outTmpE += tmp2;
                outTmpE *= nu;
                outarray[nDim_plus_one][p] = outTmpE; //store 1x

            }
        }
        else
        {
            // Loop over the points
            for (size_t p = 0; p < nPts; ++p)
            {
                // Always zero
                outarray[0][p] = 0.0; // store 1x

                // load 1/rho
                NekDouble oneOrho = 1.0 / inaverg[0][p]; // load 1x
                // get vel, vel^2, sum of vel^2
                std::array<NekDouble, 3> u{}, u2{};
                NekDouble u2sum{};
                for (size_t d = 0; d < nDim; ++d)
                {
                    size_t d_plus_one = d + 1;
                    u[d] = inaverg[d_plus_one][p] * oneOrho; // load 1x
                    u2[d] = u[d]*u[d];
                    u2sum += u2[d];
                    // Not all directions are set
                    // one could work out the one that is not set
                    outarray[d_plus_one][p] = 0.0; // store 1x
                }

                // get E - sum v^2
                NekDouble E_minus_u2sum = inaverg[nDim_plus_one][p]; // load 1x
                E_minus_u2sum *= oneOrho;
                E_minus_u2sum -= u2sum;

                // get nu = mu/rho
                NekDouble nu = mu[p] * oneOrho; // load 1x


                // ^^^^ above is almost the same for both loops

                NekDouble tmp1 = u[DerivDirection] * injumpp[0][p] -
                    injumpp[DerivDirection_plus_one][p]; // load 2x
                tmp1 *= TwoThird;
                outarray[FluxDirection_plus_one][p] = nu * tmp1; // store 1x

                tmp1 = - u[FluxDirection] * injumpp[0][p] +
                    injumpp[FluxDirection_plus_one][p];
                outarray[DerivDirection_plus_one][p] = nu * tmp1; // store 1x

                tmp1 = OneThird * u[FluxDirection] * u[DerivDirection];
                tmp1 *= injumpp[0][p];

                NekDouble tmp2 = TwoThird * u[FluxDirection] *
                    injumpp[DerivDirection_plus_one][p];

                tmp1 += tmp2;

                tmp1 = u[DerivDirection] * injumpp[FluxDirection_plus_one][p] -
                    tmp1;
                outarray[nDim_plus_one][p] = nu * tmp1; // store 1x

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


    /**
     * @brief Update viscosity
     * todo: add artificial viscosity here
     */
    void NavierStokesCFE::GetViscosityAndThermalCondFromTemp(
        const Array<OneD, NekDouble> &temperature,
        Array<OneD, NekDouble> &mu,
        Array<OneD, NekDouble> &thermalCond)
    {
        int nPts       = temperature.size();

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(temperature, mu);
        }
        else
        {
            Vmath::Fill(nPts, m_muRef, mu, 1);
        }
        NekDouble tRa = m_Cp / m_Prandtl;
        Vmath::Smul(nPts, tRa, mu, 1, thermalCond, 1);

    }

}
