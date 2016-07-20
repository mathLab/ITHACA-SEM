///////////////////////////////////////////////////////////////////////////////
//
// File EulerArtificialDiffusionCFE.cpp
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
// Description: Euler equations in conservative variables with artificial
// diffusion
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/EulerADCFE.h>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
    string EulerADCFE::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "EulerADCFE", EulerADCFE::create,
        "Euler equations in conservative variables with "
        "artificial diffusion.");

    EulerADCFE::EulerADCFE(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : CompressibleFlowSystem(pSession)
    {
    }

    void EulerADCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        m_session->LoadParameter ("FL",            m_FacL,          0.0);
        m_session->LoadParameter ("FH",            m_FacH,          0.0);
        m_session->LoadParameter ("hFactor",       m_hFactor,       1.0);
        m_session->LoadParameter ("C1",            m_C1,            3.0);
        m_session->LoadParameter ("C2",            m_C2,            5.0);
        m_session->LoadSolverInfo("ShockCaptureType",
                                  m_shockCaptureType,    "Off");

        if (m_shockCaptureType == "Smooth")
        {
            m_diffusion->SetArtificialDiffusionVector(
                        &EulerADCFE::GetSmoothArtificialViscosity, this);

            ASSERTL0(m_fields.num_elements() == m_spacedim + 3,
                     "Not enough variables for smooth shock capturing; "
                     "make sure you have added eps to variable list.");
            m_smoothDiffusion = true;
        }
        if (m_shockCaptureType=="NonSmooth")
        {
            m_diffusion->SetArtificialDiffusionVector(
                        &EulerADCFE::GetArtificialDynamicViscosity, this);
        }
    }

    EulerADCFE::~EulerADCFE()
    {

    }

    void EulerADCFE::v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD,       Array<OneD, NekDouble> > &outarray)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

        for (i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }

        m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff);

        if (m_shockCaptureType == "NonSmooth")
        {
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarray[i], 1,
                            outarrayDiff[i], 1,
                            outarray[i], 1);
            }
        }
        if(m_shockCaptureType == "Smooth")
        {
            const Array<OneD, int> ExpOrder = GetNumExpModesPerExp();

            NekDouble pOrder = Vmath::Vmax(ExpOrder.num_elements(), ExpOrder, 1);

            Array <OneD, NekDouble > a_vel  (npoints, 0.0);
            Array <OneD, NekDouble > u_abs  (npoints, 0.0);
            Array <OneD, NekDouble > pres   (npoints, 0.0);
            Array <OneD, NekDouble > wave_sp(npoints, 0.0);

            m_varConv->GetPressure(inarray, pres);
            m_varConv->GetSoundSpeed(inarray, pres, a_vel);
            m_varConv->GetAbsoluteVelocity(inarray, u_abs);

            Vmath::Vadd(npoints, a_vel, 1, u_abs, 1, wave_sp, 1);

            NekDouble max_wave_sp = Vmath::Vmax(npoints, wave_sp, 1);

            Vmath::Smul(npoints,
                        m_C2,
                        outarrayDiff[nvariables-1], 1,
                        outarrayDiff[nvariables-1], 1);

            Vmath::Smul(npoints,
                        max_wave_sp,
                        outarrayDiff[nvariables-1], 1,
                        outarrayDiff[nvariables-1], 1);

            Vmath::Smul(npoints,
                        pOrder,
                        outarrayDiff[nvariables-1], 1,
                        outarrayDiff[nvariables-1], 1);

            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarray[i], 1,
                            outarrayDiff[i], 1,
                            outarray[i], 1);
            }

            Array<OneD, Array<OneD, NekDouble> > outarrayForcing(nvariables);

            for (i = 0; i < nvariables; ++i)
            {
                outarrayForcing[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }

            GetForcingTerm(inarray, outarrayForcing);

            for (i = 0; i < nvariables; ++i)
            {
                // Add Forcing Term
                Vmath::Vadd(npoints,
                            outarray[i], 1,
                            outarrayForcing[i], 1,
                            outarray[i], 1);
            }
        }
    }

    void EulerADCFE::GetForcingTerm(
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
        Array <OneD, NekDouble > absVelocity(nPts, 0.0);

        Array<OneD,int> pOrderElmt = GetNumExpModesPerExp();
        Array<OneD, NekDouble> pOrder (nPts, 0.0);

        // Thermodynamic related quantities
        m_varConv->GetPressure(inarray, pressure);
        m_varConv->GetSoundSpeed(inarray, pressure, soundspeed);
        m_varConv->GetAbsoluteVelocity(inarray, absVelocity);
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

    void EulerADCFE::GetSmoothArtificialViscosity(
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
        m_varConv->GetPressure(physfield, pressure);
        m_varConv->GetTemperature(physfield, pressure, temperature);
        m_varConv->GetSoundSpeed(physfield, pressure, soundspeed);
        m_varConv->GetAbsoluteVelocity(physfield, absVelocity);
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

    void EulerADCFE::GetArtificialDynamicViscosity(
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

        m_varConv->GetAbsoluteVelocity(physfield, absVelocity);
        m_varConv->GetPressure        (physfield, pressure);
        m_varConv->GetSoundSpeed      (physfield, pressure, soundspeed);
        GetSensor          (physfield, Sensor, SensorKappa);

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

    void EulerADCFE::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        CompressibleFlowSystem::v_ExtraFldOutput(fieldcoeffs, variables);

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

            Array<OneD, NekDouble> smooth(nPhys);
            GetSmoothArtificialViscosity    (tmp, smooth);

            Array<OneD, NekDouble> smoothFwd(nCoeffs);
            m_fields[0]->FwdTrans(smooth, smoothFwd);

            variables.push_back  ("SmoothVisc");
            fieldcoeffs.push_back(smoothFwd);
        }
    }

}
