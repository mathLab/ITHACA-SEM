///////////////////////////////////////////////////////////////////////////////
//
// File: SmoothShockCapture.cpp
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
// Description: Smooth artificial diffusion for shock capture
//
///////////////////////////////////////////////////////////////////////////////

#include "SmoothShockCapture.h"

using namespace std;

namespace Nektar
{

std::string SmoothShockCapture::className = GetArtificialDiffusionFactory().
    RegisterCreatorFunction("Smooth",
                            SmoothShockCapture::create,
                            "Smooth artificial diffusion for shock capture.");

SmoothShockCapture::SmoothShockCapture(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const int spacedim)
    : ArtificialDiffusion(pSession, pFields, spacedim)
{
    ASSERTL0(m_fields.num_elements() == spacedim + 3,
             "Not enough variables for smooth shock capturing; "
             "make sure you have added eps to variable list.");
}

void SmoothShockCapture::v_DoArtificialDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    int i;
    int nvariables = inarray.num_elements();
    int npoints    = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

    for (i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }

    m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff);

    const Array<OneD, int> ExpOrder = m_fields[0]->EvalBasisNumModesMaxPerExp();

    NekDouble pOrder = Vmath::Vmax(ExpOrder.num_elements(), ExpOrder, 1);

    Array <OneD, NekDouble > a_vel  (npoints, 0.0);
    Array <OneD, NekDouble > u_abs  (npoints, 0.0);
    Array <OneD, NekDouble > tmp   (npoints, 0.0);

    m_varConv->GetPressure(inarray, tmp);
    m_varConv->GetSoundSpeed(inarray, tmp, a_vel);
    m_varConv->GetAbsoluteVelocity(inarray, u_abs);

    Vmath::Vadd(npoints, a_vel, 1, u_abs, 1, tmp, 1);

    NekDouble max_wave_sp = Vmath::Vmax(npoints, tmp, 1);

    NekDouble fac = m_C2*max_wave_sp*pOrder;

    Vmath::Smul(npoints,
                fac,
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

void SmoothShockCapture::v_GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu)
{
    int nvariables = physfield.num_elements();
    int nPts       = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble > sensor     (nPts, 0.0);
    Array<OneD, NekDouble > SensorKappa(nPts, 0.0);

    // Calculate sensor
    m_varConv->GetSensor(m_fields[0], physfield, sensor, SensorKappa);

    NekDouble ThetaH = m_FacH;
    NekDouble ThetaL = m_FacL;

    NekDouble Phi0     = (ThetaH+ThetaL)/2;
    NekDouble DeltaPhi = ThetaH-Phi0;

    for (int e = 0; e < mu.num_elements(); e++)
    {
        if (physfield[nvariables-1][e] <= (Phi0 - DeltaPhi))
        {
            mu[e] = 0;
        }
        else if(physfield[nvariables-1][e] >= (Phi0 + DeltaPhi))
        {
            mu[e] = m_mu0;
        }
        else if(abs(physfield[nvariables-1][e]-Phi0) < DeltaPhi)
        {
            mu[e] = m_mu0/2*(1+sin(M_PI*
                (physfield[nvariables-1][e]-Phi0)/(2*DeltaPhi)));
        }
    }
}

void SmoothShockCapture::GetForcingTerm(
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

    Array<OneD,int> pOrderElmt = m_fields[0]->EvalBasisNumModesMaxPerExp();
    Array<OneD, NekDouble> pOrder (nPts, 0.0);

    // Thermodynamic related quantities
    m_varConv->GetPressure(inarray, pressure);
    m_varConv->GetSoundSpeed(inarray, pressure, soundspeed);
    m_varConv->GetAbsoluteVelocity(inarray, absVelocity);
    m_varConv->GetSensor(m_fields[0], inarray, Sensor, SensorKappa);

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

}
