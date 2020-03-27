///////////////////////////////////////////////////////////////////////////////
//
// File: NonSmoothShockCapture.cpp
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
// Description: NonSmooth artificial diffusion for shock capture
//
///////////////////////////////////////////////////////////////////////////////

#include "NonSmoothShockCapture.h"

using namespace std;

namespace Nektar
{

std::string NonSmoothShockCapture::className = GetArtificialDiffusionFactory().
    RegisterCreatorFunction("NonSmooth",
                        NonSmoothShockCapture::create,
                        "NonSmooth artificial diffusion for shock capture.");

NonSmoothShockCapture::NonSmoothShockCapture(
           const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const int spacedim)
    : ArtificialDiffusion(pSession, pFields, spacedim)
{
    m_session->LoadParameter ("SensorOffset",  m_offset, 1);
}

void NonSmoothShockCapture::v_GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu)
{
    int nPts = m_fields[0]->GetTotPoints();
    int nElements = m_fields[0]->GetExpSize();

    // Determine the maximum wavespeed
    Array <OneD, NekDouble > Lambda(nPts, 0.0);
    Array <OneD, NekDouble > soundspeed(nPts, 0.0);
    Array <OneD, NekDouble > absVelocity(nPts, 0.0);
    m_varConv->GetSoundSpeed(physfield, soundspeed);
    m_varConv->GetAbsoluteVelocity(physfield, absVelocity);
    Vmath::Vadd(nPts, absVelocity, 1, soundspeed, 1, Lambda, 1);

    // Compute sensor based on rho
    Array<OneD, NekDouble> Sensor(nPts, 0.0);
    m_varConv->GetSensor(m_fields[0], physfield, Sensor, mu, m_offset);

    // Scale AV
    Array<OneD, NekDouble> tmp;
    for (int e = 0; e < nElements; e++)
    {
        int physOffset  = m_fields[0]->GetPhys_Offset(e);
        int nElmtPoints = m_fields[0]->GetExp(e)->GetTotPoints();

        // Get max wavespeed per element
        NekDouble LambdaElmt = 0.0;
        LambdaElmt = Vmath::Vmax(nElmtPoints, tmp = Lambda + physOffset, 1);

        // Scale viscosity by the maximum wave speed and h/p
        LambdaElmt *= m_mu0 * m_hOverP[e];
        Vmath::Smul(nElmtPoints, LambdaElmt, tmp = mu + physOffset, 1,
            tmp = mu + physOffset, 1);
    }
}
}