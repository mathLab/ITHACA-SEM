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
}

void NonSmoothShockCapture::v_GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu)
{
    const int nElements = m_fields[0]->GetExpSize();
    int PointCount      = 0;
    int nTotQuadPoints  = m_fields[0]->GetTotPoints();

    Array<OneD, NekDouble> Sensor     (nTotQuadPoints, 0.0);
    Array<OneD, NekDouble> SensorKappa(nTotQuadPoints, 0.0);

    m_varConv->GetSensor(m_fields[0], physfield, Sensor, SensorKappa);

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
                mu[n+PointCount] = 0;
            }
            else if (Sensor[n+PointCount] >= (m_Skappa-m_Kappa)
                    && Sensor[n+PointCount] <= (m_Skappa+m_Kappa))
            {
                mu[n+PointCount] = mu_0 * (0.5 * (1 + sin(
                                       M_PI * (Sensor[n+PointCount] -
                                               m_Skappa - m_Kappa) /
                                                        (2*m_Kappa))));
            }
            else if (Sensor[n+PointCount] > (m_Skappa+m_Kappa))
            {
                mu[n+PointCount] = mu_0;
            }
        }

        PointCount += nQuadPointsElement;
    }
}

}
