 ///////////////////////////////////////////////////////////////////////////////
//
// File BetaPressureArea.cpp
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
// Description: BetaPressureArea class
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/BetaPressureArea.h>

using namespace std;

namespace Nektar
{

std::string BetaPressureArea::className =
    GetPressureAreaFactory().RegisterCreatorFunction(
        "Beta", BetaPressureArea::create,
        "Beta law pressure area relationship for the arterial system");

BetaPressureArea::BetaPressureArea(
    Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
    const LibUtilities::SessionReaderSharedPtr pSession)
    : PulseWavePressureArea(pVessel, pSession)
{
}

BetaPressureArea::~BetaPressureArea()
{
}

void BetaPressureArea::v_GetPressure(NekDouble &P, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &dAUdx,
                                 const NekDouble &gamma, const NekDouble &alpha)
{
    P = m_PExt + beta * (sqrt(A) - sqrt(A0)) - gamma * dAUdx / sqrt(A); // Viscoelasticity
}

void BetaPressureArea::v_GetC(NekDouble &c, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &alpha)
{
    c = sqrt(beta / (2 * m_rho)) * sqrt(sqrt(A)); // Elastic
}

void BetaPressureArea::v_GetW1(NekDouble &W1, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble I = 0.0;
    GetCharIntegral(I, beta, A, A0);

    W1 = u + I; // Elastic and assumes u0 = 0
}

void BetaPressureArea::v_GetW2(NekDouble &W2, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble I = 0.0;
    GetCharIntegral(I, beta, A, A0);

    W2 = u - I; // Elastic and assumes u0 = 0
}

void BetaPressureArea::v_GetAFromChars(NekDouble &A, const NekDouble &W1,
                const NekDouble &W2, const NekDouble &beta, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    A = pow((W1 - W2) * sqrt(2 * m_rho / beta) / 8 + sqrt(sqrt(A0)), 4);
}

void BetaPressureArea::v_GetUFromChars(NekDouble &u, const NekDouble &W1,
                                                            const NekDouble &W2)
{
    u = (W1 + W2) / 2; // Necessarily the case for all tube laws
}

void BetaPressureArea::v_GetCharIntegral(NekDouble &I, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &alpha)
{
    NekDouble c  = 0.0;
    NekDouble c0 = 0.0;

    GetC(c,  beta, A,  A0);
    GetC(c0, beta, A0, A0);

    I = 4 * (c - c0);
}

void BetaPressureArea::v_GetJacobianInverse(NekMatrix<NekDouble> &invJ,
             const Array<OneD, NekDouble> &Au, const Array<OneD, NekDouble> &uu,
           const Array<OneD, NekDouble> &beta, const Array<OneD, NekDouble> &A0,
                   const Array<OneD, NekDouble> &alpha, const std::string &type)
{
    /*
    In the interest of speed, the inverse of the Jacobians for bifurcations,
    merges and junctions for the beta law have been calculated analytically.
    This can be done for other laws too, or the general formulation can be used
    instead.
    */

    NekDouble k = 0.0;

    if (type == "Bifurcation")
    {
        NekMatrix<NekDouble> J(6, 6);
        Array<OneD, NekDouble> c(3);
        Array<OneD, NekDouble> K(3);

        for (int i = 0; i < 3; ++i)
        {
            GetC(c[i], beta[i], Au[i], A0[i], alpha[i]);
        }

        k = c[0] * Au[1] * c[2] + Au[0] * c[2] * c[1] + Au[2] * c[0] * c[1];
        K[0] = (c[0] - uu[0]) * k;
        K[1] = (c[1] + uu[1]) * k;
        K[2] = (c[2] + uu[2]) * k;

        invJ.SetValue(0, 0, (-c[1] * uu[0] * c[2] * Au[0] + Au[2] * c[1] * c[0]
                                   * c[0] + Au[1] * c[0] * c[0] * c[2]) / K[0]);
        invJ.SetValue(0, 1, Au[1] * (c[1] - uu[1]) * c[0] * c[2] / K[0]);
        invJ.SetValue(0, 2, Au[2] * (c[2] - uu[2]) * c[0] * c[1] / K[0]);
        invJ.SetValue(0, 3, c[0] * c[1] * c[2] / K[0]);
        invJ.SetValue(0, 4, -0.5 * c[0] * Au[1] * c[2] / K[0]);
        invJ.SetValue(0, 5, -0.5 * Au[2] * c[0] * c[1] / K[0]);

        invJ.SetValue(1, 0, Au[0] * (c[0] + uu[0]) * c[1] * c[2] / K[1]);
        invJ.SetValue(1, 1, (c[0] * uu[1] * c[2] * Au[1] + Au[2] * c[0] * c[1]
                                   * c[1] + c[2] * c[1] * c[1] * Au[0]) / K[1]);
        invJ.SetValue(1, 2, -Au[2] * (c[2] - uu[2]) * c[0] * c[1] / K[1]);
        invJ.SetValue(1, 3, -c[0] * c[1] * c[2] / K[1]);
        invJ.SetValue(1, 4, -0.5 * (c[0] * Au[2] + Au[0] * c[2]) * c[1] / K[1]);
        invJ.SetValue(1, 5, 0.5 * Au[2] * c[0] * c[1] / K[1]);

        invJ.SetValue(2, 0, Au[0] * (c[0] + uu[0]) * c[1] * c[2] / K[2]);
        invJ.SetValue(2, 1, -Au[1] * (c[1] - uu[1]) * c[0] * c[2] / K[2]);
        invJ.SetValue(2, 2, (c[0] * c[1] * uu[2] * Au[2] + c[0] * Au[1] * c[2]
                                   * c[2] + c[1] * c[2] * c[2] * Au[0]) / K[2]);
        invJ.SetValue(2, 3, -c[0] * c[1] * c[2] / K[2]);
        invJ.SetValue(2, 4, 0.5 * c[0] * Au[1] * c[2] / K[2]);
        invJ.SetValue(2, 5, -0.5 * (Au[1] * c[0] + c[1] * Au[0]) * c[2] / K[2]);

        invJ.SetValue(3, 0, Au[0] * (Au[0] * c[2] * c[1] - uu[0] * c[2] * Au[1]
                                                - uu[0] * c[1] * Au[2]) / K[0]);
        invJ.SetValue(3, 1, -Au[0] * Au[1] * (c[1] - uu[1]) * c[2] / K[0]);
        invJ.SetValue(3, 2, -Au[0] * Au[2] * (c[2] - uu[2]) * c[1] / K[0]);
        invJ.SetValue(3, 3, -Au[0] * c[2] * c[1] / K[0]);
        invJ.SetValue(3, 4, 0.5 * Au[0] * Au[1] * c[2] / K[0]);
        invJ.SetValue(3, 5, 0.5 * Au[0] * c[1] * Au[2] / K[0]);

        invJ.SetValue(4, 0, Au[0] * Au[1] * (c[0] + uu[0]) * c[2] / K[1]);
        invJ.SetValue(4, 1, -Au[1] * (c[0] * Au[1] * c[2] + c[0] * uu[1] * Au[2]
                                                + c[2] * uu[1] * Au[0]) / K[1]);
        invJ.SetValue(4, 2, -Au[2] * Au[1] * (c[2] - uu[2]) * c[0] / K[1]);
        invJ.SetValue(4, 3, -c[0] * Au[1] * c[2] / K[1]);
        invJ.SetValue(4, 4, -0.5 * Au[1] * (c[0] * Au[2] + Au[0] * c[2]) / K[1]);
        invJ.SetValue(4, 5, 0.5 * Au[2] * Au[1] * c[0] / K[1]);

        invJ.SetValue(5, 0, Au[0] * Au[2] * (c[0] + uu[0]) * c[1] / K[2]);
        invJ.SetValue(5, 1, -Au[2] * Au[1] * (c[1] - uu[1]) * c[0] / K[2]);
        invJ.SetValue(5, 2, -Au[2] * (Au[2] * c[0] * c[1] + c[0] * uu[2] * Au[1]
                                                + c[1] * uu[2] * Au[0]) / K[2]);
        invJ.SetValue(5, 3, -Au[2] * c[0] * c[1] / K[2]);
        invJ.SetValue(5, 4, 0.5 * Au[2] * Au[1] * c[0] / K[2]);
        invJ.SetValue(5, 5, -0.5 * Au[2] * (Au[1] * c[0] + c[1] * Au[0]) / K[2]);
    }
    else if (type == "Merge")
    {
        NekMatrix<NekDouble> J(6, 6);
        Array<OneD, NekDouble> c(3);
        Array<OneD, NekDouble> K(3);

        for (int i = 0; i < 3; ++i)
        {
            GetC(c[i], beta[i], Au[i], A0[i], alpha[i]);
        }

        k = c[0] * Au[1] * c[2] + Au[0] * c[2] * c[1] + Au[2] * c[0] * c[1];
        K[0] = (c[0] - uu[0]) * k;
        K[1] = (c[1] + uu[1]) * k;
        K[2] = (c[2] + uu[2]) * k;

        invJ.SetValue(0, 0, (-c[1] * uu[0] * c[2] * Au[0] + Au[2] * c[1] * c[0]
                                   * c[0] + Au[1] * c[0] * c[0] * c[2]) / K[0]);
        invJ.SetValue(0, 1, Au[1] * (c[1] - uu[1]) * c[0] * c[2] / K[0]);
        invJ.SetValue(0, 2, Au[2] * (c[2] - uu[2]) * c[0] * c[1] / K[0]);
        invJ.SetValue(0, 3, c[0] * c[1] * c[2] / K[0]);
        invJ.SetValue(0, 4, -0.5 * c[0] * Au[1] * c[2] / K[0]);
        invJ.SetValue(0, 5, -0.5 * Au[2] * c[0] * c[1] / K[0]);

        invJ.SetValue(1, 0, Au[0] * (c[0] + uu[0]) * c[1] * c[2] / K[1]);
        invJ.SetValue(1, 1, (c[0] * uu[1] * c[2] * Au[1] + Au[2] * c[0] * c[1]
                                   * c[1] + c[2] * c[1] * c[1] * Au[0]) / K[1]);
        invJ.SetValue(1, 2, -Au[2] * (c[2] - uu[2]) * c[0] * c[1] / K[1]);
        invJ.SetValue(1, 3, -c[0] * c[1] * c[2] / K[1]);
        invJ.SetValue(1, 4, -0.5 * (c[0] * Au[2] + Au[0] * c[2]) * c[1] / K[1]);
        invJ.SetValue(1, 5, 0.5 * Au[2] * c[0] * c[1] / K[1]);

        invJ.SetValue(2, 0, Au[0] * (c[0] - uu[0]) * c[1] * c[2] / K[2]);
        invJ.SetValue(2, 1, -Au[1] * (c[1] + uu[1]) * c[0] * c[2] / K[2]);
        invJ.SetValue(2, 2, -(c[0] * uu[2] * c[1] * Au[2] - Au[1] * c[0] *
                              c[2] * c[2] - c[1] * c[2] * c[2] * Au[0]) / K[2]);
        invJ.SetValue(2, 3, -c[0] * c[1] * c[2] / K[2]);
        invJ.SetValue(2, 4, -0.5 * Au[1] * c[0] * c[2] / K[2]);
        invJ.SetValue(2, 5, 0.5 * (Au[1] * c[0] + Au[0] * c[1]) * c[2] / K[2]);

        invJ.SetValue(3, 0, -Au[0] * (Au[0] * c[2] * c[1] + uu[0] * c[2] *
                                         Au[1] + uu[0] * c[1] * Au[2]) / K[0]);
        invJ.SetValue(3, 1, Au[0] * Au[1] * (c[1] + uu[1]) * c[2] / K[0]);

        invJ.SetValue(3, 2, -Au[0] * Au[2] * (c[2] - uu[2]) * c[1] / K[0]);
        invJ.SetValue(3, 3, -Au[0] * c[2] * c[1] / K[0]);
        invJ.SetValue(3, 4, 0.5 * Au[0] * Au[1] * c[2] / K[0]);
        invJ.SetValue(3, 5, 0.5 * Au[0] * c[1] * Au[2] / K[0]);

        invJ.SetValue(4, 0, Au[0] * Au[1] * (c[0] + uu[0]) * c[2] / K[1]);
        invJ.SetValue(4, 1, -Au[1] * (c[0] * Au[1] * c[2] + c[0] * uu[1] * Au[2]
                                                + c[2] * uu[1] * Au[0]) / K[1]);
        invJ.SetValue(4, 2, -Au[2] * Au[1] * (c[2] - uu[2]) * c[0] / K[1]);
        invJ.SetValue(4, 3, -c[0] * Au[1] * c[2] / K[1]);
        invJ.SetValue(4, 4, -0.5 * Au[1] * (c[0] * Au[2] + Au[0] * c[2]) / K[1]);
        invJ.SetValue(4, 5, 0.5 * Au[2] * Au[1] * c[0] / K[1]);

        invJ.SetValue(5, 0, Au[0] * Au[2] * (c[0] + uu[0]) * c[1] / K[2]);
        invJ.SetValue(5, 1, -Au[2] * Au[1] * (c[1] - uu[1]) * c[0] / K[2]);
        invJ.SetValue(5, 2, -Au[2] * (Au[2] * c[0] * c[1] + c[0] * uu[2] * Au[1]
                                                + c[1] * uu[2] * Au[0]) / K[2]);
        invJ.SetValue(5, 3, -Au[2] * c[0] * c[1] / K[2]);
        invJ.SetValue(5, 4, 0.5 * Au[2] * Au[1] * c[0] / K[2]);
        invJ.SetValue(5, 5, -0.5 * Au[2] * (Au[1] * c[0] + c[1] * Au[0]) / K[2]);
    }
    else if (type == "Junction")
    {
        NekMatrix<NekDouble> J(4, 4);
        Array<OneD, NekDouble> c(2);
        Array<OneD, NekDouble> K(2);

        for (int i = 0; i < 2; ++i)
        {
            GetC(c[i], beta[i], Au[i], A0[i], alpha[i]);
        }

        k = (c[0] * Au[1] + Au[0] * c[1]);
        K[0] = (c[0] + uu[0]) * k;
        K[1] = (c[1] + uu[1]) * k;

        invJ.SetValue(0, 0, (Au[1] * c[0] * c[0] - c[1] * uu[0] * Au[0]) / K[0]);
        invJ.SetValue(0, 1, Au[1] * (c[1] - uu[1]) * c[0] / K[0]);
        invJ.SetValue(0, 2, c[0] * c[1] / K[0]);
        invJ.SetValue(0, 3, -0.5 * c[0] * Au[1] / K[0]);

        invJ.SetValue(1, 0, Au[0] * (c[0] + uu[0]) * c[1] / K[1]);
        invJ.SetValue(1, 1, (c[0] * uu[1] * Au[1] + c[1] * c[1] * Au[0]) / K[1]);
        invJ.SetValue(1, 2, -c[0] * c[1] / K[1]);
        invJ.SetValue(1, 3, -0.5 * Au[0] * c[1] / K[1]);

        invJ.SetValue(2, 0, Au[0] * (Au[0] * c[1] - uu[0] * Au[1]) / K[0]);
        invJ.SetValue(2, 1, -Au[0] * Au[1] * (c[1] - uu[1]) / K[0]);
        invJ.SetValue(2, 2, -Au[0] * c[1] / K[0]);
        invJ.SetValue(2, 3, 0.5 * Au[1] * Au[0] / K[0]);

        invJ.SetValue(3, 0, Au[0] * Au[1] * (c[0] + uu[0]) / K[1]);
        invJ.SetValue(3, 1, -Au[1] * (c[0] * Au[1] + uu[1] * Au[0]) / K[1]);
        invJ.SetValue(3, 2, -c[0] * Au[1] / K[1]);
        invJ.SetValue(3, 3, -0.5 * Au[1] * Au[0] / K[1]);
    }
}

} // namespace Nektar
