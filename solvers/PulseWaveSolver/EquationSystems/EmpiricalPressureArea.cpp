///////////////////////////////////////////////////////////////////////////////
//
// File EmpiricalPressureArea.cpp
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
// Description: EmpiricalPressureArea class
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/EmpiricalPressureArea.h>

using namespace std;

namespace Nektar
{

std::string EmpiricalPressureArea::className =
    GetPressureAreaFactory().RegisterCreatorFunction(
        "Empirical", EmpiricalPressureArea::create,
        "Empirical pressure area relationship for the arterial system");

EmpiricalPressureArea::EmpiricalPressureArea(
    Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
    const LibUtilities::SessionReaderSharedPtr pSession)
    : PulseWavePressureArea(pVessel, pSession)
{
}

EmpiricalPressureArea::~EmpiricalPressureArea()
{
}

void EmpiricalPressureArea::v_GetPressure(NekDouble &P, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &dAUdx,
                                 const NekDouble &gamma, const NekDouble &alpha)
{
    NekDouble kappa = 0.0;
    GetKappa(kappa, A, A0, alpha);

    P = m_PExt - beta * sqrt(A0) * log(kappa) / (2 * alpha)
                                                      - gamma * dAUdx / sqrt(A); //Viscoelasticity
}

void EmpiricalPressureArea::v_GetC(NekDouble &c, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &alpha)
{
    NekDouble kappa = 0.0;
    GetKappa(kappa, A, A0, alpha);

    c = sqrt(beta * sqrt(A0) / (2 * m_rho * kappa)); // Elastic
}

void EmpiricalPressureArea::v_GetW1(NekDouble &W1, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble I = 0.0;
    GetCharIntegral(I, beta, A, A0, alpha);

    W1 = u + I; // Elastic and assumes u0 = 0
}

void EmpiricalPressureArea::v_GetW2(NekDouble &W2, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble I = 0.0;
    GetCharIntegral(I, beta, A, A0, alpha);

    W2 = u - I; // Elastic and assumes u0 = 0
}

void EmpiricalPressureArea::v_GetAFromChars(NekDouble &A, const NekDouble &W1,
                const NekDouble &W2, const NekDouble &beta, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble xi = (W1 - W2) * sqrt(m_rho / (8 * beta * sqrt(A0)));

    A = A0 * exp(xi * (2 - alpha * xi));
}

void EmpiricalPressureArea::v_GetUFromChars(NekDouble &u, const NekDouble &W1,
                                                            const NekDouble &W2)
{
    u = (W1 + W2) / 2;
}

void EmpiricalPressureArea::v_GetCharIntegral(NekDouble &I,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble kappa = 0.0;
    GetKappa(kappa, A, A0, alpha);

    I = sqrt(2 * beta * sqrt(A0) / m_rho) * (1 - sqrt(kappa)) / alpha;
}

void EmpiricalPressureArea::v_GetJacobianInverse(NekMatrix<NekDouble> &invJ,
             const Array<OneD, NekDouble> &Au, const Array<OneD, NekDouble> &uu,
           const Array<OneD, NekDouble> &beta, const Array<OneD, NekDouble> &A0,
                   const Array<OneD, NekDouble> &alpha, const std::string &type)
{
    // General formulation
    if (type == "Bifurcation")
    {
        NekMatrix<NekDouble> J(6, 6);
        Array<OneD, NekDouble> c(3);

        for (int i = 0; i < 3; ++i)
        {
            GetC(c[i], beta[i], Au[i], A0[i], alpha[i]);
        }

        J.SetValue(0, 0, 1);
        J.SetValue(0, 1, 0);
        J.SetValue(0, 2, 0);
        J.SetValue(0, 3, c[0] / Au[0]);
        J.SetValue(0, 4, 0);
        J.SetValue(0, 5, 0);

        J.SetValue(1, 0, 0);
        J.SetValue(1, 1, 1);
        J.SetValue(1, 2, 0);
        J.SetValue(1, 3, 0);
        J.SetValue(1, 4, -c[1] / Au[1]);
        J.SetValue(1, 5, 0);

        J.SetValue(2, 0, 0);
        J.SetValue(2, 1, 0);
        J.SetValue(2, 2, 1);
        J.SetValue(2, 3, 0);
        J.SetValue(2, 4, 0);
        J.SetValue(2, 5, -c[2] / Au[2]);

        J.SetValue(3, 0,  Au[0]);
        J.SetValue(3, 1, -Au[1]);
        J.SetValue(3, 2, -Au[2]);
        J.SetValue(3, 3,  uu[0]);
        J.SetValue(3, 4, -uu[1]);
        J.SetValue(3, 5, -uu[2]);

        J.SetValue(4, 0,  2 * uu[0]);
        J.SetValue(4, 1, -2 * uu[1]);
        J.SetValue(4, 2, 0);
        J.SetValue(4, 3,  2 * c[0] * c[0] / Au[0]);
        J.SetValue(4, 4, -2 * c[1] * c[1] / Au[1]);
        J.SetValue(4, 5, 0);

        J.SetValue(5, 0,  2 * uu[0]);
        J.SetValue(5, 1, 0);
        J.SetValue(5, 2, -2 * uu[2]);
        J.SetValue(5, 3,  2 * c[0] * c[0] / Au[0]);
        J.SetValue(5, 4, 0);
        J.SetValue(5, 5, -2 * c[2] * c[2] / Au[2]);

        invJ = J;
        invJ.Invert();
    }
    else if (type == "Merge")
    {
        NekMatrix<NekDouble> J(6, 6);
        Array<OneD, NekDouble> c(3);

        for (int i = 0; i < 3; ++i)
        {
            GetC(c[i], beta[i], Au[i], A0[i], alpha[i]);
        }

        J.SetValue(0, 0, 1);
        J.SetValue(0, 1, 0);
        J.SetValue(0, 2, 0);
        J.SetValue(0, 3, -c[0] / Au[0]);
        J.SetValue(0, 4, 0);
        J.SetValue(0, 5, 0);

        J.SetValue(1, 0, 0);
        J.SetValue(1, 1, 1);
        J.SetValue(1, 2, 0);
        J.SetValue(1, 3, 0);
        J.SetValue(1, 4, c[1] / Au[1]);
        J.SetValue(1, 5, 0);

        J.SetValue(2, 0, 0);
        J.SetValue(2, 1, 0);
        J.SetValue(2, 2, 1);
        J.SetValue(2, 3, 0);
        J.SetValue(2, 4, 0);
        J.SetValue(2, 5, c[2] / Au[2]);

        J.SetValue(3, 0,  Au[0]);
        J.SetValue(3, 1, -Au[1]);
        J.SetValue(3, 2, -Au[2]);
        J.SetValue(3, 3,  uu[0]);
        J.SetValue(3, 4, -uu[1]);
        J.SetValue(3, 5, -uu[2]);

        J.SetValue(4, 0,  2 * uu[0]);
        J.SetValue(4, 1, -2 * uu[1]);
        J.SetValue(4, 2, 0);
        J.SetValue(4, 3,  2 * c[0] * c[0] / Au[0]);
        J.SetValue(4, 4, -2 * c[1] * c[1] / Au[1]);
        J.SetValue(4, 5, 0);

        J.SetValue(5, 0,  2 * uu[0]);
        J.SetValue(5, 1, 0);
        J.SetValue(5, 2, -2 * uu[2]);
        J.SetValue(5, 3,  2 * c[0] * c[0] / Au[0]);
        J.SetValue(5, 4, 0);
        J.SetValue(5, 5, -2 * c[2] * c[2] / Au[2]);

        invJ = J;
        invJ.Invert();
    }
    else if (type == "Junction")
    {
        NekMatrix<NekDouble> J(4, 4);
        Array<OneD, NekDouble> c(2);

        for (int i = 0; i < 2; ++i)
        {
            GetC(c[i], beta[i], Au[i], A0[i], alpha[i]);
        }

        J.SetValue(0, 0, 1);
        J.SetValue(0, 1, 0);
        J.SetValue(0, 2, c[0] / Au[0]);
        J.SetValue(0, 3, 0);

        J.SetValue(1, 0, 0);
        J.SetValue(1, 1, 1);
        J.SetValue(1, 2, 0);
        J.SetValue(1, 3, -c[1] / Au[1]);

        J.SetValue(2, 0,  Au[0]);
        J.SetValue(2, 1, -Au[1]);
        J.SetValue(2, 2,  uu[0]);
        J.SetValue(2, 3, -uu[1]);

        J.SetValue(3, 0,  2 * uu[0]);
        J.SetValue(3, 1, -2 * uu[1]);
        J.SetValue(3, 2,  2 * c[0] * c[0] / Au[0]);
        J.SetValue(3, 3, -2 * c[1] * c[1] / Au[1]);

        invJ = J;
        invJ.Invert();
    }
}

void EmpiricalPressureArea::GetKappa(NekDouble &kappa, const NekDouble &A,
                                    const NekDouble &A0, const NekDouble &alpha)
{
    kappa = 1 - alpha * log(A / A0);
}

} // namespace Nektar
