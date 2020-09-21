///////////////////////////////////////////////////////////////////////////////
//
// File PowerPressureArea.cpp
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
// Description: PowerPressureArea class
//
///////////////////////////////////////////////////////////////////////////////
#include <PulseWaveSolver/EquationSystems/PowerPressureArea.h>

using namespace std;

namespace Nektar
{

std::string PowerPressureArea::className =
    GetPressureAreaFactory().RegisterCreatorFunction(
        "Power", PowerPressureArea::create,
        "Power law pressure area relationship for the arterial system");

PowerPressureArea::PowerPressureArea(
    Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
    const LibUtilities::SessionReaderSharedPtr pSession)
    : PulseWavePressureArea(pVessel, pSession)
{
    m_session->LoadParameter("P_Collapse", P_Collapse, -13.3322); // -10mmHg converted to kg / (cm s^2)
}

PowerPressureArea::~PowerPressureArea()
{
}

void PowerPressureArea::v_GetPressure(NekDouble &P, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &dAUdx,
                                 const NekDouble &gamma, const NekDouble &alpha)
{
    NekDouble c0 = 0.0;
    GetC0(c0, beta, A0);

    NekDouble b = 0.0;
    GetB(b, c0);

    P = m_PExt + (2 * m_rho * c0 * c0 / b) * (pow((A / A0), b / 2) - 1) // Power law by Smith/Canic/Mynard
                                           - A0 * gamma * dAUdx / (A * sqrt(A)); // Viscoelasticity
}

void PowerPressureArea::v_GetC(NekDouble &c, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &alpha)
{
    NekDouble c0 = 0.0;
    GetC0(c0, beta, A0);

    NekDouble b = 0.0;
    GetB(b, c0);

    c = c0 * pow((A / A0), b / 4); // Elastic
}

void PowerPressureArea::v_GetW1(NekDouble &W1, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble I = 0.0;
    GetCharIntegral(I, beta, A, A0, alpha);

    W1 = u + I; // Elastic and assumes u0 = 0
}

void PowerPressureArea::v_GetW2(NekDouble &W2, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble I = 0.0;
    GetCharIntegral(I, beta, A, A0, alpha);

    W2 = u - I; // Elastic and assumes u0 = 0
}

void PowerPressureArea::v_GetAFromChars(NekDouble &A, const NekDouble &W1,
                const NekDouble &W2, const NekDouble &beta, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    NekDouble c0 = 0.0;
    GetC0(c0, beta, A0);

    NekDouble b = 0.0;
    GetB(b, c0);

    A = A0 * pow(((b / (8 * c0)) * (W1 - W2) + 1), 4 / b);
}

void PowerPressureArea::v_GetUFromChars(NekDouble &u, const NekDouble &W1,
                                                            const NekDouble &W2)
{
    u = (W1 + W2) / 2;
}

void PowerPressureArea::v_GetCharIntegral(NekDouble &I, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &alpha)
{
    NekDouble c = 0.0;
    NekDouble c0 = 0.0;

    GetC0(c0, beta, A0);
    GetC(c, beta, A, A0);

    NekDouble b = 0.0;
    GetB(b, c0);

    I = (4 / b) * (c - c0);
}

void PowerPressureArea::v_GetJacobianInverse(NekMatrix<NekDouble> &invJ,
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

void PowerPressureArea::GetC0(NekDouble &c0, const NekDouble &beta,
                                                            const NekDouble &A0)
{
    // Reference c0 from the beta law
    c0 = sqrt(beta * sqrt(A0) / (2 * m_rho));

    /*
    // Empirical approximation from Olufsen et al (1999)
    NekDouble k1 = 3E3;
    NekDouble k2 = -9;
    NekDouble k3 = 337;
    NekDouble PI = 3.14159265359;

    NekDouble R0 = sqrt(A0 / PI);

    c0 = sqrt(2 / (3 * m_rho) * (k1 * exp(k2 * R0) + k3));
    */
}

void PowerPressureArea::GetB(NekDouble &b, const NekDouble &c0)
{
    b = 2 * m_rho * c0 * c0 / (m_PExt - P_Collapse);
}

} // namespace Nektar
