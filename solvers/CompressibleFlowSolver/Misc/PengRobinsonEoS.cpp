///////////////////////////////////////////////////////////////////////////////
//
// File: PengRobinsonEoS.cpp
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
// Description: Peng-Robinson equation of state
//
///////////////////////////////////////////////////////////////////////////////

#include "PengRobinsonEoS.h"

using namespace std;

namespace Nektar
{

std::string PengRobinsonEoS::className =
    GetEquationOfStateFactory().RegisterCreatorFunction(
        "PengRobinson", PengRobinsonEoS::create,
        "Peng-Robinson equation of state.");

PengRobinsonEoS::PengRobinsonEoS(
    const LibUtilities::SessionReaderSharedPtr &pSession)
    : EquationOfState(pSession)
{
    pSession->LoadParameter("Tcrit", m_Tc);
    pSession->LoadParameter("Pcrit", m_Pc);
    pSession->LoadParameter("AcentricFactor", m_omega);

    m_a  = 0.45724 * m_gasConstant * m_gasConstant * m_Tc * m_Tc / m_Pc;
    m_b  = 0.0778 * m_gasConstant * m_Tc / m_Pc;
    m_fw = 0.37464 + 1.54226 * m_omega - 0.2699 * m_omega * m_omega;
}

NekDouble PengRobinsonEoS::v_GetTemperature(const NekDouble &rho,
                                            const NekDouble &e)
{
    // First we need to evaluate the log term
    //    ln[(1/rho + b - b*sqrt(2)) / (1/rho + b + b*sqrt(2))]
    NekDouble sqrt2   = sqrt(2.0);
    NekDouble logTerm = LogTerm(rho);

    // The temperature can be expressed as an equation in the form
    //      A * (T^1/2)^2 + B * T^1/2 + C = 0
    NekDouble A, B, C;

    A = m_gasConstant / (m_gamma - 1);
    B = -m_a / (m_b * 2 * sqrt2) * logTerm / sqrt(m_Tc) * m_fw * (1 + m_fw);
    C = m_a / (m_b * 2 * sqrt2) * logTerm * (1 + m_fw) * (1 + m_fw) - e;

    // Solve for T^1/2 (positive root)
    NekDouble sqrtT = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
    // Calculate the temperature
    return sqrtT * sqrtT;
}

NekDouble PengRobinsonEoS::v_GetPressure(const NekDouble &rho,
                                         const NekDouble &e)
{
    NekDouble T = GetTemperature(rho, e);

    NekDouble p =
        m_gasConstant * T / (1.0 / rho - m_b) -
        m_a * Alpha(T) / (1.0 / (rho * rho) + 2.0 * m_b / rho - m_b * m_b);

    return p;
}

NekDouble PengRobinsonEoS::v_GetEntropy(const NekDouble &rho,
                                        const NekDouble &e)
{
    NekDouble T       = GetTemperature(rho, e);
    NekDouble logTerm = LogTerm(rho);
    // Entropy for an ideal gas
    NekDouble sIg =
        m_gasConstant / (m_gamma - 1) * log(T) - m_gasConstant * log(rho);

    // First sqrt(Alpha) = 1+_w*(1-sqrt(Tr)) and sqrt(Tr)
    NekDouble sqrtA  = sqrt(Alpha(T));
    NekDouble sqrtTr = sqrt(T / m_Tc);

    NekDouble deltaS = m_gasConstant * log(1 - m_b * rho);
    deltaS += m_a * sqrtA * m_fw * logTerm * (sqrtTr / T) / (m_b * sqrt(8));

    return sIg + deltaS;
}

NekDouble PengRobinsonEoS::v_GetDPDrho_e(const NekDouble &rho,
                                         const NekDouble &e)
{
    NekDouble T    = GetTemperature(rho, e);
    NekDouble dPde = GetDPDe_rho(rho, e);

    // First calculate the denominator 1/rho^2 + 2*b/rho - b^2
    //    and alpha = [1+f_w*(1-sqrt(Tr))]^2
    NekDouble denom = 1.0 / (rho * rho) + 2.0 * m_b / rho - m_b * m_b;
    NekDouble alpha = Alpha(T);

    // Calculate dPdrho_T
    NekDouble dPdrho_T =
        m_gasConstant * T / ((1.0 - m_b * rho) * (1.0 - m_b * rho)) -
        2 * m_a * alpha * rho * (1 + m_b * rho) /
            ((denom * rho * rho) * (denom * rho * rho));

    // Calculate dedrho_T
    NekDouble dedrho_T =
        -m_a * sqrt(alpha) * (1.0 + m_fw) / (denom * rho * rho);

    // The result is dPdrho_e = dPdrho_T - dPde_rho * dedrho_T
    return dPdrho_T - dPde * dedrho_T;
}

NekDouble PengRobinsonEoS::v_GetDPDe_rho(const NekDouble &rho,
                                         const NekDouble &e)
{
    NekDouble T       = GetTemperature(rho, e);
    NekDouble logTerm = LogTerm(rho);

    // First calculate the denominator 1/rho^2 + 2*b/rho - b^2
    //    and sqrt(Alpha) = 1+f_w*(1-sqrt(Tr))
    NekDouble denom = 1.0 / (rho * rho) + 2.0 * m_b / rho - m_b * m_b;
    NekDouble sqrtA = sqrt(Alpha(T));

    // Compute cv = dedT_rho
    NekDouble cv = m_gasConstant / (m_gamma - 1) -
                   m_a / (2 * m_b * sqrt(8)) * logTerm * (m_fw * (1 + m_fw)) /
                       sqrt(T * m_Tc);

    // Now we obtain dPdT_rho
    NekDouble dPdT = m_gasConstant / (1.0 / rho - m_b) +
                     m_a / sqrt(T * m_Tc) * m_fw * sqrtA / denom;

    // The result is dPde_rho = dPdT_rho / cv
    return dPdT / cv;
}

NekDouble PengRobinsonEoS::v_GetEFromRhoP(const NekDouble &rho,
                                          const NekDouble &p)
{
    NekDouble denom   = 1.0 / (rho * rho) + 2.0 * m_b / rho - m_b * m_b;
    NekDouble logTerm = LogTerm(rho);
    // First we solve for the temperature, which can be expressed as
    //      A * (T^1/2)^2 + B * T^1/2 + C = 0
    NekDouble A, B, C;

    A = m_gasConstant / (1.0 / rho - m_b) -
        (m_a * m_fw * m_fw) / (denom * m_Tc);
    B = 2 * m_a / denom * m_fw * (1.0 + m_fw) / sqrt(m_Tc);
    C = -m_a * (1.0 + m_fw) * (1 + m_fw) / denom - p;

    // Solve for T^1/2 (positive root)
    NekDouble T = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
    T           = T * T;

    // Calculate alpha(T))
    NekDouble alpha = Alpha(T);
    // sqrt(Tr)
    NekDouble sqrtTr = sqrt(T / m_Tc);
    // Calculate internal energy
    return m_gasConstant * T / (m_gamma - 1) +
           m_a / (m_b * sqrt(8)) * logTerm *
               (alpha + sqrt(alpha) * m_fw * sqrtTr);
}

NekDouble PengRobinsonEoS::v_GetRhoFromPT(const NekDouble &p,
                                          const NekDouble &T)
{
    // First solve for the compressibility factor Z using the cubic equation
    //    Z^3 + k1 * Z^2 + k2 * Z + k3 = 0
    //    for PengRobinson:
    //        k1 = B-1, k2 = A - 2*B - 3*B^2,  k3 = - AB + B^2 + B^3
    //    where A = a*alpha(T)*P/(RT)^2, B = bP/(RT)
    NekDouble A = m_a * Alpha(T) * p / (m_gasConstant * m_gasConstant * T * T);
    NekDouble B = m_b * p / (m_gasConstant * T);

    NekDouble k1 = B - 1.0;
    NekDouble k2 = A - 2 * B - 3 * B * B;
    NekDouble k3 = -A * B + B * B + B * B * B;

    // Use ideal gas (Z=1) as starting guess for iteration
    NekDouble Z = 1.0;
    // Newton-Raphson iteration to find Z
    NekDouble tol      = 1e-6;
    NekDouble maxIter  = 100;
    NekDouble residual = 1;
    NekDouble f, df;
    unsigned int cnt = 0;
    while (abs(residual) > tol && cnt < maxIter)
    {
        f        = Z * Z * Z + k1 * Z * Z + k2 * Z + k3;
        df       = 3 * Z * Z + 2 * k1 * Z + k2;
        residual = f / df;
        Z -= residual;
        ++cnt;
    }
    if (cnt == maxIter)
    {
        cout << "Newton-Raphson in PengRobinsonEoS::v_GetRhoFromPT did not "
                "converge in "
             << maxIter << " iterations (residual = " << residual << ")"
             << endl;
    }

    // Now calculate rho = p/(ZRT)
    return p / (Z * m_gasConstant * T);
}

NekDouble PengRobinsonEoS::Alpha(const NekDouble &T)
{
    NekDouble sqrtAlpha = 1.0 + m_fw * (1.0 - sqrt(T / m_Tc));
    return sqrtAlpha * sqrtAlpha;
}

NekDouble PengRobinsonEoS::LogTerm(const NekDouble &rho)
{
    return log((1.0 / rho + m_b - m_b * sqrt(2)) /
               (1.0 / rho + m_b + m_b * sqrt(2)));
}
}
