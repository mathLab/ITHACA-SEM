///////////////////////////////////////////////////////////////////////////////
//
// File: RedlichKwongEoS.cpp
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
// Description: Redlich-Kwong equation of state
//
///////////////////////////////////////////////////////////////////////////////

#include "RedlichKwongEoS.h"

using namespace std;

namespace Nektar
{

std::string RedlichKwongEoS::className =
    GetEquationOfStateFactory().RegisterCreatorFunction(
        "RedlichKwong", RedlichKwongEoS::create,
        "Redlich-Kwong equation of state.");

RedlichKwongEoS::RedlichKwongEoS(
    const LibUtilities::SessionReaderSharedPtr &pSession)
    : EquationOfState(pSession)
{
    pSession->LoadParameter("Tcrit", m_Tc);
    pSession->LoadParameter("Pcrit", m_Pc);

    m_a = 0.42748 * m_gasConstant * m_gasConstant * m_Tc * m_Tc / m_Pc;
    m_b = 0.08664 * m_gasConstant * m_Tc / m_Pc;
}

NekDouble RedlichKwongEoS::v_GetTemperature(const NekDouble &rho,
                                            const NekDouble &e)
{
    // First we need to evaluate the log term
    //    ln[1 + b*rho]
    NekDouble logTerm = LogTerm(rho);

    // The temperature can be expressed as an equation in the form
    //      (T^1/2)^3 + A* T^1/2 + B = 0, which we solve iteratively
    NekDouble A, B;

    A = -e * (m_gamma - 1) / m_gasConstant;
    B = -3.0 * m_a / (2.0 * m_b * m_gasConstant) * (m_gamma - 1) * sqrt(m_Tc) *
        logTerm;

    // Use ideal gas solution as starting guess for iteration
    NekDouble sqrtT = sqrt(e * (m_gamma - 1) / m_gasConstant);
    // Newton-Raphson iteration to find T^(1/2)
    NekDouble tol      = 1e-6;
    NekDouble maxIter  = 100;
    NekDouble residual = 1;
    NekDouble f, df;
    unsigned int cnt = 0;
    while (abs(residual) > tol && cnt < maxIter)
    {
        f        = sqrtT * sqrtT * sqrtT + A * sqrtT + B;
        df       = 3 * sqrtT * sqrtT + A;
        residual = f / df;
        sqrtT -= residual;
        ++cnt;
    }
    if (cnt == maxIter)
    {
        cout << "Newton-Raphson in RedlichKwongEoS::v_GetTemperature did not "
                "converge in "
             << maxIter << " iterations (residual = " << residual << ")"
             << endl;
    }

    // Calculate the temperature
    return sqrtT * sqrtT;
}

NekDouble RedlichKwongEoS::v_GetPressure(const NekDouble &rho,
                                         const NekDouble &e)
{
    NekDouble T = GetTemperature(rho, e);

    NekDouble p = m_gasConstant * T / (1.0 / rho - m_b) -
                  m_a * Alpha(T) / (1.0 / (rho * rho) + m_b / rho);

    return p;
}

NekDouble RedlichKwongEoS::v_GetEntropy(const NekDouble &rho,
                                        const NekDouble &e)
{
    NekDouble T       = GetTemperature(rho, e);
    NekDouble logTerm = LogTerm(rho);
    // Entropy for an ideal gas
    NekDouble sIg =
        m_gasConstant / (m_gamma - 1) * log(T) - m_gasConstant * log(rho);

    NekDouble deltaS = m_gasConstant * log(1 - m_b * rho);
    deltaS -= m_a * Alpha(T) * logTerm / (2 * m_b * T);

    return sIg + deltaS;
}

NekDouble RedlichKwongEoS::v_GetDPDrho_e(const NekDouble &rho,
                                         const NekDouble &e)
{
    NekDouble T     = GetTemperature(rho, e);
    NekDouble alpha = Alpha(T);
    NekDouble dPde  = GetDPDe_rho(rho, e);

    // Calculate dPdrho_T
    NekDouble dPdrho_T =
        m_gasConstant * T / ((1.0 - m_b * rho) * (1.0 - m_b * rho)) -
        m_a * alpha * rho * (m_b * rho + 2) /
            ((1 + m_b * rho) * (1 + m_b * rho));

    // Calculate dedrho_T
    NekDouble dedrho_T = -3 * m_a * alpha / (2 * (1 + m_b * rho));

    // The result is dPdrho_e = dPdrho_T - dPde_rho * dedrho_T
    return dPdrho_T - dPde * dedrho_T;
}

NekDouble RedlichKwongEoS::v_GetDPDe_rho(const NekDouble &rho,
                                         const NekDouble &e)
{
    NekDouble T       = GetTemperature(rho, e);
    NekDouble alpha   = Alpha(T);
    NekDouble logTerm = LogTerm(rho);

    // First calculate the denominator 1/rho^2 + 2*b/rho - b^2
    //    and sqrt(Alpha) = 1+f_w*(1-sqrt(Tr))
    NekDouble denom = 1.0 / (rho * rho) + m_b / rho;

    // Compute cv = dedT_rho
    NekDouble cv = m_gasConstant / (m_gamma - 1) +
                   3 * m_a * alpha * logTerm / (4 * m_b * T);

    // Now we obtain dPdT_rho
    NekDouble dPdT =
        m_gasConstant / (1.0 / rho - m_b) + m_a * alpha / (denom * 2 * T);

    // The result is dPde_rho = dPdT_rho / cv
    return dPdT / cv;
}

NekDouble RedlichKwongEoS::v_GetEFromRhoP(const NekDouble &rho,
                                          const NekDouble &p)
{
    NekDouble logTerm = LogTerm(rho);
    // First calculate the temperature, which can be expressed as
    //      (T^1/2)^3 + A* T^1/2 + B = 0
    NekDouble A, B;

    A = -p * (1.0 / rho - m_b) / m_gasConstant;
    B = -m_a * sqrt(m_Tc) * (1.0 / rho - m_b) /
        (1.0 / (rho * rho) + m_b / rho) / m_gasConstant;

    // Use ideal gas solution as starting guess for iteration
    NekDouble sqrtT = sqrt(p / (rho * (m_gamma - 1)));
    // Newton-Raphson iteration to find T^(1/2)
    NekDouble tol      = 1e-6;
    NekDouble maxIter  = 100;
    NekDouble residual = 1;
    NekDouble f, df;
    unsigned int cnt = 0;
    while (abs(residual) > tol && cnt < maxIter)
    {
        f        = sqrtT * sqrtT * sqrtT + A * sqrtT + B;
        df       = 3 * sqrtT * sqrtT + A;
        residual = f / df;
        sqrtT -= residual;
        ++cnt;
    }
    if (cnt == maxIter)
    {
        cout << "Newton-Raphson in RedlichKwongEoS::v_GetEFromRhoP did not "
                "converge in "
             << maxIter << " iterations (residual = " << residual << ")"
             << endl;
    }

    // Calculate T
    NekDouble T = sqrtT * sqrtT;

    // Calculate internal energy
    return m_gasConstant * T / (m_gamma - 1) -
           3 * m_a * Alpha(T) / (2 * m_b) * logTerm;
}

NekDouble RedlichKwongEoS::v_GetRhoFromPT(const NekDouble &p,
                                          const NekDouble &T)
{
    // First solve for the compressibility factor Z using the cubic equation
    //    Z^3 + k1 * Z^2 + k2 * Z + k3 = 0
    //    for RedlichKwong:
    //        k1 = -1.0, k2 = A - B - B^2,  k3 = -AB
    //    where A = a*alpha(T)*P/(RT)^2, B = bP/(RT)
    NekDouble A = m_a * Alpha(T) * p / (m_gasConstant * m_gasConstant * T * T);
    NekDouble B = m_b * p / (m_gasConstant * T);

    NekDouble k1 = -1.0;
    NekDouble k2 = A - B - B * B;
    NekDouble k3 = -A * B;

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
        cout << "Newton-Raphson in RedlichKwongEoS::v_GetRhoFromPT did not "
                "converge in "
             << maxIter << " iterations (residual = " << residual << ")"
             << endl;
    }

    // Now calculate rho = p/(ZRT)
    return p / (Z * m_gasConstant * T);
}

NekDouble RedlichKwongEoS::Alpha(const NekDouble &T)
{
    return 1.0 / sqrt(T / m_Tc);
}

NekDouble RedlichKwongEoS::LogTerm(const NekDouble &rho)
{
    return log(1 + m_b * rho);
}
}
