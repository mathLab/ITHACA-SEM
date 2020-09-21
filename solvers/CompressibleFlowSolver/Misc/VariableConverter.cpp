///////////////////////////////////////////////////////////////////////////////
//
// File VariableConverter.cpp
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
// Description: Auxiliary functions to convert variables in
//              the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include <iostream>

#include <CompressibleFlowSolver/Misc/VariableConverter.h>

using namespace std;

namespace Nektar
{
VariableConverter::VariableConverter(
    const LibUtilities::SessionReaderSharedPtr &pSession, const int spaceDim)
    : m_session(pSession), m_spacedim(spaceDim)
{
    // Create equation of state object
    std::string eosType;
    m_session->LoadSolverInfo("EquationOfState", eosType, "IdealGas");
    m_eos = GetEquationOfStateFactory().CreateInstance(eosType, m_session);

    // Parameters for dynamic viscosity
    m_session->LoadParameter("pInf", m_pInf, 101325);
    m_session->LoadParameter("rhoInf", m_rhoInf, 1.225);
    m_session->LoadParameter("GasConstant", m_gasConstant, 287.058);
    m_session->LoadParameter("mu", m_mu, 1.78e-05);

    // Parameters for sensor
    m_session->LoadParameter("Skappa", m_Skappa, -1.0);
    m_session->LoadParameter("Kappa", m_Kappa, 0.25);

}

/**
 * @brief Destructor for VariableConverter class.
 */
VariableConverter::~VariableConverter()
{
}

/**
 * @brief Compute the dynamic energy
 *        \f$ e = rho*V^2/2 \f$.
 */
void VariableConverter::GetDynamicEnergy(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &energy)
{
    size_t nPts = physfield[m_spacedim + 1].size();
    Vmath::Zero(nPts, energy, 1);

    // tmp = (rho * u_i)^2
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nPts, physfield[i + 1], 1, physfield[i + 1], 1, energy, 1,
                     energy, 1);
    }
    // Divide by rho and multiply by 0.5 --> tmp = 0.5 * rho * u^2
    Vmath::Vdiv(nPts, energy, 1, physfield[0], 1, energy, 1);
    Vmath::Smul(nPts, 0.5, energy, 1, energy, 1);
}

/**
 * @brief Compute the specific internal energy
 *        \f$ e = (E - rho*V^2/2)/rho \f$.
 */
void VariableConverter::GetInternalEnergy(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &energy)
{
    int nPts = physfield[0].size();
    Array<OneD, NekDouble> tmp(nPts);

    GetDynamicEnergy(physfield, tmp);

    // Calculate rhoe = E - rho*V^2/2
    Vmath::Vsub(nPts, physfield[m_spacedim + 1], 1, tmp, 1, energy, 1);
    // Divide by rho
    Vmath::Vdiv(nPts, energy, 1, physfield[0], 1, energy, 1);
}

/**
 * @brief Compute the specific enthalpy \f$ h = e + p/rho \f$.
 */
void VariableConverter::GetEnthalpy(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &enthalpy)
{
    int nPts = physfield[0].size();
    Array<OneD, NekDouble> energy(nPts, 0.0);
    Array<OneD, NekDouble> pressure(nPts, 0.0);

    GetInternalEnergy(physfield, energy);
    GetPressure(physfield, pressure);

    // Calculate p/rho
    Vmath::Vdiv(nPts, pressure, 1, physfield[0], 1, enthalpy, 1);
    // Calculate h = e + p/rho
    Vmath::Vadd(nPts, energy, 1, enthalpy, 1, enthalpy, 1);
}

/**
 * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
 * \f$ \rho\mathbf{v} \f$.
 *
 * @param physfield  Momentum field.
 * @param velocity   Velocity field.
 */
void VariableConverter::GetVelocityVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    const int nPts = physfield[0].size();

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vdiv(nPts, physfield[1 + i], 1, physfield[0], 1, velocity[i], 1);
    }
}

/**
 * @brief Compute the mach number \f$ M = \| \mathbf{v} \|^2 / c \f$.
 *
 * @param physfield    Input physical field.
 * @param soundfield   The speed of sound corresponding to physfield.
 * @param mach         The resulting mach number \f$ M \f$.
 */
void VariableConverter::GetMach(Array<OneD, Array<OneD, NekDouble>> &physfield,
                                Array<OneD, NekDouble> &soundspeed,
                                Array<OneD, NekDouble> &mach)
{
    const int nPts = physfield[0].size();

    Vmath::Vmul(nPts, physfield[1], 1, physfield[1], 1, mach, 1);

    for (int i = 1; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nPts, physfield[1 + i], 1, physfield[1 + i], 1, mach, 1,
                     mach, 1);
    }

    Vmath::Vdiv(nPts, mach, 1, physfield[0], 1, mach, 1);
    Vmath::Vdiv(nPts, mach, 1, physfield[0], 1, mach, 1);
    Vmath::Vsqrt(nPts, mach, 1, mach, 1);

    Vmath::Vdiv(nPts, mach, 1, soundspeed, 1, mach, 1);
}

/**
 * @brief Compute the dynamic viscosity using the Sutherland's law
 * \f$ \mu = \mu_star * (T / T_star)^3/2 * (1 + C) / (T / T_star + C) \f$,
 * where:   \mu_star = 1.7894 * 10^-5 Kg / (m * s)
 *          T_star   = 288.15 K
 *          C        = 110. / 288.15
 *
 * WARNING, if this routine is modified the same must be done in the
 * FieldConvert utility ProcessWSS.cpp (this class should be restructured).
 *
 * @param physfield    Input physical field.
 * @param mu           The resulting dynamic viscosity.
 */
void VariableConverter::GetDynamicViscosity(
    const Array<OneD, const NekDouble> &temperature, Array<OneD, NekDouble> &mu)
{
    const int nPts    = temperature.size();
    const NekDouble C = .38175;
    NekDouble mu_star = m_mu;
    NekDouble T_star  = m_pInf / (m_rhoInf * m_gasConstant);
    NekDouble ratio;

    for (int i = 0; i < nPts; ++i)
    {
        ratio = temperature[i] / T_star;
        mu[i] = mu_star * ratio * sqrt(ratio) * (1 + C) / (ratio + C);
    }
}

void VariableConverter::GetAbsoluteVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &Vtot)
{
    const int nPts = physfield[0].size();

    // Getting the velocity vector on the 2D normal space
    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);

    Vmath::Zero(Vtot.size(), Vtot, 1);

    for (int i = 0; i < m_spacedim; ++i)
    {
        velocity[i] = Array<OneD, NekDouble>(nPts);
    }

    GetVelocityVector(physfield, velocity);

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nPts, velocity[i], 1, velocity[i], 1, Vtot, 1, Vtot, 1);
    }

    Vmath::Vsqrt(nPts, Vtot, 1, Vtot, 1);
}

void VariableConverter::GetSensor(
    const MultiRegions::ExpListSharedPtr &field,
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    Array<OneD, NekDouble> &Sensor, Array<OneD, NekDouble> &SensorKappa,
    int offset)
{
    NekDouble Skappa;
    NekDouble order;
    Array<OneD, NekDouble> tmp;
    Array<OneD, int> expOrderElement = field->EvalBasisNumModesMaxPerExp();

    for (int e = 0; e < field->GetExpSize(); e++)
    {
        int numModesElement = expOrderElement[e];
        int nElmtPoints     = field->GetExp(e)->GetTotPoints();
        int physOffset      = field->GetPhys_Offset(e);
        int nElmtCoeffs     = field->GetExp(e)->GetNcoeffs();
        int numCutOff       = numModesElement - offset;

        if (numModesElement <= offset)
        {
            Vmath::Fill(nElmtPoints, 0.0,
                    tmp = Sensor + physOffset, 1);
            Vmath::Fill(nElmtPoints, 0.0,
                    tmp = SensorKappa + physOffset, 1);
            continue;
        }

        // create vector to save the solution points per element at P = p;
        Array<OneD, NekDouble> elmtPhys(nElmtPoints,
            tmp = physarray[0] + physOffset);
        // Compute coefficients
        Array<OneD, NekDouble> elmtCoeffs(nElmtCoeffs, 0.0);
        field->GetExp(e)->FwdTrans(elmtPhys, elmtCoeffs);

        // ReduceOrderCoeffs reduces the polynomial order of the solution
        // that is represented by the coeffs given as an inarray. This is
        // done by projecting the higher order solution onto the orthogonal
        // basis and padding the higher order coefficients with zeros.
        Array<OneD, NekDouble> reducedElmtCoeffs(nElmtCoeffs, 0.0);
        field->GetExp(e)->ReduceOrderCoeffs(numCutOff, elmtCoeffs,
                                            reducedElmtCoeffs);

        Array<OneD, NekDouble> reducedElmtPhys(nElmtPoints, 0.0);
        field->GetExp(e)->BwdTrans(reducedElmtCoeffs, reducedElmtPhys);

        NekDouble numerator   = 0.0;
        NekDouble denominator = 0.0;

        // Determining the norm of the numerator of the Sensor
        Array<OneD, NekDouble> difference(nElmtPoints, 0.0);
        Vmath::Vsub(nElmtPoints, elmtPhys, 1, reducedElmtPhys, 1, difference,
                    1);

        numerator = Vmath::Dot(nElmtPoints, difference, difference);
        denominator = Vmath::Dot(nElmtPoints, elmtPhys, elmtPhys);

        NekDouble elmtSensor = sqrt(numerator / denominator);
        elmtSensor = log10(max(elmtSensor, NekConstants::kNekSqrtTol));

        Vmath::Fill(nElmtPoints, elmtSensor, tmp = Sensor + physOffset, 1);

        // Compute reference value for sensor
        order = max(numModesElement-1, 1);
        if (order > 0 )
        {
            Skappa = m_Skappa - 4.25 * log10(static_cast<NekDouble>(order));
        }
        else
        {
            Skappa = 0.0;
        }

        // Compute artificial viscosity
        NekDouble elmtSensorKappa;
        if (elmtSensor < (Skappa-m_Kappa))
        {
            elmtSensorKappa = 0;
        }
        else if (elmtSensor > (Skappa + m_Kappa))
        {
            elmtSensorKappa = 1.0;
        }
        else
        {
            elmtSensorKappa = 0.5 *
                (1 + sin(M_PI * (elmtSensor - Skappa) / (2 * m_Kappa)));
        }
        Vmath::Fill(nElmtPoints, elmtSensorKappa,
                tmp = SensorKappa + physOffset, 1);
    }
}

/**
 * @brief Calculate the pressure using the equation of state.
 *
 * @param physfield  Input momentum.
 * @param pressure   Computed pressure field.
 */
void VariableConverter::GetPressure(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    int nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetInternalEnergy(physfield, energy);

    for (int i = 0; i < nPts; ++i)
    {
        pressure[i] = m_eos->GetPressure(physfield[0][i], energy[i]);
    }
}

/**
 * @brief Compute the temperature using the equation of state.
 *
 * @param physfield    Input physical field.
 * @param temperature  The resulting temperature \f$ T \f$.
 */
void VariableConverter::GetTemperature(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &temperature)
{
    int nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetInternalEnergy(physfield, energy);

    for (int i = 0; i < nPts; ++i)
    {
        temperature[i] = m_eos->GetTemperature(physfield[0][i], energy[i]);
    }
}

/**
 * @brief Compute the sound speed using the equation of state.
 *
 * @param physfield    Input physical field
 * @param soundspeed   The resulting sound speed \f$ c \f$.
 */
void VariableConverter::GetSoundSpeed(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &soundspeed)
{
    int nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetInternalEnergy(physfield, energy);

    for (int i = 0; i < nPts; ++i)
    {
        soundspeed[i] = m_eos->GetSoundSpeed(physfield[0][i], energy[i]);
    }
}

/**
 * @brief Compute the entropy using the equation of state.
 *
 * @param physfield    Input physical field
 * @param soundspeed   The resulting sound speed \f$ c \f$.
 */
void VariableConverter::GetEntropy(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &entropy)
{
    int nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetInternalEnergy(physfield, energy);

    for (int i = 0; i < nPts; ++i)
    {
        entropy[i] = m_eos->GetEntropy(physfield[0][i], energy[i]);
    }
}

/**
 * @brief Compute \f$ e(rho,p) \f$ using the equation of state.
 *
 * @param rho          Input density
 * @param pressure     Input pressure
 * @param energy       The resulting internal energy.
 */
void VariableConverter::GetEFromRhoP(const Array<OneD, NekDouble> &rho,
                                     const Array<OneD, NekDouble> &pressure,
                                     Array<OneD, NekDouble> &energy)
{
    int nPts = rho.size();

    for (int i = 0; i < nPts; ++i)
    {
        energy[i] = m_eos->GetEFromRhoP(rho[i], pressure[i]);
    }
}

/**
 * @brief Compute \f$ rho(p,T) \f$ using the equation of state.
 *
 * @param pressure     Input pressure
 * @param temperature  Input temperature
 * @param rho          The resulting density
 */
void VariableConverter::GetRhoFromPT(const Array<OneD, NekDouble> &pressure,
                                     const Array<OneD, NekDouble> &temperature,
                                     Array<OneD, NekDouble> &rho)
{
    int nPts = pressure.size();

    for (int i = 0; i < nPts; ++i)
    {
        rho[i] = m_eos->GetRhoFromPT(pressure[i], temperature[i]);
    }
}
}
