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
// Description: Auxiliary functions to convert variables in 
//              the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>

#include <CompressibleFlowSolver/Misc/VariableConverter.h>

using namespace std;

namespace Nektar
{
    VariableConverter::VariableConverter(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const int                                   spaceDim)
        : m_session(pSession),
            m_spacedim(spaceDim)
    {
        m_session->LoadParameter ("Gamma", m_gamma, 1.4);
        m_session->LoadParameter ("pInf", m_pInf, 101325);
        m_session->LoadParameter ("rhoInf", m_rhoInf, 1.225);
        m_session->LoadParameter ("GasConstant",   m_gasConstant,   287.058);
        m_session->LoadParameter ("mu",            m_mu,            1.78e-05);

        if( m_session->DefinesParameter("thermalConductivity"))
        {
            ASSERTL0( !m_session->DefinesParameter("Pr"),
                 "Cannot define both Pr and thermalConductivity.");

            m_session->LoadParameter ("thermalConductivity",
                                        m_thermalConductivity);
        }
        else
        {
            NekDouble Pr, Cp;
            m_session->LoadParameter ("Pr",
                                        Pr, 0.72);
            Cp = m_gamma / (m_gamma - 1.0) * m_gasConstant;
            m_thermalConductivity = Cp * m_mu / Pr;
        }
    }

    /**
     * @brief Destructor for VariableConverter class.
     */
    VariableConverter::~VariableConverter()
    {

    }

    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal
     * gas law.
     *
     * @param physfield  Input momentum.
     * @param pressure   Computed pressure field.
     */
    void VariableConverter::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD,                   NekDouble>   &pressure)
    {
        int       nPts  = physfield[0].num_elements();
        NekDouble alpha = -0.5;

        // Calculate ||rho v||^2
        Vmath::Vmul(nPts, physfield[1], 1, physfield[1], 1, pressure, 1);
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nPts, physfield[1+i], 1, physfield[1+i], 1,
                               pressure,       1, pressure,       1);
        }
        // Divide by rho to get rho*||v||^2
        Vmath::Vdiv (nPts, pressure, 1, physfield[0], 1, pressure, 1);
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(nPts, alpha,
                     pressure, 1, physfield[m_spacedim+1], 1, pressure, 1);
        // Multiply by (gamma-1)
        Vmath::Smul (nPts, m_gamma-1, pressure, 1, pressure, 1);
    }

    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal
     * gas law.
     *
     * This is a slightly optimised way to calculate the pressure field which
     * avoids division by the density field if the velocity field has already
     * been calculated.
     *
     * @param physfield  Input momentum.
     * @param velocity   Velocity vector.
     * @param pressure   Computed pressure field.
     */
    void VariableConverter::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        const Array<OneD, const Array<OneD, NekDouble> > &velocity,
              Array<OneD,                   NekDouble>   &pressure)
    {
        int nPts = physfield[0].num_elements();
        NekDouble alpha = -0.5;

        // Calculate ||\rho v||^2.
        Vmath::Vmul (nPts, velocity[0], 1, physfield[1], 1, pressure, 1);
        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nPts, velocity[i], 1, physfield[1+i], 1,
                               pressure,    1, pressure,       1);
        }
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(nPts,     alpha,
                     pressure, 1, physfield[m_spacedim+1], 1, pressure, 1);
        // Multiply by (gamma-1).
        Vmath::Smul (nPts, m_gamma-1, pressure, 1, pressure, 1);
    }

    /**
     * @brief Compute the enthalpy term \f$ H = E + p/rho \$.
     */
    void VariableConverter::GetEnthalpy(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &pressure,
                  Array<OneD,                   NekDouble>   &enthalpy)
    {
        int nPts  = physfield[0].num_elements();
        Array<OneD, NekDouble> tmp(nPts, 0.0);

        // Calculate E = rhoE/rho
        Vmath::Vdiv(nPts, physfield[m_spacedim+1], 1, physfield[0], 1, tmp, 1);
        // Calculate p/rho
        Vmath::Vdiv(nPts, pressure, 1, physfield[0], 1, enthalpy, 1);
        // Calculate H = E + p/rho
        Vmath::Vadd(nPts, tmp, 1, enthalpy, 1, enthalpy, 1);
    }

    /**
     * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
     * \f$ \rho\mathbf{v} \f$.
     *
     * @param physfield  Momentum field.
     * @param velocity   Velocity field.
     */
    void VariableConverter::GetVelocityVector(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &velocity)
    {
        const int nPts = physfield[0].num_elements();

        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vdiv(nPts, physfield[1+i], 1, physfield[0], 1,
                              velocity[i],    1);
        }
    }

    /**
     * @brief Compute the temperature \f$ T = p/\rho R \f$.
     *
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param temperature  The resulting temperature \f$ T \f$.
     */
    void VariableConverter::GetTemperature(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        Array<OneD,                         NekDouble  > &pressure,
        Array<OneD,                         NekDouble  > &temperature)
    {
        const int nPts = physfield[0].num_elements();

        Vmath::Vdiv(nPts, pressure, 1, physfield[0], 1, temperature, 1);
        Vmath::Smul(nPts, 1.0/m_gasConstant, temperature, 1, temperature, 1);
    }

    /**
     * @brief Compute the sound speed \f$ c = sqrt(\gamma p/\rho) \f$.
     *
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param soundspeed   The resulting sound speed \f$ c \f$.
     */
    void VariableConverter::GetSoundSpeed(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD,             NekDouble  > &pressure,
              Array<OneD,             NekDouble  > &soundspeed)
    {
        const int nPts = physfield[0].num_elements();
        Vmath::Vdiv (nPts, pressure, 1, physfield[0], 1, soundspeed, 1);
        Vmath::Smul (nPts, m_gamma, soundspeed, 1, soundspeed, 1);
        Vmath::Vsqrt(nPts, soundspeed, 1, soundspeed, 1);
    }

    /**
     * @brief Compute the mach number \f$ M = \| \mathbf{v} \|^2 / c \f$.
     *
     * @param physfield    Input physical field.
     * @param soundfield   The speed of sound corresponding to physfield.
     * @param mach         The resulting mach number \f$ M \f$.
     */
    void VariableConverter::GetMach(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD,             NekDouble  > &soundspeed,
        Array<OneD,             NekDouble  > &mach)
    {
        const int nPts = physfield[0].num_elements();

        Vmath::Vmul(nPts, physfield[1], 1, physfield[1], 1, mach, 1);

        for (int i = 1; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nPts, physfield[1+i], 1, physfield[1+i], 1,
                             mach,           1, mach,           1);
        }

        Vmath::Vdiv(nPts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(nPts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vsqrt(nPts, mach, 1, mach, 1);

        Vmath::Vdiv(nPts, mach, 1, soundspeed,   1, mach, 1);
    }

    /**
     * @brief Compute the dynamic viscosity using the Sutherland's law
     * \f$ \mu = \mu_star * (T / T_star)^3/2 * (T_star + 110) / (T + 110) \f$,
     * where:   \mu_star = 1.7894 * 10^-5 Kg / (m * s)
     *          T_star   = 288.15 K
     *
     * @param physfield    Input physical field.
     * @param mu           The resulting dynamic viscosity.
     */
    void VariableConverter::GetDynamicViscosity(
        const Array<OneD, const NekDouble> &temperature,
              Array<OneD,       NekDouble> &mu)
    {
        const int nPts    = temperature.num_elements();
        NekDouble mu_star = m_mu;
        NekDouble T_star  = m_pInf / (m_rhoInf * m_gasConstant);
        NekDouble ratio;

        for (int i = 0; i < nPts; ++i)
        {
            ratio = temperature[i] / T_star;
            mu[i] = mu_star * ratio * sqrt(ratio) *
                    (T_star + 110.0) / (temperature[i] + 110.0);
        }
    }

    /**
     * @brief Calculate entropy.
     */
    void VariableConverter::GetEntropy(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        const Array<OneD, const NekDouble>               &pressure,
        const Array<OneD, const NekDouble>               &temperature,
              Array<OneD,       NekDouble>               &entropy)
    {
        const int nPts = physfield[0].num_elements();
        const NekDouble temp_inf = m_pInf/(m_rhoInf*m_gasConstant);;

        for (int i = 0; i < nPts; ++i)
        {
            entropy[i] = m_gamma / (m_gamma - 1.0) * m_gasConstant *
                            log(temperature[i]/temp_inf) - m_gasConstant *
                            log(pressure[i] / m_pInf);
        }
    }

    void VariableConverter::GetAbsoluteVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,                   NekDouble>   &Vtot)
    {
        const int nPts = inarray[0].num_elements();

        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);

        Vmath::Zero(Vtot.num_elements(), Vtot, 1);

        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nPts);
        }

        GetVelocityVector(inarray, velocity);

        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nPts,
                         velocity[i], 1,
                         velocity[i], 1,
                         Vtot, 1,
                         Vtot, 1);
        }

        Vmath::Vsqrt(nPts,Vtot,1,Vtot,1);
    }

}
