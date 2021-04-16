///////////////////////////////////////////////////////////////////////////////
//
// File VariableConverter.h
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_VARCONVERT_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_VARCONVERT_H

#include "EquationOfState.h"
#include <SolverUtils/UnsteadySystem.h>

namespace Nektar
{
// Forward declarations
class VariableConverter;
typedef std::shared_ptr<VariableConverter> VariableConverterSharedPtr;
/**
 *
 */
class VariableConverter
{
public:
    VariableConverter(const LibUtilities::SessionReaderSharedPtr &pSession,
                      const int spaceDim);

    ~VariableConverter();

    // Variable manipulations valid for all fluids
    void GetDynamicEnergy(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);
    void GetInternalEnergy(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);
    template <class T, typename = typename std::enable_if
    <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T GetInternalEnergy(T* physfield)
    {
        // get dynamic energy
        T oneOrho = 1.0 / physfield[0];
        T dynEne{};
        for (size_t d = 1; d < m_spacedim + 1; ++d)
        {
            T tmp = physfield[d]; //load 1x
            dynEne += tmp * tmp;
        }
        dynEne = 0.5 * dynEne * oneOrho;

        // Calculate rhoe = E - rho*V^2/2
        T energy = physfield[m_spacedim + 1] - dynEne;
        return energy * oneOrho;
    }
    void GetEnthalpy(const Array<OneD, const Array<OneD, NekDouble>> &physfield,
                     Array<OneD, NekDouble> &enthalpy);
    void GetVelocityVector(const Array<OneD, Array<OneD, NekDouble>> &physfield,
              Array<OneD, Array<OneD, NekDouble>> &velocity);
    void GetMach(Array<OneD, Array<OneD, NekDouble>> &physfield,
                 Array<OneD, NekDouble> &soundspeed,
                 Array<OneD, NekDouble> &mach);
    void GetDynamicViscosity(const Array<OneD, const NekDouble> &temperature,
                             Array<OneD, NekDouble> &mu);

    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T GetDynamicViscosity(T &temperature)
    {
        constexpr NekDouble C = .38175;
        constexpr NekDouble onePlusC = 1.0 + C;

        NekDouble mu_star = m_mu;

        T ratio = temperature * m_oneOverT_star;
        return mu_star * ratio * sqrt(ratio) * onePlusC / (ratio + C);
    }

    void GetAbsoluteVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &Vtot);
    void GetSensor(const MultiRegions::ExpListSharedPtr &field,
                   const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                   Array<OneD, NekDouble> &Sensor,
                   Array<OneD, NekDouble> &SensorKappa, int offset = 1);

    // Transformations depending on the equation of state
    void GetTemperature(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &temperature);
    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T GetTemperature(T* physfield)
    {
        T energy = GetInternalEnergy(physfield);
        return m_eos->GetTemperature(physfield[0], energy);
    }
    //
    void GetPressure(const Array<OneD, const Array<OneD, NekDouble>> &physfield,
                     Array<OneD, NekDouble> &pressure);
    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T GetPressure(T* physfield)
    {
        T energy = GetInternalEnergy(physfield);
        return m_eos->GetPressure(physfield[0], energy);
    }

    void GetSoundSpeed(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &soundspeed);
    void GetEntropy(const Array<OneD, const Array<OneD, NekDouble>> &physfield,
                    Array<OneD, NekDouble> &entropy);
    void GetEFromRhoP(const Array<OneD, NekDouble> &rho,
                      const Array<OneD, NekDouble> &pressure,
                      Array<OneD, NekDouble> &energy);
    void GetRhoFromPT(const Array<OneD, NekDouble> &pressure,
                      const Array<OneD, NekDouble> &temperature,
                      Array<OneD, NekDouble> &rho);
    void GetDmuDT(
        const Array<OneD, const NekDouble>  &temperature, 
        const Array<OneD, const NekDouble>  &mu, 
              Array<OneD, NekDouble>        &DmuDT);

    const EquationOfStateSharedPtr Geteos()
    {
        return m_eos;
    }

protected:
    LibUtilities::SessionReaderSharedPtr m_session;
    EquationOfStateSharedPtr m_eos;
    int m_spacedim;
    NekDouble m_pInf;
    NekDouble m_rhoInf;
    NekDouble m_gasConstant;
    NekDouble m_mu;
    NekDouble m_Skappa;
    NekDouble m_Kappa;
    NekDouble m_oneOverT_star;
};







}
#endif
