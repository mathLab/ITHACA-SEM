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
    void GetEnthalpy(const Array<OneD, const Array<OneD, NekDouble>> &physfield,
                     Array<OneD, NekDouble> &enthalpy);
    void GetVelocityVector(const Array<OneD, Array<OneD, NekDouble>> &physfield,
              Array<OneD, Array<OneD, NekDouble>> &velocity);
    void GetMach(Array<OneD, Array<OneD, NekDouble>> &physfield,
                 Array<OneD, NekDouble> &soundspeed,
                 Array<OneD, NekDouble> &mach);
    void GetDynamicViscosity(const Array<OneD, const NekDouble> &temperature,
                             Array<OneD, NekDouble> &mu);
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
    void GetPressure(const Array<OneD, const Array<OneD, NekDouble>> &physfield,
                     Array<OneD, NekDouble> &pressure);
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
};
}
#endif
