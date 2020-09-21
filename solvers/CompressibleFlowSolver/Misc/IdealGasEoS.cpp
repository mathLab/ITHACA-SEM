///////////////////////////////////////////////////////////////////////////////
//
// File: IdealGasEoS.cpp
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
// Description: Ideal gas equation of state
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "IdealGasEoS.h"

using namespace std;

namespace Nektar
{

std::string IdealGasEoS::className = GetEquationOfStateFactory().
    RegisterCreatorFunction("IdealGas",
                            IdealGasEoS::create,
                            "Ideal gas equation of state.");

IdealGasEoS::IdealGasEoS(const LibUtilities::SessionReaderSharedPtr& pSession)
    : EquationOfState(pSession)
{

}

NekDouble IdealGasEoS::v_GetTemperature(
    const NekDouble &rho, const NekDouble &e)
{
    boost::ignore_unused(rho);
    return e*(m_gamma-1)/m_gasConstant;
}

NekDouble IdealGasEoS::v_GetPressure(
    const NekDouble &rho, const NekDouble &e)
{
    return rho*e*(m_gamma-1);
}

NekDouble IdealGasEoS::v_GetSoundSpeed(
    const NekDouble &rho, const NekDouble &e)
{
    NekDouble T = GetTemperature(rho,e);
    return sqrt(m_gamma * m_gasConstant * T);
}

NekDouble IdealGasEoS::v_GetEntropy(
    const NekDouble &rho, const NekDouble &e)
{
    NekDouble T = GetTemperature(rho,e);
    return m_gasConstant/(m_gamma-1) * log(T) - m_gasConstant * log(rho);
}

NekDouble IdealGasEoS::v_GetDPDrho_e(
    const NekDouble &rho, const NekDouble &e)
{
    boost::ignore_unused(rho);
    return e*(m_gamma-1);
}

NekDouble IdealGasEoS::v_GetDPDe_rho(
    const NekDouble &rho, const NekDouble &e)
{
    boost::ignore_unused(e);
    return rho*(m_gamma-1);
}

NekDouble IdealGasEoS::v_GetEFromRhoP(
            const NekDouble &rho, const NekDouble &p)
{
    return p / (rho * (m_gamma-1));
}

NekDouble IdealGasEoS::v_GetRhoFromPT(
            const NekDouble &p, const NekDouble &T)
{
    return p/(m_gasConstant*T);
}

}
