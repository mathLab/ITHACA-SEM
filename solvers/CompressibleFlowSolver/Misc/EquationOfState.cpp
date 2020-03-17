///////////////////////////////////////////////////////////////////////////////
//
// File: EquationOfState.cpp
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
// Description: Abstract base class for equations of state.
//
///////////////////////////////////////////////////////////////////////////////

#include "EquationOfState.h"

using namespace std;

namespace Nektar
{
EquationOfStateFactory &GetEquationOfStateFactory()
{
    static EquationOfStateFactory instance;
    return instance;
}

EquationOfState::EquationOfState(
    const LibUtilities::SessionReaderSharedPtr &pSession)
{
    pSession->LoadParameter("Gamma", m_gamma, 1.4);
    pSession->LoadParameter("GasConstant", m_gasConstant, 287.058);
}

NekDouble EquationOfState::GetTemperature(const NekDouble &rho,
                                          const NekDouble &e)
{
    return v_GetTemperature(rho, e);
}

NekDouble EquationOfState::GetPressure(const NekDouble &rho, const NekDouble &e)
{
    return v_GetPressure(rho, e);
}

NekDouble EquationOfState::GetSoundSpeed(const NekDouble &rho,
                                         const NekDouble &e)
{
    return v_GetSoundSpeed(rho, e);
}

NekDouble EquationOfState::GetEntropy(const NekDouble &rho, const NekDouble &e)
{
    return v_GetEntropy(rho, e);
}

NekDouble EquationOfState::GetDPDrho_e(const NekDouble &rho, const NekDouble &e)
{
    return v_GetDPDrho_e(rho, e);
}

NekDouble EquationOfState::GetDPDe_rho(const NekDouble &rho, const NekDouble &e)
{
    return v_GetDPDe_rho(rho, e);
}

NekDouble EquationOfState::GetEFromRhoP(const NekDouble &rho,
                                        const NekDouble &p)
{
    return v_GetEFromRhoP(rho, p);
}

NekDouble EquationOfState::GetRhoFromPT(const NekDouble &p, const NekDouble &T)
{
    return v_GetRhoFromPT(p, T);
}

// General implementation for v_GetSoundSpeed: c^2 = xi + kappa * h
//    where xi = dpdrho - e/rho * dp/de    and  kappa = dp/de / rho
NekDouble EquationOfState::v_GetSoundSpeed(const NekDouble &rho,
                                           const NekDouble &e)
{
    NekDouble p      = GetPressure(rho, e);
    NekDouble dpde   = GetDPDe_rho(rho, e);
    NekDouble dpdrho = GetDPDrho_e(rho, e);

    NekDouble enthalpy = e + p / rho;

    NekDouble chi   = dpdrho - e / rho * dpde;
    NekDouble kappa = dpde / rho;

    return sqrt(chi + kappa * enthalpy);
}
}
