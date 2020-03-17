///////////////////////////////////////////////////////////////////////////////
//
// File: VanDerWaalsEoS.h
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
// Description: Van der Waals equation of state
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_VANDERWAALSEOS
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_VANDERWAALSEOS

#include "EquationOfState.h"

namespace Nektar
{

/**
* @brief van der Waals equation of state:
 *       p = RT/(1/rho - b) - a * rho^2
 *       with a = 27/64 * (R*Tc)^2 / Pc
 *            b = 1/8   * (R*Tc) / Pc
*/
class VanDerWaalsEoS : public EquationOfState
{
public:
    friend class MemoryManager<VanDerWaalsEoS>;

    /// Creates an instance of this class
    static EquationOfStateSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession)
    {
        EquationOfStateSharedPtr p =
            MemoryManager<VanDerWaalsEoS>::AllocateSharedPtr(pSession);
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    NekDouble m_a;
    NekDouble m_b;

    virtual NekDouble v_GetTemperature(const NekDouble &rho,
                                       const NekDouble &e);

    virtual NekDouble v_GetPressure(const NekDouble &rho, const NekDouble &e);

    virtual NekDouble v_GetEntropy(const NekDouble &rho, const NekDouble &e);

    virtual NekDouble v_GetDPDrho_e(const NekDouble &rho, const NekDouble &e);

    virtual NekDouble v_GetDPDe_rho(const NekDouble &rho, const NekDouble &e);

    virtual NekDouble v_GetEFromRhoP(const NekDouble &rho, const NekDouble &p);

    virtual NekDouble v_GetRhoFromPT(const NekDouble &rho, const NekDouble &p);

private:
    VanDerWaalsEoS(const LibUtilities::SessionReaderSharedPtr &pSession);

    virtual ~VanDerWaalsEoS(void){};
};
}

#endif
