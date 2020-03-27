///////////////////////////////////////////////////////////////////////////////
//
// File: EquationOfState.h
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_EQUATIONOFSTATE
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_EQUATIONOFSTATE

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>

namespace Nektar
{
//  Forward declaration
class EquationOfState;

/// A shared pointer to an equation of state object
typedef std::shared_ptr<EquationOfState> EquationOfStateSharedPtr;

/// Declaration of the equation of state factory
typedef LibUtilities::NekFactory<std::string, EquationOfState,
                                 const LibUtilities::SessionReaderSharedPtr &>
    EquationOfStateFactory;

/// Declaration of the equation of state factory singleton
EquationOfStateFactory &GetEquationOfStateFactory();

/**
 * @class EquationOfState
 * @brief Encapsulates equations of state allowing us to obtain thermodynamic
 *        properties: most relations are in the form X(rho,e)
 */
class EquationOfState
{
public:
    virtual ~EquationOfState()
    {
    }

    /// Calculate the temperature
    NekDouble GetTemperature(const NekDouble &rho, const NekDouble &e);

    /// Calculate the pressure
    NekDouble GetPressure(const NekDouble &rho, const NekDouble &e);

    /// Calculate the sound speed
    NekDouble GetSoundSpeed(const NekDouble &rho, const NekDouble &e);

    /// Calculate the entropy
    NekDouble GetEntropy(const NekDouble &rho, const NekDouble &e);

    /// Calculate the partial derivative of P(rho,e) with respect to rho
    NekDouble GetDPDrho_e(const NekDouble &rho, const NekDouble &e);

    /// Calculate the partial derivative of P(rho,e) with respect to e
    NekDouble GetDPDe_rho(const NekDouble &rho, const NekDouble &e);

    /// Obtain the internal energy from rho and P
    NekDouble GetEFromRhoP(const NekDouble &rho, const NekDouble &p);

    /// Obtain the density from P and T
    NekDouble GetRhoFromPT(const NekDouble &p, const NekDouble &T);

protected:
    NekDouble m_gamma;
    NekDouble m_gasConstant;

    /// Constructor
    EquationOfState(const LibUtilities::SessionReaderSharedPtr &pSession);

    virtual NekDouble v_GetTemperature(const NekDouble &rho,
                                       const NekDouble &e) = 0;

    virtual NekDouble v_GetPressure(const NekDouble &rho,
                                    const NekDouble &e) = 0;

    virtual NekDouble v_GetSoundSpeed(const NekDouble &rho, const NekDouble &e);

    virtual NekDouble v_GetEntropy(const NekDouble &rho,
                                   const NekDouble &e) = 0;

    virtual NekDouble v_GetDPDrho_e(const NekDouble &rho,
                                    const NekDouble &e) = 0;

    virtual NekDouble v_GetDPDe_rho(const NekDouble &rho,
                                    const NekDouble &e) = 0;

    virtual NekDouble v_GetEFromRhoP(const NekDouble &rho,
                                     const NekDouble &p) = 0;

    virtual NekDouble v_GetRhoFromPT(const NekDouble &rho,
                                     const NekDouble &p) = 0;
};
}

#endif
