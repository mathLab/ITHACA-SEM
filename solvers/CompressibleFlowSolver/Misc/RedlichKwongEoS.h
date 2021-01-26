///////////////////////////////////////////////////////////////////////////////
//
// File: RedlichKwongEoS.h
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_REDLICHKWONGEOS
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_REDLICHKWONGEOS

#include "EquationOfState.h"

#include <LibUtilities/SimdLib/io.hpp>

namespace Nektar
{

/**
* @brief Redlich-Kwong equation of state:
 *       p = RT/(1/rho - b) - a/( sqrt(T / Tc) * (1/rho^2 + b/rho)
 *       with a = 0.42748 * (R*Tc)^2 / Pc
 *            b = 0.08664 * (R*Tc) / Pc
*/
class RedlichKwongEoS : public EquationOfState
{
public:
    friend class MemoryManager<RedlichKwongEoS>;

    /// Creates an instance of this class
    static EquationOfStateSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession)
    {
        EquationOfStateSharedPtr p =
            MemoryManager<RedlichKwongEoS>::AllocateSharedPtr(pSession);
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    NekDouble m_a;
    NekDouble m_b;
    NekDouble m_Tc;
    NekDouble m_Pc;

    NekDouble GetTemperature(const NekDouble& rho, const NekDouble& e) final;

    vec_t GetTemperature(const vec_t& rho, const vec_t& e) final;

    NekDouble GetPressure(const NekDouble& rho, const NekDouble& e) final;

    vec_t GetPressure(const vec_t& rho, const vec_t& e) final;

    NekDouble v_GetEntropy(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetDPDrho_e(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetDPDe_rho(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetEFromRhoP(const NekDouble &rho, const NekDouble &p) final;

    NekDouble v_GetRhoFromPT(const NekDouble &rho, const NekDouble &p) final;

private:
    RedlichKwongEoS(const LibUtilities::SessionReaderSharedPtr &pSession);

    ~RedlichKwongEoS(void){};

    // Alpha term of Redlich-Kwong EoS ( 1.0/sqrt(Tr))
    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T Alpha(const T& temp)
    {
        return 1.0 / sqrt(temp / m_Tc);
    }

    // Log term term of Peng-Robinson EoS
    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T LogTerm(const T& rho)
    {
        return log(1.0 + m_b * rho);
    }

    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T GetTemperatureKernel(const T& rho, const T& e)
    {
        // First we need to evaluate the log term
        //    ln[1 + b*rho]
        T logTerm = LogTerm(rho);

        // The temperature can be expressed as an equation in the form
        //      (T^1/2)^3 + A* T^1/2 + B = 0, which we solve iteratively
        T A = e * (1.0 - m_gamma) / m_gasConstant;
        T B = -3.0 * m_a / (2.0 * m_b * m_gasConstant) * (m_gamma - 1) * sqrt(m_Tc)
            * logTerm;

        // Use ideal gas solution as starting guess for iteration
        T sqrtT = sqrt(e * (m_gamma - 1) / m_gasConstant);
        // Newton-Raphson iteration to find T^(1/2)
        T tol      = 1e-6;
        T residual = 1;
        unsigned int maxIter  = 100;
        unsigned int cnt = 0;
        while (abs(residual) > tol && cnt < maxIter)
        {
            T f = sqrtT * sqrtT * sqrtT + A * sqrtT + B;
            T df = 3 * sqrtT * sqrtT + A;
            residual = f / df;
            sqrtT -= residual;
            ++cnt;
        }
        if (cnt == maxIter)
        {
            std::cout << "Newton-Raphson in RedlichKwongEoS::v_GetTemperature did not "
                    "converge in "
                 << maxIter << " iterations (residual = " << residual << ")"
                 << std::endl;
        }

        // Calculate the temperature
        return sqrtT * sqrtT;
    }

    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T GetPressureKernel(const T& rho, const T& e)
    {
        T temp = GetTemperatureKernel(rho, e);
        T oneOrho = 1.0 / rho;
        T p = m_gasConstant * temp / (oneOrho - m_b) - m_a * Alpha(temp) /
            (oneOrho * (oneOrho + m_b));
        return p;
    }




};
}
#endif
