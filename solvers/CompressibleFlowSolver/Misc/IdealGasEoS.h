///////////////////////////////////////////////////////////////////////////////
//
// File: IdealGasEoS.h
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_IDEALGASEOS
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_MISC_IDEALGASEOS

#include "EquationOfState.h"


namespace Nektar
{

/**
* @brief Ideal gas equation of state:
 *       p = rho * R * T
*/
class IdealGasEoS : public EquationOfState
{
    public:

        friend class MemoryManager<IdealGasEoS>;

        /// Creates an instance of this class
        static EquationOfStateSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            EquationOfStateSharedPtr p = MemoryManager<IdealGasEoS>::
                                    AllocateSharedPtr(pSession);
            return p;
        }

        ///Name of the class
        static std::string className;

    protected:

        NekDouble GetTemperature(const NekDouble& rho, const NekDouble& e) final;

        vec_t GetTemperature(const vec_t& rho, const vec_t& e) final;

        NekDouble GetPressure(const NekDouble& rho, const NekDouble& e) final;

        vec_t GetPressure(const vec_t& rho, const vec_t& e) final;

        NekDouble v_GetSoundSpeed(const NekDouble &rho, const NekDouble &e) final;

        NekDouble v_GetEntropy(const NekDouble &rho, const NekDouble &e) final;

        NekDouble v_GetDPDrho_e(const NekDouble &rho, const NekDouble &e) final;

        NekDouble v_GetDPDe_rho(const NekDouble &rho, const NekDouble &e) final;

        NekDouble v_GetEFromRhoP(const NekDouble &rho, const NekDouble &p) final;

        NekDouble v_GetRhoFromPT(const NekDouble &rho, const NekDouble &p) final;

    private:
        IdealGasEoS(const LibUtilities::SessionReaderSharedPtr& pSession);

        ~IdealGasEoS(void){};


    // type agnostic kernels
    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T GetTemperatureKernel(const T& e)
    {
        return e * m_gammaMoneOgasConst;
    }

    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline T GetPressureKernel(const T& rho, const T& e)
    {
        return rho * e * m_gammaMone;
    }




};

} // namespace

#endif
