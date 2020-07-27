///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionIP.h
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
// Description: IP diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG
#define NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
namespace SolverUtils
{
class DiffusionIP : public Diffusion
{
public:
    static DiffusionSharedPtr create(std::string diffType)
    {
        boost::ignore_unused(diffType);
        return DiffusionSharedPtr(new DiffusionIP());
    }

    static std::string type;

    /// Calculate the average of conservative variables on traces
    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline void ConsVarAve(
    const size_t nConvectiveFields,
    const T& Bweight,
    const std::vector<T>& vFwd,
    const std::vector<T>& vBwd,
    std::vector<T>& aver)
    {
        constexpr int nvelst = 1;
        int nEngy  = nConvectiveFields - 1;
        int nveled = nEngy;

        T Fweight = 1.0 - Bweight;
        for (size_t v = 0; v < nEngy; ++v)
        {
            aver[v] = Fweight * vFwd[v] + Bweight * vBwd[v];
        }

        T LinternalEngy{};
        T RinternalEngy{};
        T AinternalEngy{};
        for (size_t j = nvelst; j < nveled; ++j)
        {
            LinternalEngy += vFwd[j] * vFwd[j];
            RinternalEngy += vBwd[j] * vBwd[j];
            AinternalEngy += aver[j] * aver[j];
        }
        LinternalEngy *= -0.5 / vFwd[0];
        RinternalEngy *= -0.5 / vBwd[0];
        LinternalEngy += vFwd[nEngy];
        RinternalEngy += vBwd[nEngy];

        aver[nEngy] = Fweight * LinternalEngy + Bweight * RinternalEngy;
        aver[nEngy] += AinternalEngy * (0.5 / aver[0]);

    }


protected:
    DiffusionIP();

    std::string m_shockCaptureType;
    NekDouble m_IPSymmFluxCoeff;
    NekDouble m_IP2ndDervCoeff;
    NekDouble m_IPPenaltyCoeff;

    Array<OneD, NekDouble> m_MuVarTrace;
    Array<OneD, Array<OneD, NekDouble>> m_traceNormals;
    Array<OneD, Array<OneD, NekDouble>> m_traceAver;
    Array<OneD, Array<OneD, NekDouble>> m_traceJump;
    /// Workspace for v_Diffusion
    Array<OneD, Array<OneD, NekDouble>> m_wspDiff;
    /// Workspace for CallTraceNumFlux
    TensorOfArray3D<NekDouble> m_wspNumDerivBwd;
    TensorOfArray3D<NekDouble> m_wspNumDerivFwd;

    Array<OneD, NekDouble> m_tracBwdWeightAver;
    Array<OneD, NekDouble> m_tracBwdWeightJump;
    Array<OneD, NekDouble> m_traceNormDirctnElmtLength;
    Array<OneD, NekDouble> m_traceNormDirctnElmtLengthRecip;
    LibUtilities::SessionReaderSharedPtr m_session;

    /// Get IP penalty factor based on order
    void GetPenaltyFactor(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        Array<OneD, NekDouble> &factor);
    /// Get a constant IP penalty factor
    void GetPenaltyFactor_const(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        Array<OneD, NekDouble> &factor);

    /// Add symmetric flux integration to field (in coefficient space)
    void AddSymmFluxIntegralToCoeff(
        const std::size_t nConvectiveFields, const size_t nDim,
        const size_t nPts, const size_t nTracePts,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, const int> &nonZeroIndex,
        TensorOfArray3D<NekDouble> &tracflux,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    /// Add symmetric flux integration to field (in physical space)
    void AddSymmFluxIntegralToPhys(
        const std::size_t nConvectiveFields, const size_t nDim,
        const size_t nPts, const size_t nTracePts,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, const int> &nonZeroIndex,
        TensorOfArray3D<NekDouble> &tracflux,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    /// Calculate symmetric flux on traces
    void CalTraceSymFlux(
        const std::size_t nConvectiveFields, const size_t nDim,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &solution_Aver,
        Array<OneD, Array<OneD, NekDouble>> &solution_jump,
        Array<OneD, int> &nonZeroIndexsymm,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &traceSymflux);

    /// Calculate symmetric flux on traces interface
    void DiffuseTraceSymmFlux(
        const std::size_t nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfield,
        TensorOfArray3D<NekDouble> &VolumeFlux,
        TensorOfArray3D<NekDouble> &SymmFlux,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd,
        Array<OneD, int> &nonZeroIndex);

    virtual void v_InitObject(
        LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields);

    virtual void v_Diffuse(
        const std::size_t nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    virtual void v_DiffuseCoeffs(
        const std::size_t nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    virtual void v_DiffuseVolumeFlux(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfields,
        TensorOfArray3D<NekDouble> &VolumeFlux, Array<OneD, int> &nonZeroIndex);

    virtual void v_DiffuseTraceFlux(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfields,
        TensorOfArray3D<NekDouble> &VolumeFlux,
        Array<OneD, Array<OneD, NekDouble>> &TraceFlux,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd,
        Array<OneD, int> &nonZeroIndex);
    virtual void v_AddDiffusionSymmFluxToCoeff(
        const std::size_t nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfield,
        TensorOfArray3D<NekDouble> &VolumeFlux,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    virtual void v_AddDiffusionSymmFluxToPhys(
        const std::size_t nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfield,
        TensorOfArray3D<NekDouble> &VolumeFlux,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    virtual void v_DiffuseCalculateDerivative(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfield,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    virtual void v_ConsVarAveJump(
        const std::size_t nConvectiveFields, const size_t npnts,
        const Array<OneD, const Array<OneD, NekDouble>> &vFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &vBwd,
        Array<OneD, Array<OneD, NekDouble>> &aver,
        Array<OneD, Array<OneD, NekDouble>> &jump);

    virtual const Array<OneD, const Array<OneD, NekDouble>> &v_GetTraceNormal()
    {
        return m_traceNormals;
    }

    /// Calculate numerical flux on traces
    void CalTraceNumFlux(
        const std::size_t nConvectiveFields, const size_t nDim,
        const size_t nPts, const size_t nTracePts,
        const NekDouble PenaltyFactor2,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        const TensorOfArray3D<NekDouble> &qfield,
        const Array<OneD, Array<OneD, NekDouble>> &vFwd,
        const Array<OneD, Array<OneD, NekDouble>> &vBwd,
        const Array<OneD, NekDouble> &MuVarTrace,
        Array<OneD, int> &nonZeroIndexflux,
        TensorOfArray3D<NekDouble> &traceflux,
        Array<OneD, Array<OneD, NekDouble>> &solution_Aver,
        Array<OneD, Array<OneD, NekDouble>> &solution_jump);

    /// Add second derivative term to trace jump (for DDG scheme)
    void AddSecondDerivToTrace(
        const std::size_t nConvectiveFields, const size_t nDim,
        const size_t nPts, const size_t nTracePts,
        const NekDouble PenaltyFactor2,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const TensorOfArray3D<NekDouble> &qfield,
        TensorOfArray3D<NekDouble> &numDerivFwd,
        TensorOfArray3D<NekDouble> &numDerivBwd);
};
} // namespace SolverUtils
} // namespace Nektar

#endif
