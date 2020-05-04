///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionIP.cpp
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

#include <SolverUtils/Diffusion/DiffusionIP.h>
#include <iomanip>
#include <iostream>

namespace Nektar
{
namespace SolverUtils
{
std::string DiffusionIP::type = GetDiffusionFactory().RegisterCreatorFunction(
    "InteriorPenalty", DiffusionIP::create);

DiffusionIP::DiffusionIP()
{
}

void DiffusionIP::v_InitObject(
    LibUtilities::SessionReaderSharedPtr pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
    m_session = pSession;

    m_session->LoadSolverInfo("ShockCaptureType", m_shockCaptureType, "Off");

    m_session->LoadParameter("IPSymmFluxCoeff", m_IPSymmFluxCoeff, 0.0); // 0.5

    m_session->LoadParameter("IP2ndDervCoeff", m_IP2ndDervCoeff,
                             0.0); // 1.0/12.0

    m_session->LoadParameter("IPPenaltyCoeff", m_IPPenaltyCoeff,
                             4.0); // 1.0/12.0

    // Setting up the normals
    size_t nDim      = pFields[0]->GetCoordim(0);
    size_t nVariable = pFields.size();
    size_t nTracePts = pFields[0]->GetTrace()->GetTotPoints();

    m_traceNormals = Array<OneD, Array<OneD, NekDouble>>{nDim};
    for (int i = 0; i < nDim; ++i)
    {
        m_traceNormals[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
    }
    m_traceAver = Array<OneD, Array<OneD, NekDouble>>{nVariable};
    m_traceJump = Array<OneD, Array<OneD, NekDouble>>{nVariable};
    for (int i = 0; i < nVariable; ++i)
    {
        m_traceAver[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
        m_traceJump[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
    }

    pFields[0]->GetTrace()->GetNormals(m_traceNormals);
    Array<OneD, NekDouble> lengthFwd{nTracePts, 0.0};
    Array<OneD, NekDouble> lengthBwd{nTracePts, 0.0};
    pFields[0]->GetTrace()->GetElmtNormalLength(lengthFwd, lengthBwd);

    const MultiRegions::AssemblyMapDGSharedPtr TraceMap =
        pFields[0]->GetTraceMap();
    pFields[0]->PeriodicBwdCopy(lengthFwd, lengthBwd);
    TraceMap->GetAssemblyCommDG()->PerformExchange(lengthFwd, lengthBwd);

    Vmath::Vadd(nTracePts, lengthBwd, 1, lengthFwd, 1, lengthFwd, 1);
    m_traceNormDirctnElmtLength      = lengthFwd;
    m_traceNormDirctnElmtLengthRecip = lengthBwd;
    Vmath::Sdiv(nTracePts, 1.0, m_traceNormDirctnElmtLength, 1,
                m_traceNormDirctnElmtLengthRecip, 1);

    m_tracBwdWeightAver = Array<OneD, NekDouble>{nTracePts, 0.0};
    m_tracBwdWeightJump = Array<OneD, NekDouble>{nTracePts, 0.0};
    pFields[0]->GetBwdWeight(m_tracBwdWeightAver, m_tracBwdWeightJump);
    Array<OneD, NekDouble> tmpBwdWeight{nTracePts, 0.0};
    Array<OneD, NekDouble> tmpBwdWeightJump{nTracePts, 0.0};
    for (int i = 1; i < nVariable; ++i)
    {
        pFields[i]->GetBwdWeight(tmpBwdWeight, tmpBwdWeightJump);
        Vmath::Vsub(nTracePts, tmpBwdWeight, 1, m_tracBwdWeightAver, 1,
                    tmpBwdWeight, 1);
        Vmath::Vabs(nTracePts, tmpBwdWeight, 1, tmpBwdWeight, 1);
        NekDouble norm = 0.0;
        for (int j = 0; j < nTracePts; ++j)
        {
            norm += tmpBwdWeight[j];
        }
        ASSERTL0(norm < 1.0E-11,
                 "different BWD for different variable not coded yet");
    }

    m_MuVarTrace = NullNekDouble1DArray;
    if (m_ArtificialDiffusionVector)
    {
        m_MuVarTrace = Array<OneD, NekDouble>{nTracePts, 0.0};
    }
}

void DiffusionIP::v_Diffuse(
    const std::size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{

    size_t nCoeffs = fields[0]->GetNcoeffs();
    Array<OneD, Array<OneD, NekDouble>> tmp{nConvectiveFields};
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        tmp[i] = Array<OneD, NekDouble>{nCoeffs, 0.0};
    }
    DiffusionIP::v_DiffuseCoeffs(nConvectiveFields, fields, inarray, tmp, pFwd,
                                 pBwd);
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        fields[i]->BwdTrans(tmp[i], outarray[i]);
    }
}

void DiffusionIP::v_DiffuseCoeffs(
    const std::size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    size_t nDim      = fields[0]->GetCoordim(0);
    size_t nPts      = fields[0]->GetTotPoints();
    size_t nCoeffs   = fields[0]->GetNcoeffs();
    size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> Fwd{nTracePts, 0.0};
    Array<OneD, NekDouble> Bwd{nTracePts, 0.0};

    TensorOfArray3D<NekDouble> elmtFlux{nDim};
    TensorOfArray3D<NekDouble> qfield{nDim};
    for (int j = 0; j < nDim; ++j)
    {
        qfield[j]   = Array<OneD, Array<OneD, NekDouble>>{nConvectiveFields};
        elmtFlux[j] = Array<OneD, Array<OneD, NekDouble>>{nConvectiveFields};
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            qfield[j][i] = Array<OneD, NekDouble>{nPts, 0.0};
        }
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            elmtFlux[j][i] = Array<OneD, NekDouble>{nPts, 0.0};
        }
    }

    Array<OneD, Array<OneD, NekDouble>> vFwd{nConvectiveFields};
    Array<OneD, Array<OneD, NekDouble>> vBwd{nConvectiveFields};
    if (pFwd == NullNekDoubleArrayofArray || pBwd == NullNekDoubleArrayofArray)
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            vFwd[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
            vBwd[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
        }
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            fields[i]->GetFwdBwdTracePhys(inarray[i], vFwd[i], vBwd[i]);
        }
    }
    else
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            vFwd[i] = pFwd[i];
            vBwd[i] = pBwd[i];
        }
    }

    DiffuseCalculateDerivative(fields, inarray, qfield, vFwd, vBwd);

    Array<OneD, int> nonZeroIndex;
    DiffuseVolumeFlux(fields, inarray, qfield, elmtFlux, nonZeroIndex);

    Array<OneD, Array<OneD, NekDouble>> tmpFluxIprdct{nDim};
    // volume intergration: the nonZeroIndex indicates which flux is nonzero
    for (int i = 0; i < nonZeroIndex.size(); ++i)
    {
        int j = nonZeroIndex[i];
        for (int k = 0; k < nDim; ++k)
        {
            tmpFluxIprdct[k] = elmtFlux[k][j];
        }
        fields[j]->IProductWRTDerivBase(tmpFluxIprdct, outarray[j]);
        Vmath::Neg(nCoeffs, outarray[j], 1);
    }
    // release qfield, elmtFlux and muvar;
    for (int j = 0; j < nDim; ++j)
    {
        elmtFlux[j] = NullNekDoubleArrayofArray;
    }

    Array<OneD, Array<OneD, NekDouble>> Traceflux{nConvectiveFields};
    for (int j = 0; j < nConvectiveFields; ++j)
    {
        Traceflux[j] = Array<OneD, NekDouble>{nTracePts, 0.0};
    }

    DiffuseTraceFlux(fields, inarray, qfield, elmtFlux, Traceflux, vFwd, vBwd,
                     nonZeroIndex);

    for (int i = 0; i < nonZeroIndex.size(); ++i)
    {
        int j = nonZeroIndex[i];

        fields[j]->AddTraceIntegral(Traceflux[j], outarray[j]);
        fields[j]->SetPhysState(false);
    }

    AddDiffusionSymmFluxToCoeff(nConvectiveFields, fields, inarray, qfield,
                                elmtFlux, outarray, vFwd, vBwd);

    for (int i = 0; i < nonZeroIndex.size(); ++i)
    {
        int j = nonZeroIndex[i];

        fields[j]->MultiplyByElmtInvMass(outarray[j], outarray[j]);
    }
}

void DiffusionIP::v_DiffuseCalculateDerivative(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfield,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    size_t nConvectiveFields = fields.size();
    boost::ignore_unused(pFwd, pBwd);

    size_t nDim = fields[0]->GetCoordim(0);

    Array<OneD, Array<OneD, NekDouble>> qtmp{3};
    for (int nd = 0; nd < 3; ++nd)
    {
        qtmp[nd] = NullNekDouble1DArray;
    }
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        for (int nd = 0; nd < nDim; ++nd)
        {
            qtmp[nd] = qfield[nd][i];
        }
        fields[i]->PhysDeriv(inarray[i], qtmp[0], qtmp[1], qtmp[2]);
    }
}

void DiffusionIP::v_DiffuseVolumeFlux(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfield, TensorOfArray3D<NekDouble> &VolumeFlux,
    Array<OneD, int> &nonZeroIndex)
{
    size_t nDim = fields[0]->GetCoordim(0);
    size_t nPts = fields[0]->GetTotPoints();

    Array<OneD, NekDouble> muvar = NullNekDouble1DArray;
    if (m_ArtificialDiffusionVector)
    {
        muvar = Array<OneD, NekDouble>{nPts, 0.0};
        GetAVmu(fields, inarray, muvar, m_MuVarTrace);
    }

    Array<OneD, Array<OneD, NekDouble>> tmparray2D = NullNekDoubleArrayofArray;

    m_FunctorDiffusionfluxCons(nDim, inarray, qfield, VolumeFlux, nonZeroIndex,
                               tmparray2D, muvar);
}

void DiffusionIP::v_DiffuseTraceFlux(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfield, TensorOfArray3D<NekDouble> &VolumeFlux,
    Array<OneD, Array<OneD, NekDouble>> &TraceFlux,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd,
    Array<OneD, int> &nonZeroIndex)
{
    boost::ignore_unused(VolumeFlux);
    size_t nDim      = fields[0]->GetCoordim(0);
    size_t nPts      = fields[0]->GetTotPoints();
    size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> traceflux3D(1);
    traceflux3D[0] = TraceFlux;

    size_t nConvectiveFields = fields.size();
    CalTraceNumFlux(nConvectiveFields, nDim, nPts, nTracePts, m_IP2ndDervCoeff,
                    fields, inarray, qfield, pFwd, pBwd, m_MuVarTrace,
                    nonZeroIndex, traceflux3D, m_traceAver, m_traceJump);
}

void DiffusionIP::v_AddDiffusionSymmFluxToCoeff(
    const std::size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfield, TensorOfArray3D<NekDouble> &VolumeFlux,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    if (fabs(m_IPSymmFluxCoeff) > 1.0E-12)
    {
        size_t nDim      = fields[0]->GetCoordim(0);
        size_t nPts      = fields[0]->GetTotPoints();
        size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();
        TensorOfArray3D<NekDouble> traceSymflux{nDim};
        for (int nd = 0; nd < nDim; ++nd)
        {
            traceSymflux[nd] =
                Array<OneD, Array<OneD, NekDouble>>{nConvectiveFields};
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                traceSymflux[nd][j] = Array<OneD, NekDouble>{nTracePts, 0.0};
            }
        }
        Array<OneD, int> nonZeroIndex;
        DiffuseTraceSymmFlux(nConvectiveFields, fields, inarray, qfield,
                             VolumeFlux, traceSymflux, pFwd, pBwd,
                             nonZeroIndex);

        AddSymmFluxIntegralToCoeff(nConvectiveFields, nDim, nPts, nTracePts,
                                   fields, nonZeroIndex, traceSymflux,
                                   outarray);
    }
}

void DiffusionIP::v_AddDiffusionSymmFluxToPhys(
    const std::size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfield, TensorOfArray3D<NekDouble> &VolumeFlux,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    if (fabs(m_IPSymmFluxCoeff) > 1.0E-12)
    {
        size_t nDim      = fields[0]->GetCoordim(0);
        size_t nPts      = fields[0]->GetTotPoints();
        size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();
        TensorOfArray3D<NekDouble> traceSymflux{nDim};
        for (int nd = 0; nd < nDim; ++nd)
        {
            traceSymflux[nd] =
                Array<OneD, Array<OneD, NekDouble>>{nConvectiveFields};
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                traceSymflux[nd][j] = Array<OneD, NekDouble>{nTracePts, 0.0};
            }
        }
        Array<OneD, int> nonZeroIndex;
        DiffuseTraceSymmFlux(nConvectiveFields, fields, inarray, qfield,
                             VolumeFlux, traceSymflux, pFwd, pBwd,
                             nonZeroIndex);

        AddSymmFluxIntegralToPhys(nConvectiveFields, nDim, nPts, nTracePts,
                                  fields, nonZeroIndex, traceSymflux, outarray);
    }
}

void DiffusionIP::DiffuseTraceSymmFlux(
    const std::size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfield, TensorOfArray3D<NekDouble> &VolumeFlux,
    TensorOfArray3D<NekDouble> &SymmFlux,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd,
    Array<OneD, int> &nonZeroIndex)
{
    boost::ignore_unused(inarray, qfield, VolumeFlux, pFwd, pBwd);
    size_t nDim = fields[0]->GetCoordim(0);

    CalTraceSymFlux(nConvectiveFields, nDim, fields, m_traceAver, m_traceJump,
                    nonZeroIndex, SymmFlux);
}

void DiffusionIP::CalTraceSymFlux(
    const std::size_t nConvectiveFields, const size_t nDim,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &solution_Aver,
    Array<OneD, Array<OneD, NekDouble>> &solution_jump,
    Array<OneD, int> &nonZeroIndexsymm,
    TensorOfArray3D<NekDouble> &traceSymflux)
{
    size_t nTracePts = solution_jump[nConvectiveFields - 1].size();

    for (int i = 0; i < nConvectiveFields; ++i)
    {
        Vmath::Smul(nTracePts, m_IPSymmFluxCoeff, solution_jump[i], 1,
                    solution_jump[i], 1);
    }

    m_FunctorSymmetricfluxCons(nDim, solution_Aver, solution_jump, traceSymflux,
                               nonZeroIndexsymm, m_traceNormals);

    for (int i = 0; i < nConvectiveFields; ++i)
    {
        MultiRegions::ExpListSharedPtr tracelist = fields[i]->GetTrace();
        for (int nd = 0; nd < nDim; ++nd)
        {
            tracelist->MultiplyByQuadratureMetric(traceSymflux[nd][i],
                                                  traceSymflux[nd][i]);
        }
    }
}

void DiffusionIP::AddSymmFluxIntegralToCoeff(
    const std::size_t nConvectiveFields, const size_t nDim, const size_t nPts,
    const size_t nTracePts,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, const int> &nonZeroIndex,
    TensorOfArray3D<NekDouble> &tracflux,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    boost::ignore_unused(nTracePts);

    size_t nCoeffs = outarray[nConvectiveFields - 1].size();
    Array<OneD, NekDouble> tmpCoeff{nCoeffs, 0.0};
    Array<OneD, Array<OneD, NekDouble>> tmpfield(nDim);
    for (int i = 0; i < nDim; ++i)
    {
        tmpfield[i] = Array<OneD, NekDouble>{nPts, 0.0};
    }
    int nv = 0;
    for (int j = 0; j < nonZeroIndex.size(); ++j)
    {
        nv = nonZeroIndex[j];
        for (int nd = 0; nd < nDim; ++nd)
        {
            Vmath::Zero(nPts, tmpfield[nd], 1);

            fields[nv]->AddTraceQuadPhysToField(tracflux[nd][nv],
                                                tracflux[nd][nv], tmpfield[nd]);
            fields[nv]->DivideByQuadratureMetric(tmpfield[nd], tmpfield[nd]);
        }
        fields[nv]->IProductWRTDerivBase(tmpfield, tmpCoeff);
        Vmath::Vadd(nCoeffs, tmpCoeff, 1, outarray[nv], 1, outarray[nv], 1);
    }
}

void DiffusionIP::AddSymmFluxIntegralToPhys(
    const std::size_t nConvectiveFields, const size_t nDim, const size_t nPts,
    const size_t nTracePts,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, const int> &nonZeroIndex,
    TensorOfArray3D<NekDouble> &tracflux,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    boost::ignore_unused(nTracePts);

    size_t nCoeffs = outarray[nConvectiveFields - 1].size();
    Array<OneD, NekDouble> tmpCoeff{nCoeffs, 0.0};
    Array<OneD, NekDouble> tmpPhysi{nPts, 0.0};
    Array<OneD, Array<OneD, NekDouble>> tmpfield{nDim};
    for (int i = 0; i < nDim; ++i)
    {
        tmpfield[i] = Array<OneD, NekDouble>{nPts, 0.0};
    }
    for (int j = 0; j < nonZeroIndex.size(); ++j)
    {
        int nv = nonZeroIndex[j];
        for (int nd = 0; nd < nDim; ++nd)
        {
            Vmath::Zero(nPts, tmpfield[nd], 1);

            fields[nv]->AddTraceQuadPhysToField(tracflux[nd][nv],
                                                tracflux[nd][nv], tmpfield[nd]);
            fields[nv]->DivideByQuadratureMetric(tmpfield[nd], tmpfield[nd]);
        }
        fields[nv]->IProductWRTDerivBase(tmpfield, tmpCoeff);
        fields[nv]->BwdTrans(tmpCoeff, tmpPhysi);
        Vmath::Vadd(nPts, tmpPhysi, 1, outarray[nv], 1, outarray[nv], 1);
    }
}

void DiffusionIP::GetPenaltyFactor(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    Array<OneD, NekDouble> &factor)
{
    MultiRegions::ExpListSharedPtr tracelist = fields[0]->GetTrace();
    std::shared_ptr<LocalRegions::ExpansionVector> traceExp =
        tracelist->GetExp();
    size_t ntotTrac = (*traceExp).size();
    int nTracPnt, noffset;

    const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap =
        fields[0]->GetLocTraceToTraceMap();

    const Array<OneD, const Array<OneD, int>> LRAdjExpid =
        locTraceToTraceMap->GetLeftRightAdjacentExpId();
    const Array<OneD, const Array<OneD, bool>> LRAdjflag =
        locTraceToTraceMap->GetLeftRightAdjacentExpFlag();

    std::shared_ptr<LocalRegions::ExpansionVector> fieldExp =
        fields[0]->GetExp();

    Array<OneD, NekDouble> factorFwdBwd{2, 0.0};

    NekDouble spaceDim = NekDouble(fields[0]->GetCoordim(0));

    for (int ntrace = 0; ntrace < ntotTrac; ++ntrace)
    {
        noffset  = tracelist->GetPhys_Offset(ntrace);
        nTracPnt = tracelist->GetTotPoints(ntrace);

        factorFwdBwd[0] = 0.0;
        factorFwdBwd[1] = 0.0;

        for (int nlr = 0; nlr < 2; ++nlr)
        {
            if (LRAdjflag[nlr][ntrace])
            {
                int numModes = fields[0]->GetNcoeffs(LRAdjExpid[nlr][ntrace]);
                NekDouble numModesdir =
                    pow(NekDouble(numModes), (1.0 / spaceDim));
                factorFwdBwd[nlr] = 1.0 * numModesdir * (numModesdir + 1.0);
            }
        }

        for (int np = 0; np < nTracPnt; ++np)
        {
            factor[noffset + np] = std::max(factorFwdBwd[0], factorFwdBwd[1]);
        }
    }
}

void DiffusionIP::GetPenaltyFactor_const(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    Array<OneD, NekDouble> &factor)
{
    boost::ignore_unused(fields);
    Vmath::Fill(factor.size(), m_IPPenaltyCoeff, factor, 1);
}

void DiffusionIP::v_ConsVarAveJump(
    const std::size_t nConvectiveFields, const size_t npnts,
    const Array<OneD, const Array<OneD, NekDouble>> &vFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &vBwd,
    Array<OneD, Array<OneD, NekDouble>> &aver,
    Array<OneD, Array<OneD, NekDouble>> &jump)
{
    ConsVarAve(nConvectiveFields, npnts, vFwd, vBwd, aver);

    m_SpecialBndTreat(aver);

    // note: here the jump is 2.0*(aver-vFwd)
    //       because Viscous wall use a symmetry value as the Bwd,
    //       not the target one
    Array<OneD, NekDouble> tmpF{npnts, 0.0};
    Array<OneD, NekDouble> tmpB{npnts, 0.0};

    Array<OneD, NekDouble> Fweight{npnts, 2.0};
    Array<OneD, NekDouble> Bweight;
    Bweight = m_tracBwdWeightJump;
    Vmath::Vsub(npnts, Fweight, 1, Bweight, 1, Fweight, 1);

    for (int i = 0; i < nConvectiveFields; ++i)
    {
        Vmath::Vsub(npnts, aver[i], 1, vFwd[i], 1, tmpF, 1);
        Vmath::Vsub(npnts, vBwd[i], 1, aver[i], 1, tmpB, 1);

        Vmath::Vmul(npnts, tmpF, 1, Fweight, 1, tmpF, 1);
        Vmath::Vmul(npnts, tmpB, 1, Bweight, 1, tmpB, 1);
        Vmath::Vadd(npnts, tmpF, 1, tmpB, 1, jump[i], 1);
    }
}

void DiffusionIP::ConsVarAve(
    const std::size_t nConvectiveFields, const size_t npnts,
    const Array<OneD, const Array<OneD, NekDouble>> &vFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &vBwd,
    Array<OneD, Array<OneD, NekDouble>> &aver)
{
    NekDouble LinternalEngy = 0.0;
    NekDouble RinternalEngy = 0.0;
    NekDouble AinternalEngy = 0.0;

    Array<OneD, NekDouble> Fweight{npnts, 1.0};
    Array<OneD, NekDouble> Bweight;

    Bweight = m_tracBwdWeightAver;

    Vmath::Vsub(npnts, Fweight, 1, Bweight, 1, Fweight, 1);

    for (int i = 0; i < nConvectiveFields - 1; ++i)
    {
        Vmath::Vmul(npnts, Fweight, 1, vFwd[i], 1, aver[i], 1);
        Vmath::Vvtvp(npnts, Bweight, 1, vBwd[i], 1, aver[i], 1, aver[i], 1);
    }

    int nengy  = nConvectiveFields - 1;
    int nvelst = 1;
    int nveled = nengy;
    for (int nt = 0; nt < npnts; ++nt)
    {
        LinternalEngy = 0.0;
        for (int j = nvelst; j < nveled; ++j)
        {
            LinternalEngy += vFwd[j][nt] * vFwd[j][nt];
        }
        LinternalEngy *= -0.5 / vFwd[0][nt];
        LinternalEngy += vFwd[nengy][nt];

        RinternalEngy = 0.0;
        for (int j = nvelst; j < nveled; ++j)
        {
            RinternalEngy += vBwd[j][nt] * vBwd[j][nt];
        }
        RinternalEngy *= -0.5 / vBwd[0][nt];
        RinternalEngy += vBwd[nengy][nt];

        AinternalEngy = 0.0;
        aver[nengy][nt] =
            Fweight[nt] * LinternalEngy + Bweight[nt] * RinternalEngy;
        for (int j = nvelst; j < nveled; ++j)
        {
            AinternalEngy += aver[j][nt] * aver[j][nt];
        }
        aver[nengy][nt] += AinternalEngy * (0.5 / aver[0][nt]);
    }
}

void DiffusionIP::CalTraceNumFlux(
    const std::size_t nConvectiveFields, const size_t nDim, const size_t nPts,
    const size_t nTracePts, const NekDouble PenaltyFactor2,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const TensorOfArray3D<NekDouble> &qfield,
    const Array<OneD, Array<OneD, NekDouble>> &vFwd,
    const Array<OneD, Array<OneD, NekDouble>> &vBwd,
    const Array<OneD, NekDouble> &MuVarTrace,
    Array<OneD, int> &nonZeroIndexflux, TensorOfArray3D<NekDouble> &traceflux,
    Array<OneD, Array<OneD, NekDouble>> &solution_Aver,
    Array<OneD, Array<OneD, NekDouble>> &solution_jump)
{
    boost::ignore_unused(inarray);
    const MultiRegions::AssemblyMapDGSharedPtr TraceMap =
        fields[0]->GetTraceMap();

    TensorOfArray3D<NekDouble> numDerivBwd{nDim};
    // Fwd is also used for final numerical results
    TensorOfArray3D<NekDouble> numDerivFwd{nDim};
    for (int nd = 0; nd < nDim; ++nd)
    {
        numDerivBwd[nd] =
            Array<OneD, Array<OneD, NekDouble>>{nConvectiveFields};
        numDerivFwd[nd] =
            Array<OneD, Array<OneD, NekDouble>>{nConvectiveFields};
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            numDerivBwd[nd][i] = Array<OneD, NekDouble>{nTracePts, 0.0};
            numDerivFwd[nd][i] = Array<OneD, NekDouble>{nTracePts, 0.0};
        }
    }

    Array<OneD, NekDouble> Fwd{nTracePts, 0.0};
    Array<OneD, NekDouble> Bwd{nTracePts, 0.0};

    if (fabs(PenaltyFactor2) > 1.0E-12)
    {
        AddSecondDerivToTrace(nConvectiveFields, nDim, nPts, nTracePts,
                              PenaltyFactor2, fields, qfield, numDerivFwd,
                              numDerivBwd);
    }

    for (int nd = 0; nd < nDim; ++nd)
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Zero(nTracePts, Bwd, 1);
            Vmath::Zero(nTracePts, Fwd, 1);
            fields[i]->GetFwdBwdTracePhysDerivSerial(nd, qfield[nd][i], Fwd,
                                                      Bwd);
            Vmath::Svtvp(nTracePts, 0.5, Bwd, 1, numDerivBwd[nd][i], 1,
                         numDerivBwd[nd][i], 1);
            Vmath::Svtvp(nTracePts, 0.5, Fwd, 1, numDerivFwd[nd][i], 1,
                         numDerivFwd[nd][i], 1);
            TraceMap->GetAssemblyCommDG()->PerformExchange(numDerivFwd[nd][i],
                                                           numDerivBwd[nd][i]);
            Vmath::Vadd(nTracePts, numDerivFwd[nd][i], 1, numDerivBwd[nd][i], 1,
                        numDerivFwd[nd][i], 1);
        }
    }

    for (int nd = 0; nd < nDim; ++nd)
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            numDerivBwd[nd][i] = NullNekDouble1DArray;
        }
    }

    ConsVarAveJump(nConvectiveFields, nTracePts, vFwd, vBwd, solution_Aver,
                   solution_jump);

    Array<OneD, NekDouble> jumpTmp       = Fwd;
    Array<OneD, NekDouble> PenaltyFactor = Bwd;
    GetPenaltyFactor_const(fields, PenaltyFactor);

    Vmath::Vmul(nTracePts, PenaltyFactor, 1, m_traceNormDirctnElmtLengthRecip,
                1, PenaltyFactor, 1);
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        Vmath::Vmul(nTracePts, solution_jump[i], 1, PenaltyFactor, 1, jumpTmp,
                    1);
        for (int nd = 0; nd < nDim; ++nd)
        {
            Vmath::Vvtvp(nTracePts, m_traceNormals[nd], 1, jumpTmp, 1,
                         numDerivFwd[nd][i], 1, numDerivFwd[nd][i], 1);
        }
    }
    jumpTmp       = NullNekDouble1DArray;
    PenaltyFactor = NullNekDouble1DArray;

    // Calculate normal viscous flux
    m_FunctorDiffusionfluxCons(nDim, solution_Aver, numDerivFwd, traceflux,
                               nonZeroIndexflux, m_traceNormals, MuVarTrace);
}

void DiffusionIP::AddSecondDerivToTrace(
    const std::size_t nConvectiveFields, const size_t nDim, const size_t nPts,
    const size_t nTracePts, const NekDouble PenaltyFactor2,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const TensorOfArray3D<NekDouble> &qfield,
    TensorOfArray3D<NekDouble> &numDerivFwd,
    TensorOfArray3D<NekDouble> &numDerivBwd)
{
    Array<OneD, NekDouble> Fwd{nTracePts, 0.0};
    Array<OneD, NekDouble> Bwd{nTracePts, 0.0};
    Array<OneD, NekDouble> tmp{nTracePts, 0.0};

    Array<OneD, Array<OneD, NekDouble>> elmt2ndDerv{nDim};
    for (int nd1 = 0; nd1 < nDim; ++nd1)
    {
        elmt2ndDerv[nd1] = Array<OneD, NekDouble>{nPts, 0.0};
    }

    Array<OneD, Array<OneD, NekDouble>> qtmp{3};
    for (int nd = 0; nd < 3; ++nd)
    {
        qtmp[nd] = NullNekDouble1DArray;
    }
    for (int nd2 = 0; nd2 < nDim; ++nd2)
    {
        qtmp[nd2] = elmt2ndDerv[nd2];
    }

    Vmath::Smul(nTracePts, PenaltyFactor2, m_traceNormDirctnElmtLength, 1, tmp,
                1);
    // the derivatives are assumed to be exchangable
    for (int nd1 = 0; nd1 < nDim; ++nd1)
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            fields[i]->PhysDeriv(qfield[nd1][i], qtmp[0], qtmp[1], qtmp[2]);

            for (int nd2 = nd1; nd2 < nDim; ++nd2)
            {
                Vmath::Zero(nTracePts, Bwd, 1);
                fields[i]->GetFwdBwdTracePhysDerivSerial(nd2, elmt2ndDerv[nd2],
                                                         Fwd, Bwd);
                Vmath::Vmul(nTracePts, tmp, 1, Bwd, 1, Bwd, 1);
                Vmath::Vvtvp(nTracePts, m_traceNormals[nd2], 1, Bwd, 1,
                             numDerivBwd[nd1][i], 1, numDerivBwd[nd1][i], 1);
                Vmath::Vmul(nTracePts, tmp, 1, Fwd, 1, Fwd, 1);
                Vmath::Vvtvm(nTracePts, m_traceNormals[nd2], 1, Fwd, 1,
                             numDerivFwd[nd1][i], 1, numDerivFwd[nd1][i], 1);
                Vmath::Neg(nTracePts, numDerivFwd[nd1][i], 1);

                if (nd2 != nd1)
                {
                    Vmath::Vvtvp(nTracePts, m_traceNormals[nd1], 1, Bwd, 1,
                                 numDerivBwd[nd2][i], 1, numDerivBwd[nd2][i],
                                 1);
                    Vmath::Vvtvm(nTracePts, m_traceNormals[nd1], 1, Fwd, 1,
                                 numDerivFwd[nd2][i], 1, numDerivFwd[nd2][i],
                                 1);
                    Vmath::Neg(nTracePts, numDerivFwd[nd2][i], 1);
                }
            }
        }
    }
}
} // namespace SolverUtils
} // namespace Nektar
