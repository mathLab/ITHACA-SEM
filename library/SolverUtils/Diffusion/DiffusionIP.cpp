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

#include <LibUtilities/BasicUtils/Timer.h>

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
    pFields[0]->PeriodicBwdCopy(lengthFwd, lengthBwd);

    const MultiRegions::AssemblyMapDGSharedPtr TraceMap =
        pFields[0]->GetTraceMap();
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

    // workspace for v_diffuse
    size_t nCoeffs = pFields[0]->GetNcoeffs();
    m_wspDiff = Array<OneD, Array<OneD, NekDouble>>{nVariable};
    for (int i = 0; i < nVariable; ++i)
    {
        m_wspDiff[i] = Array<OneD, NekDouble>{nCoeffs, 0.0};
    }

    // workspace for callnumtraceflux
    m_wspNumDerivBwd = TensorOfArray3D<NekDouble>{nDim};
    m_wspNumDerivFwd = TensorOfArray3D<NekDouble>{nDim};
    for (int nd = 0; nd < nDim; ++nd)
    {
        m_wspNumDerivBwd[nd] =
            Array<OneD, Array<OneD, NekDouble>>{nVariable};
        m_wspNumDerivFwd[nd] =
            Array<OneD, Array<OneD, NekDouble>>{nVariable};
        for (int i = 0; i < nVariable; ++i)
        {
            m_wspNumDerivBwd[nd][i] = Array<OneD, NekDouble>{nTracePts, 0.0};
            m_wspNumDerivFwd[nd][i] = Array<OneD, NekDouble>{nTracePts, 0.0};
        }
    }
}

void DiffusionIP::v_Diffuse(
    const size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{

    LibUtilities::Timer timer;
    timer.Start();
    DiffusionIP::v_DiffuseCoeffs(nConvectiveFields, fields, inarray, m_wspDiff,
        pFwd, pBwd);
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        fields[i]->BwdTrans(m_wspDiff[i], outarray[i]);
    }
    timer.Stop();
    timer.AccumulateRegion("Diffusion IP");
}

void DiffusionIP::v_DiffuseCoeffs(
    const size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>>       &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    LibUtilities::Timer timer;
    timer.Start();

    size_t nDim      = fields[0]->GetCoordim(0);
    size_t nPts      = fields[0]->GetTotPoints();
    size_t nCoeffs   = fields[0]->GetNcoeffs();
    size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();

    // pre-allocate this?
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

    // pre-allocate this?
    Array<OneD, Array<OneD, NekDouble>> vFwd{nConvectiveFields};
    Array<OneD, Array<OneD, NekDouble>> vBwd{nConvectiveFields};
    // when does this happen?
    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
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

    DiffuseCalcDerivative(fields, inarray, qfield, vFwd, vBwd);
    
    Array<OneD, int> nonZeroIndex;
    DiffuseVolumeFlux(fields, inarray, qfield, elmtFlux, nonZeroIndex);

    // pre-allocate this?
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
        elmtFlux[j] = NullNekDoubleArrayOfArray;
    }

    // pre-allocate this?
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
    timer.Stop();
    timer.AccumulateRegion("DiffIP:Diffusion Coeff:");
}

void DiffusionIP::v_DiffuseCoeffs(
        const std::size_t nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, Array<OneD, NekDouble>> &vFwd,
        const Array<OneD, Array<OneD, NekDouble>> &vBwd,
        TensorOfArray3D<NekDouble>                &qfield,
        Array< OneD, int >                        &nonZeroIndex)
{
    int i, j;
    int nDim      = fields[0]->GetCoordim(0);
    int nPts      = fields[0]->GetTotPoints();
    int nCoeffs   = fields[0]->GetNcoeffs();
    int nTracePts = fields[0]->GetTrace()->GetTotPoints();
    
    TensorOfArray3D<NekDouble> elmtFlux(nDim);
    for (j = 0; j < nDim; ++j)
    {
        elmtFlux[j] = Array<OneD, Array<OneD, NekDouble>>(nConvectiveFields);
        for (i = 0; i < nConvectiveFields; ++i)
        {
            elmtFlux[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
        }
    }
    
    DiffuseVolumeFlux(fields,inarray,qfield,elmtFlux,nonZeroIndex);
    
    Array<OneD, Array<OneD, NekDouble>> tmpFluxIprdct(nDim);
    // volume intergration: the nonZeroIndex indicates which flux is nonzero
    for(i = 0; i < nonZeroIndex.size(); ++i)
    {
        int j = nonZeroIndex[i];
        for (int k = 0; k < nDim; ++k)
        {
            tmpFluxIprdct[k] = elmtFlux[k][j];
        }
        fields[j]->IProductWRTDerivBase(tmpFluxIprdct,outarray[j]);
        Vmath::Neg                      (nCoeffs, outarray[j], 1);
    }
    
    Array<OneD, Array<OneD, NekDouble > > Traceflux(nConvectiveFields);
    for (int j = 0; j < nConvectiveFields; ++j)
    {
        Traceflux[j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
    }
    
    DiffuseTraceFlux(fields,inarray,qfield,elmtFlux,
                     Traceflux,vFwd,vBwd,nonZeroIndex);
    
    // release qfield, elmtFlux and muvar;
    for (j = 0; j < nDim; ++j)
    {
        elmtFlux[j]     = NullNekDoubleArrayOfArray;
    }
    
    for(i = 0; i < nonZeroIndex.size(); ++i)
    {
        int j = nonZeroIndex[i];
        
        fields[j]->AddTraceIntegral     (Traceflux[j], outarray[j]);
        fields[j]->SetPhysState         (false);
    }

    AddDiffusionSymmFluxToCoeff(nConvectiveFields, fields, inarray,qfield,
                                elmtFlux, outarray, vFwd, vBwd);    
}

void DiffusionIP::v_DiffuseCalcDerivative(
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

    Array<OneD, Array<OneD, NekDouble>> tmparray2D = NullNekDoubleArrayOfArray;

    LibUtilities::Timer timer;
    timer.Start();
    m_FunctorDiffusionfluxCons(nDim, inarray, qfield, VolumeFlux, nonZeroIndex,
                               tmparray2D, muvar);
    timer.Stop();
    timer.AccumulateRegion("DiffIP:_FunctorDiffFluxCons",1);
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

    TensorOfArray3D<NekDouble> traceflux3D(1);
    traceflux3D[0] = TraceFlux;

    size_t nConvectiveFields = fields.size();

    LibUtilities::Timer timer;
    timer.Start();
    CalcTraceNumFlux(m_IP2ndDervCoeff,
                    fields, inarray, qfield, pFwd, pBwd, m_MuVarTrace,
                    nonZeroIndex, traceflux3D, m_traceAver, m_traceJump);
    timer.Stop();
    timer.AccumulateRegion("DiffIP:_CalcTraceNumFlux",1);

    ApplyFluxBndConds(nConvectiveFields, fields, TraceFlux);
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
    const Array<OneD, Array<OneD, NekDouble>>         &inarray,
    TensorOfArray3D<NekDouble>                        &qfield,
    TensorOfArray3D<NekDouble>                        &VolumeFlux,
    Array<OneD, Array<OneD, NekDouble>>               &outarray,
    const Array<OneD, Array<OneD, NekDouble>>         &pFwd,
    const Array<OneD, Array<OneD, NekDouble>>         &pBwd)
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

    CalcTraceSymFlux(nConvectiveFields, nDim, fields, m_traceAver, m_traceJump,
                    nonZeroIndex, SymmFlux);
}

void DiffusionIP::CalcTraceSymFlux(
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

void DiffusionIP::GetPenaltyFactorConst(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    Array<OneD, NekDouble> &factor)
{
    boost::ignore_unused(fields);
    Vmath::Fill(factor.size(), m_IPPenaltyCoeff, factor, 1);
}

void DiffusionIP::v_ConsVarAveJump(
    const std::size_t nConvectiveFields, const size_t nPts,
    const Array<OneD, const Array<OneD, NekDouble>> &vFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &vBwd,
    Array<OneD, Array<OneD, NekDouble>> &aver,
    Array<OneD, Array<OneD, NekDouble>> &jump)
{
    // ConsVarAve(nConvectiveFields, nPts, vFwd, vBwd, aver);

    std::vector<NekDouble> vFwdTmp(nConvectiveFields),
        vBwdTmp(nConvectiveFields), averTmp(nConvectiveFields) ;
    for (size_t p = 0; p < nPts; ++p)
    {
        // re-arrange data
        for (size_t f = 0; f < nConvectiveFields; ++f)
        {
            vFwdTmp[f] = vFwd[f][p];
            vBwdTmp[f] = vBwd[f][p];
        }

        ConsVarAve(nConvectiveFields, m_tracBwdWeightAver[p], vFwdTmp, vBwdTmp,
            averTmp);

        // store
        for (size_t f = 0; f < nConvectiveFields; ++f)
        {
            aver[f][p] = averTmp[f];
        }
    }

    // if this can be make function of points, the nPts loop can be lifted more
    m_SpecialBndTreat(aver);

    // note: here the jump is 2.0*(aver-vFwd)
    //       because Viscous wall use a symmetry value as the Bwd,
    //       not the target one

    for (size_t f = 0; f < nConvectiveFields; ++f)
    {
        for (size_t p = 0; p < nPts; ++p)
        {
            NekDouble Bweight = m_tracBwdWeightJump[p];
            NekDouble Fweight = 2.0 - Bweight;

            NekDouble tmpF = aver[f][p] - vFwd[f][p];
            NekDouble tmpB = vBwd[f][p] - aver[f][p];
            jump[f][p] = tmpF * Fweight + tmpB * Bweight;
        }
    }
}

void DiffusionIP::CalcTraceNumFlux(
    const NekDouble PenaltyFactor2,
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
    size_t nDim      = fields[0]->GetCoordim(0);
    size_t nPts      = fields[0]->GetTotPoints();
    size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();
    size_t nConvectiveFields = fields.size();

    boost::ignore_unused(inarray);
    const MultiRegions::AssemblyMapDGSharedPtr TraceMap =
        fields[0]->GetTraceMap();


    LibUtilities::Timer timer;
    timer.Start();
    // with further restructuring this iniziatilization could be eliminated
    for (int nd = 0; nd < nDim; ++nd)
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Zero(nTracePts, m_wspNumDerivBwd[nd][i], 1);
            Vmath::Zero(nTracePts, m_wspNumDerivFwd[nd][i], 1);
        }
    }

    // could this be pre-allocated?
    Array<OneD, NekDouble> Fwd{nTracePts, 0.0};
    Array<OneD, NekDouble> Bwd{nTracePts, 0.0};

    timer.Stop();
    timer.AccumulateRegion("DiffIP:_CalcTraceNumFlux_alloc",1);


    timer.Start();
    if (fabs(PenaltyFactor2) > 1.0E-12)
    {
        AddSecondDerivToTrace(nConvectiveFields, nDim, nPts, nTracePts,
                              PenaltyFactor2, fields, qfield, m_wspNumDerivFwd,
                              m_wspNumDerivBwd);
    }

    timer.Stop();
    timer.AccumulateRegion("DiffIP:_AddSecondDerivToTrace",1);

    for (int nd = 0; nd < nDim; ++nd)
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            // this sequence of operations is really exepensive,
            // it should be done collectively for all fields together instead of
            // one by one
            timer.Start();
            fields[i]->GetFwdBwdTracePhys(qfield[nd][i], Fwd, Bwd, true, true,
                false);
            timer.Stop();
            timer.AccumulateRegion("DiffIP:_GetFwdBwdTracePhys",1);

            for (size_t p = 0; p < nTracePts; ++p)
            {
                m_wspNumDerivBwd[nd][i][p] += 0.5 * Bwd[p];
                m_wspNumDerivFwd[nd][i][p] += 0.5 * Fwd[p];
            }

            timer.Start();
            TraceMap->GetAssemblyCommDG()->PerformExchange(m_wspNumDerivFwd[nd][i],
                                                           m_wspNumDerivBwd[nd][i]);
            timer.Stop();
            timer.AccumulateRegion("DiffIP:_PerformExchange",1);

            Vmath::Vadd(nTracePts, m_wspNumDerivFwd[nd][i], 1, m_wspNumDerivBwd[nd][i], 1,
                        m_wspNumDerivFwd[nd][i], 1);
        }
    }

    timer.Start();
    ConsVarAveJump(nConvectiveFields, nTracePts, vFwd, vBwd, solution_Aver,
                   solution_jump);
    timer.Stop();
    timer.AccumulateRegion("DiffIP:_ConsVarAveJump",1);

    for (size_t p = 0; p < nTracePts; ++p)
    {
        NekDouble PenaltyFactor = m_IPPenaltyCoeff * m_traceNormDirctnElmtLengthRecip[p]; // load 1x

        for (size_t f = 0; f < nConvectiveFields; ++f)
        {
            NekDouble jumpTmp = solution_jump[f][p] * PenaltyFactor; // load 1x
            for (size_t d = 0; d < nDim; ++d)
            {
                NekDouble tmp = m_traceNormals[d][p] * jumpTmp + m_wspNumDerivFwd[d][f][p]; // load 2x
                m_wspNumDerivFwd[d][f][p] = tmp; // store 1x
            }
        }
    }

    timer.Start();
    // Calculate normal viscous flux
    m_FunctorDiffusionfluxConsTrace(nDim, solution_Aver, m_wspNumDerivFwd, traceflux,
                               nonZeroIndexflux, m_traceNormals, MuVarTrace);
    timer.Stop();
    timer.AccumulateRegion("DiffIP:_FunctorDiffFluxConsTrace",1);
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
                fields[i]->GetFwdBwdTracePhys(elmt2ndDerv[nd2], Fwd, Bwd,
                                              true, true, false);
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


/**
 * @brief aplly Neuman boundary conditions on flux 
 *        Currently only consider WallAdiabatic
 *
 **/
void DiffusionIP::ApplyFluxBndConds(
    const int                                               nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
    Array<OneD,       Array<OneD, NekDouble> >              &flux)
{            
    int nengy       = nConvectiveFields-1;
    int cnt;
    // Compute boundary conditions  for Energy
    cnt = 0;
    size_t nBndRegions = fields[nengy]->
    GetBndCondExpansions().size();
    for (size_t j = 0; j < nBndRegions; ++j)
    {
        if (fields[nengy]->GetBndConditions()[j]->
            GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        size_t nBndEdges = fields[nengy]->
        GetBndCondExpansions()[j]->GetExpSize();
        for (size_t e = 0; e < nBndEdges; ++e)
        {
            size_t nBndEdgePts = fields[nengy]->
            GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

            std::size_t id2 = fields[0]->GetTrace()->
                GetPhys_Offset(fields[0]->GetTraceMap()->
                                GetBndCondIDToGlobalTraceID(cnt++));

            if (fields[0]->GetBndConditions()[j]->
                GetUserDefined() =="WallAdiabatic")
            {
                Vmath::Zero(nBndEdgePts, &flux[nengy][id2], 1);
            }                    
        }
    }
}

} // namespace SolverUtils
} // namespace Nektar
