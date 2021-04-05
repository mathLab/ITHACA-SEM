///////////////////////////////////////////////////////////////////////////////
//
// File:  PreconCfsBRJ.cpp
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
// Description:  PreconCfsBRJ definition
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/Preconditioner/PreconCfsBRJ.h>

using namespace std;

namespace Nektar
{
/**
 * @class  PreconCfsBRJ
 *
 * Solves a linear system using iterative methods.
 */
std::string PreconCfsBRJ::className =
    GetPreconCfsOpFactory().RegisterCreatorFunction(
        "PreconCfsBRJ", PreconCfsBRJ::create,
        "Block Relaxed Jacobi Preconditioner for CFS.");

PreconCfsBRJ::PreconCfsBRJ(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm)
    : PreconCfsOp(pFields, pSession, vComm)
{
    pSession->LoadParameter("PreconItsStep", m_PreconItsStep, 7);
    pSession->LoadParameter("BRJRelaxParam", m_BRJRelaxParam, 1.0);

    int nvariables     = pFields.size();
    m_PreconMatStorage = eDiagonal;

    m_PreconMatVarsSingle = TensorOfArray2D<SNekBlkMatSharedPtr>(nvariables);
    for (int i = 0; i < nvariables; i++)
    {
        m_PreconMatVarsSingle[i] = Array<OneD, SNekBlkMatSharedPtr>(nvariables);
    }
    AllocatePreconBlkDiagCoeff(pFields, m_PreconMatVarsSingle);

    int nelmts = pFields[0]->GetNumElmts();
    int nelmtcoef;
    Array<OneD, unsigned int> nelmtmatdim(nelmts);
    for (int i = 0; i < nelmts; i++)
    {
        nelmtcoef      = pFields[0]->GetExp(i)->GetNcoeffs();
        nelmtmatdim[i] = nelmtcoef * nvariables;
    }
    AllocateNekBlkMatDig(m_PreconMatSingle, nelmtmatdim, nelmtmatdim);
}

void PreconCfsBRJ::v_InitObject()
{
    PreconCfsOp::v_InitObject();
}

void PreconCfsBRJ::v_DoPreconCfs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray,
    const bool &flag)
{
    boost::ignore_unused(flag);

    int nBRJIterTot = m_PreconItsStep;
    if (0 == nBRJIterTot)
    {
        DoNullPrecon(inarray, outarray, flag);
    }
    else
    {
        const NekDouble BRJParam   = m_BRJRelaxParam;
        const NekDouble OmBRJParam = 1.0 - BRJParam;

        unsigned int nvariables = pFields.size();
        unsigned int npoints    = pFields[0]->GetNcoeffs();
        unsigned int ntotpnt    = inarray.size();

        ASSERTL0(nvariables * npoints == ntotpnt,
                 "nvariables*npoints!=ntotpnt in PreconCoeff");

        Array<OneD, NekDouble> rhs(ntotpnt);

        Array<OneD, NekDouble> outN(ntotpnt);
        Array<OneD, NekDouble> outTmp(ntotpnt);
        Array<OneD, Array<OneD, NekDouble>> rhs2d(nvariables);
        Array<OneD, Array<OneD, NekDouble>> out_2d(nvariables);
        Array<OneD, Array<OneD, NekDouble>> outTmp_2d(nvariables);
        for (int m = 0; m < nvariables; m++)
        {
            int moffset  = m * npoints;
            rhs2d[m]     = rhs + moffset;
            out_2d[m]    = outarray + moffset;
            outTmp_2d[m] = outTmp + moffset;
            pFields[m]->MultiplyByMassMatrix(inarray + moffset, rhs2d[m]);
        }

        int nphysic   = pFields[0]->GetNpoints();
        int nTracePts = pFields[0]->GetTrace()->GetNpoints();
        TensorOfArray3D<NekDouble> qfield(m_spacedim);
        for (int i = 0; i < m_spacedim; i++)
        {
            qfield[i] = Array<OneD, Array<OneD, NekDouble>>(nvariables);
            for (int j = 0; j < nvariables; j++)
            {
                qfield[i][j] = Array<OneD, NekDouble>(nphysic);
            }
        }
        int ntmpTrace = 4 + 2 * m_spacedim;
        TensorOfArray3D<NekDouble> tmpTrace(ntmpTrace);
        for (int i = 0; i < ntmpTrace; i++)
        {
            tmpTrace[i] = Array<OneD, Array<OneD, NekDouble>>(nvariables);
            for (int j = 0; j < nvariables; j++)
            {
                tmpTrace[i][j] = Array<OneD, NekDouble>(nTracePts);
            }
        }
        Array<OneD, Array<OneD, NekDouble>> FwdFluxDeriv(nvariables);
        Array<OneD, Array<OneD, NekDouble>> BwdFluxDeriv(nvariables);
        for (int j = 0; j < nvariables; j++)
        {
            FwdFluxDeriv[j] = Array<OneD, NekDouble>(nTracePts);
            BwdFluxDeriv[j] = Array<OneD, NekDouble>(nTracePts);
        }

        bool flagUpdateDervFlux = false;

        const int nwspTraceDataType = nvariables + 1;
        NekSingle tmpSingle;
        Array<OneD, Array<OneD, NekSingle>> wspTraceDataType(nwspTraceDataType);
        for (int m = 0; m < nwspTraceDataType; m++)
        {
            wspTraceDataType[m] = Array<OneD, NekSingle>(nTracePts);
        }

        PreconBlkDiag(pFields, rhs, outarray, m_PreconMatSingle, tmpSingle);

        for (int nrelax = 0; nrelax < nBRJIterTot - 1; nrelax++)
        {
            Vmath::Smul(ntotpnt, OmBRJParam, outarray, 1, outN, 1);

            MinusOffDiag2Rhs(
                pFields, nvariables, npoints, rhs2d, out_2d, flagUpdateDervFlux,
                FwdFluxDeriv, BwdFluxDeriv, qfield, tmpTrace, wspTraceDataType,
                m_TraceJacArraySingle, m_TraceJacDerivArraySingle,
                m_TraceJacDerivSignSingle, m_TraceIPSymJacArraySingle);

            PreconBlkDiag(pFields, outarray, outTmp, m_PreconMatSingle,
                          tmpSingle);
            Vmath::Svtvp(ntotpnt, BRJParam, outTmp, 1, outN, 1, outarray, 1);
        }
    }
}

void PreconCfsBRJ::v_BuildPreconCfs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, const Array<OneD, NekDouble>> &intmp,
    const NekDouble time, const NekDouble lambda)
{
    if (0 < m_PreconItsStep)
    {
        m_operator.DoCalcPreconMatBRJCoeff(
            intmp, m_PreconMatVarsSingle, m_PreconMatSingle, m_TraceJacSingle,
            m_TraceJacDerivSingle, m_TraceJacDerivSignSingle,
            m_TraceJacArraySingle, m_TraceJacDerivArraySingle,
            m_TraceIPSymJacArraySingle);

        if (m_verbose && m_root)
        {
            cout << "     ## CalcuPreconMat " << endl;
        }
    }

    m_BndEvaluateTime   = time;
    m_DtLambdaPreconMat = lambda;

    m_CalcPreconMatFlag  = false;
    m_PreconTimesCounter = 1;
}

bool PreconCfsBRJ::UpdatePreconMatCheck(const Array<OneD, const NekDouble> &res,
                                        const NekDouble dtLambda)
{
    boost::ignore_unused(res);

    bool flag = false;

    if (m_CalcPreconMatFlag || (m_DtLambdaPreconMat != dtLambda))
    {
        flag = true;
    }

    if (m_PreconMatFreezNumb < m_PreconTimesCounter)
    {
        flag = true;
    }

    m_CalcPreconMatFlag = flag;
    return flag;
}

template <typename DataType, typename TypeNekBlkMatSharedPtr>
void PreconCfsBRJ::PreconBlkDiag(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray,
    const TypeNekBlkMatSharedPtr &PreconMatVars, const DataType &tmpDataType)
{
    boost::ignore_unused(tmpDataType);

    unsigned int nvariables = pFields.size();
    unsigned int npoints    = pFields[0]->GetNcoeffs();
    unsigned int npointsVar = nvariables * npoints;
    Array<OneD, DataType> Sinarray(npointsVar);
    Array<OneD, DataType> Soutarray(npointsVar);
    NekVector<DataType> tmpVect(npointsVar, Sinarray, eWrapper);
    NekVector<DataType> outVect(npointsVar, Soutarray, eWrapper);

    std::shared_ptr<LocalRegions::ExpansionVector> expvect =
        pFields[0]->GetExp();
    int nTotElmt = (*expvect).size();

    for (int m = 0; m < nvariables; m++)
    {
        int nVarOffset = m * npoints;
        for (int ne = 0; ne < nTotElmt; ne++)
        {
            int nCoefOffset = pFields[0]->GetCoeff_Offset(ne);
            int nElmtCoef   = pFields[0]->GetNcoeffs(ne);
            int inOffset    = nVarOffset + nCoefOffset;
            int outOffset   = nCoefOffset * nvariables + m * nElmtCoef;
            for (int i = 0; i < nElmtCoef; i++)
            {
                Sinarray[outOffset + i] = DataType(inarray[inOffset + i]);
            }
        }
    }

    outVect = (*PreconMatVars) * tmpVect;

    for (int m = 0; m < nvariables; m++)
    {
        int nVarOffset = m * npoints;
        for (int ne = 0; ne < nTotElmt; ne++)
        {
            int nCoefOffset = pFields[0]->GetCoeff_Offset(ne);
            int nElmtCoef   = pFields[0]->GetNcoeffs(ne);
            int inOffset    = nVarOffset + nCoefOffset;
            int outOffset   = nCoefOffset * nvariables + m * nElmtCoef;
            for (int i = 0; i < nElmtCoef; i++)
            {
                outarray[inOffset + i] = NekDouble(Soutarray[outOffset + i]);
            }
        }
    }
}

template <typename DataType>
void PreconCfsBRJ::MinusOffDiag2Rhs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const int nvariables, const int nCoeffs,
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, bool flagUpdateDervFlux,
    Array<OneD, Array<OneD, NekDouble>> &FwdFluxDeriv,
    Array<OneD, Array<OneD, NekDouble>> &BwdFluxDeriv,
    TensorOfArray3D<NekDouble> &qfield, TensorOfArray3D<NekDouble> &wspTrace,
    Array<OneD, Array<OneD, DataType>> &wspTraceDataType,
    const TensorOfArray4D<DataType> &TraceJacArray,
    const TensorOfArray4D<DataType> &TraceJacDerivArray,
    const Array<OneD, const Array<OneD, DataType>> &TraceJacDerivSign,
    const TensorOfArray5D<DataType> &TraceIPSymJacArray)
{
    boost::ignore_unused(flagUpdateDervFlux, qfield, TraceJacDerivArray,
                         TraceJacDerivSign, FwdFluxDeriv, BwdFluxDeriv,
                         TraceIPSymJacArray);

    int nTracePts = pFields[0]->GetTrace()->GetNpoints();
    int npoints   = pFields[0]->GetNpoints();
    int nDim      = m_spacedim;

    Array<OneD, Array<OneD, NekDouble>> outpnts(nvariables);
    for (int i = 0; i < nvariables; i++)
    {
        outpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
        pFields[i]->BwdTrans(outarray[i], outpnts[i]);
    }

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd;
    Array<OneD, Array<OneD, NekDouble>> Bwd;
    Array<OneD, Array<OneD, NekDouble>> FwdFlux;
    Array<OneD, Array<OneD, NekDouble>> BwdFlux;
    TensorOfArray3D<NekDouble> numDerivBwd(nDim);
    TensorOfArray3D<NekDouble> numDerivFwd(nDim);
    int indexwspTrace = 0;
    Fwd               = wspTrace[indexwspTrace], indexwspTrace++;
    Bwd               = wspTrace[indexwspTrace], indexwspTrace++;
    FwdFlux           = wspTrace[indexwspTrace], indexwspTrace++;
    BwdFlux           = wspTrace[indexwspTrace], indexwspTrace++;

    for (int i = 0; i < nvariables; ++i)
    {
        pFields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
    }

    int indexwspTraceDataType = 0;
    Array<OneD, Array<OneD, DataType>> Fwdarray(nvariables);
    for (int m = 0; m < nvariables; ++m)
    {
        Fwdarray[m] = wspTraceDataType[indexwspTraceDataType],
        indexwspTraceDataType++;
    }
    Array<OneD, DataType> Fwdreslt;
    Fwdreslt = wspTraceDataType[indexwspTraceDataType], indexwspTraceDataType++;

    for (int m = 0; m < nvariables; ++m)
    {
        for (int i = 0; i < nTracePts; ++i)
        {
            Fwdarray[m][i] = DataType(Fwd[m][i]);
        }
    }
    for (int m = 0; m < nvariables; ++m)
    {
        Vmath::Zero(nTracePts, Fwdreslt, 1);
        for (int n = 0; n < nvariables; ++n)
        {
            Vmath::Vvtvp(nTracePts, TraceJacArray[0][m][n], 1, Fwdarray[n], 1,
                         Fwdreslt, 1, Fwdreslt, 1);
        }

        for (int i = 0; i < nTracePts; ++i)
        {
            FwdFlux[m][i] = NekDouble(Fwdreslt[i]);
        }
    }

    for (int m = 0; m < nvariables; ++m)
    {
        for (int i = 0; i < nTracePts; ++i)
        {
            Fwdarray[m][i] = DataType(Bwd[m][i]);
        }
    }
    for (int m = 0; m < nvariables; ++m)
    {
        Vmath::Zero(nTracePts, Fwdreslt, 1);
        for (int n = 0; n < nvariables; ++n)
        {
            Vmath::Vvtvp(nTracePts, TraceJacArray[1][m][n], 1, Fwdarray[n], 1,
                         Fwdreslt, 1, Fwdreslt, 1);
        }
        for (int i = 0; i < nTracePts; ++i)
        {
            BwdFlux[m][i] = NekDouble(Fwdreslt[i]);
        }
    }

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Fill(nCoeffs, 0.0, outarray[i], 1);
        pFields[i]->AddTraceIntegralToOffDiag(FwdFlux[i], BwdFlux[i],
                                              outarray[i]);
    }

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Svtvp(nCoeffs, -m_DtLambdaPreconMat, outarray[i], 1, inarray[i],
                     1, outarray[i], 1);
    }
}

template <typename TypeNekBlkMatSharedPtr>
void PreconCfsBRJ::AllocatePreconBlkDiagCoeff(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>> &gmtxarray,
    const int &nscale)
{

    int nvars  = pFields.size();
    int nelmts = pFields[0]->GetNumElmts();
    int nelmtcoef;
    Array<OneD, unsigned int> nelmtmatdim(nelmts);
    for (int i = 0; i < nelmts; i++)
    {
        nelmtcoef      = pFields[0]->GetExp(i)->GetNcoeffs();
        nelmtmatdim[i] = nelmtcoef * nscale;
    }

    for (int i = 0; i < nvars; i++)
    {
        for (int j = 0; j < nvars; j++)
        {
            AllocateNekBlkMatDig(gmtxarray[i][j], nelmtmatdim, nelmtmatdim);
        }
    }
}

} // namespace Nektar
