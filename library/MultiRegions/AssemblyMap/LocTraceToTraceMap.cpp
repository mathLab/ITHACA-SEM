///////////////////////////////////////////////////////////////////////////////
//
// File LocTraceToTraceMap.cpp
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
// Description: Local trace to general trace mapping information
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <MultiRegions/AssemblyMap/AssemblyMap.h>
#include <MultiRegions/AssemblyMap/AssemblyMap.h>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

using namespace std;

namespace Nektar
{
namespace MultiRegions
{

/**
 * @brief Set up trace to trace mapping components.
 *
 * @param locExp         Expansion list of full dimension problem.
 * @param trace          Expansion list of one dimension lower trace.
 * @param elmtToTrace    Mapping from elemental facets to trace.
 * @param leftAdjacents  Vector of bools denoting forwards-oriented traces.
 *
 * @todo Add 1D support
 */
LocTraceToTraceMap::LocTraceToTraceMap(
    const ExpList          &locExp,
    const ExpListSharedPtr &trace,
    const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
        &elmtToTrace,
    const vector<bool>     &LeftAdjacents)
{
    const LocalRegions::ExpansionVector &locExpVector = *(locExp.GetExp());

    // Assume that all the elements have same dimension
    m_expdim = locExpVector[0]->GetShapeDimension();

    // set up interpolation details for all dimension elements.
    Setup(locExp, trace, elmtToTrace, LeftAdjacents);
}

LocTraceToTraceMap::~LocTraceToTraceMap()
{
}

/**
 * @brief Set up member variables for a two-dimensional problem.
 *
 * @param locExp         Expansion list of elements
 * @param trace          Expansion list of the trace.
 * @param elmtToTrace    Mapping from elemental trace to unique trace.
 * @param leftAdjacents  Vector of bools denoting forwards-oriented traces.
 */
void LocTraceToTraceMap::Setup(
    const ExpList &locExp,
    const ExpListSharedPtr &trace,
    const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
        &elmtToTrace,
    const vector<bool> &LeftAdjacents)
{
    m_LocTraceToTraceMap = Array<OneD, Array<OneD, int> >(2);
    m_interpTrace        = Array<OneD, Array<OneD, InterpLocTraceToTrace> >(2);
    m_interpTraceI0      = Array<OneD, Array<OneD, DNekMatSharedPtr> >(2);
    m_interpEndPtI0      = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(2);
    m_interpPoints       = Array<OneD, Array<OneD, TraceInterpPoints> >(2);
    m_interpNfaces       = Array<OneD, Array<OneD, int> >(2);

    if(m_expdim == 3)
    {
        m_interpTraceI1 = Array<OneD, Array<OneD, DNekMatSharedPtr> >(2);
        m_interpEndPtI1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(2);
    }

    m_traceCoeffsToElmtMap   = Array<OneD, Array<OneD, int> >(2);
    m_traceCoeffsToElmtTrace = Array<OneD, Array<OneD, int> >(2);
    m_traceCoeffsToElmtSign  = Array<OneD, Array<OneD, int> >(2);

    LocalRegions::ExpansionSharedPtr elmt;
    const std::shared_ptr<LocalRegions::ExpansionVector> exp =
        locExp.GetExp();

    int cnt, n, e, phys_offset;

    int nexp    = exp->size();
    m_nTracePts = trace->GetTotPoints();

    // Count number of traces and points required for maps
    int nFwdPts       = 0;
    int nBwdPts       = 0;
    int nFwdCoeffs    = 0;
    int nBwdCoeffs    = 0;
    m_nFwdLocTracePts = 0;
    m_nLocTracePts    = 0;

    for (cnt = n = 0; n < nexp; ++n)
    {
        elmt = (*exp)[n];

        for (int i = 0; i < elmt->GetNtraces(); ++i, ++cnt)
        {
            int nLocPts = elmt->GetTraceNumPoints(i);
            m_nLocTracePts += nLocPts;

            if (LeftAdjacents[cnt])
            {
                nFwdPts += elmtToTrace[n][i]->GetTotPoints();
                nFwdCoeffs += elmtToTrace[n][i]->GetNcoeffs();
                m_nFwdLocTracePts += nLocPts;
            }
            else
            {
                nBwdPts += elmtToTrace[n][i]->GetTotPoints();
                nBwdCoeffs += elmtToTrace[n][i]->GetNcoeffs();
            }
        }
    }

    m_fieldToLocTraceMap = Array<OneD, int>(m_nLocTracePts);

    m_LocTraceToTraceMap[0] = Array<OneD, int>(nFwdPts);
    m_LocTraceToTraceMap[1] = Array<OneD, int>(nBwdPts);

    m_nTraceCoeffs[0] = nFwdCoeffs;
    m_nTraceCoeffs[1] = nBwdCoeffs;

    m_traceCoeffsToElmtMap[0]   = Array<OneD, int>(nFwdCoeffs + nBwdCoeffs);
    m_traceCoeffsToElmtMap[1]   = m_traceCoeffsToElmtMap[0] + nFwdCoeffs;
    m_traceCoeffsToElmtTrace[0] = Array<OneD, int>(nFwdCoeffs + nBwdCoeffs);
    m_traceCoeffsToElmtTrace[1] = m_traceCoeffsToElmtTrace[0] + nFwdCoeffs;
    m_traceCoeffsToElmtSign[0]  = Array<OneD, int>(nFwdCoeffs + nBwdCoeffs);
    m_traceCoeffsToElmtSign[1]  = m_traceCoeffsToElmtSign[0] + nFwdCoeffs;

    // Gather information about trace interpolations
    map<TraceInterpPoints, vector<pair<int, int> >, cmpop> TraceInterpMap;

    vector<vector<int> > TraceOrder;
    TraceOrder.resize(nexp);
    int ntrace;
    int fwdcnt = 0;
    int bwdcnt = 0;

    // Generate a map of similar traces with the same
    // interpolation requirements
    for (cnt = n = 0; n < nexp; ++n)
    {
        elmt = (*exp)[n];
        ntrace = elmt->GetNtraces();
        TraceOrder[n].resize(ntrace);

        int coeffoffset = locExp.GetCoeff_Offset(n);
        for (e = 0; e < ntrace; ++e, ++cnt)
        {
            LocalRegions::ExpansionSharedPtr elmttrace = elmtToTrace[n][e];
            StdRegions::Orientation orient = elmt->GetTraceOrient(e);

            LibUtilities::PointsKey fromPointsKey0, fromPointsKey1;
            LibUtilities::PointsKey toPointsKey0,   toPointsKey1;
            Array<OneD, int> P(2,-1);

            switch(m_expdim)
            {
            case 1:
                {
                    fromPointsKey0 = elmt->GetBasis(0)->GetPointsKey();
                    fromPointsKey1 =
                        LibUtilities::PointsKey(0, LibUtilities::eNoPointsType);
                    // dummy info since no interpolation is required in this case.
                    toPointsKey0 =
                        LibUtilities::PointsKey(0, LibUtilities::eNoPointsType);
                    toPointsKey1 =
                        LibUtilities::PointsKey(0, LibUtilities::eNoPointsType);
                }
                break;
            case 2:
                {
                    int dir0 = elmt->GetGeom()->GetDir(e, 0);

                    fromPointsKey0 = elmt->GetBasis(dir0)->GetPointsKey();
                    fromPointsKey1 =
                        LibUtilities::PointsKey(0, LibUtilities::eNoPointsType);

                    toPointsKey0 = elmttrace->GetBasis(0)->GetPointsKey();
                    toPointsKey1 =
                        LibUtilities::PointsKey(0, LibUtilities::eNoPointsType);

                    P[0] = elmttrace->GetBasisNumModes(0);
                }
                break;
            case 3:
                {
                    int dir0 = elmt->GetGeom()->GetDir(e, 0);
                    int dir1 = elmt->GetGeom()->GetDir(e, 1);

                    fromPointsKey0 = elmt->GetBasis(dir0)->GetPointsKey();
                    fromPointsKey1 = elmt->GetBasis(dir1)->GetPointsKey();

                    if (orient < StdRegions::eDir1FwdDir2_Dir2FwdDir1)
                    {
                        toPointsKey0 = elmttrace->GetBasis(0)->GetPointsKey();
                        toPointsKey1 = elmttrace->GetBasis(1)->GetPointsKey();
                    }
                    else // transpose points key evaluation
                    {
                        toPointsKey0 = elmttrace->GetBasis(1)->GetPointsKey();
                        toPointsKey1 = elmttrace->GetBasis(0)->GetPointsKey();
                    }

                    P[0] = elmttrace->GetBasisNumModes(0);
                    P[1] = elmttrace->GetBasisNumModes(1);
                }
                break;
            }

            TraceInterpPoints fpoint(fromPointsKey0, fromPointsKey1,
                                      toPointsKey0,   toPointsKey1);

            pair<int, int> epf(n, e);
            TraceInterpMap[fpoint].push_back(epf);
            TraceOrder[n][e] = cnt;

            // Setup for coefficient mapping from trace normal flux
            // to elements
            Array<OneD, unsigned int> map;
            Array<OneD, int> sign;

            elmt->GetTraceToElementMap(e, map, sign, orient, P[0], P[1]);

            int order_t  = elmttrace->GetNcoeffs();
            int t_offset = trace->GetCoeff_Offset(elmttrace->GetElmtId());

            double fac = 1.0;

            if (elmt->GetTraceExp(e)->GetRightAdjacentElementExp())
            {
                if (elmttrace->GetRightAdjacentElementExp()
                    ->GetGeom()->GetGlobalID() == elmt->GetGeom()
                    ->GetGlobalID())
                {
                    fac = -1.0;
                }
            }

            if (LeftAdjacents[cnt])
            {
                for (int i = 0; i < order_t; ++i)
                {
                    m_traceCoeffsToElmtMap[0][fwdcnt]    = coeffoffset + map[i];
                    m_traceCoeffsToElmtTrace[0][fwdcnt]  = t_offset + i;
                    m_traceCoeffsToElmtSign[0][fwdcnt++] = fac * sign[i];
                }
            }
            else
            {
                for (int i = 0; i < order_t; ++i)
                {
                    m_traceCoeffsToElmtMap[1][bwdcnt]    = coeffoffset + map[i];
                    m_traceCoeffsToElmtTrace[1][bwdcnt]  = t_offset + i;
                    m_traceCoeffsToElmtSign[1][bwdcnt++] = fac * sign[i];
                }
            }
        }
    }

    int nInterpType = TraceInterpMap.size();

    // need to decide on 1D case here !!!!!
    for (int i = 0; i < 2; ++i)
    {
        m_interpTrace[i]   = Array<OneD, InterpLocTraceToTrace>(nInterpType);
        m_interpTraceI0[i] = Array<OneD, DNekMatSharedPtr>(nInterpType);
        m_interpEndPtI0[i] = Array<OneD, Array<OneD, NekDouble> >(nInterpType);
        m_interpPoints[i]  = Array<OneD, TraceInterpPoints>(nInterpType);
        m_interpNfaces[i]  = Array<OneD, int>(nInterpType, 0);
    }

    if(m_expdim > 2)
    {
        for (int i = 0; i < 2; ++i)
        {
            m_interpTraceI1[i] = Array<OneD, DNekMatSharedPtr>(nInterpType);
            m_interpEndPtI1[i] = Array<OneD, Array<OneD, NekDouble> >
                (nInterpType);
        }
    }

    int ntracepts, ntracepts1;
    int cnt1    = 0;
    int cnt2    = 0;
    int cntFwd  = 0;
    int cntBwd  = 0;
    int cntFwd1 = 0;
    int cntBwd1 = 0;
    int set;
    Array<OneD, int> traceids;
    Array<OneD, int> locTraceToTraceMap;
    cnt = 0;

    for (auto it = TraceInterpMap.begin(); it != TraceInterpMap.end();
         ++it, ++cnt1)
    {
        LibUtilities::PointsKey fromPointsKey0 = std::get<0>(it->first);
        LibUtilities::PointsKey fromPointsKey1 = std::get<1>(it->first);
        LibUtilities::PointsKey toPointsKey0   = std::get<2>(it->first);
        LibUtilities::PointsKey toPointsKey1   = std::get<3>(it->first);

        bool fwdSet = false;
        bool bwdSet = false;

        for (int f = 0; f < it->second.size(); ++f, ++cnt2)
        {
            n = it->second[f].first;
            e = it->second[f].second;

            StdRegions::StdExpansionSharedPtr elmttrace = elmtToTrace[n][e];

            elmt        = (*exp)[n];
            phys_offset = locExp.GetPhys_Offset(n);

            // Mapping of new edge order to one that loops over elmts
            // then set up mapping of faces in standard cartesian order
            elmt->GetTracePhysMap(e, traceids);

            ntracepts  = elmt->GetTraceNumPoints(e);
            ntracepts1 = elmttrace->GetTotPoints();

            StdRegions::Orientation orient = elmt->GetTraceOrient(e);

            elmt->ReOrientTracePhysMap(orient, locTraceToTraceMap,
                                       toPointsKey0.GetNumPoints(),
                                       toPointsKey1.GetNumPoints());

            int offset = trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

            if (LeftAdjacents[TraceOrder[n][e]])
            {
                for (int i = 0; i < ntracepts; ++i)
                {
                    m_fieldToLocTraceMap[cntFwd + i] =
                        phys_offset + traceids[i];
                }

                for (int i = 0; i < ntracepts1; ++i)
                {
                    m_LocTraceToTraceMap[0][cntFwd1 + i] =
                        offset + locTraceToTraceMap[i];
                }

                cntFwd += ntracepts;
                cntFwd1 += ntracepts1;
                set = 0;
            }
            else
            {
                for (int i = 0; i < ntracepts; ++i)
                {
                    m_fieldToLocTraceMap[m_nFwdLocTracePts + cntBwd + i] =
                        phys_offset + traceids[i];
                }

                for (int i = 0; i < ntracepts1; ++i)
                {
                    m_LocTraceToTraceMap[1][cntBwd1 + i] =
                        offset + locTraceToTraceMap[i];
                }

                cntBwd += ntracepts;
                cntBwd1 += ntracepts1;
                set = 1;
            }

            m_interpNfaces[set][cnt1] += 1;

            if ((fwdSet == false && set == 0) ||
                (bwdSet == false && set == 1))
            {
                m_interpPoints[set][cnt1] = it->first;

                switch(m_expdim)
                {
                case 1:
                    {
                        // Always no interplation in this case
                        m_interpTrace[set][cnt1] = eNoInterp;
                    }
                    break;
                case 2:
                    {
                        if (fromPointsKey0 == toPointsKey0)
                        {
                            m_interpTrace[set][cnt1] = eNoInterp;
                        }
                        else
                        {
                            m_interpTrace[set][cnt1] = eInterpDir0;
                            m_interpTraceI0[set][cnt1] =
                                LibUtilities::PointsManager()
                                [fromPointsKey0]->GetI(toPointsKey0);

                            // Check to see if we can
                            // just interpolate endpoint
                            if ((fromPointsKey0.GetPointsType() ==
                                 LibUtilities::eGaussRadauMAlpha1Beta0) &&
                                (toPointsKey0.GetPointsType() ==
                                 LibUtilities::eGaussLobattoLegendre))
                            {
                                if (fromPointsKey0.GetNumPoints() + 1 ==
                                    toPointsKey0.GetNumPoints())
                                {
                                    m_interpTrace[set][cnt1] = eInterpEndPtDir0;

                                    int fnp0 = fromPointsKey0.GetNumPoints();
                                    int tnp0 = toPointsKey0.GetNumPoints();

                                    m_interpEndPtI0[set][cnt1] =
                                        Array<OneD, NekDouble>(fnp0);

                                    Vmath::Vcopy
                                        (fnp0,
                                         m_interpTraceI0[set][cnt1]->
                                         GetPtr().get() +  tnp0 - 1, tnp0,
                                         &m_interpEndPtI0[set][cnt1][0],1);
                                }
                            }
                        }
                    }
                    break;
                case 3:
                    {
                        if (fromPointsKey0 == toPointsKey0)
                        {
                            if (fromPointsKey1 == toPointsKey1)
                            {
                                m_interpTrace[set][cnt1] = eNoInterp;
                            }
                            else
                            {
                                m_interpTrace[set][cnt1] = eInterpDir1;
                                m_interpTraceI1[set][cnt1] =
                                    LibUtilities::PointsManager()
                                    [fromPointsKey1]->GetI(toPointsKey1);

                                // Check to see if we can just
                                // interpolate endpoint
                                if ((fromPointsKey1.GetPointsType() ==
                                     LibUtilities::eGaussRadauMAlpha1Beta0) &&
                                    (toPointsKey1.GetPointsType() ==
                                     LibUtilities::eGaussLobattoLegendre))
                                {
                                    if (fromPointsKey1.GetNumPoints() + 1 ==
                                        toPointsKey1.GetNumPoints())
                                    {
                                        m_interpTrace[set][cnt1] = eInterpEndPtDir1;
                                        int fnp1 = fromPointsKey1.GetNumPoints();
                                        int tnp1 = toPointsKey1.GetNumPoints();
                                        m_interpEndPtI1[set][cnt1] =
                                            Array<OneD, NekDouble>(fnp1);
                                        Vmath::Vcopy
                                            (fnp1,
                                             m_interpTraceI1[set][cnt1]->GetPtr().get()+
                                             tnp1 - 1, tnp1,
                                             &m_interpEndPtI1[set][cnt1][0], 1);
                                    }
                                }
                            }
                        }
                        else
                        {
                            if (fromPointsKey1 == toPointsKey1)
                            {
                                m_interpTrace[set][cnt1] = eInterpDir0;
                                m_interpTraceI0[set][cnt1] =
                                    LibUtilities::PointsManager()
                                    [fromPointsKey0]->GetI(toPointsKey0);

                                // Check to see if we can just
                                // interpolate endpoint
                                if ((fromPointsKey0.GetPointsType() ==
                                     LibUtilities::eGaussRadauMAlpha1Beta0) &&
                                    (toPointsKey0.GetPointsType() ==
                                     LibUtilities::eGaussLobattoLegendre))
                                {
                                    if (fromPointsKey0.GetNumPoints() + 1 ==
                                        toPointsKey0.GetNumPoints())
                                    {
                                        m_interpTrace[set][cnt1] =
                                            eInterpEndPtDir0;
                                        int fnp0 = fromPointsKey0.GetNumPoints();
                                        int tnp0 = toPointsKey0.GetNumPoints();
                                        m_interpEndPtI0[set][cnt1] =
                                            Array<OneD, NekDouble>(fnp0);
                                        Vmath::Vcopy
                                            (fnp0,
                                             m_interpTraceI0[set][cnt1]->
                                             GetPtr().get()+ tnp0 - 1,tnp0,
                                             &m_interpEndPtI0[set][cnt1][0],
                                             1);
                                    }
                                }
                            }
                            else
                            {
                                m_interpTrace[set][cnt1] = eInterpBothDirs;
                                m_interpTraceI0[set][cnt1] =
                                    LibUtilities::PointsManager()[fromPointsKey0]
                                    ->GetI(toPointsKey0);
                                m_interpTraceI1[set][cnt1] =
                                    LibUtilities::PointsManager()[fromPointsKey1]
                                    ->GetI(toPointsKey1);

                                // check to see if we can just
                                // interpolate endpoint
                                if ((fromPointsKey0.GetPointsType() ==
                                     LibUtilities::eGaussRadauMAlpha1Beta0) &&
                                    (toPointsKey0.GetPointsType() ==
                                     LibUtilities::eGaussLobattoLegendre))
                                {
                                    if (fromPointsKey0.GetNumPoints() + 1 ==
                                        toPointsKey0.GetNumPoints())
                                    {
                                        m_interpTrace[set][cnt1] =
                                            eInterpEndPtDir0InterpDir1;
                                        int fnp0 = fromPointsKey0.GetNumPoints();
                                        int tnp0 = toPointsKey0.GetNumPoints();
                                        m_interpEndPtI0[set][cnt1] =
                                            Array<OneD, NekDouble>(fnp0);
                                        Vmath::Vcopy
                                            (fnp0,
                                             m_interpTraceI0[set][cnt1]->
                                             GetPtr().get()+
                                             tnp0 - 1,tnp0,
                                             &m_interpEndPtI0[set][cnt1][0],1);
                                    }
                                }
                            }
                        }
                    }
                }

                if (set == 0)
                {
                    fwdSet = true;
                }
                else
                {
                    bwdSet = true;
                }
            }
        }
    }

    TraceLocToElmtLocCoeffMap(locExp,trace);
    FindElmtNeighbors(locExp,trace);
}

void LocTraceToTraceMap::CalcLocTracePhysToTraceIDMap(
    const ExpListSharedPtr &tracelist,
    const int               ndim)
{
    switch (ndim)
    {
        case 2:
            CalcLocTracePhysToTraceIDMap_2D(tracelist);
            break;
        case 3:
            CalcLocTracePhysToTraceIDMap_3D(tracelist);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, 
                "CalcLocTracePhysToTraceIDMap not coded");
    }
}

void LocTraceToTraceMap::CalcLocTracePhysToTraceIDMap_2D(
    const ExpListSharedPtr &tracelist)
{
    std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
    int ntotTrace            = (*traceExp).size();
    int ntPnts,noffset;

    m_LocTracephysToTraceIDMap      = Array<OneD, Array<OneD, int> > (2);
    m_LocTracephysToTraceIDMap[0]   = Array<OneD, int> (m_nFwdLocTracePts,-1);
    m_LocTracephysToTraceIDMap[1]   = Array<OneD, int> (
            m_nLocTracePts-m_nFwdLocTracePts,-1);

    Array<OneD, NekDouble> tracePnts(m_nTracePts,0.0);
    for(int nt=0; nt<ntotTrace;nt++)
    {
        ntPnts  =   tracelist->GetTotPoints(nt);
        noffset =   tracelist->GetPhys_Offset(nt);
        for(int i=0;i<ntPnts;i++)
        {
            tracePnts[noffset+i]    =   NekDouble(nt);
        }
    }
       
    Array<OneD, Array<OneD, NekDouble> > loctracePntsLR(2);
    loctracePntsLR[0]   =   Array<OneD, NekDouble> (m_nFwdLocTracePts,0.0);
    loctracePntsLR[1]   =   Array<OneD, NekDouble> (
        m_nLocTracePts-m_nFwdLocTracePts,0.0);

    for(int dir = 0; dir<2;dir++)
    {
        int cnt  = 0;
        int cnt1 = 0;

        Array<OneD, NekDouble> tmp(m_nTracePts,0.0);
        Vmath::Gathr((int)m_LocTraceToTraceMap[dir].size(),
                    tracePnts.get(),
                    m_LocTraceToTraceMap[dir].get(),
                    tmp.get());

        for (int i = 0; i < m_interpTrace[dir].size(); ++i)
        {
            if (m_interpNfaces[dir][i])
            {
                LibUtilities::PointsKey fromPointsKey0 =
                    std::get<0>(m_interpPoints[dir][i]);
                LibUtilities::PointsKey toPointsKey0 =
                    std::get<2>(m_interpPoints[dir][i]);

                int fnp    = fromPointsKey0.GetNumPoints();
                int tnp    = toPointsKey0.GetNumPoints();
                int nedges = m_interpNfaces[dir][i];

                for(int ne=0;ne<nedges;ne++)
                {
                    Vmath::Fill(fnp,tmp[cnt1],&loctracePntsLR[dir][cnt],1);
                    cnt += fnp;
                    cnt1 += tnp;
                }
            }
        }    
    }
    
    NekDouble error = 0.0;
    for(int nlr = 0; nlr<2;nlr++)
    {
        for(int i=0;i<loctracePntsLR[nlr].size();i++)
        {
            m_LocTracephysToTraceIDMap[nlr][i] =   
                std::round(loctracePntsLR[nlr][i]);
            error   +=   abs(loctracePntsLR[nlr][i] - NekDouble(
                m_LocTracephysToTraceIDMap[nlr][i]));
        }
    }
    error = error/NekDouble(m_nLocTracePts);
    ASSERTL0(error<NekConstants::kNekZeroTol,
        "m_LocTracephysToTraceIDMap may not be integer !!");
}

void LocTraceToTraceMap::CalcLocTracePhysToTraceIDMap_3D(
    const ExpListSharedPtr &tracelist)
{
    std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
    int ntotTrace            = (*traceExp).size();
    int ntPnts,noffset;

    m_LocTracephysToTraceIDMap      = Array<OneD, Array<OneD, int> > (2);
    m_LocTracephysToTraceIDMap[0]   = Array<OneD, int> (m_nFwdLocTracePts,-1);
    m_LocTracephysToTraceIDMap[1]   = Array<OneD, int> (
            m_nLocTracePts-m_nFwdLocTracePts,-1);

    Array<OneD, NekDouble> tracePnts(m_nTracePts,0.0);
    for(int nt=0; nt<ntotTrace;nt++)
    {
        ntPnts  =   tracelist->GetTotPoints(nt);
        noffset =   tracelist->GetPhys_Offset(nt);
        for(int i=0;i<ntPnts;i++)
        {
            tracePnts[noffset+i]    =   NekDouble(nt);
        }
    }
       
    Array<OneD, Array<OneD, NekDouble> > loctracePntsLR(2);
    loctracePntsLR[0]   =   Array<OneD, NekDouble> (m_nFwdLocTracePts,0.0);
    loctracePntsLR[1]   =   Array<OneD, NekDouble> (
        m_nLocTracePts-m_nFwdLocTracePts,0.0);

    for(int dir = 0; dir<2;dir++)
    {
        int cnt  = 0;
        int cnt1 = 0;

        // tmp space assuming forward map is of size of trace
        Array<OneD, NekDouble> tmp(m_nTracePts,0.0);
        Vmath::Gathr((int)m_LocTraceToTraceMap[dir].size(),
                    tracePnts.get(),
                    m_LocTraceToTraceMap[dir].get(),
                    tmp.get());

        for (int i = 0; i < m_interpTrace[dir].size(); ++i)
        {
            if (m_interpNfaces[dir][i])
            {
                LibUtilities::PointsKey fromPointsKey0 =
                    std::get<0>(m_interpPoints[dir][i]);
                LibUtilities::PointsKey fromPointsKey1 =
                    std::get<1>(m_interpPoints[dir][i]);
                LibUtilities::PointsKey toPointsKey0 =
                    std::get<2>(m_interpPoints[dir][i]);
                LibUtilities::PointsKey toPointsKey1 =
                    std::get<3>(m_interpPoints[dir][i]);

                int fnp0         = fromPointsKey0.GetNumPoints();
                int fnp1         = fromPointsKey1.GetNumPoints();
                int tnp0         = toPointsKey0.GetNumPoints();
                int tnp1         = toPointsKey1.GetNumPoints();

                int nfttl        = fnp0 * fnp1;
                
                for(int ne=0;ne<m_interpNfaces[dir][i];ne++)
                {
                    Vmath::Fill(nfttl,tmp[cnt1],&loctracePntsLR[dir][cnt],1);
                    cnt += nfttl;
                    cnt1 += tnp0 * tnp1;
                }
            }
        }    
    }
    
    NekDouble error = 0.0;
    for(int nlr = 0; nlr<2;nlr++)
    {
        for(int i=0;i<loctracePntsLR[nlr].size();i++)
        {
            m_LocTracephysToTraceIDMap[nlr][i] =   
                std::round(loctracePntsLR[nlr][i]);
            error   +=   abs(loctracePntsLR[nlr][i] - NekDouble(
                m_LocTracephysToTraceIDMap[nlr][i]));
        }
    }
    error = error/NekDouble(m_nLocTracePts);
    ASSERTL0(error<NekConstants::kNekZeroTol,
        "m_LocTracephysToTraceIDMap may not be integer !!");
}

/**
 * @brief Set up maps between coefficients on trace and in cells.
 *
 * @param locExp         Expansion list in elements
 * @param trace          Expansion list on traces.
 */
void LocTraceToTraceMap::TraceLocToElmtLocCoeffMap(
    const ExpList &locExp,
    const ExpListSharedPtr &trace)
{
    const std::shared_ptr<LocalRegions::ExpansionVector> exptrac =
        trace->GetExp();
    size_t ntrace = exptrac->size();

    Array<OneD, Array<OneD, int >> LRAdjExpid{2};
    Array<OneD, Array<OneD, bool>> LRAdjflag{2};

    TensorOfArray3D<int> elmtLRMap{2};
    TensorOfArray3D<int> elmtLRSign{2};

    for (int lr = 0; lr < 2; ++lr)
    {
        LRAdjExpid[lr]  =   Array<OneD, int > {ntrace, 0};
        LRAdjflag[lr]   =   Array<OneD, bool> {ntrace, false};
        elmtLRMap[lr]   =   Array<OneD, Array<OneD, int > > {ntrace};
        elmtLRSign[lr]  =   Array<OneD, Array<OneD, int > > {ntrace};
        for (int i = 0; i < ntrace; ++i)
        {
            size_t ncoeff  =   trace->GetNcoeffs(i);
            elmtLRMap[lr][i]      =   Array<OneD, int >{ncoeff, 0};
            elmtLRSign[lr][i]     =   Array<OneD, int >{ncoeff, 0};
        }
    }

    const Array<OneD, const pair<int, int> > field_coeffToElmt  =
            locExp.GetCoeffsToElmt();
    const Array<OneD, const pair<int, int> > trace_coeffToElmt  =
            trace->GetCoeffsToElmt();

    for (int lr = 0; lr < 2; ++lr)
    {
        int ntotcoeffs = m_nTraceCoeffs[lr];
        for (int  i = 0; i < ntotcoeffs; ++i)
        {
            int ncoeffField =   m_traceCoeffsToElmtMap[lr][i];
            int ncoeffTrace =   m_traceCoeffsToElmtTrace[lr][i];
            int sign        =   m_traceCoeffsToElmtSign[lr][i];

            int ntraceelmt   = trace_coeffToElmt[ncoeffTrace].first;
            int ntracelocN   = trace_coeffToElmt[ncoeffTrace].second;

            int nfieldelmt   = field_coeffToElmt[ncoeffField].first;
            int nfieldlocN   = field_coeffToElmt[ncoeffField].second;

            LRAdjflag[lr][ntraceelmt]    =   true;
            LRAdjExpid[lr][ntraceelmt]   =   nfieldelmt;

            elmtLRMap[lr][ntraceelmt][ntracelocN]  =   nfieldlocN;
            elmtLRSign[lr][ntraceelmt][ntracelocN]  =   sign;
        }
    }
    m_leftRightAdjacentExpId                = LRAdjExpid;
    m_leftRightAdjacentExpFlag              = LRAdjflag;
    m_traceCoeffToLeftRightExpCoeffMap      = elmtLRMap;
    m_traceCoeffToLeftRightExpCoeffSign     = elmtLRSign;
}

void LocTraceToTraceMap::FindElmtNeighbors(
    const ExpList &locExp,
    const ExpListSharedPtr &trace)
{
    const std::shared_ptr<LocalRegions::ExpansionVector> exptrac =
        trace->GetExp();
    int ntrace    = exptrac->size();

    const std::shared_ptr<LocalRegions::ExpansionVector> exp =
        locExp.GetExp();
    int nexp    = exp->size();

    Array<OneD, Array<OneD, int >> LRAdjExpid(2);
    Array<OneD, Array<OneD, bool>> LRAdjflag(2);
    LRAdjExpid  =   m_leftRightAdjacentExpId  ;
    LRAdjflag   =   m_leftRightAdjacentExpFlag;

    std::set< std::pair<int, int> > neighborSet;
    int ntmp0,ntmp1;
    for(int  nt = 0; nt < ntrace; nt++)
    {
        if(LRAdjflag[0][nt]&&LRAdjflag[1][nt])
        {
            ntmp0   =   LRAdjExpid[0][nt];
            ntmp1   =   LRAdjExpid[1][nt];
            
            ASSERTL0(ntmp0!=ntmp1, " ntmp0==ntmp1, trace inside a element?? ");

            std::set< std::pair<int, int> >::iterator it = neighborSet.begin();
            neighborSet.insert(it, std::make_pair(ntmp0,ntmp1)); 
            neighborSet.insert(it, std::make_pair(ntmp1,ntmp0));
        }
    }

    Array<OneD, int > ElemIndex(nexp,0);
    for (std::set< std::pair<int, int> >::iterator it=neighborSet.begin(); 
        it!=neighborSet.end(); ++it)
    {
        int ncurrent   =  it->first;
        ElemIndex[ncurrent]++;
    }

    Array<OneD, Array<OneD, int > > ElemNeighbsId(nexp);
    Array<OneD, Array<OneD, int > > tmpId(nexp);
    Array<OneD, int > ElemNeighbsNumb(nexp,-1);
    Vmath::Vcopy(nexp,ElemIndex,1,ElemNeighbsNumb,1);
    for(int  ne = 0; ne < nexp; ne++)
    {
        int neighb  =   ElemNeighbsNumb[ne];
        ElemNeighbsId[ne]   =   Array<OneD, int >(neighb,-1);
        tmpId[ne]           =   Array<OneD, int >(neighb,-1);
    }

    for(int  ne = 0; ne < nexp; ne++)
    {
        ElemIndex[ne]   =   0;
    }
    for (std::set< std::pair<int, int> >::iterator it=neighborSet.begin(); 
        it!=neighborSet.end(); ++it)
    {
        int ncurrent   =  it->first;
        int neighbor   =  it->second;
        ElemNeighbsId[ncurrent][ ElemIndex[ncurrent] ]    = neighbor;
        ElemIndex[ncurrent]++;
    }

    // pickout repeated indexes
    for(int  ne = 0; ne < nexp; ne++)
    {
        ElemIndex[ne]   =   0;
        for(int nb =0; nb<ElemNeighbsNumb[ne]; nb++)
        {
            int neighbId =  ElemNeighbsId[ne][nb];
            bool found = false;
            for(int nc =0; nc<ElemIndex[ne]; nc++)
            {
                if(ElemNeighbsId[ne][nb]==tmpId[ne][nc])
                {
                    found = true;
                }
            }
            if(!found)
            {
                tmpId[ne][ ElemIndex[ne] ] = neighbId;
                ElemIndex[ne]++;
            }
        }
    }
    ElemNeighbsNumb = ElemIndex;
    for(int  ne = 0; ne < nexp; ne++)
    {
        int neighb = ElemNeighbsNumb[ne];
        if(neighb>0)
        {
            ElemNeighbsId[ne]   =   Array<OneD, int >(neighb,-1);
            Vmath::Vcopy(neighb,tmpId[ne],1,ElemNeighbsId[ne],1);
        }
    }
    
    // check errors
    for(int  ne = 0; ne < nexp; ne++)
    {
        for(int nb =0; nb<ElemNeighbsNumb[ne]; nb++)
        {
            ASSERTL0( (ElemNeighbsId[ne][nb]>=0)&&(ElemNeighbsId[ne][nb]<=nexp),
                "Element id <0 or >number of total elements")
        }
    }

    m_ElemNeighbsNumb = ElemNeighbsNumb;
    m_ElemNeighbsId   = ElemNeighbsId;
}

/**
 * @brief Gather the local traces in physical space from field using
 * #m_fieldToLocTraceMap.
 *
 * @param field  Solution field in physical space
 * @param faces  Resulting local traces.
 */
void LocTraceToTraceMap::LocTracesFromField(
    const Array<OneD, const NekDouble> &field, Array<OneD, NekDouble> faces)
{
    // The static cast is necessary because m_fieldToLocTraceMap should be
    // Array<OneD, size_t> ... or at least the same type as
    // m_fieldToLocTraceMap.size() ...
    Vmath::Gathr(static_cast<int>(m_fieldToLocTraceMap.size()),
                 field,
                 m_fieldToLocTraceMap,
                 faces);
}

/**
 * @brief Reverse process of LocTracesFromField()
 * Add the local traces in physical space to field using
 * #m_fieldToLocTraceMap.
 *
 * @param field  Solution field in physical space
 * @param faces  local traces.
 */
void LocTraceToTraceMap::AddLocTracesToField(
    const Array<OneD, const NekDouble>  &faces,
    Array<OneD, NekDouble>              &field)
{
    size_t nfield  =   field.size();
    Array<OneD, NekDouble> tmp {nfield, 0.0};
    Vmath::Scatr(m_fieldToLocTraceMap.size(),
                 faces,
                 m_fieldToLocTraceMap,
                 tmp);
    Vmath::Vadd(nfield, tmp, 1, field, 1, field, 1);
}

/**
 * @brief Gather the forwards-oriented local traces in physical space from field
 * using #m_fieldToLocTraceMap.
 *
 * @param field  Solution field in physical space
 * @param faces  Resulting local forwards-oriented traces.
 */
void LocTraceToTraceMap::FwdLocTracesFromField(
    const Array<OneD, const NekDouble> &field, Array<OneD, NekDouble> faces)
{
    Vmath::Gathr(m_nFwdLocTracePts, field, m_fieldToLocTraceMap, faces);
}


void LocTraceToTraceMap::InterpLocTracesToTrace(
    const int dir,
    const Array<OneD, const NekDouble> &loctraces,
    Array<OneD, NekDouble> traces)
{
    switch(m_expdim)
    {
    case 1: // Essentially do copy
        Vmath::Scatr(m_LocTraceToTraceMap[dir].size(),
                     loctraces.get(),
                     m_LocTraceToTraceMap[dir].get(),
                     traces.get());
        break;
    case 2:
        InterpLocEdgesToTrace(dir,loctraces,traces);
        break;
    case 3:
        InterpLocFacesToTrace(dir,loctraces,traces);
        break;
    default:
        NEKERROR(ErrorUtil::efatal, "Not set up");
        break;
    }
}

/**
 * @brief Interpolate local trace edges to global trace edge point distributions
 * where required.
 *
 * @param dir       Selects forwards (0) or backwards (1) direction.
 * @param locfaces  Local trace edge storage.
 * @param faces     Global trace edge storage
 */
void LocTraceToTraceMap::InterpLocEdgesToTrace(
    const int dir,
    const Array<OneD, const NekDouble> &locedges,
    Array<OneD, NekDouble> edges)
{
    ASSERTL1(dir < 2,
             "option dir out of range, "
             " dir=0 is fwd, dir=1 is bwd");

    int cnt  = 0;
    int cnt1 = 0;

    // tmp space assuming forward map is of size of trace
    Array<OneD, NekDouble> tmp(m_nTracePts);

    for (int i = 0; i < m_interpTrace[dir].size(); ++i)
    {
        // Check if there are edges to interpolate
        if (m_interpNfaces[dir][i])
        {
            // Get to/from points
            LibUtilities::PointsKey fromPointsKey0 =
                std::get<0>(m_interpPoints[dir][i]);
            LibUtilities::PointsKey toPointsKey0 =
                std::get<2>(m_interpPoints[dir][i]);

            int fnp    = fromPointsKey0.GetNumPoints();
            int tnp    = toPointsKey0.GetNumPoints();
            int nedges = m_interpNfaces[dir][i];

            // Do interpolation here if required
            switch (m_interpTrace[dir][i])
            {
                case eNoInterp: // Just copy
                {
                    Vmath::Vcopy(nedges * fnp,
                                 locedges.get() + cnt,
                                 1,
                                 tmp.get() + cnt1,
                                 1);
                }
                break;
                case eInterpDir0:
                {
                    DNekMatSharedPtr I0 = m_interpTraceI0[dir][i];
                    Blas::Dgemm('N','N', tnp, nedges,
                                fnp,1.0, I0->GetPtr().get(),
                                tnp, locedges.get() + cnt,
                                fnp, 0.0, tmp.get() + cnt1,
                                tnp);
                }
                break;
                case eInterpEndPtDir0:
                {
                    Array<OneD, NekDouble> I0 = m_interpEndPtI0[dir][i];

                    for (int k = 0; k < nedges; ++k)
                    {
                        Vmath::Vcopy(fnp,
                                     &locedges[cnt + k * fnp],
                                     1,
                                     &tmp[cnt1 + k * tnp],
                                     1);

                        tmp[cnt1 + k * tnp + tnp - 1] = Blas::Ddot(
                            fnp, locedges.get() + cnt + k * fnp, 1, &I0[0], 1);
                    }
                }
                break;
                default:
                    NEKERROR(ErrorUtil::efatal,
                             "Invalid interpolation type for 2D elements");
                    break;
            }

            cnt += nedges * fnp;
            cnt1 += nedges * tnp;
        }
    }

    Vmath::Scatr(m_LocTraceToTraceMap[dir].size(),
                 tmp.get(),
                 m_LocTraceToTraceMap[dir].get(),
                 edges.get());
}

/**
 * @brief Right inner product with localedgetoTrace Interpolation Matrix.
 *
 * @param dir       Selects forwards (0) or backwards (1) direction.
 * @param locedges  Local trace edge storage.
 * @param edges     Global trace edge storage
 */
void LocTraceToTraceMap::RightIPTWLocEdgesToTraceInterpMat(
    const int                           dir,
    const Array<OneD, const NekDouble>  &edges,
    Array<OneD, NekDouble>              &locedges)
{
    ASSERTL1(dir < 2,
             "option dir out of range, "
             " dir=0 is fwd, dir=1 is bwd");

    int cnt  = 0;
    int cnt1 = 0;

    // tmp space assuming forward map is of size of trace
    Array<OneD, NekDouble> tmp{size_t(m_nTracePts)};
    // The static cast is necessary because m_LocTraceToTraceMap should be
    // Array<OneD, size_t> ... or at least the same type as
    // m_LocTraceToTraceMap.size() ...
    Vmath::Gathr(static_cast<int>(m_LocTraceToTraceMap[dir].size()),
                 edges,
                 m_LocTraceToTraceMap[dir],
                 tmp);

    for (int i = 0; i < m_interpTrace[dir].size(); ++i)
    {
        // Check if there are edges to interpolate
        if (m_interpNfaces[dir][i])
        {
            // Get to/from points
            LibUtilities::PointsKey fromPointsKey0 =
                std::get<0>(m_interpPoints[dir][i]);
            LibUtilities::PointsKey toPointsKey0 =
                std::get<2>(m_interpPoints[dir][i]);

            int fnp    = fromPointsKey0.GetNumPoints();
            int tnp    = toPointsKey0.GetNumPoints();
            int nedges = m_interpNfaces[dir][i];

            // Do interpolation here if required
            switch (m_interpTrace[dir][i])
            {
                case eNoInterp: // Just copy
                {
                    Vmath::Vcopy(nedges * fnp,
                                 tmp.get() + cnt1,
                                 1,
                                 locedges.get() + cnt,
                                 1);
                }
                break;
                case eInterpDir0:
                {
                    DNekMatSharedPtr I0 = m_interpTraceI0[dir][i];
                    Blas::Dgemm('T',
                                'N',
                                fnp,
                                nedges,
                                tnp,
                                1.0,
                                I0->GetPtr().get(),
                                tnp,
                                tmp.get() + cnt1,
                                tnp,
                                0.0,
                                locedges.get() + cnt,
                                fnp);
                }
                break;
                case eInterpEndPtDir0:
                {
                    Array<OneD, NekDouble> I0 = m_interpEndPtI0[dir][i];

                    for (int k = 0; k < nedges; ++k)
                    {
                        Vmath::Vcopy(fnp,
                                     &tmp[cnt1 + k * tnp],
                                     1,
                                     &locedges[cnt + k * fnp],
                                     1);

                        Vmath::Svtvp(fnp,tmp[cnt1 + k * tnp + tnp - 1],
                            &I0[0], 1,locedges.get() + cnt + k * fnp, 1,
                            locedges.get() + cnt + k * fnp, 1);
                    }
                }
                break;
                default:
                    NEKERROR(ErrorUtil::efatal,
                             "Invalid interpolation type for 2D elements");
                    break;
            }

            cnt += nedges * fnp;
            cnt1 += nedges * tnp;
        }
    }
}

/**
 * @brief Interpolate local faces to trace face point distributions where
 * required.
 *
 * @param dir       Selects forwards (0) or backwards (1) direction.
 * @param locfaces  Local trace face storage.
 * @param faces     Global trace face storage
 */
void LocTraceToTraceMap::InterpLocFacesToTrace(
    const int dir,
    const Array<OneD, const NekDouble> &locfaces,
    Array<OneD, NekDouble> faces)
{
    ASSERTL1(dir < 2,
             "option dir out of range, "
             " dir=0 is fwd, dir=1 is bwd");

    int cnt1 = 0;
    int cnt  = 0;

    // tmp space assuming forward map is of size of trace
    Array<OneD, NekDouble> tmp(m_nTracePts);

    for (int i = 0; i < m_interpTrace[dir].size(); ++i)
    {
        // Check if there are faces to interpolate
        if (m_interpNfaces[dir][i])
        {
            // Get to/from points
            LibUtilities::PointsKey fromPointsKey0 =
                std::get<0>(m_interpPoints[dir][i]);
            LibUtilities::PointsKey fromPointsKey1 =
                std::get<1>(m_interpPoints[dir][i]);
            LibUtilities::PointsKey toPointsKey0 =
                std::get<2>(m_interpPoints[dir][i]);
            LibUtilities::PointsKey toPointsKey1 =
                std::get<3>(m_interpPoints[dir][i]);

            int fnp0         = fromPointsKey0.GetNumPoints();
            int fnp1         = fromPointsKey1.GetNumPoints();
            int tnp0         = toPointsKey0.GetNumPoints();
            int tnp1         = toPointsKey1.GetNumPoints();
            int nfromfacepts = m_interpNfaces[dir][i] * fnp0 * fnp1;

            // Do interpolation here if required
            switch (m_interpTrace[dir][i])
            {
                case eNoInterp: // Just copy
                {
                    Vmath::Vcopy(nfromfacepts,
                                 locfaces.get() + cnt,
                                 1,
                                 tmp.get() + cnt1,
                                 1);
                }
                break;
                case eInterpDir0:
                {
                    DNekMatSharedPtr I0 = m_interpTraceI0[dir][i];
                    Blas::Dgemm('N',
                                'N',
                                tnp0,
                                tnp1,
                                fnp0,
                                1.0,
                                I0->GetPtr().get(),
                                tnp0,
                                locfaces.get() + cnt,
                                fnp0,
                                0.0,
                                tmp.get() + cnt1,
                                tnp0);
                }
                break;
                case eInterpEndPtDir0:
                {
                    int nfaces = m_interpNfaces[dir][i];
                    for (int k = 0; k < fnp0; ++k)
                    {
                        Vmath::Vcopy(nfaces * fnp1,
                                     locfaces.get() + cnt + k,
                                     fnp0,
                                     tmp.get() + cnt1 + k,
                                     tnp0);
                    }
                    Array<OneD, NekDouble> I0 = m_interpEndPtI0[dir][i];
                    Blas::Dgemv('T',
                                fnp0,
                                tnp1 * m_interpNfaces[dir][i],
                                1.0,
                                tmp.get() + cnt1,
                                tnp0,
                                I0.get(),
                                1,
                                0.0,
                                tmp.get() + cnt1 + tnp0 - 1,
                                tnp0);
                }
                break;
                case eInterpDir1:
                {
                    DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];
                    for (int j = 0; j < m_interpNfaces[dir][i]; ++j)
                    {
                        Blas::Dgemm('N',
                                    'T',
                                    tnp0,
                                    tnp1,
                                    fnp1,
                                    1.0,
                                    locfaces.get() + cnt + j * fnp0 * fnp1,
                                    tnp0,
                                    I1->GetPtr().get(),
                                    tnp1,
                                    0.0,
                                    tmp.get() + cnt1 + j * tnp0 * tnp1,
                                    tnp0);
                    }
                }
                break;
                case eInterpEndPtDir1:
                {
                    Array<OneD, NekDouble> I1 = m_interpEndPtI1[dir][i];
                    for (int j = 0; j < m_interpNfaces[dir][i]; ++j)
                    {
                        // copy all points
                        Vmath::Vcopy(fnp0 * fnp1,
                                     locfaces.get() + cnt + j * fnp0 * fnp1,
                                     1,
                                     tmp.get() + cnt1 + j * tnp0 * tnp1,
                                     1);

                        // interpolate end points
                        for (int k = 0; k < tnp0; ++k)
                        {
                            tmp[cnt1 + k + (j + 1) * tnp0 * tnp1 - tnp0] =
                                Blas::Ddot(fnp1,
                                           locfaces.get() + cnt +
                                               j * fnp0 * fnp1 + k,
                                           fnp0,
                                           &I1[0],
                                           1);
                        }
                    }
                }
                break;
                case eInterpBothDirs:
                {
                    DNekMatSharedPtr I0 = m_interpTraceI0[dir][i];
                    DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];
                    Array<OneD, NekDouble> wsp(m_interpNfaces[dir][i] * fnp0 *
                                               tnp1 * fnp0);

                    for (int j = 0; j < m_interpNfaces[dir][i]; ++j)
                    {
                        Blas::Dgemm('N',
                                    'T',
                                    fnp0,
                                    tnp1,
                                    fnp1,
                                    1.0,
                                    locfaces.get() + cnt + j * fnp0 * fnp1,
                                    fnp0,
                                    I1->GetPtr().get(),
                                    tnp1,
                                    0.0,
                                    wsp.get() + j * fnp0 * tnp1,
                                    fnp0);
                    }
                    Blas::Dgemm('N',
                                'N',
                                tnp0,
                                tnp1 * m_interpNfaces[dir][i],
                                fnp0,
                                1.0,
                                I0->GetPtr().get(),
                                tnp0,
                                wsp.get(),
                                fnp0,
                                0.0,
                                tmp.get() + cnt1,
                                tnp0);
                }
                break;
                case eInterpEndPtDir0InterpDir1:
                {
                    DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];

                    for (int j = 0; j < m_interpNfaces[dir][i]; ++j)
                    {
                        Blas::Dgemm('N',
                                    'T',
                                    fnp0,
                                    tnp1,
                                    fnp1,
                                    1.0,
                                    locfaces.get() + cnt + j * fnp0 * fnp1,
                                    fnp0,
                                    I1->GetPtr().get(),
                                    tnp1,
                                    0.0,
                                    tmp.get() + cnt1 + j * tnp0 * tnp1,
                                    tnp0);
                    }

                    Array<OneD, NekDouble> I0 = m_interpEndPtI0[dir][i];
                    Blas::Dgemv('T',
                                fnp0,
                                tnp1 * m_interpNfaces[dir][i],
                                1.0,
                                tmp.get() + cnt1,
                                tnp0,
                                I0.get(),
                                1,
                                0.0,
                                tmp.get() + cnt1 + tnp0 - 1,
                                tnp0);
                }
                break;
            }
            cnt += nfromfacepts;
            cnt1 += m_interpNfaces[dir][i] * tnp0 * tnp1;
        }
    }

    Vmath::Scatr(m_LocTraceToTraceMap[dir].size(),
                 tmp.get(),
                 m_LocTraceToTraceMap[dir].get(),
                 faces.get());
}

/**
 * @brief Right inner product with localedgetoTrace Interpolation Matrix.
 *
 * @param dir           Selects forwards (0) or backwards (1) direction.
 * @param traces        trace .
 * @param loctraces     Local trace
 */
void LocTraceToTraceMap::RightIPTWLocFacesToTraceInterpMat(
    const int                           dir,
    const Array<OneD, const NekDouble>  &traces,
    Array<OneD, NekDouble>              &loctraces)
{
    ASSERTL1(dir < 2,
             "option dir out of range, "
             " dir=0 is fwd, dir=1 is bwd");

    int cnt  = 0;
    int cnt1 = 0;

    // tmp space assuming forward map is of size of trace
    Array<OneD, NekDouble> tmp{size_t(m_nTracePts)};
    // The static cast is necessary because m_LocTraceToTraceMap should be
    // Array<OneD, size_t> ... or at least the same type as
    // m_LocTraceToTraceMap.size() ...
    Vmath::Gathr(static_cast<int>(m_LocTraceToTraceMap[dir].size()),
                 traces,
                 m_LocTraceToTraceMap[dir],
                 tmp);

    for (int i = 0; i < m_interpTrace[dir].size(); ++i)
    {
        // Check if there are elementboundaries to interpolate
        if (m_interpNfaces[dir][i])
        {
            // Get to/from points
            LibUtilities::PointsKey fromPointsKey0 =
                std::get<0>(m_interpPoints[dir][i]);
            LibUtilities::PointsKey fromPointsKey1 =
                std::get<1>(m_interpPoints[dir][i]);
            LibUtilities::PointsKey toPointsKey0 =
                std::get<2>(m_interpPoints[dir][i]);
            LibUtilities::PointsKey toPointsKey1 =
                std::get<3>(m_interpPoints[dir][i]);
            // Here the f(from) and t(to) are chosen to be consistent with 
            // InterpLocFacesToTrace
            int fnp0         = fromPointsKey0.GetNumPoints();
            int fnp1         = fromPointsKey1.GetNumPoints();
            int tnp0         = toPointsKey0.GetNumPoints();
            int tnp1         = toPointsKey1.GetNumPoints();
            int nfromfacepts = m_interpNfaces[dir][i] * fnp0 * fnp1;

            // Do interpolation here if required
            switch (m_interpTrace[dir][i])
            {
                case eNoInterp: // Just copy
                {
                    Vmath::Vcopy(nfromfacepts,
                                 tmp.get() + cnt1,
                                 1,
                                 loctraces.get() + cnt,
                                 1);
                }
                break;
                case eInterpDir0:
                {
                    DNekMatSharedPtr I0 = m_interpTraceI0[dir][i];
                    Blas::Dgemm('T',
                                'N',
                                fnp0,
                                tnp1,
                                tnp0,
                                1.0,
                                I0->GetPtr().get(),
                                tnp0,
                                tmp.get() + cnt1,
                                tnp0,
                                0.0,
                                loctraces.get() + cnt,
                                fnp0);
                }
                break;
                case eInterpEndPtDir0:
                {
                    int nfaces = m_interpNfaces[dir][i];
                    for (int k = 0; k < fnp0; ++k)
                    {
                        Vmath::Vcopy(nfaces * fnp1,
                                     tmp.get() + cnt1 + k,
                                     tnp0,
                                     loctraces.get() + cnt + k,
                                     fnp0);
                    }
                    Array<OneD, NekDouble> I0 = m_interpEndPtI0[dir][i];
                    for(int k = 0; k< tnp1 * m_interpNfaces[dir][i]; k++)
                    {
                        Vmath::Svtvp(fnp0,tmp[cnt1 + tnp0-1+k*tnp0],
                                    &I0[0],1,&loctraces[cnt],1,
                                    &loctraces[cnt],1);
                    }
                }
                break;
                case eInterpDir1:
                {
                    DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];

                    for (int j = 0; j < m_interpNfaces[dir][i]; ++j)
                    {
                        Blas::Dgemm('N',
                                    'N',
                                    tnp0,
                                    fnp1,
                                    tnp1,
                                    1.0,
                                    tmp.get() + cnt1 + j * tnp0 * tnp1,
                                    tnp0,
                                    I1->GetPtr().get(),
                                    tnp1,
                                    0.0,
                                    loctraces.get() + cnt + j * fnp0 * fnp1,
                                    tnp0);
                    }
                }
                break;
                case eInterpEndPtDir1:
                {
                    Array<OneD, NekDouble> I1 = m_interpEndPtI1[dir][i];
                    for (int j = 0; j < m_interpNfaces[dir][i]; ++j)
                    {
                        Vmath::Vcopy(fnp0 * fnp1,
                                     tmp.get() + cnt1 + j * tnp0 * tnp1,
                                     1,
                                     loctraces.get() + cnt + j * fnp0 * fnp1,
                                     1);

                        for(int k = 0; k< fnp1; k++)
                        {
                            Vmath::Svtvp(fnp0,I1[k],
                                &tmp[cnt1 + (j + 1) * tnp0 * tnp1 - tnp0],1,
                                &loctraces[cnt+k*fnp0],1,
                                &loctraces[cnt+k*fnp0],1);
                        }
                    }

                }
                break;
                case eInterpBothDirs:
                {
                    DNekMatSharedPtr I0 = m_interpTraceI0[dir][i];
                    DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];

                    Array<OneD, NekDouble>
                        wsp{size_t(m_interpNfaces[dir][i] * fnp0 * tnp1)};

                    Blas::Dgemm('T',
                                'N',
                                fnp0,
                                tnp1 * m_interpNfaces[dir][i],
                                tnp0,
                                1.0,
                                I0->GetPtr().get(),
                                tnp0,
                                tmp.get() + cnt1,
                                tnp0,
                                0.0,
                                wsp.get(),
                                fnp0);
                    for (int j = 0; j < m_interpNfaces[dir][i]; ++j)
                    {
                        Blas::Dgemm('N',
                                    'N',
                                    fnp0,
                                    fnp1,
                                    tnp1,
                                    1.0,
                                    wsp.get() + j * fnp0 * tnp1,
                                    fnp0,
                                    I1->GetPtr().get(),
                                    tnp1,
                                    0.0,
                                    loctraces.get() + cnt + j * fnp0 * fnp1,
                                    fnp0);
                    }
                }
                break;
                case eInterpEndPtDir0InterpDir1:
                {
                    DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];

                    for (int j = 0; j < m_interpNfaces[dir][i]; ++j)
                    {
                        Blas::Dgemm('N',
                                    'N',
                                    fnp0,
                                    fnp1,
                                    tnp1,
                                    1.0,
                                    tmp.get() + cnt1 + j * tnp0 * tnp1,
                                    tnp0,
                                    I1->GetPtr().get(),
                                    tnp1,
                                    0.0,
                                    loctraces.get() + cnt + j * fnp0 * fnp1,
                                    fnp0);
                    }

                    Array<OneD, NekDouble> I0 = m_interpEndPtI0[dir][i];
                    for(int k = 0; k< tnp1 * m_interpNfaces[dir][i]; k++)
                    {
                        Vmath::Svtvp(fnp0,tmp[cnt1 + tnp0-1+k*tnp0],
                                    &I0[0],1,&loctraces[cnt],1,
                                    &loctraces[cnt],1);
                    }
                }
                break;
            }
            cnt += nfromfacepts;
            cnt1 += m_interpNfaces[dir][i] * tnp0 * tnp1;
        }
    }
}

/**
 * @brief Add contributions from trace coefficients to the elemental field
 * storage.
 *
 * @param trace  Array of global trace coefficients.
 * @param field  Array containing field coefficients storage.
 */
void LocTraceToTraceMap::AddTraceCoeffsToFieldCoeffs(
    const Array<OneD, const NekDouble> &trace, Array<OneD, NekDouble> &field)
{
    int nvals = m_nTraceCoeffs[0] + m_nTraceCoeffs[1];
    for (int i = 0; i < nvals; ++i)
    {
        field[m_traceCoeffsToElmtMap[0][i]] +=
            m_traceCoeffsToElmtSign[0][i] *
            trace[m_traceCoeffsToElmtTrace[0][i]];
    }
}

/**
 * @brief Add contributions from backwards or forwards oriented trace
 * coefficients to the elemental field storage.
 *
 * @param dir    Selects forwards (0) or backwards (1) direction
 * @param trace  Array of global trace coefficients.
 * @param field  Array containing field coefficients storage.
 */
void LocTraceToTraceMap::AddTraceCoeffsToFieldCoeffs(
    const int dir,
    const Array<OneD, const NekDouble> &trace,
    Array<OneD, NekDouble> &field)
{
    int nvals = m_nTraceCoeffs[dir];
    for (int i = 0; i < nvals; ++i)
    {
        field[m_traceCoeffsToElmtMap[dir][i]] +=
            m_traceCoeffsToElmtSign[dir][i] *
            trace[m_traceCoeffsToElmtTrace[dir][i]];
    }
}

}
}
