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
    int m_expdim = locExpVector[0]->GetShapeDimension();

    // Switch between 1D, 2D and 3D
    switch (m_expdim)
    {
        case 1:
            break;
        case 2:
            Setup2D(locExp, trace, elmtToTrace, LeftAdjacents);
            break;
        case 3:
            Setup3D(locExp, trace, elmtToTrace, LeftAdjacents);
            break;
        default:
            ASSERTL0(false, "Number of dimensions greater than 3")
            break;
    }
}

LocTraceToTraceMap::~LocTraceToTraceMap()
{
}

/**
 * @brief Set up member variables for a two-dimensional problem.
 *
 * @param locExp         Expansion list of 2D elements
 * @param trace          Expansion list of the one-dimensional trace.
 * @param elmtToTrace    Mapping from elemental edges to trace.
 * @param leftAdjacents  Vector of bools denoting forwards-oriented traces.
 */
void LocTraceToTraceMap::Setup2D(
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

    m_traceCoeffsToElmtMap   = Array<OneD, Array<OneD, int> >(2);
    m_traceCoeffsToElmtTrace = Array<OneD, Array<OneD, int> >(2);
    m_traceCoeffsToElmtSign  = Array<OneD, Array<OneD, int> >(2);

    LocalRegions::Expansion2DSharedPtr exp2d;
    const std::shared_ptr<LocalRegions::ExpansionVector> exp =
        locExp.GetExp();

    int cnt, n, e, phys_offset;

    int nexp    = exp->size();
    m_nTracePts = trace->GetTotPoints();

    // Count number of edges and points required for maps
    int nFwdPts       = 0;
    int nBwdPts       = 0;
    int nFwdCoeffs    = 0;
    int nBwdCoeffs    = 0;
    m_nFwdLocTracePts = 0;
    m_nLocTracePts    = 0;

    for (cnt = n = 0; n < nexp; ++n)
    {
        exp2d = (*exp)[n]->as<LocalRegions::Expansion2D>();

        for (int i = 0; i < exp2d->GetNedges(); ++i, ++cnt)
        {
            int nLocPts = exp2d->GetEdgeNumPoints(i);
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
    int nedge;
    int fwdcnt = 0;
    int bwdcnt = 0;

    // Generate a map of similar traces with the same
    // interpolation requirements
    for (cnt = n = 0; n < nexp; ++n)
    {
        exp2d = (*exp)[n]->as<LocalRegions::Expansion2D>();
        nedge = exp2d->GetNedges();
        TraceOrder[n].resize(nedge);

        int coeffoffset = locExp.GetCoeff_Offset(n);
        for (e = 0; e < nedge; ++e, ++cnt)
        {
            StdRegions::StdExpansionSharedPtr edge = elmtToTrace[n][e];
            StdRegions::Orientation orient         = exp2d->GetEorient(e);

            LibUtilities::PointsKey fromPointsKey0, fromPointsKey1;
            LibUtilities::PointsKey toPointsKey0, toPointsKey1;

            // For Spencer's data structure
            int basisDir = 0;
            if (exp2d->DetShapeType() == LibUtilities::eTriangle)
            {
                basisDir = e == 0 ? 0 : 1;
            }
            else
            {
                basisDir = e % 2;
            }

            fromPointsKey0 = exp2d->GetBasis(basisDir)->GetPointsKey();
            fromPointsKey1 =
                LibUtilities::PointsKey(0, LibUtilities::eNoPointsType);
            toPointsKey0 = edge->GetBasis(0)->GetPointsKey();
            toPointsKey1 =
                LibUtilities::PointsKey(0, LibUtilities::eNoPointsType);

            // Spencer's data structure
            TraceInterpPoints fpoint(
                fromPointsKey0, fromPointsKey1, toPointsKey0, toPointsKey1);

            pair<int, int> epf(n, e);
            TraceInterpMap[fpoint].push_back(epf);
            TraceOrder[n][e] = cnt;

            // Setup for coefficient mapping from trace normal flux
            // to elements
            Array<OneD, unsigned int> map;
            Array<OneD, int> sign;
            // Not sure about the call GetBasisNumModes passing (0)!
            exp2d->GetEdgeToElementMap(
                e, orient, map, sign, edge->GetBasisNumModes(0));

            int order_f = edge->GetNcoeffs();
            int foffset = trace->GetCoeff_Offset(edge->GetElmtId());

            double fac = 1.0;

            LocalRegions::Expansion1DSharedPtr locExp1d =
                elmtToTrace[n][e]->as<LocalRegions::Expansion1D>();

            if (exp2d->GetEdgeExp(e)->GetRightAdjacentElementExp())
            {
                if (locExp1d->GetRightAdjacentElementExp()
                        ->GetGeom()
                        ->GetGlobalID() == exp2d->GetGeom()->GetGlobalID())
                {
                    fac = -1.0;
                }
            }

            if (LeftAdjacents[cnt])
            {
                for (int i = 0; i < order_f; ++i)
                {
                    m_traceCoeffsToElmtMap[0][fwdcnt]    = coeffoffset + map[i];
                    m_traceCoeffsToElmtTrace[0][fwdcnt]  = foffset + i;
                    m_traceCoeffsToElmtSign[0][fwdcnt++] = fac * sign[i];
                }
            }
            else
            {
                for (int i = 0; i < order_f; ++i)
                {
                    m_traceCoeffsToElmtMap[1][bwdcnt]    = coeffoffset + map[i];
                    m_traceCoeffsToElmtTrace[1][bwdcnt]  = foffset + i;
                    m_traceCoeffsToElmtSign[1][bwdcnt++] = fac * sign[i];
                }
            }
        }
    }

    int nInterpType = TraceInterpMap.size();

    for (int i = 0; i < 2; ++i)
    {
        m_interpTrace[i]   = Array<OneD, InterpLocTraceToTrace>(nInterpType);
        m_interpTraceI0[i] = Array<OneD, DNekMatSharedPtr>(nInterpType);
        m_interpEndPtI0[i] = Array<OneD, Array<OneD, NekDouble> >(nInterpType);
        m_interpPoints[i]  = Array<OneD, TraceInterpPoints>(nInterpType);
        m_interpNfaces[i]  = Array<OneD, int>(nInterpType, 0);
    }

    int nedgepts, nedgepts1;
    int cnt1    = 0;
    int cnt2    = 0;
    int cntFwd  = 0;
    int cntBwd  = 0;
    int cntFwd1 = 0;
    int cntBwd1 = 0;
    int set;
    Array<OneD, int> edgeids;
    Array<OneD, int> locTraceToTraceMap;
    cnt = 0;

    for (auto it = TraceInterpMap.begin(); it != TraceInterpMap.end();
         ++it, ++cnt1)
    {
        LibUtilities::PointsKey fromPointsKey0 = std::get<0>(it->first);
        LibUtilities::PointsKey toPointsKey0   = std::get<2>(it->first);

        bool fwdSet = false;
        bool bwdSet = false;

        for (int f = 0; f < it->second.size(); ++f, ++cnt2)
        {
            n = it->second[f].first;
            e = it->second[f].second;

            StdRegions::StdExpansionSharedPtr edge = elmtToTrace[n][e];

            exp2d       = (*exp)[n]->as<LocalRegions::Expansion2D>();
            phys_offset = locExp.GetPhys_Offset(n);

            // Mapping of new edge order to one that loops over elmts
            // then set up mapping of faces in standard cartesian order
            exp2d->GetEdgePhysMap(e, edgeids);

            nedgepts  = exp2d->GetEdgeNumPoints(e);
            nedgepts1 = edge->GetTotPoints();

            StdRegions::Orientation orient = exp2d->GetEorient(e);

            // Account for eBackwards orientation
            exp2d->ReOrientEdgePhysMap(elmtToTrace[n][e]->GetNverts(),
                                       orient,
                                       toPointsKey0.GetNumPoints(),
                                       locTraceToTraceMap);

            int offset = trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

            if (LeftAdjacents[TraceOrder[n][e]])
            {
                for (int i = 0; i < nedgepts; ++i)
                {
                    m_fieldToLocTraceMap[cntFwd + i] = phys_offset + edgeids[i];
                }

                for (int i = 0; i < nedgepts1; ++i)
                {
                    m_LocTraceToTraceMap[0][cntFwd1 + i] =
                        offset + locTraceToTraceMap[i];
                }

                cntFwd += nedgepts;
                cntFwd1 += nedgepts1;
                set = 0;
            }
            else
            {
                for (int i = 0; i < nedgepts; ++i)
                {
                    m_fieldToLocTraceMap[m_nFwdLocTracePts + cntBwd + i] =
                        phys_offset + edgeids[i];
                }

                for (int i = 0; i < nedgepts1; ++i)
                {
                    m_LocTraceToTraceMap[1][cntBwd1 + i] =
                        offset + locTraceToTraceMap[i];
                }

                cntBwd += nedgepts;
                cntBwd1 += nedgepts1;
                set = 1;
            }

            m_interpNfaces[set][cnt1] += 1;

            if ((fwdSet == false && set == 0) || (bwdSet == false && set == 1))
            {
                m_interpPoints[set][cnt1] = it->first;

                if (fromPointsKey0 == toPointsKey0)
                {
                    m_interpTrace[set][cnt1] = eNoInterp;
                }
                else
                {
                    m_interpTrace[set][cnt1] = eInterpDir0;
                    m_interpTraceI0[set][cnt1] =
                        LibUtilities::PointsManager()[fromPointsKey0]->GetI(
                            toPointsKey0);

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

                            Vmath::Vcopy(
                                fnp0,
                                m_interpTraceI0[set][cnt1]->GetPtr().get() +
                                    tnp0 - 1,
                                tnp0,
                                &m_interpEndPtI0[set][cnt1][0],
                                1);
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
}

/**
 * @brief Set up member variables for a three-dimensional problem.
 *
 * @param locExp         Expansion list of 3D elements
 * @param trace          Expansion list of the two-dimensional trace.
 * @param elmtToTrace    Mapping from elemental faces to trace.
 * @param leftAdjacents  Vector of bools denoting forwards-oriented traces.
 */
void LocTraceToTraceMap::Setup3D(
    const ExpList &locExp,
    const ExpListSharedPtr &trace,
    const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
        &elmtToTrace,
    const vector<bool> &LeftAdjacents)
{
    m_LocTraceToTraceMap = Array<OneD, Array<OneD, int> >(2);

    m_interpTrace   = Array<OneD, Array<OneD, InterpLocTraceToTrace> >(2);
    m_interpTraceI0 = Array<OneD, Array<OneD, DNekMatSharedPtr> >(2);
    m_interpTraceI1 = Array<OneD, Array<OneD, DNekMatSharedPtr> >(2);
    m_interpEndPtI0 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(2);
    m_interpEndPtI1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(2);
    m_interpPoints  = Array<OneD, Array<OneD, TraceInterpPoints> >(2);
    m_interpNfaces  = Array<OneD, Array<OneD, int> >(2);

    m_traceCoeffsToElmtMap   = Array<OneD, Array<OneD, int> >(2);
    m_traceCoeffsToElmtTrace = Array<OneD, Array<OneD, int> >(2);
    m_traceCoeffsToElmtSign  = Array<OneD, Array<OneD, int> >(2);

    LocalRegions::Expansion3DSharedPtr exp3d;
    const std::shared_ptr<LocalRegions::ExpansionVector> exp =
        locExp.GetExp();

    int cnt, n, e, phys_offset;

    int nexp    = exp->size();
    m_nTracePts = trace->GetTotPoints();

    // Count number of faces and points required for maps
    int nFwdPts       = 0;
    int nBwdPts       = 0;
    int nFwdCoeffs    = 0;
    int nBwdCoeffs    = 0;
    m_nFwdLocTracePts = 0;
    m_nLocTracePts    = 0;

    for (cnt = n = 0; n < nexp; ++n)
    {
        exp3d = (*exp)[n]->as<LocalRegions::Expansion3D>();

        for (int i = 0; i < exp3d->GetNfaces(); ++i, ++cnt)
        {
            int nLocPts = exp3d->GetFaceNumPoints(i);
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
    int nface;
    int fwdcnt = 0;
    int bwdcnt = 0;

    // Generate a map of similar traces with the same
    // interpolation or projection requirements
    for (cnt = n = 0; n < nexp; ++n)
    {
        exp3d = (*exp)[n]->as<LocalRegions::Expansion3D>();
        nface = exp3d->GetNfaces();
        TraceOrder[n].resize(nface);

        int coeffoffset = locExp.GetCoeff_Offset(n);
        for (e = 0; e < nface; ++e, ++cnt)
        {
            StdRegions::StdExpansionSharedPtr face = elmtToTrace[n][e];
            StdRegions::Orientation orient         = exp3d->GetForient(e);

            LibUtilities::PointsKey fromPointsKey0, fromPointsKey1;
            LibUtilities::PointsKey toPointsKey0, toPointsKey1;

            // 3D specific
            int dir0 = exp3d->GetGeom3D()->GetDir(e, 0);
            int dir1 = exp3d->GetGeom3D()->GetDir(e, 1);

            fromPointsKey0 = exp3d->GetBasis(dir0)->GetPointsKey();
            fromPointsKey1 = exp3d->GetBasis(dir1)->GetPointsKey();

            if (orient < StdRegions::eDir1FwdDir2_Dir2FwdDir1)
            {
                toPointsKey0 = face->GetBasis(0)->GetPointsKey();
                toPointsKey1 = face->GetBasis(1)->GetPointsKey();
            }
            else // transpose points key evaluation
            {
                toPointsKey0 = face->GetBasis(1)->GetPointsKey();
                toPointsKey1 = face->GetBasis(0)->GetPointsKey();
            }

            TraceInterpPoints fpoint(
                fromPointsKey0, fromPointsKey1, toPointsKey0, toPointsKey1);

            pair<int, int> epf(n, e);
            TraceInterpMap[fpoint].push_back(epf);
            TraceOrder[n][e] = cnt;

            // Setup for coefficient mapping from trace normal
            // flux to elements
            Array<OneD, unsigned int> map;
            Array<OneD, int> sign;
            exp3d->GetFaceToElementMap(e,
                                       orient,
                                       map,
                                       sign,
                                       face->GetBasisNumModes(0),
                                       face->GetBasisNumModes(1));

            int order_f = face->GetNcoeffs();
            int foffset = trace->GetCoeff_Offset(face->GetElmtId());

            int fac = 1.0;

            if (exp3d->GetFaceExp(e)->GetRightAdjacentElementExp())
            {
                if (exp3d->GetFaceExp(e)
                        ->GetRightAdjacentElementExp()
                        ->GetGeom3D()
                        ->GetGlobalID() == exp3d->GetGeom3D()->GetGlobalID())
                {
                    fac = -1.0;
                }
            }

            if (LeftAdjacents[cnt])
            {
                for (int i = 0; i < order_f; ++i)
                {
                    m_traceCoeffsToElmtMap[0][fwdcnt]    = coeffoffset + map[i];
                    m_traceCoeffsToElmtTrace[0][fwdcnt]  = foffset + i;
                    m_traceCoeffsToElmtSign[0][fwdcnt++] = fac * sign[i];
                }
            }
            else
            {
                for (int i = 0; i < order_f; ++i)
                {
                    m_traceCoeffsToElmtMap[1][bwdcnt]    = coeffoffset + map[i];
                    m_traceCoeffsToElmtTrace[1][bwdcnt]  = foffset + i;
                    m_traceCoeffsToElmtSign[1][bwdcnt++] = fac * sign[i];
                }
            }
        }
    }

    int nInterpType = TraceInterpMap.size();
    for (int i = 0; i < 2; ++i)
    {
        m_interpTrace[i]   = Array<OneD, InterpLocTraceToTrace>(nInterpType);
        m_interpTraceI0[i] = Array<OneD, DNekMatSharedPtr>(nInterpType);
        m_interpTraceI1[i] = Array<OneD, DNekMatSharedPtr>(nInterpType);
        m_interpEndPtI0[i] = Array<OneD, Array<OneD, NekDouble> >(nInterpType);
        m_interpEndPtI1[i] = Array<OneD, Array<OneD, NekDouble> >(nInterpType);
        m_interpPoints[i]  = Array<OneD, TraceInterpPoints>(nInterpType);
        m_interpNfaces[i]  = Array<OneD, int>(nInterpType, 0);
    }

    int nfacepts, nfacepts1;
    int cnt1    = 0;
    int cnt2    = 0;
    int cntFwd  = 0;
    int cntBwd  = 0;
    int cntFwd1 = 0;
    int cntBwd1 = 0;
    int set;
    Array<OneD, int> faceids;
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

            StdRegions::StdExpansionSharedPtr face = elmtToTrace[n][e];

            exp3d       = (*exp)[n]->as<LocalRegions::Expansion3D>();
            phys_offset = locExp.GetPhys_Offset(n);

            // mapping of new face order to one that loops over elmts
            // then faces set up mapping of faces in standard cartesian
            // order
            exp3d->GetFacePhysMap(e, faceids);
            nfacepts  = exp3d->GetFaceNumPoints(e);
            nfacepts1 = face->GetTotPoints();

            StdRegions::Orientation orient = exp3d->GetForient(e);

            exp3d->ReOrientFacePhysMap(elmtToTrace[n][e]->GetNverts(),
                                       orient,
                                       toPointsKey0.GetNumPoints(),
                                       toPointsKey1.GetNumPoints(),
                                       locTraceToTraceMap);

            int offset = trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

            if (LeftAdjacents[TraceOrder[n][e]])
            {
                for (int i = 0; i < nfacepts; ++i)
                {
                    m_fieldToLocTraceMap[cntFwd + i] = phys_offset + faceids[i];
                }

                for (int i = 0; i < nfacepts1; ++i)
                {
                    m_LocTraceToTraceMap[0][cntFwd1 + i] =
                        offset + locTraceToTraceMap[i];
                }

                cntFwd += nfacepts;
                cntFwd1 += nfacepts1;
                set = 0;
            }
            else
            {
                for (int i = 0; i < nfacepts; ++i)
                {
                    m_fieldToLocTraceMap[m_nFwdLocTracePts + cntBwd + i] =
                        phys_offset + faceids[i];
                }

                for (int i = 0; i < nfacepts1; ++i)
                {
                    m_LocTraceToTraceMap[1][cntBwd1 + i] =
                        offset + locTraceToTraceMap[i];
                }

                cntBwd += nfacepts;
                cntBwd1 += nfacepts1;
                set = 1;
            }

            m_interpNfaces[set][cnt1] += 1;

            if (((fwdSet == false) && (set == 0)) ||
                ((bwdSet == false) && (set == 1)))
            {
                m_interpPoints[set][cnt1] = it->first;

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
                            LibUtilities::PointsManager()[fromPointsKey1]->GetI(
                                toPointsKey1);

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
                                int fnp1                 = fromPointsKey1.GetNumPoints();
                                int tnp1 = toPointsKey1.GetNumPoints();
                                m_interpEndPtI1[set][cnt1] =
                                    Array<OneD, NekDouble>(fnp1);
                                Vmath::Vcopy(
                                    fnp1,
                                    m_interpTraceI1[set][cnt1]->GetPtr().get() +
                                        tnp1 - 1,
                                    tnp1,
                                    &m_interpEndPtI1[set][cnt1][0],
                                    1);
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
                            LibUtilities::PointsManager()[fromPointsKey0]->GetI(
                                toPointsKey0);

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
                                m_interpTrace[set][cnt1] = eInterpEndPtDir0;
                                int fnp0                 = fromPointsKey0.GetNumPoints();
                                int tnp0 = toPointsKey0.GetNumPoints();
                                m_interpEndPtI0[set][cnt1] =
                                    Array<OneD, NekDouble>(fnp0);
                                Vmath::Vcopy(
                                    fnp0,
                                    m_interpTraceI0[set][cnt1]->GetPtr().get() +
                                        tnp0 - 1,
                                    tnp0,
                                    &m_interpEndPtI0[set][cnt1][0],
                                    1);
                            }
                        }
                    }
                    else
                    {
                        m_interpTrace[set][cnt1] = eInterpBothDirs;
                        m_interpTraceI0[set][cnt1] =
                            LibUtilities::PointsManager()[fromPointsKey0]->GetI(
                                toPointsKey0);
                        m_interpTraceI1[set][cnt1] =
                            LibUtilities::PointsManager()[fromPointsKey1]->GetI(
                                toPointsKey1);

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
                                Vmath::Vcopy(
                                    fnp0,
                                    m_interpTraceI0[set][cnt1]->GetPtr().get() +
                                        tnp0 - 1,
                                    tnp0,
                                    &m_interpEndPtI0[set][cnt1][0],
                                    1);
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
    Vmath::Gathr(m_fieldToLocTraceMap.size(),
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
                    Blas::Dgemm('N',
                                'N',
                                tnp,
                                nedges,
                                fnp,
                                1.0,
                                I0->GetPtr().get(),
                                tnp,
                                locedges.get() + cnt,
                                fnp,
                                0.0,
                                tmp.get() + cnt1,
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
                    ASSERTL0(false,
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
    Vmath::Gathr(m_LocTraceToTraceMap[dir].size(),
                 edges.get(),
                 m_LocTraceToTraceMap[dir].get(),
                 tmp.get());

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
                    ASSERTL0(false,
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
            ;

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
    Vmath::Gathr(m_LocTraceToTraceMap[dir].size(),
                 traces.get(),
                 m_LocTraceToTraceMap[dir].get(),
                 tmp.get());

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

                        for(int k = 0; k< tnp1; k++)
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
