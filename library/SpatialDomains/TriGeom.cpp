////////////////////////////////////////////////////////////////////////////////
//
//  File: TriGeom.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/TriGeom.h>
#include <StdRegions/StdNodalTriExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Interp.h>

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/GeomFactors.h>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

TriGeom::TriGeom()
{
    m_shapeType = LibUtilities::eTriangle;
}

TriGeom::TriGeom(const int id,
                 const SegGeomSharedPtr edges[],
                 const CurveSharedPtr curve)
    : Geometry2D(edges[0]->GetVertex(0)->GetCoordim(), curve)
{
    int j;

    m_shapeType = LibUtilities::eTriangle;
    m_globalID = id;

    // Copy the edge shared pointers.
    m_edges.insert(m_edges.begin(), edges, edges + TriGeom::kNedges);
    m_eorient.resize(kNedges);

    for (j = 0; j < kNedges; ++j)
    {
        m_eorient[j] =
            SegGeom::GetEdgeOrientation(*edges[j], *edges[(j + 1) % kNedges]);
        m_verts.push_back(
            edges[j]->GetVertex(m_eorient[j] == StdRegions::eForwards ? 0 : 1));
    }

    m_eorient[2] = m_eorient[2] == StdRegions::eBackwards ?
        StdRegions::eForwards : StdRegions::eBackwards;

    m_coordim = edges[0]->GetVertex(0)->GetCoordim();
    ASSERTL0(m_coordim > 1, "Cannot call function with dim == 1");
}

TriGeom::TriGeom(const TriGeom &in)
    : Geometry2D(in)
{
    // From Geometry
    m_shapeType = in.m_shapeType;

    // From TriFaceComponent
    m_globalID = in.m_globalID;

    // From TriGeom
    m_verts = in.m_verts;
    m_edges = in.m_edges;
    for (int i = 0; i < kNedges; i++)
    {
        m_eorient[i] = in.m_eorient[i];
    }
}

TriGeom::~TriGeom()
{
}

NekDouble TriGeom::v_GetCoord(const int i,
                              const Array<OneD, const NekDouble> &Lcoord)
{
    ASSERTL1(m_state == ePtsFilled, "Geometry is not in physical space");

    Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints());
    m_xmap->BwdTrans(m_coeffs[i], tmp);

    return m_xmap->PhysEvaluate(Lcoord, tmp);
}

StdRegions::Orientation TriGeom::GetFaceOrientation(const TriGeom &face1,
                              const TriGeom &face2, bool doRot, int dir,
                              NekDouble angle, NekDouble tol)
{
    return GetFaceOrientation(face1.m_verts, face2.m_verts,
                              doRot, dir, angle, tol);
}

StdRegions::Orientation TriGeom::GetFaceOrientation(
              const PointGeomVector &face1, const PointGeomVector &face2,
              bool doRot, int dir, NekDouble angle, NekDouble tol)
{
    int i, j, vmap[3] = {-1, -1, -1};

    if(doRot)
    {
        PointGeom rotPt;

        for (i = 0; i < 3; ++i)
        {
            rotPt.Rotate((*face1[i]), dir, angle);
            for (j = 0; j < 3; ++j)
            {
                if (rotPt.dist(*face2[j]) < tol)
                {
                    vmap[j] = i;
                    break;
                }
            }
        }
    }
    else
    {

        NekDouble x, y, z, x1, y1, z1, cx = 0.0, cy = 0.0, cz = 0.0;

        // For periodic faces, we calculate the vector between the centre
        // points of the two faces. (For connected faces this will be
        // zero). We can then use this to determine alignment later in the
        // algorithm.
        for (i = 0; i < 3; ++i)
        {
            cx += (*face2[i])(0) - (*face1[i])(0);
            cy += (*face2[i])(1) - (*face1[i])(1);
            cz += (*face2[i])(2) - (*face1[i])(2);
        }
        cx /= 3;
        cy /= 3;
        cz /= 3;

        // Now construct a mapping which takes us from the vertices of one
        // face to the other. That is, vertex j of face2 corresponds to
        // vertex vmap[j] of face1.
        for (i = 0; i < 3; ++i)
        {
            x = (*face1[i])(0);
            y = (*face1[i])(1);
            z = (*face1[i])(2);
            for (j = 0; j < 3; ++j)
            {
                x1 = (*face2[j])(0) - cx;
                y1 = (*face2[j])(1) - cy;
                z1 = (*face2[j])(2) - cz;
                if (sqrt((x1 - x) * (x1 - x) + (y1 - y) * (y1 - y) +
                         (z1 - z) * (z1 - z)) < 1e-8)
                {
                    vmap[j] = i;
                break;
                }
            }
        }
    }

    if (vmap[1] == (vmap[0] + 1) % 3)
    {
        switch (vmap[0])
        {
        case 0:
            return StdRegions::eDir1FwdDir1_Dir2FwdDir2;
            break;
        case 1:
            return StdRegions::eDir1FwdDir2_Dir2BwdDir1;
            break;
        case 2:
            return StdRegions::eDir1BwdDir1_Dir2BwdDir2;
            break;
        }
    }
    else
    {
        switch (vmap[0])
        {
        case 0:
            return StdRegions::eDir1FwdDir2_Dir2FwdDir1;
            break;
        case 1:
            return StdRegions::eDir1BwdDir1_Dir2FwdDir2;
            break;
        case 2:
            return StdRegions::eDir1BwdDir2_Dir2BwdDir1;
            break;
        }
    }

    ASSERTL0(false, "Unable to determine triangle orientation");
    return StdRegions::eNoOrientation;
}

void TriGeom::v_GenGeomFactors()
{
    if(!m_setupState)
    {
        TriGeom::v_Setup();
    }

    if (m_geomFactorsState != ePtsFilled)
    {
        GeomType Gtype = eRegular;

        TriGeom::v_FillGeom();

        // check to see if expansions are linear
        for (int i = 0; i < m_coordim; ++i)
        {
            if (m_xmap->GetBasisNumModes(0) != 2 ||
                m_xmap->GetBasisNumModes(1) != 2)
            {
                Gtype = eDeformed;
            }
        }

        m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
            Gtype, m_coordim, m_xmap, m_coeffs);

        m_geomFactorsState = ePtsFilled;
    }
}

/**
 * Note verts and edges are listed according to anticlockwise
 * convention but points in _coeffs have to be in array format from
 * left to right.
 */
void TriGeom::v_FillGeom()
{
    // check to see if geometry structure is already filled
    if (m_state == ePtsFilled)
    {
        return;
    }

    int i, j, k;
    int nEdgeCoeffs = m_xmap->GetEdgeNcoeffs(0);

    if (m_curve)
    {
        int pdim = LibUtilities::PointsManager()[LibUtilities::PointsKey(
                                                     2, m_curve->m_ptype)]
                       ->GetPointsDim();

        // Deal with 2D points type separately
        // (e.g. electrostatic or Fekete points) to 1D tensor
        // product.
        if (pdim == 2)
        {
            int N = m_curve->m_points.size();
            int nEdgePts =
                (-1 + (int)sqrt(static_cast<NekDouble>(8 * N + 1))) / 2;

            ASSERTL0(nEdgePts * (nEdgePts + 1) / 2 == N,
                     "NUMPOINTS should be a triangle number for"
                     " triangle curved face " +
                         boost::lexical_cast<string>(m_globalID));

            // Sanity check 1: are curved vertices consistent with
            // triangle vertices?
            for (i = 0; i < 3; ++i)
            {
                NekDouble dist = m_verts[i]->dist(*(m_curve->m_points[i]));
                if (dist > NekConstants::kVertexTheSameDouble)
                {
                    std::stringstream ss;
                    ss << "Curved vertex " << i << " of triangle " << m_globalID
                       << " is separated from expansion vertex by"
                       << " more than " << NekConstants::kVertexTheSameDouble
                       << " (dist = " << dist << ")";
                    NEKERROR(ErrorUtil::ewarning, ss.str().c_str());
                }
            }

            // Sanity check 2: are curved edges from the face curvature
            // consistent with curved edges?
            for (i = 0; i < kNedges; ++i)
            {
                CurveSharedPtr edgeCurve = m_edges[i]->GetCurve();

                ASSERTL0(edgeCurve->m_points.size() == nEdgePts,
                         "Number of edge points does not correspond "
                         "to number of face points in triangle " +
                             boost::lexical_cast<string>(m_globalID));

                const int offset = 3 + i * (nEdgePts - 2);
                NekDouble maxDist = 0.0;

                // Account for different ordering of nodal coordinates
                // vs. Cartesian ordering of element.
                StdRegions::Orientation orient = m_eorient[i];

                if (i == 2)
                {
                    orient = orient == StdRegions::eForwards ?
                        StdRegions::eBackwards : StdRegions::eForwards;
                }

                if (orient == StdRegions::eForwards)
                {
                    for (j = 0; j < nEdgePts - 2; ++j)
                    {
                        NekDouble dist = m_curve->m_points[offset + j]->dist(
                            *(edgeCurve->m_points[j + 1]));
                        maxDist = dist > maxDist ? dist : maxDist;
                    }
                }
                else
                {
                    for (j = 0; j < nEdgePts - 2; ++j)
                    {
                        NekDouble dist = m_curve->m_points[offset + j]->dist(
                            *(edgeCurve->m_points[nEdgePts - 2 - j]));
                        maxDist = dist > maxDist ? dist : maxDist;
                    }
                }

                if (maxDist > NekConstants::kVertexTheSameDouble)
                {
                    std::stringstream ss;
                    ss << "Curved edge " << i << " of triangle " << m_globalID
                       << " has a point separated from edge interior"
                       << " points by more than "
                       << NekConstants::kVertexTheSameDouble
                       << " (maxdist = " << maxDist << ")";
                    NEKERROR(ErrorUtil::ewarning, ss.str().c_str());
                }
            }

            const LibUtilities::PointsKey P0(
                nEdgePts, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey P1(
                nEdgePts, LibUtilities::eGaussRadauMAlpha1Beta0);
            const LibUtilities::BasisKey T0(
                LibUtilities::eOrtho_A, nEdgePts, P0);
            const LibUtilities::BasisKey T1(
                LibUtilities::eOrtho_B, nEdgePts, P1);
            Array<OneD, NekDouble> phys(
                max(nEdgePts * nEdgePts, m_xmap->GetTotPoints()));
            Array<OneD, NekDouble> tmp(nEdgePts * nEdgePts);

            for (i = 0; i < m_coordim; ++i)
            {
                // Create a StdNodalTriExp.
                StdRegions::StdNodalTriExpSharedPtr t =
                    MemoryManager<StdRegions::StdNodalTriExp>::
                        AllocateSharedPtr(T0, T1, m_curve->m_ptype);

                for (j = 0; j < N; ++j)
                {
                    phys[j] = (m_curve->m_points[j]->GetPtr())[i];
                }

                t->BwdTrans(phys, tmp);

                // Interpolate points to standard region.
                LibUtilities::Interp2D(P0,
                                       P1,
                                       tmp,
                                       m_xmap->GetBasis(0)->GetPointsKey(),
                                       m_xmap->GetBasis(1)->GetPointsKey(),
                                       phys);

                // Forwards transform to get coefficient space.
                m_xmap->FwdTrans(phys, m_coeffs[i]);
            }
        }
        else if (pdim == 1)
        {
            int npts = m_curve->m_points.size();
            int nEdgePts = (int)sqrt(static_cast<NekDouble>(npts));
            Array<OneD, NekDouble> tmp(npts);
            Array<OneD, NekDouble> phys(m_xmap->GetTotPoints());
            LibUtilities::PointsKey curveKey(nEdgePts, m_curve->m_ptype);

            // Sanity checks:
            // - Curved faces should have square number of points;
            // - Each edge should have sqrt(npts) points.
            ASSERTL0(nEdgePts * nEdgePts == npts,
                     "NUMPOINTS should be a square number for"
                     " triangle " +
                         boost::lexical_cast<string>(m_globalID));

            for (i = 0; i < kNedges; ++i)
            {
                ASSERTL0(m_edges[i]->GetXmap()->GetNcoeffs() == nEdgePts,
                         "Number of edge points does not correspond to "
                         "number of face points in triangle " +
                             boost::lexical_cast<string>(m_globalID));
            }

            for (i = 0; i < m_coordim; ++i)
            {
                for (j = 0; j < npts; ++j)
                {
                    tmp[j] = (m_curve->m_points[j]->GetPtr())[i];
                }

                // Interpolate curve points to standard triangle
                // points.
                LibUtilities::Interp2D(curveKey,
                                       curveKey,
                                       tmp,
                                       m_xmap->GetBasis(0)->GetPointsKey(),
                                       m_xmap->GetBasis(1)->GetPointsKey(),
                                       phys);

                // Forwards transform to get coefficient space.
                m_xmap->FwdTrans(phys, m_coeffs[i]);
            }
        }
        else
        {
            ASSERTL0(false,
                     "Only 1D/2D points distributions "
                     "supported.");
        }
    }

    Array<OneD, unsigned int> mapArray(nEdgeCoeffs);
    Array<OneD, int> signArray(nEdgeCoeffs);

    for (i = 0; i < kNedges; i++)
    {
        m_edges[i]->FillGeom();
        m_xmap->GetEdgeToElementMap(i, m_eorient[i], mapArray, signArray);

        nEdgeCoeffs = m_edges[i]->GetXmap()->GetNcoeffs();

        for (j = 0; j < m_coordim; j++)
        {
            for (k = 0; k < nEdgeCoeffs; k++)
            {
                m_coeffs[j][mapArray[k]] =
                    signArray[k] * m_edges[i]->GetCoeffs(j)[k];
            }
        }
    }

    m_state = ePtsFilled;
}

/**
 *
 */
NekDouble TriGeom::v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                  Array<OneD, NekDouble> &Lcoords)
{
    NekDouble resid = 0.0;
    TriGeom::v_FillGeom();

    // calculate local coordinate for coord
    if (GetMetricInfo()->GetGtype() == eRegular)
    {
        NekDouble coords2 = (m_coordim == 3) ? coords[2] : 0.0;
        PointGeom dv1, dv2, norm, orth1, orth2;
        PointGeom xin(m_coordim, 0, coords[0], coords[1], coords2);

        // Calculate edge vectors from 0-1 and 0-2 edges.
        dv1.Sub(*m_verts[1], *m_verts[0]);
        dv2.Sub(*m_verts[2], *m_verts[0]);

        // Obtain normal to plane in which dv1 and dv2 lie
        norm.Mult(dv1, dv2);

        // Obtain vector which are proportional to normal of dv1 and dv2.
        orth1.Mult(norm, dv1);
        orth2.Mult(norm, dv2);

        // Start with vector of desired points minus vertex_0
        xin -= *m_verts[0];

        // Calculate length using L/|dv1| = (x-v0).n1/(dv1.n1) for coordiante 1
        // Then rescale to [-1,1].
        Lcoords[0] = xin.dot(orth2) / dv1.dot(orth2);
        Lcoords[0] = 2 * Lcoords[0] - 1;
        Lcoords[1] = xin.dot(orth1) / dv2.dot(orth1);
        Lcoords[1] = 2 * Lcoords[1] - 1;
    }
    else
    {
        // Determine nearest point of coords  to values in m_xmap
        int npts = m_xmap->GetTotPoints();
        Array<OneD, NekDouble> ptsx(npts), ptsy(npts);
        Array<OneD, NekDouble> tmpx(npts), tmpy(npts);

        m_xmap->BwdTrans(m_coeffs[0], ptsx);
        m_xmap->BwdTrans(m_coeffs[1], ptsy);

        const Array<OneD, const NekDouble> za = m_xmap->GetPoints(0);
        const Array<OneD, const NekDouble> zb = m_xmap->GetPoints(1);

        // guess the first local coords based on nearest point
        Vmath::Sadd(npts, -coords[0], ptsx, 1, tmpx, 1);
        Vmath::Sadd(npts, -coords[1], ptsy, 1, tmpy, 1);
        Vmath::Vmul(npts, tmpx, 1, tmpx, 1, tmpx, 1);
        Vmath::Vvtvp(npts, tmpy, 1, tmpy, 1, tmpx, 1, tmpx, 1);

        int min_i = Vmath::Imin(npts, tmpx, 1);

        Lcoords[0] = za[min_i % za.size()];
        Lcoords[1] = zb[min_i / za.size()];

        // recover cartesian coordinate from collapsed coordinate.
        Lcoords[0] = (1.0 + Lcoords[0]) * (1.0 - Lcoords[1]) / 2 - 1.0;

        // Perform newton iteration to find local coordinates
        NewtonIterationForLocCoord(coords, ptsx, ptsy, Lcoords, resid);
    }
    return resid;
}

bool TriGeom::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                              Array<OneD, NekDouble> &stdCoord,
                              NekDouble tol,
                              NekDouble &resid)
{
    //Rough check if within twice min/max point
    if (GetMetricInfo()->GetGtype() != eRegular)
    {
        if (!MinMaxCheck(gloCoord))
        {
            return false;
        }
    }

    // Convert to the local (eta) coordinates.
    resid = GetLocCoords(gloCoord, stdCoord);

    if (stdCoord[0] >= -(1 + tol) && stdCoord[1] >= -(1 + tol) &&
        stdCoord[0] + stdCoord[1] <= tol)
    {
        return true;
    }

    //Clamp local coords
    ClampLocCoords(stdCoord, tol);

    return false;
}

void TriGeom::v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces)
{
    Geometry::v_Reset(curvedEdges, curvedFaces);
    CurveMap::iterator it = curvedFaces.find(m_globalID);

    if (it != curvedFaces.end())
    {
        m_curve = it->second;
    }

    for (int i = 0; i < 3; ++i)
    {
        m_edges[i]->Reset(curvedEdges, curvedFaces);
    }

    SetUpXmap();
    SetUpCoeffs(m_xmap->GetNcoeffs());
}

void TriGeom::v_Setup()
{
    if(!m_setupState)
    {
        for (int i = 0; i < 3; ++i)
        {
            m_edges[i]->Setup();
        }
        SetUpXmap();
        SetUpCoeffs(m_xmap->GetNcoeffs());
        m_setupState = true;
    }
}

void TriGeom::SetUpXmap()
{
    int order0 = m_edges[0]->GetXmap()->GetBasis(0)->GetNumModes();
    int order1 = max(order0,
                     max(m_edges[1]->GetXmap()->GetBasis(0)->GetNumModes(),
                         m_edges[2]->GetXmap()->GetBasis(0)->GetNumModes()));

    const LibUtilities::BasisKey B0(
        LibUtilities::eModified_A,
        order0,
        LibUtilities::PointsKey(order0+1, LibUtilities::eGaussLobattoLegendre));
    const LibUtilities::BasisKey B1(
        LibUtilities::eModified_B,
        order1,
        LibUtilities::PointsKey(order1,
                                LibUtilities::eGaussRadauMAlpha1Beta0));

    m_xmap = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(B0, B1);
}

}
}
