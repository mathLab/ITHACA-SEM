////////////////////////////////////////////////////////////////////////////////
//
//  File: TetGeom.cpp
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
//  Description: Tetrahedral geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/TetGeom.h>

#include <SpatialDomains/Geometry1D.h>
#include <StdRegions/StdTetExp.h>
#include <SpatialDomains/SegGeom.h>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{
const unsigned int TetGeom::VertexEdgeConnectivity[4][3] = {
    {0, 2, 3}, {0, 1, 4}, {1, 2, 5}, {3, 4, 5}};
const unsigned int TetGeom::VertexFaceConnectivity[4][3] = {
    {0, 1, 3}, {0, 1, 2}, {0, 2, 3}, {1, 2, 3}};
const unsigned int TetGeom::EdgeFaceConnectivity[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 3}, {1, 2}, {2, 3}};

TetGeom::TetGeom()
{
    m_shapeType = LibUtilities::eTetrahedron;
}

TetGeom::TetGeom(int id, const TriGeomSharedPtr faces[])
    : Geometry3D(faces[0]->GetEdge(0)->GetVertex(0)->GetCoordim())
{
    m_shapeType = LibUtilities::eTetrahedron;
    m_globalID = id;

    /// Copy the face shared pointers
    m_faces.insert(m_faces.begin(), faces, faces + TetGeom::kNfaces);

    /// Set up orientation vectors with correct amount of elements.
    m_eorient.resize(kNedges);
    m_forient.resize(kNfaces);

    SetUpLocalEdges();
    SetUpLocalVertices();
    SetUpEdgeOrientation();
    SetUpFaceOrientation();
}

TetGeom::~TetGeom()
{
}

/**
 * @brief Determines if a point specified in global coordinates is
 * located within this tetrahedral geometry and return local caretsian
 * coordinates
 */
bool TetGeom::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                              Array<OneD, NekDouble> &locCoord,
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
    resid = GetLocCoords(gloCoord, locCoord);

    // Check local coordinate is within cartesian bounds.
    if (locCoord[0] >= -(1 + tol) && locCoord[1] >= -(1 + tol) &&
        locCoord[2] >= -(1 + tol) &&
        locCoord[0] + locCoord[1] + locCoord[2] <= -1 + tol)
    {
        return true;
    }

    //Clamp local coords
    ClampLocCoords(locCoord, tol);

    return false;
}

NekDouble TetGeom::v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                  Array<OneD, NekDouble> &Lcoords)
{
    NekDouble ptdist = 1e6;

    // calculate local coordinates (eta) for coord
    if (GetMetricInfo()->GetGtype() == eRegular)
    {
        // Point inside tetrahedron
        PointGeom r(m_coordim, 0, coords[0], coords[1], coords[2]);

        // Edges
        PointGeom er0, e10, e20, e30;
        er0.Sub(r, *m_verts[0]);
        e10.Sub(*m_verts[1], *m_verts[0]);
        e20.Sub(*m_verts[2], *m_verts[0]);
        e30.Sub(*m_verts[3], *m_verts[0]);

        // Cross products (Normal times area)
        PointGeom cp1020, cp2030, cp3010;
        cp1020.Mult(e10, e20);
        cp2030.Mult(e20, e30);
        cp3010.Mult(e30, e10);

        // Barycentric coordinates (relative volume)
        NekDouble V = e30.dot(cp1020); // Tet Volume={(e30)dot(e10)x(e20)}/6
        NekDouble beta = er0.dot(cp2030) / V; // volume1={(er0)dot(e20)x(e30)}/6
        NekDouble gamma = er0.dot(cp3010) / V; // volume1={(er0)dot(e30)x(e10)}/6
        NekDouble delta = er0.dot(cp1020) / V; // volume1={(er0)dot(e10)x(e20)}/6

        // Make tet bigger
        Lcoords[0] = 2.0 * beta - 1.0;
        Lcoords[1] = 2.0 * gamma - 1.0;
        Lcoords[2] = 2.0 * delta - 1.0;

        // Set ptdist to distance to nearest vertex
        for (int i = 0; i < 4; ++i)
        {
            ptdist = min(ptdist, r.dist(*m_verts[i]));
        }
    }
    else
    {
        v_FillGeom();

        // Determine nearest point of coords  to values in m_xmap
        int npts = m_xmap->GetTotPoints();
        Array<OneD, NekDouble> ptsx(npts), ptsy(npts), ptsz(npts);
        Array<OneD, NekDouble> tmp1(npts), tmp2(npts);

        m_xmap->BwdTrans(m_coeffs[0], ptsx);
        m_xmap->BwdTrans(m_coeffs[1], ptsy);
        m_xmap->BwdTrans(m_coeffs[2], ptsz);

        const Array<OneD, const NekDouble> za = m_xmap->GetPoints(0);
        const Array<OneD, const NekDouble> zb = m_xmap->GetPoints(1);
        const Array<OneD, const NekDouble> zc = m_xmap->GetPoints(2);

        // guess the first local coords based on nearest point
        Vmath::Sadd(npts, -coords[0], ptsx, 1, tmp1, 1);
        Vmath::Vmul(npts, tmp1, 1, tmp1, 1, tmp1, 1);
        Vmath::Sadd(npts, -coords[1], ptsy, 1, tmp2, 1);
        Vmath::Vvtvp(npts, tmp2, 1, tmp2, 1, tmp1, 1, tmp1, 1);
        Vmath::Sadd(npts, -coords[2], ptsz, 1, tmp2, 1);
        Vmath::Vvtvp(npts, tmp2, 1, tmp2, 1, tmp1, 1, tmp1, 1);

        int min_i = Vmath::Imin(npts, tmp1, 1);

        // distance from coordinate to nearest point for return value.
        ptdist = sqrt(tmp1[min_i]);

        // Get collapsed coordinate
        int qa = za.size(), qb = zb.size();
        Lcoords[2] = zc[min_i / (qa * qb)];
        min_i = min_i % (qa * qb);
        Lcoords[1] = zb[min_i / qa];
        Lcoords[0] = za[min_i % qa];

        // recover cartesian coordinate from collapsed coordinate.
        Lcoords[1] = (1.0 + Lcoords[0]) * (1.0 - Lcoords[2]) / 2 - 1.0;
        Lcoords[0] = (1.0 + Lcoords[0]) * (-Lcoords[1] - Lcoords[2]) / 2 - 1.0;

        // Perform newton iteration to find local coordinates
        NekDouble resid = 0.0;
        NewtonIterationForLocCoord(coords, ptsx, ptsy, ptsz, Lcoords, resid);
    }
    return ptdist;
}

int TetGeom::v_GetDir(const int faceidx, const int facedir) const
{
    if (faceidx == 0)
    {
        return facedir;
    }
    else if (faceidx == 1)
    {
        return 2 * facedir;
    }
    else
    {
        return 1 + facedir;
    }
}

int TetGeom::v_GetVertexEdgeMap(const int i, const int j) const
{
    return VertexEdgeConnectivity[i][j];
}

int TetGeom::v_GetVertexFaceMap(const int i, const int j) const
{
    return VertexFaceConnectivity[i][j];
}

int TetGeom::v_GetEdgeFaceMap(const int i, const int j) const
{
    return EdgeFaceConnectivity[i][j];
}

void TetGeom::SetUpLocalEdges()
{

    // find edge 0
    int i, j;
    unsigned int check;

    SegGeomSharedPtr edge;

    // First set up the 3 bottom edges

    if (m_faces[0]->GetEid(0) != m_faces[1]->GetEid(0))
    {
        std::ostringstream errstrm;
        errstrm << "Local edge 0 (eid=" << m_faces[0]->GetEid(0);
        errstrm << ") on face " << m_faces[0]->GetGlobalID();
        errstrm << " must be the same as local edge 0 (eid="
                << m_faces[1]->GetEid(0);
        errstrm << ") on face " << m_faces[1]->GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }

    int faceConnected;
    for (faceConnected = 1; faceConnected < 4; faceConnected++)
    {
        check = 0;
        for (i = 0; i < 3; i++)
        {
            if ((m_faces[0])->GetEid(i) == (m_faces[faceConnected])->GetEid(0))
            {
                edge = dynamic_pointer_cast<SegGeom>(
                    (m_faces[0])->GetEdge(i));
                m_edges.push_back(edge);
                check++;
            }
        }

        if (check < 1)
        {
            std::ostringstream errstrm;
            errstrm << "Face 0 does not share an edge with first edge of "
                       "adjacent face. Faces ";
            errstrm << (m_faces[0])->GetGlobalID() << ", "
                    << (m_faces[faceConnected])->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
        else if (check > 1)
        {
            std::ostringstream errstrm;
            errstrm << "Connected faces share more than one edge. Faces ";
            errstrm << (m_faces[0])->GetGlobalID() << ", "
                    << (m_faces[faceConnected])->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
    }

    // Then, set up the 3 vertical edges
    check = 0;
    for (i = 0; i < 3; i++) // Set up the vertical edge :face(1) and face(3)
    {
        for (j = 0; j < 3; j++)
        {
            if ((m_faces[1])->GetEid(i) == (m_faces[3])->GetEid(j))
            {
                edge = dynamic_pointer_cast<SegGeom>(
                    (m_faces[1])->GetEdge(i));
                m_edges.push_back(edge);
                check++;
            }
        }
    }
    if (check < 1)
    {
        std::ostringstream errstrm;
        errstrm << "Connected faces do not share an edge. Faces ";
        errstrm << (m_faces[1])->GetGlobalID() << ", "
                << (m_faces[3])->GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }
    else if (check > 1)
    {
        std::ostringstream errstrm;
        errstrm << "Connected faces share more than one edge. Faces ";
        errstrm << (m_faces[1])->GetGlobalID() << ", "
                << (m_faces[3])->GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }
    // Set up vertical edges: face(1) through face(3)
    for (faceConnected = 1; faceConnected < 3; faceConnected++)
    {
        check = 0;
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                if ((m_faces[faceConnected])->GetEid(i) ==
                    (m_faces[faceConnected + 1])->GetEid(j))
                {
                    edge = dynamic_pointer_cast<SegGeom>(
                        (m_faces[faceConnected])->GetEdge(i));
                    m_edges.push_back(edge);
                    check++;
                }
            }
        }

        if (check < 1)
        {
            std::ostringstream errstrm;
            errstrm << "Connected faces do not share an edge. Faces ";
            errstrm << (m_faces[faceConnected])->GetGlobalID() << ", "
                    << (m_faces[faceConnected + 1])->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
        else if (check > 1)
        {
            std::ostringstream errstrm;
            errstrm << "Connected faces share more than one edge. Faces ";
            errstrm << (m_faces[faceConnected])->GetGlobalID() << ", "
                    << (m_faces[faceConnected + 1])->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
    }
}

void TetGeom::SetUpLocalVertices()
{

    // Set up the first 2 vertices (i.e. vertex 0,1)
    if ((m_edges[0]->GetVid(0) == m_edges[1]->GetVid(0)) ||
        (m_edges[0]->GetVid(0) == m_edges[1]->GetVid(1)))
    {
        m_verts.push_back(m_edges[0]->GetVertex(1));
        m_verts.push_back(m_edges[0]->GetVertex(0));
    }
    else if ((m_edges[0]->GetVid(1) == m_edges[1]->GetVid(0)) ||
             (m_edges[0]->GetVid(1) == m_edges[1]->GetVid(1)))
    {
        m_verts.push_back(m_edges[0]->GetVertex(0));
        m_verts.push_back(m_edges[0]->GetVertex(1));
    }
    else
    {
        std::ostringstream errstrm;
        errstrm << "Connected edges do not share a vertex. Edges ";
        errstrm << m_edges[0]->GetGlobalID() << ", "
                << m_edges[1]->GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }

    // set up the other bottom vertices (i.e. vertex 2)
    for (int i = 1; i < 2; i++)
    {
        if (m_edges[i]->GetVid(0) == m_verts[i]->GetGlobalID())
        {
            m_verts.push_back(m_edges[i]->GetVertex(1));
        }
        else if (m_edges[i]->GetVid(1) == m_verts[i]->GetGlobalID())
        {
            m_verts.push_back(m_edges[i]->GetVertex(0));
        }
        else
        {
            std::ostringstream errstrm;
            errstrm << "Connected edges do not share a vertex. Edges ";
            errstrm << m_edges[i]->GetGlobalID() << ", "
                    << m_edges[i - 1]->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
    }

    // set up top vertex
    if (m_edges[3]->GetVid(0) == m_verts[0]->GetGlobalID())
    {
        m_verts.push_back(m_edges[3]->GetVertex(1));
    }
    else
    {
        m_verts.push_back(m_edges[3]->GetVertex(0));
    }

    // Check the other edges match up.
    int check = 0;
    for (int i = 4; i < 6; ++i)
    {
        if ((m_edges[i]->GetVid(0) == m_verts[i - 3]->GetGlobalID() &&
             m_edges[i]->GetVid(1) == m_verts[3]->GetGlobalID()) ||
            (m_edges[i]->GetVid(1) == m_verts[i - 3]->GetGlobalID() &&
             m_edges[i]->GetVid(0) == m_verts[3]->GetGlobalID()))
        {
            check++;
        }
    }
    if (check != 2)
    {
        std::ostringstream errstrm;
        errstrm << "Connected edges do not share a vertex. Edges ";
        errstrm << m_edges[3]->GetGlobalID() << ", "
                << m_edges[2]->GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }
}

void TetGeom::SetUpEdgeOrientation()
{

    // This 2D array holds the local id's of all the vertices
    // for every edge. For every edge, they are ordered to what we
    // define as being Forwards
    const unsigned int edgeVerts[kNedges][2] = {
        {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};

    int i;
    for (i = 0; i < kNedges; i++)
    {
        if (m_edges[i]->GetVid(0) == m_verts[edgeVerts[i][0]]->GetGlobalID())
        {
            m_eorient[i] = StdRegions::eForwards;
        }
        else if (m_edges[i]->GetVid(0) == m_verts[edgeVerts[i][1]]->GetGlobalID())
        {
            m_eorient[i] = StdRegions::eBackwards;
        }
        else
        {
            ASSERTL0(false, "Could not find matching vertex for the edge");
        }
    }
}

void TetGeom::SetUpFaceOrientation()
{

    int f, i;

    // These arrays represent the vector of the A and B
    // coordinate of the local elemental coordinate system
    // where A corresponds with the coordinate direction xi_i
    // with the lowest index i (for that particular face)
    // Coordinate 'B' then corresponds to the other local
    // coordinate (i.e. with the highest index)
    Array<OneD, NekDouble> elementAaxis(m_coordim);
    Array<OneD, NekDouble> elementBaxis(m_coordim);

    // These arrays correspond to the local coordinate
    // system of the face itself (i.e. the Geometry2D)
    // faceAaxis correspond to the xi_0 axis
    // faceBaxis correspond to the xi_1 axis
    Array<OneD, NekDouble> faceAaxis(m_coordim);
    Array<OneD, NekDouble> faceBaxis(m_coordim);

    // This is the base vertex of the face (i.e. the Geometry2D)
    // This corresponds to thevertex with local ID 0 of the
    // Geometry2D
    unsigned int baseVertex;

    // The lenght of the vectors above
    NekDouble elementAaxis_length;
    NekDouble elementBaxis_length;
    NekDouble faceAaxis_length;
    NekDouble faceBaxis_length;

    // This 2D array holds the local id's of all the vertices
    // for every face. For every face, they are ordered in such
    // a way that the implementation below allows a unified approach
    // for all faces.
    const unsigned int faceVerts[kNfaces][TriGeom::kNverts] = {
        {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};

    NekDouble dotproduct1 = 0.0;
    NekDouble dotproduct2 = 0.0;

    unsigned int orientation;

    // Loop over all the faces to set up the orientation
    for (f = 0; f < kNqfaces + kNtfaces; f++)
    {
        // initialisation
        elementAaxis_length = 0.0;
        elementBaxis_length = 0.0;
        faceAaxis_length = 0.0;
        faceBaxis_length = 0.0;

        dotproduct1 = 0.0;
        dotproduct2 = 0.0;

        baseVertex = m_faces[f]->GetVid(0);

        // We are going to construct the vectors representing the A and B axis
        // of every face. These vectors will be constructed as a
        // vector-representation
        // of the edges of the face. However, for both coordinate directions, we
        // can
        // represent the vectors by two different edges. That's why we need to
        // make sure that
        // we pick the edge to which the baseVertex of the
        // Geometry2D-representation of the face
        // belongs...

        // Compute the length of edges on a base-face
        if (baseVertex == m_verts[faceVerts[f][0]]->GetGlobalID())
        {
            for (i = 0; i < m_coordim; i++)
            {
                elementAaxis[i] = (*m_verts[faceVerts[f][1]])[i] -
                                  (*m_verts[faceVerts[f][0]])[i];
                elementBaxis[i] = (*m_verts[faceVerts[f][2]])[i] -
                                  (*m_verts[faceVerts[f][0]])[i];
            }
        }
        else if (baseVertex == m_verts[faceVerts[f][1]]->GetGlobalID())
        {
            for (i = 0; i < m_coordim; i++)
            {
                elementAaxis[i] = (*m_verts[faceVerts[f][1]])[i] -
                                  (*m_verts[faceVerts[f][0]])[i];
                elementBaxis[i] = (*m_verts[faceVerts[f][2]])[i] -
                                  (*m_verts[faceVerts[f][1]])[i];
            }
        }
        else if (baseVertex == m_verts[faceVerts[f][2]]->GetGlobalID())
        {
            for (i = 0; i < m_coordim; i++)
            {
                elementAaxis[i] = (*m_verts[faceVerts[f][1]])[i] -
                                  (*m_verts[faceVerts[f][2]])[i];
                elementBaxis[i] = (*m_verts[faceVerts[f][2]])[i] -
                                  (*m_verts[faceVerts[f][0]])[i];
            }
        }
        else
        {
            ASSERTL0(false, "Could not find matching vertex for the face");
        }

        // Now, construct the edge-vectors of the local coordinates of
        // the Geometry2D-representation of the face
        for (i = 0; i < m_coordim; i++)
        {
            faceAaxis[i] =
                (*m_faces[f]->GetVertex(1))[i] - (*m_faces[f]->GetVertex(0))[i];
            faceBaxis[i] =
                (*m_faces[f]->GetVertex(2))[i] - (*m_faces[f]->GetVertex(0))[i];

            elementAaxis_length += pow(elementAaxis[i], 2);
            elementBaxis_length += pow(elementBaxis[i], 2);
            faceAaxis_length += pow(faceAaxis[i], 2);
            faceBaxis_length += pow(faceBaxis[i], 2);
        }

        elementAaxis_length = sqrt(elementAaxis_length);
        elementBaxis_length = sqrt(elementBaxis_length);
        faceAaxis_length = sqrt(faceAaxis_length);
        faceBaxis_length = sqrt(faceBaxis_length);

        // Calculate the inner product of both the A-axis
        // (i.e. Elemental A axis and face A axis)
        for (i = 0; i < m_coordim; i++)
        {
            dotproduct1 += elementAaxis[i] * faceAaxis[i];
        }

        NekDouble norm = fabs(dotproduct1) / elementAaxis_length / faceAaxis_length;
        orientation = 0;

        // if the innerproduct is equal to the (absolute value of the ) products
        // of the lengths of both vectors, then, the coordinate systems will NOT
        // be transposed
        if (fabs(norm - 1.0) < NekConstants::kNekZeroTol)
        {
            // if the inner product is negative, both A-axis point
            // in reverse direction
            if (dotproduct1 < 0.0)
            {
                orientation += 2;
            }

            // calculate the inner product of both B-axis
            for (i = 0; i < m_coordim; i++)
            {
                dotproduct2 += elementBaxis[i] * faceBaxis[i];
            }

            norm = fabs(dotproduct2) / elementBaxis_length / faceBaxis_length;

            // check that both these axis are indeed parallel
            ASSERTL1(fabs(norm - 1.0) < NekConstants::kNekZeroTol,
                     "These vectors should be parallel");

            // if the inner product is negative, both B-axis point
            // in reverse direction
            if (dotproduct2 < 0.0)
            {
                orientation++;
            }
        }
        // The coordinate systems are transposed
        else
        {
            orientation = 4;

            // Calculate the inner product between the elemental A-axis
            // and the B-axis of the face (which are now the corresponding axis)
            dotproduct1 = 0.0;
            for (i = 0; i < m_coordim; i++)
            {
                dotproduct1 += elementAaxis[i] * faceBaxis[i];
            }

            norm = fabs(dotproduct1) / elementAaxis_length / faceBaxis_length;

            // check that both these axis are indeed parallel
            ASSERTL1(fabs(norm - 1.0) < NekConstants::kNekZeroTol,
                     "These vectors should be parallel");

            // if the result is negative, both axis point in reverse
            // directions
            if (dotproduct1 < 0.0)
            {
                orientation += 2;
            }

            // Do the same for the other two corresponding axis
            dotproduct2 = 0.0;
            for (i = 0; i < m_coordim; i++)
            {
                dotproduct2 += elementBaxis[i] * faceAaxis[i];
            }

            norm = fabs(dotproduct2) / elementBaxis_length / faceAaxis_length;

            // check that both these axis are indeed parallel
            ASSERTL1(fabs(norm - 1.0) < NekConstants::kNekZeroTol,
                     "These vectors should be parallel");

            if (dotproduct2 < 0.0)
            {
                orientation++;
            }
        }

        orientation = orientation + 5;

        // Fill the m_forient array
        m_forient[f] = (StdRegions::Orientation)orientation;
    }
}

void TetGeom::v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces)
{
    Geometry::v_Reset(curvedEdges, curvedFaces);

    for (int i = 0; i < 4; ++i)
    {
        m_faces[i]->Reset(curvedEdges, curvedFaces);
    }
}

void TetGeom::v_Setup()
{
    if(!m_setupState)
    {
        for (int i = 0; i < 4; ++i)
        {
            m_faces[i]->Setup();
        }
        SetUpXmap();
        SetUpCoeffs(m_xmap->GetNcoeffs());
        m_setupState = true;
    }
}

/**
 * Generate the geometry factors for this element.
 */
void TetGeom::v_GenGeomFactors()
{
    if(!m_setupState)
    {
        TetGeom::v_Setup();
    }

    if (m_geomFactorsState != ePtsFilled)
    {
        GeomType Gtype = eRegular;

        v_FillGeom();

        // check to see if expansions are linear
        for (int i = 0; i < m_coordim; ++i)
        {
            if (m_xmap->GetBasisNumModes(0) != 2 ||
                m_xmap->GetBasisNumModes(1) != 2 ||
                m_xmap->GetBasisNumModes(2) != 2)
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
 * @brief Set up the #m_xmap object by determining the order of each
 * direction from derived faces.
 */
void TetGeom::SetUpXmap()
{
    vector<int> tmp;
    tmp.push_back(m_faces[0]->GetXmap()->GetEdgeNcoeffs(0));
    int order0 = *max_element(tmp.begin(), tmp.end());

    tmp.clear();
    tmp.push_back(order0);
    tmp.push_back(m_faces[0]->GetXmap()->GetEdgeNcoeffs(1));
    tmp.push_back(m_faces[0]->GetXmap()->GetEdgeNcoeffs(2));
    int order1 = *max_element(tmp.begin(), tmp.end());

    tmp.clear();
    tmp.push_back(order0);
    tmp.push_back(order1);
    tmp.push_back(m_faces[1]->GetXmap()->GetEdgeNcoeffs(1));
    tmp.push_back(m_faces[1]->GetXmap()->GetEdgeNcoeffs(2));
    tmp.push_back(m_faces[3]->GetXmap()->GetEdgeNcoeffs(1));
    int order2 = *max_element(tmp.begin(), tmp.end());

    const LibUtilities::BasisKey A(
        LibUtilities::eModified_A,
        order0,
        LibUtilities::PointsKey(order0+1, LibUtilities::eGaussLobattoLegendre));
    const LibUtilities::BasisKey B(
        LibUtilities::eModified_B,
        order1,
        LibUtilities::PointsKey(order1,
                                LibUtilities::eGaussRadauMAlpha1Beta0));
    const LibUtilities::BasisKey C(
        LibUtilities::eModified_C,
        order2,
        LibUtilities::PointsKey(order2,
                                LibUtilities::eGaussRadauMAlpha2Beta0));

    m_xmap = MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(A, B, C);
}

}
}
