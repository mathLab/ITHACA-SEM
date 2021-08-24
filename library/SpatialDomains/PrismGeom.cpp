////////////////////////////////////////////////////////////////////////////////
//
//  File: PrismGeom.cpp
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
//  Description: Prismatic geometry definition.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/Geometry2D.h>
#include <StdRegions/StdPrismExp.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/GeomFactors.h>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

const unsigned int PrismGeom::VertexEdgeConnectivity[6][3] = {
    {0, 3, 4}, {0, 1, 5}, {1, 2, 6}, {2, 3, 7}, {4, 5, 8}, {6, 7, 8}};
const unsigned int PrismGeom::VertexFaceConnectivity[6][3] = {
    {0, 1, 4}, {0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {1, 2, 4}, {2, 3, 4}};
const unsigned int PrismGeom::EdgeFaceConnectivity[9][2] = {
    {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 4}, {1, 2}, {2, 3}, {3, 4}, {2, 4}};

PrismGeom::PrismGeom()
{
    m_shapeType = LibUtilities::ePrism;
}

PrismGeom::PrismGeom(int id, const Geometry2DSharedPtr faces[])
    : Geometry3D(faces[0]->GetEdge(0)->GetVertex(0)->GetCoordim())
{
    m_shapeType = LibUtilities::ePrism;
    m_globalID = id;

    /// Copy the face shared pointers.
    m_faces.insert(m_faces.begin(), faces, faces + PrismGeom::kNfaces);

    /// Set up orientation vectors with correct amount of elements.
    m_eorient.resize(kNedges);
    m_forient.resize(kNfaces);

    /// Set up local objects.
    SetUpLocalEdges();
    SetUpLocalVertices();
    SetUpEdgeOrientation();
    SetUpFaceOrientation();
}

PrismGeom::~PrismGeom()
{
}

int PrismGeom::v_GetDir(const int faceidx, const int facedir) const
{
    if (faceidx == 0)
    {
        return facedir;
    }
    else if (faceidx == 1 || faceidx == 3)
    {
        return 2 * facedir;
    }
    else
    {
        return 1 + facedir;
    }
}

void PrismGeom::v_GenGeomFactors()
{
    if(!m_setupState)
    {
        PrismGeom::v_Setup();
    }

    if (m_geomFactorsState != ePtsFilled)
    {
        int i, f;
        GeomType Gtype = eRegular;

        v_FillGeom();

        // check to see if expansions are linear
        for (i = 0; i < m_coordim; ++i)
        {
            if (m_xmap->GetBasisNumModes(0) != 2 ||
                m_xmap->GetBasisNumModes(1) != 2 ||
                m_xmap->GetBasisNumModes(2) != 2)
            {
                Gtype = eDeformed;
            }
        }

        // check to see if all quadrilateral faces are parallelograms
        if (Gtype == eRegular)
        {
            // Vertex ids of quad faces
            const unsigned int faceVerts[3][4] = {
                {0, 1, 2, 3}, {1, 2, 5, 4}, {0, 3, 5, 4}};

            for (f = 0; f < 3; f++)
            {
                // Ensure each face is a parallelogram? Check this.
                for (i = 0; i < m_coordim; i++)
                {
                    if (fabs((*m_verts[faceVerts[f][0]])(i) -
                             (*m_verts[faceVerts[f][1]])(i) +
                             (*m_verts[faceVerts[f][2]])(i) -
                             (*m_verts[faceVerts[f][3]])(i)) >
                        NekConstants::kNekZeroTol)
                    {
                        Gtype = eDeformed;
                        break;
                    }
                }

                if (Gtype == eDeformed)
                {
                    break;
                }
            }
        }

        m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
            Gtype, m_coordim, m_xmap, m_coeffs);
        m_geomFactorsState = ePtsFilled;
    }
}

int PrismGeom::v_GetVertexEdgeMap(const int i, const int j) const
{
    return VertexEdgeConnectivity[i][j];
}

int PrismGeom::v_GetVertexFaceMap(const int i, const int j) const
{
    return VertexFaceConnectivity[i][j];
}

int PrismGeom::v_GetEdgeFaceMap(const int i, const int j) const
{
    return EdgeFaceConnectivity[i][j];
}

void PrismGeom::SetUpLocalEdges()
{
    // find edge 0
    int i, j;
    unsigned int check;

    SegGeomSharedPtr edge;

    // First set up the 4 bottom edges
    int f; //  Connected face index
    for (f = 1; f < 5; f++)
    {
        int nEdges = m_faces[f]->GetNumEdges();
        check = 0;
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < nEdges; j++)
            {
                if (m_faces[0]->GetEid(i) == m_faces[f]->GetEid(j))
                {
                    edge = dynamic_pointer_cast<SegGeom>(
                        (m_faces[0])->GetEdge(i));
                    m_edges.push_back(edge);
                    check++;
                }
            }
        }

        if (check < 1)
        {
            std::ostringstream errstrm;
            errstrm << "Connected faces do not share an edge. Faces ";
            errstrm << (m_faces[0])->GetGlobalID() << ", "
                    << (m_faces[f])->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
        else if (check > 1)
        {
            std::ostringstream errstrm;
            errstrm << "Connected faces share more than one edge. Faces ";
            errstrm << (m_faces[0])->GetGlobalID() << ", "
                    << (m_faces[f])->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
    }

    // Then, set up the 4 vertical edges
    check = 0;
    for (i = 0; i < 3; i++) // Set up the vertical edge :face(1) and face(4)
    {
        for (j = 0; j < 4; j++)
        {
            if ((m_faces[1])->GetEid(i) == (m_faces[4])->GetEid(j))
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
                << (m_faces[4])->GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }
    else if (check > 1)
    {
        std::ostringstream errstrm;
        errstrm << "Connected faces share more than one edge. Faces ";
        errstrm << (m_faces[1])->GetGlobalID() << ", "
                << (m_faces[4])->GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }
    // Set up vertical edges: face(1) through face(4)
    for (f = 1; f < 4; f++)
    {
        check = 0;
        for (i = 0; i < m_faces[f]->GetNumEdges(); i++)
        {
            for (j = 0; j < m_faces[f + 1]->GetNumEdges(); j++)
            {
                if ((m_faces[f])->GetEid(i) == (m_faces[f + 1])->GetEid(j))
                {
                    edge = dynamic_pointer_cast<SegGeom>(
                        (m_faces[f])->GetEdge(i));
                    m_edges.push_back(edge);
                    check++;
                }
            }
        }

        if (check < 1)
        {
            std::ostringstream errstrm;
            errstrm << "Connected faces do not share an edge. Faces ";
            errstrm << (m_faces[f])->GetGlobalID() << ", "
                    << (m_faces[f + 1])->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
        else if (check > 1)
        {
            std::ostringstream errstrm;
            errstrm << "Connected faces share more than one edge. Faces ";
            errstrm << (m_faces[f])->GetGlobalID() << ", "
                    << (m_faces[f + 1])->GetGlobalID();
            ASSERTL0(false, errstrm.str());
        }
    }

    // Finally, set up the 1 top edge
    check = 0;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            if ((m_faces[2])->GetEid(i) == (m_faces[4])->GetEid(j))
            {
                edge = dynamic_pointer_cast<SegGeom>(
                    (m_faces[2])->GetEdge(i));
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
}

void PrismGeom::SetUpLocalVertices()
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

    // set up the other bottom vertices (i.e. vertex 2,3)
    for (int i = 1; i < 3; i++)
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

    // set up top vertices
    // First, set up vertices 4,5
    if ((m_edges[8]->GetVid(0) == m_edges[4]->GetVid(0)) ||
        (m_edges[8]->GetVid(0) == m_edges[4]->GetVid(1)))
    {
        m_verts.push_back(m_edges[8]->GetVertex(0));
        m_verts.push_back(m_edges[8]->GetVertex(1));
    }
    else if ((m_edges[8]->GetVid(1) == m_edges[4]->GetVid(0)) ||
             (m_edges[8]->GetVid(1) == m_edges[4]->GetVid(1)))
    {
        m_verts.push_back(m_edges[8]->GetVertex(1));
        m_verts.push_back(m_edges[8]->GetVertex(0));
    }
    else
    {
        std::ostringstream errstrm;
        errstrm << "Connected edges do not share a vertex. Edges ";
        errstrm << m_edges[8]->GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }
}

void PrismGeom::SetUpEdgeOrientation()
{

    // This 2D array holds the local id's of all the vertices
    // for every edge. For every edge, they are ordered to what we
    // define as being Forwards
    const unsigned int edgeVerts[kNedges][2] = {
        {0, 1}, {1, 2}, {3, 2}, {0, 3}, {0, 4}, {1, 4}, {2, 5}, {3, 5}, {4, 5}};

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

void PrismGeom::SetUpFaceOrientation()
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
    const unsigned int faceVerts[kNfaces][QuadGeom::kNverts] = {
        {0, 1, 2, 3},
        {0, 1, 4, 0}, // This is triangle requires only three vertices
        {1, 2, 5, 4},
        {3, 2, 5, 0}, // This is triangle requires only three vertices
        {0, 3, 5, 4},
    };

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
        if (f == 1 || f == 3)
        { // Face is a Triangle
            if (baseVertex == m_verts[faceVerts[f][0]]->GetVid())
            {
                for (i = 0; i < m_coordim; i++)
                {
                    elementAaxis[i] = (*m_verts[faceVerts[f][1]])[i] -
                        (*m_verts[faceVerts[f][0]])[i];
                    elementBaxis[i] = (*m_verts[faceVerts[f][2]])[i] -
                        (*m_verts[faceVerts[f][0]])[i];
                }
            }
            else if (baseVertex == m_verts[faceVerts[f][1]]->GetVid())
            {
                for (i = 0; i < m_coordim; i++)
                {
                    elementAaxis[i] = (*m_verts[faceVerts[f][1]])[i] -
                        (*m_verts[faceVerts[f][0]])[i];
                    elementBaxis[i] = (*m_verts[faceVerts[f][2]])[i] -
                        (*m_verts[faceVerts[f][1]])[i];
                }
            }
            else if (baseVertex == m_verts[faceVerts[f][2]]->GetVid())
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
        }
        else
        { // Face is a Quad
            if (baseVertex == m_verts[faceVerts[f][0]]->GetGlobalID())
            {
                for (i = 0; i < m_coordim; i++)
                {
                    elementAaxis[i] = (*m_verts[faceVerts[f][1]])[i] -
                                      (*m_verts[faceVerts[f][0]])[i];
                    elementBaxis[i] = (*m_verts[faceVerts[f][3]])[i] -
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
                    elementAaxis[i] = (*m_verts[faceVerts[f][2]])[i] -
                                      (*m_verts[faceVerts[f][3]])[i];
                    elementBaxis[i] = (*m_verts[faceVerts[f][2]])[i] -
                                      (*m_verts[faceVerts[f][1]])[i];
                }
            }
            else if (baseVertex == m_verts[faceVerts[f][3]]->GetGlobalID())
            {
                for (i = 0; i < m_coordim; i++)
                {
                    elementAaxis[i] = (*m_verts[faceVerts[f][2]])[i] -
                                      (*m_verts[faceVerts[f][3]])[i];
                    elementBaxis[i] = (*m_verts[faceVerts[f][3]])[i] -
                                      (*m_verts[faceVerts[f][0]])[i];
                }
            }
            else
            {
                ASSERTL0(false, "Could not find matching vertex for the face");
            }
        }
        // Now, construct the edge-vectors of the local coordinates of
        // the Geometry2D-representation of the face
        for (i = 0; i < m_coordim; i++)
        {
            int v = m_faces[f]->GetNumVerts() - 1;
            faceAaxis[i] =
                (*m_faces[f]->GetVertex(1))[i] - (*m_faces[f]->GetVertex(0))[i];
            faceBaxis[i] =
                (*m_faces[f]->GetVertex(v))[i] - (*m_faces[f]->GetVertex(0))[i];

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

        orientation = 0;

        // if the innerproduct is equal to the (absolute value of the ) products
        // of the lengths of both vectors, then, the coordinate systems will NOT
        // be transposed
        if (fabs(elementAaxis_length * faceAaxis_length - fabs(dotproduct1)) <
            NekConstants::kNekZeroTol)
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

            ASSERTL1(
                fabs(fabs(dotproduct2 / elementBaxis_length / faceBaxis_length)
                     - 1.0) < NekConstants::kNekZeroTol,
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

            // check that both these axis are indeed parallel
            ASSERTL1(
                fabs(fabs(dotproduct1) / elementAaxis_length / faceBaxis_length
                     - 1.0) < NekConstants::kNekZeroTol,
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

            ASSERTL1(
                fabs(fabs(dotproduct2) / elementBaxis_length / faceAaxis_length
                     - 1.0) < NekConstants::kNekZeroTol,
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

void PrismGeom::v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces)
{
    Geometry::v_Reset(curvedEdges, curvedFaces);

    for (int i = 0; i < 5; ++i)
    {
        m_faces[i]->Reset(curvedEdges, curvedFaces);
    }
}

void PrismGeom::v_Setup()
{
    if(!m_setupState)
    {
        for (int i = 0; i < 5; ++i)
        {
            m_faces[i]->Setup();
        }
        SetUpXmap();
        SetUpCoeffs(m_xmap->GetNcoeffs());
        m_setupState = true;
    }
}

/**
 * @brief Set up the #m_xmap object by determining the order of each
 * direction from derived faces.
 */
void PrismGeom::SetUpXmap()
{
    vector<int> tmp;
    int order0, order1;

    if (m_forient[0] < 9)
    {
        tmp.push_back(m_faces[0]->GetXmap()->GetTraceNcoeffs(0));
        tmp.push_back(m_faces[0]->GetXmap()->GetTraceNcoeffs(2));
        order0 = *max_element(tmp.begin(), tmp.end());
    }
    else
    {
        tmp.push_back(m_faces[0]->GetXmap()->GetTraceNcoeffs(1));
        tmp.push_back(m_faces[0]->GetXmap()->GetTraceNcoeffs(3));
        order0 = *max_element(tmp.begin(), tmp.end());
    }

    if (m_forient[0] < 9)
    {
        tmp.clear();
        tmp.push_back(m_faces[0]->GetXmap()->GetTraceNcoeffs(1));
        tmp.push_back(m_faces[0]->GetXmap()->GetTraceNcoeffs(3));
        tmp.push_back(m_faces[2]->GetXmap()->GetTraceNcoeffs(2));
        order1 = *max_element(tmp.begin(), tmp.end());
    }
    else
    {
        tmp.clear();
        tmp.push_back(m_faces[0]->GetXmap()->GetTraceNcoeffs(0));
        tmp.push_back(m_faces[0]->GetXmap()->GetTraceNcoeffs(2));
        tmp.push_back(m_faces[2]->GetXmap()->GetTraceNcoeffs(2));
        order1 = *max_element(tmp.begin(), tmp.end());
    }

    tmp.clear();
    tmp.push_back(order0);
    tmp.push_back(order1);
    tmp.push_back(m_faces[1]->GetXmap()->GetTraceNcoeffs(1));
    tmp.push_back(m_faces[1]->GetXmap()->GetTraceNcoeffs(2));
    tmp.push_back(m_faces[3]->GetXmap()->GetTraceNcoeffs(1));
    tmp.push_back(m_faces[3]->GetXmap()->GetTraceNcoeffs(2));
    int order2 = *max_element(tmp.begin(), tmp.end());

    const LibUtilities::BasisKey A(
        LibUtilities::eModified_A,
        order0,
        LibUtilities::PointsKey(order0+1, LibUtilities::eGaussLobattoLegendre));
    const LibUtilities::BasisKey B(
        LibUtilities::eModified_A,
        order1,
        LibUtilities::PointsKey(order1+1, LibUtilities::eGaussLobattoLegendre));
    const LibUtilities::BasisKey C(
        LibUtilities::eModified_B,
        order2,
        LibUtilities::PointsKey(order2,
                                LibUtilities::eGaussRadauMAlpha1Beta0));

    m_xmap = MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(A, B, C);
}

}
}
