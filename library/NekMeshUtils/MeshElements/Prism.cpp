////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshElements.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description: Mesh manipulation objects.
//
////////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdNodalPrismExp.h>
#include <LocalRegions/PrismExp.h>
#include <SpatialDomains/PrismGeom.h>

#include <NekMeshUtils/MeshElements/Prism.h>
using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

LibUtilities::ShapeType Prism::m_type =
    GetElementFactory().RegisterCreatorFunction(
        LibUtilities::ePrism, Prism::create, "Prism");

/**
 * @brief Create a prism element.
 */
Prism::Prism(ElmtConfig pConf,
             vector<NodeSharedPtr> pNodeList,
             vector<int> pTagList)
    : Element(pConf, GetNumNodes(pConf), pNodeList.size())
{
    m_tag     = "R";
    m_dim     = 3;
    m_taglist = pTagList;
    int n     = m_conf.m_order - 1;

    // Create a map to relate edge nodes to a pair of vertices
    // defining an edge. This is based on the ordering produced by
    // gmsh.
    map<pair<int, int>, int> edgeNodeMap;
    map<pair<int, int>, int>::iterator it;

    // This edge-node map is based on Nektar++ ordering.
    edgeNodeMap[pair<int, int>(1, 2)] = 7;
    edgeNodeMap[pair<int, int>(2, 3)] = 7 + n;
    edgeNodeMap[pair<int, int>(4, 3)] = 7 + 2 * n;
    edgeNodeMap[pair<int, int>(1, 4)] = 7 + 3 * n;
    edgeNodeMap[pair<int, int>(1, 5)] = 7 + 4 * n;
    edgeNodeMap[pair<int, int>(2, 5)] = 7 + 5 * n;
    edgeNodeMap[pair<int, int>(3, 6)] = 7 + 6 * n;
    edgeNodeMap[pair<int, int>(4, 6)] = 7 + 7 * n;
    edgeNodeMap[pair<int, int>(5, 6)] = 7 + 8 * n;

    // Add vertices
    for (int i = 0; i < 6; ++i)
    {
        m_vertex.push_back(pNodeList[i]);
    }

    int eid = 0;
    // Create edges (with corresponding set of edge points)
    for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
    {
        vector<NodeSharedPtr> edgeNodes;
        if (m_conf.m_order > 1)
        {
            for (int j = it->second; j < it->second + n; ++j)
            {
                edgeNodes.push_back(pNodeList[j - 1]);
            }
        }
        m_edge.push_back(EdgeSharedPtr(new Edge(pNodeList[it->first.first - 1],
                                                pNodeList[it->first.second - 1],
                                                edgeNodes,
                                                m_conf.m_edgeCurveType)));
        m_edge.back()->m_id = eid++;
    }

    if (m_conf.m_reorient)
    {
        OrientPrism();
    }
    else
    {
        m_orientation = 0;
    }

    // Create faces
    int face_ids[5][4] = {
        {0, 1, 2, 3}, {0, 1, 4, -1}, {1, 2, 5, 4}, {3, 2, 5, -1}, {0, 3, 5, 4}};
    int face_edges[5][4];

    int face_offset[5];
    face_offset[0] = 6 + 9 * n;
    for (int j = 0; j < 4; ++j)
    {
        int facenodes      = j % 2 == 0 ? n * n : n * (n - 1) / 2;
        face_offset[j + 1] = face_offset[j] + facenodes;
    }

    for (int j = 0; j < 5; ++j)
    {
        vector<NodeSharedPtr> faceVertices;
        vector<EdgeSharedPtr> faceEdges;
        vector<NodeSharedPtr> faceNodes;
        int nEdge = 3 - (j % 2 - 1);

        for (int k = 0; k < nEdge; ++k)
        {
            faceVertices.push_back(m_vertex[face_ids[j][k]]);
            NodeSharedPtr a = m_vertex[face_ids[j][k]];
            NodeSharedPtr b = m_vertex[face_ids[j][(k + 1) % nEdge]];
            unsigned int i;
            for (i = 0; i < m_edge.size(); ++i)
            {
                if ((m_edge[i]->m_n1->m_id == a->m_id &&
                     m_edge[i]->m_n2->m_id == b->m_id) ||
                    (m_edge[i]->m_n1->m_id == b->m_id &&
                     m_edge[i]->m_n2->m_id == a->m_id))
                {
                    faceEdges.push_back(m_edge[i]);
                    face_edges[j][k] = i;
                    break;
                }
            }

            if (i == m_edge.size())
            {
                face_edges[j][k] = -1;
            }
        }

        if (m_conf.m_faceNodes)
        {
            int face = j, facenodes;

            if (j % 2 == 0)
            {
                facenodes = n * n;
                if (m_orientation == 1)
                {
                    face = (face + 4) % 6;
                }
                else if (m_orientation == 2)
                {
                    face = (face + 2) % 6;
                }
            }
            else
            {
                // TODO: need to rotate these too.
                facenodes = n * (n - 1) / 2;
            }

            for (int i = 0; i < facenodes; ++i)
            {
                faceNodes.push_back(pNodeList[face_offset[face] + i]);
            }
        }
        m_face.push_back(FaceSharedPtr(new Face(
            faceVertices, faceNodes, faceEdges, m_conf.m_faceCurveType)));
    }

    // Re-order edge array to be consistent with Nektar++ ordering.
    vector<EdgeSharedPtr> tmp(9);
    ASSERTL1(face_edges[0][0] != -1, "face_edges[0][0] == -1");
    tmp[0] = m_edge[face_edges[0][0]];
    ASSERTL1(face_edges[0][1] != -1, "face_edges[0][1] == -1");
    tmp[1] = m_edge[face_edges[0][1]];
    ASSERTL1(face_edges[0][2] != -1, "face_edges[0][2] == -1");
    tmp[2] = m_edge[face_edges[0][2]];
    ASSERTL1(face_edges[0][3] != -1, "face_edges[0][3] == -1");
    tmp[3] = m_edge[face_edges[0][3]];
    ASSERTL1(face_edges[1][2] != -1, "face_edges[1][2] == -1");
    tmp[4] = m_edge[face_edges[1][2]];
    ASSERTL1(face_edges[1][1] != -1, "face_edges[1][1] == -1");
    tmp[5] = m_edge[face_edges[1][1]];
    ASSERTL1(face_edges[2][1] != -1, "face_edges[2][1] == -1");
    tmp[6] = m_edge[face_edges[2][1]];
    ASSERTL1(face_edges[3][2] != -1, "face_edges[3][2] == -1");
    tmp[7] = m_edge[face_edges[3][2]];
    ASSERTL1(face_edges[4][2] != -1, "face_edges[4][2] == -1");
    tmp[8] = m_edge[face_edges[4][2]];
    m_edge = tmp;
}

/**
 * @brief Return the number of nodes defining a prism.
 */
unsigned int Prism::GetNumNodes(ElmtConfig pConf)
{
    int n = pConf.m_order;
    if (pConf.m_faceNodes && pConf.m_volumeNodes)
        return (n + 1) * (n + 1) * (n + 2) / 2;
    else if (pConf.m_faceNodes && !pConf.m_volumeNodes)
        return 3 * (n + 1) * (n + 1) + 2 * (n + 1) * (n + 2) / 2 - 9 * (n + 1) +
               6;
    else
        return 9 * (n + 1) - 12;
}

SpatialDomains::GeometrySharedPtr Prism::GetGeom(int coordDim)
{
    SpatialDomains::Geometry2DSharedPtr faces[5];
    SpatialDomains::PrismGeomSharedPtr ret;

    for (int i = 0; i < 5; ++i)
    {
        faces[i] = m_face[i]->GetGeom(coordDim);
    }

    ret = MemoryManager<SpatialDomains::PrismGeom>::AllocateSharedPtr(faces);

    return ret;
}

/**
 * @brief .
 */
void Prism::Complete(int order)
{
    int i, j, pos;

    // Create basis key for a nodal tetrahedron.
    LibUtilities::BasisKey B0(
        LibUtilities::eOrtho_A,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussLobattoLegendre));
    LibUtilities::BasisKey B1(
        LibUtilities::eOrtho_A,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussLobattoLegendre));
    LibUtilities::BasisKey B2(
        LibUtilities::eOrtho_B,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussRadauMAlpha1Beta0));

    // Create a standard nodal prism in order to get the Vandermonde
    // matrix to perform interpolation to nodal points.
    StdRegions::StdNodalPrismExpSharedPtr nodalPrism =
        MemoryManager<StdRegions::StdNodalPrismExp>::AllocateSharedPtr(
            B0, B1, B2, LibUtilities::eNodalPrismEvenlySpaced);

    Array<OneD, NekDouble> x, y, z;
    nodalPrism->GetNodalPoints(x, y, z);

    SpatialDomains::PrismGeomSharedPtr geom =
        boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
            this->GetGeom(3));

    // Create basis key for a prism.
    LibUtilities::BasisKey C0(
        LibUtilities::eOrtho_A,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussLobattoLegendre));
    LibUtilities::BasisKey C1(
        LibUtilities::eOrtho_A,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussLobattoLegendre));
    LibUtilities::BasisKey C2(
        LibUtilities::eOrtho_B,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussRadauMAlpha1Beta0));

    // Create a prism.
    LocalRegions::PrismExpSharedPtr prism =
        MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(
            C0, C1, C2, geom);

    // Get coordinate array for tetrahedron.
    int nqtot = prism->GetTotPoints();
    Array<OneD, NekDouble> alloc(6 * nqtot);
    Array<OneD, NekDouble> xi(alloc);
    Array<OneD, NekDouble> yi(alloc + nqtot);
    Array<OneD, NekDouble> zi(alloc + 2 * nqtot);
    Array<OneD, NekDouble> xo(alloc + 3 * nqtot);
    Array<OneD, NekDouble> yo(alloc + 4 * nqtot);
    Array<OneD, NekDouble> zo(alloc + 5 * nqtot);
    Array<OneD, NekDouble> tmp;

    prism->GetCoords(xi, yi, zi);

    for (i = 0; i < 3; ++i)
    {
        Array<OneD, NekDouble> coeffs(nodalPrism->GetNcoeffs());
        prism->FwdTrans(alloc + i * nqtot, coeffs);
        // Apply Vandermonde matrix to project onto nodal space.
        nodalPrism->ModalToNodal(coeffs, tmp = alloc + (i + 3) * nqtot);
    }

    // Now extract points from the co-ordinate arrays into the
    // edge/face/volume nodes. First, extract edge-interior nodes.
    for (i = 0; i < 9; ++i)
    {
        pos = 6 + i * (order - 1);
        m_edge[i]->m_edgeNodes.clear();
        for (j = 0; j < order - 1; ++j)
        {
            m_edge[i]->m_edgeNodes.push_back(NodeSharedPtr(
                new Node(0, xo[pos + j], yo[pos + j], zo[pos + j])));
        }
    }

    // Now extract face-interior nodes.
    pos = 6 + 9 * (order - 1);
    for (i = 0; i < 5; ++i)
    {
        int facesize =
            i % 2 ? (order - 2) * (order - 1) / 2 : (order - 1) * (order - 1);
        m_face[i]->m_faceNodes.clear();
        for (j = 0; j < facesize; ++j)
        {
            m_face[i]->m_faceNodes.push_back(NodeSharedPtr(
                new Node(0, xo[pos + j], yo[pos + j], zo[pos + j])));
        }
        pos += facesize;
    }

    // Finally extract volume nodes.
    for (i = pos; i < (order + 1) * (order + 1) * (order + 2) / 2; ++i)
    {
        m_volumeNodes.push_back(
            NodeSharedPtr(new Node(0, xo[i], yo[i], zo[i])));
    }

    m_conf.m_order       = order;
    m_conf.m_faceNodes   = true;
    m_conf.m_volumeNodes = true;
}

/**
 * @brief Orient prism to align degenerate vertices.
 *
 * Orientation of prismatric elements is required so that the singular
 * vertices of triangular faces (which occur as a part of the
 * collapsed co-ordinate system) align. The algorithm is based on that
 * used in T. Warburton's thesis and in the original Nektar source.
 *
 * First the points are re-ordered so that the highest global IDs
 * represent the two singular points of the prism. Then, if necessary,
 * the nodes are rotated either clockwise or counter-clockwise (w.r.t
 * to the p-r plane) to correctly align the prism. The #orientation
 * variable is set to:
 *
 * - 0 if the prism is not rotated;
 * - 1 if the prism is rotated clockwise;
 * - 2 if the prism is rotated counter-clockwise.
 *
 * This is necessary for some input modules (e.g. #InputNek) which add
 * high-order information or bounary conditions to faces.
 */
void Prism::OrientPrism()
{
    int lid[6], gid[6];

    // Re-order vertices.
    for (int i = 0; i < 6; ++i)
    {
        lid[i] = i;
        gid[i] = m_vertex[i]->m_id;
    }

    gid[0] = gid[3] = max(gid[0], gid[3]);
    gid[1] = gid[2] = max(gid[1], gid[2]);
    gid[4] = gid[5] = max(gid[4], gid[5]);

    for (int i = 1; i < 6; ++i)
    {
        if (gid[0] < gid[i])
        {
            swap(gid[i], gid[0]);
            swap(lid[i], lid[0]);
        }
    }

    if (lid[0] == 4 || lid[0] == 5)
    {
        m_orientation = 0;
    }
    else if (lid[0] == 1 || lid[0] == 2)
    {
        // Rotate prism clockwise in p-r plane
        vector<NodeSharedPtr> vertexmap(6);
        vertexmap[0]  = m_vertex[4];
        vertexmap[1]  = m_vertex[0];
        vertexmap[2]  = m_vertex[3];
        vertexmap[3]  = m_vertex[5];
        vertexmap[4]  = m_vertex[1];
        vertexmap[5]  = m_vertex[2];
        m_vertex      = vertexmap;
        m_orientation = 1;
    }
    else if (lid[0] == 0 || lid[0] == 3)
    {
        // Rotate prism counter-clockwise in p-r plane
        vector<NodeSharedPtr> vertexmap(6);
        vertexmap[0]  = m_vertex[1];
        vertexmap[1]  = m_vertex[4];
        vertexmap[2]  = m_vertex[5];
        vertexmap[3]  = m_vertex[2];
        vertexmap[4]  = m_vertex[0];
        vertexmap[5]  = m_vertex[3];
        m_vertex      = vertexmap;
        m_orientation = 2;
    }
    else
    {
        cerr << "Warning: possible prism orientation problem." << endl;
    }
}
}
}
