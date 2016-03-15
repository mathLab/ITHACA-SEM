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

#include <StdRegions/StdNodalTriExp.h>
#include <LocalRegions/TriExp.h>
#include <NekMeshUtils/MeshElements/Triangle.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

LibUtilities::ShapeType Triangle::m_type =
    GetElementFactory().RegisterCreatorFunction(
        LibUtilities::eTriangle, Triangle::create, "Triangle");

/**
 * @brief Create a triangle element.
 */
Triangle::Triangle(ElmtConfig pConf,
                   vector<NodeSharedPtr> pNodeList,
                   vector<int> pTagList)
    : Element(pConf, GetNumNodes(pConf), pNodeList.size())
{
    m_tag       = "T";
    m_dim       = 2;
    m_taglist   = pTagList;
    m_curveType = LibUtilities::eNodalTriEvenlySpaced;
    int n       = m_conf.m_order - 1;

    // Create a map to relate edge nodes to a pair of vertices
    // defining an edge. This is based on the ordering produced by
    // gmsh.
    map<pair<int, int>, int> edgeNodeMap;
    map<pair<int, int>, int>::iterator it;
    edgeNodeMap[pair<int, int>(1, 2)] = 4;
    edgeNodeMap[pair<int, int>(2, 3)] = 4 + n;
    edgeNodeMap[pair<int, int>(3, 1)] = 4 + 2 * n;

    // Add vertices. This logic will determine (in 2D) whether the
    // element is clockwise (sum > 0) or counter-clockwise (sum < 0).
    NekDouble sum = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        int o = (i + 1) % 3;
        m_vertex.push_back(pNodeList[i]);
        sum += (pNodeList[o]->m_x - pNodeList[i]->m_x) *
               (pNodeList[o]->m_y + pNodeList[i]->m_y);
    }

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
    }

    if (pConf.m_reorient)
    {
        if (sum > 0.0)
        {
            reverse(m_edge.begin(), m_edge.end());
        }
    }

    if (m_conf.m_faceNodes)
    {
        m_volumeNodes.insert(m_volumeNodes.begin(),
                             pNodeList.begin() + 3 * m_conf.m_order,
                             pNodeList.end());
    }
}

SpatialDomains::GeometrySharedPtr Triangle::GetGeom(int coordDim)
{
    SpatialDomains::SegGeomSharedPtr edges[3];
    SpatialDomains::PointGeomSharedPtr verts[3];
    SpatialDomains::TriGeomSharedPtr ret;

    for (int i = 0; i < 3; ++i)
    {
        edges[i] = m_edge[i]->GetGeom(coordDim);
        verts[i] = m_vertex[i]->GetGeom(coordDim);
    }

    StdRegions::Orientation edgeorient[3] = {
        SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
        SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
        SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[0])};

    ret = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(
        m_id, verts, edges, edgeorient);

    return ret;
}

/**
 * @brief Return the number of nodes defining a triangle.
 */
unsigned int Triangle::GetNumNodes(ElmtConfig pConf)
{
    int n = pConf.m_order;
    if (!pConf.m_faceNodes)
        return (n + 1) + 2 * (n - 1) + 1;
    else
        return (n + 1) * (n + 2) / 2;
}

void Triangle::Complete(int order)
{
    int i, j;

    // Create basis key for a nodal tetrahedron.
    LibUtilities::BasisKey B0(
        LibUtilities::eOrtho_A,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussLobattoLegendre));
    LibUtilities::BasisKey B1(
        LibUtilities::eOrtho_B,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussRadauMAlpha1Beta0));

    // Create a standard nodal triangle in order to get the
    // Vandermonde matrix to perform interpolation to nodal points.
    StdRegions::StdNodalTriExpSharedPtr nodalTri =
        MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
            B0, B1, LibUtilities::eNodalTriEvenlySpaced);

    SpatialDomains::TriGeomSharedPtr geom =
        boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(this->GetGeom(3));

    // Create basis key for a triangle.
    LibUtilities::BasisKey C0(
        LibUtilities::eOrtho_A,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussLobattoLegendre));
    LibUtilities::BasisKey C1(
        LibUtilities::eOrtho_B,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussRadauMAlpha1Beta0));

    // Create a triangle.
    LocalRegions::TriExpSharedPtr tri =
        MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(C0, C1, geom);

    // Get coordinate array for tetrahedron.
    int nqtot = tri->GetTotPoints();
    Array<OneD, NekDouble> alloc(6 * nqtot);
    Array<OneD, NekDouble> xi(alloc);
    Array<OneD, NekDouble> yi(alloc + nqtot);
    Array<OneD, NekDouble> zi(alloc + 2 * nqtot);
    Array<OneD, NekDouble> xo(alloc + 3 * nqtot);
    Array<OneD, NekDouble> yo(alloc + 4 * nqtot);
    Array<OneD, NekDouble> zo(alloc + 5 * nqtot);
    Array<OneD, NekDouble> tmp;

    tri->GetCoords(xi, yi, zi);

    for (i = 0; i < 3; ++i)
    {
        Array<OneD, NekDouble> coeffs(nodalTri->GetNcoeffs());
        tri->FwdTrans(alloc + i * nqtot, coeffs);
        // Apply Vandermonde matrix to project onto nodal space.
        nodalTri->ModalToNodal(coeffs, tmp = alloc + (i + 3) * nqtot);
    }

    // Now extract points from the co-ordinate arrays into the
    // edge/face/volume nodes. First, extract edge-interior nodes.
    for (i = 0; i < 3; ++i)
    {
        int pos = 3 + i * (order - 1);
        m_edge[i]->m_edgeNodes.clear();
        for (j = 0; j < order - 1; ++j)
        {
            m_edge[i]->m_edgeNodes.push_back(NodeSharedPtr(
                new Node(0, xo[pos + j], yo[pos + j], zo[pos + j])));
        }
    }

    // Extract face-interior nodes.
    int pos = 3 + 3 * (order - 1);
    for (i = pos; i < (order + 1) * (order + 2) / 2; ++i)
    {
        m_volumeNodes.push_back(
            NodeSharedPtr(new Node(0, xo[i], yo[i], zo[i])));
    }

    m_conf.m_order       = order;
    m_conf.m_faceNodes   = true;
    m_conf.m_volumeNodes = true;
}
}
}
