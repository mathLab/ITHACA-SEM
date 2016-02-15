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

#include <LocalRegions/QuadExp.h>
#include <NekMeshUtils/MeshElements/Quadrilateral.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

LibUtilities::ShapeType Quadrilateral::m_type =
    GetElementFactory().RegisterCreatorFunction(
        LibUtilities::eQuadrilateral, Quadrilateral::create, "Quadrilateral");

/**
 * @brief Create a quadrilateral element.
 */
Quadrilateral::Quadrilateral(ElmtConfig pConf,
                             vector<NodeSharedPtr> pNodeList,
                             vector<int> pTagList)
    : Element(pConf, GetNumNodes(pConf), pNodeList.size())
{
    m_tag     = "Q";
    m_dim     = 2;
    m_taglist = pTagList;
    int n     = m_conf.m_order - 1;

    // Create a map to relate edge nodes to a pair of vertices
    // defining an edge. This is based on the ordering produced by
    // gmsh.
    map<pair<int, int>, int> edgeNodeMap;
    map<pair<int, int>, int>::iterator it;
    edgeNodeMap[pair<int, int>(1, 2)] = 5;
    edgeNodeMap[pair<int, int>(2, 3)] = 5 + n;
    edgeNodeMap[pair<int, int>(3, 4)] = 5 + 2 * n;
    edgeNodeMap[pair<int, int>(4, 1)] = 5 + 3 * n;

    // Add vertices. This logic will determine (in 2D) whether the
    // element is clockwise (sum > 0) or counter-clockwise (sum < 0).
    NekDouble sum = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        int o = (i + 1) % 4;
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
                             pNodeList.begin() + 4 * m_conf.m_order,
                             pNodeList.end());
    }
}

void Quadrilateral::Complete(int order)
{
    LibUtilities::BasisKey C0(
        LibUtilities::eOrtho_A,
        order + 1,
        LibUtilities::PointsKey(order + 1,
                                LibUtilities::eGaussLobattoLegendre));

    SpatialDomains::QuadGeomSharedPtr geom =
        boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(this->GetGeom(3));

    // Create a quad.
    LocalRegions::QuadExpSharedPtr quad =
        MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(C0, C0, geom);

    // Get coordinate array for quadrilateral.
    int nqtot = quad->GetTotPoints();
    Array<OneD, NekDouble> alloc(3 * nqtot);
    Array<OneD, NekDouble> x(alloc);
    Array<OneD, NekDouble> y(alloc + 1 * nqtot);
    Array<OneD, NekDouble> z(alloc + 2 * nqtot);

    quad->GetCoords(x, y, z);

    // Now extract points from the co-ordinate arrays into the edge
    // and face nodes. First, extract edge-interior nodes.
    int edgeMap[4][2] = {{0, 1},
                         {order, order + 1},
                         {nqtot - 1, -1},
                         {order * (order + 1), -order - 1}};

    for (int i = 0; i < 4; ++i)
    {
        int pos = edgeMap[i][0] + edgeMap[i][1];
        m_edge[i]->m_edgeNodes.clear();
        for (int j = 1; j < order; ++j, pos += edgeMap[i][1])
        {
            m_edge[i]->m_edgeNodes.push_back(
                NodeSharedPtr(new Node(0, x[pos], y[pos], z[pos])));
        }
    }

    // Extract face-interior nodes.
    m_volumeNodes.clear();
    for (int i = 1; i < order; ++i)
    {
        int pos = i * (order + 1);
        for (int j = 1; j < order; ++j)
        {
            m_volumeNodes.push_back(
                NodeSharedPtr(new Node(0, x[pos + j], y[pos + j], z[pos + j])));
        }
    }

    m_conf.m_order       = order;
    m_conf.m_faceNodes   = true;
    m_conf.m_volumeNodes = true;
}

SpatialDomains::GeometrySharedPtr Quadrilateral::GetGeom(int coordDim)
{
    SpatialDomains::SegGeomSharedPtr edges[4];
    SpatialDomains::PointGeomSharedPtr verts[4];
    SpatialDomains::QuadGeomSharedPtr ret;

    for (int i = 0; i < 4; ++i)
    {
        edges[i] = m_edge[i]->GetGeom(coordDim);
        verts[i] = m_vertex[i]->GetGeom(coordDim);
    }

    StdRegions::Orientation edgeorient[4] = {
        SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
        SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
        SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[3]),
        SpatialDomains::SegGeom::GetEdgeOrientation(*edges[3], *edges[0])};

    ret = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(
        m_id, verts, edges, edgeorient);

    return ret;
}

/**
 * @brief Return the number of nodes defining a quadrilateral.
 */
unsigned int Quadrilateral::GetNumNodes(ElmtConfig pConf)
{
    int n = pConf.m_order;
    if (!pConf.m_faceNodes)
        return 4 * n;
    else
        return (n + 1) * (n + 1);
}
}
}
