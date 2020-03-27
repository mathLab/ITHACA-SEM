////////////////////////////////////////////////////////////////////////////////
//
//  File: Quadrilateral.cpp
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
//  Description: Mesh quad object.
//
////////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/QuadExp.h>
#include <NekMeshUtils/MeshElements/Quadrilateral.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

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
            swap(m_vertex[1], m_vertex[3]);
        }
    }

    if (m_conf.m_faceNodes)
    {
        m_volumeNodes.insert(m_volumeNodes.begin(),
                             pNodeList.begin() + 4 * m_conf.m_order,
                             pNodeList.end());
    }
}

StdRegions::Orientation Quadrilateral::GetEdgeOrient(
    int edgeId, EdgeSharedPtr edge)
{
    int locVert = edgeId;
    if (edge->m_n1 == m_vertex[locVert])
    {
        return StdRegions::eForwards;
    }
    else if (edge->m_n2 == m_vertex[locVert])
    {
        return StdRegions::eBackwards;
    }
    else
    {
        ASSERTL1(false, "Edge is not connected to this quadrilateral.");
    }

    return StdRegions::eNoOrientation;
}

void Quadrilateral::MakeOrder(int                                order,
                              SpatialDomains::GeometrySharedPtr  geom,
                              LibUtilities::PointsType           pType,
                              int                                coordDim,
                              int                               &id,
                              bool                               justConfig)
{
    m_conf.m_order       = order;
    m_curveType          = pType;
    m_conf.m_volumeNodes = false;
    m_volumeNodes.clear();

    // Quadrilaterals of order == 1 have no interior volume points.
    if (order == 1)
    {
        m_conf.m_faceNodes = false;
        return;
    }

    m_conf.m_faceNodes = true;

    if (justConfig)
    {
        return;
    }

    int nPoints = order + 1;
    StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

    Array<OneD, NekDouble> px;
    LibUtilities::PointsKey pKey(nPoints, pType);
    ASSERTL1(pKey.GetPointsDim() == 1, "Points distribution must be 1D");
    LibUtilities::PointsManager()[pKey]->GetPoints(px);

    Array<OneD, Array<OneD, NekDouble> > phys(coordDim);

    for (int i = 0; i < coordDim; ++i)
    {
        phys[i] = Array<OneD, NekDouble>(xmap->GetTotPoints());
        xmap->BwdTrans(geom->GetCoeffs(i), phys[i]);
    }

    int nQuadIntPts = (nPoints - 2) * (nPoints - 2);
    m_volumeNodes.resize(nQuadIntPts);

    for (int i = 1, cnt = 0; i < nPoints-1; ++i)
    {
        for (int j = 1; j < nPoints-1; ++j, ++cnt)
        {
            Array<OneD, NekDouble> xp(2);
            xp[0] = px[j];
            xp[1] = px[i];

            Array<OneD, NekDouble> x(3, 0.0);
            for (int k = 0; k < coordDim; ++k)
            {
                x[k] = xmap->PhysEvaluate(xp, phys[k]);
            }

            m_volumeNodes[cnt] = std::shared_ptr<Node>(
                new Node(id++, x[0], x[1], x[2]));
        }
    }
}

SpatialDomains::GeometrySharedPtr Quadrilateral::GetGeom(int coordDim)
{
    SpatialDomains::SegGeomSharedPtr edges[4];
    SpatialDomains::QuadGeomSharedPtr ret;

    for (int i = 0; i < 4; ++i)
    {
        edges[i] = m_edge[i]->GetGeom(coordDim);
    }

    ret = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(
        m_id, edges);

    ret->Setup();
    return ret;
}

void Quadrilateral::GetCurvedNodes(std::vector<NodeSharedPtr> &nodeList) const
{
    int n = m_edge[0]->GetNodeCount();
    nodeList.resize(n * n);

    // Write vertices
    nodeList[0]         = m_vertex[0];
    nodeList[n - 1]     = m_vertex[1];
    nodeList[n * n - 1] = m_vertex[2];
    nodeList[n * (n - 1)] = m_vertex[3];

    // Write edge-interior
    int skips[4][2] = {
        {0, 1}, {n - 1, n}, {n * n - 1, -1}, {n * (n - 1), -n}};
    for (int i = 0; i < 4; ++i)
    {
        bool reverseEdge = m_edge[i]->m_n1 == m_vertex[i];

        if (!reverseEdge)
        {
            for (int j = 1; j < n - 1; ++j)
            {
                nodeList[skips[i][0] + j * skips[i][1]] =
                    m_edge[i]->m_edgeNodes[n - 2 - j];
            }
        }
        else
        {
            for (int j = 1; j < n - 1; ++j)
            {
                nodeList[skips[i][0] + j * skips[i][1]] =
                    m_edge[i]->m_edgeNodes[j - 1];
            }
        }
    }

    // Write interior
    for (int i = 1; i < n - 1; ++i)
    {
        for (int j = 1; j < n - 1; ++j)
        {
            nodeList[i * n + j] =
                m_volumeNodes[(i - 1) * (n - 2) + (j - 1)];
        }
    }
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
