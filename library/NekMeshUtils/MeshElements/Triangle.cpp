////////////////////////////////////////////////////////////////////////////////
//
//  File: Triangle.cpp
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
//  Description: Mesh triangle object.
//
////////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdNodalTriExp.h>
#include <LocalRegions/TriExp.h>
#include <NekMeshUtils/MeshElements/Triangle.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

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
            swap(m_vertex[1], m_vertex[2]);
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
    SpatialDomains::TriGeomSharedPtr ret;

    for (int i = 0; i < 3; ++i)
    {
        edges[i] = m_edge[i]->GetGeom(coordDim);
    }

    ret = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(
        m_id, edges);
    ret->Setup();

    return ret;
}

StdRegions::Orientation Triangle::GetEdgeOrient(int edgeId, EdgeSharedPtr edge)
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
        ASSERTL1(false, "Edge is not connected to this triangle.");
    }

    return StdRegions::eNoOrientation;
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

void Triangle::GetCurvedNodes(std::vector<NodeSharedPtr> &nodeList) const
{
    int n = m_edge[0]->GetNodeCount();
    nodeList.resize(n * (n + 1) / 2);

    // Populate nodelist
    std::copy(m_vertex.begin(), m_vertex.end(), nodeList.begin());
    for (int i = 0; i < 3; ++i)
    {
        std::copy(m_edge[i]->m_edgeNodes.begin(),
                  m_edge[i]->m_edgeNodes.end(),
                  nodeList.begin() + 3 + i * (n - 2));
        if (m_edge[i]->m_n1 != m_vertex[i])
        {
            // If edge orientation is reversed relative to node ordering, we
            // need to reverse order of nodes.
            std::reverse(nodeList.begin() + 3 + i * (n - 2),
                         nodeList.begin() + 3 + (i + 1) * (n - 2));
        }
    }

    // Copy volume nodes.
    std::copy(m_volumeNodes.begin(),
              m_volumeNodes.end(),
              nodeList.begin() + 3 * (n - 1));
}

void Triangle::MakeOrder(int                                order,
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

    // Triangles of order < 3 have no interior volume points.
    if (order == 1 || order == 2)
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

    Array<OneD, NekDouble> px, py;
    LibUtilities::PointsKey pKey(nPoints, pType);
    ASSERTL1(pKey.GetPointsDim() == 2, "Points distribution must be 2D");
    LibUtilities::PointsManager()[pKey]->GetPoints(px, py);

    Array<OneD, Array<OneD, NekDouble> > phys(coordDim);

    for (int i = 0; i < coordDim; ++i)
    {
        phys[i] = Array<OneD, NekDouble>(xmap->GetTotPoints());
        xmap->BwdTrans(geom->GetCoeffs(i), phys[i]);
    }

    const int nTriPts = nPoints * (nPoints + 1) / 2;
    const int nTriIntPts = (nPoints - 3) * (nPoints - 2) / 2;
    m_volumeNodes.resize(nTriIntPts);

    for (int i = 3 + 3*(nPoints-2), cnt = 0; i < nTriPts; ++i, ++cnt)
    {
        Array<OneD, NekDouble> xp(2);
        xp[0] = px[i];
        xp[1] = py[i];

        Array<OneD, NekDouble> x(3, 0.0);
        for (int j = 0; j < coordDim; ++j)
        {
            x[j] = xmap->PhysEvaluate(xp, phys[j]);
        }

        m_volumeNodes[cnt] = std::shared_ptr<Node>(
            new Node(id++, x[0], x[1], x[2]));
    }

    m_conf.m_order       = order;
    m_conf.m_faceNodes   = true;
    m_conf.m_volumeNodes = false;
}

Array<OneD, NekDouble> Triangle::Normal(bool inward)
{
    Array<OneD, NekDouble> ret(3,0.0);

    ret[0] = (m_vertex[1]->m_y - m_vertex[0]->m_y) * (m_vertex[2]->m_z - m_vertex[0]->m_z) -
             (m_vertex[1]->m_z - m_vertex[0]->m_z) * (m_vertex[2]->m_y - m_vertex[0]->m_y);
    ret[1] = (m_vertex[1]->m_z - m_vertex[0]->m_z) * (m_vertex[2]->m_x - m_vertex[0]->m_x) -
             (m_vertex[1]->m_x - m_vertex[0]->m_x) * (m_vertex[2]->m_z - m_vertex[0]->m_z);
    ret[2] = (m_vertex[1]->m_x - m_vertex[0]->m_x) * (m_vertex[2]->m_y - m_vertex[0]->m_y) -
             (m_vertex[1]->m_y - m_vertex[0]->m_y) * (m_vertex[2]->m_x - m_vertex[0]->m_x);

    NekDouble mt = ret[0] * ret[0] + ret[1] * ret[1] + ret[2] * ret[2];
    mt           = sqrt(mt);

    ret[0] /= mt;
    ret[1] /= mt;
    ret[2] /= mt;

    if(m_parentCAD)
    {
        //has cad so can orientate based on that
        if(m_parentCAD->Orientation() == CADOrientation::eBackwards)
        {
            ret[0] *= -1.0;
            ret[1] *= -1.0;
            ret[2] *= -1.0;
        }

        //by default normals point outwards so if we want inward for BLs
        if(inward)
        {
            ret[0] *= -1.0;
            ret[1] *= -1.0;
            ret[2] *= -1.0;
        }
    }
    return ret;
}

}
}
