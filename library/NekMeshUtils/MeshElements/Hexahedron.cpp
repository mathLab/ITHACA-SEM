////////////////////////////////////////////////////////////////////////////////
//
//  File: Hexahedron.cpp
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
//  Description: Hexahedral elements
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/HexGeom.h>
#include <NekMeshUtils/MeshElements/Hexahedron.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

LibUtilities::ShapeType Hexahedron::m_type =
    GetElementFactory().RegisterCreatorFunction(
        LibUtilities::eHexahedron, Hexahedron::create, "Hexahedron");

/// Vertex IDs that make up hexahedron faces.
int Hexahedron::m_faceIds[6][4] = {
    {0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5},
    {3, 2, 6, 7}, {0, 3, 7, 4}, {4, 5, 6, 7}
};

Hexahedron::Hexahedron(ElmtConfig pConf,
                       vector<NodeSharedPtr> pNodeList,
                       vector<int> pTagList)
    : Element(pConf, GetNumNodes(pConf), pNodeList.size())
{
    m_tag     = "H";
    m_dim     = 3;
    m_taglist = pTagList;
    int n     = m_conf.m_order - 1;

    // Create a map to relate edge nodes to a pair of vertices defining an edge
    // This is based on the ordering produced by gmsh.
    map<pair<int, int>, int> edgeNodeMap;
    map<pair<int, int>, int>::iterator it;
    edgeNodeMap[pair<int, int>(1, 2)] = 9;
    edgeNodeMap[pair<int, int>(2, 3)] = 9 + n;
    edgeNodeMap[pair<int, int>(3, 4)] = 9 + 2 * n;
    edgeNodeMap[pair<int, int>(4, 1)] = 9 + 3 * n;
    edgeNodeMap[pair<int, int>(1, 5)] = 9 + 4 * n;
    edgeNodeMap[pair<int, int>(2, 6)] = 9 + 5 * n;
    edgeNodeMap[pair<int, int>(3, 7)] = 9 + 6 * n;
    edgeNodeMap[pair<int, int>(4, 8)] = 9 + 7 * n;
    edgeNodeMap[pair<int, int>(5, 6)] = 9 + 8 * n;
    edgeNodeMap[pair<int, int>(6, 7)] = 9 + 9 * n;
    edgeNodeMap[pair<int, int>(7, 8)] = 9 + 10 * n;
    edgeNodeMap[pair<int, int>(8, 5)] = 9 + 11 * n;

    // Add vertices
    for (int i = 0; i < 8; ++i)
    {
        m_vertex.push_back(pNodeList[i]);
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

    // Create faces
    int face_edges[6][4];
    for (int j = 0; j < 6; ++j)
    {
        vector<NodeSharedPtr> faceVertices;
        vector<EdgeSharedPtr> faceEdges;
        vector<NodeSharedPtr> faceNodes;
        for (int k = 0; k < 4; ++k)
        {
            faceVertices.push_back(m_vertex[m_faceIds[j][k]]);
            NodeSharedPtr a = m_vertex[m_faceIds[j][k]];
            NodeSharedPtr b = m_vertex[m_faceIds[j][(k + 1) % 4]];
            for (unsigned int i = 0; i < m_edge.size(); ++i)
            {
                if (((*(m_edge[i]->m_n1) == *a) &&
                     (*(m_edge[i]->m_n2) == *b)) ||
                    ((*(m_edge[i]->m_n1) == *b) && (*(m_edge[i]->m_n2) == *a)))
                {
                    face_edges[j][k] = i;
                    faceEdges.push_back(m_edge[i]);
                    break;
                }
            }
        }

        if (m_conf.m_faceNodes)
        {
            int N = 8 + 12 * n + j * n * n;
            for (int i = 0; i < n * n; ++i)
            {
                faceNodes.push_back(pNodeList[N + i]);
            }
        }
        m_face.push_back(FaceSharedPtr(new Face(
            faceVertices, faceNodes, faceEdges, m_conf.m_faceCurveType)));
    }

    // Reorder edges to be consistent with Nektar++ ordering.
    vector<EdgeSharedPtr> tmp(12);
    tmp[0]  = m_edge[face_edges[0][0]];
    tmp[1]  = m_edge[face_edges[0][1]];
    tmp[2]  = m_edge[face_edges[0][2]];
    tmp[3]  = m_edge[face_edges[0][3]];
    tmp[4]  = m_edge[face_edges[1][3]];
    tmp[5]  = m_edge[face_edges[1][1]];
    tmp[6]  = m_edge[face_edges[2][1]];
    tmp[7]  = m_edge[face_edges[4][1]];
    tmp[8]  = m_edge[face_edges[5][0]];
    tmp[9]  = m_edge[face_edges[5][1]];
    tmp[10] = m_edge[face_edges[5][2]];
    tmp[11] = m_edge[face_edges[5][3]];
    m_edge  = tmp;
}

SpatialDomains::GeometrySharedPtr Hexahedron::GetGeom(int coordDim)
{
    SpatialDomains::QuadGeomSharedPtr faces[6];
    SpatialDomains::HexGeomSharedPtr ret;

    for (int i = 0; i < 6; ++i)
    {
        faces[i] = std::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
            m_face[i]->GetGeom(coordDim));
    }

    ret = MemoryManager<SpatialDomains::HexGeom>::AllocateSharedPtr(
        m_id, faces);

    ret->Setup();
    return ret;
}

StdRegions::Orientation Hexahedron::GetEdgeOrient(
    int edgeId, EdgeSharedPtr edge)
{
    static int edgeVerts[12][2] = { {0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,5},
                                    {2,6}, {3,7}, {4,5}, {5,6}, {6,7}, {7,4} };

    if (edge->m_n1 == m_vertex[edgeVerts[edgeId][0]])
    {
        return StdRegions::eForwards;
    }
    else if (edge->m_n1 == m_vertex[edgeVerts[edgeId][1]])
    {
        return StdRegions::eBackwards;
    }
    else
    {
        ASSERTL1(false, "Edge is not connected to this hexahedron.");
    }

    return StdRegions::eNoOrientation;
}

void Hexahedron::MakeOrder(int                                order,
                           SpatialDomains::GeometrySharedPtr  geom,
                           LibUtilities::PointsType           pType,
                           int                                coordDim,
                           int                               &id,
                           bool                               justConfig)
{
    m_curveType    = pType;
    m_conf.m_order = order;
    m_volumeNodes.clear();

    if (order == 1)
    {
        m_conf.m_volumeNodes = m_conf.m_faceNodes = false;
        return;
    }

    m_conf.m_faceNodes   = true;
    m_conf.m_volumeNodes = true;

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

    int nHexIntPts = (nPoints - 2) * (nPoints - 2) * (nPoints - 2);
    m_volumeNodes.resize(nHexIntPts);

    for (int i = 1, cnt = 0; i < nPoints-1; ++i)
    {
        for (int j = 1; j < nPoints-1; ++j)
        {
            for (int k = 1; k < nPoints-1; ++k, ++cnt)
            {
                Array<OneD, NekDouble> xp(3);
                xp[0] = px[k];
                xp[1] = px[j];
                xp[2] = px[i];

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
}

/**
 * @brief Return the number of nodes defining a hexahedron.
 */
unsigned int Hexahedron::GetNumNodes(ElmtConfig pConf)
{
    int n = pConf.m_order;
    if (pConf.m_faceNodes && pConf.m_volumeNodes)
        return (n + 1) * (n + 1) * (n + 1);
    else if (pConf.m_faceNodes && !pConf.m_volumeNodes)
        return 6 * (n + 1) * (n + 1) - 12 * (n + 1) + 8;
    else
        return 12 * (n + 1) - 16;
}
}
}
