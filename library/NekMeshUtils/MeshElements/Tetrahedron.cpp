////////////////////////////////////////////////////////////////////////////////
//
//  File: Tetrahedron.cpp
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
//  Description: Mesh tet object.
//
////////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdNodalTetExp.h>
#include <LocalRegions/TetExp.h>
#include <SpatialDomains/TetGeom.h>

#include <NekMeshUtils/MeshElements/Tetrahedron.h>
#include <NekMeshUtils/MeshElements/Triangle.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

LibUtilities::ShapeType Tetrahedron::m_type =
    GetElementFactory().RegisterCreatorFunction(
        LibUtilities::eTetrahedron, Tetrahedron::create, "Tetrahedron");

/// Vertex IDs that make up tetrahedron faces.
int Tetrahedron::m_faceIds[4][3] = {
    {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}
};

/**
 * @brief Create a tetrahedron element.
 */
Tetrahedron::Tetrahedron(ElmtConfig pConf,
                         vector<NodeSharedPtr> pNodeList,
                         vector<int> pTagList)
    : Element(pConf, GetNumNodes(pConf), pNodeList.size())
{
    m_tag     = "A";
    m_dim     = 3;
    m_taglist = pTagList;
    int n     = m_conf.m_order - 1;

    // Create a map to relate edge nodes to a pair of vertices
    // defining an edge.
    map<pair<int, int>, int> edgeNodeMap;
    map<pair<int, int>, int>::iterator it;
    edgeNodeMap[pair<int, int>(1, 2)] = 5;
    edgeNodeMap[pair<int, int>(2, 3)] = 5 + n;
    edgeNodeMap[pair<int, int>(1, 3)] = 5 + 2 * n;
    edgeNodeMap[pair<int, int>(1, 4)] = 5 + 3 * n;
    edgeNodeMap[pair<int, int>(2, 4)] = 5 + 4 * n;
    edgeNodeMap[pair<int, int>(3, 4)] = 5 + 5 * n;

    // Add vertices
    for (int i = 0; i < 4; ++i)
    {
        m_vertex.push_back(pNodeList[i]);
    }

    // Create edges (with corresponding set of edge points)
    int eid = 0;
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

    // Reorient the tet to ensure collapsed coordinates align between
    // adjacent elements.
    if (m_conf.m_reorient)
    {
        OrientTet();
    }
    else
    {
        // If we didn't need to orient the tet then set up the
        // orientation map as the identity mapping.
        for (int i = 0; i < 4; ++i)
        {
            m_orientationMap[i] = i;
        }
    }

    // Face-edge IDs
    int face_edges[4][3];

    // Create faces
    for (int j = 0; j < 4; ++j)
    {
        vector<NodeSharedPtr> faceVertices;
        vector<EdgeSharedPtr> faceEdges;
        vector<NodeSharedPtr> faceNodes;

        // Extract the edges for this face. We need to do this because
        // of the reorientation which might have been applied (see the
        // additional note below).
        for (int k = 0; k < 3; ++k)
        {
            faceVertices.push_back(m_vertex[m_faceIds[j][k]]);
            NodeSharedPtr a = m_vertex[m_faceIds[j][k]];
            NodeSharedPtr b = m_vertex[m_faceIds[j][(k + 1) % 3]];
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

        // When face curvature is supplied, it may have been the case
        // that our tetrahedron was reoriented. In this case, faces have
        // different vertex IDs and so we have to rotate the face
        // curvature so that the two align appropriately.
        if (m_conf.m_faceNodes)
        {
            const int nFaceNodes = n * (n - 1) / 2;

            // Get the vertex IDs of whatever face we are processing.
            vector<int> faceVertIds(3);
            faceVertIds[0] = faceVertices[0]->m_id;
            faceVertIds[1] = faceVertices[1]->m_id;
            faceVertIds[2] = faceVertices[2]->m_id;

            // Find out the original face number as we were given it in
            // the constructor using the orientation map.
            int origFace = -1;
            for (int i = 0; i < 4; ++i)
            {
                if (m_orientationMap[i] == j)
                {
                    origFace = i;
                    break;
                }
            }

            ASSERTL0(origFace >= 0, "Couldn't find face");

            // Now get the face nodes for the original face.
            int N = 4 + 6 * n + origFace * nFaceNodes;
            for (int i = 0; i < nFaceNodes; ++i)
            {
                faceNodes.push_back(pNodeList[N + i]);
            }

            // Find the original face vertex IDs.
            vector<int> origFaceIds(3);
            origFaceIds[0] = pNodeList[m_faceIds[origFace][0]]->m_id;
            origFaceIds[1] = pNodeList[m_faceIds[origFace][1]]->m_id;
            origFaceIds[2] = pNodeList[m_faceIds[origFace][2]]->m_id;

            // Construct a HOTriangle object which performs the
            // orientation magically for us.
            HOTriangle<NodeSharedPtr> hoTri(origFaceIds, faceNodes);
            hoTri.Align(faceVertIds);

            // Copy the face nodes back again.
            faceNodes = hoTri.surfVerts;
        }

        m_face.push_back(FaceSharedPtr(new Face(
            faceVertices, faceNodes, faceEdges, m_conf.m_faceCurveType)));
    }

    if (m_conf.m_volumeNodes)
    {
        const int nFaceNodes = n * (n - 1) / 2;
        for (int i = 4 + 6 * n + 4 * nFaceNodes; i < pNodeList.size(); ++i)
        {
            m_volumeNodes.push_back(pNodeList[i]);
        }
    }

    vector<EdgeSharedPtr> tmp(6);
    tmp[0] = m_edge[face_edges[0][0]];
    tmp[1] = m_edge[face_edges[0][1]];
    tmp[2] = m_edge[face_edges[0][2]];
    tmp[3] = m_edge[face_edges[1][2]];
    tmp[4] = m_edge[face_edges[1][1]];
    tmp[5] = m_edge[face_edges[2][1]];
    m_edge = tmp;
}

SpatialDomains::GeometrySharedPtr Tetrahedron::GetGeom(int coordDim)
{
    SpatialDomains::TriGeomSharedPtr tfaces[4];
    SpatialDomains::TetGeomSharedPtr ret;

    for (int i = 0; i < 4; ++i)
    {
        tfaces[i] = boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(
            m_face[i]->GetGeom(coordDim));
    }

    ret = MemoryManager<SpatialDomains::TetGeom>::AllocateSharedPtr(tfaces);

    return ret;
}

StdRegions::Orientation Tetrahedron::GetEdgeOrient(
    int edgeId, EdgeSharedPtr edge)
{
    static int edgeVerts[6][2] = { {0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3} };

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
        ASSERTL1(false, "Edge is not connected to this quadrilateral.");
    }

    return StdRegions::eNoOrientation;
}

void Tetrahedron::MakeOrder(int                                order,
                            SpatialDomains::GeometrySharedPtr  geom,
                            LibUtilities::PointsType           pType,
                            int                                coordDim,
                            int                               &id,
                            bool                               justConfig)
{
    m_conf.m_order = order;
    m_curveType    = pType;
    m_volumeNodes.clear();

    if (order == 1 || order == 2)
    {
        m_conf.m_volumeNodes = m_conf.m_faceNodes = false;
        return;
    }
    else if (order == 2)
    {
        m_conf.m_faceNodes   = true;
        m_conf.m_volumeNodes = false;
        return;
    }
    else if (order == 3)
    {
        m_conf.m_volumeNodes = false;
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

    Array<OneD, NekDouble> px, py, pz;
    LibUtilities::PointsKey pKey(nPoints, pType);
    ASSERTL1(pKey.GetPointsDim() == 3, "Points distribution must be 3D");
    LibUtilities::PointsManager()[pKey]->GetPoints(px, py, pz);

    Array<OneD, Array<OneD, NekDouble> > phys(coordDim);

    for (int i = 0; i < coordDim; ++i)
    {
        phys[i] = Array<OneD, NekDouble>(xmap->GetTotPoints());
        xmap->BwdTrans(geom->GetCoeffs(i), phys[i]);
    }

    const int nTetPts    = nPoints * (nPoints + 1) * (nPoints + 2) / 6;
    const int nTetIntPts = (nPoints - 4) * (nPoints - 3) * (nPoints - 2) / 6;
    m_volumeNodes.resize(nTetIntPts);

    for (int i = nTetPts - nTetIntPts, cnt = 0; i < nTetPts; ++i, ++cnt)
    {
        Array<OneD, NekDouble> xp(3);
        xp[0] = px[i];
        xp[1] = py[i];
        xp[2] = pz[i];

        Array<OneD, NekDouble> x(3, 0.0);
        for (int j = 0; j < coordDim; ++j)
        {
            x[j] = xmap->PhysEvaluate(xp, phys[j]);
        }

        m_volumeNodes[cnt] = boost::shared_ptr<Node>(
            new Node(id++, x[0], x[1], x[2]));
    }
}

/**
 * @brief Return the number of nodes defining a tetrahedron.
 */
unsigned int Tetrahedron::GetNumNodes(ElmtConfig pConf)
{
    int n = pConf.m_order;
    if (pConf.m_volumeNodes && pConf.m_faceNodes)
        return (n + 1) * (n + 2) * (n + 3) / 6;
    else if (!pConf.m_volumeNodes && pConf.m_faceNodes)
        return 4 * (n + 1) * (n + 2) / 2 - 6 * (n + 1) + 4;
    else
        return 6 * (n + 1) - 8;
}

struct TetOrient
{
    TetOrient(vector<int> nid, int fid) : nid(nid), fid(fid)
    {
    }
    vector<int> nid;
    int fid;
};

struct TetOrientHash : std::unary_function<struct TetOrient, std::size_t>
{
    std::size_t operator()(struct TetOrient const &p) const
    {
        return boost::hash_range(p.nid.begin(), p.nid.end());
    }
};
typedef boost::unordered_set<struct TetOrient, TetOrientHash> TetOrientSet;

bool operator==(const struct TetOrient &a, const struct TetOrient &b)
{
    if (a.nid.size() != b.nid.size())
    {
        return false;
    }

    for (int i = 0; i < a.nid.size(); ++i)
    {
        if (a.nid[i] != b.nid[i])
        {
            return false;
        }
    }

    return true;
}

void Tetrahedron::GetCurvedNodes(std::vector<NodeSharedPtr> &nodeList) const
{
    int n = m_edge[0]->GetNodeCount();
    nodeList.resize(n*(n+1)*(n+2)/6);

    nodeList[0] = m_vertex[0];
    nodeList[1] = m_vertex[1];
    nodeList[2] = m_vertex[2];
    nodeList[3] = m_vertex[3];
    int k = 4;

    for(int i = 0; i < 6; i++)
    {
        bool reverseEdge = false;
        if(i < 3)
        {
            reverseEdge = m_edge[i]->m_n1 == m_vertex[i];
        }
        else
        {
            reverseEdge = m_edge[i]->m_n1 == m_vertex[i-3];
        }

        if (reverseEdge)
        {
            for(int j = 0; j < n-2; j++)
            {
                nodeList[k++] = m_edge[i]->m_edgeNodes[j];
            }
        }
        else
        {
            for(int j = n-3; j >= 0; j--)
            {
                nodeList[k++] = m_edge[i]->m_edgeNodes[j];
            }
        }
    }

    vector<vector<int> > ts;
    vector<int> t(3);
    t[0] = m_vertex[0]->m_id;
    t[1] = m_vertex[1]->m_id;
    t[2] = m_vertex[2]->m_id;
    ts.push_back(t);
    t[0] = m_vertex[0]->m_id;
    t[1] = m_vertex[1]->m_id;
    t[2] = m_vertex[3]->m_id;
    ts.push_back(t);
    t[0] = m_vertex[1]->m_id;
    t[1] = m_vertex[2]->m_id;
    t[2] = m_vertex[3]->m_id;
    ts.push_back(t);
    t[0] = m_vertex[0]->m_id;
    t[1] = m_vertex[2]->m_id;
    t[2] = m_vertex[3]->m_id;
    ts.push_back(t);

    for(int i = 0; i < 4; i++)
    {
        vector<int> fcid;
        fcid.push_back(m_face[i]->m_vertexList[0]->m_id);
        fcid.push_back(m_face[i]->m_vertexList[1]->m_id);
        fcid.push_back(m_face[i]->m_vertexList[2]->m_id);

        HOTriangle<NodeSharedPtr> hot(fcid, m_face[i]->m_faceNodes);

        hot.Align(ts[i]);

        std::copy(hot.surfVerts.begin(),
                  hot.surfVerts.end(),
                  nodeList.begin() + k);
        k+= hot.surfVerts.size();
    }

    std::copy(m_volumeNodes.begin(),
              m_volumeNodes.end(),
              nodeList.begin() + k);
}

/**
 * @brief Orient tetrahedron to align degenerate vertices.
 *
 * Orientation of tetrahedral elements is required so that the
 * singular vertices of triangular faces (which occur as a part of the
 * collapsed co-ordinate system) align. The algorithm is based on that
 * used in T. Warburton's thesis and in the original Nektar source.
 *
 * First the vertices are ordered with the highest global vertex at
 * the top degenerate point, and the base degenerate point has second
 * lowest ID. These vertices are swapped if the element is incorrectly
 * oriented.
 */
void Tetrahedron::OrientTet()
{
    TetOrientSet faces;

    // Create a copy of the original vertex ordering. This is used to
    // construct a mapping, #orientationMap, which maps the original
    // face ordering to the new face ordering.
    for (int i = 0; i < 4; ++i)
    {
        vector<int> nodes(3);

        nodes[0] = m_vertex[m_faceIds[i][0]]->m_id;
        nodes[1] = m_vertex[m_faceIds[i][1]]->m_id;
        nodes[2] = m_vertex[m_faceIds[i][2]]->m_id;

        sort(nodes.begin(), nodes.end());
        struct TetOrient faceNodes(nodes, i);
        faces.insert(faceNodes);
    }

    // Store a copy of the original vertex ordering so we can create a
    // permutation map later.
    vector<NodeSharedPtr> origVert = m_vertex;

    // Order vertices with highest global vertex at top degenerate
    // point. Place second highest global vertex at base degenerate
    // point.
    sort(m_vertex.begin(), m_vertex.end());

    // Calculate a.(b x c) if negative, reverse order of
    // non-degenerate points to correctly orientate the tet.

    NekDouble ax = m_vertex[1]->m_x - m_vertex[0]->m_x;
    NekDouble ay = m_vertex[1]->m_y - m_vertex[0]->m_y;
    NekDouble az = m_vertex[1]->m_z - m_vertex[0]->m_z;
    NekDouble bx = m_vertex[2]->m_x - m_vertex[0]->m_x;
    NekDouble by = m_vertex[2]->m_y - m_vertex[0]->m_y;
    NekDouble bz = m_vertex[2]->m_z - m_vertex[0]->m_z;
    NekDouble cx = m_vertex[3]->m_x - m_vertex[0]->m_x;
    NekDouble cy = m_vertex[3]->m_y - m_vertex[0]->m_y;
    NekDouble cz = m_vertex[3]->m_z - m_vertex[0]->m_z;

    NekDouble nx   = (ay * bz - az * by);
    NekDouble ny   = (az * bx - ax * bz);
    NekDouble nz   = (ax * by - ay * bx);
    NekDouble nmag = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= nmag;
    ny /= nmag;
    nz /= nmag;

    // distance of top vertex from base
    NekDouble dist = cx * nx + cy * ny + cz * nz;

    if (fabs(dist) <= 1e-4)
    {
        cerr << "Warning: degenerate tetrahedron, 3rd vertex is = " << dist
             << " from face" << endl;
    }

    if (dist < 0)
    {
        swap(m_vertex[0], m_vertex[1]);
    }

    nx   = (ay * cz - az * cy);
    ny   = (az * cx - ax * cz);
    nz   = (ax * cy - ay * cx);
    nmag = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= nmag;
    ny /= nmag;
    nz /= nmag;

    // distance of top vertex from base
    dist = bx * nx + by * ny + bz * nz;

    if (fabs(dist) <= 1e-4)
    {
        cerr << "Warning: degenerate tetrahedron, 2nd vertex is = " << dist
             << " from face" << endl;
    }

    nx   = (by * cz - bz * cy);
    ny   = (bz * cx - bx * cz);
    nz   = (bx * cy - by * cx);
    nmag = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= nmag;
    ny /= nmag;
    nz /= nmag;

    // distance of top vertex from base
    dist = ax * nx + ay * ny + az * nz;

    if (fabs(dist) <= 1e-4)
    {
        cerr << "Warning: degenerate tetrahedron, 1st vertex is = " << dist
             << " from face" << endl;
    }

    TetOrientSet::iterator it;

    // Search for the face in the original set of face nodes. Then use
    // this to construct the #orientationMap.
    for (int i = 0; i < 4; ++i)
    {
        vector<int> nodes(3);

        nodes[0] = m_vertex[m_faceIds[i][0]]->m_id;
        nodes[1] = m_vertex[m_faceIds[i][1]]->m_id;
        nodes[2] = m_vertex[m_faceIds[i][2]]->m_id;
        sort(nodes.begin(), nodes.end());

        struct TetOrient faceNodes(nodes, 0);

        it                        = faces.find(faceNodes);
        m_orientationMap[it->fid] = i;

        for (int j = 0; j < 4; ++j)
        {
            if (m_vertex[i]->m_id == origVert[j]->m_id)
            {
                m_origVertMap[j] = i;
                break;
            }
        }
    }
}
}
}
