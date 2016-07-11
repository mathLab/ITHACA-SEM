////////////////////////////////////////////////////////////////////////////////
//
//  File: TetMesh.cpp
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
//  Description: tet meshing methods
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/TetMeshing/TetMesh.h>
#include <NekMeshUtils/ExtLibInterface/TetGenInterface.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void TetMesh::Mesh()
{
    if (m_mesh->m_verbose)
        cout << endl << endl << "Tetrahdral mesh generation" << endl;

    TetGenInterfaceSharedPtr tetgen =
        MemoryManager<TetGenInterface>::AllocateSharedPtr();

    map<int, NodeSharedPtr> TetgenIdToNode;
    map<int, int> NodeIdToTetgenId;
    // at this point all nodes are in m_mesh->m_vertexset, but if there is a
    // boundary layer, we dont want all of them, also, tetgen ids must be
    // sequential so there is a map from tetgen id to real nodes

    // build sequentially ordered maps of nodes that exist and there delta value
    // in the octree
    map<int, NekDouble> TetgenIdToDelta;
    vector<Array<OneD, int> >
        surfacetris; // surface mesh connectivity based on tetgenids

    int cnt = 0;
    if (!m_pseudosurface)
    {
        // build surface mesh and node map from all surface elements
        for (int i = 0; i < m_mesh->m_element[2].size(); i++)
        {
            ASSERTL0(m_mesh->m_element[2][i]->GetConf().m_e ==
                         LibUtilities::eTriangle,
                     "quad found in surface mesh with no prism mapping");

            vector<NodeSharedPtr> n = m_mesh->m_element[2][i]->GetVertexList();
            Array<OneD, int> tri(3);
            for (int j = 0; j < n.size(); j++)
            {
                map<int, int>::iterator it;
                it = NodeIdToTetgenId.find(n[j]->m_id);
                if (it == NodeIdToTetgenId.end())
                {
                    tri[j]                       = cnt;
                    NodeIdToTetgenId[n[j]->m_id] = cnt;
                    TetgenIdToNode[cnt]          = n[j];
                    TetgenIdToDelta[cnt] = m_octree->Query(n[j]->GetLoc());
                    cnt++;
                }
                else
                {
                    tri[j] = it->second;
                }
            }
            surfacetris.push_back(tri);
        }
    }
    else
    {
        m_surftopriface = m_blmesh->GetSurfToPri();
        // surface triangles will need to be checked against surftopriface to
        // get the right face
        // all surface elements are sequentially numbered so this should be easy
        // to find in map
        for (int i = 0; i < m_mesh->m_element[2].size(); i++)
        {
            if (m_mesh->m_element[2][i]->GetConf().m_e !=
                LibUtilities::eTriangle)
                continue; // no quads for tetgen

            map<int, FaceSharedPtr>::iterator fit;
            fit = m_surftopriface.find(m_mesh->m_element[2][i]->GetId());
            if (fit == m_surftopriface.end())
            {
                // surface element does not have a correspoding prism, build
                // tetgen surface
                // tri from surface element
                vector<NodeSharedPtr> n =
                    m_mesh->m_element[2][i]->GetVertexList();
                Array<OneD, int> tri(3);
                for (int j = 0; j < n.size(); j++)
                {
                    map<int, int>::iterator it;
                    it = NodeIdToTetgenId.find(n[j]->m_id);
                    if (it == NodeIdToTetgenId.end())
                    {
                        tri[j]                       = cnt;
                        NodeIdToTetgenId[n[j]->m_id] = cnt;
                        TetgenIdToNode[cnt]          = n[j];
                        TetgenIdToDelta[cnt] = m_octree->Query(n[j]->GetLoc());
                        cnt++;
                    }
                    else
                    {
                        tri[j] = it->second;
                    }
                }
                surfacetris.push_back(tri);
            }
            else
            {
                // surface element has a prism on it, build tetgen surface
                // element from the face
                vector<NodeSharedPtr> n = fit->second->m_vertexList;
                Array<OneD, int> tri(3);
                for (int j = 0; j < n.size(); j++)
                {
                    map<int, int>::iterator it;
                    it = NodeIdToTetgenId.find(n[j]->m_id);
                    if (it == NodeIdToTetgenId.end())
                    {
                        tri[j]                       = cnt;
                        NodeIdToTetgenId[n[j]->m_id] = cnt;
                        TetgenIdToNode[cnt]          = n[j];
                        TetgenIdToDelta[cnt] = m_octree->Query(n[j]->GetLoc());
                        cnt++;
                    }
                    else
                    {
                        tri[j] = it->second;
                    }
                }
                surfacetris.push_back(tri);
            }
        }
    }

    if (m_mesh->m_verbose)
    {
        cout << "\tInital Node Count: " << TetgenIdToNode.size() << endl;
    }

    tetgen->InitialMesh(TetgenIdToNode, surfacetris);

    vector<Array<OneD, NekDouble> > newp;
    int ctbefore = TetgenIdToNode.size();
    int newpb;

    do
    {
        newpb = newp.size();
        newp.clear();
        tetgen->GetNewPoints(ctbefore, newp);
        for (int i = 0; i < newp.size(); i++)
        {
            NekDouble d                   = m_octree->Query(newp[i]);
            TetgenIdToDelta[ctbefore + i] = d;
        }
        tetgen->RefineMesh(TetgenIdToDelta);
    } while (newpb != newp.size());

    // make new map of all nodes to build tets.

    tetgen->GetNewPoints(ctbefore, newp);
    for (int i = 0; i < newp.size(); i++)
    {
        NodeSharedPtr n = boost::shared_ptr<Node>(
            new Node(m_mesh->m_numNodes++, newp[i][0], newp[i][1], newp[i][2]));
        TetgenIdToNode[ctbefore + i] = n;
    }

    m_tetconnect = tetgen->Extract();

    // tetgen->freetet();

    // create tets
    for (int i = 0; i < m_tetconnect.size(); i++)
    {
        vector<NodeSharedPtr> n;
        n.push_back(TetgenIdToNode[m_tetconnect[i][0]]);
        n.push_back(TetgenIdToNode[m_tetconnect[i][1]]);
        n.push_back(TetgenIdToNode[m_tetconnect[i][2]]);
        n.push_back(TetgenIdToNode[m_tetconnect[i][3]]);
        ElmtConfig conf(LibUtilities::eTetrahedron, 1, false, false);
        vector<int> tags;
        tags.push_back(0);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTetrahedron, conf, n, tags);

        m_mesh->m_element[3].push_back(E);
    }

    if (m_mesh->m_verbose)
        cout << "\tTets :" << m_tetconnect.size() << endl;
}
}
}
