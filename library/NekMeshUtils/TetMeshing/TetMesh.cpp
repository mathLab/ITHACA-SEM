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

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void TetMesh::Mesh()
{
    if (m_mesh->m_verbose)
        cout << endl << endl << "Tetrahdral mesh generation" << endl;

    tetgen = MemoryManager<TetGenInterface>::AllocateSharedPtr();

    map<int, NodeSharedPtr> IdToNode;
    // at this point all nodes are in m_mesh->m_vertexset, but if there is a
    // boundary layer, we dont want all of them, also, tetgen ids must be
    // sequential so there is a map from tetgen id to real nodes

    // build sequentially ordered maps of nodes that exist and there delta value
    // in the octree
    map<int, NekDouble> IdToDelta;
    vector<Array<OneD, int> > surfacetris;

    int cnt = 0;
    if(!m_usePSurface)
    {
        NodeSet::iterator it;
        for(it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); it++)
        {
            IdToNode[(*it)->m_id] = *it;
            IdToDelta[(*it)->m_id] = m_octree->Query((*it)->GetLoc());
        }
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
                tri[j] = n[j]->m_id;
            }
            surfacetris.push_back(tri);
        }
    }
    else
    {
        ASSERTL0(false,"logic needs replacing will not work currently");
        /*vector<unsigned int> blsurfs = m_blmesh->GetBLSurfs();
        vector<ElementSharedPtr> Psurf = m_blmesh->GetPsuedoSurf();
        //surface triangles will need to be checked against surftopriface to get the right face
        //all surface elements are sequentially numbered so this should be easy to find in map
        for(int i = 0; i < m_mesh->m_element[2].size(); i++)
        {
            if (m_mesh->m_element[2][i]->GetConf().m_e !=
                LibUtilities::eTriangle)
                continue; // no quads for tetgen

            vector<unsigned int> su;
            su.push_back(m_mesh->m_element[2][i]->CADSurfId);

            vector<unsigned int> inter;

            set_intersection(blsurfs.begin(), blsurfs.end(),
                             su.begin(), su.end(),
                             back_inserter(inter));

            if(inter.size() > 0)
            {
                //dont want this surface tri because its under a prism
                continue;
            }

            //surface element does not have a correspoding prism, build tetgen surface
            //tri from surface element
            vector<NodeSharedPtr> n = m_mesh->m_element[2][i]->GetVertexList();
            Array<OneD, int> tri(3);
            for(int j = 0; j < n.size(); j++)
            {
                map<int, int>::iterator it;
                it = NodeIdToTetgenId.find(n[j]->m_id);
                if(it == NodeIdToTetgenId.end())
                {
                    tri[j] = cnt;
                    NodeIdToTetgenId[n[j]->m_id] = cnt;
                    TetgenIdToNode[cnt] = n[j];
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
        for(int i = 0; i < Psurf.size(); i++)
        {
            vector<NodeSharedPtr> n = Psurf[i]->GetVertexList();
            Array<OneD, int> tri(3);
            for(int j = 0; j < n.size(); j++)
            {
                map<int, int>::iterator it;
                it = NodeIdToTetgenId.find(n[j]->m_id);
                if(it == NodeIdToTetgenId.end())
                {
                    tri[j] = cnt;
                    NodeIdToTetgenId[n[j]->m_id] = cnt;
                    TetgenIdToNode[cnt] = n[j];
                    TetgenIdToDelta[cnt] = m_octree->Query(n[j]->GetLoc());
                    cnt++;
                }
                else
                {
                    tri[j] = it->second;
                }
            }
            surfacetris.push_back(tri);
        }*/
    }

    if (m_mesh->m_verbose)
    {
        cout << "\tInital Node Count: " << IdToNode.size() << endl;
    }

    tetgen->InitialMesh(IdToNode, surfacetris);

    vector<Array<OneD, NekDouble> > newp;
    int ctbefore = IdToNode.size();
    int newpb;

    do
    {
        newpb = newp.size();
        newp.clear();
        tetgen->GetNewPoints(ctbefore, newp);
        for (int i = 0; i < newp.size(); i++)
        {
            NekDouble d                   = m_octree->Query(newp[i]);
            IdToDelta[ctbefore + i] = d;
        }
        tetgen->RefineMesh(IdToDelta);
    } while (newpb != newp.size());

    // make new map of all nodes to build tets.
    newp.clear();
    tetgen->GetNewPoints(ctbefore, newp);
    for (int i = 0; i < newp.size(); i++)
    {
        NodeSharedPtr n = boost::shared_ptr<Node>(
            new Node(ctbefore + i, newp[i][0], newp[i][1], newp[i][2]));
        IdToNode[ctbefore + i] = n;
    }

    m_tetconnect = tetgen->Extract();

    /*m_mesh->m_vertexSet.clear(); m_mesh->m_faceSet.clear(); m_mesh->m_element[2].clear();
    m_mesh->m_edgeSet.clear(); m_mesh->m_element[3].clear();

    newp.clear();
    tetgen->GetNewPoints(0,newp);
    vector<NodeSharedPtr> nodes;
    for (int i = 0; i < newp.size(); i++)
    {
        nodes.push_back(boost::shared_ptr<Node>(
            new Node(i, newp[i][0], newp[i][1], newp[i][2])));
    }*/

    // tetgen->freetet();

    // create tets
    for (int i = 0; i < m_tetconnect.size(); i++)
    {
        vector<NodeSharedPtr> n;
        n.push_back(IdToNode[m_tetconnect[i][0]]);
        n.push_back(IdToNode[m_tetconnect[i][1]]);
        n.push_back(IdToNode[m_tetconnect[i][2]]);
        n.push_back(IdToNode[m_tetconnect[i][3]]);
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
