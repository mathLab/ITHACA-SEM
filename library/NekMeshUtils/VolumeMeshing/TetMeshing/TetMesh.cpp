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

#include <NekMeshUtils/Octree/Octree.h>
#include <NekMeshUtils/VolumeMeshing/TetMeshing/TetMesh.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void TetMesh::Mesh()
{
    if (m_mesh->m_verbose)
    {
        cout << endl << endl << "Tetrahedral mesh generation" << endl;
    }

    vector<Array<OneD, NekDouble>> voidPts = m_mesh->m_cad->GetVoidPoints();
    tetgen = MemoryManager<TetGenInterface>::AllocateSharedPtr(voidPts);

    map<int, NodeSharedPtr> IdToNode;
    map<NodeSharedPtr, int> IdToNodeRev;

    // build sequentially ordered maps of nodes that exist and there delta value
    // in the octree
    map<int, NekDouble> IdToDelta;
    vector<Array<OneD, int> > surfacetris;
    NodeSet alreadyInSurface;

    if(m_surface.size() == 0)
    {
        m_surface = m_mesh->m_element[2];
    }

    int cnt = 0;
    for(int i = 0; i < m_surface.size(); i++)
    {
        vector<NodeSharedPtr> n = m_surface[i]->GetVertexList();
        Array<OneD, int> tri(3);
        for(int j = 0; j < n.size(); j++)
        {
            pair<NodeSet::iterator,bool> testIns =
                alreadyInSurface.insert(n[j]);

            if (testIns.second)
            {
                tri[j] = cnt;
                IdToNode[cnt] = n[j];
                IdToNodeRev[n[j]] = cnt;
                IdToDelta[cnt] = m_mesh->m_octree->Query(n[j]->GetLoc());
                cnt++;
            }
            else
            {
                tri[j] = IdToNodeRev[(*testIns.first)];
            }
        }
        surfacetris.push_back(tri);
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
            NekDouble d = m_mesh->m_octree->Query(newp[i]);
            IdToDelta[ctbefore + i] = d;
        }
        tetgen->RefineMesh(IdToDelta);
    } while (newpb != newp.size());

    // make new map of all nodes to build tets.
    newp.clear();
    tetgen->GetNewPoints(ctbefore, newp);
    for (int i = 0; i < newp.size(); i++)
    {
        NodeSharedPtr n = std::shared_ptr<Node>(
            new Node(ctbefore + i, newp[i][0], newp[i][1], newp[i][2]));
        IdToNode[ctbefore + i] = n;
    }

    m_tetconnect = tetgen->Extract();

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
        tags.push_back(m_id);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTetrahedron, conf, n, tags);

        m_mesh->m_element[3].push_back(E);
    }

    if (m_mesh->m_verbose)
    {
        cout << "\tTets :" << m_tetconnect.size() << endl;
    }
}
}
}
