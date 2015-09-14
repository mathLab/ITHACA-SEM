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

#include <MeshUtils/TetMeshing/TetMesh.h>
#include <MeshUtils/ExtLibInterface/TetGenInterface.h>

using namespace std;
namespace Nektar{
namespace MeshUtils{

void TetMesh::Mesh()
{
    if(m_verbose)
        cout << endl << endl << "Tetrahdral mesh generation" << endl;

    TetGenInterfaceSharedPtr tetgen =
        MemoryManager<TetGenInterface>::AllocateSharedPtr();

    //extract the surface mesh
    m_surfacemesh->Get(Nodes,Edges,Tris);

    vector<int> nodesintris;
    vector<NekDouble> nodedelta;

    map<int, MeshNodeSharedPtr>::iterator nit;
    for(nit = Nodes.begin(); nit != Nodes.end(); nit++)
    {
        vector<int> t = nit->second->GetTris();
        vector<int> e = nit->second->GetEdges();
        if(t.size() > 0 && e.size() > 1)
        {
            nodesintris.push_back(nit->first);
            nodedelta.push_back(m_octree->Query(nit->second->GetLoc()));
        }
    }

    tetgen->InitialMesh(nodesintris, nodedelta, Tris, Nodes);

    int c = 1;
    int newpb = 20;

    vector<Array<OneD, NekDouble> > newp;

    while(newpb != newp.size())
    {
        newpb = newp.size();
        newp.clear();
        tetgen->GetNewPoints(nodesintris, newp);

        vector<NekDouble> newpointdelta;
        for(int i = 0; i < newp.size(); i++)
        {
            NekDouble d = m_octree->Query(newp[i]);
            newpointdelta.push_back(d);
        }

        tetgen->RefineMesh(nodesintris, nodedelta, Tris, Nodes, newpointdelta);
        c++;
    }

    tetgen->AddNodes(nodesintris, Nodes);

    tetgen->Extract(numtet, tetconnect);

    tetgen->freetet();

    //create tets
    for(int i = 0; i < numtet; i++)
    {
        MeshTetSharedPtr t = MemoryManager<MeshTet>::AllocateSharedPtr(
            Tets.size(),tetconnect[i][0],tetconnect[i][1],tetconnect[i][2],
            tetconnect[i][3]);
        Nodes[tetconnect[i][0]]->SetTet(Tets.size());
        Nodes[tetconnect[i][1]]->SetTet(Tets.size());
        Nodes[tetconnect[i][2]]->SetTet(Tets.size());
        Nodes[tetconnect[i][3]]->SetTet(Tets.size());
        Tets[Tets.size()] = t;
    }

    if(m_verbose)
        cout << "\tTets :" << numtet << endl;
}

}
}
