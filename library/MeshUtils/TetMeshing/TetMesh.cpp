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

    vector<int> stiener;

    map<int, MeshTriSharedPtr>::iterator trit;
    for(trit = Tris.begin(); trit != Tris.end(); trit++)
    {
        int s = trit->second->Getcid();
        Array<OneD, int> ns = trit->second->GetN();
        Array<OneD, NekDouble> uva(2); uva[0] = 0.0; uva[1] = 0.0;
        for(int i = 0; i < 3; i++)
        {
            Array<OneD, NekDouble> uv = Nodes[ns[i]]->GetS(s);
            uva[0] += uv[0]/3.0;
            uva[1] += uv[1]/3.0;
        }

        Array<OneD, NekDouble> P = m_cad->GetSurf(s)->P(uva);
        Array<OneD, NekDouble> N = m_cad->GetSurf(s)->N(uva);
        Array<OneD, NekDouble> NP(3);
        NekDouble d = m_octree->Query(P);
        NP[0] = P[0] + N[0]*d*1.41;
        NP[1] = P[1] + N[1]*d*1.41;
        NP[2] = P[2] + N[2]*d*1.41;

        while(!m_cad->InsideShape(NP))
        {
            NP[0] += N[0]*d*0.5;
            NP[1] += N[1]*d*0.5;
            NP[2] += N[2]*d*0.5;
        }

        MeshNodeSharedPtr n = MemoryManager<MeshNode>::AllocateSharedPtr(
            Nodes.size(), NP[0], NP[1], NP[2]);
        stiener.push_back(Nodes.size());
        Nodes[Nodes.size()]=n;

        nodedelta.push_back(m_octree->Query(NP));
    }

    tetgen->InitialMesh(nodesintris, stiener, Tris, Nodes);

    int c = 1;
    int newpb = 20;

    vector<Array<OneD, NekDouble> > newp;

    while(newpb != newp.size())
    {
        newpb = newp.size();
        newp.clear();
        tetgen->GetNewPoints(nodesintris.size() + stiener.size(), newp);

        vector<NekDouble> newpointdelta;
        for(int i = 0; i < newp.size(); i++)
        {
            NekDouble d = m_octree->Query(newp[i]);
            newpointdelta.push_back(d);
        }

        tetgen->RefineMesh(nodesintris.size() + stiener.size(), nodedelta, Tris, Nodes, newpointdelta);
        c++;
    }

    tetgen->AddNodes(nodesintris.size() + stiener.size(), Nodes);

    tetgen->Extract(numtet, tetconnect);

    //tetgen->freetet();

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
