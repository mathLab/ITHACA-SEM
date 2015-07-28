////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <MeshUtils/TetMesh.h>
#include <MeshUtils/TetGenInterface.h>

using namespace std;
namespace Nektar{
namespace MeshUtils{

void TetMesh::Mesh()
{
    if(m_verbose)
        cout << endl << "Tetrahdral mesh generation" << endl;

    TetGenInterfaceSharedPtr tetgen =
        MemoryManager<TetGenInterface>::AllocateSharedPtr();

    m_surfacemesh->Get(Nodes,Edges,Tris);

    map<int, MeshNodeSharedPtr>::iterator nit;
    for(nit = Nodes.begin(); nit != Nodes.end(); nit++)
    {
        vector<int> t = nit->second->GetTris();
        vector<int> e = nit->second->GetEdges();
        if(t.size() > 0 && e.size() > 1)
        {
            nodesintris.push_back(nit->first);
        }
    }

    tetgen->Assign(nodesintris, Tris, Nodes, m_stienerpoints);

    if(m_verbose)
        cout << "\tMesh iteration : 1" << endl;

    tetgen->Mesh();

    tetgen->Extract(numtet, tetconnect);

    if(m_verbose)
        cout << "\tTets : " << numtet << endl << endl;

    bool repeat = true;
    int meshcounter = 1;

    while(repeat)
    {
        repeat = Validate(Nodes);
        if(!repeat)
        {
            break;
        }

        tetgen->Assign(nodesintris, Tris, Nodes, m_stienerpoints);

        if(m_verbose)
            cout << "\tMesh iteration : " << meshcounter+1 << endl;

        tetgen->Mesh();

        tetgen->Extract(numtet, tetconnect);

        if(m_verbose)
            cout << "\tTets : " << numtet << endl << endl;

        meshcounter++;
    }

    if(m_verbose)
        cout << "\tMeshing iterations: " << meshcounter << endl <<
                "\tTets :" << numtet << endl;

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
}

bool TetMesh::Validate(std::map<int, MeshNodeSharedPtr> &Nodes)
{
    int pointsbefore = m_stienerpoints.size();
    for(int i = 0; i < numtet; i++)
    {
        int tetvert[4];
        tetvert[0] = tetconnect[i][0];
        tetvert[1] = tetconnect[i][1];
        tetvert[2] = tetconnect[i][2];
        tetvert[3] = tetconnect[i][3];

        Array<OneD, NekDouble> tetdelta(4);

        Array<OneD, NekDouble> r(6);

        r[0]=Nodes[tetvert[0]]->Distance(Nodes[tetvert[1]]);
        r[1]=Nodes[tetvert[1]]->Distance(Nodes[tetvert[2]]);
        r[2]=Nodes[tetvert[2]]->Distance(Nodes[tetvert[3]]);
        r[3]=Nodes[tetvert[3]]->Distance(Nodes[tetvert[0]]);
        r[4]=Nodes[tetvert[1]]->Distance(Nodes[tetvert[3]]);
        r[5]=Nodes[tetvert[2]]->Distance(Nodes[tetvert[0]]);

        tetdelta[0] = m_octree->Query(Nodes[tetvert[0]]->GetLoc());
        tetdelta[1] = m_octree->Query(Nodes[tetvert[1]]->GetLoc());
        tetdelta[2] = m_octree->Query(Nodes[tetvert[2]]->GetLoc());
        tetdelta[3] = m_octree->Query(Nodes[tetvert[3]]->GetLoc());

        int numvalid = 0;

        if(r[0] < tetdelta[0])
            numvalid++;

        if(r[1] < tetdelta[1])
            numvalid++;

        if(r[2] < tetdelta[2])
            numvalid++;

        if(r[3] < tetdelta[3])
            numvalid++;

        if(r[4] < tetdelta[1])
            numvalid++;

        if(r[5] < tetdelta[2])
            numvalid++;

        if(numvalid != 6)
        {
            Array<OneD, NekDouble> locn(3);
            locn[0] = 0.0; locn[1] = 0.0; locn[2] = 0.0;
            for(int j = 0; j < 4; j++)
            {
                Array<OneD, NekDouble> loc = Nodes[tetvert[j]]->GetLoc();
                locn[0]+=loc[0]/6.0;
                locn[1]+=loc[1]/6.0;
                locn[2]+=loc[2]/6.0;
            }
            AddNewPoint(locn, Nodes);
        }
    }
    if(m_stienerpoints.size() == pointsbefore)
    {
        return false;
    }
    else
    {
        return true;
    }

}


void TetMesh::AddNewPoint(Array<OneD, NekDouble> loc,
                              std::map<int, MeshNodeSharedPtr> &Nodes)
{
    NekDouble npDelta = m_octree->Query(loc);

    MeshNodeSharedPtr n = boost::shared_ptr<MeshNode>(
                        new MeshNode(Nodes.size(),loc[0],loc[1],loc[2]));

    bool add = true;

    for(int i = 0; i < nodesintris.size(); i++)
    {

        NekDouble r = Nodes[nodesintris[i]]->Distance(n);

        if(r<npDelta/1.414)
        {
            add = false;
            break;
        }

    }

    if(add)
    {
        for(int i = 0; i < m_stienerpoints.size(); i++)
        {
            NekDouble r = Nodes[m_stienerpoints[i]]->Distance(n);

            if(r<npDelta/1.414)
            {
                add = false;
                break;
            }
        }
    }

    if(add)
    {
        Nodes[Nodes.size()] = n;
        m_stienerpoints.push_back(Nodes.size()-1);
    }
}

}
}
