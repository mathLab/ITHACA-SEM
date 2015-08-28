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

    //first level of high-order awarness, placing stienner points in specific
    //locations above triangles, in this iteration points are place by linear
    //interpoliation of the center of the triangle, taking the normal of the linear triangle
    map<int,int> nototri; //map the new node id to the tiangle it was created from

    map<int, MeshTriSharedPtr>::iterator trit;
    for(trit = Tris.begin(); trit != Tris.end(); trit++)
    {
        Array<OneD, int> node = trit->second->GetN();
        Array<OneD, NekDouble> l1,l2,l3;
        l1 = Nodes[node[0]]->GetLoc();
        l2 = Nodes[node[1]]->GetLoc();
        l3 = Nodes[node[2]]->GetLoc();

        Array<OneD, NekDouble> v1(3), v2(3);
        v1[0] = l2[0] - l1[0]; v1[1] = l2[1] - l1[1]; v1[2] = l2[2] - l1[2];
        v2[0] = l3[0] - l1[0]; v2[1] = l3[1] - l1[1]; v2[2] = l3[2] - l1[2];
        Array<OneD, NekDouble> N(3);
        N[0] = v1[1]*v2[2] - v1[2]*v2[1];
        N[1] = v1[0]*v2[2] - v1[2]*v2[0];
        N[2] = v1[0]*v2[1] - v1[1]*v2[0];
        Array<OneD, NekDouble> loc(3);
        loc[0] = (l1[0]+l2[0]+l3[0])/3.0;
        loc[1] = (l1[1]+l2[1]+l3[1])/3.0;
        loc[2] = (l1[2]+l2[2]+l3[2])/3.0;
        NekDouble d = m_octree->Query(loc);

        NekDouble a,b,c;
        a = N[0]*N[0] +  N[1]*N[1] + N[2]*N[2];
        b = 2.0*(N[0]*(loc[0]-l1[0]) + N[1]*(loc[1]-l1[1]) + N[2]*(loc[2]-l1[2]));
        c = (loc[0]-l1[0])*(loc[0]-l1[0]) + (loc[1]-l1[1])*(loc[1]-l1[1]) + (loc[2]-l1[2])*(loc[2]-l1[2]) + d*d;
        NekDouble t1,t2;
        cout << b*b-4.0*a*c << endl;
        t1 = (-b+sqrt(b*b-4.0*a*c))/2.0/a; t2 = (-b-sqrt(b*b-4.0*a*c))/2.0/a;
        if(t1 > 0)
        {
            loc[0] = loc[0] + N[0]*t1;
            loc[1] = loc[1] + N[1]*t1;
            loc[2] = loc[2] + N[2]*t1;
        }
        else
        {
            loc[0] = loc[0] + N[0]*t2;
            loc[1] = loc[1] + N[1]*t2;
            loc[2] = loc[2] + N[2]*t2;
        }
        cout << loc[0] << " " << loc[1] << " " << endl;

        //should look over the neigbouring tris through edge links,
        //if the associated stiener point is too close, they should be merged.. not ignored
        MeshNodeSharedPtr n = boost::shared_ptr<MeshNode>(
                            new MeshNode(Nodes.size(),loc[0],loc[1],loc[2]));
        d = m_octree->Query(loc);
        bool add = true;
        for(int i = 0; i < m_stienerpoints.size(); i++)
        {
            NekDouble dist = Nodes[m_stienerpoints[i]]->Distance(n);
            if(dist < d/1.141)
            {
                add = false;
            }
        }
        if(true)
        {
            m_stienerpoints.push_back(Nodes.size());
            Nodes[Nodes.size()] = n;
        }

    }

    tetgen->Assign(nodesintris, Tris, Nodes, m_stienerpoints);

    tetgen->Mesh();

    tetgen->Extract(numtet, tetconnect);

    if(m_verbose)
        cout << "\tMesh iteration: 1 \tTets : " << numtet << endl;

    bool repeat = true;
    int meshcounter = 1;

    /*while(repeat)
    {
        repeat = Validate(Nodes);
        if(!repeat)
        {
            break;
        }

        tetgen->Assign(nodesintris, Tris, Nodes, m_stienerpoints);

        tetgen->Mesh();

        tetgen->Extract(numtet, tetconnect);

        if(m_verbose)
            cout << "\tMesh iteration : " << meshcounter+1 << "\tTets : " << numtet << endl;

        meshcounter++;
    }*/

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
        cout << endl << "\tMeshing iterations: " << meshcounter << endl <<
                "\tTets :" << numtet << endl <<
                "\tBoundary nodes: " << nodesintris.size() << endl <<
                "\tInterior nodes: " << m_stienerpoints.size() << endl;
}

bool TetMesh::Validate(std::map<int, MeshNodeSharedPtr> &Nodes)
{
    nodesaddedinthisvalidate.clear();
    int pointsbefore = m_stienerpoints.size();
    for(int i = 0; i < numtet; i++)
    {
        vector<int> tetvert;
        tetvert.resize(4);
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
                locn[0]+=loc[0]/4.0;
                locn[1]+=loc[1]/4.0;
                locn[2]+=loc[2]/4.0;
            }
            AddNewPoint(locn, Nodes, tetvert);
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
                              map<int, MeshNodeSharedPtr> &Nodes,
                              vector<int> tetnodes)
{
    NekDouble npDelta = m_octree->Query(loc);

    MeshNodeSharedPtr n = boost::shared_ptr<MeshNode>(
                        new MeshNode(Nodes.size(),loc[0],loc[1],loc[2]));

    bool add = true;

    for(int i = 0; i < tetnodes.size(); i++)
    {

        NekDouble r = Nodes[tetnodes[i]]->Distance(n);

        if(r<npDelta/1.414)
        {
            add = false;
            break;
        }

    }

    if(add)
    {
        for(int i = 0; i < nodesaddedinthisvalidate.size(); i++)
        {
            NekDouble r = Nodes[nodesaddedinthisvalidate[i]]->Distance(n);

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
        nodesaddedinthisvalidate.push_back(Nodes.size()-1);
        m_stienerpoints.push_back(Nodes.size()-1);
    }
}

}
}
