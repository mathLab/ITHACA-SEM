////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMesh.cpp
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
//  Description: surfacemesh object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <MeshUtils/SurfaceMeshing/SurfaceMesh.h>
#include <MeshUtils/ExtLibInterface/TriangleInterface.h>

#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

void SurfaceMesh::Mesh(std::map<int, MeshNodeSharedPtr> &Nodes,
                       std::map<int, MeshEdgeSharedPtr> &Edges,
                       std::map<int, MeshTriSharedPtr> &Tris)
{
    Stretching();

    OrientateCurves(Nodes);

    //create interface to triangle thirdparty library
    TriangleInterfaceSharedPtr pplanemesh =
        MemoryManager<TriangleInterface>::AllocateSharedPtr();

    pplanemesh->Assign(orderedLoops, m_centers, Nodes,
                       m_stienerpoints, m_id, asr/pasr);

    pplanemesh->Mesh();

    pplanemesh->Extract(nump,numtri, Connec);

    bool repeat = true;
    int meshcounter = 1;

    while (repeat)
    {
        repeat = Validate(Nodes);
        if(!repeat)
        {
            break;
        }
        pplanemesh->Assign(orderedLoops, m_centers, Nodes,
                           m_stienerpoints, m_id, asr/pasr);
        pplanemesh->Mesh();
        pplanemesh->Extract(nump,numtri,Connec);
        meshcounter++;
    }

    //look through all the existing edges and tris in the complete surface mesh
    //if not a duplicate add a new one if duplicate set up data and move on
    for(int i = 0; i < numtri; i++)
    {
        int n1,n2,n3;
        n1 = Connec[i][0];
        n2 = Connec[i][1];
        n3 = Connec[i][2];

        Nodes[n1]->SetTri(Tris.size());
        Nodes[n2]->SetTri(Tris.size());
        Nodes[n3]->SetTri(Tris.size());

        int e1,e2,e3;
        e1 = Nodes[n1]->EdgeInCommon(Nodes[n2]);
        e2 = Nodes[n2]->EdgeInCommon(Nodes[n3]);
        e3 = Nodes[n3]->EdgeInCommon(Nodes[n1]);
        if(e1 == -1)
        {
            MeshEdgeSharedPtr e = MemoryManager<MeshEdge>::AllocateSharedPtr(
                Edges.size(),n1,n2);
            e->SetSurf(m_id);
            Nodes[n1]->SetEdge(Edges.size());
            Nodes[n2]->SetEdge(Edges.size());
            e1 = Edges.size();
            Edges[Edges.size()] = e;
        }
        if(e2 == -1)
        {
            MeshEdgeSharedPtr e = MemoryManager<MeshEdge>::AllocateSharedPtr(
                Edges.size(),n2,n3);
            e->SetSurf(m_id);
            Nodes[n2]->SetEdge(Edges.size());
            Nodes[n3]->SetEdge(Edges.size());
            e2 = Edges.size();
            Edges[Edges.size()] = e;
        }
        if(e3 == -1)
        {
            MeshEdgeSharedPtr e = MemoryManager<MeshEdge>::AllocateSharedPtr(
                Edges.size(),n3,n1);
            e->SetSurf(m_id);
            Nodes[n3]->SetEdge(Edges.size());
            Nodes[n1]->SetEdge(Edges.size());
            e3 = Edges.size();
            Edges[Edges.size()] = e;
        }

        Edges[e1]->SetTri(Tris.size());
        Edges[e2]->SetTri(Tris.size());
        Edges[e3]->SetTri(Tris.size());

        Tris[Tris.size()] = MemoryManager<MeshTri>::AllocateSharedPtr(
            Tris.size(),n1,n2,n3,e1,e2,e3,m_id);
    }
}

void SurfaceMesh::Report()
{
    int edgec = 0;
    for(int i = 0; i < m_cadsurf->GetEdges().size(); i++)
    {
        edgec+=m_cadsurf->GetEdges()[i].size();
    }
    cout << scientific << "\tPoints: " << nump << "\tTris: " << numtri << "\tEdges: " << edgec <<  "\tLoops: " << orderedLoops.size() << endl;
}

void SurfaceMesh::Stretching()
{
    //define a sampling and calculate the aspect ratio of the paramter plane
    asr = 0.0;
    Array<OneD, NekDouble> bnds = m_cadsurf->GetBounds();
    pasr = (bnds[1] - bnds[0])/
           (bnds[3] - bnds[2]);

    Array<TwoD, Array<OneD,NekDouble> > stretch(40,40);

    NekDouble du = (bnds[1]-bnds[0])/(40-1);
    NekDouble dv = (bnds[3]-bnds[2])/(40-1);

    for(int i = 0; i < 40; i++)
    {
        for(int j = 0; j < 40; j++)
        {
            Array<OneD, NekDouble> uv(2);
            uv[0] = bnds[0] + i*du;
            uv[1] = bnds[2] + j*dv;
            stretch[i][j]=m_cadsurf->P(uv);
        }
    }

    int ct = 0;

    for(int i = 0; i < 40-1; i++)
    {
        for(int j = 0; j < 40-1; j++)
        {
            NekDouble ru = sqrt((stretch[i][j][0]-stretch[i+1][j][0])*
                                (stretch[i][j][0]-stretch[i+1][j][0])+
                                (stretch[i][j][1]-stretch[i+1][j][1])*
                                (stretch[i][j][1]-stretch[i+1][j][1])+
                                (stretch[i][j][2]-stretch[i+1][j][2])*
                                (stretch[i][j][2]-stretch[i+1][j][2]));
            NekDouble rv = sqrt((stretch[i][j][0]-stretch[i][j+1][0])*
                                (stretch[i][j][0]-stretch[i][j+1][0])+
                                (stretch[i][j][1]-stretch[i][j+1][1])*
                                (stretch[i][j][1]-stretch[i][j+1][1])+
                                (stretch[i][j][2]-stretch[i][j+1][2])*
                                (stretch[i][j][2]-stretch[i][j+1][2]));

            if(rv < 1E-8)
                continue;

            asr += ru/rv;
            ct++;
        }
    }

    asr/=ct;
}

bool SurfaceMesh::Validate(std::map<int, MeshNodeSharedPtr> &Nodes)
{
    //check all edges in the current mesh for length against the octree
    //if the octree is not conformed to add a new point inside the triangle
    //if no new points are added meshing can stop
    int pointBefore = m_stienerpoints.size();
    for(int i = 0; i < numtri; i++)
    {
        int triVert[3];
        triVert[0]=Connec[i][0];
        triVert[1]=Connec[i][1];
        triVert[2]=Connec[i][2];

        Array<OneD, NekDouble> triDelta(3);

        Array<OneD, NekDouble> r(3);

        r[0]=Nodes[triVert[0]]->Distance(Nodes[triVert[1]]);
        r[1]=Nodes[triVert[1]]->Distance(Nodes[triVert[2]]);
        r[2]=Nodes[triVert[2]]->Distance(Nodes[triVert[0]]);

        triDelta[0] = m_octree->Query(Nodes[triVert[0]]->GetLoc());
        triDelta[1] = m_octree->Query(Nodes[triVert[1]]->GetLoc());
        triDelta[2] = m_octree->Query(Nodes[triVert[2]]->GetLoc());

        int numValid = 0;

        if(r[0] < triDelta[0])
            numValid++;

        if(r[1] < triDelta[1])
            numValid++;

        if(r[2] < triDelta[2])
            numValid++;

        if(numValid != 3)
        {
            Array<OneD,NekDouble> ainfo,binfo,cinfo;
            ainfo = Nodes[triVert[0]]->GetS(m_id);
            binfo = Nodes[triVert[1]]->GetS(m_id);
            cinfo = Nodes[triVert[2]]->GetS(m_id);

            Array<OneD, NekDouble> uvc(2);
            uvc[0] = (ainfo[0]+binfo[0]+cinfo[0])/3.0;
            uvc[1] = (ainfo[1]+binfo[1]+cinfo[1])/3.0;
            AddNewPoint(uvc,Nodes);
        }
    }

    if(m_stienerpoints.size() == pointBefore)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void SurfaceMesh::AddNewPoint(Array<OneD, NekDouble> uv,
                              std::map<int, MeshNodeSharedPtr> &Nodes)
{
    //adds a new point but checks that there are no other points nearby first
    Array<OneD, NekDouble> np = m_cadsurf->P(uv);
    NekDouble npDelta = m_octree->Query(np);

    MeshNodeSharedPtr n = boost::shared_ptr<MeshNode>(
                        new MeshNode(Nodes.size(),np[0],np[1],np[2]));

    bool add = true;

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        for(int j = 0; j < orderedLoops[i].size(); j++)
        {
            NekDouble r = Nodes[orderedLoops[i][j]]->Distance(n);

            if(r<npDelta/1.414)
            {
                add = false;
                break;
            }
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
        n->SetSurf(m_id,uv);
        Nodes[Nodes.size()] = n;
        m_stienerpoints.push_back(Nodes.size()-1);
    }
}

void SurfaceMesh::OrientateCurves(std::map<int, MeshNodeSharedPtr> &Nodes)
{
    //create integer list of bounding loop node id's.
    for(int i = 0; i < m_edges.size(); i++)
    {
        vector<int> cE;
        for(int j = 0; j < m_edges[i].size(); j++)
        {
            vector<int> edgePoints =
                    m_curvemeshes[m_edges[i][j].first]->
                        GetMeshPoints();

            int numPoints = m_curvemeshes[m_edges[i][j].first]->
                                GetNumPoints();

            if(m_edges[i][j].second == 0)
            {
                for(int k = 0; k < numPoints-1; k++)
                {
                    cE.push_back(edgePoints[k]);
                }
            }
            else
            {
                for(int k = numPoints-1; k >0; k--)
                {
                    cE.push_back(edgePoints[k]);
                }
            }
        }
        orderedLoops.push_back(cE);
    }

    //loops made need to orientate on which is biggest and define holes

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        int half = int(orderedLoops[i].size()/2) - 1;

        MeshNodeSharedPtr n1,n2,nh;

        n1 = Nodes[orderedLoops[i][0]];
        n2 = Nodes[orderedLoops[i][1]];
        nh = Nodes[orderedLoops[i][half]];

        Array<OneD,NekDouble> n1info,n2info,nhinfo;
        n1info = n1->GetS(m_id);
        n2info = n2->GetS(m_id);
        nhinfo = nh->GetS(m_id);

        NekDouble ua = (100.0*n1info[0]+
                        100.0*n2info[0]+
                        1.0* nhinfo[0])/201.0 ;
        NekDouble va = (100.0*n1info[1]+
                        100.0*n2info[1]+
                        1.0* nhinfo[1])/201.0 ;

        vector<NekDouble> tmp;
        tmp.push_back(ua);
        tmp.push_back(va);
        m_centers.push_back(tmp);
    }

    vector<NekDouble> areas;

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        NekDouble area=0.0;
        for(int j = 0; j < orderedLoops[i].size()-1; j++)
        {
            MeshNodeSharedPtr n1,n2;
            n1 = Nodes[orderedLoops[i][j]];
            n2 = Nodes[orderedLoops[i][j+1]];
            Array<OneD,NekDouble> n1info,n2info;
            n1info = n1->GetS(m_id);
            n2info = n2->GetS(m_id);

            area+=-n2info[1]*
                        (n2info[0]-n1info[0])
                  +n1info[0]*
                        (n2info[1]-n1info[1]);
        }
        area*=0.5;
        areas.push_back(area);
    }

    int ct=0;

    do
    {
        ct=0;
        for(int i = 0; i < areas.size()-1; i++)
        {
            if(abs(areas[i])<abs(areas[i+1]))
            {
                //swap
                NekDouble areatmp = areas[i];
                vector<NekDouble> centerstmp = m_centers[i];
                vector<int> orderedlooptmp = orderedLoops[i];
                vector<pair<int,int> > edgeLoopstmp = m_edges[i];

                areas[i]=areas[i+1];
                m_centers[i]=m_centers[i+1];
                orderedLoops[i]=orderedLoops[i+1];
                m_edges[i]=m_edges[i+1];

                areas[i+1]=areatmp;
                m_centers[i+1]=centerstmp;
                orderedLoops[i+1]=orderedlooptmp;
                m_edges[i+1]=edgeLoopstmp;

                ct+=1;
            }
        }

    }while(ct>0);

    if(areas[0]<0) //reverse the first uvLoop
    {
        vector<int > tmp = orderedLoops[0];
        reverse(tmp.begin(), tmp.end());
        orderedLoops[0]=tmp;
    }

    for(int i = 1; i < orderedLoops.size(); i++)
    {
        if(areas[i]>0) //reverse the loop
        {
            vector<int> tmp = orderedLoops[i];
            reverse(tmp.begin(), tmp.end());
            orderedLoops[i]=tmp;
        }
    }


}

}
}
