////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
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

#include <string>
#include <fstream>

#include <MeshUtils/SurfaceMesh.h>
#include <MeshUtils/TriangleInterface.h>

#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {
    
    void SurfaceMesh::Mesh()
    {
        Stretching();
        
        nodesused = 0;
        
        OrientateCurves();

        TriangleInterfaceSharedPtr pplanemesh =
            MemoryManager<TriangleInterface>::AllocateSharedPtr();
        
        pplanemesh->Assign(orderedLoops, m_centers, Nodes,asr/pasr);
        
        pplanemesh->Mesh();
        
        pplanemesh->Extract(numtri, Connec);
        
        bool repeat = true;
        
        while (repeat)
        {
            repeat = Validate();
            if(!repeat)
            {
                break;
            }
            pplanemesh->Assign(orderedLoops, m_centers, Nodes,asr/pasr);
            pplanemesh->Mesh();
            pplanemesh->Extract(numtri,Connec);
        }
        
        int numedges;
        Array<OneD, Array<OneD, int> > edgelist;
        pplanemesh->GetEdges(edgelist,numedges);
        
        for(int i = 0; i < numedges; i++)
        {
            MeshEdgeSharedPtr e = boost::shared_ptr<MeshEdge>(
                                        new MeshEdge(i,
                                                     Nodes[edgelist[i][0]],
                                                     Nodes[edgelist[i][0]]));
            Edges.push_back(e);
        }
        
        //LinearOptimise();
        
        //HOMesh();

    }
    
    void SurfaceMesh::LinearOptimise()
    {
        /*for(int i = 0; i < 1; i++)
        {
            EdgeSwap();
        }*/
    }
    
    void SurfaceMesh::EdgeSwap()
    {
        /*for(int i = 0; i < numtris; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                if(neigbourlist[i][j] != -1) //has a edge to swap
                {
                    int A,B,C,D;
                    
                    if(j==0)
                    {
                        A = Connec[i][0];
                        B = Connec[i][1];
                        C = Connec[i][2];
                    }
                    else if(j==1)
                    {
                        A = Connec[i][1];
                        B = Connec[i][2];
                        C = Connec[i][0];
                    }
                    else if(j==2)
                    {
                        A = Connec[i][2];
                        B = Connec[i][0];
                        C = Connec[i][1];
                        }
                    
                    
                    
                }
            }
        }*/
    }
    
    void SurfaceMesh::HOMesh()
    {
        /*LibUtilities::PointsKey pkey(m_order+1,
                                     LibUtilities::eNodalTriEvenlySpaced);
        Array<OneD, NekDouble> u,v;
        
        TotNumPoints = LibUtilities::PointsManager()[pkey]->
                                                        GetTotNumPoints();
        LibUtilities::PointsManager()[pkey]->GetPoints(u,v);
        
        HOPoints = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(numtris);
        
        
        for(int i = 0; i < numtris; i++)
        {
            DNekMat a (3,3,1.0);
            a(0,0) = Points[Connec[i][0]][0];
            a(1,0) = Points[Connec[i][0]][1];
            a(2,0) = 1.0;
            a(0,1) = Points[Connec[i][1]][0];
            a(1,1) = Points[Connec[i][1]][1];
            a(2,1) = 1.0;
            a(0,2) = Points[Connec[i][2]][0];
            a(1,2) = Points[Connec[i][2]][1];
            a(2,2) = 1.0;
            
            DNekMat c (3,3,1.0);
            c(0,0) = u[0];
            c(1,0) = v[0];
            c(2,0) = 1.0;
            c(0,1) = u[1];
            c(1,1) = v[1];
            c(2,1) = 1.0;
            c(0,2) = u[2];
            c(1,2) = v[2];
            c(2,2) = 1.0;
            
            c.Invert();
            
            DNekMat M = a*c;
            
            DNekMat p (3,TotNumPoints,1.0);
            
            for(int j = 0; j < TotNumPoints; j++)
            {
                p(0,j) = u[j];
                p(1,j) = v[j];
            }
            
            DNekMat result = M*p;
            
            HOPoints[i] = Array<OneD, Array<OneD, NekDouble> >(TotNumPoints);
            
            for(int j = 0; j < TotNumPoints; j++)
            {
                Array<OneD, NekDouble> P(2);
                P[0] = result(0,j);
                P[1] = result(1,j);
                HOPoints[i][j]=P;
            }
        }*/
    }
    
    void SurfaceMesh::Stretching()
    {
        asr = 0.0;
        pasr = (m_cadsurf->maxU() - m_cadsurf->minU())/
               (m_cadsurf->maxV() - m_cadsurf->minV());
        
        Array<TwoD, Array<OneD,NekDouble> > stretch(40,40);
        
        NekDouble du = (m_cadsurf->maxU()-m_cadsurf->minU())/(40-1);
        NekDouble dv = (m_cadsurf->maxV()-m_cadsurf->minV())/(40-1);
        
        for(int i = 0; i < 40; i++)
        {
            for(int j = 0; j < 40; j++)
            {
                stretch[i][j]=m_cadsurf->P(m_cadsurf->minU() + i*du,
                                           m_cadsurf->minV() + j*dv);
            }
        }
        
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
                asr += ru/rv;
            }
        }
        
        asr/=(40-1)*(40-1);
        
    }
    
    bool SurfaceMesh::Validate()
    {
        int pointBefore = Nodes.size();
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
                vector<pair<int, Array<OneD,NekDouble> > > ainfo,binfo,cinfo;
                ainfo = Nodes[triVert[0]]->GetS();
                binfo = Nodes[triVert[1]]->GetS();
                cinfo = Nodes[triVert[2]]->GetS();
                
                NekDouble uc = (ainfo[0].second[0]+
                                binfo[0].second[0]+
                                cinfo[0].second[0])/3.0;
                NekDouble vc = (ainfo[0].second[1]+
                                binfo[0].second[1]+
                                cinfo[0].second[1])/3.0;
                AddNewPoint(uc,vc);
            }
        }
        
        if(Nodes.size() == pointBefore)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    
    void SurfaceMesh::AddNewPoint(NekDouble u, NekDouble v)
    {
        Array<OneD, NekDouble> np = m_cadsurf->P(u,v);
        NekDouble npDelta = m_octree->Query(np);
        
        MeshNodeSharedPtr n = boost::shared_ptr<MeshNode>(
                            new MeshNode(nodesused,np[0],np[1],np[2]));
        
        bool add = true;
        
        for(int i = 0; i < Nodes.size(); i++)
        {
            NekDouble r = Nodes[i]->Distance(n);
            
            if(r<npDelta/1.414)
            {
                add = false;
                break;
            }
        }
        
        if(add)
        {
            n->SetSurf(m_id,u,v);
            Nodes.push_back(n);
            nodesused++;
        }
    }
    
    void SurfaceMesh::OrientateCurves()
    {
        //assign nodes with common curves. done
        //fetch nodes from curves. done
        //create integer list of bounding loop node id's. done
        //locuv nodes to get uv. done
        
        MeshNodeSharedPtr ce,ne;
        vector<pair<int, NekDouble> > ceinfo,neinfo;
        for(int i = 0; i < m_edges.size(); i++)
        {
            for(int j = 0; j < m_edges[i].size()-1; j++)
            {
                if(m_edges[i][j].second == 0)
                {
                    ce=m_curvemeshes[m_edges[i][j].first-1]->
                    GetLastPoint();
                    if(m_edges[i][j+1].second == 0)
                    {
                        ne=m_curvemeshes[m_edges[i][j+1].first-1]
                                         ->GetFirstPoint();
                    }
                    else
                    {
                        ne=m_curvemeshes[m_edges[i][j+1].first-1]
                                         ->GetLastPoint();
                    }
                }
                else
                {
                    ce=m_curvemeshes[m_edges[i][j].first-1]->
                    GetFirstPoint();
                    if(m_edges[i][j+1].second == 0)
                    {
                        ne=m_curvemeshes[m_edges[i][j+1].first-1]
                        ->GetFirstPoint();
                    }
                    else
                    {
                        ne=m_curvemeshes[m_edges[i][j+1].first-1]
                        ->GetLastPoint();
                    }
                }
                
                ceinfo = ce->GetC();
                neinfo = ne->GetC();
                
                ce->SetCurve(neinfo[0].first,neinfo[0].second);
                ne->SetCurve(ceinfo[0].first,ceinfo[0].second);
                //at this point know they only have 1 entry
            }
            
            if(m_edges[i].back().second == 0)
            {
                ce=m_curvemeshes[m_edges[i].back().first-1]->
                GetLastPoint();
                if(m_edges[i][0].second == 0)
                {
                    ne=m_curvemeshes[m_edges[i][0].first-1]
                    ->GetFirstPoint();
                }
                else
                {
                    ne=m_curvemeshes[m_edges[i][0].first-1]
                    ->GetLastPoint();
                }
            }
            else
            {
                ce=m_curvemeshes[m_edges[i].back().first-1]->
                GetFirstPoint();
                if(m_edges[i][0].second == 0)
                {
                    ne=m_curvemeshes[m_edges[i][0].first-1]
                    ->GetFirstPoint();
                }
                else
                {
                    ne=m_curvemeshes[m_edges[i][0].first-1]
                    ->GetLastPoint();
                }
            }
            ceinfo = ce->GetC();
            neinfo = ne->GetC();
            
            ce->SetCurve(neinfo[0].first,neinfo[0].second);
            ne->SetCurve(ceinfo[0].first,ceinfo[0].second);
        }
        
        int nodecounter=0;
        
        for(int i = 0; i < m_edges.size(); i++)
        {
            vector<int> cE;
            for(int j = 0; j < m_edges[i].size(); j++)
            {
                vector<MeshNodeSharedPtr> edgePoints =
                        m_curvemeshes[m_edges[i][j].first-1]->
                            GetMeshPoints();
                
                int numPoints = m_curvemeshes[m_edges[i][j].first-1]->
                                    GetNumPoints();
                
                if(m_edges[i][j].second == 0)
                {
                    for(int k = 0; k < numPoints-1; k++)
                    {
                        edgePoints[k]->SetID(nodesused);
                        nodesused++;
                        Nodes.push_back(edgePoints[k]);
                        cE.push_back(nodecounter);
                        nodecounter++;
                        
                    }
                }
                else
                {
                    for(int k = numPoints-1; k >0; k--)
                    {
                        edgePoints[k]->SetID(nodesused);
                        nodesused++;
                        Nodes.push_back(edgePoints[k]);
                        cE.push_back(nodecounter);
                        nodecounter++;
                    }
                }
            }
            orderedLoops.push_back(cE);
        }
        
        for(int i = 0; i < Nodes.size(); i++)
        {
            vector<NekDouble> P;
            P.resize(2);
            m_cadsurf->locuv(P[0],P[1],Nodes[i]->GetLoc());
            Nodes[i]->SetSurf(m_id,P[0],P[1]);
        }
        
        //loops made need to orientate on which is biggest and define holes
        
        for(int i = 0; i < orderedLoops.size(); i++)
        {
            int half = int(orderedLoops[i].size()/2) - 1;
            
            MeshNodeSharedPtr n1,n2,nh;
            
            n1 = Nodes[orderedLoops[i][0]];
            n2 = Nodes[orderedLoops[i][1]];
            nh = Nodes[orderedLoops[i][half]];
            
            vector<pair<int, Array<OneD,NekDouble> > > n1info,n2info,nhinfo;
            n1info = n1->GetS();
            n2info = n2->GetS();
            nhinfo = nh->GetS();
            
            NekDouble ua = (100.0*n1info[0].second[0]+
                            100.0*n2info[0].second[0]+
                            1.0* nhinfo[0].second[0])/201.0 ;
            NekDouble va = (100.0*n1info[0].second[1]+
                            100.0*n2info[0].second[1]+
                            1.0* nhinfo[0].second[1])/201.0 ;
            
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
                vector<pair<int, Array<OneD,NekDouble> > > n1info,n2info;
                n1info = n1->GetS();
                n2info = n2->GetS();
                
                area+=-n2info[0].second[1]*
                            (n2info[0].second[0]-n1info[0].second[0])
                      +n1info[0].second[0]*
                            (n2info[0].second[1]-n1info[0].second[1]);
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

