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

using namespace std;
namespace Nektar{
namespace MeshUtils {
    
    void SurfaceMesh::Mesh()
    {
        OrientateCurves();

        TriangleInterfaceSharedPtr pplanemesh =
            MemoryManager<TriangleInterface>::AllocateSharedPtr();
        
        pplanemesh->Assign(m_uvloops, m_centers, m_extrapoints);
        
        pplanemesh->Mesh();
        
        pplanemesh->Extract(numpoints, numtris, Points,Connec);
        
        bool repeat = true;
        
        while (repeat)
        {
            repeat = Validate(numpoints,numtris,Points,Connec);
            if(!repeat)
            {
                break;
            }
            pplanemesh->Assign(m_uvloops, m_centers, m_extrapoints);
            pplanemesh->Mesh();
            pplanemesh->Extract(numpoints,numtris,Points,Connec);
        }
        
        pplanemesh->Assign(m_uvloops, m_centers, m_extrapoints);
        pplanemesh->Mesh(false,true);
        pplanemesh->Extract(numpoints,numtris,Points,Connec);
        
        /*for(int i = 0; i < numpoints; i++)
        {
            Array<OneD, NekDouble> p = m_cadsurf->P(Points[i][0],Points[i][1]);
            Points[i]=p;
        }*/

    }
    
    void SurfaceMesh::HOMesh(int order)
    {
        for(int i = 0; i < numtris; i++)
        {
            DNekMat a (3,3,1.0);
            a(0,0) = Points[Connec[i][0]][0];
            a(1,0) = Points[Connec[i][0]][1];
            a(0,1) = Points[Connec[i][1]][0];
            a(1,1) = Points[Connec[i][1]][1];
            a(0,2) = Points[Connec[i][2]][0];
            a(1,2) = Points[Connec[i][2]][1];
            
            DNekMat b (3,3,1.0);
            
            DNekMat c (3,3,1.0);
            c(0,0) = -1.0;
            c(1,0) = -1.0;
            c(0,1) = 1.0;
            c(1,1) = -1.0;
            c(0,2) = -1.0;
            c(1,2) = 1.0;
            
            c.Invert();
            
            b = a*c;
            
            DNekMat C (3,3,1.0);
            c(0,0) = -1.0;
            c(1,0) = -1.0;
            c(0,1) = 1.0;
            c(1,1) = -1.0;
            c(0,2) = -1.0;
            c(1,2) = 1.0;
            
            a= b*C;
            
            cout << a(0,0) << " " << Points[Connec[i][0]][0] << endl;
        }
    }
    
    bool SurfaceMesh::Validate(int &np,
                               int &nt,
                               Array<OneD, Array<OneD, NekDouble> > &Points,
                               Array<OneD, Array<OneD, int> > &Connec)
    {
        int pointBefore = m_extrapoints.size();
        for(int i = 0; i < nt; i++)
        {
            int triVert[3];
            triVert[0]=Connec[i][0];
            triVert[1]=Connec[i][1];
            triVert[2]=Connec[i][2];
            
            NekDouble triUV[3][2];
            triUV[0][0]=Points[triVert[0]][0];
            triUV[0][1]=Points[triVert[0]][1];
            triUV[1][0]=Points[triVert[1]][0];
            triUV[1][1]=Points[triVert[1]][1];
            triUV[2][0]=Points[triVert[2]][0];
            triUV[2][1]=Points[triVert[2]][1];
            
            Array<OneD, Array<OneD, NekDouble> > locs(3);
            Array<OneD, NekDouble> triDelta(3);
            for(int i = 0; i < 3; i++)
            {
                locs[i] = m_cadsurf->P(triUV[i][0],triUV[i][1]);
                triDelta[i] = m_octree->Query(locs[i]);
            }
            
            Array<OneD, NekDouble> r(3);
            
            r[0]=sqrt((locs[0][0]-locs[1][0])*(locs[0][0]-locs[1][0]) +
                      (locs[0][1]-locs[1][1])*(locs[0][1]-locs[1][1]) +
                      (locs[0][2]-locs[1][2])*(locs[0][2]-locs[1][2]));
            r[1]=sqrt((locs[1][0]-locs[2][0])*(locs[1][0]-locs[2][0]) +
                      (locs[1][1]-locs[2][1])*(locs[1][1]-locs[2][1]) +
                      (locs[1][2]-locs[2][2])*(locs[1][2]-locs[2][2]));
            r[2]=sqrt((locs[2][0]-locs[0][0])*(locs[2][0]-locs[0][0]) +
                      (locs[2][1]-locs[0][1])*(locs[2][1]-locs[0][1]) +
                      (locs[2][2]-locs[0][2])*(locs[2][2]-locs[0][2]));
            
            int numValid = 0;
            
            if(r[0] < triDelta[0])
                numValid++;
            
            if(r[1] < triDelta[1])
                numValid++;
            
            if(r[2] < triDelta[2])
                numValid++;
            
            if(numValid != 3)
            {
                NekDouble uc = (triUV[0][0]+triUV[1][0]+triUV[2][0])/3.0;
                NekDouble vc = (triUV[0][1]+triUV[1][1]+triUV[2][1])/3.0;
                AddNewPoint(uc,vc);
            }
        }
        
        if(m_extrapoints.size() == pointBefore)
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
        Array<OneD, NekDouble> newPoint = m_cadsurf->P(u,v);
        NekDouble newPointDelta = m_octree->Query(newPoint);
        
        bool add = true;
        
        for(int i = 0; i < m_uvloops.size(); i++)
        {
            for(int j = 0; j < m_uvloops[i].size(); j++)
            {
                Array<OneD, NekDouble> testPoint =
                        m_cadsurf->P(m_uvloops[i][j][0],m_uvloops[i][j][1]);
                
                NekDouble r =
                sqrt((testPoint[0]-newPoint[0])*(testPoint[0]-newPoint[0]) +
                     (testPoint[1]-newPoint[1])*(testPoint[1]-newPoint[1]) +
                     (testPoint[2]-newPoint[2])*(testPoint[2]-newPoint[2]));
                
                if(r<newPointDelta/1.414)
                {
                    add = false;
                    break;
                }
            }
            if(add==false)
            {
                break;
            }
        }
        if(add==true)
        {
            for(int i = 0; i < m_extrapoints.size(); i++)
            {
                Array<OneD, NekDouble> testPoint =
                    m_cadsurf->P(m_extrapoints[i][0],m_extrapoints[i][1]);
                
                NekDouble r =
                sqrt((testPoint[0]-newPoint[0])*(testPoint[0]-newPoint[0]) +
                     (testPoint[1]-newPoint[1])*(testPoint[1]-newPoint[1]) +
                     (testPoint[2]-newPoint[2])*(testPoint[2]-newPoint[2]));
                
                if(r<newPointDelta/1.414)
                {
                    add = false;
                    break;
                }
            }
        }
        
        if(add)
        {
            vector<NekDouble> uv;
            uv.push_back(u);
            uv.push_back(v);
            m_extrapoints.push_back(uv);
        }
    }
    
    void SurfaceMesh::OrientateCurves()
    {
        int edgeCounter = m_numedges;
        
        while(edgeCounter>0)
        {
            int currentEdge = firstEdgeNotUsed();
            bool currentEdgeForward = true;
            bool notClosed = true;
            
            vector<int> currentLoop;
            
            vector<NekDouble> origin = m_curvemeshes[currentEdge-1]->
                                                GetFirstPoint();
            currentLoop.push_back(currentEdge);
            edgeCounter--;
            
            while(notClosed)
            {
                vector<NekDouble> loc;
                if(currentEdgeForward)
                {
                    loc = m_curvemeshes[currentEdge-1]->GetLastPoint();
                }
                else
                {
                    loc = m_curvemeshes[currentEdge-1]->GetFirstPoint();
                }
                
                if(abs(loc[0]-origin[0]) < 1E-6 &&
                   abs(loc[1]-origin[1]) < 1E-6 &&
                   abs(loc[2]-origin[2]) < 1E-6 )
                {
                    notClosed = false;
                    break;
                }
                
                vector<NekDouble> test;
                for(int i = 0; i < m_numedges; i++)
                {
                    if(m_edges[i]==currentEdge)
                        continue;
                    
                    test = m_curvemeshes[m_edges[i]-1]->GetFirstPoint();
                    if(abs(loc[0]-test[0]) < 1E-6 &&
                       abs(loc[1]-test[1]) < 1E-6 &&
                       abs(loc[2]-test[2]) < 1E-6 )
                    {
                        currentLoop.push_back(m_edges[i]);
                        edgeCounter--;
                        currentEdge=m_edges[i];
                        currentEdgeForward=true;
                        break;
                    }
                    else
                    {
                        test = m_curvemeshes[m_edges[i]-1]->GetLastPoint();
                        if(abs(loc[0]-test[0]) < 1E-6 &&
                           abs(loc[1]-test[1]) < 1E-6 &&
                           abs(loc[2]-test[2]) < 1E-6 )
                        {
                            currentLoop.push_back(-1*m_edges[i]);
                            edgeCounter--;
                            currentEdge=m_edges[i];
                            currentEdgeForward=false;
                            break;
                        }
                    }
                }
                
            }
            
            m_edgeloops.push_back(currentLoop);
        }
        
        //create list of continuous loop locs
        
        vector<vector<vector<NekDouble> > > orderedLoops;
        
        for(int i = 0; i < m_edgeloops.size(); i++)
        {
            vector<vector<NekDouble> > cE;
            for(int j = 0; j < m_edgeloops[i].size(); j++)
            {
                vector<vector<NekDouble> > edgePoints =
                        m_curvemeshes[abs(m_edgeloops[i][j])-1]->
                            GetMeshPoints();
                int numPoints = m_curvemeshes[abs(m_edgeloops[i][j])-1]->
                                    GetNumPoints();
                if(m_edgeloops[i][j]>0)
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
        
        for(int i = 0; i < orderedLoops.size(); i++)
        {
            vector<vector<NekDouble> > lP;
            for(int j = 0; j < orderedLoops[i].size(); j++)
            {
                vector<NekDouble> P;
                P.resize(2);
                m_cadsurf->locuv(P[0],P[1],orderedLoops[i][j]);
                lP.push_back(P);
            }
            m_uvloops.push_back(lP);
        }
        
        //loops made need to orientate on which is biggest and define holes
        
        for(int i = 0; i < m_uvloops.size(); i++)
        {
            NekDouble ua = 0.0;
            NekDouble va = 0.0;
            for(int j = 0; j < m_uvloops[i].size(); j++)
            {
                ua+=m_uvloops[i][j][0];
                va+=m_uvloops[i][j][1];
            }
            ua/=m_uvloops[i].size();
            va/=m_uvloops[i].size();
            vector<NekDouble> tmp;
            tmp.push_back(ua);
            tmp.push_back(va);
            m_centers.push_back(tmp);
        }
        
        vector<NekDouble> areas;
        
        for(int i = 0; i < m_uvloops.size(); i++)
        {
            NekDouble area=0.0;
            for(int j = 0; j < m_uvloops[i].size()-1; j++)
            {
                area+=-m_uvloops[i][j][1]*
                            (m_uvloops[i][j+1][0]-m_uvloops[i][j][0])
                      +m_uvloops[i][j][0]*
                            (m_uvloops[i][j+1][1]-m_uvloops[i][j][1]);
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
                    vector<vector<NekDouble> > uvLoopstmp = m_uvloops[i];
                    vector<int> edgeLoopstmp = m_edgeloops[i];
                    
                    areas[i]=areas[i+1];
                    m_centers[i]=m_centers[i+1];
                    m_uvloops[i]=m_uvloops[i+1];
                    m_edgeloops[i]=m_edgeloops[i+1];
                    
                    areas[i+1]=areatmp;
                    m_centers[i+1]=centerstmp;
                    m_uvloops[i+1]=uvLoopstmp;
                    m_edgeloops[i+1]=edgeLoopstmp;
                    
                    ct+=1;
                }
            }
            
        }while(ct>0);
        
        if(areas[0]<0) //reverse the first uvLoop
        {
            vector<vector<NekDouble> > tmp = m_uvloops[0];
            reverse(tmp.begin(), tmp.end());
            m_uvloops[0]=tmp;
        }
        
        for(int i = 1; i < m_uvloops.size(); i++)
        {
            if(areas[i]>0) //reverse the loop
            {
                vector<vector<NekDouble> > tmp = m_uvloops[i];
                reverse(tmp.begin(), tmp.end());
                m_uvloops[i]=tmp;
            }
        }
        
        
    }
    
    
    int SurfaceMesh::firstEdgeNotUsed()
    {
        for(int i = 0; i < m_numedges; i++)
        {
            bool found = false;
            for(int j = 0; j < m_edgeloops.size(); j++)
            {
                for(int k = 0; k < m_edgeloops[j].size(); k++)
                {
                    if(m_edges[i]==abs(m_edgeloops[j][k]))
                    {
                        found = true;
                    }
                }
            }
            if(found==false)
            {
                return m_edges[i];
                break;
            }
        }
        return -1;
    }
    
}
}

