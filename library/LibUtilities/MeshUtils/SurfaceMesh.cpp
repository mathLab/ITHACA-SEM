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

#include <LibUtilities/MeshUtils/SurfaceMesh.h>
#include <LibUtilities/MeshUtils/TriangleInterface.h>

using namespace std;
namespace Nektar{
namespace LibUtilities{
namespace MeshUtils {
    
    void SurfaceMesh::Mesh()
    {
        OrientateCurves();
        
        TriangleInterfaceSharedPtr pplanemesh =
            MemoryManager<TriangleInterface>::AllocateSharedPtr();
        
        pplanemesh->Assign(m_uvloops, m_centers, m_extrapoints);
        
        pplanemesh->Mesh();
        
        Array<OneD, Array<OneD, NekDouble> > Points;
        Array<OneD, Array<OneD, int> > Connec;
        
        pplanemesh->Extract(Points,Connec);
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
            for(int j = 0; i < m_uvloops[i].size(); i++)
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
}

