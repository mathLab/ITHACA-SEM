////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.cpp
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

#include <fstream>
#include <string>
#include <sstream>

#include <boost/filesystem.hpp>

#include <LibUtilities/CADSystem/CADSystem.h>

using namespace std;
namespace Nektar{
namespace LibUtilities{
    
    void CADSystem::GetParameterPlaneBounds(int i,
                                            Array<OneD, NekDouble>& out)
    {
        out = Array<OneD, NekDouble>(4);
        out[0]=m_surfs[i-1]->minU();
        out[1]=m_surfs[i-1]->maxU();
        out[2]=m_surfs[i-1]->minV();
        out[3]=m_surfs[i-1]->maxV();
    }

    string CADSystem::GetName()
    {
        return m_name;
    }
    
    void CADSystem::N(int i, NekDouble u, NekDouble v,
                      Array<OneD, NekDouble>& out)
    {
        ASSERTL0(u>=m_surfs[i-1]->minU() &&
                 u<=m_surfs[i-1]->maxU() &&
                 v>=m_surfs[i-1]->minV() &&
                 v<=m_surfs[i-1]->maxV(), "(u,v) out of bounds");
        out = m_surfs[i-1]->N(u,v);
    }
    void CADSystem::D1(int i, NekDouble u, NekDouble v,
                                 Array<OneD, NekDouble>& out)
    {
        ASSERTL0(u>=m_surfs[i-1]->minU() &&
                 u<=m_surfs[i-1]->maxU() &&
                 v>=m_surfs[i-1]->minV() &&
                 v<=m_surfs[i-1]->maxV(), "(u,v) out of bounds");
        out = m_surfs[i-1]->D1(u,v);
    }
    void CADSystem::D2(int i, NekDouble u, NekDouble v,
                                 Array<OneD, NekDouble>& out)
    {
        ASSERTL0(u>=m_surfs[i-1]->minU() &&
                 u<=m_surfs[i-1]->maxU() &&
                 v>=m_surfs[i-1]->minV() &&
                 v<=m_surfs[i-1]->maxV(), "(u,v) out of bounds");
        out = m_surfs[i-1]->D2(u,v);
    }
    
    void CADSystem::Report()
    {
        cout << "CAD has: " << m_numCurve << " curves." << endl;
        cout << "CAD has: " << m_numSurf << " surfaces." << endl;
    }
    
    void CADSystem::GetBoundingBox(Array<OneD, NekDouble>& out)
    {
        out[0]=1000000.0;
        out[1]=0.0;
        out[2]=1000000.0;
        out[3]=0.0;
        out[4]=1000000.0;
        out[5]=0.0;
        
        for(int i = 0; i < m_curves.size(); i++)
        {
            gp_Pnt start, end;
            m_curves[i]->GetMinMax(start,end);
            if(start.X()<out[0])
                out[0]=start.X();
            if(start.X()>out[1])
                out[1]=start.X();
            if(start.Y()<out[2])
                out[2]=start.Y();
            if(start.Y()>out[3])
                out[3]=start.Y();
            if(start.Z()<out[4])
                out[4]=start.Z();
            if(start.Z()>out[5])
                out[5]=start.Z();
            
            if(end.X()<out[0])
                out[0]=end.X();
            if(end.X()>out[1])
                out[1]=end.X();
            if(end.Y()<out[2])
                out[2]=end.Y();
            if(end.Y()>out[3])
                out[3]=end.Y();
            if(end.Z()<out[4])
                out[4]=end.Z();
            if(end.Z()>out[5])
                out[5]=end.Z();
        }
    }

    bool CADSystem::LoadCAD()
    {
        if ( !boost::filesystem::exists( m_name.c_str() ) )
        {
            return false;
        }
    
        string ext;
        size_t pos = m_name.find(".");
        ext = m_name.substr(pos);

        TopoDS_Shape shape;
        
        if(ext.compare(".STEP") == 0 ||
           ext.compare(".step") == 0 ||
           ext.compare(".stp")  == 0 ||
           ext.compare(".STP")  == 0 )
        {
            STEPControl_Reader reader;
            reader = STEPControl_Reader();
            reader.ReadFile(m_name.c_str());
            reader.NbRootsForTransfer();
            reader.TransferRoots();
            shape = reader.OneShape();
            if(shape.IsNull())
            {
                return false;
            }
        }
        else if(ext.compare(".IGES") == 0 ||
                ext.compare(".iges") == 0 ||
                ext.compare(".igs")  == 0 ||
                ext.compare(".IGS")  == 0 )
        {
            IGESControl_Reader reader;
            reader = IGESControl_Reader();
            reader.ReadFile(m_name.c_str());
            reader.NbRootsForTransfer();
            reader.TransferRoots();
            shape = reader.OneShape();
            if(shape.IsNull())
            {
                return false;
            }
        }
        else
        {
            return false;
        }
        
        TopTools_IndexedMapOfShape mapOfFaces;
        TopTools_IndexedMapOfShape mapOfEdges;
        TopExp::MapShapes(shape,TopAbs_FACE,mapOfFaces);
        
        for(int i = 1; i <= mapOfFaces.Extent(); i++)
        {
            TopoDS_Shape face= mapOfFaces.FindKey(i);
            
            TopTools_IndexedMapOfShape localEdges;
            TopExp::MapShapes(face,TopAbs_EDGE,localEdges);
            
            for(int j = 1; j <= localEdges.Extent(); j++)
            {
                TopoDS_Shape edge = localEdges.FindKey(j);
                BRepAdaptor_Curve curve = BRepAdaptor_Curve(TopoDS::Edge(edge));
                if(curve.GetType() != 7)
                {
                    if(!(mapOfEdges.Contains(edge)))
                    {
                        mapOfEdges.Add(edge);
                    }
                }
            }
            
        }
        
        vector<int> edgeTest;
        edgeTest.resize(mapOfEdges.Extent());
        
        for(int i=1; i<=mapOfEdges.Extent(); i++)
        {
            TopoDS_Shape edge = mapOfEdges.FindKey(i);
            
            AddCurve(i, edge);
            
            edgeTest[i-1]=0;
        }
        
        for(int i = 1; i <= mapOfFaces.Extent(); i++)
        {
            TopoDS_Shape face= mapOfFaces.FindKey(i);
            
            TopTools_IndexedMapOfShape localEdges;
            TopExp::MapShapes(face,TopAbs_EDGE,localEdges);
            
            vector<int> edges;
            
            for(int j = 1; j <= localEdges.Extent(); j++)
            {
                TopoDS_Shape edge = localEdges.FindKey(j);
                BRepAdaptor_Curve curve = BRepAdaptor_Curve(TopoDS::Edge(edge));
                if(curve.GetType() != 7)
                {
                    if(mapOfEdges.Contains(edge))
                    {
                        edges.push_back(mapOfEdges.FindIndex(edge));
                    }
                }
            }
             
            AddSurf(i, face, edges);
            
            for(int j = 0; j < edges.size(); j++)
            {
                edgeTest[edges[j]-1]++;
            }
            
        }
        int ct= 0;
        for(int i = 0; i < edgeTest.size(); i++)
        {
            if(edgeTest[i]!=2)
            {
                ct++;
            }
        }
        ASSERTL0(ct==0,"geometry will fail to mesh");
        
        m_numCurve = m_curves.size();
        m_numSurf = m_surfs.size();
        
        OrientateEdgesOnSurface();
        
        return true;
    }
    
    void CADSystem::AddCurve(int i, TopoDS_Shape in)
    {
        CADCurveSharedPtr newCurve =
            MemoryManager<CADCurve>::
                AllocateSharedPtr(i,in);
        m_curves.push_back(newCurve);
    }
    void CADSystem::AddSurf(int i, TopoDS_Shape in, std::vector<int> ein)
    {
        CADSurfSharedPtr newSurf =
            MemoryManager<CADSurf>::
                AllocateSharedPtr(i,in,ein);
        m_surfs.push_back(newSurf);
    }
    
    void CADSystem::OrientateEdgesOnSurface()
    {
        
        for(int i = 0; i < 2; i++)
        {
            vector<vector<int> > edgeloops;
            
            vector<int> edges = m_surfs[i]->GetEdges();
            
            vector<vector<int> > edgesfb;
            
            for(int j = 0; j < edges.size(); j++)
            {
                gp_Pnt start,end;
                m_curves[edges[j]-1]->GetMinMax(start,end);
                vector<int> fb;
                fb.resize(2);
                for(int k = 0 ; k < edges.size(); k++)
                {
                    if(j==k)
                        continue;
                    gp_Pnt starttest,endtest;
                    m_curves[edges[k]-1]->GetMinMax(starttest,endtest);
                    if(start.Distance(starttest)<1E-5 ||
                       start.Distance(endtest)<1E-5)
                    {
                        fb[0]=edges[k];
                    }
                    if(end.Distance(starttest)<1E-5 ||
                       end.Distance(endtest)<1E-5)
                    {
                        fb[1]=edges[k];
                    }
                }
                edgesfb.push_back(fb);
            }
            
            vector<pair<int,int> > edgelog(edgesfb.size());
            
            for(int j = 0; j < edges.size(); j++)
            {
                edgelog[j].first = edges[j];
                edgelog[j].second = 0;
            }
            
            for(int j = 0; j < edgesfb.size(); j++)
            {
                for(int k = 0; k < edges.size(); k++)
                {
                    if(edgelog[k].first == edgesfb[j][0])
                        edgelog[k].second++;
                    if(edgelog[k].first == edgesfb[j][1])
                        edgelog[k].second++;
                }
            }
            
            int ct= 0;
            for(int j = 0; j < edgelog.size(); j++)
            {
                if(edgelog[j].second !=2)
                {
                    ct++;
                }
            }
            ASSERTL0(ct==0,"error in connecting edges, cannot mesh");
            
            
        }
        
        exit(-1);
    }

}
}
