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

#include <boost/filesystem.hpp>

#include <LibUtilities/CADSystem/CADSystem.h>

using namespace std;
namespace Nektar{
namespace LibUtilities{
    
    CADCurve::CADCurve(int i, TopoDS_Shape in) : ID(i)
    {
        gp_Trsf transform;
        gp_Pnt ori(0.0,0.0,0.0);
        transform.SetScale(ori,1.0/1000.0);
        TopLoc_Location mv(transform);
        in.Move(mv);
        occCurve = BRepAdaptor_Curve(TopoDS::Edge(in));
    }
    
    void CADCurve::GetMinMax(gp_Pnt &start, gp_Pnt &end)
    {
        start = occCurve.Value(occCurve.FirstParameter());
        end  = occCurve.Value(occCurve.LastParameter());
    }
    
    CADSurf::CADSurf(int i, TopoDS_Shape in, vector<int> ein) : ID(i), edges(ein)
    {
        gp_Trsf transform;
        gp_Pnt ori(0.0,0.0,0.0);
        transform.SetScale(ori,1.0/1000.0);
        TopLoc_Location mv(transform);
        s = BRep_Tool::Surface(TopoDS::Face(in));
        in.Move(mv);
        occSurface = BRepAdaptor_Surface(TopoDS::Face(in));
    }
    
    Array<OneD, NekDouble> CADSurf::N(NekDouble u, NekDouble v)
    {
        Array<OneD, NekDouble> out(3);
        gp_Pnt Loc;
        gp_Vec D1U,D1V;
        occSurface.D1(u,v,Loc,D1U,D1V);
        gp_Vec n = D1U.Crossed(D1V);
        if(n.X()==0 && n.Y() ==0 && n.Z()==0)
        {
            out[0]=0.0;
            out[1]=0.0;
            out[2]=0.0;
        }
        else
        {
            n.Normalize();
            out[0]=n.X();
            out[1]=n.Y();
            out[2]=n.Z();
        }
        
        return out;
    }
    
    Array<OneD, NekDouble> CADSurf::D1(NekDouble u, NekDouble v)
    {
        Array<OneD, NekDouble> out(9);
        gp_Pnt Loc;
        gp_Vec D1U,D1V;
        occSurface.D1(u,v,Loc,D1U,D1V);

        out[0]=Loc.X();
        out[1]=Loc.Y();
        out[2]=Loc.Z();
        out[3]=D1U.X();
        out[4]=D1U.Y();
        out[5]=D1U.Z();
        out[6]=D1V.X();
        out[7]=D1V.Y();
        out[8]=D1V.Z();
        
        return out;
    }
    
    Array<OneD, NekDouble> CADSurf::D2(NekDouble u, NekDouble v)
    {
        Array<OneD, NekDouble> out(18);
        gp_Pnt Loc;
        gp_Vec D1U,D1V,D2U,D2V,D2UV;
        occSurface.D2(u,v,Loc,D1U,D1V,D2U,D2V,D2UV);
        
        out[0]=Loc.X();
        out[1]=Loc.Y();
        out[2]=Loc.Z();
        out[3]=D1U.X();
        out[4]=D1U.Y();
        out[5]=D1U.Z();
        out[6]=D1V.X();
        out[7]=D1V.Y();
        out[8]=D1V.Z();
        out[9]=D2U.X();
        out[10]=D2U.Y();
        out[11]=D2U.Z();
        out[12]=D2V.X();
        out[13]=D2V.Y();
        out[14]=D2V.Z();
        out[15]=D2UV.X();
        out[16]=D2UV.Y();
        out[17]=D2UV.Z();
        
        return out;
    }
    
    void CADSystem::GetParameterPlaneBounds(int i,
                                            Array<OneD, NekDouble>& out)
    {
        out = Array<OneD, NekDouble>(4);
        out[0]=m_surfs[i-1].minU();
        out[1]=m_surfs[i-1].maxU();
        out[2]=m_surfs[i-1].minV();
        out[3]=m_surfs[i-1].maxV();
    }

    string CADSystem::GetName()
    {
        return m_name;
    }
    
    void CADSystem::N(int i, NekDouble u, NekDouble v,
                      Array<OneD, NekDouble>& out)
    {
        ASSERTL0(u>=m_surfs[i-1].minU() &&
                 u<=m_surfs[i-1].maxU() &&
                 v>=m_surfs[i-1].minV() &&
                 v<=m_surfs[i-1].maxV(), "(u,v) out of bounds");
        out = m_surfs[i-1].N(u,v);
    }
    void CADSystem::D1(int i, NekDouble u, NekDouble v,
                                 Array<OneD, NekDouble>& out)
    {
        ASSERTL0(u>=m_surfs[i-1].minU() &&
                 u<=m_surfs[i-1].maxU() &&
                 v>=m_surfs[i-1].minV() &&
                 v<=m_surfs[i-1].maxV(), "(u,v) out of bounds");
        out = m_surfs[i-1].D1(u,v);
    }
    void CADSystem::D2(int i, NekDouble u, NekDouble v,
                                 Array<OneD, NekDouble>& out)
    {
        ASSERTL0(u>=m_surfs[i-1].minU() &&
                 u<=m_surfs[i-1].maxU() &&
                 v>=m_surfs[i-1].minV() &&
                 v<=m_surfs[i-1].maxV(), "(u,v) out of bounds");
        out = m_surfs[i-1].D2(u,v);
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
            m_curves[i].GetMinMax(start,end);
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
        
        TopTools_IndexedMapOfShape unfilteredEdges;
        TopTools_IndexedMapOfShape mapOfFaces;
        TopTools_IndexedMapOfShape mapOfEdges;
        TopExp::MapShapes(shape,TopAbs_FACE,mapOfFaces);
        TopExp::MapShapes(shape,TopAbs_EDGE,unfilteredEdges);
        
        for(int i=1; i<=unfilteredEdges.Extent(); i++)
        {
            TopoDS_Shape edge = unfilteredEdges.FindKey(i);
            BRepAdaptor_Curve curve = BRepAdaptor_Curve(TopoDS::Edge(edge));
            if(curve.GetType() == 7)
            {
                continue;
            }
            mapOfEdges.Add(edge);
        }

        m_numCurve = mapOfFaces.Extent();
        m_numSurf =  mapOfEdges.Extent();
        
        for(int i=1; i<=mapOfEdges.Extent(); i++)
        {
            TopoDS_Shape edge = mapOfEdges.FindKey(i);
            
            AddCurve(i, edge);
        }
        
        for(int i = 1; i <= mapOfFaces.Extent(); i++)
        {
            TopoDS_Shape face= mapOfFaces.FindKey(i);
            
            TopTools_IndexedMapOfShape localEdges;
            TopTools_IndexedMapOfShape unfilteredLocalEdges;
            TopExp::MapShapes(face,TopAbs_EDGE,unfilteredLocalEdges);

            for(int j = 1; j <= unfilteredLocalEdges.Extent(); j++)
            {
                TopoDS_Shape edge = unfilteredLocalEdges.FindKey(j);
                BRepAdaptor_Curve curve = BRepAdaptor_Curve(TopoDS::Edge(edge));
                if(curve.GetType() != 7)
                {
                    localEdges.Add(edge);
                }
            }
            
            vector<int> edges;
            edges.resize(localEdges.Extent());
            
            for(int j=0; j<localEdges.Extent(); j++)
            {
                edges[j] = mapOfEdges.FindIndex(localEdges.FindKey(j+1));
            }
            
            AddSurf(i, face, edges);
        }
        
        m_numCurve = m_curves.size();
        m_numSurf = m_surfs.size();
        
        return true;
    }
    
    

}
}
