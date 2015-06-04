
#include "occ.h"
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

    string CADSystem::GetName()
    {
        return m_name;
    }

    bool CADSystem::LoadCAD()
    {
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
        
        return true;
    }
    
    

}
}
