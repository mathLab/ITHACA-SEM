





#ifndef NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSYSTEM_H
#define NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSYSTEM_H

#include <string>

#include <boost/shared_ptr.hpp>

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>

#include <TColStd_HSequenceOfTransient.hxx>

#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>

#include <GeomAPI_ProjectPointOnSurf.hxx>

#include <BRepTools.hxx>
#include <BRep_Tool.hxx>

#include <gp_Trsf.hxx>

#include <TopLoc_Location.hxx>

namespace Nektar {
namespace LibUtilities {
    
    class CADCurve
    {
    public:
        CADCurve(int i, TopoDS_Shape in);
    
    private:
        int ID;
        BRepAdaptor_Curve occCurve;
    };
        
    class CADSurf
    {
    public:
        CADSurf(int i, TopoDS_Shape in, std::vector<int> ein);
        
    private:
        int ID;
        BRepAdaptor_Surface occSurface;
        Handle(Geom_Surface) s;
        std::vector<int> edges;
    };

	class CADSystem
	{
	public:
	    friend class MemoryManager<CADSystem>;

	    LIB_UTILITIES_EXPORT CADSystem(std::string name) : m_name(name)
	    {
	    }
	    LIB_UTILITIES_EXPORT std::string GetName();
        LIB_UTILITIES_EXPORT bool LoadCAD();

	private:
        
        void AddCurve(int i, TopoDS_Shape in)
        {
            CADCurve newCurve(i, in);
            m_curves.push_back(newCurve);
        }
        void AddSurf(int i, TopoDS_Shape in, std::vector<int> ein)
        {
            CADSurf newSurf(i, in, ein);
            m_surfs.push_back(newSurf);
        }
        
	    std::string m_name;
        int m_numCurve;
        int m_numSurf;
        std::vector<CADCurve> m_curves;
        std::vector<CADSurf> m_surfs;
	};

	typedef boost::shared_ptr<CADSystem> CADSystemSharedPtr;

}
}

#endif
