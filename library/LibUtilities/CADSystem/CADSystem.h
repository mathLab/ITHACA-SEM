////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.h
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


#ifndef NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSYSTEM_H
#define NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSYSTEM_H

#include <string>

#include <boost/shared_ptr.hpp>

#include <LibUtilities/CADSystem/CADCurve.h>
#include <LibUtilities/CADSystem/CADSurf.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
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
        
    

	class CADSystem
	{
	public:
	    friend class MemoryManager<CADSystem>;

	    LIB_UTILITIES_EXPORT CADSystem(const std::string &name) : m_name(name)
	    {
	    };
        
	    LIB_UTILITIES_EXPORT std::string GetName();
        LIB_UTILITIES_EXPORT bool LoadCAD();
        LIB_UTILITIES_EXPORT void Report();
        LIB_UTILITIES_EXPORT void GetBoundingBox(Array<OneD, NekDouble>& out);
        LIB_UTILITIES_EXPORT int GetNumSurf(){return m_numSurf;}
        LIB_UTILITIES_EXPORT int GetNumCurve(){return m_numCurve;}
        LIB_UTILITIES_EXPORT void GetParameterPlaneBounds(int i,
                                                    Array<OneD, NekDouble>& out);
        LIB_UTILITIES_EXPORT void N(int i, NekDouble u, NekDouble v,
                                    Array<OneD, NekDouble>& out);
        LIB_UTILITIES_EXPORT void D1(int i, NekDouble u, NekDouble v,
                                    Array<OneD, NekDouble>& out);
        LIB_UTILITIES_EXPORT void D2(int i, NekDouble u, NekDouble v,
                                    Array<OneD, NekDouble>& out);
        LIB_UTILITIES_EXPORT const CADCurveSharedPtr GetCurve(int i )
        {
            return m_curves[i-1];
        }
        LIB_UTILITIES_EXPORT CADSurfSharedPtr GetSurf(int i)
        {
            return m_surfs[i-1];
        }

	private:
        
        void AddCurve(int i, TopoDS_Shape in);
        void AddSurf(int i, TopoDS_Shape in, std::vector<int> ein);
        
        void OrientateEdgesOnSurface();
        
	    std::string m_name;
        int m_numCurve;
        int m_numSurf;
        std::vector<CADCurveSharedPtr> m_curves;
        std::vector<CADSurfSharedPtr> m_surfs;
	};

	typedef boost::shared_ptr<CADSystem> CADSystemSharedPtr;

}
}

#endif
