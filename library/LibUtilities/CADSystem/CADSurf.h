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


#ifndef NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSURF_H
#define NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSURF_H

#include <string>

#include <boost/shared_ptr.hpp>

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
        
        class CADSurf
        {
        public:
            friend class MemoryManager<CADSurf>;
            
            CADSurf(int i, TopoDS_Shape in, std::vector<int> ein);
            NekDouble minU(){return occSurface.FirstUParameter();}
            NekDouble maxU(){return occSurface.LastUParameter();}
            NekDouble minV(){return occSurface.FirstVParameter();}
            NekDouble maxV(){return occSurface.LastVParameter();}
            Array<OneD, NekDouble> N(NekDouble u, NekDouble v);
            Array<OneD, NekDouble> D1(NekDouble u, NekDouble v);
            Array<OneD, NekDouble> D2(NekDouble u, NekDouble v);
            
        private:
            int ID;
            BRepAdaptor_Surface occSurface;
            Handle(Geom_Surface) s;
            std::vector<int> edges;
        };
        
        typedef boost::shared_ptr<CADSurf> CADSurfSharedPtr;
    }
}

#endif