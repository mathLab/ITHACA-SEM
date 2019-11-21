////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description: Python wrapper for Geometry.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/Geometry2D.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

// Thin wrapper for ContainsPoint
bool Geometry_ContainsPoint(GeometrySharedPtr geom,
                            const Array<OneD, const NekDouble>& gloCoord)
{
    return geom->ContainsPoint(gloCoord);
}

void Geometry_GenGeomFactors(GeometrySharedPtr geom)
{
    GeomFactorsSharedPtr geomFactors = geom->GetGeomFactors();
}

void export_Geometry()
{
    py::class_<Geometry,
               std::shared_ptr<Geometry>,
               boost::noncopyable>(
                   "Geometry", py::no_init)

        .def("GetCoordim",     &Geometry::GetCoordim)
        .def("GetGlobalID",    &Geometry::GetGlobalID)

        .def("Setup",          &Geometry::Setup)
        .def("FillGeom",       &Geometry::FillGeom)
        .def("GenGeomFactors", &Geometry_GenGeomFactors)

        .def("ContainsPoint",  &Geometry_ContainsPoint)

        .def("GetVertex",      &Geometry::GetVertex)
        .def("GetEdge",        &Geometry::GetEdge)
        .def("GetFace",        &Geometry::GetFace)
        .def("GetVid",         &Geometry::GetVid)
        .def("GetEid",         &Geometry::GetEid)
        .def("GetFid",         &Geometry::GetFid)
        .def("GetTid",         &Geometry::GetTid)

        .def("GetNumVerts",    &Geometry::GetNumVerts)
        .def("GetNumEdges",    &Geometry::GetNumEdges)
        .def("GetNumFaces",    &Geometry::GetNumFaces)
        .def("GetShapeDim",    &Geometry::GetShapeDim)
        .def("GetShapeType",   &Geometry::GetShapeType)
        .def("GetEorient",     &Geometry::GetEorient)
        .def("GetForient",     &Geometry::GetForient)

        .def("GetXmap",        &Geometry::GetXmap)
        .def("GetCoeffs",      &Geometry::GetCoeffs,
             py::return_value_policy<py::copy_const_reference>())
        ;
}
