#include <SpatialDomains/Geometry.h>
#include <NekPyConfig.hpp>

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
        .def("FillGeom",       &Geometry::FillGeom)
        .def("GetXmap",        &Geometry::GetXmap)
        .def("GenGeomFactors", &Geometry_GenGeomFactors)
        .def("ContainsPoint",  &Geometry_ContainsPoint)
        .def("GetCoeffs",      &Geometry::GetCoeffs,
             py::return_value_policy<py::copy_const_reference>())
        ;
}
