#include <SpatialDomains/SegGeom.h>
#include <NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

void export_SegGeom()
{
    py::class_<SegGeom, py::bases<Geometry1D>, std::shared_ptr<SegGeom> >(
        "SegGeom", py::init<>());
}
