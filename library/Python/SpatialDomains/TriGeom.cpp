#include <SpatialDomains/TriGeom.h>
#include <NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

void export_TriGeom()
{
    py::class_<TriGeom, py::bases<Geometry2D>, std::shared_ptr<TriGeom> >(
        "TriGeom", py::init<>());
}
