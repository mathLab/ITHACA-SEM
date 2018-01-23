#include <SpatialDomains/QuadGeom.h>
#include <NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

void export_QuadGeom()
{
    py::class_<QuadGeom, py::bases<Geometry2D>, std::shared_ptr<QuadGeom> >(
        "QuadGeom", py::init<>());
}
