#include <SpatialDomains/Geometry1D.h>
#include <NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

void export_Geometry1D()
{
    py::class_<Geometry1D, py::bases<Geometry>, std::shared_ptr<Geometry1D>,
               boost::noncopyable>(
                   "Geometry1D", py::no_init);
}
