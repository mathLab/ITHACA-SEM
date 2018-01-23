#include <SpatialDomains/Geometry2D.h>
#include <NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

void export_Geometry2D()
{
    py::class_<Geometry2D, py::bases<Geometry>, std::shared_ptr<Geometry2D>,
               boost::noncopyable>(
                   "Geometry2D", py::no_init);
}
