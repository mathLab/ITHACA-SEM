#include <NekPyConfig.hpp>
#include <LocalRegions/Expansion.h>

using namespace Nektar;
using namespace Nektar::StdRegions;
using namespace Nektar::LocalRegions;
using namespace Nektar::SpatialDomains;

void export_Expansion()
{
    py::class_<Expansion,
               std::shared_ptr<Expansion>, py::bases<StdExpansion>,
               boost::noncopyable>(
                   "Expansion", py::no_init)

        .def("GetGeom", &Expansion::GetGeom)
        ;
}
