#include <NekPyConfig.hpp>
#include <LocalRegions/SegExp.h>

using namespace Nektar;
using namespace Nektar::LocalRegions;

void export_SegExp()
{
    py::class_<SegExp, py::bases<Expansion, StdRegions::StdSegExp>,
               std::shared_ptr<SegExp> >(
                   "SegExp", py::init<const LibUtilities::BasisKey&,
                   const SpatialDomains::Geometry1DSharedPtr &>());
}
