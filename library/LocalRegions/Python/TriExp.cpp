#include <NekPyConfig.hpp>
#include <LocalRegions/TriExp.h>

using namespace Nektar;
using namespace Nektar::LocalRegions;

void export_TriExp()
{
    py::class_<TriExp, py::bases<Expansion, StdRegions::StdTriExp>,
               std::shared_ptr<TriExp> >(
                   "TriExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const SpatialDomains::TriGeomSharedPtr &>());
}
