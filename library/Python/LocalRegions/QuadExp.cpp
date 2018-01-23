#include <NekPyConfig.hpp>
#include <LocalRegions/QuadExp.h>

using namespace Nektar;
using namespace Nektar::LocalRegions;

void export_QuadExp()
{
    py::class_<QuadExp, py::bases<Expansion, StdRegions::StdQuadExp>,
               std::shared_ptr<QuadExp> >(
                   "QuadExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&,
                   const SpatialDomains::QuadGeomSharedPtr &>());
}
