#include <NekPyConfig.hpp>
#include <StdRegions/StdTriExp.h>

using namespace Nektar;
using namespace Nektar::StdRegions;

void export_StdTriExp()
{
    py::class_<StdTriExp, py::bases<StdExpansion>,
               std::shared_ptr<StdTriExp> >(
                   "StdTriExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());
}
