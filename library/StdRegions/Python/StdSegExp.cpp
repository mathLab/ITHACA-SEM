#include <NekPyConfig.hpp>
#include <StdRegions/StdSegExp.h>

using namespace Nektar;
using namespace Nektar::StdRegions;

void export_StdSegExp()
{
    py::class_<StdSegExp, py::bases<StdExpansion>,
               std::shared_ptr<StdSegExp> >(
                   "StdSegExp", py::init<const LibUtilities::BasisKey&>());
}
