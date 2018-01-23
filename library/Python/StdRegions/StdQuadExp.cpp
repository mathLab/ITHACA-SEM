#include <NekPyConfig.hpp>
#include <StdRegions/StdQuadExp.h>

using namespace Nektar;
using namespace Nektar::StdRegions;

void export_StdQuadExp()
{
    py::class_<StdQuadExp, py::bases<StdExpansion>,
               std::shared_ptr<StdQuadExp> >(
                   "StdQuadExp", py::init<const LibUtilities::BasisKey&,
                   const LibUtilities::BasisKey&>());
}
