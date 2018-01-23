#include <NekPyConfig.hpp>
#include <MultiRegions/ExpList2D.h>

using namespace Nektar;
using namespace Nektar::MultiRegions;

void export_ExpList2D()
{
    py::class_<ExpList2D, py::bases<ExpList>,
               std::shared_ptr<ExpList2D> >(
                   "ExpList2D", py::init<
                   const LibUtilities::SessionReaderSharedPtr &,
                   const SpatialDomains::MeshGraphSharedPtr &>());
}
