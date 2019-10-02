#include <LibUtilities/Python/NekPyConfig.hpp>
#include <LibUtilities/Communication/Comm.h>

using namespace Nektar;
using namespace Nektar::LibUtilities;

void export_Comm()
{
	py::class_<Comm, std::shared_ptr<Comm>, boost::noncopyable>
		("Comm", py::no_init)
	;
}
