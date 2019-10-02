#include <LibUtilities/Python/NekPyConfig.hpp>
#include <FieldUtils/Module.h>
#include <boost/program_options.hpp>

using namespace Nektar;
using namespace Nektar::FieldUtils;

// Wrapper around Module::Process(&vm).
// Performs switching of m_comm if nparts > 1.
void Module_Process(ModuleSharedPtr m)
{
	if (m->m_f->m_nParts > 1)
	{
		if (m->GetModulePriority() == eOutput)
		{
			m->m_f->m_comm = m->m_f->m_partComm;
			if (m->GetModuleName() != "OutputInfo")
			{
				m->RegisterConfig("writemultiplefiles");
			}
		}
		else if (m->GetModulePriority() == eCreateGraph)
		{
			m->m_f->m_comm = m->m_f->m_partComm;
		}
		else
		{
			m->m_f->m_comm = m->m_f->m_defComm;
		}
	}
	m->SetDefaults();
	m->Process(m->m_f->m_vm);
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(RegisterConfig_overloads, RegisterConfig, 1, 2);

void export_Module()
{	
	py::class_<Module, std::shared_ptr<Module>, boost::noncopyable>
				("Module",  py::no_init)
		.def("GetModuleDescription", &Module::GetModuleDescription)
		.def("GetModulePriority", &Module::GetModulePriority)
		.def("PrintConfig", &Module::PrintConfig)
		.def("RegisterConfig", &Module::RegisterConfig, RegisterConfig_overloads())
		.def("SetDefaults", &Module::SetDefaults)
		.def("Run", &Module_Process)
	;
	
	py::class_<InputModule, py::bases<Module>, std::shared_ptr<InputModule>,
				boost::noncopyable>
				("InputModule", py::no_init)
			.def("AddFile", &InputModule::AddFile)
	;
	
	py::class_<ProcessModule, py::bases<Module>, std::shared_ptr<ProcessModule>,
				boost::noncopyable>
				("ProcessModule", py::no_init)
	;
	
	py::class_<OutputModule, py::bases<Module>, std::shared_ptr<OutputModule>,
				boost::noncopyable>
				("OutputModule", py::no_init)
	;
	NEKPY_WRAP_ENUM(ModulePriority, ModulePriorityMap);
}

