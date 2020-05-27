#include <LibUtilities/Python/NekPyConfig.hpp>
#include <FieldUtils/OutputModules/OutputFld.h>
#include <FieldUtils/OutputModules/OutputInfo.h>
#include <FieldUtils/OutputModules/OutputPts.h>
#include <FieldUtils/OutputModules/OutputTecplot.h>
#include <FieldUtils/OutputModules/OutputVtk.h>
#include <FieldUtils/OutputModules/OutputXml.h>


using namespace Nektar;
using namespace Nektar::FieldUtils;

void SetOutputFile(std::shared_ptr<OutputModule> m, std::string filename)
{
	m->RegisterConfig("outfile", filename);
}

std::shared_ptr<OutputFld> OutputFld_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<OutputFld> m = MemoryManager<OutputFld>::AllocateSharedPtr(f);
	
	if (filename.size()) 
	{
		SetOutputFile(m, filename);
	}
	
    return m;
}

std::shared_ptr<OutputInfo> OutputInfo_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<OutputInfo> m = MemoryManager<OutputInfo>::AllocateSharedPtr(f);
	std::string nparts = std::to_string(m->m_f->m_nParts);
	m->RegisterConfig("nparts", nparts);
	if (filename.size())
	{
		SetOutputFile(m, filename);
	}
	return m;
}

std::shared_ptr<OutputPts> OutputPts_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<OutputPts> m = MemoryManager<OutputPts>::AllocateSharedPtr(f);
	if (filename.size())
	{
		SetOutputFile(m, filename);
	}
	return m;
}

std::shared_ptr<OutputTecplot> OutputTecplot_Init(FieldSharedPtr f, std::string filename = "", 
												bool doubleprecision = false)
{
	std::shared_ptr<OutputTecplot> m = MemoryManager<OutputTecplot>::AllocateSharedPtr(f);
	if (filename.size())
	{
		SetOutputFile(m, filename);
	}
	if (doubleprecision)
	{
		m->RegisterConfig("double");
	}
	return m;
}

std::shared_ptr<OutputTecplotBinary> OutputTecplotBinary_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<OutputTecplotBinary> m = MemoryManager<OutputTecplotBinary>::AllocateSharedPtr(f);
	if (filename.size())
	{
		SetOutputFile(m, filename);
	}
	return m;
}

std::shared_ptr<OutputVtk> OutputVtk_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<OutputVtk> m = MemoryManager<OutputVtk>::AllocateSharedPtr(f);
	if (filename.size()) 
	{
		SetOutputFile(m, filename);
	}
    return m;
}

std::shared_ptr<OutputXml> OutputXml_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<OutputXml> m = MemoryManager<OutputXml>::AllocateSharedPtr(f);
	if (filename.size()) 
	{
		SetOutputFile(m, filename);
	}
    return m;
}

void export_OutputModules()
{
	py::class_<OutputFileBase, boost::noncopyable, py::bases<OutputModule>, 
				std::shared_ptr<OutputFld> >
				("OutputFileBase", py::no_init)
	;
	py::class_<OutputFld, boost::noncopyable, py::bases<OutputFileBase>, 
				std::shared_ptr<OutputFld> >
				("OutputFld", py::no_init)
		.def("__init__", py::make_constructor(&OutputFld_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetOutputFile", &SetOutputFile)
	;
	py::class_<OutputInfo, boost::noncopyable, py::bases<OutputModule>,
				std::shared_ptr<OutputInfo> >
				("OutputInfo", py::no_init)
		.def("__init__", py::make_constructor(&OutputInfo_Init,
			py::default_call_policies(), (py::arg("f"),
			py::arg("filename") = "")))
		.def("SetOutputFile", &SetOutputFile)
	;
	py::class_<OutputPts, boost::noncopyable, py::bases<OutputFileBase>, 
				std::shared_ptr<OutputPts> >
				("OutputPts", py::no_init)
		.def("__init__", py::make_constructor(&OutputPts_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetOutputFile", &SetOutputFile)
	;
	py::class_<OutputTecplot, boost::noncopyable, py::bases<OutputFileBase>, 
				std::shared_ptr<OutputTecplot> >
				("OutputTecplot", py::no_init)
		.def("__init__", py::make_constructor(&OutputTecplot_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "", py::arg("double") = false)))
		.def("SetOutputFile", &SetOutputFile)
	;
	py::class_<OutputTecplotBinary, boost::noncopyable, py::bases<OutputTecplot>, 
				std::shared_ptr<OutputTecplotBinary> >
				("OutputTecplotBinary", py::no_init)
		.def("__init__", py::make_constructor(&OutputTecplotBinary_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetOutputFile", &SetOutputFile)
	;
	py::class_<OutputVtk, boost::noncopyable, py::bases<OutputFileBase>, 
				std::shared_ptr<OutputVtk> >
				("OutputVtk", py::no_init)
		.def("__init__", py::make_constructor(&OutputVtk_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetOutputFile", &SetOutputFile)
	;
	py::class_<OutputXml, boost::noncopyable, py::bases<OutputModule>, 
				std::shared_ptr<OutputXml> >
				("OutputXml", py::no_init)
		.def("__init__", py::make_constructor(&OutputXml_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetOutputFile", &SetOutputFile)
	;

}
