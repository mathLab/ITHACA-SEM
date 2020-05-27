#include <LibUtilities/Python/NekPyConfig.hpp>
#include <FieldUtils/InputModules/InputDat.h>
#include <FieldUtils/InputModules/InputFld.h>
#include <FieldUtils/InputModules/InputNek5000.h>
#include <FieldUtils/InputModules/InputPts.h>
#include <FieldUtils/InputModules/InputSemtex.h>
#include <FieldUtils/InputModules/InputXml.h>


using namespace Nektar;
using namespace Nektar::FieldUtils;

std::string GetExt(std::string filename)
{
	std::string ext;
	int dot = filename.find_last_of('_.') + 1;
	ext = filename.substr(dot, filename.length() - dot);
		
	if(ext.back() == fs::path::preferred_separator) 
	{
		ext.pop_back();
	}
		
	if(ext == "gz")
	{
            std::string tmp = filename.substr(0,dot-1);
            dot = tmp.find_last_of('_.') + 1;
            ext = filename.substr(dot,filename.length()-dot);
	}
	return ext;
	
}

void SetInputDatFile(std::shared_ptr<InputDat> m, std::string filename)
{
	m->AddFile("dat", filename);
	m->RegisterConfig("infile", filename);
}

std::shared_ptr<InputDat> InputDat_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<InputDat> m = MemoryManager<InputDat>::AllocateSharedPtr(f);
	if (filename.size()) 
	{
		SetInputDatFile(m, filename);
	}
    return m;
}

void SetInputFldFile(std::shared_ptr<InputFld> m, std::string filename)
{
	m->AddFile("fld", filename);
	m->RegisterConfig("infile", filename);
}

std::shared_ptr<InputFld> InputFld_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<InputFld> m = MemoryManager<InputFld>::AllocateSharedPtr(f);
	if (filename.size()) 
	{
		SetInputFldFile(m, filename);
	}
    return m;
}

void SetInputNek5000File(std::shared_ptr<InputNek5000> m, std::string filename)
{
	m->AddFile("fld5000", filename);
	m->RegisterConfig("infile", filename);
}

std::shared_ptr<InputNek5000> InputNek5000_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<InputNek5000> m = MemoryManager<InputNek5000>::AllocateSharedPtr(f);
	if (filename.size()) 
	{
		SetInputNek5000File(m, filename);
	}
    return m;
}

void SetInputPtsFile(std::shared_ptr<InputPts> m, std::string filename)
{
	std::string ext = GetExt(filename);
	if (ext == "pts")
	{
		m->AddFile("pts", filename);
	}
	else if (ext == "csv")
	{
		m->AddFile("csv", filename);
	}
	else
	{
		std::cout << "Assuming a .pts file." << std::endl;
		m->AddFile("pts", filename);
	}
	m->RegisterConfig("infile", filename);
}

std::shared_ptr<InputPts> InputPts_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<InputPts> m = MemoryManager<InputPts>::AllocateSharedPtr(f);
	if (filename.size()) 
	{
		SetInputPtsFile(m, filename);
	}
    return m;
}

void SetInputSemtexFile(std::shared_ptr<InputSemtex> m, std::string filename)
{
	m->AddFile("fldsem", filename);
	m->RegisterConfig("infile", filename);
}

std::shared_ptr<InputSemtex> InputSemtex_Init(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<InputSemtex> m = MemoryManager<InputSemtex>::AllocateSharedPtr(f);
	if (filename.size()) 
	{
		SetInputSemtexFile(m, filename);
	}
    return m;
}

void SetInputXmlFile(std::shared_ptr<InputXml> m, std::string filename)
{
	m->AddFile("xml", filename);
	m->RegisterConfig("infile", filename);
}

std::shared_ptr<InputXml> InputXml_Init_list(FieldSharedPtr f, py::list filenames)
{
	std::shared_ptr<InputXml> m = MemoryManager<InputXml>::AllocateSharedPtr(f);
	for (int i = 0; i < py::len(filenames); ++i)
	{
		std::string filename = py::extract<std::string>(filenames[i]);
		m->AddFile("xml", filename);
	}
    return m;
}

std::shared_ptr<InputXml> InputXml_Init_string(FieldSharedPtr f, std::string filename = "")
{
	std::shared_ptr<InputXml> m = MemoryManager<InputXml>::AllocateSharedPtr(f);
	if (filename.size()) 
	{
		SetInputXmlFile(m, filename);
	}
    return m;
}

void export_InputModules()
{
	py::class_<InputDat, py::bases<InputModule>, 
				std::shared_ptr<InputDat> >
				("InputDat", py::no_init)
		.def("__init__", py::make_constructor(&InputDat_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetInputFile", &SetInputDatFile)
	;
	py::class_<InputFld, py::bases<InputModule>, 
			std::shared_ptr<InputFld> >("InputFld", py::no_init)
		.def("__init__", py::make_constructor(&InputFld_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetInputFile", &SetInputFldFile)
	;
	py::class_<InputNek5000, py::bases<InputModule>, 
				std::shared_ptr<InputNek5000> >
				("InputNek5000", py::no_init)
		.def("__init__", py::make_constructor(&InputNek5000_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetInputFile", &SetInputNek5000File)
	;
	py::class_<InputPts, py::bases<InputModule>, 
				std::shared_ptr<InputPts> >
				("InputPts", py::no_init)
		.def("__init__", py::make_constructor(&InputPts_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetInputFile", &SetInputPtsFile)
	;
	py::class_<InputSemtex, py::bases<InputModule>, 
				std::shared_ptr<InputSemtex> >
				("InputSemtex", py::no_init)
		.def("__init__", py::make_constructor(&InputSemtex_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("SetInputFile", &SetInputSemtexFile)
	;
	py::class_<InputXml, py::bases<InputModule>, 
				std::shared_ptr<InputXml> >
				("InputXml", py::no_init)
		.def("__init__", py::make_constructor(&InputXml_Init_string, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("filename") = "")))
		.def("__init__", py::make_constructor(&InputXml_Init_list, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("files"))))
		.def("SetInputFile", &SetInputXmlFile)
	;
}
