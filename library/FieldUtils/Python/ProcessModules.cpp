#include <LibUtilities/Python/NekPyConfig.hpp>
#include <FieldUtils/ProcessModules/ProcessAddCompositeID.h>
#include <FieldUtils/ProcessModules/ProcessAddFld.h>
#include <FieldUtils/ProcessModules/ProcessBoundaryExtract.h>
#include <FieldUtils/ProcessModules/ProcessC0Projection.h>
#include <FieldUtils/ProcessModules/ProcessCombineAvg.h>
#include <FieldUtils/ProcessModules/ProcessCreateExp.h>
#include <FieldUtils/ProcessModules/ProcessDeform.h>
#include <FieldUtils/ProcessModules/ProcessDisplacement.h>
#include <FieldUtils/ProcessModules/ProcessDOF.h>
#include <FieldUtils/ProcessModules/ProcessEquiSpacedOutput.h>
#include <FieldUtils/ProcessModules/ProcessFieldFromString.h>
#include <FieldUtils/ProcessModules/ProcessGrad.h>
#include <FieldUtils/ProcessModules/ProcessHomogeneousPlane.h>
#include <FieldUtils/ProcessModules/ProcessHomogeneousStretch.h>
#include <FieldUtils/ProcessModules/ProcessInnerProduct.h>
#include <FieldUtils/ProcessModules/ProcessInterpField.h>
#include <FieldUtils/ProcessModules/ProcessInterpPointDataToFld.h>
#include <FieldUtils/ProcessModules/ProcessInterpPoints.h>
#include <FieldUtils/ProcessModules/ProcessInterpPtsToPts.h>
#include <FieldUtils/ProcessModules/ProcessIsoContour.h>
#include <FieldUtils/ProcessModules/ProcessJacobianEnergy.h>
#include <FieldUtils/ProcessModules/ProcessL2Criterion.h>
#include <FieldUtils/ProcessModules/ProcessMapping.h>
#include <FieldUtils/ProcessModules/ProcessMean.h>
#include <FieldUtils/ProcessModules/ProcessMeanMode.h>
#include <FieldUtils/ProcessModules/ProcessMultiShear.h>
#include <FieldUtils/ProcessModules/ProcessNumModes.h>
#include <FieldUtils/ProcessModules/ProcessPointDataToFld.h>
#include <FieldUtils/ProcessModules/ProcessPrintFldNorms.h>
#include <FieldUtils/ProcessModules/ProcessQCriterion.h>
#include <FieldUtils/ProcessModules/ProcessQualityMetric.h>
#include <FieldUtils/ProcessModules/ProcessRemoveField.h>
#include <FieldUtils/ProcessModules/ProcessScaleInFld.h>
#include <FieldUtils/ProcessModules/ProcessScalGrad.h>
#include <FieldUtils/ProcessModules/ProcessStreamFunction.h>
#include <FieldUtils/ProcessModules/ProcessSurfDistance.h>
#include <FieldUtils/ProcessModules/ProcessVorticity.h>
#include <FieldUtils/ProcessModules/ProcessWSS.h>


using namespace Nektar;
using namespace Nektar::FieldUtils;

std::shared_ptr<ProcessAddFld> ProcessAddFld_Init(FieldSharedPtr f,
										std::string scale = "",
										std::string fromfld = "")
{
	std::shared_ptr<ProcessAddFld> m = 
			MemoryManager<ProcessAddFld>::AllocateSharedPtr(f);
	
	if (scale.size()) 
	{
		m->RegisterConfig("scale", scale);
	}
	if (fromfld.size()) 
	{
		m->RegisterConfig("fromfld", fromfld);
	}
    return m;
}

std::shared_ptr<ProcessBoundaryExtract> ProcessBoundaryExtract_Init(FieldSharedPtr f,
										std::string bnd = "",
										bool addnormals = false)
{
	std::shared_ptr<ProcessBoundaryExtract> m = 
			MemoryManager<ProcessBoundaryExtract>::AllocateSharedPtr(f);
	
	if (bnd.size()) 
	{
		m->RegisterConfig("bnd", bnd);
	}
	if (addnormals) 
	{
		m->RegisterConfig("addnormals");
	}
    return m;
}



std::shared_ptr<ProcessC0Projection> ProcessC0Projection_Init(FieldSharedPtr f,
										bool localtoglobalmap = false,
										bool usexmlbcs = false,
										std::string helmsmoothing = "")
{
	std::shared_ptr<ProcessC0Projection> m = 
			MemoryManager<ProcessC0Projection>::AllocateSharedPtr(f);
	
	if (localtoglobalmap) 
	{
		m->RegisterConfig("localtoglobalmap");
	}
	if (usexmlbcs) 
	{
		m->RegisterConfig("usexmlbcs");
	}
	if (helmsmoothing.size()) 
	{
		m->RegisterConfig("helmsmoothing", helmsmoothing);
	}
    return m;
}

std::shared_ptr<ProcessCombineAvg> ProcessCombineAvg_Init(FieldSharedPtr f,
										std::string fromfld = "")
{
	std::shared_ptr<ProcessCombineAvg> m = 
			MemoryManager<ProcessCombineAvg>::AllocateSharedPtr(f);
	
	if (fromfld.size()) 
	{
		m->RegisterConfig("fromfld", fromfld);
	}
    return m;
}

std::shared_ptr<ProcessDisplacement> ProcessDisplacement_Init(FieldSharedPtr f,
										std::string to = "",
										bool usevertexids = false)
{
	std::shared_ptr<ProcessDisplacement> m = 
			MemoryManager<ProcessDisplacement>::AllocateSharedPtr(f);
	
	if (to.size()) 
	{
		m->RegisterConfig("to", to);
	}
	if (usevertexids) 
	{
		m->RegisterConfig("usevertexids");
	}
    return m;
}

std::shared_ptr<ProcessEquiSpacedOutput> ProcessEquiSpacedOutput_Init(FieldSharedPtr f,
										bool tetonly = false,
										bool modalenergy = false)
{
	std::shared_ptr<ProcessEquiSpacedOutput> m = 
			MemoryManager<ProcessEquiSpacedOutput>::AllocateSharedPtr(f);
	
	if (tetonly) 
	{
		m->RegisterConfig("tetonly");
	}
	if (modalenergy) 
	{
		m->RegisterConfig("modalenergy");
	}
    return m;
}

std::shared_ptr<ProcessFieldFromString> ProcessFieldFromString_Init(FieldSharedPtr f,
										std::string fieldstr = "",
										std::string fieldname = "")
{
	std::shared_ptr<ProcessFieldFromString> m = 
			MemoryManager<ProcessFieldFromString>::AllocateSharedPtr(f);
	
	if (fieldstr.size()) 
	{
		m->RegisterConfig("fieldstr", fieldstr);
	}
	if (fieldname.size()) 
	{
		m->RegisterConfig("fieldname", fieldname);
	}
    return m;
}

std::shared_ptr<ProcessHomogeneousPlane> ProcessHomogeneousPlane_Init(FieldSharedPtr f,
										std::string planeid = "",
										bool wavespace = "")
{
	std::shared_ptr<ProcessHomogeneousPlane> m = 
			MemoryManager<ProcessHomogeneousPlane>::AllocateSharedPtr(f);
	
	if (planeid.size()) 
	{
		m->RegisterConfig("planeid", planeid);
	}
	if (wavespace) 
	{
		m->RegisterConfig("wavespace");
	}
    return m;
}

std::shared_ptr<ProcessHomogeneousStretch> ProcessHomogeneousStretch_Init(FieldSharedPtr f,
										std::string factor = "")
{
	std::shared_ptr<ProcessHomogeneousStretch> m = 
			MemoryManager<ProcessHomogeneousStretch>::AllocateSharedPtr(f);
	
	if (factor.size()) {
		m->RegisterConfig("factor", factor);
	}
    return m;
}

std::shared_ptr<ProcessInnerProduct> ProcessInnerProduct_Init(FieldSharedPtr f,
										std::string fromfld = "",
										std::string fields = "",
										std::string multifldids = "",
										bool allfromflds = false)
{
	std::shared_ptr<ProcessInnerProduct> m = 
			MemoryManager<ProcessInnerProduct>::AllocateSharedPtr(f);
	
	if (fromfld.size()) 
	{
		m->RegisterConfig("fromfld", fromfld);
	}
	if (fields.size()) 
	{
		m->RegisterConfig("fields", fields);
	}
	if (multifldids.size()) 
	{
		m->RegisterConfig("multifldids", multifldids);
	}
	if (allfromflds) 
	{
		m->RegisterConfig("allfromflds");
	}
    return m;
}

std::shared_ptr<ProcessInterpField> ProcessInterpField_Init(FieldSharedPtr f,
									std::string fromxml = "",
									std::string fromfld = "",
									std::string clamptolowervalue = "",
									std::string clamptouppervalue = "",
									std::string defaultvalue = "")
{
	std::shared_ptr<ProcessInterpField> m = 
			MemoryManager<ProcessInterpField>::AllocateSharedPtr(f);
	
	if (fromxml.size()) 
	{
		m->RegisterConfig("fromxml", fromxml);
	}
	if (fromfld.size()) 
	{
		m->RegisterConfig("fromfld", fromfld);
	}
	if (clamptolowervalue.size()) 
	{
		m->RegisterConfig("clamptolowervalue", clamptolowervalue);
	}
	if (clamptouppervalue.size()) 
	{
		m->RegisterConfig("clamptouppervalue", clamptouppervalue);
	}
	if (defaultvalue.size()) 
	{
		m->RegisterConfig("defaultvalue", defaultvalue);
	}
    return m;
}

std::shared_ptr<ProcessInterpPointDataToFld> ProcessInterpPointDataToFld_Init(FieldSharedPtr f,
										std::string frompts = "",
										std::string interpcoord = "")
{
	std::shared_ptr<ProcessInterpPointDataToFld> m = 
		MemoryManager<ProcessInterpPointDataToFld>::AllocateSharedPtr(f);
	
	if (frompts.size()) 
	{
		m->RegisterConfig("frompts", frompts);
	}
	if (interpcoord.size()) 
	{
		m->RegisterConfig("interpcoord", interpcoord);
	}
    return m;
}

std::shared_ptr<ProcessInterpPoints> ProcessInterpPoints_Init(FieldSharedPtr f,
									std::string fromxml = "",
									std::string fromfld = "",
									std::string topts = "",
									std::string line = "",
									std::string plane = "",
									std::string box = "",
									std::string clamptolowervalue = "",
									std::string clamptouppervalue = "",
									std::string defaultvalue = "",
									std::string cp = "")
{
	std::shared_ptr<ProcessInterpPoints> m = 
			MemoryManager<ProcessInterpPoints>::AllocateSharedPtr(f);
	
	if (fromxml.size()) 
	{
		m->RegisterConfig("fromxml", fromxml);
	}
	if (fromfld.size())
	{
		m->RegisterConfig("fromfld", fromfld);
	}
	if (topts.size()) 
	{
		m->RegisterConfig("topts", topts);
	}
	if (line.size()) 
	{
		m->RegisterConfig("line", line);
	}
	if (plane.size()) 
	{
		m->RegisterConfig("plane", plane);
	}
	if (box.size()) 
	{
		m->RegisterConfig("box", box);
	}
	if (clamptolowervalue.size()) 
	{
		m->RegisterConfig("clamptolowervalue", clamptolowervalue);
	}
	if (clamptouppervalue.size()) 
	{
		m->RegisterConfig("clamptouppervalue", clamptouppervalue);
	}
	if (defaultvalue.size()) 
	{
		m->RegisterConfig("defaultvalue", defaultvalue);
	}
	if (cp.size()) 
	{
		m->RegisterConfig("cp", cp);
	}
    return m;
}

std::shared_ptr<ProcessInterpPtsToPts> ProcessInterpPtsToPts_Init(FieldSharedPtr f,
									std::string topts = "",
									std::string line = "",
									std::string plane = "",
									std::string box = "",
									std::string clamptolowervalue = "",
									std::string clamptouppervalue = "",
									std::string defaultvalue = "",
									std::string cp = "")
{
	std::shared_ptr<ProcessInterpPtsToPts> m = 
			MemoryManager<ProcessInterpPtsToPts>::AllocateSharedPtr(f);
	
	if (topts.size()) 
	{
		m->RegisterConfig("topts", topts);
	}
	if (line.size()) 
	{
		m->RegisterConfig("line", line);
	}
	if (plane.size()) 
	{
		m->RegisterConfig("plane", plane);
	}
	if (box.size()) 
	{
		m->RegisterConfig("box", box);
	}
	if (clamptolowervalue.size()) 
	{
		m->RegisterConfig("clamptolowervalue", clamptolowervalue);
	}
	if (clamptouppervalue.size()) 
	{
		m->RegisterConfig("clamptouppervalue", clamptouppervalue);
	}
	if (defaultvalue.size()) 
	{
		m->RegisterConfig("defaultvalue", defaultvalue);
	}
	if (cp.size()) 
	{
		m->RegisterConfig("cp", cp);
	}
    return m;
}

std::shared_ptr<ProcessIsoContour> ProcessIsoContour_Init(FieldSharedPtr f,
									std::string fieldstr = "",
									std::string fieldname = "",
									std::string fieldid = "",
									std::string fieldvalue = "",
									bool globalcondense = false,
									bool smooth = false,
									std::string smoothiter = "",
									std::string smoothposdiffusion = "",
									std::string smoothnegdiffusion = "",
									std::string removesmallcontour = "")
{
	std::shared_ptr<ProcessIsoContour> m = 
			MemoryManager<ProcessIsoContour>::AllocateSharedPtr(f);
	
	if (fieldstr.size()) 
	{
		m->RegisterConfig("fieldstr", fieldstr);
	}
	if (fieldname.size()) 
	{
		m->RegisterConfig("fieldname", fieldname);
	}
	if (fieldid.size()) 
	{
		m->RegisterConfig("fieldid", fieldid);
	}
	if (fieldvalue.size()) 
	{
		m->RegisterConfig("fieldvalue", fieldvalue);
	}
	if (globalcondense) 
	{
		m->RegisterConfig("globalcondense");
	}
	if (smooth) 
	{
		m->RegisterConfig("smooth");
	}
	if (smoothiter.size()) 
	{
		m->RegisterConfig("smoothiter", smoothiter);
	}
	if (smoothposdiffusion.size()) 
	{
		m->RegisterConfig("smoothposdiffusion", smoothposdiffusion);
	}
	if (smoothnegdiffusion.size()) 
	{
		m->RegisterConfig("smoothnegdiffusion", smoothnegdiffusion);
	}
	if (removesmallcontour.size()) 
	{
		m->RegisterConfig("removesmallcontour", removesmallcontour);
	}
    return m;
}

std::shared_ptr<ProcessJacobianEnergy> ProcessJacobianEnergy_Init(FieldSharedPtr f,
										std::string topmodes = "")
{
	std::shared_ptr<ProcessJacobianEnergy> m = 
			MemoryManager<ProcessJacobianEnergy>::AllocateSharedPtr(f);
	
	if (topmodes.size())
	{
		m->RegisterConfig("topmodes", topmodes);
	}
    return m;
}

std::shared_ptr<ProcessMultiShear> ProcessMultiShear_Init(FieldSharedPtr f,
										std::string N = "",
										std::string fromfld = "")
{
	std::shared_ptr<ProcessMultiShear> m = 
			MemoryManager<ProcessMultiShear>::AllocateSharedPtr(f);
	
	if (N.size())
	{
		m->RegisterConfig("N", N);
	}
	if (fromfld.size())
	{
		m->RegisterConfig("fromfld", fromfld);
	}
    return m;
}

std::shared_ptr<ProcessPointDataToFld> ProcessPointDataToFld_Init(FieldSharedPtr f,
										std::string setnantovalue = "",
										std::string frompts = "")
{
	std::shared_ptr<ProcessPointDataToFld> m = 
			MemoryManager<ProcessPointDataToFld>::AllocateSharedPtr(f);
	
	if (setnantovalue.size())
	{
		m->RegisterConfig("setnantovalue", setnantovalue);
	}
	if (frompts.size())
	{
		m->RegisterConfig("frompts", frompts);
	}
    return m;
}

std::shared_ptr<ProcessQualityMetric> ProcessQualityMetric_Init(FieldSharedPtr f,
										bool scaled = false)
{
	std::shared_ptr<ProcessQualityMetric> m = 
			MemoryManager<ProcessQualityMetric>::AllocateSharedPtr(f);
	
	if (scaled)
	{
		m->RegisterConfig("scaled");
	}
    return m;
}

std::shared_ptr<ProcessRemoveField> ProcessRemoveField_Init(FieldSharedPtr f,
										std::string fieldname = "")
{
	std::shared_ptr<ProcessRemoveField> m = 
			MemoryManager<ProcessRemoveField>::AllocateSharedPtr(f);
	
	if (fieldname.size())
	{
		m->RegisterConfig("fieldname", fieldname);
	}
    return m;
}

std::shared_ptr<ProcessSurfDistance> ProcessSurfDistance_Init(FieldSharedPtr f,
										std::string bnd = "")
{
	std::shared_ptr<ProcessSurfDistance> m = 
			MemoryManager<ProcessSurfDistance>::AllocateSharedPtr(f);
	
	if (bnd.size())
	{
		m->RegisterConfig("bnd", bnd);
	}
    return m;
}

//std::shared_ptr<ProcessScaleInFld> ProcessScaleInFld_Init(FieldSharedPtr f,
										//std::string scale = "")
//{
	//std::shared_ptr<ProcessScaleInFld> m = 
			//MemoryManager<ProcessScaleInFld>::AllocateSharedPtr(f);
	
	//if (scale.size())
	//{
		//m->RegisterConfig("scale", scale);
	//}
    //return m;
//}

std::shared_ptr<ProcessWSS> ProcessWSS_Init(FieldSharedPtr f,
										std::string bnd = "",
										bool addnormals = false)
{
	std::shared_ptr<ProcessWSS> m = 
					MemoryManager<ProcessWSS>::AllocateSharedPtr(f);
	
	if (bnd.size()) 
	{
		m->RegisterConfig("bnd", bnd);
	}
	if (addnormals) 
	{
		m->RegisterConfig("addnormals");
	}
    return m;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(LoadFieldData_overloads, LoadFieldData, 0, 1);

void export_ProcessModules()
{
	py::class_<ProcessAddCompositeID, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessAddCompositeID> >
				("ProcessAddCompositeID", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessAddFld, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessAddFld> >
				("ProcessAddFld", py::no_init)
		.def("__init__", py::make_constructor(&ProcessAddFld_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("scale") = "", py::arg("fromfld") = ""))) 
	;
	py::class_<ProcessBoundaryExtract, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessBoundaryExtract> >
				("ProcessBoundaryExtract", py::no_init)
		.def("__init__", py::make_constructor(&ProcessBoundaryExtract_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("bnd") = "", py::arg("addnormals") = false))) 
	;
	py::class_<ProcessC0Projection, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessC0Projection> >
				("ProcessC0Projection", py::no_init)
		.def("__init__", py::make_constructor(&ProcessC0Projection_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("localtoglobalmap") = false, py::arg("usexmlbcs") = false,
			py::arg("helmsmoothing") = ""))) 
	;
	py::class_<ProcessCombineAvg, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessCombineAvg> >
				("ProcessCombineAvg", py::no_init)
		.def("__init__", py::make_constructor(&ProcessCombineAvg_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("fromfld") = ""))) 
	;
	py::class_<ProcessCreateExp, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessCreateExp> >
				("ProcessCreateExp", py::init<FieldSharedPtr>())
		.def("LoadFieldData", &ProcessCreateExp::LoadFieldData, LoadFieldData_overloads())
	;
	py::class_<ProcessDeform, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessDeform> >
				("ProcessDeform", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessDisplacement, py::bases<ProcessBoundaryExtract>, 
				std::shared_ptr<ProcessDisplacement> >
				("ProcessDisplacement", py::no_init)
		.def("__init__", py::make_constructor(&ProcessDisplacement_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("to") = "", py::arg("usevertexids") = false))) 
	;
	py::class_<ProcessDOF, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessDOF> >
				("ProcessDOF", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessEquiSpacedOutput, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessEquiSpacedOutput> >
				("ProcessEquiSpacedOutput", py::no_init)
		.def("__init__", py::make_constructor(&ProcessEquiSpacedOutput_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("tetonly") = false, py::arg("modalenergy") = false))) 
	;
	py::class_<ProcessFieldFromString, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessFieldFromString> >
				("ProcessFieldFromString", py::no_init)
		.def("__init__", py::make_constructor(&ProcessFieldFromString_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("fieldstr") = "", py::arg("fieldname") = ""))) 
	;
	py::class_<ProcessGrad, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessGrad> >
				("ProcessGrad", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessHomogeneousPlane, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessHomogeneousPlane> >
				("ProcessHomogeneousPlane", py::no_init)
		.def("__init__", py::make_constructor(&ProcessHomogeneousPlane_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("planeid") = "", py::arg("wavespace") = false))) 
	;
	py::class_<ProcessHomogeneousStretch, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessHomogeneousStretch> >
				("ProcessHomogeneousStretch", py::no_init)
		.def("__init__", py::make_constructor(&ProcessHomogeneousStretch_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("factor") = ""))) 
	;
	py::class_<ProcessInnerProduct, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessInnerProduct> >
				("ProcessInnerProduct", py::no_init)
		.def("__init__", py::make_constructor(&ProcessInnerProduct_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("fromfld") = "", py::arg("fields") = "", 
			py::arg("multifldids") = "", py::arg("allfromflds") = false))) 
	;
	py::class_<ProcessInterpField, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessInterpField> >
				("ProcessInterpField", py::no_init)
		.def("__init__", py::make_constructor(&ProcessInterpField_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("fromxml") = "", py::arg("fromfld") = "", 
			py::arg("clamptolowervalue") = "", py::arg("clamptouppervalue") = "",
			py::arg("defaultvalue") = ""))) 
	;
	py::class_<ProcessInterpPointDataToFld, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessInterpPointDataToFld> >
				("ProcessInterpPointDataToFld", py::no_init)
		.def("__init__", py::make_constructor(&ProcessInterpPointDataToFld_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("frompts") = "", py::arg("interpcoord") = "")))
	;
	py::class_<ProcessInterpPoints, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessInterpPoints> >
				("ProcessInterpPoints", py::no_init)
		.def("__init__", py::make_constructor(&ProcessInterpPoints_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("fromxml") = "", py::arg("fromfld") = "", 
			py::arg("topts") = "", py::arg("line") = "", 
			py::arg("plane") = "", py::arg("box") = "", 
			py::arg("clamptolowervalue") = "", py::arg("clamptouppervalue") = "",
			py::arg("defaultvalue") = "", py::arg("cp") = ""))) 
	;
	py::class_<ProcessInterpPtsToPts, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessInterpPtsToPts> >
				("ProcessInterpPtsToPts", py::no_init)
		.def("__init__", py::make_constructor(&ProcessInterpPtsToPts_Init, 
			py::default_call_policies(), (py::arg("f"),
			py::arg("topts") = "", py::arg("line") = "", 
			py::arg("plane") = "", py::arg("box") = "", 
			py::arg("clamptolowervalue") = "", py::arg("clamptouppervalue") = "",
			py::arg("defaultvalue") = "", py::arg("cp") = ""))) 
	;
	py::class_<ProcessIsoContour, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessIsoContour> >
				("ProcessIsoContour", py::no_init)
		.def("__init__", py::make_constructor(&ProcessIsoContour_Init, 
			py::default_call_policies(), (py::arg("f"),
			py::arg("fieldstr") = "", py::arg("fieldname") = "", 
			py::arg("fieldid") = "", py::arg("fieldvalue") = "", 
			py::arg("globalcondense") = false, py::arg("smooth") = false,
			py::arg("smoothiter") = "", py::arg("smoothposdiffusion") = "", 
			py::arg("smoothnegdiffusion") = "", py::arg("removesmallcontour") = ""))) 
	;
	py::class_<ProcessJacobianEnergy, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessJacobianEnergy> >
				("ProcessJacobianEnergy", py::no_init)
		.def("__init__", py::make_constructor(&ProcessJacobianEnergy_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("topmodes") = ""))) 
	;
	py::class_<ProcessL2Criterion, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessL2Criterion> >
				("ProcessL2Criterion", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessMapping, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessMapping> >
				("ProcessMapping", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessMean, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessMean> >
				("ProcessMean", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessMeanMode, py::bases<ProcessHomogeneousPlane>, 
				std::shared_ptr<ProcessMeanMode> >
				("ProcessMeanMode", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessMultiShear, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessMultiShear> >
				("ProcessMultiShear", py::no_init)
		.def("__init__", py::make_constructor(&ProcessMultiShear_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("N") = "", py::arg("fromfld") = ""))) 
	;
	py::class_<ProcessNumModes, py::bases<ProcessHomogeneousPlane>, 
				std::shared_ptr<ProcessNumModes> >
				("ProcessNumModes", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessPointDataToFld, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessPointDataToFld> >
				("ProcessPointDataToFld", py::no_init)
		.def("__init__", py::make_constructor(&ProcessPointDataToFld_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("setnantovalue") = "", py::arg("frompts") = ""))) 
	;
	py::class_<ProcessPrintFldNorms, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessPrintFldNorms> >
				("ProcessPrintFldNorms", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessQCriterion, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessQCriterion> >
				("ProcessQCriterion", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessQualityMetric, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessQualityMetric> >
				("ProcessQualityMetric", py::no_init)
		.def("__init__", py::make_constructor(&ProcessQualityMetric_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("scaled") = false))) 
	;
	py::class_<ProcessRemoveField, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessRemoveField> >
				("ProcessRemoveField", py::no_init)
		.def("__init__", py::make_constructor(&ProcessRemoveField_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("fieldname") = ""))) 
	;
	//py::class_<ProcessScaleInFld, py::bases<ProcessModule>, 
				//std::shared_ptr<ProcessScaleInFld> >
				//("ProcessScaleInFld", py::no_init)
		//.def("__init__", py::make_constructor(&ProcessScaleInFld_Init, 
			//py::default_call_policies(), (py::arg("f"), 
			//py::arg("scale") = ""))) 
	//;
	py::class_<ProcessScalGrad, py::bases<ProcessBoundaryExtract>, 
				std::shared_ptr<ProcessScalGrad> >
				("ProcessScalGrad", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessStreamFunction, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessStreamFunction> >
				("ProcessStreamFunction", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessSurfDistance, py::bases<ProcessBoundaryExtract>, 
				std::shared_ptr<ProcessSurfDistance> >
				("ProcessSurfDistance", py::no_init)
	.def("__init__", py::make_constructor(&ProcessSurfDistance_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("bnd") = ""))) 
	;
	py::class_<ProcessVorticity, py::bases<ProcessModule>, 
				std::shared_ptr<ProcessVorticity> >
				("ProcessVorticity", py::init<FieldSharedPtr>())
	;
	py::class_<ProcessWSS, py::bases<ProcessBoundaryExtract>, 
				std::shared_ptr<ProcessWSS> >
				("ProcessWSS", py::no_init)
		.def("__init__", py::make_constructor(&ProcessWSS_Init, 
			py::default_call_policies(), (py::arg("f"), 
			py::arg("bnd") = "", py::arg("addnormals") = false))) 
	;
}
