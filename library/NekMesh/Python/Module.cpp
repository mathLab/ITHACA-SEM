///////////////////////////////////////////////////////////////////////////////
//
// File: Module.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Python wrapper for Module.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <NekMesh/Module/Module.h>

#include <boost/algorithm/string.hpp>
#include <boost/python/raw_function.hpp>

using namespace Nektar;
using namespace Nektar::NekMesh;

/**
 * @brief Module wrapper to handle virtual function calls in @c Module and its
 * subclasses as defined by the template parameter @tparam MODTYPE.
 */
template<class MODTYPE>
struct ModuleWrap : public MODTYPE, public py::wrapper<MODTYPE>
{
    /**
     * @brief Constructor, which is identical to NekMesh::Module::Module.
     *
     * @param mesh  Input mesh.
     */
    ModuleWrap(MeshSharedPtr mesh) : MODTYPE(mesh), py::wrapper<MODTYPE>()
    {
    }

    /**
     * @brief Concrete implementation of the Module::Process function.
     */
    void Process()
    {
        this->get_override("Process")();
    }

    /**
     * @brief Defines a configuration option for this module.
     *
     * @param key     The name of the configuration option.
     * @param def     The option's default value.
     * @param desc    A text description of the option.
     * @param isBool  If true, this option is a boolean-type (true/false).
     */
    void AddConfigOption(std::string key, std::string def, std::string desc,
                         bool isBool)
    {
        ConfigOption conf(isBool, def, desc);
        this->m_config[key] = conf;
    }

    // We expose Module::m_mesh as a public member variable so that we can
    // adjust this using Python attributes.
    using MODTYPE::m_mesh;
};

template<typename T>
T Module_GetConfig(std::shared_ptr<Module> mod,
                   const std::string &key)
{
    return mod->GetConfigOption(key).as<T>();
}

/**
 * @brief Lightweight wrapper for NekMesh::Module::ProcessEdges.
 *
 * @param mod        Module to call.
 * @param reprocess  If true then edges will be reprocessed (i.e. match 1D
 *                   elements to 2D element edges to construct boundaries).
 */
void Module_ProcessEdges(std::shared_ptr<Module> mod,
                         bool reprocess = true)
{
    mod->ProcessEdges(reprocess);
}

/**
 * @brief Lightweight wrapper for NekMesh::Module::ProcessFaces.
 *
 * @param mod        Module to call.
 * @param reprocess  If true then faces will be reprocessed (i.e. match 2D
 *                   elements to 3D element faces to construct boundaries).
 */
void Module_ProcessFaces(std::shared_ptr<Module> mod,
                         bool reprocess = true)
{
    mod->ProcessFaces(reprocess);
}

template<typename MODTYPE>
struct ModuleTypeProxy
{
};

template<>
struct ModuleTypeProxy<InputModule>
{
    static const ModuleType value = eInputModule;
};

template<>
struct ModuleTypeProxy<ProcessModule>
{
    static const ModuleType value = eProcessModule;
};

template<>
struct ModuleTypeProxy<OutputModule>
{
    static const ModuleType value = eOutputModule;
};

const ModuleType ModuleTypeProxy<InputModule>::value;
const ModuleType ModuleTypeProxy<ProcessModule>::value;
const ModuleType ModuleTypeProxy<OutputModule>::value;

/**
 * @brief Lightweight wrapper for Module factory creation function.
 *
 * @param  modType  Module type (input/process/output).
 * @param  modName  Module name (typically filename extension).
 * @param  mesh     Mesh that will be passed between modules.
 * @tparam MODTYPE  Subclass of Module (e.g #InputModule, #OutputModule)
 */
template<typename MODTYPE>
ModuleSharedPtr Module_Create(py::tuple args, py::dict kwargs)
{
    std::string modName = py::extract<std::string>(args[0]);
    ModuleKey modKey = std::make_pair(ModuleTypeProxy<MODTYPE>::value, modName);
    MeshSharedPtr mesh = py::extract<MeshSharedPtr>(args[1]);
    ModuleSharedPtr mod = GetModuleFactory().CreateInstance(modKey, mesh);

    // Process keyword arguments.
    py::list items = kwargs.items();

    for (int i = 0; i < py::len(items); ++i)
    {
        std::string arg = py::extract<std::string>(items[i][0]);
        std::string val = py::extract<std::string>(items[i][1].attr("__str__")());
        mod->RegisterConfig(arg, val);
    }

    mod->SetDefaults();

    return mod;
}

/**
 * @brief Lightweight wrapper for NekMesh::Module::RegisterConfig.
 *
 * @param mod    Module to call
 * @param key    Configuration key.
 * @param value  Optional value (some configuration options are boolean).
 */
void Module_RegisterConfig(std::shared_ptr<Module> mod,
                           std::string const &key,
                           std::string const &value)
{
    mod->RegisterConfig(key, value);
}

template<typename MODTYPE>
void ModuleWrap_AddConfigOption(std::shared_ptr<ModuleWrap<MODTYPE>> mod,
                                std::string const &key,
                                std::string const &defValue,
                                std::string const &desc,
                                bool isBool)
{
    mod->AddConfigOption(key, defValue, desc, isBool);
}

class ModuleRegisterHelper
{
public:
    ModuleRegisterHelper(py::object obj) : m_obj(obj)
    {
        py::incref(obj.ptr());
    }

    ~ModuleRegisterHelper()
    {
        py::decref(m_obj.ptr());
    }

    ModuleSharedPtr create(MeshSharedPtr mesh)
    {
        py::object inst = m_obj(mesh);
        return py::extract<ModuleSharedPtr>(inst);
    }

protected:
    py::object m_obj;
};

#if PY_MAJOR_VERSION == 2
void ModuleCapsuleDestructor(void *ptr)
{
    ModuleRegisterHelper *tmp = (ModuleRegisterHelper *)ptr;
    delete tmp;
}
#else
void ModuleCapsuleDestructor(PyObject *ptr)
{
    ModuleRegisterHelper *tmp =
        (ModuleRegisterHelper *)PyCapsule_GetPointer(ptr, 0);
    delete tmp;
}
#endif

/**
 * @brief Lightweight wrapper for the Module factory RegisterCreatorFunction, to
 * support the ability for Python subclasses of Module to register themselves
 * with the Nektar++ Module factory.
 *
 * This function wraps the NekFactory RegisterCreatorFunction. This function
 * expects a function pointer to a C++ object that will construct a Module. In
 * this case we therefore need to construct a function call that will construct
 * our Python object (which is a subclass of Module), and then pass this back to
 * Boost.Python to give the Python object back.
 *
 * We have to do some indirection here to get this to work, but we can
 * achieve this with the following strategy:
 *
 * - Create a @c ModuleRegisterHelper object, which as an argument will store
 *   the Python class instance that will be instantiated from the Python side.
 * - Using std::bind, construct a function pointer to the helper's creation
 *   function, ModuleRegisterHelper::create.
 * - Create a Python capsule that will contain the @c ModuleRegisterHelper
 *   instance, and register this in the global namespace of the current
 *   module. This then ties the capsule to the lifetime of the module.
 */
void Module_Register(ModuleType const  &modType,
                     std::string const &modName,
                     py::object &obj)
{
    // Create a module register helper, which will call the C++ function to
    // create the module.
    ModuleRegisterHelper *helper = new ModuleRegisterHelper(obj);

    // Register this with the module factory using std::bind to grab a function
    // pointer to that particular object's function.
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, modName),
        std::bind(
            &ModuleRegisterHelper::create, helper, std::placeholders::_1));

    // Create a capsule that will be embedded in the __main__ namespace. So
    // deallocation will occur, but only once Python ends or the Python module
    // is deallocated.
    std::string modkey =
        "_" + std::string(ModuleTypeMap[modType]) + "_" + modName;

#if PY_MAJOR_VERSION == 2
    py::object capsule(
        py::handle<>(PyCObject_FromVoidPtr(helper, ModuleCapsuleDestructor)));
#else
    py::object capsule(
        py::handle<>(PyCapsule_New(helper, 0, ModuleCapsuleDestructor)));
#endif

    // Embed this in __main__.
    py::import("__main__").attr(modkey.c_str()) = capsule;
}

template <typename MODTYPE>
struct ModuleWrapConverter
{
    ModuleWrapConverter()
    {
        // An important bit of code which will register allow
        // shared_ptr<MODTYPE> as something that boost::python recognises,
        // otherwise modules constructed from the factory will not work from
        // Python.
        py::objects::class_value_wrapper<
            std::shared_ptr<MODTYPE>,
            py::objects::make_ptr_instance<
                MODTYPE,
                py::objects::pointer_holder<std::shared_ptr<MODTYPE>, MODTYPE>>
            >();
    }
};

template <typename MODTYPE>
struct PythonModuleClass
{
    PythonModuleClass(std::string modName)
    {
        py::class_<ModuleWrap<MODTYPE>,
                   std::shared_ptr<ModuleWrap<MODTYPE>>,
                   py::bases<Module>,
                   boost::noncopyable>(
                       modName.c_str(), py::init<MeshSharedPtr>())

            .def("AddConfigOption", ModuleWrap_AddConfigOption<MODTYPE>, (
                     py::arg("key"), py::arg("defValue"), py::arg("desc"),
                     py::arg("isBool") = false))

            // Allow direct access to mesh object through a property.
            .def_readwrite("mesh", &ModuleWrap<MODTYPE>::m_mesh)

            .def("Process", py::pure_virtual(&MODTYPE::Process))
            .def("Create", py::raw_function(Module_Create<MODTYPE>))
            .staticmethod("Create")
            ;

        ModuleWrapConverter<MODTYPE>();
    }
};

void export_Module()
{
    // Export ModuleType enum.
    NEKPY_WRAP_ENUM(ModuleType, ModuleTypeMap);

    // Define ModuleWrap to be implicitly convertible to a Module, since it
    // seems that doesn't sometimes get picked up.
    py::implicitly_convertible<std::shared_ptr<ModuleWrap<Module>>,
                               std::shared_ptr<Module>>();
    py::implicitly_convertible<std::shared_ptr<ModuleWrap<InputModule>>,
                               std::shared_ptr<Module>>();
    py::implicitly_convertible<std::shared_ptr<ModuleWrap<OutputModule>>,
                               std::shared_ptr<Module>>();
    py::implicitly_convertible<std::shared_ptr<ModuleWrap<ProcessModule>>,
                               std::shared_ptr<Module>>();

    // Wrapper for the Module class. Note that since Module contains a pure
    // virtual function, we need the ModuleWrap helper class to handle this for
    // us. In the lightweight wrappers above, we therefore need to ensure we're
    // passing std::shared_ptr<Module> as the first argument, otherwise they
    // won't accept objects constructed from Python.
    py::class_<ModuleWrap<Module>,
               std::shared_ptr<ModuleWrap<Module>>,
               boost::noncopyable>(
                   "Module", py::init<MeshSharedPtr>())

        // Process function for this module.
        .def("Process", py::pure_virtual(&Module::Process))

        // Configuration options.
        .def("RegisterConfig", Module_RegisterConfig, (
                 py::arg("key"), py::arg("value") = ""))
        .def("PrintConfig", &Module::PrintConfig)
        .def("SetDefaults", &Module::SetDefaults)
        .def("GetStringConfig", Module_GetConfig<std::string>)
        .def("GetFloatConfig", Module_GetConfig<double>)
        .def("GetIntConfig", Module_GetConfig<int>)
        .def("GetBoolConfig", Module_GetConfig<bool>)
        .def("AddConfigOption", ModuleWrap_AddConfigOption<Module>, (
                 py::arg("key"), py::arg("defValue"), py::arg("desc"),
                 py::arg("isBool") = false))

        // Mesh accessor method.
        .def("GetMesh", &Module::GetMesh)

        // Mesh processing functions.
        .def("ProcessVertices", &Module::ProcessVertices)
        .def("ProcessEdges", Module_ProcessEdges, (
                 py::arg("reprocessEdges") = true))
        .def("ProcessFaces", Module_ProcessFaces, (
                 py::arg("reprocessFaces") = true))
        .def("ProcessElements", &Module::ProcessElements)
        .def("ProcessComposites", &Module::ProcessComposites)
        .def("ClearElementLinks", &Module::ClearElementLinks)

        // Allow direct access to mesh object through a property.
        .def_readwrite("mesh", &ModuleWrap<Module>::m_mesh)

        // Factory functions.
        .def("Register", &Module_Register)
        .staticmethod("Register")
        ;

    ModuleWrapConverter<Module>();

    PythonModuleClass<InputModule>("InputModule");
    PythonModuleClass<ProcessModule>("ProcessModule");
    PythonModuleClass<OutputModule>("OutputModule");
}
