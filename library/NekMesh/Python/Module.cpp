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
 * @brief A basic Python logger object, which writes log messages to sys.stdout.
 *
 * This log streamer writes log messages to Python's sys.stdout. Although this
 * may seem superfluous, in some cases certain Python applications (in
 * particular, Jupyter notebooks), override sys.stdout to redirect output to
 * e.g. a browser or other IOStream. This therefore enables C++ log messages to
 * be redirected accordingly.
 */
class PythonStream : public LogOutput
{
public:
    /// Default constructor.
    PythonStream() : LogOutput()
    {
    }

    /**
     * @brief Write a log message @p msg.
     *
     * The logger runs the Python command `sys.stdout.write(msg)`.
     *
     * @param msg   The message to be printed.
     */
    void Log(const std::string &msg) override
    {
        // Create a Python string with the message.
        py::str pyMsg(msg);
        // Write to the Python stdout.
        py::import("sys").attr("stdout").attr("write")(msg);
    }

private:
    /// Finalise function. Nothing required in this case.
    void Finalise() override
    {
    }
};

/// A horrible but necessary static variable to enable/disable verbose output by
/// default.
static bool default_verbose = false;

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
    ModuleType modType = ModuleTypeProxy<MODTYPE>::value;

    if (modType == eProcessModule && py::len(args) != 2)
    {
        throw NekMeshError("ProcessModule.Create() requires two arguments: "
                           "module name and a Mesh object.");
    }
    else if (modType != eProcessModule && py::len(args) != 3)
    {
        throw NekMeshError(ModuleTypeMap[modType] + "Module.Create() requires "
                           "three arguments: module name, a Mesh object, and a "
                           "filename");
    }

    std::string modName = py::extract<std::string>(args[0]);
    ModuleKey modKey = std::make_pair(modType, modName);

    if (!py::extract<MeshSharedPtr>(args[1]).check())
    {
        throw NekMeshError("Second argument to Create() should be a mesh "
                           "object.");
    }

    MeshSharedPtr mesh = py::extract<MeshSharedPtr>(args[1]);

    if (!GetModuleFactory().ModuleExists(modKey))
    {
        throw NekMeshError("Module '" + modName + "' does not exist.");
    }

    ModuleSharedPtr mod = GetModuleFactory().CreateInstance(modKey, mesh);

    // First argument for input/output module should be the filename.
    if (modKey.first == eInputModule)
    {
        mod->RegisterConfig("infile", py::extract<std::string>(args[2]));
    }
    else if (modKey.first == eOutputModule)
    {
        mod->RegisterConfig("outfile", py::extract<std::string>(args[2]));
    }

    // Process keyword arguments.
    py::list items = kwargs.items();

    // Set default verbosity.
    bool verbose = default_verbose;

    for (int i = 0; i < py::len(items); ++i)
    {
        std::string arg = py::extract<std::string>(items[i][0]), val;

        // Enable or disable verbose for this module accordingly.
        if (arg == "verbose")
        {
            verbose = py::extract<bool>(items[i][1]);
            continue;
        }

        val = py::extract<std::string>(items[i][1].attr("__str__")());
        mod->RegisterConfig(arg, val);
    }

    // Set a logger for this module.
    auto pythonLog = std::make_shared<PythonStream>();
    Logger log = Logger(pythonLog, verbose ? VERBOSE : INFO);
    mod->SetLogger(log);

    // Set other default arguments.
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

/**
 * @brief Add a configuration option for the module.
 *
 * @tparam MODTYPE   Module type (e.g. #ProcessModule).
 * @param  mod       Module object.
 * @param  key       Name of the configuration option.
 * @param  defValue  Default value.
 * @param  desc      Description of the option.
 * @param  isBool    If true, denotes that the option will be bool-type.
 */
template<typename MODTYPE>
void ModuleWrap_AddConfigOption(std::shared_ptr<ModuleWrap<MODTYPE>> mod,
                                std::string const &key,
                                std::string const &defValue,
                                std::string const &desc,
                                bool isBool)
{
    mod->AddConfigOption(key, defValue, desc, isBool);
}

/**
 * @brief Helper class to handle module registration.
 *
 * This class is used in combination with the #Module_Register function to
 * handle module registration. In particular, given a Python object that
 * represents the user-supplied class (which inherits from Module), we use the
 * functions in here to handle reference management, as well as the creation of
 * the module that can then be passed back to Python.
 */
class ModuleRegisterHelper
{
public:
    /**
     * @brief Constructor.
     *
     * @param obj  Python class object that inherits from Module.
     */
    ModuleRegisterHelper(py::object obj) : m_obj(obj)
    {
        py::incref(obj.ptr());
    }

    /**
     * @brief Destructor.
     */
    ~ModuleRegisterHelper()
    {
        py::decref(m_obj.ptr());
    }

    /**
     * @brief Constructs a module from a given mesh, using the Python object
     * stored in #m_obj.
     */
    ModuleSharedPtr create(MeshSharedPtr mesh)
    {
        py::object inst = m_obj(mesh);
        return py::extract<ModuleSharedPtr>(inst);
    }

protected:
    /// Python object that represents a subclass of Module.
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
        ModuleKey(modType, modName),
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

/**
 * @brief Enables or disables verbose output by default.
 */
void Module_Verbose(bool verbose)
{
    default_verbose = verbose;
}

template <typename MODTYPE>
struct ModuleWrapConverter
{
    ModuleWrapConverter()
    {
        // An important bit of code which will allow shared_ptr<MODTYPE> to be
        // interpreted as something that boost::python recognises, otherwise
        // modules constructed from the factory will not work from Python.
        py::objects::class_value_wrapper<
            std::shared_ptr<MODTYPE>,
            py::objects::make_ptr_instance<
                MODTYPE,
                py::objects::pointer_holder<std::shared_ptr<MODTYPE>, MODTYPE>>
            >();
    }
};

/**
 * @brief Wrapper for subclasses of the Module class, e.g. #InputModule, which
 * can then be inhereted from inside Python.
 */
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
    NEKPY_WRAP_ENUM_STRING(ModuleType, ModuleTypeMap);

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

        // Enable verbose output (or not).
        .def("Verbose", &Module_Verbose)
        .staticmethod("Verbose")
        ;

    ModuleWrapConverter<Module>();

    // Wrap the Module subclasses.
    PythonModuleClass<InputModule>("InputModule");
    PythonModuleClass<ProcessModule>("ProcessModule");
    PythonModuleClass<OutputModule>("OutputModule");
}
