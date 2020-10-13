////////////////////////////////////////////////////////////////////////////////
//
//  File: Module.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Mesh converter module base classes.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESH_MODULE
#define NEKMESH_MODULE

#include <NekMesh/Module/Log.hpp>

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <set>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <NekMesh/MeshElements/Mesh.h>
#include <NekMesh/NekMeshDeclspec.h>

namespace io = boost::iostreams;

namespace Nektar
{
namespace NekMesh
{
    /**
     * Denotes different types of mesh converter modules: so far only
     * input, output and process modules are defined.
     */
    enum ModuleType {
        eInputModule,
        eProcessModule,
        eOutputModule,
        SIZE_ModuleType
    };

    const char* const ModuleTypeMap[] =
    {
        "Input",
        "Process",
        "Output"
    };

    typedef std::map<int, std::pair<FaceSharedPtr, std::vector<int> > > PerMap;

    /**
     * @brief Represents a command-line configuration option.
     */
    struct ConfigOption
    {
        /**
         * @brief Construct a new configuration option.
         *
         * @param isBool    True if the option is boolean type.
         * @param defValue  Default value of the option.
         * @param desc      Description of the option.
         */
        ConfigOption(bool isBool, std::string defValue, std::string desc) :
            isBool(isBool), beenSet(false), value(), defValue(defValue),
            desc(desc) {}
        ConfigOption() :
            isBool(false), beenSet(false), value(), defValue(), desc() {}

        /**
         * @brief Re-interpret the value stored in #value as some type using
         * boost::lexical_cast.
         */
        template<typename T>
        T as() const
        {
            try
            {
                return boost::lexical_cast<T>(value);
            }
            catch(const std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                abort();
            }
        }

        /**
         * @brief Interpret the value stored in #value as some type using
         * boost::lexical_cast and return true of false depending on cast
         */
        template<typename T>
        bool isType() const
        {
            bool returnval = true;
            try
            {
                boost::lexical_cast<T>(value);
            }
            catch(const std::exception &e)
            {
                returnval = false;
            }

            return returnval;
        }


        /// True if the configuration option is a boolean (thus does not
        /// need additional arguments).
        bool   isBool;
        /// True if the configuration option has been set at command
        /// line. If false, the default value will be put into #value.
        bool   beenSet;
        /// The value of the configuration option.
        std::string value;
        /// Default value of the configuration option.
        std::string defValue;
        /// Description of the configuration option.
        std::string desc;
    };


    /**
     * Abstract base class for mesh converter modules. Each subclass
     * implements the Process() function, which in some way alters the
     * mesh #m.
     */
    class Module
    {
    public:
        NEKMESH_EXPORT Module(MeshSharedPtr p_m) : m_mesh(p_m), m_log()
        {
        }

        NEKMESH_EXPORT virtual void Process() = 0;
        NEKMESH_EXPORT virtual ~Module() = default;

        NEKMESH_EXPORT void SetLogger(Logger &log)
        {
            m_log = log;
            m_log.SetPrefix(GetModuleName());
        }

        NEKMESH_EXPORT Logger &GetLogger()
        {
            return m_log;
        }

        NEKMESH_EXPORT void RegisterConfig(std::string key,
                                           std::string value = std::string());
        NEKMESH_EXPORT void PrintConfig();
        NEKMESH_EXPORT void SetDefaults();

        MeshSharedPtr GetMesh()
        {
            return m_mesh;
        }

        const ConfigOption &GetConfigOption(std::string key) const
        {
            auto it = m_config.find(key);
            ASSERTL0(it != m_config.end(), "Configuration key not found!");
            return it->second;
        }

        virtual std::string GetModuleName()
        {
            return "";
        }

        /// Extract element vertices
        NEKMESH_EXPORT void ProcessVertices();
        /// Extract element edges
        NEKMESH_EXPORT void ProcessEdges(bool ReprocessEdges = true);
        /// Extract element faces
        NEKMESH_EXPORT void ProcessFaces(bool ReprocessFaces = true);
        /// Generate element IDs
        NEKMESH_EXPORT void ProcessElements();
        /// Generate composites
        NEKMESH_EXPORT void ProcessComposites();

        NEKMESH_EXPORT void ClearElementLinks();

    protected:
        /// Mesh object
        MeshSharedPtr m_mesh;
        /// List of configuration values.
        std::map<std::string, ConfigOption> m_config;
        /// Logger object.
        Logger m_log;

        NEKMESH_EXPORT void ReorderPrisms(PerMap &perFaces);
        NEKMESH_EXPORT void PrismLines(int                     prism,
                                       PerMap                  &perFaces,
                                       std::set<int>           &prismsDone,
                                       std::vector<ElementSharedPtr> &line);
    };

    /**
     * @brief Abstract base class for input modules.
     *
     * Input modules should read the contents of #m_mshFile in the Process()
     * function and populate the members of #m. Typically any given module
     * should populate Mesh::expDim, Mesh::spaceDim, Mesh::node and
     * Mesh::element, then call the protected ProcessX functions to
     * generate edges, faces, etc.
     */
    class InputModule : public Module
    {
    public:
        NEKMESH_EXPORT InputModule(MeshSharedPtr p_m);
        NEKMESH_EXPORT void OpenStream();

    protected:
        /// Print summary of elements.
        NEKMESH_EXPORT void PrintSummary();
        /// Input stream
        io::filtering_istream m_mshFile;
        /// Input stream
        std::ifstream m_mshFileStream;
    };

    /**
     * @brief Abstract base class for processing modules.
     *
     * Processing modules take a populated %Mesh object and process it in
     * some fashion; for example the %ProcessJac module calculates the
     * Jacobian of each element and prints warnings for non-positive
     * elements.
     */
    class ProcessModule : public Module
    {
    public:
        NEKMESH_EXPORT ProcessModule(MeshSharedPtr p_m) : Module(p_m) {}
    };

    /**
     * @brief Abstract base class for output modules.
     *
     * Output modules take the mesh #m and write to the file specified by
     * the stream #m_mshFile.
     */
    class OutputModule : public Module
    {
    public:
        NEKMESH_EXPORT OutputModule(MeshSharedPtr p_m);
        NEKMESH_EXPORT void OpenStream();


    protected:
        /// Output stream
        io::filtering_ostream m_mshFile;
        /// Input stream
        std::ofstream m_mshFileStream;
    };

    typedef std::pair<ModuleType,std::string> ModuleKey;
    NEKMESH_EXPORT std::ostream& operator<<(std::ostream& os, const ModuleKey& rhs);

    typedef std::shared_ptr<Module> ModuleSharedPtr;
    typedef LibUtilities::NekFactory< ModuleKey, Module, MeshSharedPtr > ModuleFactory;

    NEKMESH_EXPORT ModuleFactory& GetModuleFactory();
}
}

#endif
