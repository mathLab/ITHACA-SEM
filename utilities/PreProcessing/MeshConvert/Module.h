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
//  License for the specific language governing rights and limitations under
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

#ifndef UTILITIES_PREPROCESSING_MESHCONVERT_MODULE
#define UTILITIES_PREPROCESSING_MESHCONVERT_MODULE

#include <string>
#include <iostream>
#include <LibUtilities/BasicUtils/NekFactory.hpp>

#include "MeshElements.h"

namespace Nektar
{
    namespace Utilities
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
            ConfigOption(bool isBool, string defValue, string desc) :
                isBool(isBool), beenSet(false), value(), defValue(defValue), 
                desc(desc) {}
            ConfigOption() :
                isBool(false), beenSet(false), value(), defValue(), desc() {}
            
            /**
             * @brief Re-interpret the value stored in #value as some type using
             * boost::lexical_cast.
             */
            template<typename T>
            T as()
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
            
            /// True if the configuration option is a boolean (thus does not
            /// need additional arguments).
            bool   isBool;
            /// True if the configuration option has been set at command
            /// line. If false, the default value will be put into #value.
            bool   beenSet;
            /// The value of the configuration option.
            string value;
            /// Default value of the configuration option.
            string defValue;
            /// Description of the configuration option.
            string desc;
        };
        
        /**
         * Abstract base class for mesh converter modules. Each subclass
         * implements the Process() function, which in some way alters the
         * mesh #m.
         */
        class Module
        {
        public:
            Module(MeshSharedPtr p_m) : m(p_m) {}
            virtual void Process() = 0;
            
            void RegisterConfig(string key, string value);
            void PrintConfig();
            void SetDefaults();
            
        protected:
            /// Mesh object
            MeshSharedPtr m;
            /// List of configuration values.
            map<string, ConfigOption> config;
            
            /// Extract element vertices
            virtual void ProcessVertices();
            /// Extract element edges
            virtual void ProcessEdges();
            /// Extract element faces
            virtual void ProcessFaces();
            /// Generate element IDs
            virtual void ProcessElements();
            /// Generate composites
            virtual void ProcessComposites();
        };
        
        /**
         * @brief Abstract base class for input modules.
         *
         * Input modules should read the contents of #mshFile in the Process()
         * function and populate the members of #m. Typically any given module
         * should populate Mesh::expDim, Mesh::spaceDim, Mesh::node and
         * Mesh::element, then call the protected ProcessX functions to
         * generate edges, faces, etc.
         */
        class InputModule : public Module
        {
        public:
            InputModule(MeshSharedPtr p_m);
            void OpenStream();
            
        protected:
            /// Print summary of elements.
            void         PrintSummary();
            /// Input stream
            std::ifstream mshFile;
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
            ProcessModule(MeshSharedPtr p_m) : Module(p_m) {}
        };
        
        /**
         * @brief Abstract base class for output modules.
         *
         * Output modules take the mesh #m and write to the file specified by
         * the stream #mshFile.
         */
        class OutputModule : public Module
        {
        public:
            OutputModule(MeshSharedPtr p_m);
            void OpenStream();
            
        protected:
            /// Output stream
            std::ofstream mshFile;
        };
        
        typedef std::pair<ModuleType,std::string> ModuleKey;
        std::ostream& operator<<(std::ostream& os, const ModuleKey& rhs);

        typedef boost::shared_ptr<Module> ModuleSharedPtr;
        typedef LibUtilities::NekFactory< ModuleKey, Module, MeshSharedPtr > ModuleFactory;
        
        ModuleFactory& GetModuleFactory();
    }
}

#endif
