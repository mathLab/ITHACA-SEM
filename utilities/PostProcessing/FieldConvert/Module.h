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
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Communication/CommSerial.h>

#include "Field.hpp"

namespace po = boost::program_options;

namespace Nektar
{
    namespace Utilities
    {
        using namespace std;

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
            m_isBool(isBool), m_beenSet(false), m_value(), 
                m_defValue(defValue), m_desc(desc) {}
        ConfigOption() :
            m_isBool(false), m_beenSet(false), m_value(), 
                m_defValue(), m_desc() {}
            
            /**
             * @brief Re-interpret the value stored in #value as some type using
             * boost::lexical_cast.
             */
            template<typename T>
            T as()
            {
                try
                {
                    return boost::lexical_cast<T>(m_value);
                }
                catch(const exception &e)
                {
                    cerr << e.what() << endl;
                    abort();
                }
            }
            
            /// True if the configuration option is a boolean (thus does not
            /// need additional arguments).
            bool   m_isBool;
            /// True if the configuration option has been set at command
            /// line. If false, the default value will be put into #value.
            bool   m_beenSet;
            /// The value of the configuration option.
            string m_value;
            /// Default value of the configuration option.
            string m_defValue;
            /// Description of the configuration option.
            string m_desc;
        };
        
        /**
         * Abstract base class for mesh converter modules. Each subclass
         * implements the Process() function, which in some way alters the
         * mesh #m.
         */
        class Module
        {
        public:
            Module(FieldSharedPtr p_f) : m_f(p_f), m_requireEquiSpaced(false) {}
            virtual void Process(po::variables_map &vm) = 0;
            
            void RegisterConfig(string key, string value);
            void PrintConfig();
            void SetDefaults();
            
            bool GetRequireEquiSpaced(void)
            {
                return m_requireEquiSpaced;
            }

            void SetRequireEquiSpaced(bool pVal)
            {
                m_requireEquiSpaced = pVal;
            }

        protected:
            /// Field object
            FieldSharedPtr m_f;
            /// List of configuration values.
            map<string, ConfigOption> m_config;
            bool m_requireEquiSpaced;
            
        };
        
        /**
         * @brief Abstract base class for input modules.
         *
         * Input modules should read the contents of #fldFile in the Process()
         * function and populate the members of #m. Typically any given module
         * should populate Mesh::expDim, Mesh::spaceDim, Mesh::node and
         * Mesh::element, then call the protected ProcessX functions to
         * generate edges, faces, etc.
         */
        class InputModule : public Module
        {
        public:
            InputModule(FieldSharedPtr p_m);
            void AddFile(string fileType, string fileName);
            
        protected:
            /// Print summary of elements.
            void PrintSummary();
            set<string> m_allowedFiles;

        };
        
        typedef boost::shared_ptr<InputModule> InputModuleSharedPtr;

        /**
         * @brief Abstract base class for processing modules.
         * 
         */
        class ProcessModule : public Module
        {
        public:
            ProcessModule(FieldSharedPtr p_f) : Module(p_f) {}

        };
        
        /**
         * @brief Abstract base class for output modules.
         *
         * Output modules take the mesh #m and write to the file specified by
         * the stream.
         */
        class OutputModule : public Module
        {
        public:
            OutputModule(FieldSharedPtr p_f);
            void OpenStream();
            
        protected:
            /// Output stream
            ofstream m_fldFile;
        
        };
        
        typedef pair<ModuleType,string> ModuleKey;
        ostream& operator<<(ostream& os, const ModuleKey& rhs);
        
        typedef boost::shared_ptr<Module> ModuleSharedPtr;
        typedef LibUtilities::NekFactory<ModuleKey, Module, FieldSharedPtr>
            ModuleFactory;
        
        ModuleFactory& GetModuleFactory();
        
        class FieldConvertComm : public  LibUtilities::CommSerial
        {
        public:
            FieldConvertComm(int argc, char* argv[], int size, int rank) : CommSerial(argc, argv)
            {
                m_size = size;
                m_rank = rank;
                m_type = "FieldConvert parallel";
            }
            FieldConvertComm(int size, int rank) : CommSerial(0, NULL)
            {
                m_size = size;
                m_rank = rank;
                m_type = "FieldConvert parallel";
            }
            virtual ~FieldConvertComm() {}
            void v_SplitComm(int pRows, int pColumns)
            {
            // Compute row and column in grid.
                m_commRow    = boost::shared_ptr<FieldConvertComm>(new FieldConvertComm(pColumns,m_rank));
                m_commColumn = boost::shared_ptr<FieldConvertComm>(new FieldConvertComm(pRows,0));
            }

        protected:
            int v_GetRank(void)
            {
                return m_rank;
            }

            bool v_TreatAsRankZero(void)
            {
                return true;
            }
            
            bool v_RemoveExistingFiles(void)
            {
                return false;
            }
        private:
            int m_rank;
        };
    }
}

#endif
