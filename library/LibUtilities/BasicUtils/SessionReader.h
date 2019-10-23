///////////////////////////////////////////////////////////////////////////////
//
// File SessionReader.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_SESSIONREADER_H
#define NEKTAR_LIB_UTILITIES_SESSIONREADER_H

#include <map>
#include <string>
#include <memory>

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/program_options/variables_map.hpp>

class TiXmlElement;
class TiXmlDocument;

namespace Nektar
{
    namespace LibUtilities
    {
        typedef std::map<std::string, std::string>  SolverInfoMap;
        typedef std::map<std::string, NekDouble>    ParameterMap;
        typedef std::map<std::string, std::string>  GeometricInfoMap;
        typedef std::map<std::string, std::string>  ExpressionMap;
        typedef std::vector<std::string>            VariableList;
        typedef std::map<std::string, std::string>  TagMap;
        typedef std::map<std::string, std::string>  FilterParams;
        typedef std::vector<
            std::pair<std::string, FilterParams> >   FilterMap;

        struct CmdLineArg
        {
            std::string shortName;
            std::string description;
            bool        isFlag;
        };

        typedef std::map<std::string, CmdLineArg>    CmdLineArgMap;

        typedef std::map<std::string, int>           EnumMap;
        typedef std::map<std::string, EnumMap>       EnumMapList;

        typedef std::map<std::string, std::string>   GloSysInfoMap;
        typedef std::map<std::string, GloSysInfoMap> GloSysSolnInfoList;

        typedef std::map<std::string, std::string>   GloSysInfoMap;
        typedef std::map<std::string, GloSysInfoMap> GloSysSolnInfoList;

        enum FunctionType
        {
            eFunctionTypeExpression,
            eFunctionTypeFile,
            eFunctionTypeTransientFile,
            eSIZE_FunctionType,
        };
        const char* const FunctionTypeMap[] =
        {
            "No Function type",
            "Expression",
            "File"
        };

        struct FunctionVariableDefinition
        {
            enum FunctionType m_type;
            std::string       m_filename;
            EquationSharedPtr m_expression;
            std::string       m_fileVariable;
        };

        typedef std::map<std::pair<std::string,int>, FunctionVariableDefinition>
            FunctionVariableMap;
        typedef std::map<std::string, FunctionVariableMap >
            FunctionMap;

        class SessionReader;
        typedef std::shared_ptr<SessionReader> SessionReaderSharedPtr;

        /// Reads and parses information from a Nektar++ XML session file.
        class SessionReader
        {
        public:
            /// Support creation through MemoryManager.
            friend class MemoryManager<SessionReader>;

            /**
             * @brief Creates an instance of the SessionReader class.
             *
             * This function should be used by an application to instantiate the
             * session reader. It should be called at the very beginning of the
             * application before any other processing of command-line
             * arguments. After instantiating the class and setting up any
             * parallel communication, it also calls the main initialisation
             * of the object.
             */
            LIB_UTILITIES_EXPORT static SessionReaderSharedPtr CreateInstance(
                int argc, char *argv[])
            {
                SessionReaderSharedPtr p = MemoryManager<
                    LibUtilities::SessionReader>::AllocateSharedPtr(argc, argv);
                return p;
            }

            /**
             * @brief Creates an instance of the SessionReader class initialised
             *        using a separate list of XML documents.
             *
             * This function should be used by an application to instantiate the
             * session reader. It may be called after processing of command-line
             * arguments. After instantiating the class and setting up any
             * parallel communication, it also calls the main initialisation
             * of the object.
             */
            LIB_UTILITIES_EXPORT static SessionReaderSharedPtr CreateInstance(
                int                       argc,
                char                     *argv[],
                std::vector<std::string> &pFilenames,
                const CommSharedPtr      &pComm = CommSharedPtr())
            {
                SessionReaderSharedPtr p = MemoryManager<
                    LibUtilities::SessionReader>
                        ::AllocateSharedPtr(argc, argv, pFilenames, pComm);
                return p;
            }

            LIB_UTILITIES_EXPORT SessionReader(
                int                             argc,
                char                           *argv[],
                const std::vector<std::string> &pFilenames,
                const CommSharedPtr            &pComm);

            /// Destructor
            LIB_UTILITIES_EXPORT ~SessionReader();

            LIB_UTILITIES_EXPORT void InitSession(
                const std::vector<std::string> &filenames =
                    std::vector<std::string>());

            /// Provides direct access to the TiXmlDocument object.
            LIB_UTILITIES_EXPORT TiXmlDocument &GetDocument();
            /// Provides direct access to the TiXmlElement specified.
            LIB_UTILITIES_EXPORT TiXmlElement *GetElement(
                const std::string& pPath);
            /// Tests if a specified element is defined in the XML document.
            LIB_UTILITIES_EXPORT bool DefinesElement(
                const std::string& pPath) const;
            /// Returns the filename of the loaded XML document.
            LIB_UTILITIES_EXPORT const std::vector<std::string>
                                                    &GetFilenames() const;
            /// Returns the session name of the loaded XML document.
            LIB_UTILITIES_EXPORT const std::string &GetSessionName() const;
            /// Returns the session name with process rank
            LIB_UTILITIES_EXPORT const std::string  GetSessionNameRank() const;
            /// Returns the communication object.
            LIB_UTILITIES_EXPORT CommSharedPtr GetComm();
            /// Returns if file system shared
            LIB_UTILITIES_EXPORT bool GetSharedFilesystem();
            /// Finalises the session.
            LIB_UTILITIES_EXPORT void Finalise();

            /* ------ PARAMETERS --------*/
            /// Checks if a parameter is specified in the XML document.
            LIB_UTILITIES_EXPORT bool DefinesParameter(
                const std::string &name) const;
            /// Returns the value of the specified parameter.
            LIB_UTILITIES_EXPORT const NekDouble &GetParameter(
                const std::string &pName) const;
            /// Load an integer parameter
            LIB_UTILITIES_EXPORT void LoadParameter(
                const std::string &name,
                int               &var) const;
            /// Check for and load an integer parameter.
            LIB_UTILITIES_EXPORT void LoadParameter(
                const std::string &name,
                int               &var,
                const int         &def) const;
            /// Load a double precision parameter
            LIB_UTILITIES_EXPORT void LoadParameter(
                const std::string &name,
                NekDouble         &var) const;
            /// Check for and load a double-precision parameter.
            LIB_UTILITIES_EXPORT void LoadParameter(
                const std::string &name,
                NekDouble         &var,
                const NekDouble   &def) const;
            /// Set an integer parameter
            LIB_UTILITIES_EXPORT void SetParameter(
                const std::string &name,
                int               &var);
            /// Set a double precision parameter
            LIB_UTILITIES_EXPORT void SetParameter(
                const std::string &name,
                NekDouble         &var);


            /* ------ SOLVER INFO ------ */
            /// Checks if a solver info property is specified.
            LIB_UTILITIES_EXPORT bool DefinesSolverInfo(
                const std::string &name) const;
            /// Returns the value of the specified solver info property.
            LIB_UTILITIES_EXPORT const std::string& GetSolverInfo(
                const std::string &pProperty) const;
            /// Sets the value of the specified solver info property.
            LIB_UTILITIES_EXPORT void SetSolverInfo(
                const std::string &pProperty, const std::string &pValue);
            /// Returns the value of the specified solver info property as enum
            template<typename T>
            inline const T GetSolverInfoAsEnum(const std::string &pName) const;
            /// Returns the value of the specified property and value as enum
            template<typename T>
            inline const T GetValueAsEnum(const std::string &pName,
                                          const std::string &vValue) const;
            /// Check for and load a solver info property.
            LIB_UTILITIES_EXPORT void LoadSolverInfo(
                const std::string &name,
                      std::string &var,
                const std::string &def = "") const;
            /// Check if the value of a solver info property matches.
            LIB_UTILITIES_EXPORT void MatchSolverInfo(
                const std::string &name,
                const std::string &trueval,
                      bool        &var,
                const bool        &def = false) const;
            /// Check if the value of a solver info property matches.
            LIB_UTILITIES_EXPORT bool MatchSolverInfo(
                const std::string &name,
                const std::string &trueval) const;
            /// Check if the value of a solver info property matches.
            template<typename T>
            inline bool MatchSolverInfoAsEnum(
                const std::string &name,
                const T           &trueval) const;
            /// Registers an enumeration value.
            LIB_UTILITIES_EXPORT inline static std::string RegisterEnumValue(
                std::string pEnum,
                std::string pString,
                int         pEnumValue);
            /// Registers the default string value of a solver info property.
            LIB_UTILITIES_EXPORT inline static std::string
                RegisterDefaultSolverInfo(
                const std::string &pName,
                const std::string &pValue);

            /* ----GlobalSysSolnInfo ----- */

            LIB_UTILITIES_EXPORT bool DefinesGlobalSysSolnInfo(
                const std::string &variable,
                const std::string &property) const;

            LIB_UTILITIES_EXPORT const std::string& GetGlobalSysSolnInfo(
                const std::string &variable,
                const std::string &property) const;


            /* ------ GEOMETRIC INFO ------ */
            LIB_UTILITIES_EXPORT std::string GetGeometryType() const;
            /// Checks if a geometric info property is defined.
            LIB_UTILITIES_EXPORT bool DefinesGeometricInfo(
                const std::string &name) const;
            /// Checks for and load a geometric info string property.
            LIB_UTILITIES_EXPORT void LoadGeometricInfo(
                const std::string &name,
                      std::string &var,
                const std::string &def = "") const;
            /// Checks for and loads a geometric info boolean property.
            LIB_UTILITIES_EXPORT void LoadGeometricInfo(
                const std::string &name,
                      bool        &var,
                const bool        &def = false) const;
            /// Checks for and loads a geometric info double-precision property.
            LIB_UTILITIES_EXPORT void LoadGeometricInfo(
                const std::string &name,
                      NekDouble   &var,
                const NekDouble   &def = 0.0) const;
            /// Check if the value of a geometric info string property matches.
            LIB_UTILITIES_EXPORT void MatchGeometricInfo(
                const std::string &name,
                const std::string &trueval,
                      bool        &var,
                const bool        &def = false) const;

            /* ------ VARIABLES ------ */
            /// Returns the name of the variable specified by the given index.
            LIB_UTILITIES_EXPORT const std::string& GetVariable(
                const unsigned int &idx) const;
            LIB_UTILITIES_EXPORT void SetVariable(
                const unsigned int &idx,
                      std::string newname);

            /// Returns the names of all variables.
            LIB_UTILITIES_EXPORT std::vector<std::string> GetVariables() const;

            /* ------ FUNCTIONS ------*/
            /// Checks if a specified function is defined in the XML document.
            LIB_UTILITIES_EXPORT bool DefinesFunction(
                const std::string &name) const;
            /// Checks if a specified function has a given variable defined.
            LIB_UTILITIES_EXPORT bool DefinesFunction(
                const std::string &name,
                const std::string &variable,
                const int pDomain = 0) const;
            /// Returns an EquationSharedPtr to a given function variable.
            LIB_UTILITIES_EXPORT EquationSharedPtr GetFunction(
                const std::string &name,
                const std::string &variable,
                const int pDomain = 0) const;
            /// Returns an EquationSharedPtr to a given function variable index.
            LIB_UTILITIES_EXPORT EquationSharedPtr GetFunction(
                const std::string  &name,
                const unsigned int &var,
                const int pDomain = 0) const;
            /// Returns the type of a given function variable.
            LIB_UTILITIES_EXPORT enum FunctionType GetFunctionType(
                const std::string &name,
                const std::string &variable,
                const int pDomain = 0) const;
            /// Returns the type of a given function variable index.
            LIB_UTILITIES_EXPORT enum FunctionType GetFunctionType(
                const std::string  &pName,
                const unsigned int &pVar,
                const int pDomain = 0) const;
            /// Returns the filename to be loaded for a given variable.
            LIB_UTILITIES_EXPORT std::string GetFunctionFilename(
                const std::string &name,
                const std::string &variable,
                const int pDomain = 0) const;
            /// Returns the filename to be loaded for a given variable index.
            LIB_UTILITIES_EXPORT std::string GetFunctionFilename(
                const std::string  &name,
                const unsigned int &var,
                const int pDomain = 0) const;
            /// Returns the filename variable to be loaded for a given variable
            /// index.
            LIB_UTILITIES_EXPORT std::string GetFunctionFilenameVariable(
                const std::string  &name,
                const std::string &variable,
                const int pDomain = 0) const;

            /// Returns the instance of the Interpreter specific to this
            /// session.
            LIB_UTILITIES_EXPORT InterpreterSharedPtr GetInterpreter();

            /* ------ TAGS ------ */
            /// Checks if a specified tag is defined.
            LIB_UTILITIES_EXPORT bool DefinesTag(
                const std::string& pName) const;
            /// Sets a specified tag.
            LIB_UTILITIES_EXPORT void SetTag(
                const std::string& pName,
                const std::string& pValue);
            /// Returns the value of a specified tag.
            LIB_UTILITIES_EXPORT const std::string &GetTag(
                const std::string& pName) const;

            /* ------ FILTERS ------ */
            LIB_UTILITIES_EXPORT const FilterMap& GetFilters() const;

            /* ------ CMDLINE ARGUMENTS ------- */
            /// Checks if a specified cmdline argument has been given.
            LIB_UTILITIES_EXPORT bool DefinesCmdLineArgument(
                const std::string& pName) const;
            /// Retrieves a command-line argument value.
            template <typename T>
            T GetCmdLineArgument(
                const std::string& pName) const
            {
                return m_cmdLineOptions.find(pName)->second.as<T>();
            }
            /// Registers a command-line argument with the session reader.
            LIB_UTILITIES_EXPORT inline static std::string
              RegisterCmdLineArgument(
                const std::string &pName,
                const std::string &pShortName,
                const std::string &pDescription);
            /// Registers a command-line flag with the session reader.
            LIB_UTILITIES_EXPORT inline static std::string
              RegisterCmdLineFlag(
                const std::string &pName,
                const std::string &pShortName,
                const std::string &pDescription);

            /// Substitutes expressions defined in the XML document.
            LIB_UTILITIES_EXPORT void SubstituteExpressions(std::string &expr);

            LIB_UTILITIES_EXPORT void SetUpXmlDoc();

        private:
            boost::program_options::variables_map m_cmdLineOptions;

            /// Communication object.
            CommSharedPtr                             m_comm;
            /// Filenames
            std::vector<std::string>                  m_filenames;
            /// Session name of the loaded XML document (filename minus ext).
            std::string                               m_sessionName;
            /// Pointer to the loaded XML document.
            TiXmlDocument*                            m_xmlDoc;
            /// Parameters.
            ParameterMap                              m_parameters;
            /// Solver information properties.
            SolverInfoMap                             m_solverInfo;
            /// Geometric information properties.
            GeometricInfoMap                          m_geometricInfo;
            /// Expressions.
            ExpressionMap                             m_expressions;
            /// Interpreter instance.
            InterpreterSharedPtr                      m_interpreter;
            /// Functions.
            FunctionMap                               m_functions;
            /// Variables.
            VariableList                              m_variables;
            /// Custom tags.
            TagMap                                    m_tags;
            /// Filters map.
            FilterMap                                 m_filters;
            /// Be verbose
            bool                                      m_verbose;
            /// Running on a shared filesystem
            bool                                      m_sharedFilesystem;
            /// String to enumeration map for Solver Info parameters.
            LIB_UTILITIES_EXPORT static EnumMapList&  GetSolverInfoEnums();
            /// Default solver info options.
            LIB_UTILITIES_EXPORT static SolverInfoMap& GetSolverInfoDefaults();
            /// GlobalSysSoln Info map.
            LIB_UTILITIES_EXPORT static GloSysSolnInfoList& GetGloSysSolnList();
            /// CmdLine argument map.
            LIB_UTILITIES_EXPORT static CmdLineArgMap& GetCmdLineArgMap();

            /// Main constructor
            LIB_UTILITIES_EXPORT SessionReader(
                int                             argc,
                char                           *argv[]);

            LIB_UTILITIES_EXPORT void TestSharedFilesystem();

            /// Parse the program arguments and fill #m_cmdLineOptions
            std::vector<std::string> ParseCommandLineArguments(
                int argc, char *argv[]);

            /// Parse the session name.
            std::string ParseSessionName(std::vector<std::string> &filenames);

            /// Loads an xml file into a tinyxml doc and decompresses if needed
            LIB_UTILITIES_EXPORT void LoadDoc(
                const std::string &pFilename,
                TiXmlDocument* pDoc) const;
            /// Creates an XML document from a list of input files.
            LIB_UTILITIES_EXPORT TiXmlDocument *MergeDoc(
                const std::vector<std::string> &pFilenames) const;
            /// Loads and parses the specified file.
            LIB_UTILITIES_EXPORT void ParseDocument();
            /// Loads the given XML document and instantiates an appropriate
            /// communication object.
            LIB_UTILITIES_EXPORT void CreateComm(
                int               &argc,
                char*              argv[]);

            /// Partitions the comm object based on session parameters.
            LIB_UTILITIES_EXPORT void PartitionComm();

            /// Reads the PARAMETERS section of the XML document.
            LIB_UTILITIES_EXPORT void ReadParameters(TiXmlElement *conditions);
            /// Reads the SOLVERINFO section of the XML document.
            LIB_UTILITIES_EXPORT void ReadSolverInfo(TiXmlElement *conditions);
            /// Reads the GLOBALSYSSOLNINFO section of the XML document.
            LIB_UTILITIES_EXPORT void ReadGlobalSysSolnInfo(
                    TiXmlElement *conditions);
            /// Reads the EXPRESSIONS section of the XML document.
            LIB_UTILITIES_EXPORT void ReadExpressions(TiXmlElement *conditions);
            /// Reads the VARIABLES section of the XML document.
            LIB_UTILITIES_EXPORT void ReadVariables(TiXmlElement *conditions);
            /// Reads the FUNCTIONS section of the XML document.
            LIB_UTILITIES_EXPORT void ReadFunctions(TiXmlElement *conditions);
            /// Reads the FILTERS section of the XML document.
            LIB_UTILITIES_EXPORT void ReadFilters(TiXmlElement *filters);
            /// Enforce parameters from command line arguments.
            LIB_UTILITIES_EXPORT void CmdLineOverride();
            /// Check values of solver info options are valid.
            LIB_UTILITIES_EXPORT void VerifySolverInfo();

            /// Parse a string in the form lhs = rhs.
            LIB_UTILITIES_EXPORT void ParseEquals(
                const std::string &line,
                      std::string &lhs,
                      std::string &rhs);
        };


        /**
         *
         */
        template<typename T>
        inline bool SessionReader::MatchSolverInfoAsEnum(
            const std::string &name, const T &trueval) const
        {
            return (GetSolverInfoAsEnum<T>(name) == trueval);
        }


        /**
         *
         */
        template<typename T>
        inline const T SessionReader::GetSolverInfoAsEnum(
            const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ASSERTL0(DefinesSolverInfo(vName),
                     "Solver info '" + pName + "' not defined.");

            std::string vValue = GetSolverInfo(vName);
            auto x = GetSolverInfoEnums().find(vName);
            ASSERTL0(x != GetSolverInfoEnums().end(),
                     "Enum for SolverInfo property '" + pName + "' not found.");

            auto y = x->second.find(vValue);
            ASSERTL0(y != x->second.end(),
                     "Value of SolverInfo property '" + pName +
                     "' is invalid.");

            return T(y->second);
        }



        /**
         *
         */
        template<typename T>
        inline const T SessionReader::GetValueAsEnum(
            const std::string &pName,
            const std::string &pValue) const
        {
            std::string vName  = boost::to_upper_copy(pName);

            auto x = GetSolverInfoEnums().find(vName);
            ASSERTL0(x != GetSolverInfoEnums().end(),
                     "Enum for property '" + pName + "' not found.");

            auto y = x->second.find(pValue);
            ASSERTL0(y != x->second.end(),
                     "Value of property '" + pValue + "' is invalid.");
            return T(y->second);
        }



        /**
         * A set of valid values for a given solver info property may be
         * registered using this function. It must be called statically during
         * the initialisation of a static variable. For example:
         *
         * @code
         * std::string GlobalLinSys::lookupIds[2] = {
         *     LibUtilities::SessionReader::RegisterEnumValue(
         *                                  "GlobalSysSoln",
         *                                  "DirectFull",
         *                                  MultiRegions::eDirectFullMatrix),
         *     LibUtilities::SessionReader::RegisterEnumValue(
         *                                  "GlobalSysSoln",
         *                                  "DirectStaticCond",
         *                                  MultiRegions::eDirectStaticCond)
         * }
         * @endcode
         *
         * @param   pEnum       The name of the property.
         * @param   pString     A valid value for the property.
         * @param   pEnumValue  An enumeration value corresponding to this
         *                      value.
         *
         * @return The value for the property provided by #pString.
         */
        inline std::string SessionReader::RegisterEnumValue(
            std::string pEnum, std::string pString, int pEnumValue)
        {
            std::string vEnum = boost::to_upper_copy(pEnum);
            auto x = GetSolverInfoEnums().find(vEnum);

            if (x == GetSolverInfoEnums().end())
            {
                GetSolverInfoEnums()[vEnum] = EnumMap();
                x = GetSolverInfoEnums().find(vEnum);
            }

            x->second[pString] = pEnumValue;
            return pString;
        }


        /**
         * A default value for a given solver info property may be registered
         * using this function. The property will take this value until it is
         * overwritten by a value specified in the XML document, or specified
         * as a command-line argument. Usage has the form:
         *
         * @code
         * std::string GlobalLinSys::def
         *     = LibUtilities::SessionReader::RegisterDefaultSolverInfo(
         *         "GlobalSysSoln","DirectMultiLevelStaticCond");
         * @endcode
         *
         * @param   pName       The name of the property.
         * @param   pValue      The default value of the property.
         *
         * @return The default value of the property provided by #pValue.
         */
        inline std::string SessionReader::RegisterDefaultSolverInfo(
            const std::string &pName,
            const std::string &pValue)
        {
            std::string vName = boost::to_upper_copy(pName);
            GetSolverInfoDefaults()[vName] = pValue;
            return pValue;
        }


        /**
         *
         */
        inline std::string SessionReader::RegisterCmdLineArgument(
            const std::string &pName,
            const std::string &pShortName,
            const std::string &pDescription)
        {
            ASSERTL0(!pName.empty(), "Empty name for cmdline argument.");
            CmdLineArg x;
            x.shortName = pShortName;
            x.description = pDescription;
            x.isFlag = false;
            GetCmdLineArgMap()[pName] = x;
            return pName;
        }


        /**
         *
         */
        inline std::string SessionReader::RegisterCmdLineFlag(
            const std::string &pName,
            const std::string &pShortName,
            const std::string &pDescription)
        {
            ASSERTL0(!pName.empty(), "Empty name for cmdline argument.");
            CmdLineArg x;
            x.shortName = pShortName;
            x.description = pDescription;
            x.isFlag = true;
            GetCmdLineArgMap()[pName] = x;
            return pName;
        }
    }
}

#endif

