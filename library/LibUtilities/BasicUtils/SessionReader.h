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
// License for the specific language governing rights and limitations under
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

#include <iostream>
#include <map>
#include <string>

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

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
        typedef std::map<std::string, EquationSharedPtr>  EquationMap;
        typedef std::map<std::string, std::string>  TagMap;

        enum FunctionType
        {
            eFunctionTypeNone,
            eFunctionTypeExpression,
            eFunctionTypeFile,
            eSIZE_FunctionType
        };
        const char* const FunctionTypeMap[] =
        {
            "No Function type",
            "Expression",
            "File"
        };
        struct FunctionDefinition
        {
            enum FunctionType m_type;
            std::string       m_filename;
            EquationMap       m_expressions;
        };
        typedef std::map<std::string, FunctionDefinition > FunctionMap;


        class SessionReader
        {
        public:
            LIB_UTILITIES_EXPORT SessionReader(std::string& pFilename);
            LIB_UTILITIES_EXPORT SessionReader(const SessionReader& pSrc);
            LIB_UTILITIES_EXPORT ~SessionReader();

            LIB_UTILITIES_EXPORT TiXmlDocument& GetDocument();
            LIB_UTILITIES_EXPORT TiXmlElement* GetElement(const std::string& pPath);
            LIB_UTILITIES_EXPORT bool DefinesElement(const std::string& pPath);
            LIB_UTILITIES_EXPORT const std::string& GetFilename();
            LIB_UTILITIES_EXPORT const std::string& GetSolverInfo(const std::string &pProperty);

            LIB_UTILITIES_EXPORT NekDouble GetParameter(std::string pName);

            /// Check for and load an integer parameter
            /// Check for and load a double precision parameter
            LIB_UTILITIES_EXPORT void LoadParameter(const std::string name, int &var, int def = 0);
            LIB_UTILITIES_EXPORT void LoadParameter(const std::string name, NekDouble& var, const NekDouble def= 0.0);
            LIB_UTILITIES_EXPORT bool DefinesParameter(const std::string name);

            LIB_UTILITIES_EXPORT void LoadSolverInfo(const std::string name, std::string& var, const std::string def = "");
            LIB_UTILITIES_EXPORT void MatchSolverInfo(const std::string name, const std::string trueval, bool& var, const bool def = false);
            LIB_UTILITIES_EXPORT bool DefinesSolverInfo(const std::string name);

            LIB_UTILITIES_EXPORT void LoadGeometricInfo(const std::string name, std::string& var, const std::string def = "");
            LIB_UTILITIES_EXPORT void LoadGeometricInfo(const std::string name, bool& var, const bool def = false);
            LIB_UTILITIES_EXPORT void MatchGeometricInfo(const std::string name, const std::string trueval, bool& var, const bool def = false);
            LIB_UTILITIES_EXPORT bool DefinesGeometricInfo(const std::string name);

            LIB_UTILITIES_EXPORT std::string GetVariable(const unsigned int idx) const;

            LIB_UTILITIES_EXPORT EquationSharedPtr GetFunction(const std::string& name, const std::string& variable) const;
            LIB_UTILITIES_EXPORT EquationSharedPtr GetFunction(const std::string& name, unsigned int var) const;
            LIB_UTILITIES_EXPORT enum FunctionType GetFunctionType(const std::string& name) const;
            LIB_UTILITIES_EXPORT std::string GetFunctionFilename(const std::string& name) const;
            LIB_UTILITIES_EXPORT bool DefinesFunction(const std::string& name) const;
            LIB_UTILITIES_EXPORT bool DefinesFunction(const std::string& name, const std::string& variable) const;

            LIB_UTILITIES_EXPORT bool DefinesTag(const std::string& pName);
            LIB_UTILITIES_EXPORT void SetTag(const std::string& pName, const std::string& pValue);
            LIB_UTILITIES_EXPORT const std::string GetTag(const std::string& pName);

        private:
            std::string                 m_filename;
            TiXmlDocument*              m_xmlDoc;

            SolverInfoMap               m_solverInfo;
            ParameterMap                m_parameters;
            GeometricInfoMap            m_geometricInfo;
            ExpressionMap               m_expressions;
            FunctionMap                 m_functions;
            VariableList                m_variables;
            TagMap                      m_tags;

            void ReadParameters(TiXmlElement *conditions);
            void ReadSolverInfo(TiXmlElement *conditions);
            void ReadGeometricInfo(TiXmlElement *geometry);
            void ReadExpressions(TiXmlElement *conditions);
            void ReadVariables(TiXmlElement *conditions);
            void ReadFunctions(TiXmlElement *conditions);

            /// Perform a case-insensitive string comparison.
            int NoCaseStringCompare(const std::string & s1, const std::string& s2);
            void SubstituteExpressions(std::string &expr);
        };

        typedef boost::shared_ptr<SessionReader> SessionReaderSharedPtr;
    }
}

#endif

