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
#ifndef NEKTAR_SOLVERS_ADRSOLVER_SESSIONREADER_H
#define NEKTAR_SOLVERS_ADRSOLVER_SESSIONREADER_H

#include <iostream>
#include <map>
#include <string>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

class TiXmlElement;
class TiXmlDocument;

namespace Nektar
{
    typedef std::map<std::string, std::string>  SolverInfoMap;
    typedef std::map<std::string, NekDouble>    ParameterMap;

    class SessionReader
    {
    public:
        SessionReader(std::string& pFilename);
        SessionReader(const SessionReader& pSrc);
        ~SessionReader();

        const std::string& GetFilename();
        const std::string& GetSolverInfo(const std::string &pProperty);

        NekDouble GetParameter(std::string pName);

        /// Check for and load an integer parameter
        /// Check for and load a double precision parameter
        void LoadParameter(const std::string name, int &var, int def = 0);
        void LoadParameter(const std::string name, NekDouble& var, const NekDouble def= 0.0);
        bool DefinesParameter(const std::string name);

        void LoadSolverInfo(const std::string name, std::string& var, const std::string def = "");
        void MatchSolverInfo(const std::string name, const std::string trueval, bool& var, const bool def = false);
        bool DefinesSolverInfo(const std::string name);

    private:
        std::string                 m_filename;
        TiXmlDocument*              m_xmlDoc;

        SolverInfoMap               m_solverInfo;
        ParameterMap                m_parameters;

        void ReadParameters(TiXmlElement *conditions);
        void ReadSolverInfo(TiXmlElement *conditions);

        /// Perform a case-insensitive string comparison.
        int NoCaseStringCompare(const std::string & s1, const std::string& s2);

    };

    typedef boost::shared_ptr<SessionReader> SessionReaderSharedPtr;
}

#endif

