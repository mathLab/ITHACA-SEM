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

        const std::string& getFilename();
        const std::string& getSolverInfo(const std::string &pProperty);

        NekDouble getParameter(std::string pName);

        /// Check for and load an integer parameter
        /// Check for and load a double precision parameter
        void loadParameter(const std::string name, int &var, int def = 0);
        void loadParameter(const std::string name, NekDouble& var, const NekDouble def= 0.0);
        bool definesParameter(const std::string name);

        void loadSolverInfo(const std::string name, std::string& var, const std::string def = "");
        void matchSolverInfo(const std::string name, const std::string trueval, bool& var, const bool def = false);
        bool definesSolverInfo(const std::string name);

    private:
        std::string                 mFilename;
        TiXmlDocument*              mXmlDoc;

        SolverInfoMap               mSolverInfo;
        ParameterMap                mParameters;

        void ReadParameters(TiXmlElement *conditions);
        void ReadSolverInfo(TiXmlElement *conditions);

        /// Perform a case-insensitive string comparison.
        int NoCaseStringCompare(const std::string & s1, const std::string& s2);

    };

    typedef boost::shared_ptr<SessionReader> SessionReaderSharedPtr;
}

#endif
