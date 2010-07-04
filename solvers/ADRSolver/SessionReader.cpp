// Use the stl version, primarily for string.
#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <iostream>

#include <tinyxml/tinyxml.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <SpatialDomains/Equation.h>

#include <ADRSolver/SessionReader.h>

namespace Nektar
{

    SessionReader::SessionReader(std::string& pFilename)
        : mFilename(pFilename)
    {
        mXmlDoc = new TiXmlDocument(pFilename);
        ASSERTL0(mXmlDoc, "Failed to create XML document object.");

        bool loadOkay = mXmlDoc->LoadFile();
        ASSERTL0(loadOkay, (std::string("Unable to load file: ") +
            pFilename).c_str());

        TiXmlHandle docHandle(mXmlDoc);

        TiXmlNode* n = NULL;
        TiXmlElement* e = NULL;

        /// Look for all data in CONDITIONS block.
        e = docHandle.FirstChildElement("NEKTAR").FirstChildElement("CONDITIONS").Element();
        ASSERTL0(e, "Unable to find CONDITIONS tag in file.");

        ReadParameters(e);
        ReadSolverInfo(e);
    }

    SessionReader::SessionReader(const SessionReader& pSrc)
    {
        mFilename = pSrc.mFilename;
        mXmlDoc   = pSrc.mXmlDoc;
        mSolverInfo = pSrc.mSolverInfo;
        mParameters = pSrc.mParameters;
    }

    SessionReader::~SessionReader()
    {

    }

    const std::string& SessionReader::getFilename()
    {
        return mFilename;
    }

    const std::string& SessionReader::getSolverInfo(const std::string &pProperty)
    {
        SolverInfoMap::iterator slvIter = mSolverInfo.find(pProperty);

        ASSERTL1(slvIter != mSolverInfo.end(),
            (std::string("Unable to find requested property: ") + pProperty).c_str());

        return slvIter->second;
    }

    NekDouble SessionReader::getParameter(std::string pName)
    {
        ParameterMap::iterator paramMapIter = mParameters.find(pName);

        ASSERTL0(paramMapIter != mParameters.end(),
            (std::string("Unable to find requested parameter: ") + pName).c_str());

        return paramMapIter->second;
    }

    void SessionReader::loadParameter(const std::string pName, int &pVar, int pDefault)
    {
        ParameterMap::iterator paramMapIter = mParameters.find(pName);
        if(paramMapIter != mParameters.end())
        {
            pVar = paramMapIter->second;
        }
        else
        {
            pVar  = pDefault;
        }
    }

    void SessionReader::loadParameter(const std::string pName, NekDouble& pVar, const NekDouble pDefault)
    {
        ParameterMap::iterator paramMapIter = mParameters.find(pName);
        if(paramMapIter != mParameters.end())
        {
            pVar = paramMapIter->second;
        }
        else
        {
            pVar  = pDefault;
        }
    }

    bool SessionReader::definesParameter(const std::string pName)
    {
        ParameterMap::iterator paramMapIter = mParameters.find(pName);
        return (paramMapIter != mParameters.end());
    }

    void SessionReader::loadSolverInfo(const std::string pName, std::string& pVar, const std::string pDefault)
    {
        SolverInfoMap::iterator solverInfoMapIter = mSolverInfo.find(pName);
        if(solverInfoMapIter != mSolverInfo.end())
        {
            pVar = solverInfoMapIter->second;
        }
        else
        {
            pVar  = pDefault;
        }
    }

    void SessionReader::matchSolverInfo(const std::string pName, const std::string pTrueVal, bool& pVar, const bool pDefault)
    {
        SolverInfoMap::iterator solverInfoMapIter = mSolverInfo.find(pName);
        if(solverInfoMapIter != mSolverInfo.end())
        {
            pVar = (NoCaseStringCompare(solverInfoMapIter->second, pTrueVal) == 0);
        }
        else
        {
            pVar  = pDefault;
        }
    }

    bool SessionReader::definesSolverInfo(const std::string pName)
    {
        SolverInfoMap::iterator solverInfoMapIter = mSolverInfo.find(pName);
        return (solverInfoMapIter != mSolverInfo.end());
    }

    void SessionReader::ReadSolverInfo(TiXmlElement *conditions)
    {
        TiXmlElement *solverInfoElement = conditions->FirstChildElement("SOLVERINFO");

        if (solverInfoElement)
        {
            TiXmlElement *solverInfo = solverInfoElement->FirstChildElement("I");

            while (solverInfo)
            {
                std::string solverProperty = solverInfo->Attribute("PROPERTY");
                // make sure that solver property is capitalised
                transform(solverProperty.begin(), solverProperty.end(), solverProperty.begin(), (int(*)(int))std::toupper);
                ASSERTL0(!solverProperty.empty(), "Unable to find PROPERTY value.");

                std::string solverValue    = solverInfo->Attribute("VALUE");
                ASSERTL0(!solverValue.empty(),"Unable to find VALUE string");

                SolverInfoMap::iterator solverInfoIter = mSolverInfo.find(solverProperty);

                ASSERTL0(solverInfoIter == mSolverInfo.end(),
                         (std::string("SolverInfo value: ") + solverProperty
                          + std::string(" already specified.")).c_str());

                // Set Variable
                mSolverInfo[solverProperty] = solverValue;
                solverInfo = solverInfo->NextSiblingElement("I");
            }
        }
    }


    void SessionReader::ReadParameters(TiXmlElement *conditions)
    {
        TiXmlElement *parametersElement = conditions->FirstChildElement("PARAMETERS");

        // See if we have parameters defined.  They are optional so we go on if not.
        if (parametersElement)
        {
            TiXmlElement *parameter = parametersElement->FirstChildElement("P");
            LibUtilities::ExpressionEvaluator expEvaluator;

            // Multiple nodes will only occur if there is a comment in between
            // definitions.
            while (parameter)
            {
                TiXmlNode *node = parameter->FirstChild();

                while (node && node->Type() != TiXmlNode::TEXT)
                {
                    node = node->NextSibling();
                }

                if (node)
                {
                    // Format is "paramName = value"
                    std::string line = node->ToText()->Value();

                    /// Pull out lhs and rhs and eliminate any spaces.
                    int beg=line.find_first_not_of(" ");
                    int end=line.find_first_of("=");
                    std::string lhs = line.substr(line.find_first_not_of(" "), end-beg-1);
                    lhs = lhs.substr(0, lhs.find_last_not_of(" ")+1);

                    std::string rhs = line.substr(line.find_last_of("=")+1);
                    rhs = rhs.substr(rhs.find_first_not_of(" "));
                    rhs = rhs.substr(0, rhs.find_last_not_of(" ")+1);

                    /// We want the list of parameters to have their RHS evaluated,
                    /// so we use the expression evaluator to do the dirty work.
                    if (!lhs.empty() && !rhs.empty())
                    {
                        NekDouble value=0.0;
                        expEvaluator.DefineFunction("", rhs);
                        value =  expEvaluator.Evaluate();
                        mParameters[lhs] = value;
                        expEvaluator.SetParameter(lhs, value);
                    }
                }

                parameter = parameter->NextSiblingElement();
            }

            // Set ourselves up for evaluation later.
            //SpatialDomains::Equation::SetConstParameters(mParameters);
        }
    }

    int SessionReader::NoCaseStringCompare(const std::string & s1, const std::string& s2)
    {
        //if (s1.size() < s2.size()) return -1;
        //if (s1.size() > s2.size()) return 1;

        std::string::const_iterator it1=s1.begin();
        std::string::const_iterator it2=s2.begin();

        //stop when either string's end has been reached
        while ( (it1!=s1.end()) && (it2!=s2.end()) )
        {
            if(::toupper(*it1) != ::toupper(*it2)) //letters differ?
            {
                // return -1 to indicate smaller than, 1 otherwise
                return (::toupper(*it1)  < ::toupper(*it2)) ? -1 : 1;
            }

            //proceed to the next character in each string
            ++it1;
            ++it2;
        }

        size_t size1=s1.size();
        size_t size2=s2.size();// cache lengths

        //return -1,0 or 1 according to strings' lengths
        if (size1==size2)
        {
            return 0;
        }

        return (size1 < size2) ? -1 : 1;
    }
}
