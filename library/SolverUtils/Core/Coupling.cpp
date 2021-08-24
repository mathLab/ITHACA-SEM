////////////////////////////////////////////////////////////////////////////////
//
// File: Coupling.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Coupling
//
////////////////////////////////////////////////////////////////////////////////

#include "Coupling.h"

#include <LibUtilities/BasicUtils/ParseUtils.h>

namespace Nektar
{
namespace SolverUtils
{

using namespace std;

CouplingFactory &GetCouplingFactory()
{
    static CouplingFactory instance;
    return instance;
}

Coupling::Coupling(MultiRegions::ExpListSharedPtr field)
    : m_couplingName(""), m_evalField(field), m_nSendVars(0), m_sendSteps(0),
      m_nRecvVars(0), m_recvSteps(0)
{
    m_config["RECEIVESTEPS"]     = "0";
    m_config["RECEIVEVARIABLES"] = "";

    m_config["SENDSTEPS"]     = "0";
    m_config["SENDVARIABLES"] = "";
}

void Coupling::v_Init()
{
    LibUtilities::SessionReaderSharedPtr session = m_evalField->GetSession();

    TiXmlElement *vCoupling = session->GetElement("Nektar/Coupling");
    ASSERTL0(vCoupling, "Invalid Coupling config");

    vCoupling->QueryStringAttribute("NAME", &m_couplingName);
    ASSERTL0(m_couplingName.size(), "No Coupling NAME attribute set");

    TiXmlElement *element = vCoupling->FirstChildElement("I");
    while (element)
    {
        std::stringstream tagcontent;
        tagcontent << *element;
        // read the property name
        ASSERTL0(element->Attribute("PROPERTY"),
                 "Missing PROPERTY attribute in Coupling section "
                 "XML element: \n\t'" +
                     tagcontent.str() + "'");
        std::string property = element->Attribute("PROPERTY");
        ASSERTL0(!property.empty(),
                 "PROPERTY attribute must be non-empty in XML "
                 "element: \n\t'" +
                     tagcontent.str() + "'");

        // make sure that solver property is capitalised
        std::string propertyUpper = boost::to_upper_copy(property);

        CouplingConfigMap::const_iterator x = m_config.find(propertyUpper);
        ASSERTL0(x != m_config.end(),
                 "Invalid PROPERTY attribute in Coupling section "
                 "XML element: \n\t'" +
                     tagcontent.str() + "'");

        // read the value
        ASSERTL0(element->Attribute("VALUE"),
                 "Missing VALUE attribute in Coupling section "
                 "XML element: \n\t'" +
                     tagcontent.str() + "'");
        std::string value = element->Attribute("VALUE");
        ASSERTL0(!value.empty(),
                 "VALUE attribute must be non-empty in XML "
                 "element: \n\t'" +
                     tagcontent.str() + "'");

        // Set Variable
        m_config[propertyUpper] = value;

        element = element->NextSiblingElement("I");
    }

    // mangle config into variables. This is ugly
    ParseUtils::GenerateVector(m_config["RECEIVEVARIABLES"], m_recvFieldNames);
    m_nRecvVars = m_recvFieldNames.size();

    ParseUtils::GenerateVector(m_config["SENDVARIABLES"], m_sendFieldNames);
    m_nSendVars = m_sendFieldNames.size();

    m_recvSteps = boost::lexical_cast<int>(m_config["RECEIVESTEPS"]);
    m_sendSteps = boost::lexical_cast<int>(m_config["SENDSTEPS"]);

    if (session->GetComm()->GetRank() == 0 &&
        session->DefinesCmdLineArgument("verbose") && m_config.size() > 0)
    {
        cout << "Coupling Config:" << endl;
        CouplingConfigMap::iterator x;
        for (x = m_config.begin(); x != m_config.end(); ++x)
        {
            cout << "\t" << x->first << " = '" << x->second << "'" << endl;
        }
    }
}

vector<int> Coupling::GenerateVariableMapping(vector<string> &vars,
                                              vector<string> &transVars)
{
    vector<int> transToVars;
    Array<OneD, Array<OneD, NekDouble> > sendField(transVars.size());
    for (int i = 0; i < transVars.size(); ++i)
    {
        auto it2 = find(vars.begin(), vars.end(), transVars[i]);
        ASSERTL0(it2 != vars.end(),
                 "send variable " + transVars[i] + " not found");
        int id = distance(vars.begin(), it2);

        transToVars.push_back(id);
    }

    return transToVars;
}
}
}
