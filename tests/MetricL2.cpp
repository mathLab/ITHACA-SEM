///////////////////////////////////////////////////////////////////////////////
//
// File: MetricL2.cpp
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
// Description: Implementation of the L2 metric.
//
///////////////////////////////////////////////////////////////////////////////

#include <MetricL2.h>

namespace Nektar
{
    std::string MetricL2::type = GetMetricFactory().
        RegisterCreatorFunction("L2", MetricL2::create);
    
    // Guess default tolerance for generation routine.
    std::string MetricL2::defaultTolerance = "1e-12";
    
    MetricL2::MetricL2(TiXmlElement *metric, bool generate) : 
        MetricRegex(metric, generate)
    {
        // Set up the regular expression. This (optionally) matches a variable
        // name if it exists: first field is variable name, second field is L2
        // error.
        m_regex = "^L 2 error\\s*(?:\\(variable "
                  "(\\w+)\\))?\\s*:\\s*([+-]?\\d.+\\d|-?0).*";

        // Find the L2 error to match against.
        TiXmlElement *value = metric->FirstChildElement("value");
        ASSERTL0(value || m_generate, "Missing value tag for L2 metric!");
        while (value)
        {
            // Set up a match with two fields which correspond with the
            // subexpression above. The first is the variable name, second is
            // the L2 error.
            ASSERTL0(value->Attribute("tolerance"),
                     "Missing tolerance in L2 metric");
            ASSERTL0(value->GetText() || value->GetText() == "",
                     "Missing value in L2 metric.");

            MetricRegexFieldValue var;
            if (value->Attribute("variable"))
            {
                var.m_value = value->Attribute("variable");
            }

            MetricRegexFieldValue val;
            val.m_value = value->GetText();
            val.m_useTolerance = true;
            val.m_tolerance = atof(value->Attribute("tolerance"));

            if (!m_generate)
            {
                std::vector<MetricRegexFieldValue> tmp(2);
                tmp[0] = var;
                tmp[1] = val;
                m_matches.push_back(tmp);
            }
            else
            {
                m_varTolerance[var.m_value] = value->Attribute("tolerance");
            }

            value = value->NextSiblingElement("value");
        }
    }
    
    void MetricL2::v_Generate(std::istream& pStdout, std::istream& pStderr)
    {
        // Run MetricRegex to generate matches.
        MetricRegex::v_Generate(pStdout, pStderr);

        // First remove all existing values.
        m_metric->Clear();

        // Now create new values.
        for (int i = 0; i < m_matches.size(); ++i)
        {
            ASSERTL0(m_matches[i].size() == 2,
                     "Wrong number of matches for regular expression.");

            bool          tolSet = false;
            std::string   tol    = MetricL2::defaultTolerance;
            TiXmlElement *value  = new TiXmlElement("value");
            
            // See if there is a tolerance found already for this variable
            // (including empty variables).
            std::map<std::string,std::string>::iterator it = 
                m_varTolerance.find(m_matches[i][0].m_value);

            if (it != m_varTolerance.end())
            {
                tol    = it->second;
                tolSet = true;
            }                
            
            if (m_matches[i][0].m_value.size() > 0)
            {
                value->SetAttribute("variable", m_matches[i][0].m_value);
                
                if (m_matches[i][0].m_value == "p" && !tolSet)
                {
                    // Set lower tolerance for pressure fields automatically if
                    // we haven't already got a tolerance from the existing
                    // file.
                    tol = "1e-8";
                }
            }
            
            value->SetAttribute("tolerance", tol);
            value->LinkEndChild(new TiXmlText(m_matches[i][1].m_value));
            m_metric->LinkEndChild(value);
        }
    }
}
