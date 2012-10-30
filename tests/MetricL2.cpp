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
    
    MetricL2::MetricL2(TiXmlElement *metric) : MetricRegex(metric)
    {
        // Set up the regular expression. This (optionally) matches a variable
        // name if it exists: first field is variable name, second field is L2
        // error.
        m_regex = "^L 2 error\\s*(?:\\(variable "
                  "(\\w+)\\))?\\s*:\\s*([+-]?\\d.+\\d|0).*";
        
        // Find the L2 error to match against.
        TiXmlElement *value = metric->FirstChildElement("value");
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

            std::vector<MetricRegexFieldValue> tmp(2);
            tmp[0] = var;
            tmp[1] = val;
            m_matches.push_back(tmp);
            
            value = value->NextSiblingElement("value");
        }
    }

}
