///////////////////////////////////////////////////////////////////////////////
//
// File: MetricLInf.cpp
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
// Description: Implementation of the LInf metric.
//
///////////////////////////////////////////////////////////////////////////////

#include <MetricLInf.h>

namespace Nektar
{
    std::string MetricLInf::type = GetMetricFactory().
        RegisterCreatorFunction("LINF", MetricLInf::create);

    MetricLInf::MetricLInf(TiXmlElement *metric) : MetricRegex(metric)
    {
        // Set up the regular expression. This (optionally) matches a variable
        // name if it exists: first field is variable name, second field is L2
        // error.
        m_regex =
           "^L inf error\\s*(?:\\(variable (\\w+)\\))?\\s*:\\s*([+-]?\\d.+\\d|0).*";

        // Find the L2 error to match against.
        TiXmlElement *value = metric->FirstChildElement("value");
        while (value)
        {
            // Find name of field.
            std::string variable = value->Attribute("variable");

            // Set up a match with two fields which correspond with the
            // subexpression above. The first is the variable name, second is
            // the LInf error.
            std::vector<std::string> tmp(2);
            tmp[0] = variable;
            tmp[1] = value->GetText();
            m_matches.push_back(tmp);

            // Indicate that the L2 error needs tolerance testing.
            m_tolerance.insert(1);

            value = value->NextSiblingElement("value");
        }
    }

}
