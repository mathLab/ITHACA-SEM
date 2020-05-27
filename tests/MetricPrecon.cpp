///////////////////////////////////////////////////////////////////////////////
//
// File: MetricPrecon.cpp
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
// Description: Implementation of the preconditioner metric.
//
///////////////////////////////////////////////////////////////////////////////

#include <MetricPrecon.h>

namespace Nektar
{
    std::string MetricPrecon::type = GetMetricFactory().
        RegisterCreatorFunction("PRECON", MetricPrecon::create);

    // Guess default tolerance for generation routine.
    std::string MetricPrecon::defaultTolerance = "2";

    MetricPrecon::MetricPrecon(TiXmlElement *metric, bool generate) :
        MetricRegex(metric, generate)
    {
        // Set up the regular expression. This (optionally) matches a variable
        // name if it exists: first field is variable name, second field is L2
        // error.
        m_regex = "^CG iterations made = (\\d+) .*";

        // Find the number of iterations to match against.
        TiXmlElement *value = metric->FirstChildElement("value");
        ASSERTL0(value || m_generate, "Missing value tag for precon metric!");

        while (value)
        {
            // Set up a match with two fields which correspond with the
            // subexpression above. The first is the variable name, second is
            // the precon iteration count.
            ASSERTL0(value->Attribute("tolerance"),
                     "Missing tolerance in preconditioner metric");
            ASSERTL0(!EmptyString(value->GetText()),
                     "Missing value in preconditioner metric.");

            MetricRegexFieldValue val;
            val.m_value = value->GetText();
            val.m_useIntTolerance = true;
            val.m_intTolerance = atoi(value->Attribute("tolerance"));

            if (!m_generate)
            {
                std::vector<MetricRegexFieldValue> tmp(1);
                tmp[0] = val;
                m_matches.push_back(tmp);
            }
            else
            {
                m_varTolerance = value->Attribute("tolerance");
            }

            value = value->NextSiblingElement("value");
        }
    }

    void MetricPrecon::v_Generate(std::istream& pStdout, std::istream& pStderr)
    {
        // Run MetricRegex to generate matches.
        MetricRegex::v_Generate(pStdout, pStderr);

        // First remove all existing values.
        m_metric->Clear();

        // Now create new values.
        for (int i = 0; i < m_matches.size(); ++i)
        {
            ASSERTL0(m_matches[i].size() == 1,
                     "Wrong number of matches for regular expression.");

            std::string   tol    = MetricPrecon::defaultTolerance;
            TiXmlElement *value  = new TiXmlElement("value");

            if (m_varTolerance != "")
            {
                tol = m_varTolerance;
            }

            value->SetAttribute("tolerance", tol);
            value->LinkEndChild(new TiXmlText(m_matches[i][0].m_value));
            m_metric->LinkEndChild(value);
        }
    }
}
