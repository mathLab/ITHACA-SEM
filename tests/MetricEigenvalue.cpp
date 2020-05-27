///////////////////////////////////////////////////////////////////////////////
//
// File: MetricEigenvalue.cpp
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
// Description: Implementation of the eigenvalue metric.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/algorithm/string.hpp>

#include <MetricEigenvalue.h>

namespace Nektar
{
    std::string MetricEigenvalue::type = GetMetricFactory().
        RegisterCreatorFunction("EIGENVALUE", MetricEigenvalue::create);

    // Guess default tolerance for generation routine.
    std::string MetricEigenvalue::defaultTolerance = "1e-03";

    MetricEigenvalue::MetricEigenvalue(TiXmlElement *metric, bool generate) :
        MetricRegex(metric, generate)
    {
        // We do not mind which order the converged eigenvalues are listed.
        m_unordered = true;

        // Regex for FP numbers of forms: 120, -23, 4.345, 2.4563e-01, -nan
        std::string fp = "-?\\d+\\.?\\d*(?:e[+-]\\d+)?|-?nan";

        // Set up the regular expression. This matches lines beginning with EV:
        // followed by an eigenvalue index and then at least 2 floating-point
        // values comprising real and imaginary components of complex evals.
        // Comparison is made only on the captured eigenvalue components.
        m_regex = "^EV:\\s+\\d+\\s+(" + fp + ")\\s+(" + fp + ").*";

        // Find the number of iterations to match against.
        TiXmlElement *value = metric->FirstChildElement("value");
        ASSERTL0(value || m_generate, "Missing value tag for eigenvalue metric!");

        while (value)
        {
            ASSERTL0(value->Attribute("tolerance"),
                     "Missing tolerance in eigenvalue metric");
            ASSERTL0(!EmptyString(value->GetText()),
                     "Missing value in preconditioner metric.");

            MetricRegexFieldValue mag, angle;

            // Read valute as comma-separate mag,angle parts
            std::string cmplx = value->GetText();
            std::vector<std::string> cmpts;
            boost::split(cmpts, cmplx, boost::is_any_of(","));
            ASSERTL0(cmpts.size() == 2,
                     "Value should be magnitude and angle, separated by comma");

            mag.m_value = cmpts[0];
            mag.m_useTolerance = true;
            mag.m_tolerance = atof(value->Attribute("tolerance"));

            angle.m_value = cmpts[1];
            angle.m_useTolerance = true;
            angle.m_tolerance = atof(value->Attribute("tolerance"));

            if (!m_generate)
            {
                std::vector<MetricRegexFieldValue> tmp(2);
                tmp[0] = mag;
                tmp[1] = angle;
                m_matches.push_back(tmp);
            }
            else
            {
                m_varTolerance = value->Attribute("tolerance");
            }

            value = value->NextSiblingElement("value");
        }
    }

    void MetricEigenvalue::v_Generate(std::istream& pStdout, std::istream& pStderr)
    {
        // Run MetricRegex to generate matches.
        MetricRegex::v_Generate(pStdout, pStderr);

        // First remove all existing values.
        m_metric->Clear();

        // Now create new values.
        for (int i = 0; i < m_matches.size(); ++i)
        {
            ASSERTL0(m_matches[i].size() == 3,
                     "Wrong number of matches for regular expression.");

            std::string   tol    = MetricEigenvalue::defaultTolerance;
            TiXmlElement *value  = new TiXmlElement("value");

            if (m_varTolerance != "")
            {
                tol = m_varTolerance;
            }

            value->SetAttribute("index", m_matches[i][0].m_value);
            value->SetAttribute("tolerance", tol);
            value->LinkEndChild(new TiXmlText(m_matches[i][1].m_value + ","
                                    + m_matches[i][2].m_value));
            m_metric->LinkEndChild(value);
        }
    }
}
