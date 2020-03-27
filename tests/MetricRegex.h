///////////////////////////////////////////////////////////////////////////////
//
// File: MetricRegex.h
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
// Description: Definition of the regular-expression metric.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_TESTS_METRICREGEX_H
#define NEKTAR_TESTS_METRICREGEX_H

#include <Metric.h>
#include <boost/regex.hpp>
#include <vector>

namespace Nektar
{
    /**
     * @brief Data structure for a Regex value to match.
     */
    struct MetricRegexFieldValue
    {
        MetricRegexFieldValue()
            : m_value(""), m_useTolerance(false), m_tolerance(0.0),
              m_useIntTolerance(false), m_intTolerance(0)
        {
        }

        std::string m_value;
        bool m_useTolerance;
        double m_tolerance;

        bool m_useIntTolerance;
        int m_intTolerance;
    };

    class MetricRegex : public Metric
    {
    public:
        virtual ~MetricRegex() {}

        static MetricSharedPtr create(TiXmlElement *metric, bool generate)
        {
            return MetricSharedPtr(new MetricRegex(metric, generate));
        }

        static std::string type;

    protected:
        /// Storage for the boost regex.
        boost::regex                                     m_regex;
        /// Stores the multiple matches defined in each <MATCH> tag.
        std::vector<std::vector<MetricRegexFieldValue> > m_matches;
        /// If true, regex matches may be in any order in output
        bool m_unordered;

        MetricRegex(TiXmlElement *metric, bool generate);

        virtual bool v_Test    (std::istream& pStdout, std::istream& pStderr);
        virtual void v_Generate(std::istream& pStdout, std::istream& pStderr);
    };
}

#endif
