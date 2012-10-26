///////////////////////////////////////////////////////////////////////////////
//
// File: MetricRegex.cpp
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
// Description: Implementation of the regular-expression metric.
//
///////////////////////////////////////////////////////////////////////////////

#include <MetricRegex.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace Nektar
{
    std::string MetricRegex::type = GetMetricFactory().
        RegisterCreatorFunction("regex", MetricRegex::create);
    
    /**
     * @brief Constructor.
     */
    MetricRegex::MetricRegex(TiXmlElement *metric) : Metric(metric)
    {
    
    }

    void MetricRegex::v_Parse(TiXmlElement *metric)
    {
        // Parse a regex <METRIC> tag. This would populate m_regex, m_matches
        // and m_tolerance below.
    }

    bool MetricRegex::v_TestLine(std::string line)
    {
        // If we have matched everything, nothing to do.
        if (m_matches.size() == 0)
        {
            return true;
        }
        
        boost::cmatch             matches;
        std::vector<std::string> &okValues = m_matches[0];

        // Test to see if we have a match.
        if (boost::regex_match(line.c_str(), matches, m_regex))
        {
            for (int i = 0; i < matches.size(); ++i)
            {
                std::string match(matches[i].first, matches[i].second);

                if (m_tolerance.count(i) > 0)
                {
                    // If the okValues are not within tolerance, failed the
                    // test.
                    if (fabs(boost::lexical_cast<int>(okValues[i]) -
                             boost::lexical_cast<int>(match)) > 1e-6)
                    {
                        return false;
                    }
                }
                else
                {
                    // Case insensitive match.
                    if (!boost::iequals(match, okValues[i]))
                    {
                        return false;
                    }
                }
            }
            
            // Remove this match from the list of matches.
            m_matches.erase(m_matches.begin());
        }
    }
}
