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

using namespace std;

namespace Nektar
{
    std::string MetricRegex::type = GetMetricFactory().
        RegisterCreatorFunction("REGEX", MetricRegex::create);
    
    /**
     * @brief Constructor.
     */
    MetricRegex::MetricRegex(TiXmlElement *metric) : Metric(metric)
    {
        // Parse a regex <METRIC> tag. This would populate m_regex, m_matches
        // and m_tolerance below.
    }

    bool MetricRegex::v_Test(std::istream& pStdout, std::istream& pStderr)
    {
        // If we have matched everything, nothing to do.
        if (m_matches.size() == 0)
        {
            return true;
        }
        
        boost::cmatch             matches;
        std::vector<std::string> &okValues = m_matches[0];

        // Process output file line by line searching for regex matches
        std::string line;
        while (getline(pStdout, line))
        {
            // Test to see if we have a match on this line.
            if (boost::regex_match(line.c_str(), matches, m_regex))
            {
                // Error if no fields in regex then throw an error.
                if (matches.size() == 1)
                {
                    cout << "No test sections in regex!" << endl;
                    return false;
                }

                // Check each field in turn
                for (int i = 1; i < matches.size(); ++i)
                {
                    std::string match(matches[i].first, matches[i].second);

                    if (m_tolerance.count(i-1) > 0)
                    {
                        // If the okValues are not within tolerance, failed the
                        // test.
                        if (fabs(boost::lexical_cast<double>(okValues[i-1]) -
                                 boost::lexical_cast<double>(match)) > 1e-6)
                        {
                            cout << "Failed tolerance match." << endl;
                            cout << "  Expected: " << okValues[i-1] << endl;
                            cout << "  Result:   " << match << endl;
                            return false;
                        }
                    }
                    else
                    {
                        // Case insensitive match.
                        if (!boost::iequals(match, okValues[i-1]))
                        {
                            cout << "Failed case-insensitive match." << endl;
                            cout << "  Expected: " << okValues[i-1] << endl;
                            cout << "  Result:   " << match << endl;
                            return false;
                        }
                    }
                }

                // Remove this match from the list of matches.
                m_matches.erase(m_matches.begin());
            }
        }

        return true;
    }
}
