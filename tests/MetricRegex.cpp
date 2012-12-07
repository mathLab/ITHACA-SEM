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
    MetricRegex::MetricRegex(TiXmlElement *metric, bool generate) : 
        Metric(metric, generate)
    {
        // If we are a derived class, do nothing
        if (m_type != "REGEX")
        {
            return;
        }

        // Parse Regex expression
        TiXmlElement *regex = metric->FirstChildElement("regex");
        ASSERTL0(regex, "No Regex defined.");
        ASSERTL0(regex->GetText(), "Failed to get text");
        m_regex = regex->GetText();

        // Parse matching values if not generating.
        if (m_generate)
        {
            return;
        }
        
        TiXmlElement *matches = metric->FirstChildElement("matches");
        ASSERTL0(matches, "No matches defined.");
        TiXmlElement *match = matches->FirstChildElement("match");
        while (match)
        {
            std::vector<MetricRegexFieldValue> tmp;
            TiXmlElement *field = match->FirstChildElement("field");
            while (field)
            {
                MetricRegexFieldValue v;
                v.m_value = field->GetText();
                v.m_useTolerance = false;

                const char * tol = field->Attribute("tolerance");
                if (tol)
                {
                    v.m_useTolerance = true;
                    v.m_tolerance = atof(tol);
                }
                tmp.push_back(v);
                field = field->NextSiblingElement("field");
            }
            m_matches.push_back(tmp);

            match = match->NextSiblingElement("match");
        }
    }

    /**
     * @brief Test output against a regex expression and set of matches.
     */
    bool MetricRegex::v_Test(std::istream& pStdout, std::istream& pStderr)
    {
        ASSERTL0(m_matches.size(), "No test conditions defined for Regex.");

        std::vector<MetricRegexFieldValue> &okValues = m_matches[0];
        int                                 nMatch   = m_matches.size();
        bool                                success  = true;
        boost::cmatch                       matches;

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
                    cerr << "No test sections in regex!" << endl;
                    return false;
                }

                // Check each field in turn
                for (int i = 1; i < matches.size(); ++i)
                {
                    std::string match(matches[i].first, matches[i].second);
                    
                    if (okValues[i-1].m_useTolerance)
                    {
                        double val;
                        try 
                        {
                            val = fabs(fabs(boost::lexical_cast<double>(
                                           okValues[i-1].m_value))
                                     - fabs(boost::lexical_cast<double>(match))
                                      );
                        }
                        catch(boost::bad_lexical_cast &e)
                        {
                            cerr << "Could not convert one of " << match 
                                 << " (match) or " << okValues[i-1].m_value
                                 << " (comparison value) to double" << endl;
                            success = false;
                            continue;
                        }
                        
                        // If the okValues are not within tolerance, failed the
                        // test.
                        if (val > okValues[i-1].m_tolerance)
                        {
                            cerr << "Failed tolerance match." << endl;
                            cerr << "  Expected: " << okValues[i-1].m_value
                                 << " +/- " << okValues[i-1].m_tolerance
                                 << endl;
                            cerr << "  Result:   " << match << endl;
                            success = false;
                        }
                    }
                    else
                    {
                        // Case insensitive match.
                        if (!boost::iequals(match, okValues[i-1].m_value))
                        {
                            cerr << "Failed case-insensitive match." << endl;
                            cerr << "  Expected: " << okValues[i-1].m_value
                                 << endl;
                            cerr << "  Result:   " << match << endl;
                            success = false;
                        }
                    }
                }

                // Remove this match from the list of matches.
                m_matches.erase(m_matches.begin());
            }
        }

        if (m_matches.size() != 0)
        {
            cerr << "Expected " << nMatch << " matches but only found "
                 << (nMatch - m_matches.size()) << "!" << endl;
            success = false;
        }

        return success;
    }

    /**
     * @brief Test output against a regex expression and set of matches.
     */
    void MetricRegex::v_Generate(std::istream& pStdout, std::istream& pStderr)
    {
        boost::cmatch matches;

        // Process output file line by line searching for regex matches
        std::string line;
        while (getline(pStdout, line))
        {
            // Test to see if we have a match on this line.
            if (boost::regex_match(line.c_str(), matches, m_regex))
            {
                // Error if no fields in regex then throw an error.
                ASSERTL0(matches.size() != 1, "No test sections in regex!");
                
                vector<MetricRegexFieldValue> okValues;

                for (int i = 1; i < matches.size(); ++i)
                {
                    // Create new field.
                    MetricRegexFieldValue okValue;
                    okValue.m_useTolerance = false;
                    okValue.m_value        = std::string(matches[i].first, 
                                                         matches[i].second);
                    okValues.push_back(okValue);
                }
                
                m_matches.push_back(okValues);
            }
        }

        // If we are not a derived class then create a new structure.
        if (m_type == "REGEX")
        {
            // Remove matches if they already exist.
            TiXmlElement *matches = m_metric->FirstChildElement("matches");
            if (matches)
            {
                ASSERTL0(m_metric->RemoveChild(matches), 
                         "Couldn't remove matches from metric!");
            }

            // Create new matches element.
            matches = new TiXmlElement("matches");
            m_metric->LinkEndChild(matches);
            
            for (int i = 0; i < m_matches.size(); ++i)
            {
                TiXmlElement *match = new TiXmlElement("match");
                matches->LinkEndChild(match);
                
                for (int j = 0; j < m_matches[i].size(); ++j)
                {
                    TiXmlElement *field = new TiXmlElement("field");
                    match->LinkEndChild(field);
                    
                    field->SetAttribute(
                        "id", boost::lexical_cast<std::string>(j));
                    
                    if (m_matches[i][j].m_useTolerance)
                    {
                        field->SetAttribute(
                            "tolerance", boost::lexical_cast<
                                std::string>(m_matches[i][j].m_tolerance));
                    }
                    
                    field->LinkEndChild(new TiXmlText(m_matches[i][j].m_value));
                }
            }
        }
    }
}
