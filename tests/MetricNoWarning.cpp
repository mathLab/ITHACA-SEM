///////////////////////////////////////////////////////////////////////////////
//
// File: MetricNoWarning.cpp
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
// Description: 
// Implementation of the no-warning metric. This test metric can be used
// in either of two ways: 
//
// 1. Default mode: Test fails if "WARNING" appears in output OR error stream
//      
//      <metric type="nowarning" id="1">
//
// 2. Advanced mode: Test fails if the specified regex expression appears in
//    the output OR the error stream, AND the regex groups (caught with the
//    parenthesises) match all the expressions specified within ONE of the 
//    match tags (one expression = one field tag). The example here contains 
//    two regex groups, and two possible matches. The test will therefore fail
//    if, e.g.,
//              "WARNING: Invalid value"
//    or
//              "WARNING: Invalid input value" 
//    appears in the output OR the error stream. A simple "WARNING" will 
//    however not fail the test. The advanced mode option is therefore 
//    intended to target very  specific warnings/errors. Note that a regex
//    without specific regex groups can also be specified.
//
//      <metric type="nowarning" id="1">
//          <regex>.*WARNING:.(\w+).(\w+).*</regex>
//          <matches>
//              <match>
//                  <field id="1">Invalid</field>
//                  <field id="2">value</field>
//              </match>
//              <match>
//                  <field id="1">Invalid</field>
//                  <field id="2">input</field>
//              </match>
//          </matches>
//      </metric>
//
//
///////////////////////////////////////////////////////////////////////////////

#include <MetricNoWarning.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>

namespace Nektar
{
    std::string MetricNoWarning::type = GetMetricFactory().
        RegisterCreatorFunction("NOWARNING", MetricNoWarning::create);
    
    /**
     * @brief Constructor
     */
    MetricNoWarning::MetricNoWarning(TiXmlElement *metric, bool generate) : 
        Metric(metric, generate)
    {
        // Check if user has provided custom warning regex, which can contain
        // several matching groups to target specific Warning messages
        TiXmlElement *regex = metric->FirstChildElement("regex");
        if (regex)
        {
            // Check that user provided a regex for this xml tag
            const char *tmp = regex->GetText();
            ASSERTL0(tmp, "No text found in regex tag");
            m_regexWarning = tmp;
        }

        // If we generate from output, we don't need to parse the xml file
        if (m_generate)
        {
            return;
        }

        // Check if user has provided custom matching groups for regex
        // Note that this is optional, a regex without matching groups is
        // also permissible
        TiXmlElement *matches = metric->FirstChildElement("matches");
        if (matches)
        {
            // Go through all "match" tags
            for (
                TiXmlElement *match = matches->FirstChildElement("match");
                match;
                match = match->NextSiblingElement("match")
                )
            {
                // Add vector with regex groups
                m_matches.push_back(std::vector<std::string>());

                // Go though all "field" tags
                for (
                    TiXmlElement *field = match->FirstChildElement("field");
                    field; 
                    field = field->NextSiblingElement("field"))
                    {
                        // Extract field text (i.e. regex group)
                        const char *tmp = field->GetText();
                        // Check that user provided a text
                        ASSERTL0(tmp, "No text that specifies regex group "
                            "found in field tag");
                        m_matches.back().push_back(tmp);
                    }
            }

            // Check if user provided any match tags
            ASSERTL0(m_matches.size(), "No match tag that specifies "
                "regex groups was found inside matches tag.");

            // Check that all set of regex groups are the same size,
            // i.e. contains the same number of "field" tags
            int size_min = 1E3;
            int size_max = 0;
            for (const auto &match : m_matches)
            {
                size_min = std::min((int) match.size(), size_min);
                size_max = std::max((int) match.size(), size_max);
            }
            ASSERTL0(size_min!=0, "No valid field tags found "
                "for one of the match tags");
            ASSERTL0(size_min == size_max, "Number of valid field tags "
                "not the same for all match tags");
        }
    }

    /**
     * @brief Test if the output contains the warning message.
     *        If so, the test fails, otherwise, it passes
     */
    bool MetricNoWarning::v_Test(std::istream& pStdout, std::istream& pStderr)
    {
        boost::ignore_unused(pStdout, pStderr);

        // Loop over both standard output and error output
        for (
            std::string line;
            getline(pStdout, line)||getline(pStderr, line);
            )
        {
            boost::smatch matches;

            // Check if regex match against given output
            if (boost::regex_search(line, matches, m_regexWarning))
            {
                // Test fails if regex matches line and 
                // no matching groups have been specified
                if (!m_matches.size())
                {
                    return false;
                }
                // If regex groups are specified, the test fails if any of
                // them matches the groups in "matches"
                for (const auto &match : m_matches)
                {
                    bool all_match = false;
                    for (
                        int i=1;
                        i<matches.size() && matches.size() == match.size()+1;
                        ++i)
                    {
                        // Construct string object that contains submatch
                        std::string submatch(matches[i].first,
                            matches[i].second);
                        
                        // Compare to specified pattern
                        if(submatch == match[i-1])
                        {
                            all_match = true;
                            continue;
                        }
                        else
                        {
                            all_match = false;
                            break;
                        }
                    }
                    if (all_match)
                    {
                        return false;
                    }
                    else
                    {
                        continue;
                    }
                    
                }
            }
        }

        // If we arrived here, the test passed
        return true;

    }

    /**
     * @brief Test if the output contains the warning message.
     *        If so, generate the .tst file
     */
    void MetricNoWarning::v_Generate(std::istream& pStdout, std::istream& pStderr)
    {
        boost::ignore_unused(pStderr);

        // Check both standard output and error output
        for (
            std::string line;
            getline(pStdout, line)||getline(pStderr, line);
            )
        {
            boost::smatch matches;

            // Check if regex match against given output
            if (boost::regex_search(line, matches, m_regexWarning))
            {
                // If regex groups were not found, continue
                if(matches.size() == 1)
                {
                    continue;
                }
                
                // Add vector with regex groups
                m_matches.push_back(std::vector<std::string>());

                // Save all regex groups
                for (int i=1; i<matches.size(); ++i)
                {
                    // Construct string object that contains submatch
                    std::string submatch(matches[i].first,matches[i].second);

                    m_matches.back().push_back(submatch);
                }
            }
        }

        // Remove matches if they already exist.
        TiXmlElement *matches = m_metric->FirstChildElement("matches");
        if (matches)
        {
            ASSERTL0(m_metric->RemoveChild(matches), "Couldn't remove matches "
                "from metric")
        }

        // Add new "matches" tag
        matches = new TiXmlElement("matches");
        m_metric->LinkEndChild(matches);

        // Add all "match" tags under "matches" tag
        for (const auto &match_it : m_matches)
        {
            TiXmlElement *match = new TiXmlElement("match");
            matches->LinkEndChild(match);

            int j = 0;
            for (const auto &field_it : match_it)
            {
                TiXmlElement *field = new TiXmlElement("field");
                match->LinkEndChild(field);

                field->SetAttribute(
                    "id", boost::lexical_cast<std::string>(j++));

                field->LinkEndChild(new TiXmlText(field_it));
            }
        }
    }

} // namespace Nektar
