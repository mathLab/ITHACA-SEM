///////////////////////////////////////////////////////////////////////////////
//
// File: TestData.cpp
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
// Description: Encapsulation of test XML file.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <TestData.h>
#include <TestException.hpp>

using namespace std;

namespace Nektar
{
    TestData::TestData(const fs::path& pFilename, po::variables_map& pVm)
            : m_cmdoptions(pVm)
    {
        // Process test file format.
        m_doc = new TiXmlDocument(pFilename.string().c_str());

        bool loadOkay = m_doc->LoadFile();

        ASSERTL0(loadOkay, "Failed to load test definition file: "
                 + pFilename.string() + "\n"
                 + string(m_doc->ErrorDesc()));

        Parse(m_doc);
    }

    TestData::TestData(const TestData& pSrc)
    {
        boost::ignore_unused(pSrc);
    }

    const std::string& TestData::GetDescription() const
    {
        return m_description;
    }

    const Command& TestData::GetCommand(unsigned int pId) const
    {
        ASSERTL0(pId < m_commands.size(),
                 "Command ID '" + std::to_string(pId) + "' not found");
        return m_commands[pId];
    }

    unsigned int TestData::GetNumCommands() const
    {
        return m_commands.size();
    }

    std::string TestData::GetMetricType(unsigned int pId) const
    {
        ASSERTL0(pId < m_metrics.size(), "Metric ID out of range.");

        // read the property name
        ASSERTL0(m_metrics[pId]->Attribute("type"),
                 "Missing 'type' attribute in metric "
                 + boost::lexical_cast<string>(pId) + ").");
        return boost::to_upper_copy(string(m_metrics[pId]->Attribute("type")));
    }

    unsigned int TestData::GetNumMetrics() const
    {
        return m_metrics.size();
    }

    TiXmlElement* TestData::GetMetric(unsigned int pId)
    {
        ASSERTL0(pId < m_metrics.size(), "Metric index out of range.");
        return m_metrics[pId];
    }

    unsigned int TestData::GetMetricId(unsigned int pId)
    {
        ASSERTL0(pId < m_metrics.size(), "Metric index out of range.");
        const char *id = m_metrics[pId]->Attribute("id");
        ASSERTL0(id, "No ID found for metric!");
        return boost::lexical_cast<unsigned int>(id);
    }

    DependentFile TestData::GetDependentFile(unsigned int pId) const
    {
        ASSERTL0(pId < m_files.size(), "File index out of range.");
        return m_files[pId];
    }

    unsigned int TestData::GetNumDependentFiles() const
    {
        return m_files.size();
    }

    Command TestData::ParseCommand(TiXmlElement *elmt) const
    {
        Command cmd;
        TiXmlElement *tmp;

        cmd.m_pythonTest = false;

        // Parse executable tag. Do not enforce a check because this might be
        // overridden on the command line.
        if (elmt->FirstChildElement("executable"))
        {
            tmp = elmt->FirstChildElement("executable");
            cmd.m_executable = fs::path(tmp->GetText());

            // Test to see if this test requires Python
            std::string needsPython;
            tmp->QueryStringAttribute("python", &needsPython);
            cmd.m_pythonTest = needsPython == "true";

#if defined(RELWITHDEBINFO)
            cmd.m_executable += cmd.m_pythonTest ? "" : "-rg";
#elif !defined(NDEBUG)
            cmd.m_executable += cmd.m_pythonTest ? "" : "-g";
#endif
        }

        // Find associated parameters.
        tmp = elmt->FirstChildElement("parameters");
        ASSERTL0(tmp, "Cannot find 'parameters' for test.");
        if (tmp->GetText())
        {
            cmd.m_parameters = string(tmp->GetText());
        }

        // Find parallel processes tah.
        tmp = elmt->FirstChildElement("processes");
        if (tmp)
        {
            cmd.m_processes = atoi(tmp->GetText());
        }
        else
        {
            cmd.m_processes = 1;
        }

        return cmd;
    }

    void TestData::Parse(TiXmlDocument* pDoc)
    {
        TiXmlHandle handle(pDoc);
        TiXmlElement *testElement, *tmp, *metrics, *files;
        testElement = handle.FirstChildElement("test").Element();
        ASSERTL0(testElement, "Cannot find 'test' root element.");

        // Find description tag.
        tmp = testElement->FirstChildElement("description");
        ASSERTL0(tmp, "Cannot find 'description' for test.");
        m_description = string(tmp->GetText());

        // Find command(s) to run.
        if (m_cmdoptions.count("executable"))
        {
            m_commands.push_back(ParseCommand(testElement));
            m_commands.back().m_executable = fs::path(
                m_cmdoptions["executable"].as<std::string>());
        }
        else if (testElement->FirstChildElement("executable"))
        {
            m_commands.push_back(ParseCommand(testElement));
        }
        else if ((tmp = testElement->FirstChildElement("segment")))
        {
            ASSERTL0(m_cmdoptions.count("executable") == 0,
                     "Test files defining more than one command in segment "
                     "blocks cannot use --executable.");

            while (tmp)
            {
                m_commands.push_back(ParseCommand(tmp));
                tmp = tmp->NextSiblingElement("segment");
            }
        }

        ASSERTL0(m_commands.size() > 0,
                 "No executable / command segments specified for test.");

        // Extract metric tags
        metrics = testElement->FirstChildElement("metrics");
        ASSERTL0(metrics, "No metrics defined for test.");

        tmp = metrics->FirstChildElement("metric");
        while (tmp)
        {
            m_metrics.push_back(tmp);
            tmp = tmp->NextSiblingElement("metric");
        }

        // Extract list of dependent files
        files = testElement->FirstChildElement("files");
        if (files)
        {
            tmp = files->FirstChildElement("file");
            while (tmp)
            {
                DependentFile f;
                f.m_filename = string(tmp->GetText());
                if (tmp->Attribute("description"))
                {
                    f.m_description = string(tmp->Attribute("description"));
                }
                m_files.push_back(f);
                tmp = tmp->NextSiblingElement("file");
            }
        }
    }

    void TestData::SaveFile()
    {
        m_doc->SaveFile();
    }
}
