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
// Description: Encapsulation of test XML file.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/version.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#include <TestData.h>

using namespace std;

namespace Nektar
{
    TestData::TestData(const fs::path& pFilename)
    {
        // Process test file format.
#if BOOST_VERSION > 104200
        m_doc = new TiXmlDocument(pFilename.c_str());
#else
        m_doc = new TiXmlDocument(pFilename.file_string().c_str());
#endif
        bool loadOkay = m_doc->LoadFile();

        ASSERTL0(loadOkay, "Failed to load test definition file: "
                 + pFilename.string() + "\n"
                 + string(m_doc->ErrorDesc()));

        Parse(m_doc);
    }

    TestData::TestData(const TestData& pSrc)
    {

    }

    const std::string& TestData::GetDescription() const
    {
        return m_description;
    }

    const std::string TestData::GetExecutable() const
    {
        std::string execname = m_executable;
    #if defined(RELWITHDEBINFO)
        execname += "-rg";
    #elif !defined(NDEBUG)
        execname += "-g";
    #endif
        return execname;
    }

    const std::string& TestData::GetParameters() const
    {
        return m_parameters;
    }

    const unsigned int& TestData::GetNProcesses() const
    {
        return m_processes;
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

    void TestData::Parse(TiXmlDocument* pDoc)
    {
        TiXmlHandle handle(m_doc);
        TiXmlElement *testElement, *tmp, *metrics, *files;
        testElement = handle.FirstChildElement("test").Element();
        ASSERTL0(testElement, "Cannot find 'test' root element.");

        // Find description tag.
        tmp = testElement->FirstChildElement("description");
        ASSERTL0(tmp, "Cannot find 'description' for test.");
        m_description = string(tmp->GetText());

        // Find executable tag.
        tmp = testElement->FirstChildElement("executable");
        ASSERTL0(tmp, "Cannot find 'executable' for test.");
        m_executable = string(tmp->GetText());

        // Find parameters tag.
        tmp = testElement->FirstChildElement("parameters");
        ASSERTL0(tmp, "Cannot find 'parameters' for test.");
        m_parameters = string(tmp->GetText());

        // Find parallel processes tah.
        tmp = testElement->FirstChildElement("processes");
        if (tmp)
        {
            m_processes = atoi(tmp->GetText());
        }
        else
        {
            m_processes = 1;
        }

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
