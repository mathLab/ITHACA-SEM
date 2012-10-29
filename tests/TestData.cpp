/*
 * TestData.cpp
 *
 *  Created on: 18 Oct 2012
 *      Author: cc
 */

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include "TestData.h"

using namespace std;

namespace Nektar
{
namespace Test
{

TestData::TestData(const fs::path& pFilename)
{
    // Process test file format.
    m_doc = new TiXmlDocument(pFilename.c_str());
    bool loadOkay = m_doc->LoadFile();

    ASSERTL0(loadOkay, "Failed to load file: " + string(m_doc->ErrorDesc()));

    Parse(m_doc);
}

TestData::TestData(const TestData& pSrc)
{

}

const std::string& TestData::GetDescription() const
{
    return m_description;
}

const std::string& TestData::GetExecutable() const
{
    return m_executable;
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

}
}
