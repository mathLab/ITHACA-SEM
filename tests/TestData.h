/*
 * TestData.h
 *
 *  Created on: 18 Oct 2012
 *      Author: cc
 */

#ifndef TESTDATA_H_
#define TESTDATA_H_

#include <boost/filesystem.hpp>

#include <string>
#include <vector>

#include <tinyxml/tinyxml.h>

namespace fs = boost::filesystem;

namespace Nektar
{
namespace Test
{

struct DependentFile
{
    std::string m_description;
    std::string m_filename;
};

class TestData
{
public:
    TestData(const fs::path& pFilename);
    TestData(const TestData& pSrc);

    const std::string& GetDescription() const;
    const std::string& GetExecutable() const;
    const std::string& GetParameters() const;

    std::string GetMetricType(unsigned int pId) const;
    unsigned int GetNumMetrics() const;
    TiXmlElement* GetMetric(unsigned int pId);

    DependentFile GetDependentFile(unsigned int pId) const;
    unsigned int GetNumDependentFiles() const;

private:
    std::string                     m_filename;
    std::string                     m_description;
    std::string                     m_executable;
    std::string                     m_parameters;
    TiXmlDocument*                  m_doc;
    std::vector<TiXmlElement*>      m_metrics;
    std::vector<DependentFile>      m_files;

    void Parse(TiXmlDocument* pDoc);

};

}
}

#endif /* TESTDATA_H_ */
