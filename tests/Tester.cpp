///////////////////////////////////////////////////////////////////////////////
//
// File: Tester.cpp
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
// Description: Tester executible.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>

#include <MetricL2.h>
#include <MetricRegex.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/version.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#ifdef _WINDOWS
#define COPY_COMMAND "copy "
#else
#define COPY_COMMAND "cp "
#endif

using namespace std;
using namespace Nektar;
namespace po = boost::program_options;

std::string PortablePath(const boost::filesystem::path& path)
{
    boost::filesystem::path temp = path;
#if BOOST_VERSION > 104200
    temp.make_preferred();
    return temp.string();
#else
    return temp.file_string();
#endif
}

int main(int argc, char *argv[])
{
    // Process test file format.
    TiXmlDocument *file = new TiXmlDocument(argv[1]);
    bool loadOkay = file->LoadFile();
    
    if (!loadOkay)
    {
        cout << file->ErrorDesc() << endl;
        return 1;
    }
    
    TiXmlHandle handle(file);
    TiXmlElement *testElement;
    testElement = handle.FirstChildElement("test").Element();
    
    // Find description tag.
    TiXmlElement *descElement = handle.FirstChildElement("test").
        FirstChildElement("description").Element();
    string description(descElement->GetText());
    
    // Find solver tag.
    TiXmlElement *solverElement = handle.FirstChildElement("test").
        FirstChildElement("executable").Element();
    string executable(solverElement->GetText());
    
    // Find input tag.
    TiXmlElement *inputElement = handle.FirstChildElement("test").
        FirstChildElement("parameters").Element();
    string input(inputElement->GetText());

    /*
     * Process and set up metric tags.
     */
    TiXmlElement *metrics = testElement->FirstChildElement("metrics");
    
    if (!metrics)
    {
        cout << "Couldn't find metrics tag." << endl;
        return 1;
    }
    
    // Construct vector of metrics.
    TiXmlElement *metric = metrics->FirstChildElement("metric");
    vector<MetricSharedPtr> metricList;
    
    while (metric)
    {
        string  metricType = metric->Attribute("type");
        int     metricId   = boost::lexical_cast<int>(metric->Attribute("id"));

        boost::algorithm::to_lower(metricType);
        
        MetricSharedPtr m = GetMetricFactory().CreateInstance(
            metricType, metricId);
        metricList.push_back(m);
        
        cout << "Metric type: " << metricType << ", id = " << metricId << endl;
        
        m->Parse(metric);
        
        metric = metric->NextSiblingElement("metric");
    }

    /*
     * Copy required files for this test.
     */
    TiXmlElement *required = handle.FirstChildElement("test").
        FirstChildElement("required").Element();
    TiXmlElement *reqfile = required->FirstChildElement("file");
    struct stat vFileInfo;
    
    while (reqfile)
    {
        string fname   = reqfile->GetText();
        string source  = PortablePath(std::string(TEST_PATH)) + fname;
        string command = std::string(COPY_COMMAND) + source + " .";
        int vNotPresent = stat(source.c_str(), &vFileInfo);
        
        if (vNotPresent)
        {
            cerr << "Required file " << fname << " not found." << endl;
            return 1;
        }
        
        int status = system(command.c_str());
        
        if (status)
        {
            cerr << "Unable to copy file:" << source 
                 << " to current location" << endl;
            return 1;
        }
        
        reqfile = metric->NextSiblingElement("required");
    }
    
    /*
     * Run executable and perform tests.
     */
    string cmd = executable + " " + input;
    
    
    /*
    ifstream file(filename);
        
    while (getline(file, line))
    {
        for (int i = 0; i < metricList.size(); ++i)
        {
            if (!metricList[i]->TestLine(line))
            {
                return 1;
            }
        }
    }

    for (int i = 0; i < metricList.size(); ++i)
    {
        if (!metricList[i]->FinishTest())
        {
            return 1;
        }
    }
    */
    
    return 0;
}
