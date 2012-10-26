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

#include <TestData.h>
#include <MetricL2.h>
#include <MetricRegex.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/version.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#ifdef _WINDOWS
#define COPY_COMMAND "copy "
#define MKDIR_COMMAND "mkdir "
#define DEL_COMMAND "del "
#else
#define COPY_COMMAND "cp "
#define MKDIR_COMMAND "mkdir "
#define DEL_COMMAND "rm -rf "
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
    // Parse the test file
    Test::TestData file("example.xml");

    // Generate the metric objects
    vector<MetricSharedPtr> metrics;
    for (unsigned int i = 0; i < file.GetNumMetrics(); ++i)
    {
        metrics.push_back(GetMetricFactory().CreateInstance(file.GetMetricType(i), file.GetMetric(i)));
    }

    // Copy required files for this test.
    struct stat vFileInfo;
    for (unsigned int i = 0; i < file.GetNumDependentFiles(); ++i)
    {
        Test::DependentFile f = file.GetDependentFile(i);
        string fname = file.GetDependentFile(i).m_filename;
        string source  = PortablePath(std::string(SOURCE_PATH) + "/" + fname);
        string command = std::string(COPY_COMMAND)
                        + source + " .";
        int vNotPresent = stat(source.c_str(), &vFileInfo);
        
        if (vNotPresent)
        {
            cerr << "Required file " << source << " not found." << endl;
            return 1;
        }
        
        int status = system(command.c_str());
        
        if (status)
        {
            cerr << "Unable to copy file:" << source
                 << " to current location" << endl;
            return 1;
        }
    }

    // Construct command to run
    std::string command;
    command += PortablePath(std::string(BUILD_PATH) + "/" + file.GetExecutable());
//#if defined(NDEBUG)
    command += "-g";
//#endif
    command += " " + file.GetParameters();
    command += " 1>output.out 2>output.err";

    // Run executable and perform tests.
    cout << "EXECUTING: " << command << endl;
    int status=system(command.c_str());
    
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
