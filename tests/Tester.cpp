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
// Description: Tester executable.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <TestData.h>
#include <Metric.h>

#include <boost/filesystem.hpp>
#include <boost/version.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace Nektar;

// Define some namespace aliases
namespace po = boost::program_options;
namespace fs = boost::filesystem;

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
    int status = 0;
    std::string command;

    if (argc != 2)
    {
        cerr << "Error: incorrect number of arguments" << endl;
        cerr << "Usage: " << argv[0] << " [test-definition.tst]" << endl;
        return 2;
    }

    // Path to test definition file
    const fs::path specFile(argv[1]);

    // Parent path of test definition file containing dependent files
    const fs::path specPath = specFile.parent_path();

    // Temporary directory to create and in which to conduct test
    const fs::path tmpDir = fs::current_path()
                            / fs::path("tmp_" + specFile.stem().string());

    // The current directory
    const fs::path startDir = fs::current_path();

    try
    {
        // Parse the test file
        Test::TestData file(specFile);

        // Generate the metric objects
        vector<MetricSharedPtr> metrics;
        for (unsigned int i = 0; i < file.GetNumMetrics(); ++i)
        {
            metrics.push_back( GetMetricFactory().CreateInstance(
                                                    file.GetMetricType(i),
                                                    file.GetMetric(i)
                                                  ));
        }

        // Remove the temporary directory if left from a previous test
        if (fs::exists(tmpDir))
        {
            fs::remove_all(tmpDir);
        }

        // Create temporary directory
        fs::create_directory(tmpDir);

        // Change working directory to the temporary directory
        fs::current_path(tmpDir);

        // Copy required files for this test from the test definition directory
        // to the temporary directory.
        for (unsigned int i = 0; i < file.GetNumDependentFiles(); ++i)
        {
            fs::path source_file(file.GetDependentFile(i).m_filename);

            fs::path source = specPath / source_file;
            fs::path dest   = tmpDir   / source_file;
            fs::copy_file(source, dest);
        }

        // Construct test command to run. If in debug mode, append "-g"
        // Output from stdout and stderr are directed to the files output.out
        // and output.err, respectively.
        if (file.GetNProcesses() > 1)
        {
            command += "mpirun -np "
                    + boost::lexical_cast<std::string>(file.GetNProcesses())
                    + " ";
        }
        fs::path execPath = startDir / fs::path(file.GetExecutable());
        command += PortablePath(execPath);
    #if defined(NDEBUG)
        command += " ";
    #else
        command += "-g ";
    #endif
        command += file.GetParameters();
        command += " 1>output.out 2>output.err";

        // Run executable to perform test.
        if (system(command.c_str()))
        {
            cerr << "Error occurred running test:" << endl;
            cerr << "Command: " << command << endl;
            throw 1;
        }

        // Check output files exist
        if (!(fs::exists("output.out") && fs::exists("output.err")))
        {
            cerr << "One or more test output files are missing." << endl;
            throw 1;
        }

        // Open output files and check they are readable
        ifstream vStdout("output.out");
        ifstream vStderr("output.err");
        if (vStdout.bad() || vStderr.bad())
        {
            cerr << "One or more test output files are unreadable." << endl;
            throw 1;
        }

        // Test against all metrics
        status = 0;
        for (int i = 0; i < metrics.size(); ++i)
        {
            vStdout.seekg(0, ios::beg);
            vStderr.seekg(0, ios::beg);
            if (!metrics[i]->Test(vStdout, vStderr))
            {
                status = 1;
            }
        }

        // Change back to the original path and delete temporary directory
        fs::current_path(startDir);
        fs::remove_all(tmpDir);

        // Return status of test. 0 = PASS, 1 = FAIL
        return status;
    }
    catch (const fs::filesystem_error& e)
    {
        cerr << "Filesystem operation error occurred:" << endl;
        cerr << "  " << e.what() << endl;
        cerr << "  Files left in " << tmpDir.string() << endl;
    }
    catch (...)
    {
        cerr << "  Files left in " << tmpDir.string() << endl;
    }

    // If a system error, return 2
    return 2;
}
