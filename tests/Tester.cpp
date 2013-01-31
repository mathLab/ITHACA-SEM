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
#include <boost/thread.hpp>

using namespace std;
using namespace Nektar;

// Define some namespace aliases
namespace po = boost::program_options;
namespace fs = boost::filesystem;

string PortablePath(const boost::filesystem::path& path)
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
    string command;

    // Set up command line options.
    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",                 "Produce this help message.")
        ("verbose,v",              "Turn on verbosity.")
        ("generate-metric,g",      po::value<vector<int> >(), 
                                   "Generate a single metric.")
        ("generate-all-metrics,a", "Generate all metrics.");
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file",   po::value<string>(), "Input filename");
    
    po::options_description cmdline_options;
    cmdline_options.add(hidden).add(desc);
    
    po::options_description visible("Allowed options");
    visible.add(desc);
    
    po::positional_options_description p;
    p.add("input-file", -1);
    
    po::variables_map vm;
    
    try
    {
        po::store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);
    } 
    catch (const exception& e)
    {
        cerr << e.what() << endl;
        cerr << desc;
        return 1;
    }

    if (vm.count("help") || vm.count("input-file") != 1) {
        cerr << "Usage: Tester [options] input-file.tst" << endl;
        cout << desc;
        return 1;
    }
    
    // Set up set containing metrics to be generated.
    vector<int> metricGenVec;
    if (vm.count("generate-metric"))
    {
        metricGenVec = vm["generate-metric"].as<vector<int> >();
        cout << "SIZE = " << metricGenVec.size() << endl;
    }
    set<int> metricGen(metricGenVec.begin(), metricGenVec.end());

    // Path to test definition file
    const fs::path specFile(vm["input-file"].as<string>());

    // Parent path of test definition file containing dependent files
    fs::path specPath = specFile.parent_path();

    if (specPath.empty())
    {
        specPath = fs::current_path();
    }

#if BOOST_VERSION > 104200
    string specFileStem = specFile.stem().string();
#else
    string specFileStem = specFile.stem();
#endif

    // Temporary directory to create and in which to conduct test
    const fs::path tmpDir = fs::current_path()
        / fs::path("tmp_" + specFileStem);

    // The current directory
    const fs::path startDir = fs::current_path();

    try
    {
        // Parse the test file
        TestData file(specFile);

        // Generate the metric objects
        vector<MetricSharedPtr> metrics;
        for (unsigned int i = 0; i < file.GetNumMetrics(); ++i)
        {
            set<int>::iterator it = metricGen.find(file.GetMetricId(i));
            bool genMetric = it != metricGen.end() || 
                             (vm.count("generate-all-metrics") > 0);
            
            metrics.push_back( GetMetricFactory().CreateInstance(
                                                    file.GetMetricType(i),
                                                    file.GetMetric(i),
                                                    genMetric
                                                  ));
            
            if (it != metricGen.end())
            {
                metricGen.erase(it);
            }
        }

        if (metricGen.size() != 0)
        {
            string s = metricGen.size() == 1 ? "s" : "";
            set<int>::iterator it;
            cerr << "Unable to find metric"+s+" with ID"+s+" ";
            for (it = metricGen.begin(); it != metricGen.end(); ++it)
            {
                cerr << *it << " ";
            }
            cerr << endl;
            return 1;
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
                    + boost::lexical_cast<string>(file.GetNProcesses())
                    + " ";
        }

        // If executable doesn't exist in path then hope that it is in the
        // user's PATH environment variable.
        fs::path execPath = startDir / fs::path(file.GetExecutable());
        if (!fs::exists(execPath))
        {
            execPath = fs::path(file.GetExecutable());
        }

        command += PortablePath(execPath);
        command += " ";
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
        string line;
        for (int i = 0; i < metrics.size(); ++i)
        {
            vStdout.clear();
            vStderr.clear();
            vStdout.seekg(0, ios::beg);
            vStderr.seekg(0, ios::beg);
            if (!metrics[i]->Test(vStdout, vStderr))
            {
                status = 1;
            }
        }

        // Dump output files to terminal for debugging purposes on fail.
        if (status == 1)
        {
            vStdout.clear();
            vStderr.clear();
            vStdout.seekg(0, ios::beg);
            vStderr.seekg(0, ios::beg);

            cout << "=== Output ===" << endl;
            while(vStdout.good())
            {
                getline(vStdout, line);
                cout << line << endl;
            }
            cout << "=== Errors ===" << endl;
            while(vStderr.good())
            {
                getline(vStderr, line);
                cout << line << endl;
            }
        }

        // Close output files.
        vStdout.close();
        vStderr.close();

        // Change back to the original path and delete temporary directory.
        fs::current_path(startDir);

        // Repeatedly try deleting directory with sleep for filesystems which
        // work asynchronously. This allows time for the filesystem to register
        // the output files are closed so they can be deleted and not cause a 
        // filesystem failure. Attempts made for 1 second.
        int i = 1000;
        while (i > 0)
        {
            try
            {
                // If delete successful, stop trying.
                fs::remove_all(tmpDir);
                break;
            }
            catch (const fs::filesystem_error& e)
            {
                //usleep(1000);
                boost::this_thread::sleep(boost::posix_time::milliseconds(1));
                i--;
                if (i > 0)
                {
                    cout << "Locked files encountered. "
                         << "Retring after 1ms..." << endl;
                }
                else
                {
                    // If still failing after 1sec, we consider it a permanent
                    // filesystem error and abort.
                    throw e;
                }
            }
        }
        
        // Save any changes.
        if (vm.count("generate-metric")      > 0 || 
            vm.count("generate-all-metrics") > 0)
        {
            file.SaveFile();
        }
        
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
