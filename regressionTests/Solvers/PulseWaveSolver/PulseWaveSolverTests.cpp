///////////////////////////////////////////////////////////////////////////////
//
// File PulseWaveSolverTests.cpp
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
// Description: Run a series of tests on the PulseWaveSolver
//
///////////////////////////////////////////////////////////////////////////////
#include "../../Auxiliary/RegressBase.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <boost/filesystem.hpp>

void RunL2RegressionTest(std::string demo, std::string input, std::string info);
void MakeOkFile(std::string demo, std::string input, std::string info);

#ifdef MAKE_OK_FILE
#define Execute($1,$2,$3)  MakeOkFile($1,$2,$3)
#else
#define Execute($1,$2,$3)  RunL2RegressionTest($1,$2,$3)
#endif

#ifdef _WINDOWS
#define COPY_COMMAND "copy "

#else
#define COPY_COMMAND "cp "
#endif

static int tests_total = 0;
static int tests_passed = 0;
static int tests_failed = 0;
static bool verbose = false;
static bool quiet = false;

int main(int argc, char* argv[])
{
    if (argc > 2)
    {
        std::cout << "Usage: " << argv[0] << " [-v|-q]" << std::endl;
        return -1;
    }
    if (argc == 2 && std::string(argv[1]) == "-v")
    {
        verbose = true;
    }
    if (argc == 2 && std::string(argv[1]) == "-q")
    {
        quiet = true;
    }
   
    // Test of linear SWE 
    Execute("PulseWaveSolver","TestPaper1_DG.xml","1D discontinuous material properties, P=5");
	Execute("PulseWaveSolver","TestPaper2_DG.xml","1D straight vessel, P=5");
	Execute("PulseWaveSolver","TestPaper3_DG.xml","1D tapering vessel, P=5");
	Execute("PulseWaveSolver","Test_Q_inflow.xml","Q-inflow Riemann solver, P=5");
	Execute("PulseWaveSolver","Test_Junction.xml","Junction Riemann solver, P=5");
	Execute("PulseWaveSolver","Test_Bifurcation.xml","Bifurcation Riemann solver, P=5");
	Execute("PulseWaveSolver","Test_TwoBifurcations.xml","Double Bifurcation, P=11");
	Execute("PulseWaveSolver","TestPaper4_55Network.xml","55 Artery Network, P=5");

    
    if (tests_failed && !quiet)
    {
        std::cout << "WARNING: " << tests_failed << " test(s) failed." << std::endl;
    }
    else if (verbose)
    {
        std::cout << "All tests passed successfully!" << std::endl;
    }
    return tests_failed;
}

void RunL2RegressionTest(std::string Demo, std::string input, std::string info)
{
    tests_total++;
    if (!quiet)
    {
        std::cout << "TESTING: " << Demo << std::flush;
    }
	std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "/bin/";
    RegressBase Test(NektarSolverDir.c_str(),Demo,input,"Solvers/PulseWaveSolver/OkFiles/");
    int fail;

    // Copy input file to current location
	boost::filesystem::path sourceFile(std::string(REG_PATH) + "Solvers/PulseWaveSolver/InputFiles/" + input);
    std::string syscommand = std::string(COPY_COMMAND) + PortablePath(sourceFile) + " .";
    int status = system(syscommand.c_str());
    if(status)
    {
        std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
        exit(2);
    }

    if(fail = Test.TestL2()) // test failed
    {
        if (!quiet)
        {
            std::cout << "\rFAILED: " << Demo << " " << input << " (" << info << ")" << std::endl;
            std::cout << "===========================================================\n";
            // Explain cause of error if available
            Test.PrintTestError(fail);
            std::cout << "===========================================================\n";
        }
        tests_failed++;
    }
    else
    {
        if (verbose)
        {
            std::cout << "\rPASSED: " << Demo << " " << input << " (" << info << ")" << std::endl;
        }
        else if (quiet)
        {
            // print nothing
        }
        else {
            std::cout << "\rPASSED: " << Demo << " (" << info << ")" << std::endl;
        }
        tests_passed++;
    }

#ifdef _WINDOWS
	std::string cleanup = "del /Q *.xml *.fld *.chk *.his";
#else
    std::string cleanup = "rm -f *.xml *.fld *.chk *.his";
#endif
    system(cleanup.c_str());
};

void MakeOkFile(std::string Demo, std::string input, std::string info)
{
    tests_total++;
    std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "/bin/";
    RegressBase Test(NektarSolverDir,Demo,input,"Solvers/PulseWaveSolver/OkFiles/");
    int fail;

    // Copy input file to current location
	boost::filesystem::path sourceFile(std::string(REG_PATH) + "Solvers/PulseWaveSolver/InputFiles/" + input);
    std::string syscommand = std::string(COPY_COMMAND) + PortablePath(sourceFile) +" .";
    int status = system(syscommand.c_str());
    if(status)
    {
        std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
        exit(2);
    }

    if(fail = Test.MakeOkFile())
    {
        std::cout << "Failed to make OK file\n";
        // Explain cause of error if available
        std::cout << "===========================================================\n";
        Test.PrintTestError(fail);
        std::cout << "===========================================================\n";
        tests_failed++;
	}
#ifdef _WINDOWS
	std::string cleanup = "del /Q *.xml *.fld *.chk *.his";
#else
    std::string cleanup = "rm -f *.xml *.fld *.chk *.his";
#endif
    system(cleanup.c_str());

}




