///////////////////////////////////////////////////////////////////////////////
//
// File MultiRegionsDemoTests.cpp
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
// Description: Run a series of tests on the Incompressible Navier Stokes Solver
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
    
    //Test Channel Flow
    Execute("IncNavierStokesSolver","Test_ChanFlow_m3.xml","Channel Flow P=3");
    Execute("IncNavierStokesSolver","Test_ChanFlow_m8.xml","Channel Flow P=8");
    Execute("IncNavierStokesSolver","Test_ChanFlow_m8_BodyForce.xml","Channel Flow P=8 BodyForce");    
    Execute("IncNavierStokesSolver","Test_ChanFlow_m8_singular.xml","Channel Flow P=8 Singularity Check");
    Execute("IncNavierStokesSolver","Test_ChanFlow2D_bcsfromfiles.xml","Channel Flow P=5 Boundary Conditions from files");  
    
    //Test Kovasznay Flow
    Execute("IncNavierStokesSolver","Test_KovaFlow_m3.xml","Kovasznay Flow P=3");
    Execute("IncNavierStokesSolver","Test_KovaFlow_m8.xml","Kovasznay Flow P=8");
    
    //Test Decaying Vortex
    Execute("IncNavierStokesSolver","Test_TaylorVor_dt1.xml","Convergence: Taylor Vortex IMEXOrder2 dt=0.01");
    Execute("IncNavierStokesSolver","Test_TaylorVor_dt2.xml","Convergence: Taylor Vortex IMEXOrder2 dt=0.001");


    //Test Coupled LinearNS Kovasznay Flow
    Execute("IncNavierStokesSolver","Test_KovaFlow_m8.xml","Steady Oseen Kovasznay flow P=6");
    Execute("IncNavierStokesSolver","Test_Kovas_Quad6_Tri4_mixedbcs.xml","Steady Oseen Kovasznay flow, mixed elements and bcs P=7");


    //Test Coupled LinearNS unsteady Channel Flow
    Execute("IncNavierStokesSolver","Test_ChanFlow_LinNS_m8.xml","Unsteady channel flow with coupled solve , P=8");

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
	std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "/dist/bin/";
    RegressBase Test(NektarSolverDir.c_str(),Demo,input,"Solvers/IncNavierStokesSolver/OkFiles/");
    int fail;
    
	std::string input1=input;
	input1.erase(input1.end()-4,input1.end());
        
    // Copy input file to current location
	input.erase(input.end()-3,input.end());
	boost::filesystem::path filePath(std::string(REG_PATH) + "Solvers/IncNavierStokesSolver/InputFiles/" + input);
	boost::filesystem::path filePath1(std::string(REG_PATH) + "Solvers/IncNavierStokesSolver/InputFiles/" + input1);	
	std::string syscommand = std::string(COPY_COMMAND) + filePath.file_string() + "xml .";
	std::string syscommand2 = std::string(COPY_COMMAND) + filePath.file_string() + "rst .";
	std::string syscommand3 = std::string(COPY_COMMAND) + filePath1.file_string() + "_u_1.bc .";	
	std::string syscommand4 = std::string(COPY_COMMAND) + filePath1.file_string() + "_u_3.bc .";
	std::string syscommand5 = std::string(COPY_COMMAND) + filePath1.file_string() + "_v_3.bc .";	

	int status = system(syscommand.c_str());
    if(status)
    {
        std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
        exit(2);
    }

	//Restart files just needed for the Kovasznay Flow so far
	if(input=="Test_KovaFlow_m3." || input=="Test_KovaFlow_m8.")
	{
	   int status2 = system(syscommand2.c_str());
		if(status2)
		{
			std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
			exit(2);
		}

	}	
	//files for the Test_ChanFlow2D_bcsfromfiles.xml
	if(input=="Test_ChanFlow2D_bcsfromfiles.")
	{
	   int status3 = system(syscommand3.c_str());
	   int status4 = system(syscommand4.c_str());
	   int status5 = system(syscommand5.c_str());	   
		if((status3)||(status4)||(status5))
		{
			std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
			exit(2);
		}
	}	

	input = input+"xml";

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
	std::string cleanup = "del /Q *.xml *.fld *.chk *.rst *.bc";
#else
    std::string cleanup = "rm -f *.xml *.fld *.chk *.rst *.bc";
#endif

    system(cleanup.c_str());
};

void MakeOkFile(std::string Demo, std::string input,				std::string info)
{
    tests_total++;
	std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "/dist/bin/";
    RegressBase Test(NektarSolverDir.c_str(),Demo,input,"Solvers/IncNavierStokesSolver/OkFiles/");
    int fail;
    
	std::string input1=input;
	input1.erase(input1.end()-4,input1.end());    


    // Copy input file to current location
	input.erase(input.end()-3,input.end());
	boost::filesystem::path filePath(std::string(REG_PATH) + "Solvers/IncNavierStokesSolver/InputFiles/" + input);	
	boost::filesystem::path filePath1(std::string(REG_PATH) + "Solvers/IncNavierStokesSolver/InputFiles/" + input1);	
	std::string syscommand1 = std::string(COPY_COMMAND) + filePath.file_string() + "xml .";
	std::string syscommand2 = std::string(COPY_COMMAND) + filePath.file_string() + "rst .";
	std::string syscommand3 = std::string(COPY_COMMAND) + filePath1.file_string() + "_u_1.bc .";	
	std::string syscommand4 = std::string(COPY_COMMAND) + filePath1.file_string() + "_u_3.bc .";
	std::string syscommand5 = std::string(COPY_COMMAND) + filePath1.file_string() + "_v_3.bc .";	

    int status1 = system(syscommand1.c_str());
    if(status1)
    {
        std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
        exit(2);
    }

	//Restart files just needed for the Kovasznay Flow so far
	if((input == "Test_KovaFlow_m3.") || (input == "Test_KovaFlow_m8."))
	{
		int status2 = system(syscommand2.c_str());
		if(status2)
		{
			std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
			exit(2);
		}

	}
	//files for the Test_ChanFlow2D_bcsfromfiles.xml
	if(input=="Test_ChanFlow2D_bcsfromfiles.")
	{
	   int status3 = system(syscommand3.c_str());
	   int status4 = system(syscommand4.c_str());
	   int status5 = system(syscommand5.c_str());	   
		if((status3)||(status4)||(status5))
		{
			std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
			exit(2);
		}
	}

	input = input+"xml";

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
	std::string cleanup = "del /Q *.xml *.fld *.chk *.rst *.bc";
#else
    std::string cleanup = "rm -f *.xml *.fld *.chk *.rst *.bc";
#endif

    system(cleanup.c_str());
}




