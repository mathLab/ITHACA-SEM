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

int main(int argc, char* argv[]) 
{ 
    //Test Channel Flow
	Execute("IncNavierStokesSolver","Test_ChanFlow_m3.xml","Testing Channel Flow modes=3");
	Execute("IncNavierStokesSolver","Test_ChanFlow_m8.xml","Testing Channel Flow modes=8");
	
	//Test Kovasznay Flow
	Execute("IncNavierStokesSolver","Test_KovaFlow_m3.xml","Testing Kovasznay Flow modes=3");
	Execute("IncNavierStokesSolver","Test_KovaFlow_m8.xml","Testing Kovasznay Flow modes=8");
	
	//Test Decaying Vortex
	Execute("IncNavierStokesSolver","Test_TaylorVor_m3.xml","Testing Taylor Vortex modes=3");
	Execute("IncNavierStokesSolver","Test_TaylorVor_m8.xml","Testing Taylor Vortex modes=8");
	
    return 0;
}

void RunL2RegressionTest(std::string Demo, std::string input, std::string info)
{
	std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "/dist/bin/";
    RegressBase Test(NektarSolverDir.c_str(),Demo,input,"Solvers/IncNavierStokesSolver/OkFiles/");
    int fail;

    // Copy input file to current location
	input.erase(input.end()-3,input.end());
	boost::filesystem::path filePath(std::string(REG_PATH) + "Solvers/IncNavierStokesSolver/InputFiles/" + input);
    std::string syscommand = std::string(COPY_COMMAND) + filePath.file_string() + "xml .";
	std::string syscommand2 = std::string(COPY_COMMAND) + filePath.file_string() + "rst .";
	
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
	
	input = input+"xml";

    std::cout << Demo << ":  info = \"" << info <<"\""<<std::endl;
    if(fail = Test.TestL2()) // test failed
    {
        std::cout <<" status: FAILED" << std::endl;
        std::cout << "===========================================================\n";   
        // Explain cause of error if available
        Test.PrintTestError(fail);
        std::cout << "===========================================================\n";
    }
    else
    {
        std:: cout <<" status: PASSED" << std::endl;
    }            
    
#ifdef _WINDOWS
	std::string cleanup = "del /Q *.xml *.fld *.chk *.rst";
#else
    std::string cleanup = "rm -f *.xml *.fld *.chk *.rst";
#endif

    system(cleanup.c_str());
};

void MakeOkFile(std::string Demo, std::string input,				std::string info)
{
	std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "/dist/bin/";
    RegressBase Test(NektarSolverDir.c_str(),Demo,input,"Solvers/IncNavierStokesSolver/OkFiles/");
    int fail;


    // Copy input file to current location
	input.erase(input.end()-3,input.end());
	boost::filesystem::path filePath(std::string(REG_PATH) + "Solvers/IncNavierStokesSolver/InputFiles/" + input);
	std::string syscommand1 = std::string(COPY_COMMAND) + filePath.file_string() + "xml .";
	std::string syscommand2 = std::string(COPY_COMMAND) + filePath.file_string() + "rst .";
	
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
	
	input = input+"xml";

    if(fail = Test.MakeOkFile())
    {
        std::cout << "Failed to make OK file\n";
        // Explain cause of error if available
        std::cout << "===========================================================\n";   
        Test.PrintTestError(fail);
        std::cout << "===========================================================\n";    
	}
#ifdef _WINDOWS
	std::string cleanup = "del /Q *.xml *.fld *.chk *.rst";
#else
    std::string cleanup = "rm -f *.xml *.fld *.chk *.rst";
#endif
	
    system(cleanup.c_str());
}



	
