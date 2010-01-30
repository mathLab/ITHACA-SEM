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

void RunL2RegressionTest(std::string demo, std::string input, std::string info);
void MakeOkFile(std::string demo, std::string input, std::string info);

#ifdef MAKE_OK_FILE
#define Execute($1,$2,$3)  MakeOkFile($1,$2,$3)
#else
#define Execute($1,$2,$3)  RunL2RegressionTest($1,$2,$3)
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
	Execute("IncNavierStokesSolver","Test_DecVortex_m3.xml","Testing Decaying Vortex modes=3");
	Execute("IncNavierStokesSolver","Test_DecVortex_m8.xml","Testing Decaying Vortex modes=8");
	
    return 0;
}

void RunL2RegressionTest(std::string Demo, std::string input, std::string info)
{
    RegressBase Test("../solvers/builds/IncNavierStokesSolver/",Demo,input,"Solvers/IncNavierStokesSolver/OkFiles/");
    int fail;

    // Copy input file to current location
    std::string syscommand = "cp ../../../Solvers/IncNavierStokesSolver/InputFiles/"+input +" .";
    int status = system(syscommand.c_str());
    if(status)
    {
        std::cerr << "Unable to copy file:" << input << " to current location" << std::endl;
        exit(2);
    }

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
    
    std::string cleanup = "rm -f *.xml *.fld";
    system(cleanup.c_str());
};

void MakeOkFile(std::string Demo, std::string input, std::string info)
{

    RegressBase Test("../solvers/builds/IncNavierStokesSolver/",Demo,input,"Solvers/IncNavierStokesSolver/OkFiles/");
    int fail;


    // Copy input file to current location
    std::string syscommand = "cp ../../../Solvers/IncNavierStokesSolver/InputFiles/"+input +" .";
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
}
}



	
