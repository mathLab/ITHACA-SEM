///////////////////////////////////////////////////////////////////////////////
//
// File ADRSolverTests.cpp
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
// Description: Run a series of tests on the ADR Solver 
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
    // //Test Steady Diffusion Advection
    Execute("ADRSolver","Test_Helmholtz1D_8modes.xml","Testing 1D Helmholtz/Steady Diffusion Reaction Modes=8");

    Execute("ADRSolver","Test_Helmholtz1D_8modes_DG.xml","Testing 1D Helmholtz/Steady Diffusion Reaction (DG) Modes=8");

    Execute("ADRSolver","Test_Helmholtz1D_8nodes.xml","Testing 1D Helmholtz/Steady Diffusion Reaction Modes=8");

    Execute("ADRSolver","Test_Helmholtz2D_modal.xml","Testing 2D Helmholtz/Steady Diffusion Reaction Modes=7");

    Execute("ADRSolver","Test_Helmholtz2D_modal_DG.xml","Testing 2D Helmholtz/Steady Diffusion Reaction (DG) Modes=7");

    Execute("ADRSolver","Test_Helmholtz2D_nodal.xml","Testing 2D Helmholtz/Steady Diffusion Reaction Modes=7");

    // Test Steady Advection Diffusion Reaction 
    Execute("ADRSolver","Test_SteadyAdvDiffReact2D_modal.xml","Testing 2D Steady Advection Diffusion Reaction Modes=9");

    // //Test Advection
    Execute("ADRSolver","Test_Advection1D_m8_Order1.xml","Testing 1D unsteady advection, order 1, modes=8");
    Execute("ADRSolver","Test_Advection1D_m12_Order2.xml","Testing 1D unsteady advection, order 2, modes=12");
    Execute("ADRSolver","Test_Advection1D_m14_Order4.xml","Testing 1D unsteady advection, order 4, modes=14");

    Execute("ADRSolver","Test_Advection_m12_Order1.xml","Testing 2D unsteady advection, order 1, modes=12");
    Execute("ADRSolver","Test_Advection_m12_Order2.xml","Testing 2D unsteady advection, order 2, modes=12");
    Execute("ADRSolver","Test_Advection_m14_Order4.xml","Testing 2D unsteady advection, order 4, modes=14");
    Execute("ADRSolver","Test_Advection_m12_DG_Order1.xml","Testing 2D unsteady DG advection, order 1, modes=12");
    Execute("ADRSolver","Test_Advection_m12_DG_Order2.xml","Testing 2D unsteady DG advection, order 2, modes=12");
    Execute("ADRSolver","Test_Advection_m14_DG_Order4.xml","Testing 2D unsteady DG advection, order 4, modes=14");
    
    return 0;
}

void RunL2RegressionTest(std::string Demo, std::string input, std::string info)
{
    std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "ADRSolver/";
    RegressBase Test(NektarSolverDir.c_str(),Demo,input,"Solvers/ADRSolver/OkFiles/");
    int fail;

    // Copy input file to current location
    std::string syscommand = std::string("cp ") + REG_PATH + "Solvers/ADRSolver/InputFiles/"+input +" .";
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
    std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "ADRSolver/";
    RegressBase Test(NektarSolverDir,Demo,input,"Solvers/ADRSolver/OkFiles/");
    int fail;

    // Copy input file to current location
    std::string syscommand = std::string("cp ") + REG_PATH + "Solvers/ADRSolver/InputFiles/"+input +" .";
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



	
