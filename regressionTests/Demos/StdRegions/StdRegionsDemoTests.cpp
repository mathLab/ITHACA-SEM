///////////////////////////////////////////////////////////////////////////////
//
// File StdRegionsDemoTests.cpp
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
// Description: Run a series of tests on StdRegions Demos using
// StdRegionsDemo class
//
///////////////////////////////////////////////////////////////////////////////
#include "../../Auxiliary/RegressBase.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

void RunL2RegressionTest(std::string demo, std::string input, std::string info);
void MakeOkFile(std::string demo, std::string input, std::string info);

#ifdef MAKE_OK_FILE
#define Execute($1,$2,$3)  MakeOkFile($1,$2,$3)
#else
#define Execute($1,$2,$3)  RunL2RegressionTest($1,$2,$3)
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
	
    // 1D Projection Tests
    Execute("StdProject1D", "1 6 7","Seg Ortho. Basis, P=6, Q=7");

    Execute("StdProject1D", "4 6 7","Seg Mod. Basis, P=6, Q=7");
	
	Execute("StdProject1D", "12 2 2","Single Mode Fourier Basis, P=2, Q=8");


    // 1D Differentiation Projection Tests
    Execute("StdProject_Diff1D", "1 6 7","Seg Ortho. Basis, P=6, Q=7");

    Execute("StdProject_Diff1D", "4 6 7","Seg Mod. Basis, P=6, Q=7");

    // 2D Projection Tests
    Execute("StdProject2D", "2 1 2 6 6 7 7", "Tri Ortho. Basis, P=6, Q=7");

    Execute("StdProject2D", "2 4 5 6 6 7 7", "Tri Mod. Basis P=6, Q=7");

    Execute("StdProject2D", "2 12 11 6 6 7 7", "Tri Nodal Basis P=6, Q=7");

    Execute("StdProject2D", "3 1 1 6 6 7 7", "Quad Ortho Basis, P=4, Q=5");

    Execute("StdProject2D", "3 4 4 6 6 7 7", "Quad Mod. Basis P=6, Q=7");

    Execute("StdProject2D", "3 8 8 6 6 7 7", "Quad Lagrange Basis P=6, Q=7");

    Execute("StdProject2D", "3 7 7 6 6 8 8 ", "Quad Fourier Basis P=6, Q=8");
	
	Execute("StdProject2D", "3 11 11 2 2 6 6", "Quad Fourier Single mode Basis P=2, Q=6");


    // 2D Differntiation and Projection Tests
    Execute("StdProject_Diff2D", "2 1 2 6 6 7 7 ", "Tri Ortho Basis, P=6, Q=7");

    Execute("StdProject_Diff2D", "2 4 5 6 6 7 7  ", "Tri Mod. Basis P=6, Q=7");

    Execute("StdProject_Diff2D", "2 11 11 6 6 7 7  ", "Tri Nodal Basis P=6, Q=7");

    Execute("StdProject_Diff2D", "3 1 1 6 6 7 7 ", "Quad Ortho Basis, P=6, Q=7");

    Execute("StdProject_Diff2D", "3 4 4 6 6 7 7  ", "Quad Modified Basis P=6, Q=7");

    Execute("StdProject_Diff2D", "3 8 8 6 6 7 7  ", "Quad Lagrange Basis P=6, Q=7");

    Execute("StdProject_Diff2D", "3 7 7 6 6 8 8 ", "Quad Fourier Basis P=6, Q=8");


    // 3D Projection Tests
    Execute("StdProject3D", "4 1 2 3 6 6 6 7 7 7", "Tet Ortho Basis, P=6, Q=7");

    Execute("StdProject3D", "4 4 5 6 6 6 6 7 7 7", "Tet Modified Basis P=6, Q=7");

    Execute("StdProject3D", "6 1 1 2 6 6 6 7 7 7", "Prism Ortho Basis P=6, Q=7");

    Execute("StdProject3D", "6 4 4 5 6 6 6 7 7 7", "Prism Modified Basis P=6 Q=7");

    Execute("StdProject3D", "7 1 1 1 6 6 6 7 7 7", "Hex Ortho Basis, P=4, Q=7");

    Execute("StdProject3D", "7 4 4 4 6 6 6 7 7 7", "Hex Modified Basis P=6, Q=7");

    // Fourier not yet working
//    Execute("StdProject3D", "7 7 7 7 6 6 6 8 8 8 ", "Hexahedral Fourier Basis modes=6, quad=8");

    Execute("StdProject3D", "7 8 8 8 6 6 6 7 7 7", "Hex Lagrange Basis P=6, Q=7");

    Execute("StdProject3D", "7 9 9 9 6 6 6 7 7 7", "Hex Legendre Basis P=6, Q=7");

    Execute("StdProject3D", "7 10 10 10 6 6 6 7 7 7", "Hex Chebyshev Basis P=6, Q=7");


    // 3D Differentiation and Projection Tests
    Execute("StdProject_Diff3D", "4 1 2 3 6 6 6 7 7 7", "Tet Ortho Basis, P=6, Q=7");

    Execute("StdProject_Diff3D", "4 4 5 6 6 6 6 7 7 7", "Tet Modified Basis P=6, Q=7");

    Execute("StdProject_Diff3D", "6 1 1 2 6 6 6 7 7 7", "Prism Ortho Basis P=6, Q=7");

    Execute("StdProject_Diff3D", "6 4 4 5 6 6 6 7 7 7", "Prism Modified Basis P=6 Q=7");

    Execute("StdProject_Diff3D", "7 1 1 1 6 6 6 7 7 7", "Hex Ortho Basis, P=4, Q=7");

    Execute("StdProject_Diff3D", "7 4 4 4 6 6 6 7 7 7", "Hex Modified Basis P=6, Q=7");

    Execute("StdProject_Diff3D", "7 8 8 8 6 6 6 7 7 7", "Hex Lagrange Basis P=6, Q=7");

    Execute("StdProject_Diff3D", "7 9 9 9 6 6 6 7 7 7", "Hex Legendre Basis P=6, Q=7");

    Execute("StdProject_Diff3D", "7 10 10 10 6 6 6 7 7 7", "Hex Cheby. Basis P=6, Q=7");

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
    RegressBase Test(NEKTAR_BIN_DIR,Demo,input,"Demos/StdRegions/OkFiles/");
    int fail;

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
};

void MakeOkFile(std::string Demo, std::string input, std::string info)
{
    tests_total++;
    RegressBase Test(NEKTAR_BIN_DIR,Demo,input,"Demos/StdRegions/OkFiles/");
    int fail;

    if(fail = Test.MakeOkFile())
    {
        std::cout << "Failed to make OK file\n";
        // Explain cause of error if available
        std::cout << "===========================================================\n";
        Test.PrintTestError(fail);
        std::cout << "===========================================================\n";
        tests_failed++;
    }
}




