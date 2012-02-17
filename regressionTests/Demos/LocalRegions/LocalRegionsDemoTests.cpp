///////////////////////////////////////////////////////////////////////////////
//
// File LocalRegionsDemoTests.cpp
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
// Description: Run a series of tests on LocalRegions Demos using
// LocalRegionsDemo class
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
    Execute("LocProject1D", "1 6 7 0.0 1.0","Seg [0,1] Ortho Basis, P=6, Q=7");

    Execute("LocProject1D", "4 6 7 0.0 1.0","Seg [0,1] Mod. Basis, P=6, Q=7");

    // 1D Differentiation Projection Tests
    Execute("LocProject_Diff1D", "1 6 7 0.0 1.0","Seg [0,1] Ortho Basis, P=6, Q=7");

    Execute("LocProject_Diff1D", "4 6 7 0.0 1.0","Seg [0.1] Mod. Basis, P=6, Q=7");

    // 2D Projection Tests
    Execute("LocProject2D", "2 1 2 6 6 7 7 0.0 0.0 1.0 1.0 0.5 1.0", "Tri Ortho Basis, P=6, Q=7");

    Execute("LocProject2D", "2 4 5 6 6 7 7 0.0 0.0 1.0 1.0 0.5 1.0", "Tri Modified Basis P=6, Q=7");

    Execute("LocProject2D", "2 11 11 6 6 7 7  0.0 0.0 1.0 1.0 0.5 1.0", "Tri Nodal Basis P=6, Q=7");

    Execute("LocProject2D", "3 1 1 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Reg. Quad Ortho Basis, P=4, Q=5");

    Execute("LocProject2D", "3 4 4 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Reg. Quad Mod. Basis P=6, Q=7");

    Execute("LocProject2D", "3 8 8 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Reg. Quad Lagrange Basis P=6, Q=7");

    Execute("LocProject2D", "3 1 1 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Lin. Deformed Quad Ortho Basis, P=4, Q=5");

    Execute("LocProject2D", "3 4 4 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Lin. Deformed Quad Modified Basis P=6, Q=7");

    Execute("LocProject2D", "3 8 8 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Lin. Deformed Quad Lagrange Basis P=6, Q=7");

    // 2D Differntiation and Projection Tests
    Execute("LocProject_Diff2D", "2 1 2 6 6 7 7 0.0 0.0 1.0 1.0 0.5 1.0", "Tri Ortho Basis, P=6, Q=7");

    Execute("LocProject_Diff2D", "2 4 5 6 6 7 7 0.0 0.0 1.0 1.0 0.5 1.0", "Tri Mod. Basis P=6, Q=7");

    Execute("LocProject_Diff2D", "2 11 11 6 6 7 7  0.0 0.0 1.0 1.0 0.5 1.0", "Tri Nodal Basis P=6, Q=7");

    Execute("LocProject_Diff2D", "3 1 1 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Reg. Quad Ortho Basis, P=4, Q=5");

    Execute("LocProject_Diff2D", "3 4 4 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Reg. Quad Mod. Basis P=6, Q=7");

    Execute("LocProject_Diff2D", "3 8 8 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Reg. Quad Lagrange Basis P=6, Q=7");

    Execute("LocProject_Diff2D", "3 1 1 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Lin. Deformed Quad Ortho Basis, P=4, Q=5");

    Execute("LocProject_Diff2D", "3 4 4 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Lin. Deformed Quad Mod. Basis P=6, Q=7");

    Execute("LocProject_Diff2D", "3 8 8 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Lin. Deformed Quad Lagrange Basis P=6, Q=7");


    // 3D Projection Tests
    Execute("LocProject3D", "4 1 2 3 6 6 6 7 7 7 0 0 0  1 0 0  0 1 0  0 0 1", "Reg. Tet Ortho Basis, P=6, Q=7");

    Execute("LocProject3D", "4 4 5 6 6 6 6 7 7 7 0 0 0  1 0 0  0 1 0  0 0 1", "Reg. Tet Mod. Basis, P=6, Q=7");
    
    Execute("LocProject3D", "6 1 1 2 6 6 6 7 7 6 0 0 0  1 0 0  1 1 0  0 1 0  0.5 0 1  0.5 1 1", "Reg. Prism Ortho Basis, P=6, Q=7");

    Execute("LocProject3D", "6 4 4 5 6 6 6 7 7 6 0 0 0  1 0 0  1 1 0  0 1 0  0.5 0 1  0.5 1 1", "Reg. Prism Mod. Basis, P=6, Q=7");

    Execute("LocProject3D", "7 1 1 1 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Reg. Hex Ortho Basis, P=6, Q=7");

    Execute("LocProject3D", "7 4 4 4 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Reg. Hex Mod. Basis, P=6, Q=7");

    Execute("LocProject3D", "7 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Reg. Hex Lagrange Basis, P=6, Q=7");

    Execute("LocProject3D", "7 1 1 1 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Lin. Deformed Hex Ortho Basis, P=6, Q=7");

    Execute("LocProject3D", "7 4 4 4 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Lin. Deformed Hex Mod. Basis, P=6, Q=7");

    Execute("LocProject3D", "7 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Lin. Deformed Hex Lagrange Basis, P=6, Q=7");

    // 3D Differentiation and Projection Tests
    Execute("LocProject_Diff3D", "4 1 2 3 6 6 6 7 7 7 0 0 0  1 0 0  0 1 0  0 0 1", "Reg. Tet Ortho Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "4 4 5 6 6 6 6 7 7 7 0 0 0  1 0 0  0 1 0  0 0 1", "Reg. Tet Mod. Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "6 1 1 2 6 6 6 7 7 6 0 0 0  1 0 0  1 1 0  0 1 0  0.5 0 1  0.5 1 1", "Reg. Prism Ortho Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "6 4 4 5 6 6 6 7 7 6 0 0 0  1 0 0  1 1 0  0 1 0  0.5 0 1  0.5 1 1", "Reg. Prism Mod. Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "7 1 1 1 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Reg. Hex Ortho Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "7 4 4 4 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Reg. Hex Mod. Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "7 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Reg. Hex Lagrange Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "7 1 1 1 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Lin. Deformed Hex Ortho Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "7 4 4 4 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Lin. Deformed Hex Mod. Basis, P=6, Q=7");

    Execute("LocProject_Diff3D", "7 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Lin. Deformed Hex Lagrange Basis, P=6, Q=7");

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
    RegressBase Test(NEKTAR_BIN_DIR,Demo,input,"Demos/LocalRegions/OkFiles/");
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
    RegressBase Test(NEKTAR_BIN_DIR,Demo,input,"Demos/LocalRegions/OkFiles/");
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




