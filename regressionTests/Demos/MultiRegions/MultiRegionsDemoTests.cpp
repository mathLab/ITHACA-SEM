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
// Description: Run a series of tests on MultiRegions Demos
//
///////////////////////////////////////////////////////////////////////////////
#include "../../Auxiliary/RegressBase.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <boost/filesystem.hpp>

void RunL2RegressionTest(std::string demo, std::string input, std::string info, unsigned int np = 1);
void MakeOkFile(std::string demo, std::string input, std::string info);

#ifdef MAKE_OK_FILE
#define Execute($1,$2,$3)  MakeOkFile($1,$2,$3)
#define ExecuteParallel($1,$2,$3,$4) MakeOkFile($1,$2,$3)
#else
#define Execute($1,$2,$3)  RunL2RegressionTest($1,$2,$3)
#define ExecuteParallel($1,$2,$3,$4)  RunL2RegressionTest($1,$2,$3,$4)
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
    // 1D Demos
    Execute("Helmholtz1D", "helmholtz1D_8modes.xml","CG Helmholtz1D  P=8");

    Execute("Helmholtz1D", "helmholtz1D_8modes_RBC.xml","CG Helmholtz1D  P=8 with Robin BC");

    Execute("HDGHelmholtz1D", "helmholtz1D_8modes.xml","HDG Helmholtz1D  P=8");

    Execute("HDGHelmholtz1D", "helmholtz1D_8modes_RBC.xml","HDG Helmholtz1D  P=8 with Robin BC");

    // 2D Demos
    Execute("Helmholtz2D", "helmholtz2D_7modes.xml","CG Helmholtz2D  P=7");

    Execute("Helmholtz2D", "helmholtz2D_7nodes.xml","CG Helmholtz2D  P=7");

    Execute("Helmholtz2D", "helmholtz2D_7modes_AllBCs_mlsc.xml","CG Helmholtz2D  P=7 All BCs, Direct ML Static Condensation");

    Execute("Helmholtz2D", "helmholtz2D_7modes_AllBCs_sc.xml","CG Helmholtz2D  P=7 All BCs, Direct Static Condensation");

    Execute("Helmholtz2D", "helmholtz2D_7modes_AllBCs.xml","CG Helmholtz2D  P=7 All BCs, Direct Full Matrix");
    
    //Execute("Helmholtz2D", "helmholtz2D_7modes_AllBCs_iter.xml","CG Helmholtz2D  P=7 All BCs, Iterative Full Matrix");

    Execute("Helmholtz2D", "helmholtz2D_7modes_AllBCs_iter_ml.xml","CG Helmholtz2D  P=7 All BCs, Iterative ML Static Condensation");

    Execute("Helmholtz2D", "helmholtz2D_9modes_varcoeff.xml","CG Helmholtz2D P=9 Var Coeffs");

    Execute("Helmholtz2D", "helmholtz2D_curved_quad.xml","CG Helmholtz2D P=7 curved quads");

    Execute("Helmholtz2D", "helmholtz2D_curved_tri.xml","CG Helmholtz2D P=7 curved tris");

    Execute("HDGHelmholtz2D", "helmholtz2D_7modes.xml","HDG Helmholtz2D  P=7");

    Execute("HDGHelmholtz2D", "helmholtz2D_7modes_AllBCs.xml","HDG Helmholtz2D  P=7 All BCs, ML Static Condensation");

    Execute("SteadyAdvectionDiffusionReaction2D", "linearadvdiffreact2D_7modes.xml","Steady ADR 2D  P=7");

    // 3D Demos
    Execute("Helmholtz3D", "helmholtz3D_hex.xml","CG Helmholtz3D Hex");

    Execute("Helmholtz3D", "helmholtz3D_tet.xml","CG Helmholtz3D Tet");

    Execute("Helmholtz3D", "helmholtz3D_prism.xml","CG Helmholtz3D Prism");

    Execute("Helmholtz3D", "helmholtz3D_prism_deformed.xml","CG Helmholtz3D Prism (deformed)");

    // 3D Demos Homogeneous1D
    Execute("Helmholtz3DHomo1D", "helmholtz3D_homo1D_7modes_8nz.xml","CG Helmholtz3D Homogeneous 1D");

    Execute("HDGHelmholtz3DHomo1D", "helmholtz3D_homo1D_7modes_8nz.xml","HDG Helmholtz3D Homogeneous 1D");
	
	Execute("Deriv3DHomo1D", "derivatives3Dhomo1D.xml","Testing 3D homogeneous 1D derivatives");
	
	Execute("Deriv3DHomo2D", "derivatives3Dhomo2D.xml","Testing 3D homogeneous 2D derivatives");

#ifdef NEKTAR_USE_MPI
    ExecuteParallel("Helmholtz1D", "helmholtz1D_8modes.xml","Par(2) CG Helmholtz1D  P=8", 2);
//    ExecuteParallel("Helmholtz2D", "helmholtz2D_7modes.xml","Par(2) CG Helmholtz2D  P=7", 2);
//    ExecuteParallel("Helmholtz2D", "helmholtz2D_7modes.xml","Par(5) CG Helmholtz2D  P=7", 5);
    ExecuteParallel("Helmholtz3D", "helmholtz3D_hex.xml",   "Par(2) CG Helmholtz3D Hex",  2);
#endif

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

void RunL2RegressionTest(std::string Demo, std::string input, std::string info, unsigned int np)
{
    tests_total++;
    if (!quiet)
    {
        std::cout << "TESTING: " << Demo << std::flush;
    }
#ifdef USE_PATH
    RegressBase Test(NEKTAR_BIN_DIR,Demo,input,"Demos/MultiRegions/OkFiles/",np);
#else
    RegressBase Test("",Demo,input,"Demos/MultiRegions/OkFiles/",np);
#endif
    int fail;

    // Copy input file to current location
	boost::filesystem::path filePath(std::string(REG_PATH) + "Demos/MultiRegions/InputFiles/" + input);
    std::string syscommand = std::string(COPY_COMMAND) + PortablePath(filePath) + " .";
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
	std::string cleanup = "del /Q *.xml *.fld";
#else
    std::string cleanup = "rm -f *.xml *.fld";
#endif

    system(cleanup.c_str());
};

void MakeOkFile(std::string Demo, std::string input, std::string info)
{
    tests_total++;
    RegressBase Test(NEKTAR_BIN_DIR,Demo,input,"Demos/MultiRegions/OkFiles/");
    int fail;


    // Copy input file to current location
	boost::filesystem::path filePath(std::string(REG_PATH) + "Demos/MultiRegions/InputFiles/" + input);
    std::string syscommand = std::string(COPY_COMMAND) + PortablePath(filePath) + " .";
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
}




