////////////////////////////////////////////////////////////////////OS///////////
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

    // Test Steady Diffusion Advection
    Execute("ADRSolver","Test_Helmholtz1D_8modes.xml","1D Helmholtz/Steady Diffusion Reaction P=8");

    Execute("ADRSolver","Test_Helmholtz1D_8modes_DG.xml","1D Helmholtz/Steady Diffusion Reaction (DG) P=8");

    Execute("ADRSolver","Test_Helmholtz1D_8nodes.xml","1D Helmholtz/Steady Diffusion Reaction P=8");

    Execute("ADRSolver","Test_Helmholtz2D_modal.xml","2D Helmholtz/Steady Diffusion Reaction P=7");

    Execute("ADRSolver","Test_Helmholtz2D_modal_DG.xml","2D Helmholtz/Steady Diffusion Reaction (DG) P=7");

    Execute("ADRSolver","Test_Helmholtz2D_nodal.xml","2D Helmholtz/Steady Diffusion Reaction P=7");

    Execute("ADRSolver","Test_Helmholtz3D_modal.xml","3D Helmholtz/Steady Diffusion Reaction P=7");

    Execute("ADRSolver","Test_Helmholtz3D_nodal.xml","3D Helmholtz/Steady Diffusion Reaction P=7");
	
	Execute("ADRSolver","Test_Helmholtz_3DHomo1D_MVM.xml","3D-Homogeneous-1D Helmholtz/Steady Diffusion (MVM)");
	
	Execute("ADRSolver","Test_Helmholtz_3DHomo2D_MVM.xml","3D-Homogeneous-2D Helmholtz/Steady Diffusion (MVM)");
	
#ifdef NEKTAR_USING_FFTW
	Execute("ADRSolver","Test_Helmholtz_3DHomo1D_FFT.xml","3D-Homogeneous-1D Helmholtz/Steady Diffusion (FFT)");
	
	Execute("ADRSolver","Test_Helmholtz_3DHomo2D_FFT.xml","3D-Homogeneous-2D Helmholtz/Steady Diffusion (FFT)");
#endif
	
    // Test Steady Advection Diffusion Reaction
    Execute("ADRSolver","Test_SteadyAdvDiffReact2D_modal.xml","2D Steady Advection Diffusion Reaction P=9");

    Execute("ADRSolver","Test_Advection1D_m7_DG.xml","1D unsteady DG advection, P=7");

    Execute("ADRSolver","Test_Advection1D_m8_Order1.xml","1D unsteady advection, order 1, P=8");

    Execute("ADRSolver","Test_Advection1D_m12_Order2.xml","1D unsteady advection, order 2, P=12");

    Execute("ADRSolver","Test_Advection1D_m14_Order4.xml","1D unsteady advection, order 4, P=14");

    Execute("ADRSolver","Test_Advection_m12_Order1.xml","2D unsteady advection, order 1, P=12");

    Execute("ADRSolver","Test_Advection_m12_Order2.xml","2D unsteady advection, order 2, P=12");

    Execute("ADRSolver","Test_Advection_m14_Order4.xml","2D unsteady advection, order 4, P=14");

    Execute("ADRSolver","Test_Advection_m12_DG_Order1.xml","2D unsteady DG advection, order 1, P=12");

    Execute("ADRSolver","Test_Advection_m12_DG_Order2.xml","2D unsteady DG advection, order 2, P=12");

    Execute("ADRSolver","Test_Advection_m14_DG_Order4.xml","2D unsteady DG advection, order 4, P=14");

    Execute("ADRSolver","Test_Advection_m12_DG_periodic.xml","2D unsteady DG advection, order 1, P=12, periodic bcs");

    Execute("ADRSolver","Test_Advection3D_m12_DG_hex_periodic.xml","3D unsteady DG advection, hexahedra, order 1, P=12, periodic bcs");

    Execute("ADRSolver","Test_Advection3D_m12_DG_hex.xml","3D unsteady DG advection, hexahedra, order 4, P=12");

    Execute("ADRSolver","Test_Advection3D_m12_DG_tet.xml","3D unsteady DG advection, tetrahedra, order 4, P=12");

    Execute("ADRSolver","Test_Advection3D_m12_DG_prism.xml","3D unsteady DG advection, prisms, order 4, P=14");

    Execute("ADRSolver","Test_ImDiffusion_m6.xml","2D unsteady DG implicit diffusion, order 3, P=6");

    Execute("ADRSolver","Test_ImDiffusion_m12.xml","2D unsteady DG implicit diffusion, order 3, P=12");

    Execute("ADRSolver","Test_ImDiffusion_VarCoeff.xml","2D unsteady CG implicit diffusion, variable coeffs.");

    Execute("ADRSolver","Test_ExDiffusion_m3.xml","2D unsteady DG explicit diffusion, order 4, P=3");

    Execute("ADRSolver","Test_ExDiffusion_m8.xml","2D unsteady DG explicit diffusion, order 4, P=8");

	Execute("ADRSolver","Test_UnsteadyAdvectionDiffusion_Order1_001.xml","2D unsteady advection-diffusion, IMEXOrder1, P=9, dt=0.001");

	Execute("ADRSolver","Test_UnsteadyAdvectionDiffusion_Order1_0001.xml","2D unsteady advection-diffusion, IMEXOrder1, P=9, dt=0.0001");

	Execute("ADRSolver","Test_UnsteadyAdvectionDiffusion_Order2_001.xml","2D unsteady advection-diffusion, IMEXOrder2, P=9, dt=0.001");

	Execute("ADRSolver","Test_UnsteadyAdvectionDiffusion_Order2_0001.xml","2D unsteady advection-diffusion, IMEXOrder2, P=9, dt=0.0001");
	
	Execute("ADRSolver","Test_UnsteadyAdvectionDiffusion_3DHomo1D_MVM.xml","3D-Homogeneous-1D unsteady advection-diffusion (MVM)");
	
	Execute("ADRSolver","Test_UnsteadyAdvectionDiffusion_3DHomo2D_MVM.xml","3D-Homogeneous-2D unsteady advection-diffusion (MVM)");
    
    // Test inviscid Burger equation in 1D for DG and FR
    Execute("ADRSolver","Test_InviscidBurger1D_WeakDG_GLL_LAGRANGE.xml","1D unsteady WeakDG inviscidBurger GLL_LAGRANGE, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_WeakDG_MODIFIED.xml","1D unsteady WeakDG inviscidBurger MODIFIED, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRDG_GLL_LAGRANGE.xml","1D unsteady FRDG inviscidBurger GLL_LAGRANGE, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRDG_MODIFIED.xml","1D unsteady FRDG inviscidBurger MODIFIED, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRSD_GLL_LAGRANGE.xml","1D unsteady FRSD inviscidBurger GLL_LAGRANGE, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRSD_MODIFIED.xml","1D unsteady FRSD inviscidBurger MODIFIED, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRHU_GLL_LAGRANGE.xml","1D unsteady FRHU inviscidBurger GLL_LAGRANGE, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRHU_MODIFIED.xml","1D unsteady FRHU inviscidBurger MODIFIED, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRDG_GLL_LAGRANGE_SEM.xml","1D unsteady FRDG inviscidBurger GLL_LAGRANGE_SEM, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRSD_GLL_LAGRANGE_SEM.xml","1D unsteady FRSD inviscidBurger GLL_LAGRANGE_SEM, P=10");
    
    Execute("ADRSolver","Test_InviscidBurger1D_FRHU_GLL_LAGRANGE_SEM.xml","1D unsteady FRHU inviscidBurger GLL_LAGRANGE_SEM, P=10");
    
    
    // Test linear advection equation in 1D for DG and FR
    Execute("ADRSolver","Test_Advection1D_WeakDG_GLL_LAGRANGE.xml","1D unsteady WeakDG advection GLL_LAGRANGE, P=3");
    
    Execute("ADRSolver","Test_Advection1D_WeakDG_MODIFIED.xml","1D unsteady WeakDG advection MODIFIED, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRDG_GLL_LAGRANGE.xml","1D unsteady FRDG advection GLL_LAGRANGE, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRDG_MODIFIED.xml","1D unsteady FRDG advection MODIFIED, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRSD_GLL_LAGRANGE.xml","1D unsteady FRSD advection GLL_LAGRANGE, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRSD_MODIFIED.xml","1D unsteady FRSD advection MODIFIED, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRHU_GLL_LAGRANGE.xml","1D unsteady FRHU advection GLL_LAGRANGE, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRHU_MODIFIED.xml","1D unsteady FRHU advection MODIFIED, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRDG_GLL_LAGRANGE_SEM.xml","1D unsteady FRDG advection GLL_LAGRANGE_SEM, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRSD_GLL_LAGRANGE_SEM.xml","1D unsteady FRSD advection GLL_LAGRANGE_SEM, P=3");
    
    Execute("ADRSolver","Test_Advection1D_FRHU_GLL_LAGRANGE_SEM.xml","1D unsteady FRHU advection GLL_LAGRANGE_SEM, P=3");
	
#ifdef NEKTAR_USING_FFTW
	Execute("ADRSolver","Test_UnsteadyAdvectionDiffusion_3DHomo1D_FFT.xml","3D-Homogeneous-1D unsteady advection-diffusion (FFT)");
	
	Execute("ADRSolver","Test_UnsteadyAdvectionDiffusion_3DHomo2D_FFT.xml","3D-Homogeneous-2D unsteady advection-diffusion (FFT)");
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

void RunL2RegressionTest(std::string Demo, std::string input, std::string info)
{
    tests_total++;
    if (!quiet)
    {
        std::cout << "TESTING: " << Demo << std::flush;
    }
	std::string NektarSolverDir =std::string("") +  NEKTAR_SOLVER_DIR + "/bin/";
    RegressBase Test(NektarSolverDir.c_str(),Demo,input,"Solvers/ADRSolver/OkFiles/");
    int fail;

    // Copy input file to current location
	boost::filesystem::path sourceFile(std::string(REG_PATH) + "Solvers/ADRSolver/InputFiles/" + input);
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
    RegressBase Test(NektarSolverDir,Demo,input,"Solvers/ADRSolver/OkFiles/");
    int fail;

    // Copy input file to current location
	boost::filesystem::path sourceFile(std::string(REG_PATH) + "Solvers/ADRSolver/InputFiles/" + input);
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




