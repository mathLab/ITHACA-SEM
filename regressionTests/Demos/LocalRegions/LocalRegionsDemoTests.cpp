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

int main(int argc, char* argv[]) 
{ 

    // 1D Projection Tests
    Execute("LocProject1D", "1 6 7 0.0 1.0","Segment [0,1] Ortho Basis, modes=6, quad=7");

    Execute("LocProject1D", "4 6 7 0.0 1.0","Segment [0,1] Modified Basis, modes=6, quad=7");

    // 1D Differentiation Projection Tests
    Execute("LocProject_Diff1D", "1 6 7 0.0 1.0","Segment [0,1] Ortho Basis, modes=6, quad=7");

    Execute("LocProject_Diff1D", "4 6 7 0.0 1.0","Segment [0.1] Modified Basis, modes=6, quad=7");

    // 2D Projection Tests
    Execute("LocProject2D", "2 1 2 6 6 7 7 0.0 0.0 1.0 1.0 0.5 1.0", "Triangular Ortho Basis, modes=6, quad=7");

    Execute("LocProject2D", "2 4 5 6 6 7 7 0.0 0.0 1.0 1.0 0.5 1.0", "Triangular Modified Basis modes=6, quad=7");

    Execute("LocProject2D", "2 11 11 6 6 7 7  0.0 0.0 1.0 1.0 0.5 1.0", "Triangular Nodal Basis modes=6, quad=7");

    Execute("LocProject2D", "3 1 1 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Regular Quadrilateral Ortho Basis, modes=4, quad=5");

    Execute("LocProject2D", "3 4 4 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Regular Quadrilateral Modified Basis modes=6, quad=7");

    Execute("LocProject2D", "3 8 8 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Regular Quadrilateral Nodal/Lagrange Basis modes=6, quad=7");

    Execute("LocProject2D", "3 1 1 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Linearly Deformed Quadrilateral Ortho Basis, modes=4, quad=5");

    Execute("LocProject2D", "3 4 4 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Linearly Deformed Quadrilateral Modified Basis modes=6, quad=7");

    Execute("LocProject2D", "3 8 8 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Linearly Deformed Quadrilateral Nodal/Lagrange Basis modes=6, quad=7");

    // 2D Differntiation and Projection Tests
    Execute("LocProject_Diff2D", "2 1 2 6 6 7 7 0.0 0.0 1.0 1.0 0.5 1.0", "Triangular Ortho Basis, modes=6, quad=7");

    Execute("LocProject_Diff2D", "2 4 5 6 6 7 7 0.0 0.0 1.0 1.0 0.5 1.0", "Triangular Modified Basis modes=6, quad=7");

    Execute("LocProject_Diff2D", "2 11 11 6 6 7 7  0.0 0.0 1.0 1.0 0.5 1.0", "Triangular Nodal Basis modes=6, quad=7");

    Execute("LocProject_Diff2D", "3 1 1 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Regular Quadrilateral Ortho Basis, modes=4, quad=5");

    Execute("LocProject_Diff2D", "3 4 4 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Regular Quadrilateral Modified Basis modes=6, quad=7");

    Execute("LocProject_Diff2D", "3 8 8 6 6 7 7 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0", "Regular Quadrilateral Nodal/Lagrange Basis modes=6, quad=7");

    Execute("LocProject_Diff2D", "3 1 1 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Linearly Deformed Quadrilateral Ortho Basis, modes=4, quad=5");

    Execute("LocProject_Diff2D", "3 4 4 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Linearly Deformed Quadrilateral Modified Basis modes=6, quad=7");

    Execute("LocProject_Diff2D", "3 8 8 6 6 7 7 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0", "Linearly Deformed Quadrilateral Nodal/Lagrange Basis modes=6, quad=7");

    
    // 3D Projection Tests
    Execute("LocProject3D", "4 1 2 3 6 6 6 7 7 7 0 0 0  1 0 0  0 1 0  0 0 1", "Regular Tetrahedron Ortho Basis, modes=6, quad=7");
    
    Execute("LocProject3D", "4 4 5 6 6 6 6 7 7 7 0 0 0  1 0 0  0 1 0  0 0 1", "Regular Tetrahedron Modified Basis, modes=6, quad=7");

    Execute("LocProject3D", "7 1 1 1 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Regular Hexahedron Ortho Basis, modes=6, quad=7");
    
    Execute("LocProject3D", "7 4 4 4 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Regular Hexahedron Modified Basis, modes=6, quad=7");

    Execute("LocProject3D", "7 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Regular Hexahedron Nodal/Lagrange Basis, modes=6, quad=7");
    
    Execute("LocProject3D", "7 1 1 1 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Linearly Deformed Hexahedron Ortho Basis, modes=6, quad=7");
    
    Execute("LocProject3D", "7 4 4 4 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Linearly Deformed Hexahedron Modified Basis, modes=6, quad=7");

    Execute("LocProject3D", "7 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Linearly Deformed Hexahedron Nodal/Lagrange Basis, modes=6, quad=7");

    // 3D Differentiation and Projection Tests
    Execute("LocProject_Diff3D", "4 1 2 3 6 6 6 7 7 7 0 0 0  1 0 0  0 1 0  0 0 1", "Regular Tetrahedron Ortho Basis, modes=6, quad=7");
    
    Execute("LocProject_Diff3D", "4 4 5 6 6 6 6 7 7 7 0 0 0  1 0 0  0 1 0  0 0 1", "Regular Tetrahedron Modified Basis, modes=6, quad=7");

    Execute("LocProject_Diff3D", "7 1 1 1 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Regular Hexahedron Ortho Basis, modes=6, quad=7");
    
    Execute("LocProject_Diff3D", "7 4 4 4 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Regular Hexahedron Modified Basis, modes=6, quad=7");

    Execute("LocProject_Diff3D", "7 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1 0  0 1 0  0 0 1  1 0 1  1 1 1  0 1 1", "Regular Hexahedron Nodal/Lagrange Basis, modes=6, quad=7");
    
    Execute("LocProject_Diff3D", "7 1 1 1 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Linearly Deformed Hexahedron Ortho Basis, modes=6, quad=7");
    
    Execute("LocProject_Diff3D", "7 4 4 4 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Linearly Deformed Hexahedron Modified Basis, modes=6, quad=7");

    Execute("LocProject_Diff3D", "7 8 8 8 6 6 6 7 7 7 0 0 0  1 0 0  1 1.5 0  0 1 0  0 0 1  1.5 0 1  1 1 1  0 1 1.5", "Linearly Deformed Hexahedron Nodal/Lagrange Basis, modes=6, quad=7");
    
    return 0;
}

void RunL2RegressionTest(std::string Demo, std::string input, std::string info)
{
    RegressBase Test(NEKTAR_BIN_DIR,Demo,input,"Demos/LocalRegions/OkFiles/");
    int fail;

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
};

void MakeOkFile(std::string Demo, std::string input, std::string info)
{
    RegressBase Test(NEKTAR_BIN_DIR,Demo,input,"Demos/LocalRegions/OkFiles/");
    int fail;

    if(fail = Test.MakeOkFile())
    {
        std::cout << "Failed to make OK file\n";
        // Explain cause of error if available
        std::cout << "===========================================================\n";   
        Test.PrintTestError(fail);
        std::cout << "===========================================================\n";    
}
}



	
