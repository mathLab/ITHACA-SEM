#include "Demo/regDemo.h"
#include "Solver/regSolver.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>
#include <sys/stat.h>

using namespace std;

int main(int argc, char* argv[]) { 
    
    regDemo TestH2D;
    string solver="AdvectionDiffusionReactionSolver",
        mesh="Test_Advection.xml";
    regSolver TestADR(solver,mesh);
    string result;
    int fail;   
    
    // TEST HELMHOLTZ2D Demo
    cout << "\nTesting Helmholtz2D Demo..." ;
    fail=TestH2D.TestL2();
    fail? result="FAIL":result="PASS";
    cout << result << endl; 

    // TEST ADVECTION - Constructed by default
    cout << "\nTesting "<< solver <<":\n\n\t"<<mesh<<"..." ;
    // TEST L2 ERROR FOR ALL VARIABLES
    fail=TestADR.TestL2();
    fail? result="FAIL":result="PASS";
    cout << result << endl;

/*
    // TEST DIFFUSIONREACTION
    mesh="Test_DiffusionReaction.xml";
    cout <<"\n\t"<<mesh<<"..." ;
    TestADR.setMesh(mesh);
    // TEST L2 ERROR FOR ALL VARIABLES
    fail=TestADR.TestL2();
    fail? result="FAIL":result="PASS";
    cout << result << endl;
    
    // TEST IMPLICIT DIFFUSION
    mesh="Test_ImDiffusion.xml";
    cout <<"\n\t"<<mesh<<"..." ;
    TestADR.setMesh(mesh);
    // TEST L2 ERROR FOR ALL VARIABLES
    fail=TestADR.TestL2();
    fail? result="FAIL":result="PASS";
    cout << result << endl;
*/

    // TEST EXPLICIT DIFFUSION
    mesh="Test_ExDiffusion.xml";
    cout <<"\n\t"<<mesh<<"..." ;
    TestADR.setMesh(mesh);
    // TEST L2 ERROR FOR ALL VARIABLES
    fail=TestADR.TestL2();
    fail? result="FAIL":result="PASS";
    cout << result << endl;
	
	//////////////////////////////////
	// INCOMPRESSIBLE NAVIER-STOKES SOLVER
	
	solver = "IncNavierStokesSolver";
	mesh   = "Test_ChaFlow.xml";
    
	regSolver TestINS(solver,mesh);
	
	// TEST Channal Flow
	cout << "\nTesting "<< solver <<":\n\n\t"<<mesh<<"..." ;
    // TEST L2 ERROR FOR ALL VARIABLES
    fail=TestINS.TestL2();
    fail? result="FAIL":result="PASS";
    cout << result << endl;
	
	// TEST Kovasznay Flow 
    mesh = "Test_KovFlow.xml";
    cout <<"\n\t"<<mesh<<"..." ;
    TestINS.setMesh(mesh);
    // TEST L2 ERROR FOR ALL VARIABLES
    fail=TestINS.TestL2();
    fail? result="FAIL":result="PASS";
    cout << result << endl;
	
	//////////////////////////////////
	// EULER SOLVER
	
	solver = "EulerSolver";
	mesh   = "Test_IsentropicVortex1.xml";
	
	regSolver TestEUL(solver,mesh);
	
	// TEST Isentropic vortex
	cout << "\nTesting "<< solver <<":\n\n\t"<<mesh<<"..." ;
	// TEST L2 ERROR FOR ALL VARIABLES
	fail=TestEUL.TestL2();
	fail? result="FAIL":result="PASS";
	cout << result << endl;
	
    return 0; 
};



    
