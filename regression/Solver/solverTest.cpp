#include "regSolver.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>

using namespace std;

int main(int argc, char* argv[]) { 
	string solver="AdvectionDiffusionReactionSolver",
		mesh="Test_Advection.xml",
		comm="";
	regSolver TestADR(solver,mesh);
	string result;
	int fail;	
		
	TestADR.print();

	cout<<"\nTesting "<<solver<<"\nPlease wait..."<<endl;

	// Test Advection - Constructed by default
    fail=TestADR.TestL2();
	fail? result="FAIL":result="PASS";
	cout << "\nTest Advection: " << result << endl;
	
	// EXPLAIN CAUSES FOR ERROR, IF ANY
	TestADR.printTestError(fail);

	// Test DiffusionReaction
	mesh="Test_DiffusionReaction.xml";
	TestADR.setMesh(mesh);
	fail=TestADR.TestL2();
	fail? result="FAIL":result="PASS";
	cout << "\nTest DiffusionReaction: " << result << endl;
	
	// EXPLAIN CAUSES FOR ERROR, IF ANY
	TestADR.printTestError(fail);
	
	// Test Implicit Diffusion
	mesh="Test_ImDiffusion.xml";
	TestADR.setMesh(mesh);
	fail=TestADR.TestL2();
	fail? result="FAIL":result="PASS";
	cout << "\nTest Implicit Diffusion: " << result << endl;

	// EXPLAIN CAUSES FOR ERROR, IF ANY
	TestADR.printTestError(fail);

	// Test Explicit Diffusion
	mesh="Test_ExDiffusion.xml";
	TestADR.setMesh(mesh);
	fail=TestADR.TestL2();
	fail? result="FAIL":result="PASS";
	cout << "\nTest Explicit Diffusion: " << result << endl;

	// EXPLAIN CAUSES FOR ERROR, IF ANY
	TestADR.printTestError(fail);

	return 0; 
};



	
