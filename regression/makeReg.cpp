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
	
	// Test Helmholtz2D;
	TestH2D.makeRegress();
	
	// Test Advection - Constructed by default
	TestADR.makeRegress();

	// Test DiffusionReaction
	mesh="Test_DiffusionReaction.xml";
	TestADR.setMesh(mesh);
	TestADR.makeRegress();
	
	// Test Implicit Diffusion
	mesh="Test_ImDiffusion.xml";
	TestADR.setMesh(mesh);
	TestADR.makeRegress();

	// Test Explicit Diffusion
	mesh="Test_ExDiffusion.xml";
	TestADR.setMesh(mesh);
	TestADR.makeRegress();

	return 0; 
};



	
