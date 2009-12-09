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
	
	// IncNavierStokesSolver tests
	string solverINS = "IncNavierStokesSolver",
	       meshINS1  = "Test_ChaFlow.xml",
	       meshINS2  = "Test_KovFlow.xml";
	
	regSolver TestINS(solverINS,meshINS1);
	
	//Test IncNavierStokesSolver: Channal Flow
	TestINS.makeRegress();
	
	//Test IncNavierStokesSolver: Kovasznay Flow
	TestINS.setMesh(meshINS2);
	TestINS.makeRegress();
	
	// EulerSolver tests
	string solverEUL = "EulerSolver",
	       meshEUL1  = "Test_IsentropicVortex1.xml";
	
	regSolver TestEUL(solverEUL,meshEUL1);
	
	//Test EulerSolver: IsentropicVortex with quads
	TestEUL.makeRegress();

	
	return 0; 
};



	
