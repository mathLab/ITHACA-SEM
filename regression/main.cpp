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

    return 0; 
};



    
