#include "regDemo.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) { 
	string demo="Helmholtz2D";
	string mesh="helmholtz2D.xml";
	regDemo TestH2D(demo,mesh);
	string result;
	int fail;	
	
	TestH2D.print();
	fail=TestH2D.TestL2();
	fail? result="FAIL":result="PASS";
	cout << "\nTest Helmholtz2D: " << result << endl<<endl;
	
	// EXPLAIN CAUSES FOR ERROR, IF ANY
	TestH2D.printTestError(fail);
	return 0; 
};



	
