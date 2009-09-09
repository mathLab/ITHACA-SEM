#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "../Auxiliary/regBase.h"

#ifndef REGDEMO_H
#define REGDEMO_H

#define DEMO_MODES_STR "No. modes"
#define DEMO_L2ERROR_STR "L 2 error"
#define DEMO_LINF_ERROR_STR "L infinity error"

using namespace std;
//header of class regDemo
class regDemo: public regBase{
	
public:
	/*-----------------------------------------------------------
	/ CLASS CONSTRUCTORS
	/-----------------------------------------------------------*/
	regDemo();
	regDemo(std::string D);
	regDemo(std::string D, std::string bM);
	regDemo(std::string D, std::string Geo[2]);

	/*---------------------------------------------------------/
	/ Test Errors
	/---------------------------------------------------------*/
	int Test(std::string);
	int TestL2(void);
	int TestLInf(void);
	
	/*---------------------------------------------------------/
	/ PROCESS .OK FILE
	/---------------------------------------------------------*/
	int makeRegress(void);
	int readRegress(std::string mesh, std::string errStr, std::string &L2e,int &modes);
	// Scalar version
	
	std::string demoPath(void);
	// Generates path to MultiRegion Demos
	std::string mPath(void);

protected:
	std::string makeCommand(int);
	
	/*---------------------------------------------------------/
	/ PROCESS .ERR FILE
	/---------------------------------------------------------*/
	int readOutput(std::string errStr, std::string &err);
};

#endif
