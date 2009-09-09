#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "../Auxiliary/regBase.h"

#ifndef regSolver_H
#define regSolver_H

#define SOLVER_MAX_VARS 4
#define SOLVER_MODES_STR "Max Exp. Order"
#define SOLVER_L2ERROR_STR "L2 Error"
#define SOLVER_TIME_STR "Time Step"

using namespace std;
//header of class regSolver
class regSolver: public regBase{
public:
	/*-----------------------------------------------------------
	/ CLASS CONSTRUCTORS
	/-----------------------------------------------------------*/
	regSolver(void);
	regSolver(std::string D);
	regSolver(std::string D, std::string bM);
	regSolver(std::string D, std::string Geo[2]);
	
	/*---------------------------------------------------------/
	/ Test Errors
	/---------------------------------------------------------*/
	int Test(std::string);
	int TestL2(void);
	
	/*---------------------------------------------------------/
	/ Process .OK File
	/---------------------------------------------------------*/
	int makeRegress(void);
	int readRegress(std::string mesh, std::string errStr, std::string *, std::string *, int &);
	/*---------------------------------------------------------/
	/ arguments: meshfile, error string, vars[], errors[], modes
	/---------------------------------------------------------*/

	std::string solverPath(void);
	std::string mPath(void);

protected:
	std::string makeCommand(int);
	
	/*---------------------------------------------------------/
	/ PROCESS OUTPUT .ERR FILE
	/----------------------------------------------------------*/
	int readOutput(std::string errStr, std::string *,std::string *);
	/*---------------------------------------------------------/
	/ arguments: vars[], errors[]
	/---------------------------------------------------------*/
};
#endif // regSolver_H
