#include <cstdio>
#include <cstdlib>
#include <string>
using namespace std;

#ifndef regBase_H
#define regBase_H

#ifndef USR_WIN_32
#define USR_WIN_32 false
#endif //USR_WIN_32

#ifndef MODES_HIGH
#define MODES_HIGH 8
#define MODES_LOW 3
#endif //MODES_HIGH

//header of class regBase
class regBase{
protected:
	// Name of Prog and Path
	std::string Prog;
	
	// number of tests (2: high + low)
	int numTest;
	
	// Base Meshfile
	std::string baseMesh;
	// Low modes and High modes mesh
	std::string Geometry[2];
	int Modes[2];
	
	// Path to Prog, ok and mesh
	std::string Path;
	std::string meshPath;
	std::string okPath;
	
	// Name of .OK File and Output of Prog .ERR File
	std::string Regress;
	std::string Output;
	

public:
	/*-----------------------------------------------------------
	/ CLASS CONSTRUCTORS
	/-----------------------------------------------------------*/
	regBase(void);

	/*-----------------------------------------------------------
	/ PRINT ATTRIBUTES
	/-----------------------------------------------------------*/
	void print(void);
	
	/*-----------------------------------------------------------
	/ SET ATTRIBUTES
	/-----------------------------------------------------------*/
	int setProg(std::string);
	int setPath(std::string);
	int setNumTest(int);
	int setMesh(std::string);
	int setGeometry(std::string[]);
	int setModes(int[]);
	int setOutput(std::string);
	int setRegress(std::string);
	int setOkPath(std::string);
	int setMeshPath(std::string);

	/*-----------------------------------------------------------
	/ Get Attributes
	/-----------------------------------------------------------*/
	std::string getProg(void);
	std::string getPath(void);
	int getNumTest(void);
	std::string getMesh(void);
	std::string getGeometry(int);
	/*-----------------------------------------------------------/
	/ Input: 0, to get meshfile with low number of modes
	/ any non-zero int, to get meshfile with high number of modes.
	/-----------------------------------------------------------*/
	std::string getModes(int);
	// Returns modes[int] as a string
	std::string getOutput(void);
	std::string getRegress(void);
	std::string getOkPath(void);
	std::string getMeshPath(void);

	/*-----------------------------------------------------------
	/ COMPARE ERRORS TO MASTER ERRORS
	/-----------------------------------------------------------*/
	virtual int Test(std::string);
	/*-----------------------------------------------------------
	/ Input: Error String (e.g. "L 2 error")
	/ Returns: 0, if Prog Error is equal to Master Error for both 
	/			  mehsfiles
	/          1, otherwise. 
	/-----------------------------------------------------------*/
	
	virtual int TestL2(void);
	/*-----------------------------------------------------------
	/ Compares the Master Error stored in .ok file to error produced
	/ when running Prog, for both Meshfiles stored in Geometry. 
	/ Returns: 0, if Prog Error is equal to Master Error for both 
	/			  mehsfiles
	/          1, otherwise. 
	/-----------------------------------------------------------*/
	int printTestError(int);
	/*-----------------------------------------------------------
	/ Prints to screen the causes of the error reference by int
	/ as output from Test();
	/-----------------------------------------------------------*/


	/*-----------------------------------------------------------
	/ PROCESS .OK FILE
	/-----------------------------------------------------------*/
	virtual int makeRegress(void);
	/*-----------------------------------------------------------
	/ Generate Master Errors for the Meshfiles in Geometry using 
	/ a known working version of the code.
	/-----------------------------------------------------------*/
	virtual int readRegress(std::string mesh, std::string errStr, std::string &, std::string &);
	/*-----------------------------------------------------------
	/ Inputs: Meshfile, Error String, Errors, Modes
	/ If Data for Meshfile exists in .ok file, reads Error (e.g. L2) and 
	/ stores its value in Errors, and Number of modes in Modes
	/ Returns: 0, if Error and Modes are retirived successfully
	/          1, otherwise. 
	/-----------------------------------------------------------*/

protected:
	virtual std::string makeCommand(int);

	virtual int checkFiles(void);	
	/*-----------------------------------------------------------
	/ Returns: 0, if all files needed are in the right place
	/          1, if Prog.exe is not found
	/		   2, if baseMesh!="" and baseMesh.xml not found
	/          3, if Geometry files not in place or not defined
	/		   4, if .ok file not in place or not defined
	/		   7, if both Geometry and .ok files missing
	/-----------------------------------------------------------*/

	int RunProg(string);
	/*-----------------------------------------------------------
	/ Runs the Prog to be tested, takes command as input.
	/-----------------------------------------------------------*/

	/*-----------------------------------------------------------
	/ EXTRACT FROM FILES
	/-----------------------------------------------------------*/
	std::string getTline(std::string, std::string, bool);
	/*-----------------------------------------------------------
	/ Inputs: Name of File (string)
	/         String contained in line (e.g. "L 2 error")
	/         bool true: use path
	/ Returns: Full line containing input string.
	/-----------------------------------------------------------*/
	std::string getError(std::string);
	/*-----------------------------------------------------------
	/ Inputs:  String that defines error (e.g. "L 2 error")
	/ Returns: Value of error as string (e.g. "0.1234e-004")
	/-----------------------------------------------------------*/
	int getVal(std::string, std::string, std::string &);
	/*-----------------------------------------------------------
	/ Inputs:  String containing Line of file
	/          String with code to look for (e.g. "L 2 error")
	/          String where value is stored
	/ Returns: 0 if Value is extracted successfully
	/          1 otherwise
	/-----------------------------------------------------------*/

	
	/*-----------------------------------------------------------
	/ PROCESS ERRORS
	/-----------------------------------------------------------*/
	double errorToDouble(std::string);
	/*-----------------------------------------------------------
	/ Inputs:  Error formated string as given by getError
	/ Returns: Double Precission containing error 
	/-----------------------------------------------------------*/
	std::string getDig(std::string);
	/*-----------------------------------------------------------
	/ Inputs:  Error formated string as given by getError
	/ Returns: First 4 significant figures 
	/-----------------------------------------------------------*/
	int getExp(std::string, std::string &);
	/*-----------------------------------------------------------
	/ Inputs:  Error formated string as given by getError
	/ Returns: string conatining exponent info (e.g. -004)
	/-----------------------------------------------------------*/
	int compareExp(std::string, std::string);
	/*-----------------------------------------------------------
	/ Inputs:  Error formated string as given by getError
	/ Returns: 0 if exponents are equal, 1 otherwise
	/-----------------------------------------------------------*/
	int isDoubChar(std::string);
	/*-----------------------------------------------------------
	/ Input: Any char
	/ Returns: 1 if char is a number 0-9
	/          2 if char is "e","E","." or "-"
	/          0 otherwise
	/-----------------------------------------------------------*/
	int compareError(std::string, std::string);
	/*-----------------------------------------------------------
	/ Inputs:  Two strings containing errors formated as in getError
	/ Returns: 0 if they are equal
    /          1 otherwise
	/-----------------------------------------------------------*/

	/*-----------------------------------------------------------
	/ PROCESS OUTPUT .ERR FILE
	/-----------------------------------------------------------*/
	virtual int readOutput(std::string errStr, std::string&);

	/*-----------------------------------------------------------
	/ PROCESS GEOMETRY FILES
	/-----------------------------------------------------------*/
	int createGeometry();
	/*-----------------------------------------------------------
	/ Check Mesh and Generate Geometry[]
	/ Output: 0 if successfull
	/		  >0 if unsuccessfull
	/-----------------------------------------------------------*/
	std::string getMeshProperty(std::string,std::string);
	/*-----------------------------------------------------------
	/ Inputs: Meshfile 
	/	      Parameter (e.g. "Projection")
	/ Output: Value (e.g. "Continuous")
	/-----------------------------------------------------------*/
	
	void clearBuffer(char*,int);
};
#endif // regBase_H
