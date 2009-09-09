#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <sys/stat.h>
#include "regBase.h"

using namespace std;

regBase::regBase(){
		Prog="";
		Path="";
		numTest=2;
		baseMesh="";
		Geometry[0]=""; //  LOW ORDER EXPANSION 
		Geometry[1]=""; //  HIGH ORDER EXPANSION
		Modes[0]=MODES_LOW; //  Number of Modes {High, Low}
		Modes[1]=MODES_HIGH;  
		Regress=Prog+".ok";
		Output=Prog+".err";
	}

void regBase::print(){
	std::cout << "===============================================\n";
	std::cout << "Program: " << regBase::Prog <<std::endl;
	std::cout << "Path: \n" << regBase::Path <<std::endl;
	std::cout << "okPath: \n" << regBase::okPath <<std::endl;
	std::cout << "meshPath: \n" << regBase::meshPath <<std::endl;
	std::cout << "Geometry: " <<std::endl;
	for(int i=0;i<regBase::numTest;++i){
		if(regBase::Geometry[i]!=""){
			std::cout << "\t" << regBase::Geometry[i];
			std::cout <<" (Modes=" << regBase::Modes[i] << ")" << std::endl;
		}
	}
	std::cout << "Master and Output: " <<std::endl;
	std::cout << "\t" << regBase::Output <<std::endl;
	std::cout << "\t" << regBase::Regress <<std::endl;
	std::cout << "===============================================\n";
};

/*---------------------------------------------------------/
/ Set Attributes
/---------------------------------------------------------*/
int regBase::setProg(std::string D){
	regBase::Prog=D;
	return 0;
};

int regBase::setPath(std::string P){
	regBase::Path=P;
	return 0;
};

int regBase::setMesh(std::string M){
	int pos=M.find(".");
	std::string msh=M.substr(0,pos);
	regBase::baseMesh=M;
	regBase::Geometry[0]=msh+"_low.xml"; 
	regBase::Geometry[1]=msh+"_high.xml";
	return regBase::createGeometry();
};

int regBase::setNumTest(int n){
	regBase::numTest=n;
	return 0;
};

int regBase::setGeometry(std::string G[]){
	for (int i=0;i<regBase::numTest;++i){
		regBase::Geometry[i]=G[i];
	}
	return 0;
};			

int regBase::setModes(int m[]){
	for(int i=0;i<regBase::numTest;++i){
		regBase::Modes[i]=m[i];
	}
	// UPDATE GEOMETRY FILES AND RETURN OUTCOME
	return regBase::createGeometry();
};

int regBase::setOutput(std::string out){
	regBase::Output=out;
	return 0;
};

int regBase::setRegress(std::string reg){
	regBase::Regress=reg;
	return 0;
};

int regBase::setOkPath(std::string op){
	regBase::okPath=op;
	return 0;
};

int regBase::setMeshPath(std::string mp){
	regBase::meshPath=mp;
	return 0;
};

/*---------------------------------------------------------/
/ Get Attributes
/---------------------------------------------------------*/
std::string regBase::getProg(void){
	return regBase::Prog;
};

std::string regBase::getPath(void){
	return regBase::Path;
};

int regBase::getNumTest(void){
	return regBase::numTest;
};

std::string regBase::getGeometry(int i){
	// Check that i is between 0 and num_test
	if(i>=regBase::numTest){i=regBase::numTest-1;}
	if(i<0){i=0;}
	return regBase::Geometry[i];
};

std::string regBase::getModes(int i){
	std::string num="";
	std::stringstream out;
	if (i!=0){i=1;}
	out << regBase::Modes[i];
	num = out.str();
	return num;
};

std::string regBase::getOutput(void){
	return regBase::Output;
};

std::string regBase::getRegress(void){
	return regBase::Regress;
};

std::string regBase::getOkPath(void){
	return regBase::okPath;
};

std::string regBase::getMeshPath(void){
	return regBase::meshPath;
};

/*---------------------------------------------------------/
/ Test Error
/---------------------------------------------------------*/
int regBase::Test(std::string errStr){
	return 1;
}
int regBase::TestL2(void){	
	return 1;
};
int regBase::printTestError(int fail){
	switch(fail){
		case 0: break;
		case 1:{
			cout << "\n* The file "<< regBase::Prog << "could not be found!"<<endl;
			cout << "* Suggestions:\n\tCheck Path:\n\t"<<regBase::Path<<endl;
			break;
		   };
		case 2:{
			cout << "\n* The file "<< regBase::baseMesh << "could not be found!"<<endl;
			cout << "* Suggestions:\n\tCheck Path:\n\t"<<regBase::Path<<endl;
			break;
		   };
		case 3:{
			cout << "\n* The file "<< regBase::Regress << "could not be found!"<<endl;
			cout << "* Suggestions:\n\tCheck Path:\n\t"<<regBase::Path<<endl;
			break;
		   };
		case 4:{
			cout << "\n* The number of Expansion Modes (High AND/OR Low) disagree with info stored in .OK file!"<<endl;
			break;
		   };
		case 5:{
		   cout << "\n* There is no info for this mesh in .OK file! (High AND/OR Low Modes)"<<endl;
		   break;
		   };
		case 6:{
		   cout << "\n* The file "<< regBase::Output << "could not be found!"<<endl;
			cout << "* Suggestions:\n\tCheck Path:\n\t"<<regBase::Path<<endl;
		   break;
		   };
		case 7:{
		   cout << "\n* The error obtained for one of the Input Files tested does not agree with the value stored on the .OK file!"<<endl;
		   break;
		   };
	};
	return 0;
};
/*---------------------------------------------------------/
/ Process .ok File
/---------------------------------------------------------*/
int regBase::makeRegress(){
	// Not implemented
	return 1;
};

int regBase::readRegress(std::string mesh, std::string errStr, std::string &L2e,std::string &modes){
	// Not implemented
	return 1;
};



// PRIVATE METHODS
/*---------------------------------------------------------/
/ Check Files and Run Program
/---------------------------------------------------------*/
int regBase::checkFiles(){
	// CHECK EXISTENCE OF FILES	//
	int fail=0;
	struct stat buf;
	std::string fname="";

	// Check Prog
	fname=regBase::Path+regBase::Prog; 
	if (USR_WIN_32){
		fname=fname+".exe";
	};
	if (stat(fname.c_str(),&buf)){
		// Prog is not in folder
		cout << fname << " is not available!" << endl;
		std::cout << "Please check Path is correct: \n\t"<<regBase::Path<<std::endl;
		return 1;
	}

	// Check meshfiles
	if(regBase::baseMesh!=""){
		fname=regBase::meshPath+regBase::baseMesh;
		if (stat(fname.c_str(),&buf)){
			// baseMesh is not in folder
			cout << fname << " is not available!" << endl;
			std::cout << "Please check Path is correct: \n\t"<<regBase::meshPath<<std::endl;
			return 2;
		}
	}
	for (int i=0; i<regBase::numTest;++i){
		fname=regBase::Path+regBase::Geometry[i];
		if (stat(fname.c_str(),&buf)){
			// Geometry is not in folder or has not been generated
			fail=3;
		}
	}
	// Check Regress
	fname=regBase::okPath+regBase::Regress;
	if (stat(fname.c_str(),&buf)){
		// Regress is not in folder
		fail+=4;
	}
	return fail;
};		

std::string regBase::makeCommand(int i){
	std::string prog="",mesh="",comm="";
	prog=regBase::Path+regBase::Prog;
	mesh=regBase::Path+regBase::Geometry[i];
	comm=prog+" "+mesh+" "+mesh;
	comm+=" > "+regBase::Output;
	return comm;
}

int regBase::RunProg(std::string comm){			
// Run Prog with Geometry files, storing output in Output file.	//
	int status;
	status=system(comm.c_str());
	if (status) {
		std::cout << "Unable to run command " << comm << std::endl;
		return 1;
	}
	return 0;
};
/*---------------------------------------------------------/
/ Auxiliary methods to process errors
/---------------------------------------------------------*/
std::string regBase::getTline(std::string source, std::string T, bool usePath){
	// Extract the line containing T Error from Output of Prog. //	
	int lineNo=0;
	FILE *res;
	char buffer[40];
	std::string line="",fname=source;
	bool gotit=false;
	//Open file Source
	if(usePath){fname=regBase::Path+source;}
	res = fopen(fname.c_str(), "r");
	if (!res) {
		std::cout << "Could not open output file (" << source << ")!" << std::endl;
		return line;
	}
	// Process the output of Prog.
	while(fgets(buffer,sizeof(buffer),res) && !gotit){
		++lineNo;
		line.assign(buffer,sizeof(buffer)-1);
				
		// Clear buffer
		regBase::clearBuffer(buffer,sizeof(buffer));

		// Retrieve line containing L2 Error
		if (line.find(T)<std::string::npos){
			gotit=true;
		}
	} // while
	fclose(res);
	if(!gotit){line="";}
	return line;
};

std::string regBase::getError(std::string l2){
	int pos=15; //initial position
	std::string err="";
	bool stop=false,got1=false;
	if (l2.find("L 2 error")<std::string::npos){
		while(!stop){
			if(regBase::isDoubChar(l2.substr(pos,1))){
				err.append(l2.substr(pos,1));
				got1=true;
			}
			else{
				if (got1){
					stop=true;
				}
			}
			++pos;
		}
		//std::cout << "getError: " << err << std::endl;
		//err=l2.substr(18,12);
	}
	return err;
};

double regBase::errorToDouble(std::string error){
	double err_d=0.0;
	err_d=std::strtod(error.c_str(),NULL);
	return err_d;
};
std::string regBase::getDig(std::string error){
	return error.substr(0,6);
};

int regBase::getExp(std::string error, std::string &exp){
	// Assuming error is formated using getError or getVal
	int pos=0;
	exp="";
	if (error.find("e")<std::string::npos){
		pos=error.find("e")+1;
		exp=error.substr(pos,std::string::npos);
		return 0;
	}
	if (error.find("E")<std::string::npos){
		pos=error.find("E")+1;
		exp=error.substr(pos,std::string::npos);
		return 0;
	}
	return 1;
};

int regBase::isDoubChar(std::string c){
	bool isd=false;	
	isd=(c=="0")||(c=="1")||(c=="2")||(c=="3")||(c=="4")||(c=="5");
	isd=isd||(c=="6")||(c=="7")||(c=="8")||(c=="9");
	if (isd){
		return 1;
	}
	isd=isd||(c==".")||(c=="e")||(c=="E")||(c=="-");
	if (isd){
		return 2;
	}
	else{
		return 0;
	};
};
int regBase::getVal(std::string l, std::string code, std::string &val){
	int pos=0; //initial position
	std::string line="";
	bool stop=false,got1=false;
	val="";
	if (l.find(code)<std::string::npos){
		pos=l.find(":");
		line=l.substr(pos,std::string::npos);
		pos=0;
		while(!stop && !(pos>line.length())){
			if(regBase::isDoubChar(line.substr(pos,1))){
				val.append(line.substr(pos,1));
				if ((line.substr(pos,1)!=".")&&(line.substr(pos,1)!="e")){
					got1=true;
				}
			}
			else{
				if (got1){
					stop=true;
				}
			}
			++pos;
		}
	}
	if (val==""){return 1;}
	return 0;
}
int regBase::compareExp(std::string e1, std::string e2){
	int intE1,intE2;
	if (e1!=""){
		if (e2==""){
			// e1!=e2
			return 1;
		}
		intE1=std::atoi(e1.c_str());
		intE2=std::atoi(e2.c_str());
		if (intE1==intE2){
			return 0;
		}		
		else{
			return 1;
		}
	}
	else{
		if (e2==""){
			return 0;
		}		
		else{
			return 1;
		}
	}
};

/*---------------------------------------------------------/
/ Compare Errors
/---------------------------------------------------------*/
int regBase::compareError(std::string e1, std::string e2){
	double err, errM;

	/////// CONVERT ERRORS FROM STRING TO DOUBLE /////////
	errM=regBase::errorToDouble(e1);
	err=regBase::errorToDouble(e2);
	if(errM!=err){
		return 1;
	}
	return 0;
	/*
		// std::string f6="",exp=""; //first 5 digits "0.1234" and exp "e-004"
		// std::string f6M="",expM=""; //first 5 digits and exp of Master Error
		
		// Get first 6 digits and exp from Master errors
		f6M=regBase::getDig(L2);
		regBase::getExp(L2,expM);

		int pos=error.find(f6M,0);
		if (pos!=0){
			std::cout << "First 5 Digits of Prog Error and Master Error don't match!"<<std::endl;
			++fail;
			break;
		}
		
		//////// COMPARE EXPONENTS  ///////////
		regBase::getExp(error,exp);
		if(regBase::compareExp(expM, exp)){
			std::cout << "Exponents do not correspond!"<<std::endl;
			++fail;
		}

*/	
};




/*---------------------------------------------------------/
/ PROCESS .ERR FILE
/---------------------------------------------------------*/
int regBase::readOutput(std::string errStr, std::string &err){
	// Not implemented
	return 1;	
};

/*-----------------------------------------------------------
/ PROCESS GEOMETRY FILES
/-----------------------------------------------------------*/
int regBase::createGeometry(){
	std::string low="",high="";
	FILE *lowF,*highF,*bm;
	std::string fname="",line="",lowL="",highL="";
	std::string before="", after="";
	int p1=0,p2=0;
	bool inExp=false;
	char buffer[50];
	struct stat buf;

	// Retrieve Number or Modes for Low and High
	low=regBase::getModes(0);	
	high=regBase::getModes(1);
	
	// Check if Geometry Files need to be Created
	if(regBase::baseMesh==""){
		return 1;
	}
	
	// Check if baseMesh File exists
	fname=regBase::meshPath+regBase::baseMesh;
	if (stat(fname.c_str(),&buf)){
		// Prog is not in folder
		cout << fname << " is not available!" << endl;
		return 2;
	}

	// Open baseMesh File
	bm=fopen(fname.c_str(),"r");
	if(!bm){
		// Unable to open file baseMesh
		return 2;
	}
	
	// Open Geometry Files (Create them if they don't exist)
	fname=regBase::Geometry[0];	
	lowF=fopen(fname.c_str(),"w");
	fname=regBase::Geometry[1];
	highF=fopen(fname.c_str(),"w");

	// Read from bm
	while(fgets(buffer, sizeof(buffer), bm)){
		line.assign(buffer,sizeof(buffer)-1);

		if(!inExp){ // Copy line to lowF and highF
			fputs(line.c_str(),lowF);
			fputs(line.c_str(),highF);
		}
		else{
			// Modify number of modes in line
			p1=line.find("NUMMODES=")+sizeof("NUMMODES=");
			
			if (p1<std::string::npos){
				p2=line.find("\"",p1);
				before=line.substr(0,p1);
				after=line.substr(p2,std::string::npos);

				// Copy new line to lowF and highF
				lowL=before+low+after;
				fputs(lowL.c_str(),lowF);
				
				highL=before+high+after;
				fputs(highL.c_str(),highF);
			}
		}
		if (line.find("EXPANSIONS")<std::string::npos){
			// We are reading expansions
			inExp=true;
		}
		if (line.find("/EXPANSIONS")<std::string::npos){
			// There are no more lines to be changed
			inExp=false;
		}	
	}//while
	
	// Close Files
	fclose(bm);
	fclose(lowF);
	fclose(highF);
	return 0;
};

std::string regBase::getMeshProperty(std::string filename,std::string prop){
	FILE *bm;
	std::string line="",VAL="";
	int p1=0,p2=0;
	char buffer[70];
	struct stat buf;
	std::string fname=filename;

	// Check if filename exists
	if (stat(fname.c_str(),&buf)){
		// Prog is not in folder
		cout << fname << " is not available!" << endl;
		return VAL;
	}

	// Open baseMesh File
	bm=fopen(fname.c_str(),"r");
	if(!bm){
		// Unable to open file baseMesh
		return VAL;
	}

	// Read from bm
	while(fgets(buffer, sizeof(buffer), bm)){
		line.assign(buffer,sizeof(buffer)-1);
		if(line.find(prop)<std::string::npos){ 
			// Get Value
			p1=line.find("VALUE=")+sizeof("VALUE=");
		
			if (p1<std::string::npos){
				p2=line.find("\"",p1);
				VAL=line.substr(p1,p2-p1);
			}
		}
	}//while
	
	// Close Files
	fclose(bm);
	return VAL;
};

void regBase::clearBuffer(char* buf, int bSize){
	int i=0; 
	char empty=0;
	char *bufArray=buf;
	for(int i=0;i<bSize;i++){
		*bufArray=empty;
		++bufArray;
	}
};
