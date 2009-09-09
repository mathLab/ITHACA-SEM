#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "regSolver.h"

using namespace std;

/*-----------------------------------------------------------
/ CLASS CONSTRUCTORS
/-----------------------------------------------------------*/
regSolver::regSolver(){
	Prog="AdvectionDiffusionReactionSolver";
	baseMesh="";
	numTest=2;
	// Default : Test Advection
	regSolver::Geometry[0]="Test_Advection_low.xml"; 
	regSolver::Geometry[1]="Test_Advection_high.xml";
	Modes[0]=MODES_LOW;  //  Number of Modes {High, Low}
	Modes[1]=MODES_HIGH;
	Output=Prog+".err";
	Regress=Prog+".ok";
	// SOLVER and .OK PATH
	Path=regSolver::solverPath();
	okPath.assign(REG_PATH);
	okPath+="/Solver/";
	meshPath=regSolver::mPath();
};

regSolver::regSolver(std::string D, std::string bM){
	std::string msh="";
	int pos;

	Prog=D;
	numTest=2;
	baseMesh=bM;
	pos=bM.find(".");
	msh=bM.substr(0,pos);
	Geometry[0]=msh+"_low.xml"; 
	Geometry[1]=msh+"_high.xml";
	Modes[0]=MODES_LOW;  //  Number of Modes {High, Low}
	Modes[1]=MODES_HIGH;
	Output=D+".err";
	Regress=D+".ok";
	// SOLVER and .OK PATH
	Path=regSolver::solverPath();
	okPath.assign(REG_PATH);
	okPath+="/Solver/";
	meshPath=regSolver::mPath();
}

regSolver::regSolver(std::string D, std::string Geo[2]){
	Prog=D;
	numTest=2;
	baseMesh="";
	Geometry[0]=Geo[0]; 
	Geometry[1]=Geo[1];
	Modes[0]=MODES_LOW;  //  Number of Modes {High, Low}
	Modes[1]=MODES_HIGH;
	Output=D+".err";
	Regress=D+".ok";

	// SOLVER and .OK PATH
	Path=regSolver::solverPath();
	okPath.assign(REG_PATH);
	okPath+="/Solver/";
	meshPath=regSolver::mPath();
};
/*---------------------------------------------------------/
/ Test Errors
/---------------------------------------------------------*/
int regSolver::Test(std::string errStr){
	std::string outErrors[SOLVER_MAX_VARS];
	std::string outVars[SOLVER_MAX_VARS];
	std::string regErrors[SOLVER_MAX_VARS];
	std::string regVars[SOLVER_MAX_VARS];
	int modes=0, fail=0;
	int i=0,n=0;

	// CHECKFILES 
	int chk=regSolver::checkFiles();
	switch (chk){
		case 1:case 2: {
			// Prog or baseMesh (if defined) not in place
			return chk;
				  }
		case 4:case 7:{
			// .OK file Does not exist or is not in place.
			return 3;
			   }
		default:{
			if(regSolver::baseMesh!=""){
				// CREATE GEOMETRY FILES	
				regSolver::createGeometry();
			}
			// 0 -> No problem, overwrite Geometry			
			// 3 -> Geometry files don't exist
				}
	}

	for(i=0;i<regSolver::numTest;++i){
		// READ REGRESS .OK FILE
		if(regSolver::readRegress(regSolver::Geometry[i],errStr,regVars,regErrors,modes)){
			fail=4;
			//cout << regSolver::Regress << " not available!" << endl;
			continue;
		}
		
		// CHECK MODES
		if(modes!=regSolver::Modes[i]){
			cout <<"\nModes in " << regSolver::Geometry[i];
			cout << ": " << modes << " :: Modes in .ok: " << regSolver::Modes[i] << endl;
			fail=5;
			continue;
		}

		// RUN SOLVER WITH GEOMETRY
		regSolver::RunProg(regSolver::makeCommand(i));

		// READ OUTPUT .ERR FILE
		if(regSolver::readOutput(errStr,outVars,outErrors)){
			fail=6;
			continue;
		};

		for(n=0;n<SOLVER_MAX_VARS;++n){
			if(outVars[n].length()){
				// COMPARE ERRORS FROM .ERR AND .OK		
				if(regSolver::compareError(outErrors[n],regErrors[n])){
					cout <<"\n variable "<<outVars[n]<<" (modes="<<regSolver::Modes[i];
					cout << ") error: " << outErrors[n] << " :: .ok error: " << regErrors[n] << endl;
					fail=7;
				};		
			}
		} // for Vars
	}// for Geo

	return fail;
};

int regSolver::TestL2(void){
	return regSolver::Test(SOLVER_L2ERROR_STR);
};

/*---------------------------------------------------------/
/ Process .OK File
/---------------------------------------------------------*/
int regSolver::makeRegress(void){
	std::string command;
	std::string lineT="",line="",fname="";
	std::string prop="Projection", val="";
	std::string tName=regSolver::okPath+"temp.ok";	
	FILE *mRes,*out,*temp;	
	char buffer[40];
	bool errLine=false,gotMesh=false,gotProj=false,newGeo=false,stop=false;
	
	// CHECKFILES 
	int chk=regSolver::checkFiles();
	switch (chk){
		case 1:case 2: {
			// Prog or baseMesh (if defined) not in place
			return chk;
				  }
		default:{
			if(regSolver::baseMesh!=""){
				// CREATE GEOMETRY FILES			
				regSolver::createGeometry();
			}
			// 0,4 -> No problem, overwrite Geometry			
			// 3,7 -> Geometry files don't exist
				}
	}
	
	// INFORM USER THAT MAKEREGRESS IS RUNNING
	std::cout << "\nGenerating Regression Test Resutls for\n\t" << regSolver::Prog << std::endl;
	for(int i=0;i<regSolver::numTest;++i){
		std::cout << "\t" << regSolver::Geometry[i];
		std::cout <<" (Modes=" << regBase::Modes[i] << ")" << std::endl;
	}
	std::cout<< "Please Wait..."<<std::endl;

	// CREATE REGRESS.OK (IF NEEDED)
	fname=regSolver::okPath+regSolver::Regress;
	if ((chk==4)||(chk==7)){
		// Regress is not in folder
		// Create Regress.OK and add solver name
		mRes=fopen(fname.c_str(),"w+");
		if(!mRes){
			// Unable to open file
			return 3;
		}
		fputs(regSolver::Prog.c_str(),mRes);
		fputs("\n==========================================\n",mRes);
		fclose(mRes);
	}
	
	// FOR EACH GEOMETRY EXTRACT AND UPDATE OR INSERT INFO INTO OK FILE
	for (int i=0;i<regSolver::numTest;++i){
		
		// OPEN TEMPORARY .OK FILE
		temp=fopen(tName.c_str(),"w");
		
		// OPEN .OK FILE TO READ 
		mRes=fopen(fname.c_str(),"r");

		if(!temp || !mRes){
			// Problems opening files
			return 3;
		}

		// COPY .OK FILE TO TEMP UNTIL GEOMETRY IS FOUND
		gotMesh=false;
		while(fgets(buffer,sizeof(buffer),mRes)){
			lineT.assign(buffer,sizeof(buffer)-1);
			
			// Clear buffer
			regSolver::clearBuffer(buffer,sizeof(buffer));
			
			// FIND IF GEOMETRY IS THERE
			if (lineT.find(regSolver::Geometry[i])<std::string::npos){
				gotMesh=true;
				break;
			}
			else{
				// COPY LINE TO TEMP
				fputs(lineT.c_str(),temp);
			}
		} // while

		// WRITE INFO CORRESPONDING TO GEOMETRY INTO TEMP
		fputs(regSolver::Geometry[i].c_str(),temp);
		fputs("\n",temp);
		
		// GET PROJECTION TYPE (CG or DG)
		val=regSolver::getMeshProperty(regSolver::Geometry[i],prop);

		if(gotMesh){
			// FIND OUT IF THERE IS INFO ON THIS PROJECTION
			gotProj=false;
			stop=false;
			while(!(stop||gotProj)){				
				if(!fgets(buffer,sizeof(buffer),mRes)){
					stop=true; continue;
				}
				lineT.assign(buffer,sizeof(buffer)-1);
				
				// Clear buffer
				regSolver::clearBuffer(buffer,sizeof(buffer));

				if(lineT.find("="+val)<std::string::npos){
					// There is info on this projection
					gotProj=true;
					continue;
				}
				if(lineT.find("-------------------------")<std::string::npos){
					// There's no info on this projection 
					stop=true;
					continue;
				}
				// IF STILL IN THE LOOP COPY LINE TO TEMP
				fputs(lineT.c_str(),temp);
			} // while
		};

		// ADD PROJECTION TYPE TO TEMP
		lineT=prop+"="+val;
		fputs(lineT.c_str(),temp);
		fputs("\n",temp);
	
		// RUN SOLVER WITH GEOMETRY
		command=regSolver::makeCommand(i);
		if(regSolver::RunProg(command)){
			// cout << "Failed to run command: \n\t"<<command<<endl;
			return 2;
		}
		
		// GET TIME STEPS 
		lineT=regSolver::getTline(regSolver::Output.c_str(),SOLVER_TIME_STR,false);
		fputs(lineT.c_str(),temp);

		// GET NUMBER OF MODES
		lineT=regSolver::getTline(regSolver::Output.c_str(),SOLVER_MODES_STR,false);
		fputs(lineT.c_str(),temp);

		// OPEN FILE CONTAINING OUTPUT OF PROG
		out = fopen(regSolver::Output.c_str(), "r");
		if (!out) {
			// std::cout << "Could not open output file (" << regSolver::Output << ")!" << std::endl;
			return 3;
		}

		// PROCESS THE OUTPUT OF PROG
		while(fgets(buffer,sizeof(buffer),out)){
			lineT.assign(buffer,sizeof(buffer)-1);

			// Clear buffer
			regSolver::clearBuffer(buffer,sizeof(buffer));
	
			// RETRIEVE THE LINES CONTAINING ERROR INFO
			errLine=(lineT.find("error")<std::string::npos);
			errLine=errLine||(lineT.find("Error")<std::string::npos);
			errLine=errLine||(lineT.find("ERROR")<std::string::npos);

			if (errLine){
				// ADD ALL THE LINES CONATINING ERRORS TO THE .OK FILE
				fputs(lineT.c_str(),temp);
			}
		} // while
		fclose(out);

		// ADD A LINE THAT SEPARATES INFO FOR DIFFERENT GEOMETRIES
		fputs("----------------------------------------\n",temp);
		
		newGeo=false;

		if(gotMesh){
			// FINISH COPYING FROM ORIGINAL TO TEMP
			while(fgets(buffer,sizeof(buffer),mRes)){
				lineT.assign(buffer,sizeof(buffer)-1);
				
				// Clear buffer
				regSolver::clearBuffer(buffer,sizeof(buffer));

				if(newGeo){
					// COPY TO TEMP
					fputs(lineT.c_str(),temp);
				}
				else{
					if (lineT.find(".xml")<std::string::npos){
						fputs(lineT.c_str(),temp);
						newGeo=true;
					}
				}
			} // while
		} // if gotmesh
		
		// Close Files
		fclose(mRes);
		fclose(temp);
		
		// Update Regress
		remove(fname.c_str());
		rename(tName.c_str(),fname.c_str());

	} // for	
	return 0;
};

int regSolver::readRegress(std::string mesh, std::string errStr, std::string *vars, std::string *errors, int &modes){
	std::string var_str="(variable ";
	FILE *mRes;
	char buffer[40];
	std::string line="",fname="";
	bool nextGeo=false, gotProj=false, gotMesh=false, gotErr=false;
	std::string err="", timestep="", var="", mm="";
	std::string prop="Projection", val="";
	modes=0;

	// REGRESS.OK FILE IS IN okPath
	fname=regSolver::okPath+regSolver::Regress;
	// vars and errors have to be initialized arrays
	std::string *varArray=vars;
	std::string *errArray=errors;
	int pos=0;

	//	OPEN .OK FILE
	mRes=fopen(fname.c_str(),"r");
	if(!mRes){
		// Unable to open file
		return 1;
	}

	// FIND MESH IN .OK FILE 
	while(fgets(buffer, sizeof(buffer), mRes)&&!gotMesh){
		line.assign(buffer,sizeof(buffer)-1);
		
		// Clear buffer
		regSolver::clearBuffer(buffer,sizeof(buffer));

		if (line.find(mesh)<std::string::npos){
			// CORRECT MESH FILE IN REGRESS
			gotMesh=true;
			
			// GET PROJECTION TYPE (CG or DG)
			val=regSolver::getMeshProperty(mesh,prop);
			gotProj=false;
			while(fgets(buffer, sizeof(buffer), mRes)&&!gotProj){
				line.assign(buffer,sizeof(buffer)-1);
				if(line.find(val)<std::string::npos){
					gotProj=true;
				}
			};
			
			if(!gotProj){
				// THERE IS NO INFO ON THIS PROJECTION FOR THE MESH
				return 3;
			}

			// EXTRACT INFO CORRESPONDING TO MESH
			while(fgets(buffer, sizeof(buffer), mRes)&&!nextGeo){
				line.assign(buffer,sizeof(buffer)-1);		
				// Clear buffer
				// regSolver::clearBuffer(buffer,sizeof(buffer));

				if(line.find(".xml")<std::string::npos){
					// NEXT MESH INFO
					nextGeo=true;
					continue;
				};
				
				// GET TIME STEP
				if(line.find(SOLVER_TIME_STR)<std::string::npos){
					getVal(line,SOLVER_TIME_STR,timestep);
					// cout << "Time: " << timestep <<endl;
				};

				// GET MODES
				if(line.find(SOLVER_MODES_STR)<std::string::npos){
					getVal(line,SOLVER_MODES_STR,mm);
					modes=atoi(mm.c_str());
					// cout<<"Modes: "<<mm<<modes;
				}

				// GET VARS AND ERRORS
				if(line.find(errStr)<std::string::npos){
					gotErr=regSolver::getVal(line,errStr, err);
					*errArray=err;
					++errArray;

					// GET VAR 
					pos=line.find(var_str)+var_str.length();
					var=line.substr(pos,1);
					*varArray=var;
					++varArray;

					// cout<<"Error: "<<var<<"->"<<err<<endl;
				}
			} // nextGeo
		}// if mesh
	}// while find mesh
	
	if(!gotMesh){
		// Failed to find mesh info in .OK file
		return 2;
	};

	fclose(mRes);
	
	return 0;
};


/*---------------------------------------------------------/
/ Process .ERR File
/---------------------------------------------------------*/
int regSolver::readOutput(std::string errStr, std::string *vars,std::string *errors){
	FILE *mRes;
	char buffer[50];
	bool gotErr=false,
		 stop=false;	
	// vars and errors have to be initialized arrays
	std::string *varArray=vars;
	std::string *errArray=errors;
	//// OUTPUT.ERR IS LOCATED IN WORKING DIR
	std::string fname=regSolver::Output;

	std::string line="",time="",project="",err="";
	std::string var_str="(variable ",var="";
	int pos=0;
	
	// OPEN FILE .ERR
	mRes=fopen(fname.c_str(),"r");
	if(!mRes){
		// Unable to open file
		return 1;
	}
	while(fgets(buffer, sizeof(buffer), mRes)){	
		line.assign(buffer,sizeof(buffer)-1);

		// GET TIME STEP
		if (line.find(SOLVER_TIME_STR)<std::string::npos){
			time=line;
		}

		// GET PROJECTION
		if (line.find("Projection=")<std::string::npos){
			project=line;
		}

		if (line.find(errStr)<std::string::npos){
			// GET ERROR AND CONVERT TO DOUBLE 
			gotErr=regSolver::getVal(line,errStr, err);
			*errArray=err;
			++errArray;

			// GET VAR 
			pos=line.find(var_str)+var_str.length();
			var=line.substr(pos,1);
			*varArray=var;
			++varArray;
		}		
	}
	fclose(mRes);
	return 0;
};

/*---------------------------------------------------------/
/ Create Command to be Run as required by solver
/---------------------------------------------------------*/
std::string regSolver::makeCommand(int i){
	std::string prog="",mesh="",comm="";
	prog=regBase::Path+regBase::Prog;
	mesh=regBase::Geometry[i];
	comm=prog+" "+mesh;
	comm+=" > "+regSolver::Output;
	return comm;
};
/*---------------------------------------------------------/
/ Create Solver and Mesh Paths
/---------------------------------------------------------*/
std::string regSolver::solverPath(void){
	std::string p="";
	p.assign(REG_PATH);
	p=p+"/../solvers/builds/";
	p=p+regSolver::getProg()+"/";
	if(USR_WIN_32){
		p=p+"Release/";
	}
	return p;
};
std::string regSolver::mPath(void){
	std::string p="";
	p.assign(REG_PATH);
	p=p+"/../solvers/"+regSolver::Prog+"/";
	return p;
};
