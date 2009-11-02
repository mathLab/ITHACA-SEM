///////////////////////////////////////////////////////////////////////////////
//
// File regDemo.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Regression test for Helmholtz 2D demo
//
///////////////////////////////////////////////////////////////////////////////

#include "regDemo.h"
using namespace std;

/*-----------------------------------------------------------
/ CLASS CONSTRUCTORS
/-----------------------------------------------------------*/
regDemo::regDemo(){
    // Default Constructur uses Helmholtz2D Demo
#ifdef NEKTAR_DEBUG
    Prog="Helmholtz2D-g";
#else
    Prog="Helmholtz2D";
#endif
    numTest=2;
    baseMesh="helmholtz2D.xml";
    Geometry[0]="helmholtz2D_low.xml"; //  LOW ORDER EXPANSION - 4 MODES
    Geometry[1]="helmholtz2D_high.xml"; //  HIGH ORDER EXPANSION - 10 MODES
    Regress=Prog+".ok";
    Output=Prog+".err";
    // DEMO and .OK PATH
    Path=regDemo::demoPath();
    okPath.assign(REG_PATH);
    okPath+="/Demo/";
    meshPath=regDemo::mPath();
}


regDemo::regDemo(std::string D){
    Prog=D;
    numTest=2;
    baseMesh=D+".xml";
    Geometry[0]=D+"_low.xml";
    Geometry[1]=D+"_high.xml";
    Regress=D+".ok";
    Output=D+".err";
    // DEMO and .OK PATH
    Path=regDemo::demoPath();
    okPath.assign(REG_PATH);
    okPath+="/Demo/";
    meshPath=regDemo::mPath();
}

regDemo::regDemo(std::string D, std::string bM){
    Prog=D;
    baseMesh=bM;
    numTest=2;
    Geometry[0]=D+"_low.xml";
    Geometry[1]=D+"_high.xml";
    Regress=D+".ok";
    Output=D+".err";
    // DEMO and .OK PATH
    Path=regDemo::demoPath();
    okPath.assign(REG_PATH);
    okPath+="/Demo/";
    meshPath=regDemo::mPath();
}

/*---------------------------------------------------------/
/  Test Errors
/---------------------------------------------------------*/
int regDemo::Test(std::string errStr){
    int modes=0,fail=0;
    int i=0,n=0;
    std::string outErr="", regErr="", comm="";

    // CHECKFILES
    int chk=regDemo::checkFiles();
    switch (chk){
        case 1:case 2: {
            // Prog or baseMesh (if defined) not in place
            return chk;
        }
        case 4:case 7: {
            // .OK FILE DOES NOT EXIST OR IS NOT IN PLACE
            return 3;
        }
        default:{
            if(regDemo::baseMesh!=""){
                // CREATE GEOMETRY FILES
                regDemo::createGeometry();
            }
            // 0,3 -> No problem, create or overwrite Geometry
        }
    } // switch

    for(i=0;i<regDemo::numTest;++i){
        // READ REGRESS .OK FILE
        if(regDemo::readRegress(regDemo::Geometry[i],errStr,regErr,modes)){
            fail=4;
            continue;
        }
        // CHECK MODES
        if(modes!=regDemo::Modes[i]){
            cout << regDemo::Geometry[i] << " -> " << modes << "::"
                 << regDemo::Modes[i] << endl;
            fail=5;
            continue;
        }

        // RUN SOLVER WITH GEOMETRY
        if(regDemo::RunProg(regDemo::makeCommand(i))){
            fail=6;
            continue;
        };

        // READ OUTPUT .ERR FILE
        regDemo::readOutput(errStr,outErr);
        // COMPARE ERRORS FROM .ERR AND .OK
        if(regDemo::compareError(outErr,regErr)){
            fail=7;
            cout << regDemo::baseMesh<<" (modes="<<regDemo::Modes[i];
            cout << ") errors: " << outErr << " :: " << regErr << endl;
        };
    }// for Geo

    return fail;
};

int regDemo::TestL2(void){
    return regDemo::Test(DEMO_L2ERROR_STR);
};

int regDemo::TestLInf(void){
    return regDemo::Test(DEMO_LINF_ERROR_STR);
};
/*---------------------------------------------------------/
/  PROCESS .OK FILE
/---------------------------------------------------------*/
int regDemo::makeRegress(){
    std::string command;
    std::string lineT="",line="",fname="";
    std::string tName=regDemo::okPath+"temp.ok";
    FILE *mRes,*out,*temp;
    char buffer[40];
    bool errLine=false,gotMesh=false,gotProj=false,newGeo=false;

    // CHECKFILES
    int chk=regDemo::checkFiles();
    switch (chk){
        case 1:case 2: {
            // Prog or baseMesh (if defined) not in place
            return chk;
        }
        default:{
            if(regDemo::baseMesh!=""){
                // CREATE GEOMETRY FILES
                regDemo::createGeometry();
            }
            // 0,4 -> No problem, overwrite Geometry
            // 3,7 -> Geometry files don't exist
        }
    }

    // INFORM USER THAT MAKEREGRESS IS RUNNING
    std::cout << "\nGenerating Regression Test Resutls for\n\t"
              << regDemo::Prog << std::endl;
    for(int i=0;i<regDemo::numTest;++i){
        std::cout << "\t" << regDemo::Geometry[i];
        std::cout <<" (Modes=" << regBase::Modes[i] << ")" << std::endl;
    }
    std::cout<< "Please Wait..."<<std::endl;

    // CREATE REGRESS.OK (IF NEEDED)
    fname=regDemo::okPath+regDemo::Regress;
    if ((chk==4)||(chk==7)){
        // Regress is not in folder
        // Create Regress.OK and add solver name
        mRes=fopen(fname.c_str(),"w+");
        if(!mRes){
            // Unable to open file
            return 3;
        }
        fputs(regDemo::Prog.c_str(),mRes);
        fputs("\n==========================================\n",mRes);
        fclose(mRes);
    }

    // FOR EACH GEOMETRY EXTRACT AND UPDATE OR INSERT INFO INTO OK FILE
    for (int i=0;i<regDemo::numTest;++i){

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
        while(fgets(buffer,sizeof(buffer),mRes) && !gotMesh){
            lineT.assign(buffer,sizeof(buffer)-1);

            // Clear buffer
            regDemo::clearBuffer(buffer,sizeof(buffer));

            // RETRIEVE LINE AND FIND IF GEOMETRY IS THERE
            if (lineT.find(regDemo::Geometry[i])<std::string::npos){
                gotMesh=true;
            }
            else{
                // COPY TO TEMP
                fputs(lineT.c_str(),temp);
            }
        } // while

        // WRITE GEOMETRY NAME INTO TEMP
        fputs(regDemo::Geometry[i].c_str(),temp);
        fputs("\n",temp);

        // RUN SOLVER WITH GEOMETRY
        command=regDemo::makeCommand(i);
        if(regDemo::RunProg(command)){
            // cout << "Failed to run command: \n\t"<<command<<endl;
            return 2;
        }

        // GET NUMBER OF MODES
        lineT=regDemo::getTline(regDemo::Output.c_str(),DEMO_MODES_STR,false);
        fputs(lineT.c_str(),temp);

        // OPEN FILE CONTAINING OUTPUT OF PROG
        out = fopen(regDemo::Output.c_str(), "r");
        if (!out) {
            // std::cout << "Could not open output file ("
            //           << regDemo::Output << ")!" << std::endl;
            return 3;
        }

        // PROCESS THE OUTPUT OF PROG
        while(fgets(buffer,sizeof(buffer),out)){
            // n gives number of variables {u,v,w,p}
            line.assign(buffer,sizeof(buffer)-1);

            // Clear buffer
            regDemo::clearBuffer(buffer,sizeof(buffer));

            // Retrieve line containing Error
            errLine=(line.find("error")<std::string::npos);
            errLine=errLine||(line.find("Error")<std::string::npos);
            errLine=errLine||(line.find("ERROR")<std::string::npos);

            if (errLine){
                fputs(line.c_str(),temp);
            }
        } // while

        fclose(out);

        // ADD A LINE THAT SEPARATES INFO FOR DIFFERENT GEOMETRIES
        fputs("----------------------------------------\n",temp);

        if(gotMesh){
            newGeo=false;
            // FINISH COPYING FROM ORIGINAL TO TEMP
            while(fgets(buffer,sizeof(buffer),mRes)){
                lineT.assign(buffer,sizeof(buffer)-1);

                // Clear buffer
                regDemo::clearBuffer(buffer,sizeof(buffer));

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
}

int regDemo::readRegress(std::string mesh, std::string errStr,
                         std::string &L2e, int &modes){
    FILE *mRes;
    char buffer[50];
    std::string line="",fname="";
    bool gotErr=false,
        nextGeo=false,
        stop=false;
    std::string mm="";
    L2e="";
    modes=0;

    fname=regDemo::okPath+regDemo::Regress;
    mRes=fopen(fname.c_str(),"r");
    if(!mRes){
        // Unable to open file
        return 1;
    }
    while(fgets(buffer, sizeof(buffer), mRes)&&!stop){
        line.assign(buffer,sizeof(buffer)-1);

        // Clear buffer
        regDemo::clearBuffer(buffer,sizeof(buffer));

        if (line.find(mesh)<std::string::npos){
            // CORRECT MESH FILE IN REGRESS
            stop=true;

            // EXTRACT INFO CORRESPONDING TO MESH
            while(fgets(buffer, sizeof(buffer), mRes)&&!nextGeo){
                line.assign(buffer,sizeof(buffer)-1);

                // Clear buffer
                regDemo::clearBuffer(buffer,sizeof(buffer));

                if(line.find(".xml")<std::string::npos){
                    // NEXT MESH INFO
                    nextGeo=true;
                    continue;
                };

                // GET MODES
                if(line.find(DEMO_MODES_STR)<std::string::npos){
                    getVal(line,DEMO_MODES_STR,mm);
                    modes=atoi(mm.c_str());
                }

                // GET ERROR
                if(line.find(errStr)<std::string::npos){
                    gotErr=getVal(line,errStr, L2e);
                }
            } // While in Geo.xml
        } // if
    } // while
    fclose(mRes);
    if(!gotErr){return 0;}
    else {return 1;}
};

/*---------------------------------------------------------/
/  PROCESS COMMAND TO RUN
/---------------------------------------------------------*/
std::string regDemo::makeCommand(int i){
    std::string prog="",mesh="",comm="";
    prog=regDemo::Path+regDemo::Prog;
    mesh=regDemo::Geometry[i];
    comm=prog+" "+mesh+" "+mesh;
    comm+=" > "+regDemo::Output;
    return comm;
};

/*---------------------------------------------------------/
/ PROCESS .ERR FILE
/---------------------------------------------------------*/
int regDemo::readOutput(std::string errStr, std::string &err){
    FILE *mRes;
    char buffer[50];
    bool gotErr=false,
         stop=false;
    //// OUTPUT.ERR IS LOCATED IN WORKING DIR
    std::string fname=regDemo::Output;
    std::string line="";
    int pos=0;
    err="";

    // OPEN FILE .ERR
    mRes=fopen(fname.c_str(),"r");
    if(!mRes){
        // Unable to open file
        return 1;
    }
    while(fgets(buffer, sizeof(buffer), mRes)){
        line.assign(buffer,sizeof(buffer)-1);

        if (line.find(errStr)<std::string::npos){
            // GET ERROR AND CONVERT TO DOUBLE
            gotErr=regDemo::getVal(line,errStr, err);
        }
    }
    fclose(mRes);
    return 0;
};

/*---------------------------------------------------------/
/ Create Demo Path
/---------------------------------------------------------*/
std::string regDemo::demoPath(void){
    std::string p="";
    p.assign(REG_PATH);
    p=p+"/../builds/Demos/MultiRegions/";
    if(USR_WIN_32){
        p=p+"Release/";
    }
    return p;
};

std::string regDemo::mPath(void){
    std::string p="";
    p.assign(REG_PATH);
    p=p+"/../library/Demos/MultiRegions/";
    return p;
};
