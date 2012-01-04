///////////////////////////////////////////////////////////////////////////////
//
// File RegressBase.cpp
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
// Description: Base class for regression tests 
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <sys/stat.h>
#include "RegressBase.h"

#include <boost/lexical_cast.hpp>
#include <boost/version.hpp>

RegressBase::RegressBase(){
    m_prog="";
    m_progPath="";
    m_input  ="";
    m_okFile=m_prog+".ok";
    m_output=m_prog+".err";
}


RegressBase::RegressBase(std::string relative_demo_path, std::string demo_name, std::string input, std::string relative_okfile_path, unsigned int np)
{
    m_prog     = demo_name;
    m_okFile   = demo_name+".ok";
    m_output   = demo_name+".err";
    m_input    = input;
    m_np       = np;

    // DEMO and .OK PATH
    //m_progPath.assign(REG_PATH);
    m_progPath+=relative_demo_path;
    if(USR_WIN_32)
    {
        //m_progPath=m_progPath+"Release/";
    }
    m_okPath.assign(REG_PATH);
    m_okPath+=relative_okfile_path;
}

void RegressBase::Print(){
    std::cout << "===============================================\n";
    std::cout << "Program: " << m_prog <<std::endl;
    std::cout << "Path:    " << m_progPath <<std::endl;
    std::cout << "Output:  " << m_output <<std::endl;
    std::cout << "OKFile:  " << m_okFile <<std::endl;
    std::cout << "OkPath:  " << m_okPath <<std::endl;
    std::cout << "===============================================\n";
};

/*---------------------------------------------------------/
/ Set Attributes  routines
/---------------------------------------------------------*/
int RegressBase::SetProg(std::string D)
{
    m_prog=D;
    return 0;
};

int RegressBase::SetProgPath(std::string P){
    m_progPath=P;
    return 0;
};

int RegressBase::SetOutput(std::string out)
{
    m_output=out;
    return 0;
};

int RegressBase::SetOkFile(std::string reg)
{
    m_okFile=reg;
    return 0;
};

int RegressBase::SetOkPath(std::string op)
{
    m_okPath=op;
    return 0;
};


/*---------------------------------------------------------/
  / Get Attributes
  /---------------------------------------------------------*/
inline std::string RegressBase::GetProg(void)
{
    return m_prog;
};

inline std::string RegressBase::GetProgPath(void)
{
    return m_progPath;
};


inline std::string RegressBase::GetOutput(void)
{
    return m_output;
};

inline std::string RegressBase::GetOkFile(void)
{
    return m_okFile;
};

inline std::string RegressBase::GetOkPath(void)
{
    return m_okPath;
};

/*---------------------------------------------------------/
  / Test Error
  /---------------------------------------------------------*/
int RegressBase::Test(std::string errStr)
{
    int modes=0,fail=0;
    int i=0,n=0;
    std::string outErr, regErr, comm;
    
    // CHECKFILES 
    int chk = CheckFiles();

    switch (chk)
    {
    case 1:
    case 2: 
        {
            // Prog or baseMesh (if defined) not in place
            return chk;
        }
    case 4:
    case 7: 
        {
            // .OK FILE DOES NOT EXIST OR IS NOT IN PLACE
            return 3;
        }
    default:
        {
        }
    } 
    
    // READ  .OK FILE
    if(ReadOkFile(errStr,regErr))
    {
        fail=4;
    }

    // RUN SOLVER WITH GEOMETRY
    if(RunProg(MakeCommand()))
    {
        fail=6;
    };
    
    // READ OUTPUT .ERR FILE
    ReadOutput(errStr,outErr);
    
    // COMPARE ERRORS FROM .ERR AND .OK		
    if(CompareError(outErr,regErr))
    {
        std::cout << std::endl;
        std::cout << "Test result:     " << outErr << std::endl;
        std::cout << "Expected result: " << regErr << std::endl;
        fail=7;
    }
    
    return fail;
};

int RegressBase::TestL2(void)
{
    return Test(DEMO_L2ERROR_STR);
};

int RegressBase::TestLinf(void)
{
    return Test(DEMO_LINF_ERROR_STR);
};

int RegressBase::PrintTestError(int fail)
{
    switch(fail)
    {
    case 0: 
        break;
    case 1:
        {
            std::cout << "- The file "<< m_prog << " could not be found!"<<std::endl;
            std::cout << "- Suggestions:\n\tCheck Path:\n\t"<<m_progPath<<std::endl;
            break;
        };
    case 3:
        {
            std::cout << "- The file "<< m_okFile << " could not be found/opened!"<<std::endl;
            std::cout << "- Suggestions:\n\tCheck Path:\n\t"<<m_okPath<<std::endl;
            break;
        };
    case 4:
        {
            std::cout << "- The file "<< m_okFile << " was found but no values matched this set of parameters"<< std::endl;
            std::cout << "- Check a reference value has been generated for:\n\t"
                      << "Program    : " << m_prog << "\n\tParameters: " << m_input << std::endl;
        }
    case 5:
        {
            std::cout << "- There is no info for this mesh in .OK file!"<<std::endl;
            break;
        };
    case 6:
        {
            std::cout << "- The file "<< m_output << " could not be found!"<<std::endl;
            std::cout << "- Suggestions:\n\tCheck Path:\n\t"<<m_progPath<<std::endl;
            break;
        };
    case 7:
        {
            std::cout << "- The error obtained the demo test does not agree with the\n  value stored in the .OK file"<<std::endl;
            break;
        };
    };
    return 0;
};

/*---------------------------------------------------------/
  / Process .ok File
  /---------------------------------------------------------*/
int RegressBase::MakeOkFile()
{
    std::string command;
    std::string lineT = "",line = "",fname = "";
    std::string tName = m_okPath+"temp.ok";	
    FILE *OkFile,*out,*temp;	
    char buffer[BUFSIZ];
    bool errLine=false,gotInput=false,newNextInput=false;
    
    // CHECKFILES 
    int chk=CheckFiles();
    switch (chk)
    {
    case 1:case 2: 
        {
            // Prog or baseMesh (if defined) not in place
            return chk;
        }
    default:
        {
        }
    }
    
    // INFORM USER THAT MAKEREGRESS IS RUNNING
    std::cout << "\nGenerating Regression Test Results for\n\t" << m_prog <<" "<< m_input << std::endl;
    
    // CREATE REGRESS.OK (IF NEEDED)
    fname=m_okPath+m_okFile;
    if ((chk==4)||(chk==7))
    {
        // Regress is not in folder
        // Create Regress.OK and add solver name
        OkFile=fopen(fname.c_str(),"w+");
        if(!OkFile)
        {
            // Unable to open file
            return 3;
        }
        fputs(m_prog.c_str(),OkFile);
        fputs("\n==========================================\n",OkFile);
        fclose(OkFile);
    }
    
    // OPEN TEMPORARY .OK FILE
    temp=fopen(tName.c_str(),"w");
    
    // OPEN .OK FILE TO READ 
    OkFile=fopen(fname.c_str(),"r");
    
    if(!temp || !OkFile)
    {
        // Problems opening files
        return 3;
    }
    
    // Copy .ok file to temp until input is found
    gotInput=false;
    while(fgets(buffer,sizeof(buffer),OkFile) && !gotInput)
    {
        lineT.assign(buffer,sizeof(buffer)-1);
        
        // Clear buffer
        ClearBuffer(buffer,sizeof(buffer));
        
        // Retrieve line and find if input is there 
        if (lineT.find(m_input) < std::string::npos)
        {
            gotInput = true;
        }
        else
        {
            // COPY TO TEMP
            fputs(lineT.c_str(),temp);
        }
    } // while
    
    // Write input  into  temp file
    fputs(m_input.c_str(),temp);
    fputs("\n",temp);
    
    // RUN SOLVER 
    if(RunProg(MakeCommand()))
    {
        // Failed to run command
        return 2;
    }
    
    // OPEN FILE CONTAINING OUTPUT OF PROG
    out = fopen(m_output.c_str(), "r");
    if (!out) 
    {
        return 3;
    }
    
    // Process the output of prog
    while(fgets(buffer,sizeof(buffer),out))
    {
        line.assign(buffer,sizeof(buffer)-1);
        
        // Clear buffer
        ClearBuffer(buffer,sizeof(buffer));
        
        // Retrieve line containing Error
        errLine=(line.find("error")<std::string::npos);
        errLine=errLine||(line.find("Error")<std::string::npos);
        errLine=errLine||(line.find("ERROR")<std::string::npos);
        
        if (errLine)
        {
            fputs(line.c_str(),temp);
        }
    } // while
    
    fclose(out);
    
    // ADD A LINE THAT SEPARATES INFO FOR DIFFERENT GEOMETRIES
    fputs("----------------------------------------\n",temp);
    
    if(gotInput)
    {
        newNextInput=false;
        // FINISH COPYING FROM ORIGINAL TO TEMP
        while(fgets(buffer,sizeof(buffer),OkFile))
        {
            lineT.assign(buffer,sizeof(buffer)-1);
            
            // Clear buffer
            ClearBuffer(buffer,sizeof(buffer));
            
            if(newNextInput)
            {
                // COPY TO TEMP
                fputs(lineT.c_str(),temp);
            }
            else // Find separator and
            {
                if (lineT.find("------")<std::string::npos)
                {
                    ClearBuffer(buffer,sizeof(buffer));
                    newNextInput=true;
                }
            }
        } // while
    } // if gotmesh
    
    // Close Files
    fclose(OkFile);
    fclose(temp);
    
    // Update Regress
    remove(fname.c_str());
    rename(tName.c_str(),fname.c_str());

    return 0;
}

int RegressBase::ReadOkFile(std::string &errStr, std::string &L2e)
{
    FILE *OkFile;
    char buffer[BUFSIZ];
    std::string line = "", fname = "";
    bool gotErr = false, nextInfo = true, stop = false;
    std::string mm = "";
    L2e = "";
    
    fname = m_okPath+m_okFile;
    OkFile  = fopen(fname.c_str(),"r");
    
    if(!OkFile)
    {
        // Unable to open file
        return 1;
    }

    while(fgets(buffer, sizeof(buffer), OkFile)&&!stop)
    {
        line.assign(buffer,sizeof(buffer)-1);
        
        // Clear buffer
        ClearBuffer(buffer,sizeof(buffer));
	
        if (line.find(m_prog)<std::string::npos)
        {
            // Correct prog
            stop=true;
            
            // Extract info corresponding to input 
            while(fgets(buffer, sizeof(buffer), OkFile)&&nextInfo)
            {
                line.assign(buffer,sizeof(buffer)-1);		
                
                // Clear buffer
                ClearBuffer(buffer,sizeof(buffer));			
		
                if(line.find(m_input.c_str())<std::string::npos)
                {
                    // Next info 
                    nextInfo=false;
                
                    // Get Error
                    while(fgets(buffer, sizeof(buffer), OkFile))
                    {
                        line.assign(buffer,sizeof(buffer)-1);		

                        ClearBuffer(buffer,sizeof(buffer));
                        
                        if(line.find(errStr)<std::string::npos)
                        {
                            gotErr=!GetVal(line,errStr, L2e);
                        }
                        else if(line.find("-------")< std::string::npos)
                        {
                            break;
                        }
                    }
                };
            } 
        } // if 
    } // while

    fclose(OkFile);

    if(!gotErr)
    {
        return 1;
    }
    else 
    {
        return 0;
    }
};

// PRIVATE METHODS
/*---------------------------------------------------------/
  / Check Files and Run Program
  /---------------------------------------------------------*/
int RegressBase::CheckFiles()
{
    // CHECK EXISTENCE OF FILES	//
    int fail=0;
    struct stat buf;
    std::string fname;
    
    // Check Prog
#if defined(NDEBUG)
    fname=m_progPath+m_prog; 
#else
    fname=m_progPath+m_prog+"-g"; 
#endif
    if (USR_WIN_32)
    {
        fname=fname+".exe";
    };

    if (stat(fname.c_str(),&buf))
    {
        // Prog is not in folder
        std::cout << fname << " is not available!" << std::endl;
        std::cout << "Please check Path is correct: \n\t"<<m_progPath<<std::endl;
        return 1;
    }

    // Check Regress
    fname=m_okPath+m_okFile;
    if (stat(fname.c_str(),&buf))
    {
        // Regress is not in folder
        fail+=4;
    }

    return fail;
};		

std::string RegressBase::MakeCommand(void)
{
    std::string prog="",mesh="",comm="";
#ifdef NEKTAR_USE_MPI
    if (m_np > 1)
    {
        prog="mpirun -np " + boost::lexical_cast<std::string>(m_np) + " ";
    }
#endif
    prog+="\"" + m_progPath+m_prog;
#if defined(NDEBUG) 
    comm=prog+"\" "+m_input;
#else
    comm=prog+"-g\" "+m_input;
#endif
    comm+=" > "+m_output;
    return comm;
}

int RegressBase::RunProg(std::string comm)
{			
    // Run Prog with Geometry files, storing output in Output file.	//
    int status;
    status=system(comm.c_str());
    if (status) 
    {
        std::cout << "Unable to run command " << comm << std::endl;
        return 1;
    }
    return 0;
};

/*---------------------------------------------------------/
  / Auxiliary methods to process errors
  /---------------------------------------------------------*/
std::string RegressBase::GetTline(std::string source, std::string T, bool usePath)
{
    // Extract the line containing T Error from Output of Prog. //	
    int lineNo=0;
    FILE *res;
    char buffer[40];
    std::string line="",fname=source;
    bool gotit=false;
    //Open file Source
    if(usePath){fname=m_progPath+source;}
    res = fopen(fname.c_str(), "r");
    if (!res) 
    {
        std::cout << "Could not open output file (" << source << ")!" << std::endl;
        return line;
    }
    // Process the output of Prog.
    while(fgets(buffer,sizeof(buffer),res) && !gotit)
    {
        ++lineNo;
        line.assign(buffer,sizeof(buffer)-1);
	
        // Clear buffer
        ClearBuffer(buffer,sizeof(buffer));
        
        // Retrieve line containing L2 Error
        if (line.find(T)<std::string::npos)
        {
            gotit=true;
        }
    } // while
    fclose(res);
    if(!gotit){line="";}
    return line;
};

std::string RegressBase::GetError(std::string code)
{
    int pos=15; //initial position
    std::string err="";
    bool stop=false,got1=false;
    if (code.find("L 2 error")<std::string::npos)
    {
        while(!stop)
        {
            if(IsDoubleChar(code.substr(pos,1)))
            {
                err.append(code.substr(pos,1));
                got1=true;
            }
            else{
                if (got1)
                {
                    stop=true;
                }
            }
            ++pos;
        }
    }
    return err;
};

int RegressBase::IsDoubleChar(std::string c)
{
    bool isd=false;	
    isd=(c=="0")||(c=="1")||(c=="2")||(c=="3")||(c=="4")||(c=="5");
    isd=isd||(c=="6")||(c=="7")||(c=="8")||(c=="9");
    if (isd){
        return 1;
    }
    isd=isd||(c==".")||(c=="e")||(c=="E")||(c=="-");
    if (isd)
    {
        return 2;
    }
    else
    {
        return 0;
    };
};

int RegressBase::GetVal(std::string l, std::string code, std::string &val)
{
    unsigned int pos=0; //initial position
    std::string line="";
    bool stop=false, got1=false;
    val="";
    if (l.find(code)<std::string::npos)
    {
        pos=l.find(":");
        line=l.substr(pos,std::string::npos);
        pos=0;
        while(!stop && !(pos>line.length()))
        {
            if(IsDoubleChar(line.substr(pos,1)))
            {
                val.append(line.substr(pos,1));
                if ((line.substr(pos,1)!=".")&&(line.substr(pos,1)!="e"))
                {
                    got1=true;
                }
            }
            else
            {
                if (got1)
                {
                    stop=true;
                }
            }
            ++pos;
        }
    }
    if (val==""){return 1;}

    return 0;
}
/*---------------------------------------------------------/
/ Compare Errors
/---------------------------------------------------------*/
int RegressBase::CompareError(std::string e1, std::string e2)
{
    double err1, err2, diff;

    // Convert error string to double and compare if different is less
    // that REGRESS_DOUBLE_TOL

    err1 = std::strtod(e1.c_str(),NULL);
    err2 = std::strtod(e2.c_str(),NULL);

    diff = err1-err2;
    diff = (diff < 0)? -diff:diff;

    if(diff < REGRESS_DOUBLE_TOL)
    {
        return 0;
    }
    return 1;
};


/*---------------------------------------------------------/
  / PROCESS .ERR FILE
  /---------------------------------------------------------*/
/** \brief Process the .err file 
 *  
 */ 
int RegressBase::ReadOutput(std::string &errStr, std::string &err)
{
    FILE *OkFile;
    char buffer[50];
    int gotErr=0;
    bool stop=false;
    //// OUTPUT.ERR IS LOCATED IN WORKING DIR
    std::string fname=m_output;
    std::string line="";
    int pos=0;
    err="";
    
    // OPEN FILE .ERR
    OkFile=fopen(fname.c_str(),"r");
    if(!OkFile)
    {
        // Unable to open file
        return 1;
    }

    while(fgets(buffer, sizeof(buffer), OkFile))
    {
        line.assign(buffer,sizeof(buffer)-1);
        
        if (line.find(errStr)<std::string::npos)
        {
            // GET ERROR
            gotErr=GetVal(line,errStr, err);
        }		
    }
    fclose(OkFile);
    return 0;	
};


void RegressBase::ClearBuffer(char* buf, int bSize)
{
    int i=0; 
    char empty=0;
    char *bufArray=buf;
    for(int i=0;i<bSize;i++)
    {
        *bufArray=empty;
        ++bufArray;
    }
}

std::string PortablePath(const boost::filesystem::path& path)
{
    boost::filesystem::path temp = path;
    #if BOOST_VERSION > 104200
    temp.make_preferred();
    return temp.string();
    #else
    return temp.file_string();
    #endif
    
}
