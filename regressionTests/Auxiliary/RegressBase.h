///////////////////////////////////////////////////////////////////////////////
//
// File RegressBase.h
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
// Description: header file for Base class for regression tests 
//
///////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>

#ifndef RegressBase_H
#define RegressBase_H

#ifndef USR_WIN_32
#define USR_WIN_32 false
#endif //USR_WIN_32

#define DEMO_L2ERROR_STR "L 2 error"
#define DEMO_LINF_ERROR_STR "L infinity error"

#define SOLVER_L2ERROR_STR "L 2 error"

#define REGRESS_DOUBLE_TOL 1e-8

#include <boost/filesystem.hpp>

//header of class RegressBase
class RegressBase
{
protected:

    std::string m_prog;     /**< Name of Prog/Demo        */
    std::string m_input;    /**< input string to demo     */

    std::string m_progPath; /**< Path to demo             */
    std::string m_okPath;   /**< Path to ok file          */ 
    
    // Name of .OK File and Output of Prog .ERR File
    std::string m_okFile;   /**< name of .ok file         */ 
    std::string m_output;   /**< name of demo output file */ 
	
    unsigned int m_np;      /**< Number of processors     */

public:
    /*-----------------------------------------------------------
      / CLASS CONSTRUCTORS
      /-----------------------------------------------------------*/
    // default constructor
    RegressBase(void);
    
    RegressBase(std::string relative_demo_path, std::string demo_name, std::string input, std::string relative_okfile_path, unsigned int np = 1);
    
    /* \brief Print Attributes
     */
    void Print(void);
    
    /*-----------------------------------------------------------
                SET ATTRIBUTES
     -----------------------------------------------------------*/
    int SetProg(std::string);
    int SetOutput(std::string);
    int SetOkFile(std::string);
    int SetProgPath(std::string);
    int SetOkPath(std::string);
    
    /*-----------------------------------------------------------
                Get Attributes
      -----------------------------------------------------------*/
    std::string GetProg(void);
    std::string GetOutput(void);
    std::string GetOkFile(void);
    std::string GetProgPath(void);
    std::string GetOkPath(void);
    
    /** \brief compare errors to master errors 
     */
    int Test(std::string);
    
    /**  \brief Input: Error String (e.g. "L 2 error")
     *    Returns: 0, if Prog Error is equal to Error in ok file 
     */ 
    int TestL2(void);

    int TestLinf(void);
                
    /** \brief Prints to screen the suggested causes of the error
     *  reference by int as output from Test();
     */
    int PrintTestError(int fail);

    /*-----------------------------------------------------------
              PROCESS .OK FILE
      -----------------------------------------------------------*/

    /** \brief 
     *   Generate Master Errors for the Meshfiles in Geometry using 
     *   a known working version of the code.
     */
    int MakeOkFile(void);
                
    /** \brief If Data for Meshfile exists in .ok file, reads Error
     * (e.g. L2) and  stores its value in Errors, and Number of modes
     * in Modes  Returns: 0, if Error and Modes are retirived
     * successfully 1, otherwise.
     */ 
    int ReadOkFile(std::string &errStr, std::string &);
    
protected:
    std::string MakeCommand(void);
    
    /** \brief 
     *  \return  0, if all files needed are in the right place
     *           1, if Prog.exe is not found
     *           4, if .ok file not in place or not defined
     *           7, if .ok files missing
     */
    int CheckFiles(void);	
    
    
    /** \brief  Runs the Prog to be tested, takes command as input.
     */
    int RunProg(std::string comm);
    
    /*-----------------------------------------------------------
                EXTRACT FROM FILES
     -----------------------------------------------------------*/

    /** \brief
     * 
     * \param source : Name of File (string),
     * \param T is an identity string contained in line (e.g. "L 2 error"),
     * \param usePath is a bool which when true indicated to use path,
     * 
     * \return Full line containing input string.
     */ 
    std::string GetTline(std::string source, std::string T, bool usePath);

    /** \brief 
     *  \param code String that defines error type (e.g. "L 2 error")
     *  \return  Value of error as string (e.g. "0.1234e-004")
     */ 
    std::string GetError(std::string code);

    /** \brief 
     *
     * \param l is a string containing a line of a file,
     * \param code is a string with code to look for (e.g. "L 2 error"),
     * \param val is a string where the value is stored.
     *
     * \return 0 if Value is extracted successfully,  1 otherwise
     */ 
    int GetVal(std::string l, std::string code, std::string &val);
        
    /*-----------------------------------------------------------
                 PROCESS ERRORS
       -----------------------------------------------------------*/

    /** \brief 
     *
     * \param c: Any char
     *
     * \return 1 if char is a number 0-9, 
     *         2 if char is "e","E","." or "-",
     *         0 otherwise
     */ 
    int IsDoubleChar(std::string c);

    /** \brief 
     *
     * \param e1, \param e2 Two strings containing errors formated as
     * in getError
     *
     * \return  0 if they are equal,   1 otherwise
     */
    int CompareError(std::string e1, std::string e2);
    
    /* \brief  Read and process output file. 
     */
    int ReadOutput(std::string &errStr, std::string&);
    
    void ClearBuffer(char*,int);
};

std::string PortablePath(const boost::filesystem::path& path);

#endif // RegressBase_H
