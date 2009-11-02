///////////////////////////////////////////////////////////////////////////////
//
// File regDemo.h
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
// Description: Regression test for Helmholtz 2D Demo.
//
///////////////////////////////////////////////////////////////////////////////

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
class regDemo: public regBase
{
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
    int readRegress(std::string mesh, std::string errStr, std::string &L2e,
                    int &modes);

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
