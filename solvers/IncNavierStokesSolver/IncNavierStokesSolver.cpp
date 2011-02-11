///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokesSolver.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <IncNavierStokesSolver/IncNavierStokes.h>
#include <IncNavierStokesSolver/VelocityCorrectionScheme.h>


using namespace Nektar;

int main(int argc, char *argv[])
{
#ifdef TIMING
    timeval timer1, timer2;
    NekDouble time1, time2;
    NekDouble exeTime;
#endif
    
    if(argc != 2 && (argc != 4))
    {
        cerr << "\n \t Usage: IncNavierStokes  input.xml [GlobalOptimizationFile ElementalOptimizationFile] \n" << endl;
        exit(1);
    }

    string fileNameString(argv[1]);
    string globoptfile;

    try
    {
        //----------------------------------------------------------------
        // Read the mesh and construct container class
        if(argc == 2)
        {
            globoptfile = NekNullString;
        }
        else
        {
            string eloptfile  (argv[3]);
            NekOptimize::LoadElementalOptimizationParameters(eloptfile);

            globoptfile = argv[2];
        }

        VelocityCorrectionScheme dom(fileNameString,globoptfile);
        dom.Summary(cout);
        
        switch(dom.GetEquationType())
        {
        case eUnsteadyStokes:
        case eUnsteadyNavierStokes:
            {

                int nsteps = dom.GetSteps();

                // Set initial condition using time t=0
                dom.SetInitialConditions(0.0);

#ifdef TIMING
                dom.AdvanceInTime(1);
                gettimeofday(&timer1, NULL);
#endif
                // Integrate from start time to end time
                dom.AdvanceInTime(nsteps);
#ifdef TIMING
                gettimeofday(&timer2, NULL);
                time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
                time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
                exeTime = (time2-time1)/1000000.0;
                cout << "Execution Time: " << exeTime << endl;
#endif
                break;
            }
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }

        // Dump output
        dom.Output();

        // Evaluate L2 Error
        cout << endl;
        for(int i = 0; i < dom.GetNvariables(); ++i)
        {
            cout << "L 2 error (variable " << dom.GetVariable(i) << ") : " << dom.L2Error(i,true) << endl;
            cout << "L inf error (variable " << dom.GetVariable(i) << ") : " << dom.LinfError(i) << endl;
        }

        return 0;
    }
    catch (const std::runtime_error& e)
    {
        return 1;
    }
}

/**
* $Log $
**/
