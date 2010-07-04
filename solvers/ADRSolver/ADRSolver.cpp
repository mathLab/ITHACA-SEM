///////////////////////////////////////////////////////////////////////////////
//
// File ADRSolver.cpp
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
// Description: Advection Diffusion Reaction solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath>

//#include <ADRSolver/ADR.h>
#include <ADRSolver/EquationSystem.h>
#include <ADRSolver/SessionReader.h>
#include <ADRSolver/EquationSystems/AlievPanfilov.h>
using namespace Nektar;


int main(int argc, char *argv[])
{

    ASSERTL0(argc == 2,"\n \t Usage: ADRSolver  sessionfile \n");

    string filename(argv[1]);
    time_t starttime, endtime;
    NekDouble CPUtime;

    SessionReaderSharedPtr session;
    EquationSystemSharedPtr equ;

    // Time integration function object for unsteady equations
    LibUtilities::TimeIntegrationSchemeOperators ode;

    time(&starttime);

    // Create session reader.

    session = MemoryManager<SessionReader>::AllocateSharedPtr(filename);

    // Create specific equation as specified in the session file.
    try
    {
        equ = EquationSystemFactory::CreateInstance(
                session->getSolverInfo("EQTYPE"), session, ode);
        equ->doInitialise();
    }
    catch (int e)
    {
        ASSERTL0(e == -1, "No such class defined.");
    }

    equ->printSummary(cout);

    // Solve steady-state equations or perform time integration.
    if (equ->isSteady())
    {
        equ->doSolveHelmholtz();
    }
    else
    {
        equ->SetInitialConditions();
        equ->GeneralTimeIntegration(equ->GetSteps(), equ->getTimeIntMethod(), ode);
    }

    time(&endtime);
    CPUtime = (1.0/60.0/60.0)*difftime(endtime,starttime);

    // Write  output to .fld file
    equ->Output();

    cout << "-------------------------------------------" << endl;
    cout << "Total Computation Time = " << CPUtime << " hr." << endl;

    // Evaluate L2 Error
    for(int i = 0; i < equ->GetNvariables(); ++i)
    {
        cout << "Variable " <<  equ->GetVariable(i)  << ": L 2 error = " << equ->L2Error(i) << endl;
        cout << "Variable " <<  equ->GetVariable(i)  << ": L infinity error = " << equ->LinfError(i) << endl;
    }
}

/**
* $Log $
**/
