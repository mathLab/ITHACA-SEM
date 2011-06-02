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
// Description: Advection Diffusion Reaction framework solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <ADRSolver/EquationSystem.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
using namespace Nektar;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout << "\nUsage: ADRSolver  sessionfile" << endl;
        GetEquationSystemFactory().PrintAvailableClasses();
        exit(1);
    }

    string filename(argv[1]);
    string vCommModule("Serial");
    time_t starttime, endtime;
    NekDouble CPUtime;

    LibUtilities::CommSharedPtr vComm;
    LibUtilities::SessionReaderSharedPtr session;
    EquationSystemSharedPtr equ;

    // Record start time.
    time(&starttime);

    try
    {
        // Create session reader.
        session = MemoryManager<LibUtilities::SessionReader>::AllocateSharedPtr(filename);

        // Create communicator
        if (session->DefinesSolverInfo("Communication"))
        {
            vCommModule = session->GetSolverInfo("Communication");
        }
        else if (LibUtilities::GetCommFactory().ModuleExists("ParallelMPI"))
        {
            vCommModule = "ParallelMPI";
        }
        vComm = LibUtilities::GetCommFactory().CreateInstance(vCommModule, argc, argv);

        // Create instance of module to solve the equation specified in the session.
        equ = GetEquationSystemFactory().CreateInstance(
                                    session->GetSolverInfo("EQTYPE"), vComm, session);

        // Print a summary of solver and problem parameters and initialise the
        // solver.
        equ->PrintSummary(cout);
        equ->DoInitialise();

        // Solve the problem.
        equ->DoSolve();

        // Record end time.
        time(&endtime);
        CPUtime = (1.0/60.0/60.0)*difftime(endtime,starttime);

        // Write output to .fld file
        equ->Output();

        // Evaluate and output computation time and solution accuracy.
        // The specific format of the error output is essential for the
        // regression tests to work.
        // Evaluate L2 Error
        for(int i = 0; i < equ->GetNvariables(); ++i)
        {
            NekDouble vL2Error = equ->L2Error(i,false);
            NekDouble vLinfError = equ->LinfError(i);
            if (vComm->GetRank() == 0)
            {
                cout << "L 2 error (variable " << equ->GetVariable(i) << ") : " << vL2Error << endl;
                cout << "L inf error (variable " << equ->GetVariable(i) << ") : " << vLinfError << endl;
            }
        }

        vComm->Finalise();
    }
    catch (const std::runtime_error& e)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }

    return 0;
}
