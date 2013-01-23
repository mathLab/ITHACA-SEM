///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSolver.cpp
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
// Description: Compressible Flow Equations framework solver
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/EquationSystem.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    // Create session reader
    LibUtilities::SessionReaderSharedPtr session;
    session = LibUtilities::SessionReader::CreateInstance(argc, argv);

    time_t starttime, endtime;
    NekDouble CPUtime;

    EquationSystemSharedPtr equ;

    // Record start time
    time(&starttime);

    // Create instance of module to solve the equation specified in the session
    try
    {
        equ = GetEquationSystemFactory().CreateInstance(
                                session->GetSolverInfo("EQTYPE"), session);
    }
    catch (int e)
    {
        ASSERTL0(e == -1, "No such solver class defined.");
    }

    // Print a summary of solver/problem parameters
    equ->PrintSummary(cout);
    
    // Initialise the problem
    equ->DoInitialise();

    // Solve the problem
    equ->DoSolve();

    // Record end time
    time(&endtime);
    
    // Compute the computational time in hours
    CPUtime = difftime(endtime, starttime);

    // Write output to .fld file
    equ->Output();

    // Evaluate and output computation time and solution accuracy
    if (session->GetComm()->GetRank() == 0)
    {
        cout << "-------------------------------------------" << endl;
        cout << "Total Computation Time = " << CPUtime << "s" << endl;
        cout << "-------------------------------------------" << endl;
    }
    
    for(int i = 0; i < equ->GetNvariables(); ++i)
    {
        // Get exact solution
        Array<OneD, NekDouble> exactsoln(equ->GetTotPoints(), 0.0);
        equ->EvaluateExactSolution(i, exactsoln, equ->GetFinalTime());
        
        NekDouble l2 = equ->L2Error  (i, exactsoln);
        NekDouble li = equ->LinfError(i, exactsoln);
        
        if (session->GetComm()->GetRank() == 0)
        {
            cout << "L 2 error (variable " << equ->GetVariable(i)  << "): " 
                 << l2 << endl;
            cout << "L inf error (variable " << equ->GetVariable(i)  << "): " 
                 << li << endl;
        }
    }
    
    session->Finalise();
}
