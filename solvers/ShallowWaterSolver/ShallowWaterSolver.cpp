///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterSolver.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Shallow Water Equations framework solver
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/EquationSystem.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    // Create session reader and MeshGraph.
    LibUtilities::SessionReaderSharedPtr session;
    session = LibUtilities::SessionReader::CreateInstance(argc, argv);

    SpatialDomains::MeshGraphSharedPtr graph;
    graph = SpatialDomains::MeshGraph::Read(session);

    time_t starttime, endtime;
    NekDouble CPUtime;
    EquationSystemSharedPtr equ;

    // Record start time.
    time(&starttime);

    // Create instance of module to solve the equation specified in the session.
    try
    {
        equ = GetEquationSystemFactory().CreateInstance(
            session->GetSolverInfo("EQTYPE"), session, graph);
    }
    catch (int e)
    {
        ASSERTL0(e == -1, "No such solver class defined.");
    }

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
    cout << "-------------------------------------------" << endl;
    cout << "Total Computation Time = " << CPUtime << " hr." << endl;

    
    for(int i = 0; i < equ->GetNvariables(); ++i)
    {
        cout << "L 2 error (variable " << equ->GetVariable(i)  << "): " << equ->L2Error(i,true) << endl;
        cout << "L inf error (variable " << equ->GetVariable(i)  << "): " << equ->LinfError(i) << endl;
    }

    session->Finalise();
}
