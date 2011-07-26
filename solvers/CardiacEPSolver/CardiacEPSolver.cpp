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

#include <Auxiliary/Driver.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
using namespace Nektar;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout << "\nUsage: CardiacEPSolver  sessionfile" << endl;
        GetEquationSystemFactory().PrintAvailableClasses();
        exit(1);
    }

    string filename(argv[1]);
    string vCommModule("Serial");
    string vDriverModule("Standard");

    LibUtilities::CommSharedPtr vComm;
    LibUtilities::SessionReaderSharedPtr session;
    DriverSharedPtr drv;

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

        // Create driver
        if (session->DefinesSolverInfo("Driver"))
        {
            vDriverModule = session->GetSolverInfo("Driver");
        }
        drv = GetDriverFactory().CreateInstance(vDriverModule, vComm, session);

        // Execute driver
        drv->Execute();

        // Finalise communications
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
