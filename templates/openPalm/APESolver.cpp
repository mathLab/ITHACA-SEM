//PALM_UNIT -name APESolver\
//          -functions {C APESolver}\
//          -parallel mpi\
//          -minproc 1\
//          -maxproc 100000\
//          -object_files {APESolver.o APESolver/EquationSystems/APE.o APESolver/RiemannSolvers/APESolver.o APESolver/RiemannSolvers/UpwindSolver.o APESolver/RiemannSolvers/LaxFriedrichsSolver.o OpenPalmExchange.o}\
//          -comment {APESolver openPALM unit}
//
//PALM_CWIPI_COUPLING -name PRECISE-Nektar
//
//PALM_CWIPI_OBJECT -name basefields\
//                  -coupling PRECISE-Nektar\
//                  -intent IN
//
///////////////////////////////////////////////////////////////////////////////
//
// File APESolver.cpp
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
// Description: APE Equations framework solver
//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>

#include <SolverUtils/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include "OpenPalmExchange.h"

#include <palmlibc.h>
#include <cwipi.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

extern "C" int APESolver()
{
    // redirect stdout to a file

    FILE *logFile;
    logFile = fopen ("nektar.log","wt");

    dup2(fileno(logFile), STDOUT_FILENO);
    dup2(fileno(logFile), STDERR_FILENO);

    // set casename
    int   argc = 1;
    char *argv[] = {
        "APESolverExt",
    };

    LibUtilities::SessionReaderSharedPtr session;
    string vDriverModule;
    DriverSharedPtr drv;

    try
    {
        ASSERTL0(LibUtilities::GetCommFactory().ModuleExists("OpenPalm"), "OpenPalm comm module not found");
        LibUtilities::CommSharedPtr comm = LibUtilities::GetCommFactory().CreateInstance("OpenPalm", argc, argv);

        std::vector<std::string> filenames;
        filenames.push_back("APE_3DPulseWall_FRDG_MODIFIED.xml");

        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv, filenames, comm);

        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session);
        // this initializes the APE system

        // Execute driver
        drv->Execute();
        // This runs UnsteadySystem::v_DoSolve() exchange should happen in APE:v_preIntegrate()

        // Finalise session
        session->Finalise();
    }
    catch (const std::runtime_error& e)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }

    fclose(logFile);

    return 0;
}
