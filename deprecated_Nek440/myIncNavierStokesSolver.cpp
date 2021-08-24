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

#include <SolverUtils/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include "EquationSystems/VelocityCorrectionScheme.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "EquationSystems/CoupledLinearNS_TT.h"
#include "EquationSystems/CoupledLinearNS_trafoP.h"

//#include <MultiRegions/ExpList2D.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    string vDriverModule;
    DriverSharedPtr drv;

    /*cout << "argc " << argc << endl;
    cout << "*argv " << *argv << endl;
    cout << "*(argv[1]) " << *(argv[1]) << endl;
    cout << "(*argv)[1] " << (*argv)[1] << endl;
    std::string s(argv[1]);
    cout << "s " << s << endl;*/

    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        // Create driver
        drv = GetDriverFactory().CreateInstance(vDriverModule, session); 

	// check the Solver type
	cout << "SolverType " << session->GetSolverInfo("SolverType") << endl;

	if (session->GetSolverInfo("SolverType") == "CoupledLinearisedNS_TT")
	{
		CoupledLinearNS_TT CLNS(session);
		CLNS.offline_phase();
		CLNS.online_phase();
	
		cout<<"Online phase finished"<<endl;
	}
	else
	{
		drv->Execute();
	}

        session->Finalise();
    }
    catch (const std::runtime_error&)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }
    
    return 0;
}
