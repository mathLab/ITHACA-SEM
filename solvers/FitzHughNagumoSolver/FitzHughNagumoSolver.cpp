///////////////////////////////////////////////////////////////////////////////
//
// File FitzHughNagumoSolver.cpp
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

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <FitzHughNagumoSolver/FitzHughNagumo.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    
    ASSERTL0(argc == 2,"\n \t Usage: FitzHughNagumoSolver  meshfile \n");

    // Create session reader.
    LibUtilities::SessionReaderSharedPtr session;
    session = LibUtilities::SessionReader::CreateInstance(argc, argv);

    string filename(argv[1]);
    time_t starttime, endtime;
    NekDouble CPUtime;
    
    time(&starttime);    

    //----------------------------------------------------------------
    // Read the mesh and construct container class
    FitzHughNagumo EAD(session);
    
    // Time integration function object for unsteady equations
    LibUtilities::TimeIntegrationSchemeOperators ode;

    int nsteps = EAD.GetSteps();
    int initialwavetype = EAD.initialwavetype();;

    EAD.Summary(cout);
    
    EAD.ZeroPhysFields(); // Zero phys field so that following switch is consistent

    switch(EAD.GetEquationType())
    {
        // IMEX test for diffusion reaction test problem
    case eIMEXtest:
        {
            EAD.SetInitialConditions();

	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolvetest,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEeReactionIMEXtest,EAD);	
          
	  // General Linear Time Integration
	  EAD.GeneralTimeIntegration(nsteps, EAD.GetTimeIntMethod(), ode);
        }
        break;

        // du/dt = \varepsilon \nabla + (1/\varepsilon) u*( 1- u) * { 2*u + \varepsilon  }
        // the exact solution = 1.0 / ( 1 + exp((|x| - t - 1.0)/varepsilon) )
    case eFHNtest1:
        {
            EAD.SetInitialConditions();

	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolvetest,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEeReactiontest1,EAD);	
          
	  // General Linear Time Integration
	  EAD.GeneralTimeIntegration(nsteps, EAD.GetTimeIntMethod(), ode);
        }
        break;

        // du/dt = \varepsilon \nabla + (2/\varepsilon) ( 1- u) u*u
        // the exact solution = 1.0 / ( 1 + exp((x - t - 1.5)/varepsilon) )
    case eFHNtest2:
        {
            EAD.SetInitialConditions();
            
	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolvetest,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEeReactiontest2,EAD);	
          
	  // General Linear Time Integration
	  EAD.GeneralTimeIntegration(nsteps, EAD.GetTimeIntMethod(), ode);
        }
        break;

    case eFitzHughNagumo:
        {
            EAD.SetFHNInitialConditions(initialwavetype);

	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolvemono,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEFitzHughNagumo,EAD);	
          
	  // General Linear Time Integration
	  EAD.GeneralTimeIntegration(nsteps, EAD.GetTimeIntMethod(), ode);
        }
        break;

    case eFHNRogers:
        {
            EAD.SetFHNInitialConditions(initialwavetype);

	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolvemono,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEFHNRogers,EAD);	
          
	  // General Linear Time Integration
	  EAD.GeneralTimeIntegration(nsteps, EAD.GetTimeIntMethod(), ode);
        }
        break;

    case eFHNmonohetero:
        {
            EAD.SetFHNInitialConditions(initialwavetype);

            EAD.Setdiffusivity();

	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolvehetero,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEFHNRogers,EAD);	
          
	  // General Linear Time Integration
	  EAD.GeneralTimeIntegration(nsteps, EAD.GetTimeIntMethod(), ode);
        }
        break;

    case eNoEquationType:
    default:
        ASSERTL0(false,"Unknown or undefined equation type");
    }

    time(&endtime);
    CPUtime = (1.0/60.0/60.0)*difftime(endtime,starttime);

    // Write  output to .fld file
    EAD.Output();

    cout << "-------------------------------------------" << endl;
    cout << "Total Computation Time = " << CPUtime << " hr." << endl;

    // Evaluate L2 Error
    for(int i = 0; i < EAD.GetNvariables(); ++i)
    {
        cout << "L2 Error (variable "<< EAD.GetVariable(i) <<"): " << EAD.L2Error(i) << endl;
    }

    session->Finalise();
}

/**
* $Log $
**/
