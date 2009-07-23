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

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <FitzHughNagumoSolver/FitzHughNagumo.h>


using namespace Nektar;


int main(int argc, char *argv[])
{
    
    ASSERTL0(argc == 2,"\n \t Usage: FitzHughNagumoSolver  meshfile \n");

    string fileNameString(argv[1]);
    time_t starttime, endtime;
    NekDouble CPUtime;
    
    time(&starttime);    

    //----------------------------------------------------------------
    // Read the mesh and construct container class
    FitzHughNagumo EAD(fileNameString);
    
    // Time integration function object for unsteady equations
    LibUtilities::TimeIntegrationSchemeOperators ode;

    int nsteps = EAD.GetSteps();
    NekDouble lambda = 0.0;

    EAD.Summary(cout);
    
    EAD.ZeroPhysFields(); // Zero phys field so that following switch is consistent
    
    // if not steady state equation set up initial conditions 
    EAD.SetInitialConditions();

    switch(EAD.GetEquationType())
    {
    case eFHNtest1:
        {
	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolve,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEeReactiontest1,EAD);	
          
	  // General Linear Time Integration
	  EAD.GeneralTimeIntegration(nsteps, EAD.GetTimeIntMethod(), ode);
        }

        break;

    case eFHNtest2:
        {
	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolve,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEeReactiontest2,EAD);	
          
	  // General Linear Time Integration
	  EAD.GeneralTimeIntegration(nsteps, EAD.GetTimeIntMethod(), ode);
        }

        break;

    case eFHNmonoplane:
        {
	  // Choose the method of deriving forcing functions
	  ode.DefineImplicitSolve (&FitzHughNagumo::ODEhelmSolve,EAD);	
	  ode.DefineOdeRhs        (&FitzHughNagumo::ODEeReactionmono,EAD);	
          
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
}

/**
* $Log $
**/
