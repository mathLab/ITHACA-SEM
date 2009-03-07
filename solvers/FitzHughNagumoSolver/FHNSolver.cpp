///////////////////////////////////////////////////////////////////////////////
//
// File FHNSolver.cpp
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
// Description: FitzHugh-Nagumo model solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <FitzHugh-Nagumo/FHN.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    
    ASSERTL0(argc == 2,"\n \t Usage: FHNSolver  meshfile \n");

    string fileNameString(argv[1]);
    
    //----------------------------------------------------------------
    // Read the mesh and construct container class
    FHN heart(fileNameString);
    
    int nsteps = heart.GetSteps();

    heart.Summary(cout);
    
    heart.ZeroPhysFields(); // Zero phys field so that following switch is consistent
    
    // Set up the intial conditions 
    heart.SetInitialConditions();
	
    // Create forcing function object
    LibUtilities::TimeIntegrationSchemeOperators ode;
	
    switch(heart.GetEquationType())
    {			
	case eTestmodel:
	 {
	     // Choose time integration method
	      LibUtilities::TimeIntegrationMethod IntMethod = LibUtilities::eIMEXdirk_3_4_3;	
		
	     // Choose the method of deriving forcing functions
	      ode.DefineImplicitSolve       (&FHN::ODEhelmSolve,heart);	
	      ode.DefineOdeRhs              (&FHN::ODETest_Reaction,heart);	
				
	     // General Linear Time Integration
	      heart.GeneralTimeIntegration(nsteps, IntMethod, ode);
	  }
        break;

	case eFHN1961:
	 {
	     // Choose time integration method
	      LibUtilities::TimeIntegrationMethod IntMethod = LibUtilities::eIMEXdirk_3_4_3;	
		
	     // Choose the method of deriving forcing functions
	      ode.DefineImplicitSolve       (&FHN::ODEFHN_helmSolve,heart);	
	      ode.DefineOdeRhs              (&FHN::ODEFHN_Reaction,heart);	
				
	     // General Linear Time Integration
	      heart.GeneralTimeIntegration(nsteps, IntMethod, ode);
	  }
        break;

    default:
        ASSERTL0(false,"Unknown or undefined equation type");
    }
     // Dump output
    heart.Output();

    // Evaluate L2 Error
    for(int i = 0; i < heart.GetNvariables(); ++i)
    {
        cout << "L2 Error (variable "<< heart.GetVariable(i) <<"): " << heart.L2Error(i) << endl;
    }

}

/**
* $Log $
**/
