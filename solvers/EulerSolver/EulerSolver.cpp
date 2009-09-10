///////////////////////////////////////////////////////////////////////////////
//
// File EulerSolver.cpp
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
// Description: Euler solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <EulerSolver/EulerEquations.h>


using namespace Nektar;


int main(int argc, char *argv[])
{
    
    ASSERTL0(argc == 2,"\n \t Usage: EulerSolver  meshfile \n");

    string fileNameString(argv[1]);
    
    //----------------------------------------
    // Read the mesh and construct container class
    
    EulerEquations dom(fileNameString);
    //----------------------------------------

    
    //----------------------------------------
    // print session info 
    
    dom.Summary(cout);
    //----------------------------------------

    
    //----------------------------------------   
    dom.ZeroPhysFields(); 
       
    // Set up the intial conditions 
    switch(dom.m_problemType)
      {
      case eGeneral:
	{
	  dom.SetInitialConditions();
	}
	break;
      case eIsentropicVortex:
	{
	  dom.SetIsenTropicVortex();
	}
	break;
      case eSubsonicCylinder:
	{
	  dom.SetInitialConditions();
	}
	break;
      case eRinglebFlow:
	{
	  dom.SetInitialRinglebFlow();
	}
	break;
      }

    //----------------------------------------

    
    //----------------------------------------
    // Integrate from start time to end time
    
    int nsteps = dom.GetSteps();

    // Create forcing function object
    LibUtilities::TimeIntegrationSchemeOperators ode;
    
    // Choose time integration method
    LibUtilities::TimeIntegrationMethod IntMethod = LibUtilities::eClassicalRungeKutta4;		
    
    // Choose the method of deriving forcing functions
    ode.DefineOdeRhs       (&EulerEquations::ODErhs,&dom);		
		
    // General Linear Time Integration
    dom.GeneralTimeIntegration(nsteps, IntMethod, ode);
    //----------------------------------------


    //----------------------------------------
    // Dump output
     
    dom.Output();
    //----------------------------------------

    
    //----------------------------------------
    // print error
    
    // Evaluate L2 Error   
    Array<OneD, NekDouble> exactsolution(dom.GetTotPoints(),0.0);
    switch(dom.m_problemType)
      {
      case eGeneral:
	{
	  dom.EvaluateExactSolution(0,exactsolution,0,0);
	}
	break;
      case eIsentropicVortex:
	{
	  dom.GetExactIsenTropicVortex(exactsolution, 0);
	}
	break;
      case eSubsonicCylinder:
	{
	  dom.EvaluateExactSolution(0,exactsolution,0,0);
	}
	break;
      case eRinglebFlow:
	{
	  dom.GetExactRinglebFlow(exactsolution, 0);
	}
	break;
      }	

    cout << "L2 Error (variable "<< dom.GetVariable(0) <<"): " << dom.L2Error(0,0,exactsolution) << endl;
    cout << "L2 Error (variable "<< dom.GetVariable(0) <<"): " << dom.L2Error(0) << endl;

    dom.GetExactRinglebFlow(exactsolution, 4);

    //---------------------------------------
    
   
}

/**
* $Log $
**/
