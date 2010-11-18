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

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <ADRSolver/EquationSystem.h>
#include <ADRSolver/SessionReader.h>
using namespace Nektar;

int main(int argc, char *argv[])
{
  if(argc != 2)
    {
      cout << "\n \t Usage: ShallowWaterSolver  sessionfile \n" << endl;
      exit(1);
    }
  
  string filename(argv[1]);
  time_t starttime, endtime;
  NekDouble CPUtime;
  
  SessionReaderSharedPtr session;
  EquationSystemSharedPtr equ;
  
  // Record start time.
  time(&starttime);
  
  // Create session reader.
  session = MemoryManager<SessionReader>::AllocateSharedPtr(filename);
  
  // Create instance of module to solve the equation specified in the session.
  try
    {
      equ = EquationSystemFactory::CreateInstance(session->GetSolverInfo("EQTYPE"), session);
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
}


//     //----------------------------------------
//     // Read the mesh and construct container class
    
//     EulerEquations dom(fileNameString);
//     //----------------------------------------

    
//     //----------------------------------------
//     // print session info 
    
//     dom.Summary(cout);
//     //----------------------------------------

    
//     //----------------------------------------   
//     dom.ZeroPhysFields(); 
       
//     // Set up the initial conditions 
//     switch(dom.m_problemType)
//       {
//       case eGeneral:
// 	{
// 	  dom.SetInitialConditions();
// 	}
// 	break;
//       case eIsentropicVortex:
// 	{
// 	  dom.SetIsenTropicVortex();
// 	}
// 	break;
//       case eSubsonicCylinder:
// 	{
// 	  dom.SetInitialConditions();
// 	}
// 	break;
//       case eRinglebFlow:
// 	{
// 	  dom.SetInitialRinglebFlow();
// 	}
// 	break;
//       }
//     //----------------------------------------

    
//     //----------------------------------------
//     // Integrate from start time to end time
    
//     int nsteps = dom.GetSteps();

//     // Create forcing function object
//     LibUtilities::TimeIntegrationSchemeOperators ode;
    
//     // Choose time integration method
//     LibUtilities::TimeIntegrationMethod IntMethod = LibUtilities::eClassicalRungeKutta4;		
    
//     // Choose the method of deriving forcing functions
//     ode.DefineOdeRhs       (&EulerEquations::ODErhs,&dom);		
		
//     // General Linear Time Integration
//     dom.GeneralTimeIntegration(nsteps, IntMethod, ode);
//     //----------------------------------------


//     //----------------------------------------
//     // Dump output
     
//     dom.Output();
//     //----------------------------------------

    
//     //----------------------------------------
//     // print error
    
//     // Evaluate L2 Error   
//     Array<OneD, NekDouble> exactsolution(dom.GetTotPoints(),0.0);
//     switch(dom.m_problemType)
//       {
//       case eGeneral:
// 	{
// 	  dom.EvaluateExactSolution(0,exactsolution,0);
// 	  // Evaluate L2 Error
// 	  cout << "Error:" << endl;
// 	  for(int i = 0; i < dom.GetNvariables(); ++i)
// 	    {
// 	      cout << "\t"<< dom.GetVariable(i) << ": "
// 		   << dom.LinfError(i) << " (Linf), "
// 		   << dom.L2Error(i) << " (L2) " << endl;
// 	    }
// 	}
// 	break;
//       case eIsentropicVortex:
// 	{
// 	  // Evaluate L2 Error
// 	  cout << "Error:" << endl;

// 	  // Evaluate L2 Error
// 	  for(int i = 0; i < dom.GetNvariables(); ++i)
// 	    {
// 	      dom.GetExactIsenTropicVortex(exactsolution, i);
// 	      cout << "L2 Error (variable "<< dom.GetVariable(i) <<"): " << dom.L2Error(i, exactsolution) << endl;
// 	    }
// 	}
// 	break;
//       case eSubsonicCylinder:
// 	{
// 	  dom.EvaluateExactSolution(0,exactsolution,0);
// 	}
// 	break;
//       case eRinglebFlow:
// 	{
// 	  dom.GetExactRinglebFlow(exactsolution, 0);
// 	}
// 	break;
//       }

//     // Writing exactsolution file
//     switch(dom.m_problemType)
//       {
//       case eRinglebFlow:
// 	{
// 	  dom.GetExactRinglebFlow(exactsolution, 4);
// 	}
// 	break;
//       }

//     //---------------------------------------
    
   
// }

// /**
// * $Log $
// **/
