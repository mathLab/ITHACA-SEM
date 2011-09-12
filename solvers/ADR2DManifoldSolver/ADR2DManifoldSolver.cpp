///////////////////////////////////////////////////////////////////////////////
//
// File ADR2DManifoldSolver.cpp
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
// Description: Advection Diffusion Reaction on 2D manifold solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <ADR2DManifoldSolver/ADR2DManifold.h>

using namespace Nektar;


int main(int argc, char *argv[])
{
    
    ASSERTL0(argc == 2,"\n \t Usage: ADR2DManifoldSolver  meshfile \n");

    LibUtilities::SessionReaderSharedPtr session;
    time_t starttime, endtime;
    NekDouble CPUtime;
   
    // Create session reader.
    session = LibUtilities::SessionReader::CreateInstance(argc, argv);

    time(&starttime);
    //----------------------------------------------------------------
    // Read the mesh and construct container class
    ADR2DManifold dom(session);
    
    // Time integration function object for unsteady equations
    LibUtilities::TimeIntegrationSchemeOperators ode;

    int nsteps = dom.GetSteps();
    int initialwavetype = dom.Getinitialwavetype();

    // Set up the intial conditions 
    dom.SetUSERDEFINEDInitialConditions(initialwavetype);

    dom.Summary(cout);
    dom.ZeroPhysFields(); // Zero phys field so that following switch is consistent
	
    switch(dom.GetEquationType())
    {   
    case eUnsteadyAdvection:
        {
            if(dom.GetExplicitAdvection())
            {
                // Choose the method of deriving forcing functions
                ode.DefineOdeRhs    (&ADR2DManifold::ODErhs,dom);
                
                // General Linear Time Integration
                dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);
            }
            else
            {
                ASSERTL0(false,"Implicit unsteady Advection Equation is not set up");
            }
        }
        break;

    case eUnsteadyDiffusion:
        {               
            if(dom.GetExplicitDiffusion())
            {                              
                // Choose the method of deriving forcing functions
                ode.DefineOdeRhs       (&ADR2DManifold::ODErhs,dom);	
                
                // General Linear Time Integration
                dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);
            }

            else
            {
                // Choose the method of deriving forcing functions
                ode.DefineImplicitSolve       (&ADR2DManifold::ODEhelmSolve,dom);		
                
                // General Linear Time Integration
                dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);
            }
        }
        break;
	
    case eUnsteadyDiffusionReaction:
        {
            if(dom.GetExplicitDiffusion())
            {
                // Choose the method of deriving forcing functions
                ode.DefineOdeRhs       (&ADR2DManifold::ODErhs,dom);	
                
                // General Linear Time Integration
                dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);
            }
            else
            {
                // Choose the method of deriving forcing functions
                ode.DefineImplicitSolve       (&ADR2DManifold::ODEhelmSolve,dom);		
                ode.DefineOdeRhs              (&ADR2DManifold::ODEeLinearMGReaction,dom);	
                
                // General Linear Time Integration
                dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);
            }
        }
        break; 

    case eNonlinearMorphogenensis:
        {
            // Choose the method of deriving forcing functions
            ode.DefineImplicitSolve       (&ADR2DManifold::ODEhelmSolve,dom);		
            ode.DefineOdeRhs              (&ADR2DManifold::ODEeNonlinearMorphoReaction,dom);	
            
            // General Linear Time Integration
            dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);
        }
        break; 

    case eFHNMONO:
        {
            // Choose the method of deriving forcing functions
            ode.DefineImplicitSolve (&ADR2DManifold::ODEhelmSolveFHNmono,dom);	
            ode.DefineOdeRhs        (&ADR2DManifold::ODEeReactionFHNmono,dom);	
          
            // General Linear Time Integration
            dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);
        }
        break;
    case eAlievPanfilov:
        {
            // Choose the method of deriving forcing functions
            ode.DefineImplicitSolve (&ADR2DManifold::ODEhelmSolveFHNmono,dom);  
            ode.DefineOdeRhs        (&ADR2DManifold::ODEeReactionAP,dom);  
          
            // General Linear Time Integration
            dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);            
        }
        break;
    case eBarkley:
        {
            // Choose the method of deriving forcing functions
            ode.DefineImplicitSolve (&ADR2DManifold::ODEhelmSolveFHNmono,dom);  
            ode.DefineOdeRhs        (&ADR2DManifold::ODEeReactionBarkley,dom);  
          
            // General Linear Time Integration
            dom.GeneralTimeIntegration(nsteps, dom.GetTimeIntMethod(), ode);            
        }
        break;

    case eNoEquationType:
    default:
        ASSERTL0(false,"Unknown or undefined equation type");
    }

    time(&endtime);
    CPUtime = (1.0/60.0/60.0)*difftime(endtime,starttime);

    // Write  output to .fld file
    dom.Output();

    cout << "-------------------------------------------" << endl;
    cout << "Total Computation Time = " << CPUtime << " hr." << endl << endl;

    // Evaluate L2 Error
    for(int i = 0; i < dom.GetNvariables(); ++i)
    {
        cout << "For variable "<< dom.GetVariable(i) <<", " << endl;
        cout << "Lint Error: " << dom.USERDEFINEDError(i,0,initialwavetype) << endl;
        cout << "H1 Error: " << dom.USERDEFINEDError(i,1,initialwavetype) << endl;
        cout << "L2 Error: " << dom.USERDEFINEDError(i,2,initialwavetype) << endl << endl;
    }

    session->Finalise();
 }

/**
* $Log $
**/
