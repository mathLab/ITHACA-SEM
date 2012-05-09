///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvection1DFRSolver.cpp
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
// Description: Demo to test the Flux Reconstruction 1D
// Unsteady advection (linear or spatially non linear advection) 
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <FluxReconstruction/Advection1DFR.h>

using namespace Nektar;

int main(int argc, char *argv[])
{

	/// Reading in the session file
	LibUtilities::SessionReaderSharedPtr pSession = LibUtilities::SessionReader::CreateInstance(argc, argv);
	
	/// Checking usage
	if((argc != 2)&&(argc != 3))
    {
        fprintf(stderr,"Usage: ./UnsteadyAdvection1DFRSolver meshfile [SysSolnType]\n");
        exit(1);
    }
	
	/// Instantiating the object to hold the spatial discretisation
	Advection1DFR solver(pSession);
	
	/// The time integration method to use.
	LibUtilities::TimeIntegrationMethod  TIME_SCHEME;
	/// The time integration scheme operators to use.
	LibUtilities::TimeIntegrationSchemeOperators ODE;
	
	/// Defining the RHS operator, i.e. the advection oprator
	/// the projection operator is compulsory for CG because we need to reinforce BC
	/// for DG is just a "copy" of the degrees of freedom
	ODE.DefineOdeRhs(&Advection1DFR::EvaluateAdvectionTerm,solver);
	ODE.DefineProjection(&Advection1DFR::Projection,solver);
	
	// We may need an array of time-integration schemes, in case we are using
	// high order schemes which need other schemes to be initialised
	Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
	
	/// Solution of the ODE system 
	LibUtilities::TimeIntegrationSolutionSharedPtr ode_solution;
	
	/// Creating a vector in which we strore our solution step by step.
	// At the really beginning we put in equal to the physical value of the object Domain
	// inside Advection1DFR, where we have stored the initial solution
	Array<OneD, Array<OneD, NekDouble> > U(1);
	U[0] = Array<OneD, NekDouble>(solver.GetNumPoints());
	
	U[0] = solver.GetDomain()->GetPhys();
	
	/// Determining TimeIntegrationMethod to use, specified in the session file
	for (int i = 0; i < (int)LibUtilities::SIZE_TimeIntegrationMethod; ++i)
	{
		bool match;
		pSession->MatchSolverInfo("TIMEINTEGRATIONMETHOD",LibUtilities::TimeIntegrationMethodMap[i], match, false);
		if (match)
		{
			TIME_SCHEME = (LibUtilities::TimeIntegrationMethod) i;
			break;
		}
	}
	
	// So far just one-step schemes are taken in account for simplicity,
	// later on we can add other schemes, it is quite fast
	int numMultiSteps;
	switch(TIME_SCHEME)
	{
        case LibUtilities::eForwardEuler:
        case LibUtilities::eClassicalRungeKutta4:
		case LibUtilities::eRungeKutta2_ImprovedEuler:
		case LibUtilities::eAdamsBashforthOrder1:
		{
			numMultiSteps = 1;
			
			IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
			
			LibUtilities::TimeIntegrationSchemeKey IntKey(TIME_SCHEME);
			IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
			
			ode_solution = IntScheme[0]->InitializeScheme(solver.GetTimeStep(),U,solver.GetTime(),ODE);
			break;
		}
		default:
		{
			ASSERTL0(false,"This Time-Integration scheme has not been set up for this problem");
		}
	}
	
    /// ------------------------------------------------------
	/// Performing the time integration
    /// ------------------------------------------------------
	for(int n = 0; n < solver.GetNumSteps(); ++n)
	{
		U = IntScheme[0]->TimeIntegrate(solver.GetTimeStep(),ode_solution,ODE);
		solver.UpdateTime();
		
		/// Coping the solution in the member variable m_phys of the object Domain
	    solver.GetDomain()->UpdatePhys() = U[0];
		
        /// Updating the variable m_coeffs of the object Domain (the coefficinet space)
		solver.GetDomain()->FwdTrans(U[0],solver.GetDomain()->UpdateCoeffs());
	}

	solver.SolutionPrint();
	
	pSession->Finalise();
}
