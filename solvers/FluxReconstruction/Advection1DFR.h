///////////////////////////////////////////////////////////////////////////////
//
// File Advection1DFR.h
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

#ifndef NEKTAR_SOLVERS_FluxReconstruction_Advection1DFR_H
#define NEKTAR_SOLVERS_FluxReconstruction_Advection1DFR_H

#include <cstdio>
#include <fstream>
#include <fstream>
#include <iomanip>
#include <MultiRegions/DisContField1D.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/Polylib/Polylib.h>

namespace Nektar
{
	class Advection1DFR
	{
	public:
		
		/// Constructor
		Advection1DFR(LibUtilities::SessionReaderSharedPtr &vSession);
		
        /// Destructor
		~Advection1DFR(void);
	
		/// Function which evaluate the advection term
		void EvaluateAdvectionTerm(const Array<OneD, Array<OneD, NekDouble> > & inarray,
								   Array<OneD, Array<OneD, NekDouble> > & outarray,
								   const NekDouble time);
	
		/// Projection
		void Projection(const Array<OneD, Array<OneD, NekDouble> > & inarray,
						Array<OneD, Array<OneD, NekDouble> > & outarray,
						const NekDouble time) const;
		
		/// Riemann solver to work out the interface fluxes
		void RiemannSolver(const Array<OneD, NekDouble> & inarray, 
						   Array<OneD, NekDouble> & F);
		
		/// Function which calculate the jumps and return the final flux = df/dx
		/// corrected with dGL and dGR derivatives
		void FluxesReconstruction(const Array<OneD, NekDouble> & fd,
								  const Array<OneD, NekDouble> & fdi,
								  const Array<OneD, NekDouble> & fi,
								  Array<OneD, NekDouble> & outarray);
		
		/// Function which calculate the correction function derivatives at the beginning and store the value
		void GFunctionsGrad(Array<OneD, NekDouble> & dGL,
				            Array<OneD, NekDouble> & dGR);
		
		/// Function which interpolate the total value (total) to the interface value (interface)
		void InterpToInterface(const Array<OneD, NekDouble> & total,
							   Array<OneD, NekDouble> & interfaceValue);
		
		/// Function to print out a summary and the results (the solution to file and the error to screen)
		void SolutionPrint();
		
		void AppendOutput(const Array<OneD, NekDouble> & approx,
						  const Array<OneD, NekDouble> & exact) const;
		
		void GenerateGnuplotScript() const;
		
        void GenerateMatlabScript() const;
        
        
		/// Functions to access some variable inside the problem
		MultiRegions::DisContField1DSharedPtr GetDomain() const;
		NekDouble GetTime() const;
		void UpdateTime();
		NekDouble GetTimeStep() const;
		int GetNumSteps() const;
		int GetNumPoints() const;
		
		// Variables
		MultiRegions::DisContField1DSharedPtr   Domain;     // Object in which the mesh, the projection type and the expansion are joined
		SpatialDomains::MeshGraphSharedPtr      graph1D;    // Object containg the mesh
		Array<OneD,NekDouble>                   x;          // X-coordinates of the quadrature points
		Array<OneD,NekDouble>					y;          // Y-coordinates of the quadrature points
		Array<OneD,NekDouble>					z;          // Z-coordinates of the quadrature points
		Array<OneD,NekDouble>					xi;         // X-coordinates of the interface points
		Array<OneD,NekDouble>					yi;         // Y-coordinates of the interface points
		Array<OneD,NekDouble>					zi;         // Z-coordinates of the interface points
		int										nq;         // Number of quadrature points
		int										ni;         // Number of interface points
		int										ne;         // Number of elements
		NekDouble								Time;       // Time-level
		NekDouble								TimeStep;   // Time-step
		int                                     NumSteps;   // Number of time-step
		LibUtilities::EquationSharedPtr         AdveFunc;   // Equation holding the advection definition
		LibUtilities::EquationSharedPtr         InitCond;   // Equation holding the initial condition definition
		LibUtilities::EquationSharedPtr         ExSol;      // Equation holding the exact solution definition
		Array<OneD, Array<OneD, NekDouble> >    Adv;        // Vector containg the values of the advection coefficient at the quadrature
		Array<OneD, Array<OneD, NekDouble> >    Ini;        // Vector containg the values of the initial condition at the quadrature
		Array<OneD, Array<OneD, NekDouble> >    Exac;       // Vector containg the values of the exact solution at the quadrature
		std::string                             RiemSol;    // String which defines the Riemann solver 
		std::string                             GFtype;     // Type of correction functions
		Array<OneD, NekDouble>                  K;          // Array of K which defines which type of Riemann solver
		Array<OneD, NekDouble>                  dgL;        // Left correction function derived by x
		Array<OneD, NekDouble>                  dgR;        // Right correction function derived by x
		
	private:
		
	};
}
#endif
