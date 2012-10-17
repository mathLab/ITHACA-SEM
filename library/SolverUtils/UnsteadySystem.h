///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadySystem.h
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
// Description: Generic timestepping for Unsteady solvers
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_UNSTEADYSYSTEM_H
#define NEKTAR_SOLVERUTILS_UNSTEADYSYSTEM_H

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Filters/Filter.h>
#include <time.h>

namespace Nektar
{
    namespace SolverUtils
    {
        /// Base class for unsteady solvers.
        class UnsteadySystem : public EquationSystem
        {
        public:
            /// Destructor
            SOLVER_UTILS_EXPORT virtual ~UnsteadySystem();
		
            /// Calculate the largest time-step for maintaining a stable
            /// solution (using CFL).
            SOLVER_UTILS_EXPORT NekDouble GetTimeStep(
                const Array<OneD,int> ExpOrder, 
                const Array<OneD,NekDouble> CFL, 
                NekDouble timeCFL);
            
            /// Calculate the larger time-step mantaining the problem stable, just for CFLTester equation
            SOLVER_UTILS_EXPORT NekDouble GetTimeStep(
                int ExpOrder, 
                NekDouble CFL, 
                NekDouble TimeStability);
		
            /// CFL number
            NekDouble m_cfl;
		
            // Mapping of the real convective field on the standard element.
            // This function gives back the convective filed in the standard
            // element to calculate the stability region of the problem in a
            // unique way.
            SOLVER_UTILS_EXPORT Array<OneD,NekDouble> GetStdVelocity(
                const Array<OneD, Array<OneD,NekDouble> > inarray);
            
            /// Function to calculate the stability limit for DG/CG
            SOLVER_UTILS_EXPORT NekDouble GetStabilityLimit(int n);
            
	        /// Function to calculate the stability limit for DG/CG (a vector
	        /// of them)
            SOLVER_UTILS_EXPORT Array<OneD,NekDouble> 
                GetStabilityLimitVector(
                    const Array<OneD,int> &ExpOrder);
            
        protected:
            /// Number of time steps between outputting status information.
            int                                             m_infosteps;
            /// The time integration method to use.
            LibUtilities::TimeIntegrationMethod             m_timeIntMethod;
            /// The time integration scheme operators to use.
            LibUtilities::TimeIntegrationSchemeOperators    m_ode;
            ///
            NekDouble                                       m_epsilon;
            /// Indicates if explicit or implicit treatment of diffusion is used.
            bool                                            m_explicitDiffusion;
            /// Indicates if explicit or implicit treatment of advection is used.
            bool                                            m_explicitAdvection;
            /// Indicates if explicit or implicit treatment of reaction is used.
            bool                                            m_explicitReaction;

            std::vector<int>                                m_intVariables;

            std::vector<FilterSharedPtr>                    m_filters;

            /// Initialises UnsteadySystem class members.
            SOLVER_UTILS_EXPORT UnsteadySystem(const LibUtilities::SessionReaderSharedPtr& pSession);

            SOLVER_UTILS_EXPORT virtual void v_InitObject();

            /// Solves an unsteady problem.
            SOLVER_UTILS_EXPORT virtual void v_DoSolve();

            /// Sets up initial conditions.
            SOLVER_UTILS_EXPORT virtual void v_DoInitialise();

            /// Print a summary of time stepping parameters.
            SOLVER_UTILS_EXPORT virtual void v_PrintSummary(std::ostream &out);
            
            /// Print the solution at each solution point in a txt file
            SOLVER_UTILS_EXPORT virtual void v_AppendOutput1D(Array<OneD, Array<OneD, NekDouble> > &solution1D);

            ///
            SOLVER_UTILS_EXPORT virtual void v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numflux);

            ///
            SOLVER_UTILS_EXPORT virtual void v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                Array<OneD, Array<OneD, NekDouble> > &numfluxY );

            ///
            SOLVER_UTILS_EXPORT virtual void v_NumFluxforScalar(
                Array<OneD, Array<OneD, NekDouble> > &ufield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux);

            ///
            SOLVER_UTILS_EXPORT virtual void v_NumFluxforVector(
                Array<OneD, Array<OneD, NekDouble> > &ufield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                Array<OneD, Array<OneD, NekDouble> > &qflux);

            /// Evaulate flux = m_fields*ivel for i th component of Vu for
            /// direction j
            SOLVER_UTILS_EXPORT virtual void v_GetFluxVector(
                const int i, 
                const int j,
                      Array<OneD, Array<OneD, NekDouble> > &physfield,
                      Array<OneD, Array<OneD, NekDouble> > &flux);
		
            /// Virtual function to get the time step
            SOLVER_UTILS_EXPORT virtual NekDouble v_GetTimeStep(
                const Array<OneD,int> ExpOrder, 
                const Array<OneD,NekDouble> CFL, NekDouble timeCFL);
		
            /// Virtual function to get the time step
			SOLVER_UTILS_EXPORT virtual NekDouble v_GetTimeStep(
                int ExpOrder, 
                NekDouble CFL, 
                NekDouble TimeStability);

        private:
            ///
            void WeakPenaltyforScalar(
                const int var,
                const Array<OneD, const NekDouble> &physfield,
                      Array<OneD,       NekDouble> &penaltyflux,
                      NekDouble time=0.0);

            ///
            void WeakPenaltyforVector(
                const int var,
                const int dir,
                const Array<OneD, const NekDouble> &physfield,
                Array<OneD,       NekDouble> &penaltyflux,
                NekDouble C11,
                NekDouble time=0.0);
        };
    }
}

#endif
