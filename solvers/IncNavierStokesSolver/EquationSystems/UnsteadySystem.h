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

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYSYSTEM_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYSYSTEM_H

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    /// Base class for unsteady solvers.
    class UnsteadySystem : public EquationSystem
    {
    public:
        /// Destructor
        virtual ~UnsteadySystem();

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

        /// Initialises UnsteadySystem class members.
        UnsteadySystem(SessionReaderSharedPtr& pSession);

        /// Solves an unsteady problem.
        virtual void v_DoSolve();

        /// Sets up initial conditions.
        virtual void v_DoInitialise();

        /// Print a summary of time stepping parameters.
        virtual void v_PrintSummary(std::ostream &out);

        ///
        virtual void v_NumericalFlux(
                    Array<OneD, Array<OneD, NekDouble> > &physfield,
                    Array<OneD, Array<OneD, NekDouble> > &numflux);

        ///
        virtual void v_NumericalFlux(
                    Array<OneD, Array<OneD, NekDouble> > &physfield,
                    Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                    Array<OneD, Array<OneD, NekDouble> > &numfluxY );

        ///
        virtual void v_NumFluxforScalar(
                    Array<OneD, Array<OneD, NekDouble> > &ufield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux);

        ///
        virtual void v_NumFluxforVector(
                    Array<OneD, Array<OneD, NekDouble> > &ufield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                    Array<OneD, Array<OneD, NekDouble> > &qflux);

        /// Evaulate flux = m_fields*ivel for i th component of Vu for
        /// direction j
        virtual void v_GetFluxVector(const int i, const int j,
                    Array<OneD, Array<OneD, NekDouble> > &physfield,
                    Array<OneD, Array<OneD, NekDouble> > &flux);

    private:
        ///
        void WeakPenaltyforScalar(const int var,
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

#endif
