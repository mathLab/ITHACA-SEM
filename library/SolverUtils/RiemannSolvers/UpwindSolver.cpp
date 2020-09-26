///////////////////////////////////////////////////////////////////////////////
//
// File: UpwindSolver.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Upwind Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/RiemannSolvers/UpwindSolver.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string UpwindSolver::solverName = GetRiemannSolverFactory().
            RegisterCreatorFunction("Upwind", UpwindSolver::create, "Upwind solver");

        /**
         * @class UpwindSolver
         *
         * @brief Upwind scheme Riemann solver.
         *
         * The upwind solver determines the flux based upon an advection
         * velocity \f$\mathbf{V}\f$ and trace normal \f$\mathbf{n}\f$. In
         * particular, the flux for each component of the velocity field is
         * deterined by:
         *
         * \f[ \mathbf{f}(u^+,u^-) = \begin{cases} \mathbf{V}u^+, &
         * \mathbf{V}\cdot\mathbf{n} \geq 0,\\ \mathbf{V}u^-, &
         * \mathbf{V}\cdot\mathbf{n} < 0.\end{cases} \f]
         *
         * Here the superscript + and - denotes forwards and backwards spaces
         * respectively.
         */

        /**
         * @brief Default constructor.
         */
        UpwindSolver::UpwindSolver(
            const LibUtilities::SessionReaderSharedPtr& pSession)
            : RiemannSolver(pSession)
        {
        }

        /**
         * @brief Implementation of the upwind solver.
         *
         * The upwind solver assumes that a scalar field Vn is defined, which
         * corresponds with the dot product \f$\mathbf{V}\cdot\mathbf{n}\f$,
         * where \f$\mathbf{V}\f$ is the advection velocity and \f$\mathbf{n}\f$
         * defines the normal of a vertex, edge or face at each quadrature point
         * of the trace space.
         *
         * @param Fwd   Forwards trace space.
         * @param Bwd   Backwards trace space.
         * @param flux  Resulting flux.
         */
        void UpwindSolver::v_Solve(
            const int                                         nDim,
            const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
            const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                  Array<OneD,       Array<OneD, NekDouble> > &flux)
        {
            boost::ignore_unused(nDim);

            ASSERTL1(CheckScalars("Vn"), "Vn not defined.");
            const Array<OneD, NekDouble> &traceVel = m_scalars["Vn"]();

            for (int j = 0; j < traceVel.size(); ++j)
            {
                const Array<OneD, const Array<OneD, NekDouble> > &tmp =
                    traceVel[j] >= 0 ? Fwd : Bwd;
                for (int i = 0; i < Fwd.size(); ++i)
                {
                    flux[i][j] = traceVel[j]*tmp[i][j];
                }
            }
        }
    }
}
