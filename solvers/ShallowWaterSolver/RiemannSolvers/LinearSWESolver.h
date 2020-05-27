///////////////////////////////////////////////////////////////////////////////
//
// File: LinearSWESolver.h
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
// Description: Linearized Shallow Water Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_RIEMANNSOLVER_LINEARSWESOLVER
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_RIEMANNSOLVER_LINEARSWESOLVER

#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{
    class LinearSWESolver : public RiemannSolver
    {
    protected:
        bool m_pointSolve;

        LinearSWESolver(const LibUtilities::SessionReaderSharedPtr& pSession);
        
        virtual void v_Solve(
            const int                                         nDim,
            const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
            const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                  Array<OneD,       Array<OneD, NekDouble> > &flux);
        
        virtual void v_ArraySolve(
            const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
            const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                  Array<OneD,       Array<OneD, NekDouble> > &flux)
        {
            boost::ignore_unused(Fwd, Bwd, flux);
            NEKERROR(ErrorUtil::efatal,
                     "This function should be defined by subclasses.");
        }

        virtual void v_PointSolve(
            NekDouble  etaL, NekDouble  uL, NekDouble  vL, NekDouble dL,
            NekDouble  etaR, NekDouble  uR, NekDouble  vR, NekDouble dR,
            NekDouble &etaf, NekDouble &uf, NekDouble &vf)
        {
            boost::ignore_unused(etaL, uL, vL, dL, etaR, uR, vR, dR, etaf,
                                 uf, vf);
            NEKERROR(ErrorUtil::efatal,
                     "This function should be defined by subclasses.");
        }
    };
}

#endif
