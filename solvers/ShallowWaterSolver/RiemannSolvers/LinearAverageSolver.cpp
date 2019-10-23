///////////////////////////////////////////////////////////////////////////////
//
// File: LinearAverageSolver.cpp
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
// Description: Simple mean value solver for the Linear Shallow Water 
//              Equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <ShallowWaterSolver/RiemannSolvers/LinearAverageSolver.h>

namespace Nektar
{
    std::string LinearAverageSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "LinearAverage",
            LinearAverageSolver::create,
            "Average Value Riemann solver");

    LinearAverageSolver::LinearAverageSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : LinearSWESolver(pSession)
    {

    }

    /**
     * @brief Average Value Riemann solver for the Linear Shallow 
     * Water Equations
     *
     * @param etaL   Free surface elevation left state.
     * @param etaR   Free surface elevation right state.  
     * @param uL     x-velocity  left state.  
     * @param uR     x-velocity  right state.  
     * @param vL     y-velocity  left state.  
     * @param vR     y-velocity  right state. 
     * @param dL     still water depth component left state.  
     * @param dR     still water depth component right state. 
     * @param etaf   Computed Riemann flux for continuity.
     * @param uf     Computed Riemann flux for x-momentum component 
     * @param vf     Computed Riemann flux for y-momentum component 
     */
    void LinearAverageSolver::v_PointSolve(
        NekDouble  etaL, NekDouble  uL, NekDouble  vL, NekDouble dL,
        NekDouble  etaR, NekDouble  uR, NekDouble  vR, NekDouble dR,
        NekDouble &etaf, NekDouble &uf, NekDouble &vf)
    {
        boost::ignore_unused(vL, vR);

        static NekDouble g = m_params["gravity"]();

        etaf  = 0.5 * (dL * uL + dR * uR);
        uf    = 0.5 * (g * etaL + g * etaR);
        vf    = 0.0;
    }
}
