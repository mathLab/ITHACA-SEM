///////////////////////////////////////////////////////////////////////////////
//
// File: AverageSolver.cpp
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
// Description: Simple mean value solver for the Nonlinear Shallow Water 
//              Equations.
//              Do not use for anything else than testing and comparisons
//
///////////////////////////////////////////////////////////////////////////////

#include <ShallowWaterSolver/RiemannSolvers/AverageSolver.h>

namespace Nektar
{
    std::string AverageSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "Average",
            AverageSolver::create,
            "Average Value Riemann solver");

    AverageSolver::AverageSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : NonlinearSWESolver(pSession)
    {

    }

    /**
     * @brief Average Value Riemann solver for the Nonlinear Shallow 
     * Water Equations
     *
     * @param hL      Water depth left state.
     * @param hR      Water depth right state.  
     * @param huL     x-momentum component left state.  
     * @param huR     x-momentum component right state.  
     * @param hvL     y-momentum component left state.  
     * @param hvR     y-momentum component right state.  
     * @param hf      Computed Riemann flux for density.
     * @param huf     Computed Riemann flux for x-momentum component 
     * @param hvf     Computed Riemann flux for y-momentum component 
     */
    void AverageSolver::v_PointSolve(
        NekDouble  hL, NekDouble  huL, NekDouble  hvL,
        NekDouble  hR, NekDouble  huR, NekDouble  hvR,
        NekDouble &hf, NekDouble &huf, NekDouble &hvf)
    {        
        static NekDouble g = m_params["gravity"]();
        
        // Left and Right velocities
        NekDouble uL = huL / hL;
        NekDouble vL = hvL / hL;
	NekDouble uR = huR / hR;
        NekDouble vR = hvR / hR;
      
	// Compute fluxes
	hf  = 0.5 * (hL * uL + hR * uR);
	huf = 0.5 * ((uL * uL * hL + 0.5 * g * hL * hL) +
		     (uR * uR * hR + 0.5 * g * hR * hR));
	hvf = 0.5 * (hL * uL * vL + hR * uR * vR);
    }
}
