///////////////////////////////////////////////////////////////////////////////
//
// File: LaxFriedrichsSolver.cpp
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
// Description: Lax-Friedrichs Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <ShallowWaterSolver/RiemannSolvers/LaxFriedrichsSolver.h>

namespace Nektar
{
    std::string LaxFriedrichsSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "LaxFriedrichs",
            LaxFriedrichsSolver::create,
            "LaxFriedrichs Riemann solver");

    LaxFriedrichsSolver::LaxFriedrichsSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : NonlinearSWESolver(pSession)
    {

    }

    /**
     * @brief Lax-Friedrichs Riemann solver
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
    void LaxFriedrichsSolver::v_PointSolve(
        NekDouble  hL, NekDouble  huL, NekDouble  hvL, 
        NekDouble  hR, NekDouble  huR, NekDouble  hvR, 
        NekDouble &hf, NekDouble &huf, NekDouble &hvf)
    {        
        static NekDouble g = m_params["gravity"]();
        
        // Left and right velocities
        NekDouble uL = huL / hL;
        NekDouble vL = hvL / hL;
        NekDouble uR = huR / hR;
        NekDouble vR = hvR / hR;
                
        // Left and right wave speeds
        NekDouble cL = sqrt(g * hL);
        NekDouble cR = sqrt(g * hR);
        
          // Square root of hL and hR.
        NekDouble srL  = sqrt(hL);
        NekDouble srR  = sqrt(hR);
        NekDouble srLR = srL + srR;
        
        // Velocity Roe averages
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble cRoe   = sqrt(0.5 * (cL * cL + cR * cR));
        
        // Minimum and maximum wave speeds
        NekDouble S    = std::max(uRoe+cRoe, std::max(uR+cR, -(uL-cL)));
        NekDouble sign = 1.0;
        
        if(S == -(uL-cL))
        {
            sign = -1.0;
        }
        
        // Lax-Friedrichs Riemann h flux
        hf  = 0.5 * ((huL + huR) - sign * S * (hR -hL));
        
        // Lax-Friedrichs Riemann hu flux
        huf = 0.5 * ((hL * uL * uL + 0.5 * g * hL * hL  + 
                      hR * uR * uR + 0.5 * g * hR * hR) - 
                       sign * S * (huR - huL));
        
        // Lax-Friedrichs Riemann hv flux
        hvf = 0.5 * ((hL * uL * vL + hR * uR * vR) - 
                        sign * S * (hvR - hvL));
        
    }
}
