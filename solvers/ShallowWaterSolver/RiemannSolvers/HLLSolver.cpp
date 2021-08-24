///////////////////////////////////////////////////////////////////////////////
//
// File: HLLSolver.cpp
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
// Description: HLL Riemann solver for the Nonlinear Shallow Water Equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <ShallowWaterSolver/RiemannSolvers/HLLSolver.h>

namespace Nektar
{
    std::string HLLSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "HLL",
            HLLSolver::create,
            "HLL Riemann solver");

    HLLSolver::HLLSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : NonlinearSWESolver(pSession)
    {

    }

    /**
     * @brief HLL Riemann solver for the Nonlinear Shallow Water Equations
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
    void HLLSolver::v_PointSolve(
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
      

        // Left and right wave speeds
        NekDouble cL = sqrt(g * hL);
        NekDouble cR = sqrt(g * hR);
    
        // the two-rarefaction wave assumption
        NekDouble hstar,fL,fR;
        hstar = 0.5 * (cL + cR) + 0.25 * (uL - uR);
        hstar *= hstar;
        hstar *= (1.0/g);

    
        // Compute SL
        NekDouble SL;
        if (hstar > hL)
          SL = uL - cL * sqrt(0.5*((hstar*hstar + hstar*hL)/(hL*hL)));
        else
          SL = uL - cL;
    
        // Compute SR
        NekDouble SR;
        if (hstar > hR)
          SR = uR + cR * sqrt(0.5*((hstar*hstar + hstar*hR)/(hR*hR)));
        else
          SR = uR + cR;
    
        if (SL >= 0)
          {
            hf  = hL * uL;
            huf  = uL * uL * hL + 0.5 * g * hL * hL;
            hvf  = hL * uL * vL;
          }
        else if (SR <= 0)
          {
            hf  = hR * uR;
            huf  = uR * uR * hR + 0.5 * g * hR * hR;
            hvf  = hR * uR *vR;
          }
        else 
          {
            hf = (SR * hL * uL - SL * hR * uR + 
                  SL * SR * (hR - hL)) / (SR - SL);
            fL =  uL * uL * hL + 0.5 * g * hL * hL;
            fR =  uR * uR * hR + 0.5 * g * hR * hR;
            huf =(SR * fL - SL * fR + 
                  SL * SR * (hR * uR - hL * uL)) / (SR - SL);
            fL =  uL * vL * hL;
            fR =  uR * vR * hR;
            hvf =(SR * fL - SL * fR + 
                  SL * SR * (hR * vR - hL * vL)) / (SR - SL);
          }
    }
}
