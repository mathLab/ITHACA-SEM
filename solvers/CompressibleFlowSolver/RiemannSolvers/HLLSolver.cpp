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
// Description: HLL Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/HLLSolver.h>

namespace Nektar
{
    std::string HLLSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "HLL",
            HLLSolver::create,
            "HLL Riemann solver");

    HLLSolver::HLLSolver() : CompressibleSolver()
    {

    }

    /**
     * @brief HLL Riemann solver
     *
     * @param rhoL      Density left state.
     * @param rhoR      Density right state.  
     * @param rhouL     x-velocity component left state.  
     * @param rhouR     x-velocity component right state.  
     * @param rhovL     y-velocity component left state.  
     * @param rhovR     y-velocity component right state.  
     * @param rhowL     z-velocity component left state.  
     * @param rhowR     z-velocity component right state.
     * @param EL        Energy left state.  
     * @param ER        Energy right state. 
     * @param rhof      Riemann flux for density (i.e. first equation).  
     * @param rhouf     Riemann flux for x-velocity component 
     *                  (i.e. second equation).
     * @param rhovf     Riemann flux for y-velocity component 
     *                  (i.e. third equation).
     * @param rhowf     Riemann flux for z-velocity component 
     *                  (i.e. fourth equation).
     * @param Ef        Riemann flux for energy (i.e. fifth equation).  
     *
     */
    void HLLSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {        
        NekDouble gamma = m_params["gamma"]();
        
        // Left and Right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;

        // Left and right pressures
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        NekDouble hL = (EL + pL) / rhoL;
        NekDouble hR = (ER + pR) / rhoR;
        
        // Density Roe average
        NekDouble rhoRoe = sqrt(rhoL) * sqrt(rhoR);
        
        // Velocity Roe averages
        NekDouble uRoe   = (sqrt(rhoL) * uL + 
                            sqrt(rhoR) * uR) / (sqrt(rhoL) + sqrt(rhoR));
        NekDouble vRoe   = (sqrt(rhoL) * vL + 
                            sqrt(rhoR) * vR) / (sqrt(rhoL) + sqrt(rhoR));
        NekDouble wRoe   = (sqrt(rhoL) * wL + 
                            sqrt(rhoR) * wR) / (sqrt(rhoL) + sqrt(rhoR));

        // Entalpy Roe average
        NekDouble hRoe   = (sqrt(rhoL) * hL + 
                            sqrt(rhoR) * hR) / (sqrt(rhoL) + sqrt(rhoR));
        NekDouble cRoe   = sqrt((gamma - 1.0) *
                                (hRoe - 0.5 * (uRoe * uRoe + vRoe * vRoe +
                                               wRoe * wRoe)));
        
        // Maximum wave speeds
        NekDouble SL = std::min(uL-cL, uRoe-cRoe);
        NekDouble SR = std::max(uR+cR, uRoe+cRoe);
        
        // HLL Riemann fluxes (positive case)
        if (SL >= 0)
        {
            rhof  = rhoL * uL;
            rhouf = rhoL * uL * uL + pL;
            rhovf = rhoL * uL * vL;
            rhowf = rhoL * uL * wL;
            Ef    = uL * (EL + pL);
        }
        // HLL Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            rhof  = rhoR * uR;
            rhouf = rhoR * uR * uR + pR;
            rhovf = rhoR * uR * vR;
            rhowf = rhoR * uR * wR;
            Ef    = uR * (ER + pR);
        }
        // HLL Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            rhof  = (SR * (rhoL * uL) - 
                     SL * (rhoR * uR) + 
                     SR * SL * (rhoR - rhoL)) / (SR - SL);
            
            rhouf = (SR * (rhoL * uL * uL + pL) - 
                     SL * (rhoR * uR * uR + pR) + 
                     SR * SL * (rhouR - rhouL)) / (SR - SL);
            
            rhovf = (SR * (rhoL * uL * vL) - 
                     SL * (rhoR * uR * vR) + 
                     SR * SL * (rhovR - rhovL)) / (SR - SL);
            
            Ef    = (SR * (uL * EL + uL * pL) - 
                     SL * (uR * ER + uR * pR) + 
                     SR * SL * (ER - EL)) / (SR - SL);
        }
    }
}
