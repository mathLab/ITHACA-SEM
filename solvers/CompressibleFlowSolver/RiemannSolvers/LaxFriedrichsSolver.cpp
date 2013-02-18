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
// Description: Lax-Friedrichs Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/LaxFriedrichsSolver.h>

namespace Nektar
{
    std::string LaxFriedrichsSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "LaxFriedrichs",
            LaxFriedrichsSolver::create,
            "LaxFriedrichs Riemann solver");

    LaxFriedrichsSolver::LaxFriedrichsSolver() : CompressibleSolver()
    {

    }

    /**
     * @brief Lax-Friedrichs Riemann solver
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
    void LaxFriedrichsSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {        
        NekDouble gamma = m_params["gamma"]();
        
        // Left and right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;
        
        // Left and right pressures
        NekDouble pL = (gamma - 1.0) * 
                       (EL - 0.5 * (rhouL*uL + rhovL*vL + rhowL*wL));
        
        NekDouble pR = (gamma - 1.0) * 
                       (ER - 0.5 * (rhouR*uR + rhovR*vR + rhowR*wR));
        
        // Left and right speeds of sound
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        
        // Left and right entalpies
        NekDouble hL = (EL + pL) / rhoL;
        NekDouble hR = (ER + pR) / rhoR;
        
        // Density Roe average
        NekDouble rhoRoe = sqrt(rhoL) * sqrt(rhoR);
        
        // Velocity Roe averages
        NekDouble uRoe   = (sqrt(rhoL) * uL + 
                            sqrt(rhoR) * uR) / (sqrt(rhoL) + sqrt(rhoR));
        
        NekDouble vRoe   = (sqrt(rhoL) * vL +
                            sqrt(rhoR) * vR) / (sqrt(rhoL) + sqrt(rhoR));
        
        NekDouble wRoe   = (sqrt(rhoL) * vL + 
                            sqrt(rhoR) * vR) / (sqrt(rhoL) + sqrt(rhoR));
        
        // Entalpy Roe average
        NekDouble hRoe   = (sqrt(rhoL) * hL + 
                            sqrt(rhoR) * hR) / (sqrt(rhoL) + sqrt(rhoR));
        
        // Speed of sound Roe average
        NekDouble cRoe   = sqrt((gamma - 1.0) * (hRoe - 0.5 * 
                                                 (uRoe * uRoe + vRoe * vRoe)));
        
        // Minimum and maximum wave speeds
        NekDouble S    = std::max(uRoe+cRoe, std::max(uR+cR, -uL+cL));
        NekDouble sign = 1.0;
        
        if(S == -uL+cL)
        {
            sign = -1.0;
        }
        
        // Lax-Friedrichs Riemann rho flux
        rhof  = 0.5 * ((rhouL + rhouR) - sign * S * (rhoR -rhoL));
        
        // Lax-Friedrichs Riemann rhou flux
        rhouf = 0.5 * ((rhoL * uL * uL + pL + rhoR * uR * uR + pR) - 
                       sign * S * (rhouR - rhouL));
        
        // Lax-Friedrichs Riemann rhov flux
        rhovf = 0.5 * ((rhoL * uL * vL + rhoR * uR * vR) - 
                        sign * S * (rhovR - rhovL));
        
        // Lax-Friedrichs Riemann rhow flux
        rhowf = 0.5 * ((rhoL * uL * wL + rhoR * uR * wR) - 
                        sign * S * (rhowR - rhowL));
        
        // Lax-Friedrichs Riemann E flux
        Ef    = 0.5 * ((uL * (EL + pL) + uR * (ER + pR)) - 
                        sign * S * (ER - EL));
    }
}
