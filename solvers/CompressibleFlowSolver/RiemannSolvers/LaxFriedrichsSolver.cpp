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
     * @param rhouL     x-momentum component left state.  
     * @param rhouR     x-momentum component right state.  
     * @param rhovL     y-momentum component left state.  
     * @param rhovR     y-momentum component right state.  
     * @param rhowL     z-momentum component left state.  
     * @param rhowR     z-momentum component right state.
     * @param EL        Energy left state.  
     * @param ER        Energy right state. 
     * @param rhof      Computed Riemann flux for density.
     * @param rhouf     Computed Riemann flux for x-momentum component 
     * @param rhovf     Computed Riemann flux for y-momentum component 
     * @param rhowf     Computed Riemann flux for z-momentum component 
     * @param Ef        Computed Riemann flux for energy.
     */
    void LaxFriedrichsSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {        
        static NekDouble gamma = m_params["gamma"]();
        
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
        
        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;
        
        // Velocity Roe averages
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble hRoe   = (srL * hL + srR * hR) / srLR;
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * 
            (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe)));
        
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
