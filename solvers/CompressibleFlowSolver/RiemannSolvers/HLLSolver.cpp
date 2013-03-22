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
    void HLLSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {        
        static NekDouble gamma = m_params["gamma"]();
        
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

        // Maximum wave speeds
        NekDouble SL = std::min(uL-cL, uRoe-cRoe);
        NekDouble SR = std::max(uR+cR, uRoe+cRoe);

        // HLL Riemann fluxes (positive case)
        if (SL >= 0)
        {
            rhof  = rhouL;
            rhouf = rhouL * uL + pL;
            rhovf = rhouL * vL;
            rhowf = rhouL * wL;
            Ef    = uL * (EL + pL);
        }
        // HLL Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            rhof  = rhouR;
            rhouf = rhouR * uR + pR;
            rhovf = rhouR * vR;
            rhowf = rhouR * wR;
            Ef    = uR * (ER + pR);
        }
        // HLL Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            NekDouble tmp1 = 1.0 / (SR - SL);
            NekDouble tmp2 = SR * SL;
            rhof  = (SR * rhouL - SL * rhouR + tmp2 * (rhoR - rhoL)) * tmp1;
            rhouf = (SR * (rhouL * uL + pL) - SL * (rhouR * uR + pR) +
                     tmp2 * (rhouR - rhouL)) * tmp1;
            rhovf = (SR * rhouL * vL - SL * rhouR * vR +
                     tmp2 * (rhovR - rhovL)) * tmp1;
            rhowf = (SR * rhouL * wL - SL * rhouR * wR +
                     tmp2 * (rhowR - rhowL)) * tmp1;
            Ef    = (SR * uL * (EL + pL) - SL * uR * (ER + pR) + 
                     tmp2 * (ER - EL)) * tmp1;
        }
    }
}
