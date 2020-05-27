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

    HLLSolver::HLLSolver(const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleSolver(pSession)
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
        // Left and Right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;

        // Internal energy (per unit mass)
        NekDouble eL =
                (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL)) / rhoL;
        NekDouble eR =
                (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR)) / rhoR;
        // Pressure
        NekDouble pL = m_eos->GetPressure(rhoL, eL);
        NekDouble pR = m_eos->GetPressure(rhoR, eR);
        // Speed of sound
        NekDouble cL = m_eos->GetSoundSpeed(rhoL, eL);
        NekDouble cR = m_eos->GetSoundSpeed(rhoR, eR);

        // Left and right total enthalpy
        NekDouble HL = (EL + pL) / rhoL;
        NekDouble HR = (ER + pR) / rhoR;

        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;
        
        // Roe average state
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble URoe2  = uRoe*uRoe + vRoe*vRoe + wRoe*wRoe;
        NekDouble HRoe   = (srL * HL + srR * HR) / srLR;
        NekDouble cRoe   = GetRoeSoundSpeed(
                                rhoL, pL, eL, HL, srL,
                                rhoR, pR, eR, HR, srR,
                                HRoe, URoe2, srLR);

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
