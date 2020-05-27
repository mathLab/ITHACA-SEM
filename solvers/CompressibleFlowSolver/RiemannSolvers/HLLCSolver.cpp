///////////////////////////////////////////////////////////////////////////////
//
// File: HLLCSolver.cpp
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
// Description: HLLC Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/HLLCSolver.h>

namespace Nektar
{
    std::string HLLCSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "HLLC", HLLCSolver::create, "HLLC Riemann solver");
    
    HLLCSolver::HLLCSolver(const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleSolver(pSession)
    {
        
    }
    
    /**
     * @brief HLLC Riemann solver
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
    void HLLCSolver::v_PointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef)
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
        
        // HLLC Riemann fluxes (positive case)
        if (SL >= 0)
        {
            rhof  = rhouL;
            rhouf = rhouL * uL + pL;
            rhovf = rhouL * vL;
            rhowf = rhouL * wL;
            Ef    = uL * (EL + pL);
        }
        // HLLC Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            rhof  = rhouR;
            rhouf = rhouR * uR + pR;
            rhovf = rhouR * vR;
            rhowf = rhouR * wR;
            Ef    = uR * (ER + pR);
        }
        // HLLC Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            NekDouble SM = (pR - pL + rhouL * (SL - uL) - rhouR * (SR - uR)) /
            (rhoL * (SL - uL) - rhoR * (SR - uR));
            NekDouble rhoML  = rhoL * (SL - uL) / (SL - SM);
            NekDouble rhouML = rhoML * SM;
            NekDouble rhovML = rhoML * vL;
            NekDouble rhowML = rhoML * wL;
            NekDouble EML    = rhoML * (EL / rhoL +
                                        (SM - uL) * (SM + pL / (rhoL * (SL - uL))));
            
            NekDouble rhoMR  = rhoR * (SR - uR) / (SR - SM);
            NekDouble rhouMR = rhoMR * SM;
            NekDouble rhovMR = rhoMR * vR;
            NekDouble rhowMR = rhoMR * wR;
            NekDouble EMR    = rhoMR * (ER / rhoR +
                                        (SM - uR) * (SM + pR / (rhoR * (SR - uR))));
            
            if (SL < 0.0 && SM >= 0.0)
            {
                rhof  = rhouL + SL * (rhoML - rhoL);
                rhouf = rhouL * uL + pL + SL * (rhouML - rhouL);
                rhovf = rhouL * vL + SL * (rhovML - rhovL);
                rhowf = rhouL * wL + SL * (rhowML - rhowL);
                Ef    = uL * (EL + pL) + SL * (EML - EL);
            }
            else if(SM < 0.0 && SR > 0.0)
            {
                rhof  = rhouR + SR * (rhoMR - rhoR);
                rhouf = rhouR * uR + pR + SR * (rhouMR - rhouR);
                rhovf = rhouR * vR + SR * (rhovMR - rhovR);
                rhowf = rhouR * wR + SR * (rhowMR - rhowR);
                Ef    = uR * (ER + pR) + SR * (EMR - ER);
            }
        }
    }

    void HLLCSolver::v_PointSolveVisc(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL, NekDouble  EpsL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER, NekDouble  EpsR,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef, NekDouble &Epsf)
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
        
        // HLLC Riemann fluxes (positive case)
        if (SL >= 0)
        {
            rhof  = rhouL;
            rhouf = rhouL * uL + pL;
            rhovf = rhouL * vL;
            rhowf = rhouL * wL;
            Ef    = uL * (EL + pL);
            Epsf  = 0.0;
        }
        // HLLC Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            rhof  = rhouR;
            rhouf = rhouR * uR + pR;
            rhovf = rhouR * vR;
            rhowf = rhouR * wR;
            Ef    = uR * (ER + pR);
            Epsf  = 0.0;
        }
        // HLLC Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            NekDouble SM = (pR - pL + rhouL * (SL - uL) - rhouR * (SR - uR)) /
            (rhoL * (SL - uL) - rhoR * (SR - uR));
            NekDouble rhoML  = rhoL * (SL - uL) / (SL - SM);
            NekDouble rhouML = rhoML * SM;
            NekDouble rhovML = rhoML * vL;
            NekDouble rhowML = rhoML * wL;
            NekDouble EML    = rhoML * (EL / rhoL +
                                        (SM - uL) * (SM + pL / (rhoL * (SL - uL))));
            NekDouble EpsML  = EpsL * (SL - uL) / (SL - SM);
            
            NekDouble rhoMR  = rhoR * (SR - uR) / (SR - SM);
            NekDouble rhouMR = rhoMR * SM;
            NekDouble rhovMR = rhoMR * vR;
            NekDouble rhowMR = rhoMR * wR;
            NekDouble EMR    = rhoMR * (ER / rhoR +
                                        (SM - uR) * (SM + pR / (rhoR * (SR - uR))));
            NekDouble EpsMR    = EpsR * (SL - uR) / (SL - SM);
            
            if (SL < 0.0 && SM >= 0.0)
            {
                rhof  = rhouL + SL * (rhoML - rhoL);
                rhouf = rhouL * uL + pL + SL * (rhouML - rhouL);
                rhovf = rhouL * vL + SL * (rhovML - rhovL);
                rhowf = rhouL * wL + SL * (rhowML - rhowL);
                Ef    = uL * (EL + pL) + SL * (EML - EL);
                Epsf  = 0.0 + SL * (EpsML - EpsL);
            }
            else if(SM < 0.0 && SR > 0.0)
            {
                rhof  = rhouR + SR * (rhoMR - rhoR);
                rhouf = rhouR * uR + pR + SR * (rhouMR - rhouR);
                rhovf = rhouR * vR + SR * (rhovMR - rhovR);
                rhowf = rhouR * wR + SR * (rhowMR - rhowR);
                Ef    = uR * (ER + pR) + SR * (EMR - ER);
                Epsf  = 0.0 + SR * (EpsMR - EpsR);
            }
        }
    }
}
