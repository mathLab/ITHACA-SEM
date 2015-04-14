///////////////////////////////////////////////////////////////////////////////
//
// File: ExactSolver.cpp
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
// Description: Exact Riemann solver (Toro 2009).
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/ExactSolverToro.h>

#define TOL 1e-6
#define NRITER 20

namespace Nektar
{
    std::string ExactSolverToro::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "ExactToro",
            ExactSolverToro::create,
            "Exact Riemann solver");

    ExactSolverToro::ExactSolverToro() : CompressibleSolver()
    {

    }
    
    /**
     * @brief Use either PVRS, two-rarefaction or two-shock Riemann solvers to
     * calculate an initial pressure for the Newton-Raphson scheme.
     *
     * @param g      Array of calculated gamma values.
     * @param rhoL   Density left state.
     * @param rhoR   Density right state.  
     * @param uL     x-velocity component left state.
     * @param uR     x-velocity component right state.
     * @param pL     Pressure component left state.  
     * @param pR     Pressure component right state.  
     * @param cL     Sound speed component left state.  
     * @param cR     Sound speed component right state.
     * @return       Computed initial guess for the Newton-Raphson scheme.
     */
    inline NekDouble guessp(
        NekDouble g[],
        NekDouble rhoL, NekDouble uL, NekDouble pL, NekDouble cL,
        NekDouble rhoR, NekDouble uR, NekDouble pR, NekDouble cR)
    {
        const NekDouble quser = 2.0;
        NekDouble cup, ppv, pmin, pmax, qmax;
        
        cup  = 0.25*(rhoL + rhoR)*(cL + cR);
        ppv  = 0.5 *(pL + pR) + 0.5*(uL - uR)*cup;
        ppv  = std::max(0.0, ppv);
        pmin = std::min(pL, pR);
        pmax = std::max(pL, pR);
        qmax = pmax/pmin;
        
        if (qmax <= quser && pmin <= ppv && ppv <= pmax)
        {
            // Select PVRS Riemann solver.
            return ppv;
        }
        else if (ppv < pmin)
        {
            // Select two-rarefaction Riemann solver.
            NekDouble pq = pow(pL/pR, g[1]);
            NekDouble um = (pq*uL/cL + uR/cR + g[4]*(pq - 1.0))/(pq/cL + 1.0/cR);
            NekDouble ptL = 1.0 + g[7]*(uL - um)/cL;
            NekDouble ptR = 1.0 + g[7]*(um - uR)/cR;
            return 0.5*(pL*pow(ptL, g[3]) + pR*pow(ptR, g[3]));
        }
        else
        {
            // Select two-shock Riemann solver with PVRS as estimate.
            NekDouble geL = sqrt((g[5]/rhoL)/(g[6]*pL + ppv));
            NekDouble geR = sqrt((g[5]/rhoR)/(g[6]*pR + ppv));
            return (geL*pL + geR*pR - (uR - uL))/(geL + geR);
        }
    }
    
    /**
     * @brief Evaluate pressure functions fL and fR in Newton iteration of
     * Riemann solver (see equation 4.85 and 4.86 of Toro 2009).
     *
     * @param g   Array of gamma parameters.
     * @param p   Pressure at current iteration.
     * @param dk  Density (left or right state)
     * @param pk  Pressure (left or right state)
     * @param ck  Sound speed (left or right state)
     * @param f   Computed pressure (shock).
     * @param fd  Computed pressure (rarefaction).
     */
    inline void prefun(NekDouble *g, NekDouble  p, NekDouble  dk, NekDouble pk,
                       NekDouble ck, NekDouble &f, NekDouble &fd)
    {
        if (p <= pk)
        {
            // Rarefaction wave
            NekDouble prat = p/pk;
            f = g[4]*ck*(pow(prat, g[1])-1.0);
            fd = pow(prat, -g[2])/(dk*ck);
        }
        else
        {
            // Shock wave
            NekDouble ak  = g[5]/dk;
            NekDouble bk  = g[6]*pk;
            NekDouble qrt = sqrt(ak/(bk+p));
            f             = (p-pk)*qrt;
            fd            = (1.0-0.5*(p-pk)/(bk+p))*qrt;
        }
    }

    /**
     * @brief Exact Riemann solver for the Euler equations.
     *
     * This algorithm is transcribed from:
     * 
     *   "Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical
     *   Introduction", E. F. Toro (3rd edition, 2009).
     *
     * The full Fortran 77 routine can be found at the end of chapter 4 (section
     * 4.9). This transcription is essentially the functions STARPU and SAMPLE
     * glued together, and variable names are kept mostly the same. See the
     * preceding chapter which explains the derivation of the solver. The
     * routines PREFUN and GUESSP are kept separate and are reproduced above.
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
    void ExactSolverToro::v_PointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef)
    {
        static NekDouble gamma = m_params["gamma"]();

        // Left and right variables.
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;
        
        // Left and right pressure.
        NekDouble pL = (gamma - 1.0) * 
                       (EL - 0.5 * (rhouL*uL + rhovL*vL + rhowL*wL));
        NekDouble pR = (gamma - 1.0) * 
                       (ER - 0.5 * (rhouR*uR + rhovR*vR + rhowR*wR));
        
        // Compute gammas.
        NekDouble g[] = { gamma,
                          (gamma-1.0)/(2.0*gamma),
                          (gamma+1.0)/(2.0*gamma),
                          2.0*gamma/(gamma-1.0),
                          2.0/(gamma-1.0),
                          2.0/(gamma+1.0),
                          (gamma-1.0)/(gamma+1.0),
                          0.5*(gamma-1.0),
                          gamma-1.0 };

        // Compute sound speeds.
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);

        ASSERTL0(g[4]*(cL+cR) > (uR-uL), "Vacuum is generated by given data.");

        // Guess initial pressure.
        NekDouble pOld  = guessp(g, rhoL, uL, pL, cL, rhoR, uR, pR, cR);
        NekDouble uDiff = uR - uL;
        NekDouble p, fL, fR, fLd, fRd, change;
        int k;

        // Newton-Raphson iteration for pressure in star region.
        for (k = 0; k < NRITER; ++k)
        {
            prefun(g, pOld, rhoL, pL, cL, fL, fLd);
            prefun(g, pOld, rhoR, pR, cR, fR, fRd);
            p = pOld - (fL+fR+uDiff) / (fLd+fRd);
            change = 2 * fabs((p-pOld)/(p+pOld));
            
            if (change <= TOL)
            {
                break;
            }
            
            if (p < 0.0)
            {
                p = TOL;
            }
            
            pOld = p;
        }
        
        ASSERTL0(k < NRITER, "Divergence in Newton-Raphson scheme");

        // Compute velocity in star region.
        NekDouble u = 0.5*(uL+uR+fR-fL);
        
        // -- SAMPLE ROUTINE --
        // The variable S aligns with the Fortran parameter S of the SAMPLE
        // routine, but is hard-coded as 0 (and should be optimised out by the
        // compiler). Since we are using a Godunov scheme we pick this as 0 (see
        // chapter 6 of Toro 2009).
        const NekDouble S = 0.0;

        // Computed primitive variables.
        NekDouble outRho, outU, outV, outW, outP;
        
        if (S <= u)
        {
            if (p <= pL)
            {
                // Left rarefaction
                NekDouble shL = uL - cL;
                if (S <= shL)
                {
                    // Sampled point is left data state.
                    outRho = rhoL;
                    outU   = uL;
                    outV   = vL;
                    outW   = wL;
                    outP   = pL;
                }
                else
                {
                    NekDouble cmL = cL*pow(p/pL, g[1]);
                    NekDouble stL = u - cmL;
                    
                    if (S > stL)
                    {
                        // Sampled point is star left state
                        outRho = rhoL*pow(p/pL, 1.0/gamma);
                        outU   = u;
                        outV   = vL;
                        outW   = wL;
                        outP   = p;
                    }
                    else
                    {
                        // Sampled point is inside left fan
                        NekDouble c = g[5]*(cL + g[7]*(uL - S));
                        outRho = rhoL*pow(c/cL, g[4]);
                        outU   = g[5]*(cL + g[7]*uL + S);
                        outV   = vL;
                        outW   = wL;
                        outP   = pL*pow(c/cL, g[3]);
                    }
                }
            }
            else
            {
                // Left shock
                NekDouble pmL = p / pL;
                NekDouble SL  = uL - cL*sqrt(g[2]*pmL + g[1]);
                if (S <= SL)
                {
                    // Sampled point is left data state
                    outRho = rhoL;
                    outU   = uL;
                    outV   = vL;
                    outW   = wL;
                    outP   = pL;
                }
                else
                {
                    // Sampled point is star left state
                    outRho = rhoL*(pmL + g[6])/(pmL*g[6] + 1.0);
                    outU   = u;
                    outV   = vL;
                    outW   = wL;
                    outP   = p;
                }
            }
        }
        else
        {
            if (p > pR)
            {
                // Right shock
                NekDouble pmR = p/pR;
                NekDouble SR = uR + cR*sqrt(g[2]*pmR + g[1]);
                if (S >= SR)
                {
                    // Sampled point is right data state
                    outRho = rhoR;
                    outU   = uR;
                    outV   = vR;
                    outW   = wR;
                    outP   = pR;
                }
                else
                {
                    // Sampled point is star right state
                    outRho = rhoR*(pmR + g[6])/(pmR*g[6] + 1.0);
                    outU   = u;
                    outV   = vR;
                    outW   = wR;
                    outP   = p;
                }
            }
            else
            {
                // Right rarefaction
                NekDouble shR = uR + cR;
                
                if (S >= shR)
                {
                    // Sampled point is right data state
                    outRho = rhoR;
                    outU   = uR;
                    outV   = vR;
                    outW   = wR;
                    outP   = pR;
                }
                else
                {
                    NekDouble cmR = cR*pow(p/pR, g[1]);
                    NekDouble stR = u + cmR;
                    
                    if (S <= stR)
                    {
                        // Sampled point is star right state
                        outRho = rhoR*pow(p/pR, 1.0/gamma);
                        outU   = u;
                        outV   = vR;
                        outW   = wR;
                        outP   = p;
                    }
                    else
                    {
                        // Sampled point is inside left fan
                        NekDouble c = g[5]*(cR - g[7]*(uR-S));
                        outRho = rhoR*pow(c/cR, g[4]);
                        outU   = g[5]*(-cR + g[7]*uR + S);
                        outV   = vR;
                        outW   = wR;
                        outP   = pR*pow(c/cR, g[3]);
                    }
                }
            }
        }
        
        // Transform computed primitive variables to fluxes.
        rhof  = outRho * outU;
        rhouf = outP + outRho*outU*outU;
        rhovf = outRho * outU * outV;
        rhowf = outRho * outU * outW;
        Ef    = outU*(outP/(gamma-1.0) + 0.5*outRho*
                      (outU*outU + outV*outV + outW*outW) + outP);
    }
}
