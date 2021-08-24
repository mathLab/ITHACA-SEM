///////////////////////////////////////////////////////////////////////////////
//
// File: LEEUpwindSolver.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Lax-Friedrichs solver for the LEE equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <AcousticSolver/RiemannSolvers/LEEUpwindSolver.h>

using namespace std;

namespace Nektar
{

std::string LEEUpwindSolver::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "LEEUpwind", LEEUpwindSolver::create, "Upwind Solver for LEE");

/**
 *
 */
LEEUpwindSolver::LEEUpwindSolver(
    const LibUtilities::SessionReaderSharedPtr &pSession)
    : LEESolver(pSession)
{
}

/**
 * @brief Lax-Friedrichs Riemann solver
 *
 * @param pL     Perturbation pressure left state
 * @param rhoL   Perturbation density left state
 * @param pR     Perturbation pressure right state
 * @param rhoR   Perturbation density right state
 * @param rhouL  x perturbation velocity component left state
 * @param rhouR  x perturbation velocity component right state
 * @param rhovL  y perturbation velocity component left state
 * @param rhovR  y perturbation velocity component right state
 * @param rhowL  z perturbation velocity component left state
 * @param rhowR  z perturbation velocity component right state
 * @param c0sqL  Base pressure left state
 * @param c0sqR  Base pressure right state
 * @param rho0L  Base density left state
 * @param rho0R  Base density right state
 * @param u0L    Base x velocity component left state
 * @param u0R    Base x velocity component right state
 * @param v0L    Base y velocity component left state
 * @param v0R    Base y velocity component right state
 * @param w0L    Base z velocity component left state
 * @param w0R    Base z velocity component right state
 * @param pF     Computed Riemann flux for perturbation pressure
 * @param rhoF   Computed Riemann flux for perturbation density
 * @param rhouF  Computed Riemann flux for x perturbation velocity component
 * @param rhovF  Computed Riemann flux for y perturbation velocity component
 * @param rhowF  Computed Riemann flux for z perturbation velocity component
 */
void LEEUpwindSolver::v_PointSolve(
    NekDouble  pL,    NekDouble  rhoL,  NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL,
    NekDouble  pR,    NekDouble  rhoR,  NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR,
    NekDouble  c0sqL, NekDouble  rho0L, NekDouble  u0L,   NekDouble  v0L,   NekDouble  w0L,
    NekDouble  c0sqR, NekDouble  rho0R, NekDouble  u0R,   NekDouble  v0R,   NekDouble  w0R,
    NekDouble &pF,    NekDouble &rhoF,  NekDouble &rhouF, NekDouble &rhovF, NekDouble &rhowF)
{
    boost::ignore_unused(rho0L, v0L, w0L, rho0R, v0R, w0R);

    // Speed of sound
    NekDouble c0L = sqrt(c0sqL);
    NekDouble c0R = sqrt(c0sqR);
    NekDouble c0M = (c0L + c0R) / 2.0;

    NekDouble u0M = (u0L + u0R) / 2.0;

    pF    = 0.0;
    rhoF  = 0.0;
    rhouF = 0.0;
    rhovF = 0.0;
    rhowF = 0.0;

    // lambda_1,2,3
    if (u0M > 0)
    {
        rhoF  = rhoF + u0L * (c0sqL * rhoL - pL) / c0sqL;
        rhovF = rhovF + rhovL * u0L;
        rhowF = rhowF + rhowL * u0L;
    }
    else
    {
        rhoF  = rhoF + u0R * (c0sqR * rhoR - pR) / c0sqR;
        rhovF = rhovF + rhovR * u0R;
        rhowF = rhowF + rhowR * u0R;
    }

    // lambda_4
    if (u0M - c0M > 0)
    {
        pF   = pF + 0.5 * (c0L - u0L) * (rhouL * c0L - pL);
        rhoF = rhoF + 0.5 * (c0L - u0L) * (rhouL * c0sqL - c0L * pL) /
                          pow(c0sqL, 3.0 / 2.0);
        rhouF = rhouF + 0.5 * (c0L - u0L) * (-rhouL * c0L + pL) / c0L;
    }
    else
    {
        pF   = pF + 0.5 * (c0R - u0R) * (rhouR * c0R - pR);
        rhoF = rhoF + 0.5 * (c0R - u0R) * (rhouR * c0sqR - c0R * pR) /
                          pow(c0sqR, 3.0 / 2.0);
        rhouF = rhouF + 0.5 * (c0R - u0R) * (-rhouR * c0R + pR) / c0R;
    }

    // lambda_5
    if (u0M + c0M > 0)
    {
        pF   = pF + 0.5 * (c0L + u0L) * (rhouL * c0L + pL);
        rhoF = rhoF + 0.5 * (c0L + u0L) * (rhouL * c0sqL + c0L * pL) /
                          pow(c0sqL, 3.0 / 2.0);
        rhouF = rhouF + 0.5 * (c0L + u0L) * (rhouL * c0L + pL) / c0L;
    }
    else
    {
        pF   = pF + 0.5 * (c0R + u0R) * (rhouR * c0R + pR);
        rhoF = rhoF + 0.5 * (c0R + u0R) * (rhouR * c0sqR + c0R * pR) /
                          pow(c0sqR, 3.0 / 2.0);
        rhouF = rhouF + 0.5 * (c0R + u0R) * (rhouR * c0R + pR) / c0R;
    }
}
} // namespace Nektar
