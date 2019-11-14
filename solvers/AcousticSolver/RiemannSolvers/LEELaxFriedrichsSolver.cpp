///////////////////////////////////////////////////////////////////////////////
//
// File: LEELaxFriedrichsSolver.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
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

#include <AcousticSolver/RiemannSolvers/LEELaxFriedrichsSolver.h>

using namespace std;

namespace Nektar
{

std::string LEELaxFriedrichsSolver::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "LEELaxFriedrichs", LEELaxFriedrichsSolver::create,
        "Lax-Friedrichs Solver for LEE");

/**
 *
 */
LEELaxFriedrichsSolver::LEELaxFriedrichsSolver(
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
void LEELaxFriedrichsSolver::v_PointSolve(
    NekDouble  pL,    NekDouble  rhoL,  NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL,
    NekDouble  pR,    NekDouble  rhoR,  NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR,
    NekDouble  c0sqL, NekDouble  rho0L, NekDouble  u0L,   NekDouble  v0L,   NekDouble  w0L,
    NekDouble  c0sqR, NekDouble  rho0R, NekDouble  u0R,   NekDouble  v0R,   NekDouble  w0R,
    NekDouble &pF,    NekDouble &rhoF,  NekDouble &rhouF, NekDouble &rhovF, NekDouble &rhowF)
{
    boost::ignore_unused(rho0L, v0L, w0L, rho0R, v0R, w0R);

    // Speed of sound
    NekDouble cL = sqrt(c0sqL);
    NekDouble cR = sqrt(c0sqR);

    // max absolute eigenvalue of the jacobian of F_n1
    NekDouble a_1_max = 0;
    a_1_max           = std::max(a_1_max, std::abs(u0L - cL));
    a_1_max           = std::max(a_1_max, std::abs(u0R - cR));
    a_1_max           = std::max(a_1_max, std::abs(u0L + cL));
    a_1_max           = std::max(a_1_max, std::abs(u0R + cR));

    NekDouble pFL    = rhouL * c0sqL + u0L * pL;
    NekDouble rhoFL  = rhouL + rhoL * u0L;
    NekDouble rhouFL = pL + rhouL * u0L;
    NekDouble rhovFL = 0;
    NekDouble rhowFL = 0;

    NekDouble pFR    = rhouR * c0sqR + u0R * pR;
    NekDouble rhoFR  = rhouR + rhoR * u0R;
    NekDouble rhouFR = pR + rhouR * u0R;
    NekDouble rhovFR = 0;
    NekDouble rhowFR = 0;

    // assemble the face-normal fluxes
    pF    = 0.5 * (pFL + pFR - a_1_max * (pR - pL));
    rhoF  = 0.5 * (rhoFL + rhoFR - a_1_max * (rhoR - rhoL));
    rhouF = 0.5 * (rhouFL + rhouFR - a_1_max * (rhouR - rhouL));
    rhovF = 0.5 * (rhovFL + rhovFR - a_1_max * (rhovR - rhovL));
    rhowF = 0.5 * (rhowFL + rhowFR - a_1_max * (rhowR - rhowL));
}

} // namespace Nektar
