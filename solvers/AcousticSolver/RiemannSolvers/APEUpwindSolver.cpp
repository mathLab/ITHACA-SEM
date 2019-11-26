///////////////////////////////////////////////////////////////////////////////
//
// File: APEUpwindSolver.cpp
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
// Description: Upwind Riemann solver for the APE equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <AcousticSolver/RiemannSolvers/APEUpwindSolver.h>

using namespace std;

namespace Nektar
{

std::string APEUpwindSolver::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "APEUpwind", APEUpwindSolver::create,
        "Upwind solver for the APE equation");

APEUpwindSolver::APEUpwindSolver(
    const LibUtilities::SessionReaderSharedPtr &pSession)
    : AcousticSolver(pSession)
{
}

/**
 * @brief Upwind Riemann solver
 *
 * The fluxes are straight out of sympy, so lets just hope the compiler
 * optimizes them for us.
 *
 * @param pL     Perturbation pressure left state
 * @param rhoL   Perturbation density left state
 * @param pR     Perturbation pressure right state
 * @param rhoR   Perturbation density right state
 * @param uL     x perturbation velocity component left state
 * @param uR     x perturbation velocity component right state
 * @param vL     y perturbation velocity component left state
 * @param vR     y perturbation velocity component right state
 * @param wL     z perturbation velocity component left state
 * @param wR     z perturbation velocity component right state
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
 * @param uF     Computed Riemann flux for x perturbation velocity component
 * @param vF     Computed Riemann flux for y perturbation velocity component
 * @param wF     Computed Riemann flux for z perturbation velocity component
 */
void APEUpwindSolver::v_PointSolve(
    NekDouble  pL,    NekDouble  rhoL,  NekDouble  uL,  NekDouble  vL,  NekDouble  wL,
    NekDouble  pR,    NekDouble  rhoR,  NekDouble  uR,  NekDouble  vR,  NekDouble  wR,
    NekDouble  c0sqL, NekDouble  rho0L, NekDouble  u0L, NekDouble  v0L, NekDouble  w0L,
    NekDouble  c0sqR, NekDouble  rho0R, NekDouble  u0R, NekDouble  v0R, NekDouble  w0R,
    NekDouble &pF,    NekDouble &rhoF,  NekDouble &uF,  NekDouble &vF,  NekDouble &wF)
{
    boost::ignore_unused(rhoL, rhoR, rhoF);

    // Speed of sound
    NekDouble c0L = sqrt(c0sqL);
    NekDouble c0R = sqrt(c0sqR);
    NekDouble c0M = (c0L + c0R) / 2.0;

    NekDouble u0M = (u0L + u0R) / 2.0;

    pF = 0.0;
    uF = 0.0;
    vF = 0.0;
    wF = 0.0;

    // lambda_3
    if (u0M - c0M > 0)
    {
        pF = pF + 0.5 * (c0L - u0L) *
                      (-rho0L * v0L * vL * (c0L * u0L + c0sqL) -
                       rho0L * w0L * wL * (c0L * u0L + c0sqL) +
                       (c0sqL - pow(u0L, 2)) * (rho0L * c0L * uL - pL)) /
                      (c0sqL - pow(u0L, 2));

        uF = uF +
             0.5 * (c0L - u0L) *
                 (rho0L * c0L * (c0L * u0L + c0sqL) * (v0L * vL + w0L * wL) -
                  rho0L * c0sqL * uL * (c0sqL - pow(u0L, 2)) +
                  c0L * pL * (c0sqL - pow(u0L, 2))) /
                 (rho0L * c0sqL * (c0sqL - pow(u0L, 2)));
    }
    else
    {
        pF = pF + 0.5 * (c0R - u0R) *
                      (-rho0R * v0R * vR * (c0R * u0R + c0sqR) -
                       rho0R * w0R * wR * (c0R * u0R + c0sqR) +
                       (c0sqR - pow(u0R, 2)) * (rho0R * c0R * uR - pR)) /
                      (c0sqR - pow(u0R, 2));

        uF = uF +
             0.5 * (c0R - u0R) *
                 (rho0R * c0R * (c0R * u0R + c0sqR) * (v0R * vR + w0R * wR) -
                  rho0R * c0sqR * uR * (c0sqR - pow(u0R, 2)) +
                  c0R * pR * (c0sqR - pow(u0R, 2))) /
                 (rho0R * c0sqR * (c0sqR - pow(u0R, 2)));
    }

    // lambda_4
    if (u0M + c0M > 0)
    {
        pF = pF + 0.5 * (c0L + u0L) *
                      (-rho0L * v0L * vL * (c0L * u0L - c0sqL) -
                       rho0L * w0L * wL * (c0L * u0L - c0sqL) +
                       (c0sqL - pow(u0L, 2)) * (rho0L * c0L * uL + pL)) /
                      (c0sqL - pow(u0L, 2));

        uF = uF +
             0.5 * (c0L + u0L) *
                 (-rho0L * c0L * (c0L * u0L - c0sqL) * (v0L * vL + w0L * wL) +
                  rho0L * c0sqL * uL * (c0sqL - pow(u0L, 2)) +
                  c0L * pL * (c0sqL - pow(u0L, 2))) /
                 (rho0L * c0sqL * (c0sqL - pow(u0L, 2)));
    }
    else
    {
        pF = pF + 0.5 * (c0R + u0R) *
                      (-rho0R * v0R * vR * (c0R * u0R - c0sqR) -
                       rho0R * w0R * wR * (c0R * u0R - c0sqR) +
                       (c0sqR - pow(u0R, 2)) * (rho0R * c0R * uR + pR)) /
                      (c0sqR - pow(u0R, 2));

        uF = uF +
             0.5 * (c0R + u0R) *
                 (-rho0R * c0R * (c0R * u0R - c0sqR) * (v0R * vR + w0R * wR) +
                  rho0R * c0sqR * uR * (c0sqR - pow(u0R, 2)) +
                  c0R * pR * (c0sqR - pow(u0R, 2))) /
                 (rho0R * c0sqR * (c0sqR - pow(u0R, 2)));
    }
}

} // namespace Nektar
