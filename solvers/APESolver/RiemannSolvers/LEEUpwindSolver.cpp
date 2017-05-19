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
// Description: Lax-Friedrichs solver for the LEE equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <APESolver/RiemannSolvers/LEEUpwindSolver.h>

using namespace std;

namespace Nektar
{

std::string LEEUpwindSolver::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "LEEUpwind",
        LEEUpwindSolver::create,
        "Upwind Solver for LEE");

/**
*
*/
LEEUpwindSolver::LEEUpwindSolver(const LibUtilities::SessionReaderSharedPtr& pSession) : LEESolver(pSession)
{
}

/**
 * @brief Lax-Friedrichs Riemann solver
 *
 * @param pL     Perturbation pressure left state
 * @param rhoL   Perturbation density left state
 * @param pR     Perturbation pressure right state
 * @param rhoR   Perturbation density right state
 * @param ruL    x perturbation velocity component left state
 * @param ruR    x perturbation velocity component right state
 * @param rvL     y perturbation velocity component left state
 * @param rvR     y perturbation velocity component right state
 * @param rwL    z perturbation velocity component left state
 * @param rwR    z perturbation velocity component right state
 * @param p0     Base pressure
 * @param rho0   Base density
 * @param u0     Base x velocity component
 * @param v0     Base y velocity component
 * @param w0     Base z velocity component
 * @param pF     Computed Riemann flux for perturbation pressure
 * @param rhoF   Computed Riemann flux for perturbation density
 * @param ruF    Computed Riemann flux for x perturbation velocity component
 * @param rvF    Computed Riemann flux for y perturbation velocity component
 * @param rwF    Computed Riemann flux for z perturbation velocity component
 */
void LEEUpwindSolver::v_PointSolve(
    NekDouble  pL,  NekDouble  rhoL,  NekDouble  ruL, NekDouble  rvL, NekDouble  rwL,
    NekDouble  pR,  NekDouble  rhoR,  NekDouble  ruR, NekDouble  rvR, NekDouble  rwR,
    NekDouble  p0L, NekDouble  rho0L, NekDouble  u0L, NekDouble  v0L, NekDouble  w0L,
    NekDouble  p0R, NekDouble  rho0R, NekDouble  u0R, NekDouble  v0R, NekDouble  w0R,
    NekDouble &pF,  NekDouble &rhoF,  NekDouble &ruF, NekDouble &rvF, NekDouble &rwF)
{
    ASSERTL1(CheckParams("Gamma"), "Gamma not defined.");
    const NekDouble &gamma = m_params["Gamma"]();

    // Speed of sound
    NekDouble cL = sqrt(gamma * p0L / rho0L);
    NekDouble cR = sqrt(gamma * p0R / rho0R);

    NekDouble cM  = (cL + cR) / 2;
    NekDouble u0M = (u0L + u0R) / 2;

    NekDouble h0, h1, h2;

    if (u0M > 0)
    {
        // rho - p / c^2
        h0 = rhoL - pL / (cM * cM);
    }
    else
    {
        h0 = rhoR - pR / (cM * cM);
    }

    if (u0M - cM > 0)
    {
        // ru / 2 - p / (2*c)
        h1 = ruL / 2 - pL / (2 * cM);
    }
    else
    {
        h1 = ruR / 2 - pR / (2 * cM);
    }

    if (u0M + cM > 0)
    {
        // ru / 2 + p / (2*c)
        h2 = ruL / 2 + pL / (2 * cM);
    }
    else
    {
        h2 = ruR / 2 + pR / (2 * cM);
    }

    // p = c0*(h2-h1)
    // rho = h0 + c0*(h2-h1)
    // ru = h1+h2
    NekDouble p   = cM * (h2 - h1);
    NekDouble rho = h0 + (h2 - h1) / cM;
    NekDouble ru  = h1 + h2;

    pF   = ru * cM * cM + u0M * p;
    rhoF = ru + rho * u0M;
    ruF  = p + ru * u0M;
    rvF  = (rvL + rvR) / 2 * u0M;
    rwF  = (rwL + rwR) / 2* u0M;

}
}
