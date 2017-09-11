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

#include <AcousticSolver/RiemannSolvers/LEEUpwindSolver.h>

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
 * @param ruF    Computed Riemann flux for x perturbation velocity component
 * @param rvF    Computed Riemann flux for y perturbation velocity component
 * @param rwF    Computed Riemann flux for z perturbation velocity component
 */
void LEEUpwindSolver::v_PointSolve(
    NekDouble  pL,    NekDouble  rhoL,  NekDouble  ruL, NekDouble  rvL, NekDouble  rwL,
    NekDouble  pR,    NekDouble  rhoR,  NekDouble  ruR, NekDouble  rvR, NekDouble  rwR,
    NekDouble  c0sqL, NekDouble  rho0L, NekDouble  u0L, NekDouble  v0L, NekDouble  w0L,
    NekDouble  c0sqR, NekDouble  rho0R, NekDouble  u0R, NekDouble  v0R, NekDouble  w0R,
    NekDouble &pF,    NekDouble &rhoF,  NekDouble &ruF, NekDouble &rvF, NekDouble &rwF)
{
    ASSERTL1(CheckParams("Gamma"), "Gamma not defined.");
    const NekDouble &gamma = m_params["Gamma"]();

    // Speed of sound
    NekDouble cL = sqrt(gamma * p0L / rho0L);
    NekDouble cR = sqrt(gamma * p0R / rho0R);

    NekDouble cM  = (cL + cR) / 2;
    NekDouble u0M = (u0L + u0R) / 2;

    NekDouble h1, h2, h3, h4, h5;

    if (u0M > 0)
    {
        h1 = rhoL - pL / pow(cL, 2);
        h2 = rvL;
        h3 = rwL;
    }
    else
    {
        h1 = rhoR - pR / pow(cR, 2);
        h2 = rvR;
        h3 = rwR;
    }

    if (u0M - cM > 0)
    {
        h4 = ruL / 2.0 - pL / (2 * cL);
    }
    else
    {
        h4 = ruR / 2.0 - pR / (2 * cR);
    }

    if (u0M + cM > 0)
    {
        h5 = ruL / 2.0 + pL / (2 * cL);
    }
    else
    {
        h5 = ruR / 2.0 + pR / (2 * cR);
    }

    NekDouble p   = cM * (h5 - h4);
    NekDouble rho = h1 - cM * (h4 + h5);
    NekDouble ru  = h4 + h5;
    NekDouble rv  = h2;
    NekDouble rw  = h3;

    pF   = ru * pow(cM, 2) + u0M * p;
    rhoF = ru + rho * u0M;
    ruF  = ru * u0M + p;
    rvF  = rv * u0M;
    rwF  = rw * u0M;
}
}
