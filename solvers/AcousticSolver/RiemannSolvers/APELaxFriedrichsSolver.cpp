///////////////////////////////////////////////////////////////////////////////
//
// File: APELaxFriedrichsSolver.cpp
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
// Description: Lax-Friedrichs solver for the APE equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <AcousticSolver/RiemannSolvers/APELaxFriedrichsSolver.h>

using namespace std;

namespace Nektar
{

std::string APELaxFriedrichsSolver::solverName =
    SolverUtils::GetRiemannSolverFactory().
    RegisterCreatorFunction("APELaxFriedrichs", APELaxFriedrichsSolver::create,
                            "Lax-Friedrichs Solver");


/**
*
*/
APELaxFriedrichsSolver::APELaxFriedrichsSolver(const LibUtilities::SessionReaderSharedPtr& pSession) :
    AcousticSolver(pSession)
{

}

/**
 * @brief Lax-Friedrichs Riemann solver
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
 * @param p0     Base pressure
 * @param rho0   Base density
 * @param u0     Base x velocity component
 * @param v0     Base y velocity component
 * @param w0     Base z velocity component
 * @param pF     Computed Riemann flux for perturbation pressure
 * @param rhoF   Computed Riemann flux for perturbation density
 * @param uF     Computed Riemann flux for x perturbation velocity component
 * @param vF     Computed Riemann flux for y perturbation velocity component
 * @param wF     Computed Riemann flux for z perturbation velocity component
 */
void APELaxFriedrichsSolver::v_PointSolve(
    NekDouble  pL,  NekDouble  rhoL,  NekDouble  uL,  NekDouble  vL,  NekDouble  wL,
    NekDouble  pR,  NekDouble  rhoR,  NekDouble  uR,  NekDouble  vR,  NekDouble  wR,
    NekDouble  p0L, NekDouble  rho0L, NekDouble  u0L, NekDouble  v0L, NekDouble  w0L,
    NekDouble  p0R, NekDouble  rho0R, NekDouble  u0R, NekDouble  v0R, NekDouble  w0R,
    NekDouble &pF,  NekDouble &rhoF,  NekDouble &uF,  NekDouble &vF,  NekDouble &wF)
{
    ASSERTL1(CheckParams("Gamma"), "Gamma not defined.");
    const NekDouble &gamma = m_params["Gamma"]();

    // Speed of sound
    NekDouble cL = sqrt(gamma * p0L / rho0L);
    NekDouble cR = sqrt(gamma * p0R / rho0R);

    // max absolute eigenvalue of the jacobian of F_n1
    NekDouble a_1_max = 0;
    a_1_max = std::max(a_1_max, std::abs(u0L - cL));
    a_1_max = std::max(a_1_max, std::abs(u0R - cR));
    a_1_max = std::max(a_1_max, std::abs(u0L + cL));
    a_1_max = std::max(a_1_max, std::abs(u0R + cR));

    NekDouble pFL = gamma * p0L * uL + pL * u0L;
    NekDouble uFL = pL / rho0L + u0L * uL + v0L * vL + w0L * wL;
    NekDouble vFL = 0;
    NekDouble wFL = 0;

    NekDouble pFR = gamma * p0R * uR + pR * u0R;
    NekDouble uFR = pR / rho0R + u0R * uR + v0R * vR + w0R * wR;
    NekDouble vFR = 0;
    NekDouble wFR = 0;

    // assemble the face-normal fluxes
    pF = 0.5 * (pFL + pFR - a_1_max * (pR - pL));
    uF = 0.5 * (uFL + uFR - a_1_max * (uR - uL));
    vF = 0.5 * (vFL + vFR - a_1_max * (vR - vL));
    wF = 0.5 * (wFL + wFR - a_1_max * (wR - wL));
}

}

