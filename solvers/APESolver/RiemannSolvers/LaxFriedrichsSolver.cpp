///////////////////////////////////////////////////////////////////////////////
//
// File: LaxFriedrichsSolver.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2014 Kilian Lackhove
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

#include <APESolver/RiemannSolvers/LaxFriedrichsSolver.h>

using namespace std;

namespace Nektar
{

std::string LaxFriedrichsSolver::solverName = SolverUtils::GetRiemannSolverFactory().
        RegisterCreatorFunction("LaxFriedrichs", LaxFriedrichsSolver::create, "Lax-Friedrichs Solver");


/**
*
*/
LaxFriedrichsSolver::LaxFriedrichsSolver() :
    APESolver()
{

}

/**
 * @brief Lax-Friedrichs Riemann solver
 *
 * @param pL     Perturbation pressure left state.
 * @param pR     Perturbation pressure right state.
 * @param uL     x perturbation verlocity component left state.
 * @param uR     x perturbation verlocity component right state.
 * @param vL     y perturbation verlocity component left state.
 * @param vR     y perturbation verlocity component right state.
 * @param wL     z perturbation verlocity component left state.
 * @param wR     z perturbation verlocity component right state.
 * @param p0     Base pressure.
 * @param u0     Base x verlocity component
 * @param v0     Base y verlocity component
 * @param w0     Base z verlocity component
 * @param pF     Computed Riemann flux for perturbation pressure.
 * @param uF     Computed Riemann flux for x perturbation verlocity component
 * @param vF     Computed Riemann flux for y perturbation verlocity component
 * @param wF     Computed Riemann flux for z perturbation verlocity component
 */
void LaxFriedrichsSolver::v_PointSolve(
    NekDouble  pL, NekDouble  uL, NekDouble  vL, NekDouble  wL,
    NekDouble  pR, NekDouble  uR, NekDouble  vR, NekDouble  wR,
    NekDouble  p0, NekDouble  u0, NekDouble  v0, NekDouble  w0,
    NekDouble &pF, NekDouble &uF, NekDouble &vF, NekDouble &wF)
{
    ASSERTL1(CheckParams("Gamma"), "Gamma not defined.");
    ASSERTL1(CheckParams("Rho"), "Rho not defined.");
    const NekDouble &gamma= m_params["Gamma"]();
    const NekDouble &rho = m_params["Rho"]();

    // Speed of sound
    NekDouble c = sqrt(gamma*p0/rho);

    // max absolute eigenvalue of the jacobian of F_n1
    NekDouble a_1_max = std::max(std::abs(u0 - c), std::abs(u0 + c));

    NekDouble pFL = gamma*p0*uL + pL*u0;
    NekDouble uFL = pL/rho + u0*uL + v0*vL + w0*wL;
    NekDouble vFL = 0;
    NekDouble wFL = 0;

    NekDouble pFR = gamma*p0*uR + pR*u0;
    NekDouble uFR = pR/rho + u0*uR + v0*vR + w0*wR;
    NekDouble vFR = 0;
    NekDouble wFR = 0;

    // assemble the face-normal fluxes
    pF = 0.5*(pFL + pFR - a_1_max*(pR - pL));
    uF = 0.5*(uFL + uFR - a_1_max*(uR - uL));
    vF = 0.5*(vFL + vFR - a_1_max*(vR - vL));
    wF = 0.5*(wFL + wFR - a_1_max*(wR - wL));
}

}
