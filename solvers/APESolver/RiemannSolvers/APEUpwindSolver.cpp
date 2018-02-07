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
// Description: Upwind Riemann solver for the APE equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <APESolver/RiemannSolvers/APEUpwindSolver.h>

using namespace std;

namespace Nektar
{

std::string APEUpwindSolver::solverName = SolverUtils::GetRiemannSolverFactory().
                                       RegisterCreatorFunction("APEUpwind", APEUpwindSolver::create,
                                               "Upwind solver for the APE equation");

APEUpwindSolver::APEUpwindSolver(const LibUtilities::SessionReaderSharedPtr& pSession) :
    APESolver(pSession)
{

}


/**
 * @brief Upwind Riemann solver
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
void APEUpwindSolver::v_PointSolve(
    NekDouble  pL,  NekDouble  rhoL,  NekDouble  uL,  NekDouble  vL,  NekDouble  wL,
    NekDouble  pR,  NekDouble  rhoR,  NekDouble  uR,  NekDouble  vR,  NekDouble  wR,
    NekDouble  p0L, NekDouble  rho0L, NekDouble  u0L, NekDouble  v0L, NekDouble  w0L,
    NekDouble  p0R, NekDouble  rho0R, NekDouble  u0R, NekDouble  v0R, NekDouble  w0R,
    NekDouble &pF,  NekDouble &rhoF,  NekDouble &uF,  NekDouble &vF,  NekDouble &wF)
{
    // fetch params
    ASSERTL1(CheckParams("Gamma"), "Gamma not defined.");
    const NekDouble &gamma = m_params["Gamma"]();

    // Speed of sound
    NekDouble cL = sqrt(gamma * p0L / rho0L);
    NekDouble cR = sqrt(gamma * p0R / rho0R);

    Array<OneD, NekDouble> characteristic(4);
    Array<OneD, NekDouble> W(2);
    Array<OneD, NekDouble> lambda(2);

    // compute the wave speeds
    lambda[0] = (u0L + u0R) / 2 + (cL + cR) / 2;
    lambda[1] = (u0L + u0R) / 2 - (cL + cR) / 2;

    // calculate the caracteristic variables
    // left characteristics
    characteristic[0] = pL / 2 + uL * cL * rho0L / 2;
    characteristic[1] = pL / 2 - uL * cL * rho0L / 2;
    // right characteristics
    characteristic[2] = pR / 2 + uR * cR * rho0R / 2;
    characteristic[3] = pR / 2 - uR * cR * rho0R / 2;

    // take left or right value of characteristic variable
    for (int j = 0; j < 2; j++)
    {
        if (lambda[j] >= 0)
        {
            W[j] = characteristic[j];
        }
        else
        {
            W[j] = characteristic[j + 2];
        }
    }

    // calculate conservative variables from characteristics
    NekDouble p = W[0] + W[1];
    NekDouble u = (W[0] - W[1]) / (cL * rho0L); // TODO

    // assemble the fluxes
    pF = gamma * p0L * u + u0L * p; // TODO
    uF =
        p / rho0L + u0L * u + v0L * (vL + vR) / 2 + w0L * (wL + wR) / 2; // TODO
    vF = 0.0;
    wF = 0.0;
}

}


