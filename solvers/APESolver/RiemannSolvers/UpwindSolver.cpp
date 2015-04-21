///////////////////////////////////////////////////////////////////////////////
//
// File: UpwindSolver.cpp
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

#include <APESolver/RiemannSolvers/UpwindSolver.h>

using namespace std;

namespace Nektar
{

std::string UpwindSolver::solverName = SolverUtils::GetRiemannSolverFactory().
                                       RegisterCreatorFunction("APEUpwind", UpwindSolver::create,
                                               "Upwind solver for the APE equation");

UpwindSolver::UpwindSolver() :
    APESolver()
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
void UpwindSolver::v_PointSolve(
    NekDouble  pL, NekDouble  rhoL, NekDouble  uL, NekDouble  vL, NekDouble  wL,
    NekDouble  pR, NekDouble  rhoR, NekDouble  uR, NekDouble  vR, NekDouble  wR,
    NekDouble  p0, NekDouble  rho0, NekDouble  u0, NekDouble  v0, NekDouble  w0,
    NekDouble &pF, NekDouble &rhoF, NekDouble &uF, NekDouble &vF, NekDouble &wF)
{
    // fetch params
    ASSERTL1(CheckParams("Gamma"), "Gamma not defined.");
    const NekDouble &gamma = m_params["Gamma"]();

    // Speed of sound
    NekDouble c = sqrt(gamma * p0 / rho0);

    Array<OneD, NekDouble> characteristic(4);
    Array<OneD, NekDouble> W(2);
    Array<OneD, NekDouble> lambda(2);

    // compute the wave speeds
    lambda[0] = u0 + c;
    lambda[1] = u0 - c;

    // calculate the caracteristic variables
    //left characteristics
    characteristic[0] = pL/2 + uL*c*rho0/2;
    characteristic[1] = pL/2 - uL*c*rho0/2;
    //right characteristics
    characteristic[2] = pR/2 + uR*c*rho0/2;
    characteristic[3] = pR/2 - uR*c*rho0/2;

    //take left or right value of characteristic variable
    for (int j = 0; j < 2; j++)
    {
        if (lambda[j] >= 0)
        {
            W[j] = characteristic[j];
        }
        if (lambda[j] < 0)
        {
            W[j] = characteristic[j+2];
        }
    }

    //calculate conservative variables from characteristics
    NekDouble p = W[0] + W[1];
    NekDouble u = (W[0] - W[1])/(c*rho0);

    // assemble the fluxes
    pF = rho0*u + u0*p/(c*c);
    uF = p/rho0 + u0*u + v0*vL + w0*wL;
    vF = 0.0;
    wF = 0.0;
}

}

