///////////////////////////////////////////////////////////////////////////////
//
// File: UpwindPulseSolver.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Upwind pulse Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/RiemannSolvers/UpwindPulseSolver.h>

namespace Nektar
{
std::string UpwindPulseSolver::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "UpwindPulse", UpwindPulseSolver::create, "UpwindPulseSolver");

UpwindPulseSolver::UpwindPulseSolver(
    const LibUtilities::SessionReaderSharedPtr& pSession)
    : RiemannSolver(pSession)
{
}

/**
 *  Calculates the third term of the weak form (1): numerical flux
 *  at boundary \f$ \left[ \mathbf{\psi}^{\delta} \cdot \{
 *  \mathbf{F}^u - \mathbf{F}(\mathbf{U}^{\delta}) \}
 *  \right]_{x_e^l}^{x_eû} \f$
 */
void UpwindPulseSolver::v_Solve(
    const int nDim, const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    int i;
    int nTracePts = Fwd[0].size();

    ASSERTL1(CheckScalars("A0"), "A0 not defined.");
    const Array<OneD, NekDouble> &A0 = m_scalars["A0"]();

    ASSERTL1(CheckScalars("beta"), "beta not defined.");
    const Array<OneD, NekDouble> &beta = m_scalars["beta"]();

    ASSERTL1(CheckScalars("N"), "N not defined.");
    const Array<OneD, NekDouble> &N = m_scalars["N"]();

    for (i = 0; i < nTracePts; ++i)
    {
        RiemannSolverUpwind(Fwd[0][i], Fwd[1][i], Bwd[0][i], Bwd[1][i],
                            flux[0][i], flux[1][i], A0[i], beta[i], N[i]);
    }
}

/**
 *  Riemann solver for upwinding at an interface between two
 *  elements. Uses the characteristic variables for calculating
 *  the upwinded state \f$(A_u,u_u)\f$ from the left
 *  \f$(A_L,u_L)\f$ and right state \f$(A_R,u_R)\f$.  Returns the
 *  upwinded flux $\mathbf{F}^u$ needed for the weak formulation
 *  (1). Details can be found in "Pulse wave propagation in the
 *  human vascular system", section 3.3
 *
 */
void UpwindPulseSolver::RiemannSolverUpwind(NekDouble AL, NekDouble uL,
                                            NekDouble AR, NekDouble uR,
                                            NekDouble &Aflux, NekDouble &uflux,
                                            NekDouble A_0, NekDouble beta,
                                            NekDouble n)
{
    Array<OneD, NekDouble> W(2);
    Array<OneD, NekDouble> upwindedphysfield(2);
    NekDouble cL  = 0.0;
    NekDouble cR  = 0.0;
    NekDouble p   = 0.0;
    NekDouble p_t = 0.0;

    ASSERTL1(CheckParams("rho"), "rho not defined.");
    NekDouble rho = m_params["rho"]();

    ASSERTL1(CheckParams("pext"), "pext not defined.");
    NekDouble pext = m_params["pext"]();

    // Compute the wave speeds. The use of the normal here allows
    // for the definition of the characteristics to be inverted
    // (and hence the left and right state) if n is in the -ve
    // x-direction. This means we end up with the positive
    // defintion of the flux which has to therefore be multiplied
    // by the normal at the end of the methods This is a bit of a
    // mind twister but is efficient from a coding perspective.
    cL = sqrt(beta * sqrt(AL) / (2 * rho)) * n;
    cR = sqrt(beta * sqrt(AR) / (2 * rho)) * n;

    ASSERTL1(fabs(cL + cR) > fabs(uL + uR), "Conditions are not sub-sonic");

    // If upwinding from left and right for subsonic domain
    // then know characteristics immediately
    W[0] = uL + 4 * cL;
    W[1] = uR - 4 * cR;

    // Calculate conservative variables from characteristics
    NekDouble w0mw1 = 0.25 * (W[0] - W[1]);
    NekDouble fac   = rho / (2 * beta);
    w0mw1 *= w0mw1; // squared
    w0mw1 *= w0mw1; // fourth power
    fac *= fac;     // squared
    upwindedphysfield[0] = w0mw1 * fac;
    upwindedphysfield[1] = 0.5 * (W[0] + W[1]);

    // Compute the fluxes multipled by the normal.
    Aflux = upwindedphysfield[0] * upwindedphysfield[1] * n;
    p     = pext + beta * (sqrt(upwindedphysfield[0]) - sqrt(A_0));
    p_t   = 0.5 * (upwindedphysfield[1] * upwindedphysfield[1]) + p / rho;
    uflux = p_t * n;
}
}
