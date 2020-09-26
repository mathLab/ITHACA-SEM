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
    : RiemannSolver(pSession), m_session(pSession)
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

    const Array<OneD, NekDouble> &alpha = m_scalars["alpha"]();

    ASSERTL1(CheckScalars("N"), "N not defined.");
    const Array<OneD, NekDouble> &N = m_scalars["N"]();

    for (i = 0; i < nTracePts; ++i)
    {
        RiemannSolverUpwind(Fwd[0][i], Fwd[1][i], Bwd[0][i], Bwd[1][i],
                            flux[0][i], flux[1][i], A0[i], beta[i], N[i], alpha[i]);
    }
}

/**
 *  Riemann solver for upwinding at an interface between two
 *  elements. Uses the characteristic variables for calculating
 *  the upwinded state \f$(A_u, u_u)\f$ from the left
 *  \f$(A_L, u_L)\f$ and right state \f$(A_R, u_R)\f$.  Returns the
 *  upwinded flux $\mathbf{F}^u$ needed for the weak formulation
 *  (1). Details can be found in "Pulse wave propagation in the
 *  human vascular system", section 3.3
 *
 */
void UpwindPulseSolver::RiemannSolverUpwind(NekDouble AL, NekDouble uL,
                                            NekDouble AR, NekDouble uR,
                                            NekDouble &Aflux, NekDouble &uflux,
                                            NekDouble A0, NekDouble beta,
                                            NekDouble n, NekDouble alpha)
{
    NekDouble W1 = 0.0;
    NekDouble W2 = 0.0;
    NekDouble IL = 0.0;
    NekDouble IR = 0.0;
    NekDouble Au = 0.0;
    NekDouble uu = 0.0;
    NekDouble cL = 0.0;
    NekDouble cR = 0.0;
    NekDouble P  = 0.0;

    ASSERTL1(CheckParams("rho"), "rho not defined.");
    NekDouble rho      = m_params["rho"]();
    NekDouble nDomains = m_params["domains"]();

    m_nVariables = m_session->GetVariables().size();

    m_vessels =
        Array<OneD, MultiRegions::ExpListSharedPtr>(m_nVariables * nDomains);

    if (m_session->DefinesSolverInfo("PressureArea"))
    {
        m_pressureArea = GetPressureAreaFactory().CreateInstance(
                m_session->GetSolverInfo("PressureArea"), m_vessels, m_session);
    }
    else
    {
        m_pressureArea = GetPressureAreaFactory().CreateInstance("Beta",
                                                          m_vessels, m_session);
    }

    // Compute the wave speeds to check dynamics are sub-sonic
    m_pressureArea->GetC(cL, beta, AL, A0, alpha);
    m_pressureArea->GetC(cR, beta, AR, A0, alpha);
    ASSERTL1(fabs(cL + cR) > fabs(uL + uR), "Conditions are not sub-sonic");

    /*
    Calculate the characteristics. The use of the normal here allows
    for the definition of the characteristics (and hence the left
    and right state) to be inverted if n is in the -ve
    x-direction. This means we end up with the positive
    defintion of the flux which has to therefore be multiplied
    by the normal at the end of the method. This is a bit of a
    mind twister but is efficient from a coding perspective.
    */
    m_pressureArea->GetCharIntegral(IL, beta, AL, A0, alpha);
    m_pressureArea->GetCharIntegral(IR, beta, AR, A0, alpha);
    W1 = uL + n * IL;
    W2 = uR - n * IR;

    // Calculate conservative variables from characteristics
    m_pressureArea->GetAFromChars(Au, n * W1, n * W2, beta, A0, alpha);
    m_pressureArea->GetUFromChars(uu, W1, W2);

    // Pressure for the energy flux
    m_pressureArea->GetPressure(P, beta, Au, A0, 0, 0, alpha);

    // Compute the fluxes multiplied by the normal
    Aflux = Au * uu * n;
    uflux = (uu * uu / 2 + P / rho) * n;
}

} // namespace Nektar
