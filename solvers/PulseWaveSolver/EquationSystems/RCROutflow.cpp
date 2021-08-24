///////////////////////////////////////////////////////////////////////////////
//
// File RCROutflow.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/RCROutflow.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;

namespace Nektar
{

std::string RCROutflow::className =
    GetBoundaryFactory().RegisterCreatorFunction(
        "RCR-terminal", RCROutflow::create, "RCR outflow boundary condition");

RCROutflow::RCROutflow(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
                       const LibUtilities::SessionReaderSharedPtr pSession,
                       PulseWavePressureAreaSharedPtr pressureArea)
    : PulseWaveBoundary(pVessel, pSession, pressureArea)
{
    m_session->LoadParameter("TimeStep", m_timestep);
}

RCROutflow::~RCROutflow()
{
}

void RCROutflow::v_DoBoundary(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    Array<OneD, Array<OneD, NekDouble> > &A_0,
    Array<OneD, Array<OneD, NekDouble> > &beta,
    Array<OneD, Array<OneD, NekDouble> > &alpha, const NekDouble time, int omega,
    int offset, int n)
{
    NekDouble A_r  = 0.0;
    NekDouble u_r  = 0.0;
    NekDouble A_u  = 0.0;
    NekDouble u_u  = 0.0;
    NekDouble A_l  = 0.0;
    NekDouble u_l  = 0.0;
    NekDouble c_0  = 0.0;
    NekDouble R1   = 0.0;
    NekDouble R2   = 0.0;
    NekDouble POut = m_pout;
    NekDouble rho  = m_rho;

    Array<OneD, MultiRegions::ExpListSharedPtr> vessel(2);

    // Pointers to the domains
    vessel[0] = m_vessels[2 * omega];
    vessel[1] = m_vessels[2 * omega + 1];

    /* Find the terminal RCR boundary condition and calculates
       the updated velocity and area as well as the updated
       boundary conditions */

    /* Load terminal resistance, capacitance, outflow pressure,
       and number of points from the input file */
    NekDouble RT = ((vessel[0]->GetBndCondExpansions())[n])->GetCoeffs()[0];
    NekDouble C  = ((vessel[1]->GetBndCondExpansions())[n])->GetCoeffs()[0];
    int nq       = vessel[0]->GetTotPoints();

    // Get the values of all variables needed for the Riemann problem
    A_l = inarray[0][offset + nq - 1];
    u_l = inarray[1][offset + nq - 1];

    // Goes through the first resistance; calculates c_0
    m_pressureArea->GetC(c_0, beta[omega][nq - 1], A_0[omega][nq - 1], A_0[omega][nq - 1], alpha[omega][nq - 1]);

    /* Calculate R1 and R2, R1 being calculated so as
       to eliminate reflections in the vessel */
    R1 = rho * c_0 / A_0[omega][nq - 1];

    if (R1 > 0.9 * RT)
    {
        // In case the resistance is lower than the characteristic impedance.
        R1 = 0.9 * RT;
    }

    R2 = RT - R1;

    // Call the R RiemannSolver
    R_RiemannSolver(R1, A_l, u_l, A_0[omega][nq - 1], beta[omega][nq - 1],
                                          alpha[omega][nq - 1], m_pc, A_u, u_u);

    /* Fix the boundary conditions in the virtual region to ensure
       upwind state matches the boundary condition at the next time step */
    A_r = A_l;
    u_r = 2 * u_u - u_l;

    /* Goes through the CR system, which
       just updates the pressure pc */
    m_pc += (m_timestep / C) * (A_u * u_u - (m_pc - POut) / R2);

    // Store the updated values
    (vessel[0]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = A_r;
    (vessel[1]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = u_r;
}

void RCROutflow::R_RiemannSolver(NekDouble R, NekDouble A_l, NekDouble u_l,
                                 NekDouble A_0, NekDouble beta, NekDouble alpha,
                                 NekDouble POut, NekDouble &A_u, NekDouble &u_u)
{
    NekDouble W1           = 0.0;
    NekDouble c            = 0.0;
    NekDouble cL           = 0.0;
    NekDouble I            = 0.0;
    NekDouble A_calc       = 0.0;
    NekDouble FA           = 0.0;
    NekDouble dFdA         = 0.0;
    NekDouble delta_A_calc = 0.0;
    NekDouble P            = 0.0;
    NekDouble rho          = m_rho;

    int proceed  = 1;
    int iter     = 0;
    int MAX_ITER = 100;

    // Tolerances for the algorithm
    NekDouble Tol = 1.0E-10;

    // Calculate the wave speed
    m_pressureArea->GetC(cL, beta, A_l, A_0, alpha);

    // Riemann invariant \f$W_1(Al,ul)\f$
    m_pressureArea->GetW1(W1, u_l, beta, A_l, A_0, alpha);

    // Newton Iteration (Area only)
    A_calc = A_l;
    while ((proceed) && (iter < MAX_ITER))
    {
        iter += 1;

        m_pressureArea->GetPressure(P, beta, A_calc, A_0, 0, 0, alpha);
        m_pressureArea->GetC(c, beta, A_calc, A_0, alpha);
        m_pressureArea->GetCharIntegral(I, beta, A_calc, A_0, alpha);

        FA = R * A_calc * (W1 - I) - P + POut;
        dFdA = R * (W1 - I - c) - c * c * rho / A_calc;
        delta_A_calc = FA / dFdA;
        A_calc -= delta_A_calc;

        if (sqrt(delta_A_calc * delta_A_calc) < Tol)
        {
            proceed = 0;
        }
    }

    m_pressureArea->GetPressure(P, beta, A_calc, A_0, 0, 0, alpha);

    // Obtain u_u and A_u
    u_u = (P - POut) / (R * A_calc);
    A_u = A_calc;
}

} // namespace Nektar
