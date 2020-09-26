///////////////////////////////////////////////////////////////////////////////
//
// File QInflow.cpp
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
// Description: QInflow class
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/QInflow.h>

using namespace std;

namespace Nektar
{

std::string QInflow::className = GetBoundaryFactory().RegisterCreatorFunction(
    "Q-inflow", QInflow::create, "Inflow boundary condition");

QInflow::QInflow(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
                 const LibUtilities::SessionReaderSharedPtr pSession,
                 PulseWavePressureAreaSharedPtr pressureArea)
    : PulseWaveBoundary(pVessel, pSession, pressureArea)
{
}

QInflow::~QInflow()
{
}

void QInflow::v_DoBoundary(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    Array<OneD, Array<OneD, NekDouble> > &A_0,
    Array<OneD, Array<OneD, NekDouble> > &beta,
    Array<OneD, Array<OneD, NekDouble> > &alpha,
    const NekDouble time, int omega, int offset, int n)
{
    NekDouble Q;
    NekDouble A_u;
    NekDouble u_u;
    NekDouble A_r;
    NekDouble u_r;
    NekDouble A_l;
    NekDouble u_l;

    Array<OneD, MultiRegions::ExpListSharedPtr> vessel(2);

    // Pointers to the domains
    vessel[0] = m_vessels[2 * omega];
    vessel[1] = m_vessels[2 * omega + 1];

    // Evaluates the time-dependent Q
    vessel[0]->EvaluateBoundaryConditions(time);

    // Q is contained as A in the input file
    Q = (vessel[0]->UpdateBndCondExpansion(n))->GetCoeffs()[0];

    // Initial conditions in the inlet vessel
    A_r = inarray[0][offset];
    u_r = inarray[1][offset];

    // Call the Q-inflow Riemann solver
    Q_RiemannSolver(Q, A_r, u_r, A_0[omega][0], beta[omega][0], alpha[omega][0], A_u, u_u);

    /* Fix the boundary conditions in the virtual region to ensure
       upwind state matches the boundary condition at the next time step */
    A_l = A_r;
    u_l = 2 * u_u - u_r;

    // Store the updated values in the boundary condition
    (vessel[0]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = A_l;
    (vessel[1]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = u_l;
}

/**
 *  Q-inflow Riemann solver for pulse wave propagation. This Riemann solver is
 * called by DoOdeProjection in case of a Q-inflow boundary condition. It is
 * based on the conservation of mass and total pressure and on the
 * characteristic information. For further details see "Pulse wave propagation
 * in the human vascular system", section 3.4.1 Returns the upwinded quantities
 * \f$(A_u,u_u)\f$ and stores them into the boundary values
 */
void QInflow::Q_RiemannSolver(NekDouble Q, NekDouble A_r, NekDouble u_r,
                              NekDouble A_0, NekDouble beta, NekDouble alpha,
                              NekDouble &Au, NekDouble &uu)
{
    NekDouble W2           = 0.0;
    NekDouble A_calc       = 0.0;
    NekDouble FA           = 0.0;
    NekDouble dFdA         = 0.0;
    NekDouble delta_A_calc = 0.0;
    NekDouble I            = 0.0;
    NekDouble c            = 0.0;

    int proceed  = 1;
    int iter     = 0;
    int MAX_ITER = 200;

    // Tolerances for the algorithm
    NekDouble Tol = 1.0e-10;

    // Riemann invariant \f$W_2(Ar,ur)\f$
    m_pressureArea->GetW2(W2, u_r, beta, A_r, A_0, alpha);

    // Newton Iteration (Area only)
    A_calc = A_r;
    while ((proceed) && (iter < MAX_ITER))
    {
        iter += 1;

        m_pressureArea->GetCharIntegral(I, beta, A_calc, A_0, alpha);
        m_pressureArea->GetC(c, beta, A_calc, A_0, alpha);

        FA           = Q - A_calc * (W2 + I);
        dFdA         = -W2 - I - c;
        delta_A_calc = FA / dFdA;
        A_calc      -= delta_A_calc;

        if (sqrt(delta_A_calc * delta_A_calc) < Tol)
        {
            proceed = 0;
        }
    }

    m_pressureArea->GetCharIntegral(I, beta, A_calc, A_0, alpha);

    // Obtain u_u and A_u
    uu = W2 + I;
    Au = A_calc;
}

} // namespace Nektar
