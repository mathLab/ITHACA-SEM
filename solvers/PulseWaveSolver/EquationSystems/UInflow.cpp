///////////////////////////////////////////////////////////////////////////////
//
// File UInflow.cpp
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
// Description: UInflow class
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/UInflow.h>

namespace Nektar
{

string UInflow::className = GetBoundaryFactory().RegisterCreatorFunction(
    "U-inflow", UInflow::create, "Velocity inflow boundary condition");

/**
 *
 */
// ExpListSharedPtr = 'Expansion list shared pointer'
UInflow::UInflow(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
                 const LibUtilities::SessionReaderSharedPtr pSession,
                 PulseWavePressureAreaSharedPtr pressureArea)
    : PulseWaveBoundary(pVessel, pSession, pressureArea)
{
    // Constructor
}

/**
 *
 */
UInflow::~UInflow()
{
    // Destructor
}

void UInflow::v_DoBoundary(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &A_0,
    Array<OneD, Array<OneD, NekDouble>> &beta, const NekDouble time, int omega,
    int offset, int n)
{
    NekDouble A;
    NekDouble u;
    NekDouble A_r;
    NekDouble u_r;
    NekDouble A_l;
    NekDouble u_l;
    Array<OneD, MultiRegions::ExpListSharedPtr> vessel(2);

    // Pointers to the domains
    vessel[0] = m_vessels[2 * omega];
    vessel[1] = m_vessels[2 * omega + 1];

    // Evaluates the time-dependent U
    vessel[1]->EvaluateBoundaryConditions(time);

    // Read the BC values from the input file
    A  = (vessel[0]->UpdateBndCondExpansion(n))->GetCoeffs()[0];
    u   = (vessel[1]->UpdateBndCondExpansion(n))->GetCoeffs()[0];

    // Initial conditions in inlet vessel
    A_r = inarray[0][offset];
    u_r = inarray[1][offset];

    /* Fix the boundary conditions in the virtual vessel to ensure
       upwind state matches the boundary condition at the next time step */
    A_l = A_r;
    u_l = 2 * u - u_r;

    // Store the updated values in the boundary condition
    (vessel[0]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = A_l;
    (vessel[1]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = u_l;
}

} // namespace Nektar
