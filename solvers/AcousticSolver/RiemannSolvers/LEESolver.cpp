///////////////////////////////////////////////////////////////////////////////
//
// File: LEESolver.cpp
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
// Description: Riemann solver base classs for the LEE equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <AcousticSolver/RiemannSolvers/LEESolver.h>

using namespace std;

namespace Nektar
{

/**
 *
 */
LEESolver::LEESolver(const LibUtilities::SessionReaderSharedPtr &pSession)
    : AcousticSolver(pSession)
{
    m_requiresRotation = true;
}

/**
 *
 */
void LEESolver::v_Solve(const int nDim,
                        const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                        const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
                        Array<OneD, Array<OneD, NekDouble>> &flux)
{
    int nTracePts = Fwd[0].size();

    Array<OneD, Array<OneD, NekDouble>> bfFwd(nDim + 3);
    Array<OneD, Array<OneD, NekDouble>> bfBwd(nDim + 3);
    for (int i = 0; i < nDim + 3; i++)
    {
        bfFwd[i] = Array<OneD, NekDouble>(nTracePts);
        bfBwd[i] = Array<OneD, NekDouble>(nTracePts);
    }

    GetRotBasefield(bfFwd, bfBwd);

    int expDim = nDim;
    NekDouble vF, wF;

    if (expDim == 1)
    {
        for (int i = 0; i < nTracePts; ++i)
        {
            v_PointSolve(
                  Fwd[0][i],   Fwd[1][i],   Fwd[2][i], 0.0,  0.0,
                  Bwd[0][i],   Bwd[1][i],   Bwd[2][i], 0.0,  0.0,
                bfFwd[0][i], bfFwd[1][i], bfFwd[2][i], 0.0,  0.0,
                bfBwd[0][i], bfBwd[1][i], bfBwd[2][i], 0.0,  0.0,
                 flux[0][i],  flux[1][i],  flux[2][i],  vF,   wF);
        }
    }
    else if (expDim == 2)
    {
        for (int i = 0; i < nTracePts; ++i)
        {
            v_PointSolve(
                  Fwd[0][i],   Fwd[1][i],   Fwd[2][i],   Fwd[3][i],  0.0,
                  Bwd[0][i],   Bwd[1][i],   Bwd[2][i],   Bwd[3][i],  0.0,
                bfFwd[0][i], bfFwd[1][i], bfFwd[2][i], bfFwd[3][i],  0.0,
                bfBwd[0][i], bfBwd[1][i], bfBwd[2][i], bfBwd[3][i],  0.0,
                 flux[0][i],  flux[1][i],  flux[2][i],  flux[3][i],   wF);
        }
    }
    else if (expDim == 3)
    {
        for (int i = 0; i < nTracePts; ++i)
        {
            v_PointSolve(
                  Fwd[0][i],   Fwd[1][i],   Fwd[2][i],   Fwd[3][i],   Fwd[4][i],
                  Bwd[0][i],   Bwd[1][i],   Bwd[2][i],   Bwd[3][i],   Bwd[4][i],
                bfFwd[0][i], bfFwd[1][i], bfFwd[2][i], bfFwd[3][i], bfFwd[4][i],
                bfBwd[0][i], bfBwd[1][i], bfBwd[2][i], bfBwd[3][i], bfBwd[4][i],
                 flux[0][i],  flux[1][i],  flux[2][i],  flux[3][i],  flux[4][i]);
        }
    }
}

} // namespace Nektar
