///////////////////////////////////////////////////////////////////////////////
//
// File: CompressibleSolver.cpp
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
// Description: Compressible Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/CompressibleSolver.h>

namespace Nektar
{
    CompressibleSolver::CompressibleSolver() : RiemannSolver(),
                                               m_pointSolve(true)
    {
        m_requiresRotation = true;
    }
    
    void CompressibleSolver::v_Solve(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        if (m_pointSolve)
        {
            int expDim      = nDim;
            int nvariables  = Fwd.num_elements();
            
            NekDouble rhouf, rhovf;
            
            // Check if PDE-based SC is used
            if (expDim == 1)
            {
                for (int i = 0; i < Fwd[0].num_elements(); ++i)
                {
                    v_PointSolve(
                        Fwd [0][i], Fwd [1][i], 0.0,   0.0,   Fwd [2][i],
                        Bwd [0][i], Bwd [1][i], 0.0,   0.0,   Bwd [2][i],
                        flux[0][i], flux[1][i], rhouf, rhovf, flux[2][i]);
                }
            }
            else if (expDim == 2)
            {
                if (nvariables == expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolve(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], 0.0,   Fwd [3][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], 0.0,   Bwd [3][i],
                            flux[0][i], flux[1][i], flux[2][i], rhovf, flux[3][i]);
                    }
                }
                
                if (nvariables > expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolveVisc(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], 0.0, Fwd [3][i], Fwd [4][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], 0.0, Bwd [3][i], Bwd [4][i],
                            flux[0][i], flux[1][i], flux[2][i], rhovf, flux[3][i], flux[4][i]);
                    }
                }
                
            }
            else if (expDim == 3)
            {
                for (int i = 0; i < Fwd[0].num_elements(); ++i)
                {
                    v_PointSolve(
                        Fwd [0][i], Fwd [1][i], Fwd [2][i], Fwd [3][i], Fwd [4][i],
                        Bwd [0][i], Bwd [1][i], Bwd [2][i], Bwd [3][i], Bwd [4][i],
                        flux[0][i], flux[1][i], flux[2][i], flux[3][i], flux[4][i]);
                }
                if (nvariables > expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolveVisc(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], Fwd [3][i], Fwd [4][i], Fwd [5][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], Bwd [3][i], Bwd [4][i], Bwd [5][i],
                            flux[0][i], flux[1][i], flux[2][i], flux[3][i], flux[4][i], flux[5][i]);
                    }
                }
            }
        }
        else
        {
            v_ArraySolve(Fwd, Bwd, flux);
        }
    }
}
