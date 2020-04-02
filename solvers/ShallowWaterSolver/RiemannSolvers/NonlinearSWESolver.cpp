///////////////////////////////////////////////////////////////////////////////
//
// File: NonlinearSWESolver.cpp
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
// Description: Shallow water Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <ShallowWaterSolver/RiemannSolvers/NonlinearSWESolver.h>

namespace Nektar
{
    NonlinearSWESolver::NonlinearSWESolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : RiemannSolver(pSession), m_pointSolve(true)
    {
        m_requiresRotation = true;
    }

    void  NonlinearSWESolver::v_Solve(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        boost::ignore_unused(nDim);

        if (m_pointSolve)
        {
            int expDim = Fwd.size()-1;
            NekDouble hvf;

            if (expDim == 1)
            {
                for (int i = 0; i < Fwd[0].size(); ++i)
                {
                    v_PointSolve(
                        Fwd [0][i], Fwd [1][i], 0.0,
                        Bwd [0][i], Bwd [1][i], 0.0,
                        flux[0][i], flux[1][i], hvf);
                }
            }
            else if (expDim == 2)
            {
                for (int i = 0; i < Fwd[0].size(); ++i)
                {
                    v_PointSolve(
                        Fwd [0][i], Fwd [1][i], Fwd [2][i],
                        Bwd [0][i], Bwd [1][i], Bwd [2][i],
                        flux[0][i], flux[1][i], flux[2][i]);
                }
            }
            else if (expDim == 3)
            {
	      ASSERTL0(false, "No 3D Shallow water supported.");
	    }
        }
        else
        {
            v_ArraySolve(Fwd, Bwd, flux);
        }
    }
}
