///////////////////////////////////////////////////////////////////////////////
//
// File: AUSM0Solver.h
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
// Description: AUSM0 Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_AUSM0SOLVER
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_AUSM0SOLVER

#include <CompressibleFlowSolver/RiemannSolvers/CompressibleSolver.h>

namespace Nektar
{
    class AUSM0Solver : public CompressibleSolver
    {
    public:
        static RiemannSolverSharedPtr create()
        {
            return RiemannSolverSharedPtr(
                new AUSM0Solver());
        }
        
        static std::string solverName;
        
    protected:
        AUSM0Solver();
        
        virtual void v_PointSolve(
            double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
            double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
            double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef);
        
        double M1Function(int A, double M);
        double M2Function(int A, double M);
        double M4Function(int A, double beta,  double M);
        double P5Function(int A, double alpha, double M);
    };
}

#endif
