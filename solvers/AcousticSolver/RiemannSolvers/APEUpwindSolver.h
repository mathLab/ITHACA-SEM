///////////////////////////////////////////////////////////////////////////////
//
// File: APEUpwindSolver.h
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
// Description: Upwind Riemann solver for the APE equations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ACOUSTICSOLVER_RIEMANNSOLVERS_UPWINDSOLVER
#define NEKTAR_SOLVERS_ACOUSTICSOLVER_RIEMANNSOLVERS_UPWINDSOLVER

#include <AcousticSolver/RiemannSolvers/AcousticSolver.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

class APEUpwindSolver : public AcousticSolver
{
public:
    static RiemannSolverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession)
    {
        return RiemannSolverSharedPtr(new APEUpwindSolver(pSession));
    }

    static std::string solverName;

protected:
    APEUpwindSolver(const LibUtilities::SessionReaderSharedPtr &pSession);

    virtual void v_PointSolve(
        NekDouble  pL,    NekDouble  rhoL,  NekDouble  uL,  NekDouble  vL,  NekDouble  wL,
        NekDouble  pR,    NekDouble  rhoR,  NekDouble  uR,  NekDouble  vR,  NekDouble  wR,
        NekDouble  c0sqL, NekDouble  rho0L, NekDouble  u0L, NekDouble  v0L, NekDouble  w0L,
        NekDouble  c0sqR, NekDouble  rho0R, NekDouble  u0R, NekDouble  v0R, NekDouble  w0R,
        NekDouble &pF,    NekDouble &rhoF,  NekDouble &uF,  NekDouble &vF,  NekDouble &wF);
};

} // namespace Nektar

#endif
