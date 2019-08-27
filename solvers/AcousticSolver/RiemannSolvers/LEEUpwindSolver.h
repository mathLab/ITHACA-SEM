///////////////////////////////////////////////////////////////////////////////
//
// File: LEEUpwindSolver.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Lax-Friedrichs solver for the LEE equations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ACOUSTICSOLVER_RIEMANNSOLVERS_LEELAXFRIEDRICHSSOLVER
#define NEKTAR_SOLVERS_ACOUSTICSOLVER_RIEMANNSOLVERS_LEELAXFRIEDRICHSSOLVER

#include <AcousticSolver/RiemannSolvers/LEESolver.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

class LEEUpwindSolver : public LEESolver
{
public:
    static RiemannSolverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession)
    {
        return RiemannSolverSharedPtr(new LEEUpwindSolver(pSession));
    }

    static std::string solverName;

protected:
    LEEUpwindSolver(const LibUtilities::SessionReaderSharedPtr &pSession);

    virtual void v_PointSolve(
        NekDouble  pL,    NekDouble  rhoL,  NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL,
        NekDouble  pR,    NekDouble  rhoR,  NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR,
        NekDouble  c0sqL, NekDouble  rho0L, NekDouble  u0L,   NekDouble  v0L,   NekDouble  w0L,
        NekDouble  c0sqR, NekDouble  rho0R, NekDouble  u0R,   NekDouble  v0R,   NekDouble  w0R,
        NekDouble &pF,    NekDouble &rhoF,  NekDouble &rhouF, NekDouble &rhovF, NekDouble &rhowF);
};

} // namespace Nektar
#endif
