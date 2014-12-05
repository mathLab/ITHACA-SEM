///////////////////////////////////////////////////////////////////////////////
//
// File: LaxFriedrichsSolver.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2014 Kilian Lackhove
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
// Description: Lax-Friedrichs solver for the APE equations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_APESOLVER_RIEMANNSOLVERS_LAXFRIEDRICHSSOLVER
#define NEKTAR_SOLVERS_APESOLVER_RIEMANNSOLVERS_LAXFRIEDRICHSSOLVER

#include <SolverUtils/SolverUtilsDeclspec.h>
#include <APESolver/RiemannSolvers/APESolver.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

class LaxFriedrichsSolver : public APESolver
{
public:
    static RiemannSolverSharedPtr create()
    {
        return RiemannSolverSharedPtr(new LaxFriedrichsSolver());
    }

    static std::string solverName;

protected:
    LaxFriedrichsSolver();

    virtual void v_PointSolve(
        NekDouble  pL, NekDouble  uL, NekDouble  vL, NekDouble  wL,
        NekDouble  pR, NekDouble  uR, NekDouble  vR, NekDouble  wR,
        NekDouble  p0, NekDouble  u0, NekDouble  v0, NekDouble  w0,
        NekDouble &pF, NekDouble &uF, NekDouble &vF, NekDouble &wF);
};

}
#endif
