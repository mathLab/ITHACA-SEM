///////////////////////////////////////////////////////////////////////////////
//
// File: UpwindPulseSolver.h
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
// Description: Upwind pulse Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PULSEWAVESOLVER_RIEMANNSOLVER_UPWINDPULSE
#define NEKTAR_SOLVERS_PULSEWAVESOLVER_RIEMANNSOLVER_UPWINDPULSE

#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>
#include <PulseWaveSolver/EquationSystems/PulseWavePressureArea.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{
class UpwindPulseSolver : public RiemannSolver
{
public:
    static RiemannSolverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    {
        return RiemannSolverSharedPtr(new UpwindPulseSolver(pSession));
    }

    static std::string solverName;

protected:
    UpwindPulseSolver(const LibUtilities::SessionReaderSharedPtr& pSession);

    LibUtilities::SessionReaderSharedPtr            m_session;
    int                                             m_nVariables;
    Array<OneD, MultiRegions::ExpListSharedPtr>     m_vessels;
    PulseWavePressureAreaSharedPtr                  m_pressureArea;

    virtual void v_Solve(const int nDim,
                         const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                         const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
                         Array<OneD, Array<OneD, NekDouble>> &flux);

    void RiemannSolverUpwind(NekDouble AL, NekDouble uL, NekDouble AR,
                             NekDouble uR, NekDouble &Aflux, NekDouble &uflux,
                             NekDouble A0, NekDouble beta, NekDouble n, NekDouble alpha = 0.5);
};
}

#endif
