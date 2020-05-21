///////////////////////////////////////////////////////////////////////////////
//
// File: RoeSolver.cpp
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
// Description: Roe Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/RoeSolver.h>

#include <AVXOperators/AVXUtil.hpp>

namespace Nektar
{
std::string RoeSolver::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "Roe",
        RoeSolver::create,
        "Roe Riemann solver");

RoeSolver::RoeSolver(const LibUtilities::SessionReaderSharedPtr& pSession)
    : CompressibleSolver(pSession)
{
    m_pointSolve = false;
}

/// programmatic ctor
RoeSolver::RoeSolver(): CompressibleSolver()
{
    m_pointSolve = false;
}

/**
 * @brief Roe Riemann solver.
 *
 * Stated equations numbers are from:
 *
 *   "Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical
 *   Introduction", E. F. Toro (3rd edition, 2009).
 *
 * We follow the algorithm prescribed following equation 11.70.
 *
 * @param rhoL      Density left state.
 * @param rhoR      Density right state.
 * @param rhouL     x-momentum component left state.
 * @param rhouR     x-momentum component right state.
 * @param rhovL     y-momentum component left state.
 * @param rhovR     y-momentum component right state.
 * @param rhowL     z-momentum component left state.
 * @param rhowR     z-momentum component right state.
 * @param EL        Energy left state.
 * @param ER        Energy right state.
 * @param rhof      Computed Riemann flux for density.
 * @param rhouf     Computed Riemann flux for x-momentum component
 * @param rhovf     Computed Riemann flux for y-momentum component
 * @param rhowf     Computed Riemann flux for z-momentum component
 * @param Ef        Computed Riemann flux for energy.
 */
void RoeSolver::v_PointSolve(
    double rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
    double rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
    double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
{
    static NekDouble gamma = m_params["gamma"]();

    RoeKernel(
        rhoL, rhouL, rhovL, rhowL, EL,
        rhoR, rhouR, rhovR, rhowR, ER,
        rhof, rhouf, rhovf, rhowf, Ef,
        gamma);
}


void RoeSolver::v_ArraySolve(
    const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
    const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
          Array<OneD,       Array<OneD, NekDouble> > &flux)
{
    static auto gamma = m_params["gamma"]();
    static auto nVars = Fwd.num_elements();
    static auto spaceDim = nVars-2;

    using vec_t = AVX::VecData<NekDouble, AVX::SIMD_WIDTH_SIZE>;

    // copy to aligned vectors
    size_t sizeScalar = Fwd[0].num_elements();
    size_t pad = sizeScalar % AVX::SIMD_WIDTH_SIZE;
    size_t sizeVec = sizeScalar / AVX::SIMD_WIDTH_SIZE + (pad == 0 ? 0 : 1);

    std::vector<AVX::AlignedVector<vec_t>> alignedFwd(nVars),
        alignedBwd(nVars), alignedFlux(nVars);

    for (size_t i = 0; i < nVars; ++i)
    {
        alignedFwd[i] = AVX::AlignedVector<vec_t>(sizeVec);
        alignedBwd[i] = AVX::AlignedVector<vec_t>(sizeVec);
        alignedFlux[i] = AVX::AlignedVector<vec_t>(sizeVec);

        AVX::CopyToAlignedVector(Fwd[i], alignedFwd[i]);
        AVX::CopyToAlignedVector(Bwd[i], alignedBwd[i]);
    }

    // AVX loop
    for (size_t i = 0; i < sizeVec; ++i)
    {
        vec_t rhoL{}, rhouL{}, rhovL{}, rhowL{}, EL{};
        vec_t rhoR{}, rhouR{}, rhovR{}, rhowR{}, ER{};

        // load
        rhoL  = alignedFwd[0][i];
        rhouL = alignedFwd[1][i];
        EL    = alignedFwd[spaceDim+1][i];
        rhoR  = alignedBwd[0][i];
        rhouR = alignedBwd[1][i];
        ER    = alignedBwd[spaceDim+1][i];

        if (spaceDim == 2)
        {
            rhovL = alignedFwd[2][i];
            rhovR = alignedBwd[2][i];
        }
        else if (spaceDim == 3)
        {
            rhovL = alignedFwd[2][i];
            rhowL = alignedFwd[3][i];
            rhovR = alignedBwd[2][i];
            rhowR = alignedBwd[3][i];
        }

        vec_t rhof{}, rhouf{}, rhovf{}, rhowf{}, Ef{};

        RoeKernel(
            rhoL, rhouL, rhovL, rhowL, EL,
            rhoR, rhouR, rhovR, rhowR, ER,
            rhof, rhouf, rhovf, rhowf, Ef,
            gamma);

        // store
        alignedFlux[0][i] = rhof;
        alignedFlux[1][i] = rhouf;
        alignedFlux[nVars-1][i] = Ef;
        if (spaceDim == 2)
        {
            alignedFlux[2][i] = rhovf;
        }
        else if (spaceDim == 3)
        {
            alignedFlux[2][i] = rhovf;
            alignedFlux[3][i] = rhowf;
        }

    } // avx loop

    for (size_t i = 0; i < nVars; ++i)
    {
        AVX::CopyFromAlignedVector(alignedFlux[i], flux[i]);
    }
}

}
