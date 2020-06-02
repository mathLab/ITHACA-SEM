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
    // m_pointSolve = false;
}

/// programmatic ctor
RoeSolver::RoeSolver(): CompressibleSolver()
{
    // m_pointSolve = false;
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
    const Array<OneD, const Array<OneD, NekDouble> > &fwd,
    const Array<OneD, const Array<OneD, NekDouble> > &bwd,
          Array<OneD,       Array<OneD, NekDouble> > &flux)
{
    static auto gamma = m_params["gamma"]();
    static auto nVars = fwd.num_elements();
    static auto spaceDim = nVars-2;

    using vec_t = AVX::VecData<NekDouble, AVX::SIMD_WIDTH_SIZE>;

    // get limit of vectorizable chunk
    size_t sizeScalar = fwd[0].num_elements();
    size_t sizeVec = (sizeScalar / AVX::SIMD_WIDTH_SIZE) * AVX::SIMD_WIDTH_SIZE;

    // SIMD loop
    size_t i = 0;
    for (; i < sizeVec; i+=AVX::SIMD_WIDTH_SIZE)
    {
        vec_t rhoL{}, rhouL{}, rhovL{}, rhowL{}, EL{};
        vec_t rhoR{}, rhouR{}, rhovR{}, rhowR{}, ER{};

        // load
        rhoL  = &(fwd[0][i]);
        rhouL = &(fwd[1][i]);
        EL    = &(fwd[spaceDim+1][i]);
        rhoR  = &(bwd[0][i]);
        rhouR = &(bwd[1][i]);
        ER    = &(bwd[spaceDim+1][i]);

        if (spaceDim == 2)
        {
            rhovL = &(fwd[2][i]);
            rhovR = &(bwd[2][i]);
        }
        else if (spaceDim == 3)
        {
            rhovL = &(fwd[2][i]);
            rhowL = &(fwd[3][i]);
            rhovR = &(bwd[2][i]);
            rhowR = &(bwd[3][i]);
        }

        vec_t rhof{}, rhouf{}, rhovf{}, rhowf{}, Ef{};

        RoeKernel(
            rhoL, rhouL, rhovL, rhowL, EL,
            rhoR, rhouR, rhovR, rhowR, ER,
            rhof, rhouf, rhovf, rhowf, Ef,
            gamma);

        // store
        rhof.store_nts(&(flux[0][i]));
        rhouf.store_nts(&(flux[1][i]));
        Ef.store_nts(&(flux[nVars-1][i]));
        if (spaceDim == 2)
        {
            rhovf.store_nts(&(flux[2][i]));
        }
        else if (spaceDim == 3)
        {
            rhovf.store_nts(&(flux[2][i]));
            rhowf.store_nts(&(flux[3][i]));
        }

    } // avx loop


    // spillover loop
    for (; i < sizeScalar; ++i)
    {
        NekDouble rhoL{}, rhouL{}, rhovL{}, rhowL{}, EL{};
        NekDouble rhoR{}, rhouR{}, rhovR{}, rhowR{}, ER{};

        // load
        rhoL  = fwd[0][i];
        rhouL = fwd[1][i];
        EL    = fwd[spaceDim+1][i];
        rhoR  = bwd[0][i];
        rhouR = bwd[1][i];
        ER    = bwd[spaceDim+1][i];

        if (spaceDim == 2)
        {
            rhovL = fwd[2][i];
            rhovR = bwd[2][i];
        }
        else if (spaceDim == 3)
        {
            rhovL = fwd[2][i];
            rhowL = fwd[3][i];
            rhovR = bwd[2][i];
            rhowR = bwd[3][i];
        }

        NekDouble rhof{}, rhouf{}, rhovf{}, rhowf{}, Ef{};

        RoeKernel(
            rhoL, rhouL, rhovL, rhowL, EL,
            rhoR, rhouR, rhovR, rhowR, ER,
            rhof, rhouf, rhovf, rhowf, Ef,
            gamma);

        // store
        flux[0][i] = rhof;
        flux[1][i] = rhouf;
        flux[nVars-1][i] = Ef;
        if (spaceDim == 2)
        {
            flux[2][i] = rhovf;
        }
        else if (spaceDim == 3)
        {
            flux[2][i] = rhovf;
            flux[3][i] = rhowf;
        }

    } // loop

}

} // namespace Nektar
