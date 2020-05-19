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

#include <AVXOperators/VecData.hpp>

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
    double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
    double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
    double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
{
    static NekDouble gamma = m_params["gamma"]();

    // Left and right velocities
    NekDouble uL = rhouL / rhoL;
    NekDouble vL = rhovL / rhoL;
    NekDouble wL = rhowL / rhoL;
    NekDouble uR = rhouR / rhoR;
    NekDouble vR = rhovR / rhoR;
    NekDouble wR = rhowR / rhoR;

    // Left and right pressures
    NekDouble pL = (gamma - 1.0) *
        (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
    NekDouble pR = (gamma - 1.0) *
        (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));

    // Left and right enthalpy
    NekDouble hL = (EL + pL) / rhoL;
    NekDouble hR = (ER + pR) / rhoR;

    // Square root of rhoL and rhoR.
    NekDouble srL  = sqrt(rhoL);
    NekDouble srR  = sqrt(rhoR);
    NekDouble srLR = srL + srR;

    // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
    NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
    NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
    NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
    NekDouble hRoe   = (srL * hL + srR * hR) / srLR;
    NekDouble URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
    NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));

    // Compute eigenvectors (equation 11.59).
    NekDouble k[5][5] = {
        {1, uRoe - cRoe, vRoe, wRoe, hRoe - uRoe * cRoe},
        {1, uRoe,        vRoe, wRoe, 0.5 * URoe},
        {0, 0,           1,    0,    vRoe},
        {0, 0,           0,    1,    wRoe},
        {1, uRoe+cRoe,  vRoe,  wRoe, hRoe + uRoe*cRoe}
    };

    // Calculate jumps \Delta u_i (defined preceding equation 11.67).
    NekDouble jump[5] = {
        rhoR  - rhoL,
        rhouR - rhouL,
        rhovR - rhovL,
        rhowR - rhowL,
        ER    - EL
    };

    // Define \Delta u_5 (equation 11.70).
    NekDouble jumpbar = jump[4] - (jump[2]-vRoe*jump[0])*vRoe -
        (jump[3]-wRoe*jump[0])*wRoe;

    // Compute wave amplitudes (equations 11.68, 11.69).
    NekDouble alpha[5];
    alpha[1] = (gamma-1.0)*(jump[0]*(hRoe - uRoe*uRoe) + uRoe*jump[1] -
                            jumpbar)/(cRoe*cRoe);
    alpha[0] = (jump[0]*(uRoe + cRoe) - jump[1] - cRoe*alpha[1])/(2.0*cRoe);
    alpha[4] = jump[0] - (alpha[0] + alpha[1]);
    alpha[2] = jump[2] - vRoe * jump[0];
    alpha[3] = jump[3] - wRoe * jump[0];

    // Compute average of left and right fluxes needed for equation 11.29.
    rhof  = 0.5*(rhoL*uL + rhoR*uR);
    rhouf = 0.5*(pL + rhoL*uL*uL + pR + rhoR*uR*uR);
    rhovf = 0.5*(rhoL*uL*vL + rhoR*uR*vR);
    rhowf = 0.5*(rhoL*uL*wL + rhoR*uR*wR);
    Ef    = 0.5*(uL*(EL + pL) + uR*(ER + pR));

    // Compute eigenvalues \lambda_i (equation 11.58).
    NekDouble uRoeAbs = std::abs(uRoe);
    NekDouble lambda[5] = {
        std::abs(uRoe - cRoe),
        uRoeAbs,
        uRoeAbs,
        uRoeAbs,
        std::abs(uRoe + cRoe)
    };

    // Finally perform summation (11.29).
    for (int i = 0; i < 5; ++i)
    {
        uRoeAbs = 0.5*alpha[i]*lambda[i];

        rhof  -= uRoeAbs*k[i][0];
        rhouf -= uRoeAbs*k[i][1];
        rhovf -= uRoeAbs*k[i][2];
        rhowf -= uRoeAbs*k[i][3];
        Ef    -= uRoeAbs*k[i][4];
    }
}


void RoeSolver::v_ArraySolve(
    const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
    const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
          Array<OneD,       Array<OneD, NekDouble> > &flux)
{
    static NekDouble gamma = m_params["gamma"]();
    static auto spaceDim = Fwd.num_elements()-2;

    using vec_t = AVX::VecData<NekDouble, AVX::SIMD_WIDTH_SIZE>;

    // Scalar Step to alignement
    size_t cnt = 0;
    size_t nPts = Fwd[0].num_elements();
    for (; cnt < nPts && reinterpret_cast<size_t>(&(Fwd[0][cnt])) % AVX::SIMD_WIDTH_BYTES != 0; ++cnt)
    {
        std::cout << "single step\n"
            << reinterpret_cast<size_t>(&(Fwd[0][cnt])) % AVX::SIMD_WIDTH_BYTES << '\n'
            << reinterpret_cast<size_t>(&(Fwd[1][cnt])) % AVX::SIMD_WIDTH_BYTES << '\n'
            << reinterpret_cast<size_t>(&(Fwd[2][cnt])) % AVX::SIMD_WIDTH_BYTES << std::endl;

        NekDouble  rhoL{};
        NekDouble  rhouL{};
        NekDouble  rhovL{};
        NekDouble  rhowL{};
        NekDouble  EL{};
        NekDouble  rhoR{};
        NekDouble  rhouR{};
        NekDouble  rhovR{};
        NekDouble  rhowR{};
        NekDouble  ER{};

        rhoL  = Fwd[0][cnt];
        rhouL = Fwd[1][cnt];
        EL    = Fwd[spaceDim+1][cnt];
        rhoR  = Bwd[0][cnt];
        rhouR = Bwd[1][cnt];
        ER    = Bwd[spaceDim+1][cnt];

        if (spaceDim == 2)
        {
            rhovL = Fwd[2][cnt];
            rhovR = Bwd[2][cnt];
        }
        else if (spaceDim == 3)
        {
            rhovL = Fwd[2][cnt];
            rhowL = Fwd[3][cnt];
            rhovR = Bwd[2][cnt];
            rhowR = Bwd[3][cnt];
        }


        // Left and right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;

        // Left and right pressures
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));

        // Left and right enthalpy
        NekDouble hL = (EL + pL) / rhoL;
        NekDouble hR = (ER + pR) / rhoR;

        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;

        // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble hRoe   = (srL * hL + srR * hR) / srLR;
        NekDouble URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));

        // Compute eigenvectors (equation 11.59).
        NekDouble k[5][5] = {
            {1, uRoe - cRoe, vRoe, wRoe, hRoe - uRoe * cRoe},
            {1, uRoe,        vRoe, wRoe, 0.5 * URoe},
            {0, 0,           1,    0,    vRoe},
            {0, 0,           0,    1,    wRoe},
            {1, uRoe+cRoe,  vRoe,  wRoe, hRoe + uRoe*cRoe}
        };

        // Calculate jumps \Delta u_i (defined preceding equation 11.67).
        NekDouble jump[5] = {
            rhoR  - rhoL,
            rhouR - rhouL,
            rhovR - rhovL,
            rhowR - rhowL,
            ER    - EL
        };

        // Define \Delta u_5 (equation 11.70).
        NekDouble jumpbar = jump[4] - (jump[2]-vRoe*jump[0])*vRoe -
            (jump[3]-wRoe*jump[0])*wRoe;

        // Compute wave amplitudes (equations 11.68, 11.69).
        NekDouble alpha[5];
        alpha[1] = (gamma-1.0)*(jump[0]*(hRoe - uRoe*uRoe) + uRoe*jump[1] -
                                jumpbar)/(cRoe*cRoe);
        alpha[0] = (jump[0]*(uRoe + cRoe) - jump[1] - cRoe*alpha[1])/(2.0*cRoe);
        alpha[4] = jump[0] - (alpha[0] + alpha[1]);
        alpha[2] = jump[2] - vRoe * jump[0];
        alpha[3] = jump[3] - wRoe * jump[0];

        // Compute average of left and right fluxes needed for equation 11.29.
        NekDouble rhof  = 0.5*(rhoL*uL + rhoR*uR);
        NekDouble rhouf = 0.5*(pL + rhoL*uL*uL + pR + rhoR*uR*uR);
        NekDouble rhovf = 0.5*(rhoL*uL*vL + rhoR*uR*vR);
        NekDouble rhowf = 0.5*(rhoL*uL*wL + rhoR*uR*wR);
        NekDouble Ef    = 0.5*(uL*(EL + pL) + uR*(ER + pR));

        // Compute eigenvalues \lambda_i (equation 11.58).
        NekDouble uRoeAbs = std::abs(uRoe);
        NekDouble lambda[5] = {
            std::abs(uRoe - cRoe),
            uRoeAbs,
            uRoeAbs,
            uRoeAbs,
            std::abs(uRoe + cRoe)
        };

        // Finally perform summation (11.29).
        for (size_t i = 0; i < 5; ++i)
        {
            uRoeAbs = 0.5*alpha[cnt]*lambda[cnt];

            rhof  -= uRoeAbs*k[cnt][0];
            rhouf -= uRoeAbs*k[cnt][1];
            rhovf -= uRoeAbs*k[cnt][2];
            rhowf -= uRoeAbs*k[cnt][3];
            Ef    -= uRoeAbs*k[cnt][4];
        }

        // store
        flux[0][cnt] = rhof;
        flux[1][cnt] = rhouf;
        flux[spaceDim+1][cnt] = Ef;
        if (spaceDim == 2)
        {
            flux[2][cnt] = rhovf;
        }
        else if (spaceDim == 3)
        {
            flux[2][cnt] = rhovf;
            flux[3][cnt] = rhowf;
        }

    }

    // AVX loop
    for (; (nPts - cnt) >= AVX::SIMD_WIDTH_SIZE; cnt += AVX::SIMD_WIDTH_SIZE)
    {
        std::cout << "vector step\n"
            << reinterpret_cast<size_t>(&(Fwd[0][cnt])) % AVX::SIMD_WIDTH_BYTES << '\n'
            << reinterpret_cast<size_t>(&(Fwd[1][cnt])) % AVX::SIMD_WIDTH_BYTES << '\n'
            << reinterpret_cast<size_t>(&(Fwd[2][cnt])) % AVX::SIMD_WIDTH_BYTES << std::endl;
        vec_t rhoL{};
        vec_t rhouL{};
        vec_t rhovL{};
        vec_t rhowL{};
        vec_t EL{};
        vec_t rhoR{};
        vec_t rhouR{};
        vec_t rhovR{};
        vec_t rhowR{};
        vec_t ER{};

        // load
        rhoL  = Fwd[0][cnt];
        rhouL = Fwd[1][cnt];
        EL    = Fwd[spaceDim+1][cnt];
        rhoR  = Bwd[0][cnt];
        rhouR = Bwd[1][cnt];
        ER    = Bwd[spaceDim+1][cnt];

        if (spaceDim == 2)
        {
            rhovL = Fwd[2][cnt];
            rhovR = Bwd[2][cnt];
        }
        else if (spaceDim == 3)
        {
            rhovL = Fwd[2][cnt];
            rhowL = Fwd[3][cnt];
            rhovR = Bwd[2][cnt];
            rhowR = Bwd[3][cnt];
        }


        // Left and right velocities
        vec_t uL = rhouL / rhoL;
        vec_t vL = rhovL / rhoL;
        vec_t wL = rhowL / rhoL;
        vec_t uR = rhouR / rhoR;
        vec_t vR = rhovR / rhoR;
        vec_t wR = rhowR / rhoR;

        // Left and right pressures
        vec_t pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        vec_t pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));

        // Left and right enthalpy
        vec_t hL = (EL + pL) / rhoL;
        vec_t hR = (ER + pR) / rhoR;

        // Square root of rhoL and rhoR.
        vec_t srL  = sqrt(rhoL);
        vec_t srR  = sqrt(rhoR);
        vec_t srLR = srL + srR;

        // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
        vec_t uRoe   = (srL * uL + srR * uR) / srLR;
        vec_t vRoe   = (srL * vL + srR * vR) / srLR;
        vec_t wRoe   = (srL * wL + srR * wR) / srLR;
        vec_t hRoe   = (srL * hL + srR * hR) / srLR;
        vec_t URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
        vec_t cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));

        // Compute eigenvectors (equation 11.59).
        vec_t k[5][5] = {
            {1., uRoe - cRoe, vRoe, wRoe, hRoe - uRoe * cRoe},
            {1., uRoe,        vRoe, wRoe, 0.5 * URoe},
            {0., 0.,           1.,    0.,    vRoe},
            {0., 0.,           0.,    1.,    wRoe},
            {1., uRoe+cRoe,  vRoe,  wRoe, hRoe + uRoe*cRoe}
        };

        // Calculate jumps \Delta u_i (defined preceding equation 11.67).
        vec_t jump[5] = {
            rhoR  - rhoL,
            rhouR - rhouL,
            rhovR - rhovL,
            rhowR - rhowL,
            ER    - EL
        };

        // Define \Delta u_5 (equation 11.70).
        vec_t jumpbar = jump[4] - (jump[2]-vRoe*jump[0])*vRoe -
            (jump[3]-wRoe*jump[0])*wRoe;

        // Compute wave amplitudes (equations 11.68, 11.69).
        vec_t alpha[5];
        alpha[1] = (gamma-1.0)*(jump[0]*(hRoe - uRoe*uRoe) + uRoe*jump[1] -
                                jumpbar)/(cRoe*cRoe);
        alpha[0] = (jump[0]*(uRoe + cRoe) - jump[1] - cRoe*alpha[1])/(2.0*cRoe);
        alpha[4] = jump[0] - (alpha[0] + alpha[1]);
        alpha[2] = jump[2] - vRoe * jump[0];
        alpha[3] = jump[3] - wRoe * jump[0];

        // Compute average of left and right fluxes needed for equation 11.29.
        vec_t rhof  = 0.5*(rhoL*uL + rhoR*uR);
        vec_t rhouf = 0.5*(pL + rhoL*uL*uL + pR + rhoR*uR*uR);
        vec_t rhovf = 0.5*(rhoL*uL*vL + rhoR*uR*vR);
        vec_t rhowf = 0.5*(rhoL*uL*wL + rhoR*uR*wR);
        vec_t Ef    = 0.5*(uL*(EL + pL) + uR*(ER + pR));

        // Compute eigenvalues \lambda_i (equation 11.58).
        vec_t uRoeAbs = abs(uRoe);
        vec_t lambda[5] = {
            abs(uRoe - cRoe),
            uRoeAbs,
            uRoeAbs,
            uRoeAbs,
            abs(uRoe + cRoe)
        };

        // Finally perform summation (11.29).
        for (size_t i = 0; i < 5; ++i)
        {
            uRoeAbs = 0.5*alpha[cnt]*lambda[cnt];

            rhof  = rhof  - uRoeAbs*k[cnt][0];
            rhouf = rhouf - uRoeAbs*k[cnt][1];
            rhovf = rhovf - uRoeAbs*k[cnt][2];
            rhowf = rhowf - uRoeAbs*k[cnt][3];
            Ef    = Ef    - uRoeAbs*k[cnt][4];
        }


        // store
        std::cout << "before store" << std::endl;
        std::cout << &(Fwd[0][cnt]) << '\n';
        std::cout << reinterpret_cast<size_t>(&(Fwd[0][cnt])) % AVX::SIMD_WIDTH_BYTES << '\n';
        NekDouble* tmp = &(flux[0][cnt]);
        rhof.store(tmp);
        std::cout << "after rhof store" << std::endl;
        rhouf.store(&(flux[1][cnt]));
        Ef.store(&(flux[spaceDim+1][cnt]));
        if (spaceDim == 2)
        {
            rhovf.store(&(flux[2][cnt]));
        }
        else if (spaceDim == 3)
        {
            rhovf.store(&(flux[2][cnt]));
            rhowf.store(&(flux[3][cnt]));
        }

    }


    // spill over loop
    for (; cnt < nPts; ++cnt)
    {
        std::cout << "vector step\n"
            << reinterpret_cast<size_t>(&(Fwd[0][cnt])) % AVX::SIMD_WIDTH_BYTES << '\n'
            << reinterpret_cast<size_t>(&(Fwd[1][cnt])) % AVX::SIMD_WIDTH_BYTES << '\n'
            << reinterpret_cast<size_t>(&(Fwd[2][cnt])) % AVX::SIMD_WIDTH_BYTES << std::endl;
        NekDouble  rhoL{};
        NekDouble  rhouL{};
        NekDouble  rhovL{};
        NekDouble  rhowL{};
        NekDouble  EL{};
        NekDouble  rhoR{};
        NekDouble  rhouR{};
        NekDouble  rhovR{};
        NekDouble  rhowR{};
        NekDouble  ER{};

        rhoL  = Fwd[0][cnt];
        rhouL = Fwd[1][cnt];
        EL    = Fwd[spaceDim+1][cnt];
        rhoR  = Bwd[0][cnt];
        rhouR = Bwd[1][cnt];
        ER    = Bwd[spaceDim+1][cnt];

        if (spaceDim == 2)
        {
            rhovL = Fwd[2][cnt];
            rhovR = Bwd[2][cnt];
        }
        else if (spaceDim == 3)
        {
            rhovL = Fwd[2][cnt];
            rhowL = Fwd[3][cnt];
            rhovR = Bwd[2][cnt];
            rhowR = Bwd[3][cnt];
        }


        // Left and right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;

        // Left and right pressures
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));

        // Left and right enthalpy
        NekDouble hL = (EL + pL) / rhoL;
        NekDouble hR = (ER + pR) / rhoR;

        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;

        // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble hRoe   = (srL * hL + srR * hR) / srLR;
        NekDouble URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));

        // Compute eigenvectors (equation 11.59).
        NekDouble k[5][5] = {
            {1, uRoe - cRoe, vRoe, wRoe, hRoe - uRoe * cRoe},
            {1, uRoe,        vRoe, wRoe, 0.5 * URoe},
            {0, 0,           1,    0,    vRoe},
            {0, 0,           0,    1,    wRoe},
            {1, uRoe+cRoe,  vRoe,  wRoe, hRoe + uRoe*cRoe}
        };

        // Calculate jumps \Delta u_i (defined preceding equation 11.67).
        NekDouble jump[5] = {
            rhoR  - rhoL,
            rhouR - rhouL,
            rhovR - rhovL,
            rhowR - rhowL,
            ER    - EL
        };

        // Define \Delta u_5 (equation 11.70).
        NekDouble jumpbar = jump[4] - (jump[2]-vRoe*jump[0])*vRoe -
            (jump[3]-wRoe*jump[0])*wRoe;

        // Compute wave amplitudes (equations 11.68, 11.69).
        NekDouble alpha[5];
        alpha[1] = (gamma-1.0)*(jump[0]*(hRoe - uRoe*uRoe) + uRoe*jump[1] -
                                jumpbar)/(cRoe*cRoe);
        alpha[0] = (jump[0]*(uRoe + cRoe) - jump[1] - cRoe*alpha[1])/(2.0*cRoe);
        alpha[4] = jump[0] - (alpha[0] + alpha[1]);
        alpha[2] = jump[2] - vRoe * jump[0];
        alpha[3] = jump[3] - wRoe * jump[0];

        // Compute average of left and right fluxes needed for equation 11.29.
        NekDouble rhof  = 0.5*(rhoL*uL + rhoR*uR);
        NekDouble rhouf = 0.5*(pL + rhoL*uL*uL + pR + rhoR*uR*uR);
        NekDouble rhovf = 0.5*(rhoL*uL*vL + rhoR*uR*vR);
        NekDouble rhowf = 0.5*(rhoL*uL*wL + rhoR*uR*wR);
        NekDouble Ef    = 0.5*(uL*(EL + pL) + uR*(ER + pR));

        // Compute eigenvalues \lambda_i (equation 11.58).
        NekDouble uRoeAbs = std::abs(uRoe);
        NekDouble lambda[5] = {
            std::abs(uRoe - cRoe),
            uRoeAbs,
            uRoeAbs,
            uRoeAbs,
            std::abs(uRoe + cRoe)
        };

        // Finally perform summation (11.29).
        for (size_t i = 0; i < 5; ++i)
        {
            uRoeAbs = 0.5*alpha[cnt]*lambda[cnt];

            rhof  -= uRoeAbs*k[cnt][0];
            rhouf -= uRoeAbs*k[cnt][1];
            rhovf -= uRoeAbs*k[cnt][2];
            rhowf -= uRoeAbs*k[cnt][3];
            Ef    -= uRoeAbs*k[cnt][4];
        }

        // store
        flux[0][cnt] = rhof;
        flux[1][cnt] = rhouf;
        flux[spaceDim+1][cnt] = Ef;
        if (spaceDim == 2)
        {
            flux[2][cnt] = rhovf;
        }
        else if (spaceDim == 3)
        {
            flux[2][cnt] = rhovf;
            flux[3][cnt] = rhowf;
        }

    }
}
}
