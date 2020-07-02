///////////////////////////////////////////////////////////////////////////////
//
// File: RoeSolver.h
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_ROESOLVER
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_ROESOLVER

#include <CompressibleFlowSolver/RiemannSolvers/CompressibleSolver.h>

namespace Nektar
{
class RoeSolver : public CompressibleSolver
{
public:
    static RiemannSolverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    {
        return RiemannSolverSharedPtr(
            new RoeSolver(pSession));
    }

    static std::string solverName;

    /// programmatic ctor
    RoeSolver();

protected:
    RoeSolver(const LibUtilities::SessionReaderSharedPtr& pSession);

    using ND = NekDouble;

    void v_PointSolve(
        ND  rhoL, ND  rhouL, ND  rhovL, ND  rhowL, ND  EL,
        ND  rhoR, ND  rhouR, ND  rhovR, ND  rhowR, ND  ER,
        ND &rhof, ND &rhouf, ND &rhovf, ND &rhowf, ND &Ef) final;

    void v_ArraySolve(
        const Array<OneD, const Array<OneD, ND> > &Fwd,
        const Array<OneD, const Array<OneD, ND> > &Bwd,
              Array<OneD,       Array<OneD, ND> > &flux) final;
};

template<class T>
inline void RoeKernel(
    T& rhoL, T& rhouL, T& rhovL, T& rhowL, T& EL,
    T& rhoR, T& rhouR, T& rhovR, T& rhowR, T& ER,
    T& rhof, T& rhouf, T& rhovf, T& rhowf, T& Ef,
    NekDouble gamma)
{
    // Left and right velocities
    T uL = rhouL / rhoL;
    T vL = rhovL / rhoL;
    T wL = rhowL / rhoL;
    T uR = rhouR / rhoR;
    T vR = rhovR / rhoR;
    T wR = rhowR / rhoR;

    // Left and right pressures
    T pL = (gamma - 1.0) *
        (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
    T pR = (gamma - 1.0) *
        (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));

    // Left and right enthalpy
    T hL = (EL + pL) / rhoL;
    T hR = (ER + pR) / rhoR;

    // Square root of rhoL and rhoR.
    T srL  = sqrt(rhoL);
    T srR  = sqrt(rhoR);
    T srLR = srL + srR;

    // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
    T uRoe   = (srL * uL + srR * uR) / srLR;
    T vRoe   = (srL * vL + srR * vR) / srLR;
    T wRoe   = (srL * wL + srR * wR) / srLR;
    T hRoe   = (srL * hL + srR * hR) / srLR;
    T URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
    T cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));

    // Compute eigenvectors (equation 11.59).
    T k[5][5] = {
        {1., uRoe - cRoe, vRoe, wRoe, hRoe - uRoe * cRoe},
        {1., uRoe,        vRoe, wRoe, 0.5 * URoe},
        {0., 0.,           1.,    0.,    vRoe},
        {0., 0.,           0.,    1.,    wRoe},
        {1., uRoe+cRoe,  vRoe,  wRoe, hRoe + uRoe*cRoe}
    };

    // Calculate jumps \Delta u_i (defined preceding equation 11.67).
    T jump[5] = {
        rhoR  - rhoL,
        rhouR - rhouL,
        rhovR - rhovL,
        rhowR - rhowL,
        ER    - EL
    };

    // Define \Delta u_5 (equation 11.70).
    T jumpbar = jump[4] - (jump[2]-vRoe*jump[0])*vRoe -
        (jump[3]-wRoe*jump[0])*wRoe;

    // Compute wave amplitudes (equations 11.68, 11.69).
    T alpha[5];
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
    T uRoeAbs = abs(uRoe);
    T lambda[5] = {
        abs(uRoe - cRoe),
        uRoeAbs,
        uRoeAbs,
        uRoeAbs,
        abs(uRoe + cRoe)
    };

    // Finally perform summation (11.29).
    for (size_t i = 0; i < 5; ++i)
    {
        uRoeAbs = 0.5*alpha[i]*lambda[i];

        rhof  -= uRoeAbs*k[i][0];
        rhouf -= uRoeAbs*k[i][1];
        rhovf -= uRoeAbs*k[i][2];
        rhowf -= uRoeAbs*k[i][3];
        Ef    -= uRoeAbs*k[i][4];
    }
}


} // namespace

#endif
