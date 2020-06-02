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
#include <CompressibleFlowSolver/RiemannSolvers/RoeSolverOpt.h>

#include <AVXOperators/AVXUtil.hpp>

namespace Nektar
{
std::string RoeSolverOpt::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "RoeOpt",
        RoeSolverOpt::create,
        "Roe Riemann solver opt");

RoeSolverOpt::RoeSolverOpt(const LibUtilities::SessionReaderSharedPtr& pSession)
    :
    CompressibleSolver(pSession)
{
    m_requiresRotation = false;
}

/// programmatic ctor
RoeSolverOpt::RoeSolverOpt():
CompressibleSolver()
{
    m_requiresRotation = false;
}

void RoeSolverOpt::v_Solve(
    const int                                 nDim,
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

    // get normal, vellocs
    ASSERTL1(CheckVectors("N"), "N not defined.");
    ASSERTL1(CheckAuxVec("vecLocs"), "vecLocs not defined.");
    const Array<OneD, const Array<OneD, NekDouble> > normals =
        m_vectors["N"]();
    const Array<OneD, const Array<OneD, NekDouble> > vecLocs =
        m_auxVec["vecLocs"]();

    const unsigned int vx = (int)vecLocs[0][0];
    const unsigned int vy = (int)vecLocs[0][1];
    const unsigned int vz = (int)vecLocs[0][2];

    // Generate matrices if they don't already exist.
    if (m_rotMat.num_elements() == 0)
    {
        GenerateRotationMatrices(normals);
    }

    // SIMD loop
    size_t i = 0;
    for (; i < sizeVec; i+=AVX::SIMD_WIDTH_SIZE)
    {

        // load scalars
        vec_t rhoL = &(fwd[0][i]);
        vec_t rhoR = &(bwd[0][i]);
        vec_t ER   = &(bwd[spaceDim+1][i]);
        vec_t EL   = &(fwd[spaceDim+1][i]);

        // 3D case only
        // load vectors left
        vec_t tmpIn[3], tmpOut[3];
        tmpIn[0] = &(fwd[1][i]);
        tmpIn[1] = &(fwd[2][i]);
        tmpIn[2] = &(fwd[3][i]);

        // load rotation matrix
        vec_t rotMat[9];
        for (size_t j = 0; j < 9; ++j)
        {
            rotMat[j] = &(m_rotMat[j][i]);
        }
        // rotateTo kernel Fwd
        rotateToNormalKernel(tmpIn, rotMat, tmpOut);

        vec_t rhouL = tmpOut[0];
        vec_t rhovL = tmpOut[1];
        vec_t rhowL = tmpOut[2];

        // load vectors right
        tmpIn[0] = &(bwd[1][i]);
        tmpIn[1] = &(bwd[2][i]);
        tmpIn[2] = &(bwd[3][i]);

        // rotateTo kernel Bwd
        rotateToNormalKernel(tmpIn, rotMat, tmpOut);

        vec_t rhouR = tmpOut[0];
        vec_t rhovR = tmpOut[1];
        vec_t rhowR = tmpOut[2];

        // Roe kernel
        vec_t rhof{}, rhouf{}, rhovf{}, rhowf{}, Ef{};
        RoeKernel(
            rhoL, rhouL, rhovL, rhowL, EL,
            rhoR, rhouR, rhovR, rhowR, ER,
            rhof, rhouf, rhovf, rhowf, Ef,
            gamma);

        // rotateFrom kernel
        tmpIn[0] = rhouf;
        tmpIn[1] = rhovf;
        tmpIn[2] = rhowf;
        rotateToNormalKernel(tmpIn, rotMat, tmpOut);

        // store scalar
        rhof.store_nts(&(flux[0][i]));
        Ef.store_nts(&(flux[nVars-1][i]));

        // store vector
        tmpOut[0].store_nts(&(flux[1][i]));
        tmpOut[1].store_nts(&(flux[2][i]));
        tmpOut[2].store_nts(&(flux[3][i]));

    }

}

} // namespace Nektar
