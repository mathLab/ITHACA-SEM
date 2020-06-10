///////////////////////////////////////////////////////////////////////////////
//
// File: RoeSolverSIMD.cpp
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
// Description: Roe Riemann solver using simd types.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/RoeSolver.h>
#include <CompressibleFlowSolver/RiemannSolvers/RoeSolverSIMD.h>

#include <LibUtilities/SimdLib/tinysimd.hpp>

namespace Nektar
{
std::string RoeSolverSIMD::solverName =
    SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
        "RoeOpt",
        RoeSolverSIMD::create,
        "Roe Riemann solver opt");

RoeSolverSIMD::RoeSolverSIMD(const LibUtilities::SessionReaderSharedPtr& pSession)
    :
    CompressibleSolver(pSession)
{
    m_requiresRotation = false;
}

/// programmatic ctor
RoeSolverSIMD::RoeSolverSIMD():
CompressibleSolver()
{
    m_requiresRotation = false;
}

void RoeSolverSIMD::v_Solve(
    const int                                        nDim,
    const Array<OneD, const Array<OneD, NekDouble> > &fwd,
    const Array<OneD, const Array<OneD, NekDouble> > &bwd,
          Array<OneD,       Array<OneD, NekDouble> > &flux)
{
    static auto gamma = m_params["gamma"]();
    static auto nVars = fwd.num_elements();
    static auto spaceDim = nVars-2;

    // 3D case only so far
    ASSERTL0(spaceDim == 3, "SIMD Roe implemented only for 3D case...");

    using namespace tinysimd;
    using vec_t = simd<NekDouble>;
    // using vec_t = typename tinysimd::abi::scalar<NekDouble>::type;

    // get limit of vectorizable chunk
    size_t sizeScalar = fwd[0].num_elements();
    size_t sizeVec = (sizeScalar / vec_t::width) * vec_t::width;

    // get normal, vellocs
    ASSERTL1(CheckVectors("N"), "N not defined.");
    // ASSERTL1(CheckAuxVec("vecLocs"), "vecLocs not defined.");
    const Array<OneD, const Array<OneD, NekDouble> > normals =
        m_vectors["N"]();
    // const Array<OneD, const Array<OneD, NekDouble> > vecLocs =
        // m_auxVec["vecLocs"]();

    // const unsigned int vx = (int)vecLocs[0][0];
    // const unsigned int vy = (int)vecLocs[0][1];
    // const unsigned int vz = (int)vecLocs[0][2];

    // Generate matrices if they don't already exist.
    if (m_rotMat.num_elements() == 0)
    {
        GenerateRotationMatrices(normals);
    }

    // SIMD loop
    size_t i = 0;
    for (; i < sizeVec; i+=vec_t::width)
    {
        // load scalars
        vec_t rhoL, rhoR, ER, EL;
        rhoL.load(&(fwd[0][i]), is_not_aligned);
        rhoR.load(&(bwd[0][i]), is_not_aligned);
        ER.load(&(bwd[spaceDim+1][i]), is_not_aligned);
        EL.load(&(fwd[spaceDim+1][i]), is_not_aligned);

        // load vectors left
        vec_t tmpIn[3], tmpOut[3];
        tmpIn[0].load(&(fwd[1][i]), is_not_aligned);
        tmpIn[1].load(&(fwd[2][i]), is_not_aligned);
        tmpIn[2].load(&(fwd[3][i]), is_not_aligned);

        // load rotation matrix
        vec_t rotMat[9];
        for (size_t j = 0; j < 9; ++j)
        {
            rotMat[j].load(&(m_rotMat[j][i]), is_not_aligned);
        }

        // rotateTo kernel Fwd
        rotateToNormalKernel(tmpIn, rotMat, tmpOut);

        vec_t rhouL = tmpOut[0];
        vec_t rhovL = tmpOut[1];
        vec_t rhowL = tmpOut[2];

        // load vectors right
        tmpIn[0].load(&(bwd[1][i]), is_not_aligned);
        tmpIn[1].load(&(bwd[2][i]), is_not_aligned);
        tmpIn[2].load(&(bwd[3][i]), is_not_aligned);

        // rotateTo kernel Bwd
        rotateToNormalKernel(tmpIn, rotMat, tmpOut);

        vec_t rhouR = tmpOut[0];
        vec_t rhovR = tmpOut[1];
        vec_t rhowR = tmpOut[2];

        // Roe kernel
        vec_t rhof{}, Ef{};
        RoeKernel(
            rhoL, rhouL, rhovL, rhowL, EL,
            rhoR, rhouR, rhovR, rhowR, ER,
            rhof, tmpIn[0], tmpIn[1], tmpIn[2], Ef,
            gamma);

        // rotateFrom kernel
        rotateFromNormalKernel(tmpIn, rotMat, tmpOut);

        // store scalar
        #if 0
        // nts
        rhof.store_nts(&(flux[0][i]), is_not_reused);
        Ef.store_nts(&(flux[nVars-1][i]), is_not_reused);
        #else
        //unaligned
        rhof.store(&(flux[0][i]), is_not_aligned);
        Ef.store(&(flux[nVars-1][i]), is_not_aligned);
        #endif

        // store vector 3D only
        #if 0
        // nts
        // tmpOut[0].store_nts(&(flux[1][i]), is_not_reused);
        // tmpOut[1].store_nts(&(flux[2][i]), is_not_reused);
        // tmpOut[2].store_nts(&(flux[3][i]), is_not_reused);
        #else
        // unaligned
        tmpOut[0].store(&(flux[1][i]), is_not_aligned);
        tmpOut[1].store(&(flux[2][i]), is_not_aligned);
        tmpOut[2].store(&(flux[3][i]), is_not_aligned);
        #endif
    }

    // spillover loop
    for (; i < sizeScalar; ++i)
    {
        // load scalars
        NekDouble rhoL = fwd[0][i];
        NekDouble rhoR = bwd[0][i];
        NekDouble EL   = fwd[spaceDim+1][i];
        NekDouble ER   = bwd[spaceDim+1][i];

        // 3D case only
        // load vectors left
        NekDouble tmpIn[3], tmpOut[3];
        tmpIn[0] = fwd[1][i];
        tmpIn[1] = fwd[2][i];
        tmpIn[2] = fwd[3][i];

        // load rotation matrix
        NekDouble rotMat[9];
        for (size_t j = 0; j < 9; ++j)
        {
            rotMat[j] = m_rotMat[j][i];
        }

        // rotateTo kernel Fwd
        rotateToNormalKernel(tmpIn, rotMat, tmpOut);

        NekDouble rhouL = tmpOut[0];
        NekDouble rhovL = tmpOut[1];
        NekDouble rhowL = tmpOut[2];

        // load vectors right
        tmpIn[0] = bwd[1][i];
        tmpIn[1] = bwd[2][i];
        tmpIn[2] = bwd[3][i];

        // rotateTo kernel Bwd
        rotateToNormalKernel(tmpIn, rotMat, tmpOut);

        NekDouble rhouR = tmpOut[0];
        NekDouble rhovR = tmpOut[1];
        NekDouble rhowR = tmpOut[2];

        // Roe kernel
        NekDouble rhof{}, Ef{};
        RoeKernel(
            rhoL, rhouL, rhovL, rhowL, EL,
            rhoR, rhouR, rhovR, rhowR, ER,
            rhof, tmpIn[0], tmpIn[1], tmpIn[2], Ef,
            gamma);

        // rotateFrom kernel
        rotateFromNormalKernel(tmpIn, rotMat, tmpOut);

        // store scalar
        flux[0][i] = rhof;
        flux[nVars-1][i] = Ef;

        // store vector 3D only
        flux[1][i] = tmpOut[0];
        flux[2][i] = tmpOut[1];
        flux[3][i] = tmpOut[2];

    }
}

} // namespace Nektar
