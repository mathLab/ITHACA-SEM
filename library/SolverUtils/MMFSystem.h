///////////////////////////////////////////////////////////////////////////////
//
// File MMFSystem.h
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
// Description: MMF system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_MMFSYSTEM_H
#define NEKTAR_SOLVERUTILS_MMFSYSTEM_H

#include <SolverUtils/UnsteadySystem.h>
// #include <SolverUtils/Advection/Advection.h>

namespace Nektar
{
namespace SolverUtils
{

enum SurfaceType
{
    ePlane,
    eSphere,
    eTRSphere,
    eIrregular,
    eNonconvex,
    eCube,
    SIZE_SurfaceType
};

const char *const SurfaceTypeMap[] = {
    "Plane", "Sphere", "TRSphere", "Irregular", "Nonconvex", "Cube",
};

enum BoundaryCopyType
{
    eDirichlet,
    eNeumann,
    eFwdEQBwd,
    eFwdEQNegBwd,
    SIZE_BoundaryCopyType ///< Length of enum list
};

const char *const BoundaryCopyTypeMap[] = {
    "Dirichlet", "Neumann", "FwdEQBwd", "FwdEQNegBwd",
};

enum UpwindType
{
    eNotSet,        ///< flux not defined
    eAverage,       ///< averaged (or centred) flux
    eLaxFriedrich,  ///< Lax-Friedrich flux
    eUpwind,        ///  Upwind
    eRusanov,       ///< Rusanov flux
    eHLL,           ///< Harten-Lax-Leer flux
    eHLLC,          ///< Harten-Lax-Leer Contact wave flux
    SIZE_UpwindType ///< Length of enum list
};

const char *const UpwindTypeMap[] = {
    "NoSet", "Average", "LaxFriedrich", "Upwind", "Rusanov", "HLL", "HLLC",
};

enum TestMaxwellType
{
    eMaxwell1D,
    eTestMaxwell2DPEC,
    eTestMaxwell2DPECAVGFLUX,
    eTestMaxwell2DPMC,
    eMaxwell3D,
    eScatField1D,
    eScatField2D,
    eScatField3D,
    eTotField1D,
    eTotField2D,
    eTotField3D,
    eMaxwellSphere,
    eELF2DSurface,
    SIZE_TestMaxwellType ///< Length of enum list
};

const char *const TestMaxwellTypeMap[] = {
    "Maxwell1D",        "TestMaxwell2DPEC", "TestMaxwell2DPECAVGFLUX",
    "TestMaxwell2DPMC", "Maxwell3D",        "ScatField1D",
    "ScatField2D",      "ScatField3D",      "TotField1D",
    "TotField2D",       "TotField3D",       "MaxwellSphere",
    "ELF2DSurface",
};

enum PolType
{
    eTransMagnetic,
    eTransElectric,
    SIZE_PolType
};

const char *const PolTypeMap[] = {
    "TransMagnetic", "TransElectric",
};

enum IncType
{
    ePlaneWave,
    ePlaneWaveImag,
    eCylindricalWave,
    SIZE_IncType
};

const char *const IncTypeMap[] = {
    "PlaneWave", "PlaneWaveImag", "CylindricalWave",
};

/// A base class for PDEs which include an advection component
class MMFSystem : virtual public UnsteadySystem
{
public:
    NekDouble m_pi;

    int m_shapedim;

    SurfaceType m_surfaceType;
    UpwindType m_upwindType;

    TestMaxwellType m_TestMaxwellType;
    PolType m_PolType;
    IncType m_IncType;

    Array<OneD, NekDouble> m_MMFfactors;

    SOLVER_UTILS_EXPORT MMFSystem(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph);

    SOLVER_UTILS_EXPORT virtual ~MMFSystem();

    SOLVER_UTILS_EXPORT virtual void v_GenerateSummary(SummaryList &s);

    SOLVER_UTILS_EXPORT void MMFInitObject(
        const Array<OneD, const Array<OneD, NekDouble>> &Anisotropy,
        const int TangentXelem = -1);

    SOLVER_UTILS_EXPORT void CopyBoundaryTrace(
        const Array<OneD, const NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd,
        const BoundaryCopyType BDCopyType, const int var = 0,
        const std::string btype = "NoUserDefined");

protected:
    NekDouble m_alpha;

    NekDouble m_Incfreq;
    int m_SmoothFactor;
    NekDouble m_SFinit;

    Array<OneD, Array<OneD, NekDouble>> m_movingframes;
    Array<OneD, Array<OneD, NekDouble>> m_surfaceNormal;

    Array<OneD, Array<OneD, NekDouble>> m_ncdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> m_ncdotMFBwd;

    Array<OneD, Array<OneD, NekDouble>> m_nperpcdotMFFwd;
    Array<OneD, Array<OneD, NekDouble>> m_nperpcdotMFBwd;

    Array<OneD, Array<OneD, NekDouble>> m_DivMF;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_CurlMF;

    // MFdim \times spacedim \times npts
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_MFtraceFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_MFtraceBwd;

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ntimesMFFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ntimesMFBwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ntimes_ntimesMFFwd;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_ntimes_ntimesMFBwd;

    Array<OneD, Array<OneD, NekDouble>> m_ZimFwd;
    Array<OneD, Array<OneD, NekDouble>> m_ZimBwd;
    Array<OneD, Array<OneD, NekDouble>> m_YimFwd;
    Array<OneD, Array<OneD, NekDouble>> m_YimBwd;

    Array<OneD, Array<OneD, NekDouble>> m_epsvec;
    Array<OneD, Array<OneD, NekDouble>> m_muvec;

    Array<OneD, Array<OneD, NekDouble>> m_negepsvecminus1;
    Array<OneD, Array<OneD, NekDouble>> m_negmuvecminus1;

    // m_dedxi_cdot_e[m][j][n][] = de^m / d \xi^j \cdot e^n
    Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
        m_dedxi_cdot_e;

    SpatialDomains::GeomMMF m_MMFdir;

    Array<OneD, NekDouble> m_MFlength;

    void SetUpMovingFrames(
        const Array<OneD, const Array<OneD, NekDouble>> &Anisotropy,
        const int TangentXelem);

    void CheckMovingFrames(
        const Array<OneD, const Array<OneD, NekDouble>> &movingframes);

    SOLVER_UTILS_EXPORT void ComputencdotMF();

    SOLVER_UTILS_EXPORT void ComputeDivCurlMF();

    SOLVER_UTILS_EXPORT void ComputeMFtrace();

    SOLVER_UTILS_EXPORT void VectorDotProd(
        const Array<OneD, const Array<OneD, NekDouble>> &v1,
        const Array<OneD, const Array<OneD, NekDouble>> &v2,
        Array<OneD, NekDouble> &v3);

    SOLVER_UTILS_EXPORT void VectorCrossProd(
        const Array<OneD, const Array<OneD, NekDouble>> &v1,
        const Array<OneD, const Array<OneD, NekDouble>> &v2,
        Array<OneD, Array<OneD, NekDouble>> &v3);

    SOLVER_UTILS_EXPORT void VectorCrossProd(const Array<OneD, NekDouble> &v1,
                                             const Array<OneD, NekDouble> &v2,
                                             Array<OneD, NekDouble> &v3);

    SOLVER_UTILS_EXPORT void ComputeCurl(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> CartesianToMovingframes(
        const Array<OneD, const Array<OneD, NekDouble>> &uvec,
        unsigned int field);

    SOLVER_UTILS_EXPORT void DeriveCrossProductMF(
        Array<OneD, Array<OneD, NekDouble>> &CrossProductMF);
    SOLVER_UTILS_EXPORT void ComputeNtimesMF();

    SOLVER_UTILS_EXPORT void ComputeNtimesFz(
        const int dir, const Array<OneD, Array<OneD, NekDouble>> &Fwd,
        const Array<OneD, Array<OneD, NekDouble>> &Bwd,
        const Array<OneD, const NekDouble> &imFwd,
        const Array<OneD, const NekDouble> &imBwd,
        Array<OneD, NekDouble> &outarrayFwd,
        Array<OneD, NekDouble> &outarrayBwd);

    SOLVER_UTILS_EXPORT void ComputeNtimesF12(
        const Array<OneD, Array<OneD, NekDouble>> &Fwd,
        const Array<OneD, Array<OneD, NekDouble>> &Bwd,
        const Array<OneD, const NekDouble> &im1Fwd,
        const Array<OneD, const NekDouble> &im1Bwd,
        const Array<OneD, const NekDouble> &im2Fwd,
        const Array<OneD, const NekDouble> &im2Bwd,
        Array<OneD, NekDouble> &outarrayFwd,
        Array<OneD, NekDouble> &outarrayBwd);

    SOLVER_UTILS_EXPORT void ComputeNtimestimesdFz(
        const int dir, const Array<OneD, Array<OneD, NekDouble>> &Fwd,
        const Array<OneD, Array<OneD, NekDouble>> &Bwd,
        const Array<OneD, const NekDouble> &imFwd,
        const Array<OneD, const NekDouble> &imBwd,
        Array<OneD, NekDouble> &outarrayFwd,
        Array<OneD, NekDouble> &outarrayBwd);

    SOLVER_UTILS_EXPORT void ComputeNtimestimesdF12(
        const Array<OneD, Array<OneD, NekDouble>> &Fwd,
        const Array<OneD, Array<OneD, NekDouble>> &Bwd,
        const Array<OneD, const NekDouble> &im1Fwd,
        const Array<OneD, const NekDouble> &im1Bwd,
        const Array<OneD, const NekDouble> &im2Fwd,
        const Array<OneD, const NekDouble> &im2Bwd,
        Array<OneD, NekDouble> &outarrayFwd,
        Array<OneD, NekDouble> &outarrayBwd);

    SOLVER_UTILS_EXPORT void CartesianToSpherical(
        const NekDouble x0j, const NekDouble x1j, const NekDouble x2j,
        NekDouble &sin_varphi, NekDouble &cos_varphi, NekDouble &sin_theta,
        NekDouble &cos_theta);

    SOLVER_UTILS_EXPORT void ComputeZimYim(
        Array<OneD, Array<OneD, NekDouble>> &epsvec,
        Array<OneD, Array<OneD, NekDouble>> &muvec);

    SOLVER_UTILS_EXPORT void AdddedtMaxwell(
        const Array<OneD, const Array<OneD, NekDouble>> &physarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    SOLVER_UTILS_EXPORT void GetMaxwellFluxVector(
        const int var,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &flux);

    SOLVER_UTILS_EXPORT void GetMaxwellFlux1D(
        const int var,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &flux);

    SOLVER_UTILS_EXPORT void GetMaxwellFlux2D(
        const int var,
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &flux);

    SOLVER_UTILS_EXPORT void LaxFriedrichMaxwellFlux1D(
        Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
        Array<OneD, Array<OneD, NekDouble>> &numfluxBwd);

    SOLVER_UTILS_EXPORT void UpwindMaxwellFlux1D(
        Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
        Array<OneD, Array<OneD, NekDouble>> &numfluxBwd);

    SOLVER_UTILS_EXPORT void AverageMaxwellFlux1D(
        Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
        Array<OneD, Array<OneD, NekDouble>> &numfluxBwd);

    SOLVER_UTILS_EXPORT void NumericalMaxwellFlux(
        Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
        Array<OneD, Array<OneD, NekDouble>> &numfluxBwd,
        const NekDouble time = 0.0);

    SOLVER_UTILS_EXPORT void NumericalMaxwellFluxTM(
        Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
        Array<OneD, Array<OneD, NekDouble>> &numfluxBwd, const NekDouble time);

    SOLVER_UTILS_EXPORT void NumericalMaxwellFluxTE(
        Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
        Array<OneD, Array<OneD, NekDouble>> &numfluxBwd, const NekDouble time);

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> GetIncidentField(
        const int var, const NekDouble time);
    SOLVER_UTILS_EXPORT void Computedemdxicdote();

    SOLVER_UTILS_EXPORT NekDouble
    AvgInt(const Array<OneD, const NekDouble> &inarray);
    SOLVER_UTILS_EXPORT NekDouble
    AvgAbsInt(const Array<OneD, const NekDouble> &inarray);
    SOLVER_UTILS_EXPORT NekDouble
    AbsIntegral(const Array<OneD, const NekDouble> &inarray);
    SOLVER_UTILS_EXPORT NekDouble
    RootMeanSquare(const Array<OneD, const NekDouble> &inarray);
    SOLVER_UTILS_EXPORT NekDouble VectorAvgMagnitude(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray);

    SOLVER_UTILS_EXPORT void GramSchumitz(
        const Array<OneD, const Array<OneD, NekDouble>> &v1,
        const Array<OneD, const Array<OneD, NekDouble>> &v2,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        bool KeepTheMagnitude = true);

    SOLVER_UTILS_EXPORT void BubbleSort(Array<OneD, NekDouble> &refarray,
                                        Array<OneD, NekDouble> &sortarray);
};

// Shared pointer to an MMFSystem class
typedef std::shared_ptr<MMFSystem> MMFSystemSharedPtr;
}
}

#endif
