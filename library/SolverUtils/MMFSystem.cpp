///////////////////////////////////////////////////////////////////////////////
//
// File: MMFSystem.cpp
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
// Description: Base class for MMF systems.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/MMFSystem.h>

namespace Nektar
{
namespace SolverUtils
{

MMFSystem::MMFSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr& pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}

MMFSystem::~MMFSystem()
{
}

void MMFSystem::MMFInitObject(
    const Array<OneD, const Array<OneD, NekDouble>> &Anisotropy,
    const int TangentXelem)
{
    m_pi       = 3.14159265358979323846;
    m_shapedim = m_fields[0]->GetShapeDimension();

    ASSERTL0(m_spacedim == 3, "Space Dimension should be 3.");

    // Factor for Numerical Flux
    m_session->LoadParameter("alpha", m_alpha, 1.0);

    // Factor for Numerical Flux
    m_session->LoadParameter("Incfreq", m_Incfreq, 1.0);

    // SmoothFactor
    m_session->LoadParameter("SmoothFactor", m_SmoothFactor, 1);
    m_session->LoadParameter("SFinit", m_SFinit, 0.0);

    // Define SurfaceType
    if (m_session->DefinesSolverInfo("SURFACETYPE"))
    {
        std::string SurfaceTypeStr = m_session->GetSolverInfo("SURFACETYPE");
        int i;
        for (i = 0; i < (int)SIZE_SurfaceType; ++i)
        {
            if (boost::iequals(SurfaceTypeMap[i], SurfaceTypeStr))
            {
                m_surfaceType = (SurfaceType)i;
                break;
            }
        }
    }
    else
    {
        m_surfaceType = (SurfaceType)0;
    }

    // if discontinuous Galerkin determine numerical flux to use
    for (int i = 0; i < (int)SIZE_UpwindType; ++i)
    {
        bool match;
        m_session->MatchSolverInfo("UPWINDTYPE", UpwindTypeMap[i], match,
                                   false);
        if (match)
        {
            m_upwindType = (UpwindType)i;
            break;
        }
    }

    // SetUpMovingFrames: To generate m_movingframes
    SetUpMovingFrames(Anisotropy, TangentXelem);

    // Derive movingfrmaes at interfaces
    // Get: m_MFtraceFwd and m_MFtraceBwd
    ComputeMFtrace();

    // Check Movingframes and surfraceNormal
    // Get: m_ncdotMFFwd,m_ncdotMFBwd,m_nperpcdotMFFwd,m_nperpcdotMFBwd
    ComputencdotMF();

    // Compute nabla cdot m_movingframes
    // Get: m_DivMF, m_CurlMF
    ComputeDivCurlMF();
}

void MMFSystem::SetUpMovingFrames(
    const Array<OneD, const Array<OneD, NekDouble>> &Anisotropy,
    const int TangentXelem)
{

    int nq    = m_fields[0]->GetNpoints();
    int MFdim = 3;

    // Construct The Moving Frames
    m_movingframes = Array<OneD, Array<OneD, NekDouble>>(MFdim);

    for (int j = 0; j < MFdim; ++j)
    {
        m_movingframes[j] = Array<OneD, NekDouble>(m_spacedim * nq, 0.0);
    }

    // Read MMF Geom Info
    std::string conn = "TangentX";
    m_MMFfactors     = Array<OneD, NekDouble>(4);

    m_session->LoadSolverInfo("MMFDir", conn, "LOCAL");

    // (x-x_0)^2/a^2 + (y-y_0)^2/b^2 = 1
    // factors[0] = a
    // factors[1] = b
    // factors[2] = x_0
    // factors[3] = y_0
    m_session->LoadParameter("MMFCircAxisX", m_MMFfactors[0], 1.0);
    m_session->LoadParameter("MMFCircAxisY", m_MMFfactors[1], 1.0);
    m_session->LoadParameter("MMFCircCentreX", m_MMFfactors[2], 0.0);
    m_session->LoadParameter("MMFCircCentreY", m_MMFfactors[3], 0.0);

    if (conn == "TangentX")
        m_MMFdir = SpatialDomains::eTangentX;
    if (conn == "TangentY")
        m_MMFdir = SpatialDomains::eTangentY;
    if (conn == "TangentXY")
        m_MMFdir = SpatialDomains::eTangentXY;
    if (conn == "TangentZ")
        m_MMFdir = SpatialDomains::eTangentZ;
    if (conn == "TangentCircular")
        m_MMFdir = SpatialDomains::eTangentCircular;
    if (conn == "TangentIrregular")
        m_MMFdir = SpatialDomains::eTangentIrregular;
    if (conn == "TangentNonconvex")
        m_MMFdir = SpatialDomains::eTangentNonconvex;
    if (conn == "LOCAL")
        m_MMFdir = SpatialDomains::eLOCAL;

    /*for (int j=0; j<m_shapedim; j++)
    {
        for (int k=0; k<m_spacedim; k++)
        {
            std::cout << "before ehsan " << nq << "\t"<< m_shapedim << "\t" <<
    m_spacedim << "\t" <<"m_movingframes " << m_movingframes[j][k*nq] <<
    std::endl;
        }
    }*/
    // Get Tangetn vectors from GeomFactors2D, Orthonormalized = true
    m_fields[0]->GetMovingFrames(m_MMFdir, m_MMFfactors, m_movingframes);
    /*  for (int j=0; j<m_shapedim; j++)
      {
          for (int k=0; k<m_spacedim; k++)
          {
                std::cout << "ehsan " << nq << "\t"<< m_shapedim << "\t" <<
      m_spacedim << "\t" <<"m_movingframes " << m_movingframes[j][k*nq] <<
      std::endl;
          }
      }*/
    // Align the tangentX direction after TangentXelem
    if (TangentXelem > 0)
    {
        Array<OneD, NekDouble> tmp(nq);
        m_fields[0]->GenerateElementVector(TangentXelem, 1.0, 0.0, tmp);
        for (int j = 0; j < m_shapedim; j++)
        {
            for (int k = 0; k < m_spacedim; k++)
            {
                Vmath::Vmul(nq, &tmp[0], 1, &m_movingframes[j][k * nq], 1,
                            &m_movingframes[j][k * nq], 1);
            }
        }

        m_fields[0]->GenerateElementVector(TangentXelem, 0.0, 1.0, tmp);
        Vmath::Vadd(nq, &tmp[0], 1, &m_movingframes[0][0 * nq], 1,
                    &m_movingframes[0][0 * nq], 1);
        Vmath::Vadd(nq, &tmp[0], 1, &m_movingframes[1][1 * nq], 1,
                    &m_movingframes[1][1 * nq], 1);

        int indxtmp = Vmath::Imax(nq, tmp, 1);
        std::cout << "*** MF in PML Region is aligned as MF1 = ( "
                  << m_movingframes[0][indxtmp] << " , "
                  << m_movingframes[0][nq + indxtmp] << " , "
                  << m_movingframes[0][2 * nq + indxtmp] << " ) "
                  << ", MF2 = ( " << m_movingframes[1][indxtmp] << " , "
                  << m_movingframes[1][nq + indxtmp] << " , "
                  << m_movingframes[1][2 * nq + indxtmp] << " ) " << std::endl;
    }

    // Multiply Anisotropy to movingframes
    for (int j = 0; j < m_shapedim; ++j)
    {
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vmul(nq, &Anisotropy[j][0], 1, &m_movingframes[j][k * nq], 1,
                        &m_movingframes[j][k * nq], 1);
        }
    }

    // Test the moving frames
    CheckMovingFrames(m_movingframes);
}

void MMFSystem::CheckMovingFrames(
    const Array<OneD, const Array<OneD, NekDouble>> &movingframes)
{
    NekDouble t1x, t1y, t1z, t2x, t2y, t2z, t3x, t3y, t3z;
    NekDouble dot12 = 0.0, dot23 = 0.0, dot31 = 0.0;
    NekDouble Tol = 0.0001;

    int nq = m_fields[0]->GetNpoints();

    for (int i = 0; i < nq; ++i)
    {
        t1x = movingframes[0][i];
        t1y = movingframes[0][i + nq];
        t1z = movingframes[0][i + 2 * nq];

        t2x = movingframes[1][i];
        t2y = movingframes[1][i + nq];
        t2z = movingframes[1][i + 2 * nq];

        t3x = movingframes[2][i];
        t3y = movingframes[2][i + nq];
        t3z = movingframes[2][i + 2 * nq];

        dot12 = t1x * t2x + t1y * t2y + t1z * t2z;
        dot23 = t2x * t3x + t2y * t3y + t2z * t3z;
        dot31 = t3x * t1x + t3y * t1y + t3z * t1z;
    }

    std::cout << "======================================================"
              << std::endl;
    std::cout << "======================================================"
              << std::endl;
    std::cout << "*** The first moving frame is alinged along"
              << SpatialDomains::GeomMMFMap[m_MMFdir] << std::endl;

    Array<OneD, NekDouble> tmpx(nq), tmpy(nq), tmpz(nq);

    Vmath::Vcopy(nq, &movingframes[0][0], 1, &tmpx[0], 1);
    Vmath::Vcopy(nq, &movingframes[0][nq], 1, &tmpy[0], 1);
    Vmath::Vcopy(nq, &movingframes[0][2 * nq], 1, &tmpz[0], 1);
    std::cout << nq << " , "
              << "*** Avg MF1 = ( " << AvgAbsInt(tmpx) << " , "
              << AvgAbsInt(tmpy) << " , " << AvgAbsInt(tmpz) << " ) "
              << std::endl;

    Vmath::Vcopy(nq, &movingframes[1][0], 1, &tmpx[0], 1);
    Vmath::Vcopy(nq, &movingframes[1][nq], 1, &tmpy[0], 1);
    Vmath::Vcopy(nq, &movingframes[1][2 * nq], 1, &tmpz[0], 1);
    std::cout << "*** Avg MF2 = ( " << AvgAbsInt(tmpx) << " , "
              << AvgAbsInt(tmpy) << " , " << AvgAbsInt(tmpz) << " ) "
              << std::endl;

    if (m_shapedim == 3)
    {
        Vmath::Vcopy(nq, &movingframes[2][0], 1, &tmpx[0], 1);
        Vmath::Vcopy(nq, &movingframes[2][nq], 1, &tmpy[0], 1);
        Vmath::Vcopy(nq, &movingframes[2][2 * nq], 1, &tmpz[0], 1);
        std::cout << "*** Avg MF3 = ( " << AvgAbsInt(tmpx) << " , "
                  << AvgAbsInt(tmpy) << " , " << AvgAbsInt(tmpz) << " ) "
                  << std::endl;
    }

    if ((fabs(dot12) + fabs(dot23) + fabs(dot31)) < Tol)
    {
        std::cout << "*** Moving frames are Orthogonal" << std::endl;
    }

    else
    {
        std::cout << "*** Moving frames are NOT Orthogonal" << std::endl;
    }

    Array<OneD, NekDouble> tmp;
    for (int j = 0; j < m_shapedim; ++j)
    {
        tmp = Array<OneD, NekDouble>(nq, 0.0);
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nq, &movingframes[j][k * nq], 1,
                         &movingframes[j][k * nq], 1, &tmp[0], 1, &tmp[0], 1);
        }
        Vmath::Vsqrt(nq, tmp, 1, tmp, 1);
        std::cout << "*** Avg. Magnitude of MF" << j << " = " << AvgAbsInt(tmp)
                  << std::endl;
    }
}

void MMFSystem::ComputencdotMF()
{
    int nq              = m_fields[0]->GetNpoints();
    int nTracePointsTot = GetTraceNpoints();

    // Compute MFjFwd and MFjBwd
    Array<OneD, NekDouble> tmp(nq);

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceFwd(m_shapedim);
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtraceBwd(m_shapedim);
    Array<OneD, Array<OneD, NekDouble>> SurfaceNormalFwd;
    Array<OneD, Array<OneD, NekDouble>> SurfaceNormalBwd;

    for (int j = 0; j < m_shapedim; ++j)
    {
        MFtraceFwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        MFtraceBwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

        SurfaceNormalFwd = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        SurfaceNormalBwd = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

        for (int k = 0; k < m_spacedim; ++k)
        {
            MFtraceFwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);
            MFtraceBwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);

            SurfaceNormalFwd[k] = Array<OneD, NekDouble>(nTracePointsTot);
            SurfaceNormalBwd[k] = Array<OneD, NekDouble>(nTracePointsTot);

            Vmath::Vcopy(nq, &m_movingframes[j][k * nq], 1, &tmp[0], 1);

            m_fields[0]->GetFwdBwdTracePhys(tmp, MFtraceFwd[j][k],
                                            MFtraceBwd[j][k]);

            CopyBoundaryTrace(MFtraceFwd[j][k], MFtraceBwd[j][k], eFwdEQBwd);
        }
    }

    VectorCrossProd(MFtraceFwd[0], MFtraceFwd[1], SurfaceNormalFwd);
    VectorCrossProd(MFtraceBwd[0], MFtraceBwd[1], SurfaceNormalBwd);

    // Compute n \times e^i
    m_ncdotMFFwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    m_ncdotMFBwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        m_ncdotMFFwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        m_ncdotMFBwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);

        VectorDotProd(m_traceNormals, MFtraceFwd[j], m_ncdotMFFwd[j]);
        VectorDotProd(m_traceNormals, MFtraceBwd[j], m_ncdotMFBwd[j]);
    }

    // Compute n^{\perp} \times e^i
    Array<OneD, Array<OneD, NekDouble>> SurfaceNormal(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> Tracevector(m_spacedim);
    for (int k = 0; k < m_spacedim; k++)
    {
        SurfaceNormal[k] = Array<OneD, NekDouble>(nTracePointsTot);
        Tracevector[k]   = Array<OneD, NekDouble>(nTracePointsTot);

        Vmath::Vcopy(nTracePointsTot, &m_MFtraceFwd[2][k][0], 1,
                     &SurfaceNormal[k][0], 1);
    }

    VectorCrossProd(m_traceNormals, SurfaceNormal, Tracevector);

    m_nperpcdotMFFwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    m_nperpcdotMFBwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
    for (int j = 0; j < m_shapedim; ++j)
    {
        m_nperpcdotMFFwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
        m_nperpcdotMFBwd[j] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);

        VectorDotProd(Tracevector, MFtraceFwd[j], m_nperpcdotMFFwd[j]);
        VectorDotProd(Tracevector, MFtraceBwd[j], m_nperpcdotMFBwd[j]);
    }

    if (m_shapedim == 2)
    {
        std::cout << "*** m_ncdotMFFwd = ( " << RootMeanSquare(m_ncdotMFFwd[0])
                  << " , " << RootMeanSquare(m_ncdotMFFwd[1]) << " ) "
                  << std::endl;
        std::cout << "*** m_ncdotMFBwd = ( " << RootMeanSquare(m_ncdotMFBwd[0])
                  << " , " << RootMeanSquare(m_ncdotMFBwd[1]) << " ) "
                  << std::endl;

        std::cout << "*** m_nperpcdotMFFwd = ( "
                  << RootMeanSquare(m_nperpcdotMFFwd[0]) << " , "
                  << RootMeanSquare(m_nperpcdotMFFwd[1]) << " ) " << std::endl;
        std::cout << "*** m_nperpcdotMFBwd = ( "
                  << RootMeanSquare(m_nperpcdotMFBwd[0]) << " , "
                  << RootMeanSquare(m_nperpcdotMFBwd[1]) << " ) " << std::endl;
    }

    else if (m_shapedim == 3)
    {
        std::cout << "*** m_ncdotMFFwd = ( "
                  << Vmath::Vsum(nTracePointsTot, m_ncdotMFFwd[0], 1) << " , "
                  << Vmath::Vsum(nTracePointsTot, m_ncdotMFFwd[1], 1) << " , "
                  << Vmath::Vsum(nTracePointsTot, m_ncdotMFFwd[2], 1) << " ) "
                  << std::endl;
        std::cout << "*** m_ncdotMFBwd = ( "
                  << Vmath::Vsum(nTracePointsTot, m_ncdotMFBwd[0], 1) << " , "
                  << Vmath::Vsum(nTracePointsTot, m_ncdotMFBwd[1], 1) << " , "
                  << Vmath::Vsum(nTracePointsTot, m_ncdotMFBwd[2], 1) << " ) "
                  << std::endl;

        std::cout << "*** m_nperpcdotMFFwd = ( "
                  << Vmath::Vsum(nTracePointsTot, m_nperpcdotMFFwd[0], 1)
                  << " , "
                  << Vmath::Vsum(nTracePointsTot, m_nperpcdotMFFwd[1], 1)
                  << " , "
                  << Vmath::Vsum(nTracePointsTot, m_nperpcdotMFFwd[2], 1)
                  << " ) " << std::endl;
        std::cout << "*** m_nperpcdotMFBwd = ( "
                  << Vmath::Vsum(nTracePointsTot, m_nperpcdotMFBwd[0], 1)
                  << " , "
                  << Vmath::Vsum(nTracePointsTot, m_nperpcdotMFBwd[1], 1)
                  << " , "
                  << Vmath::Vsum(nTracePointsTot, m_nperpcdotMFBwd[2], 1)
                  << " ) " << std::endl;
    }
}

void MMFSystem::ComputeDivCurlMF()
{
    int nq     = m_fields[0]->GetNpoints();
    int MMFdim = 3;

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Dtmp(nq);

    m_DivMF = Array<OneD, Array<OneD, NekDouble>>(MMFdim);
    for (int j = 0; j < MMFdim; ++j)
    {
        m_DivMF[j] = Array<OneD, NekDouble>(nq, 0.0);
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vcopy(nq, &m_movingframes[j][k * nq], 1, &tmp[0], 1);

            m_fields[0]->PhysDeriv(k, tmp, Dtmp);
            Vmath::Vadd(nq, &Dtmp[0], 1, &m_DivMF[j][0], 1, &m_DivMF[j][0], 1);
        }
    }

    std::cout << "*** Divergence of MF1 = " << AvgInt(m_DivMF[0])
              << ", MF2 = " << AvgInt(m_DivMF[1])
              << ", MF3 = " << AvgInt(m_DivMF[2]) << std::endl;

    // Compute Curl of MF: CurlMF[i][j] = (\nabla \times e^i) cdot e^j
    m_CurlMF = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(MMFdim);
    for (int i = 0; i < MMFdim; ++i)
    {
        m_CurlMF[i] = Array<OneD, Array<OneD, NekDouble>>(MMFdim);

        for (int j = 0; j < MMFdim; ++j)
        {
            m_CurlMF[i][j] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    Array<OneD, Array<OneD, NekDouble>> MFtmp(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> CurlMFtmp(m_spacedim);

    for (int i = 0; i < m_spacedim; ++i)
    {
        MFtmp[i]     = Array<OneD, NekDouble>(nq);
        CurlMFtmp[i] = Array<OneD, NekDouble>(nq);
    }

    for (int dir = 0; dir < MMFdim; dir++)
    {
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vcopy(nq, &m_movingframes[dir][i * nq], 1, &MFtmp[i][0], 1);
        }

        ComputeCurl(MFtmp, CurlMFtmp);

        for (int j = 0; j < MMFdim; ++j)
        {
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nq, &m_movingframes[j][i * nq], 1,
                             &CurlMFtmp[i][0], 1, &m_CurlMF[dir][j][0], 1,
                             &m_CurlMF[dir][j][0], 1);
            }
        }
    }

    std::cout << "*** Curl of MF1 = ( " << AvgInt(m_CurlMF[0][0]) << " , "
              << AvgInt(m_CurlMF[0][1]) << " , " << AvgInt(m_CurlMF[0][2])
              << " ) " << std::endl;
    std::cout << "*** Curl of MF2 = ( " << AvgInt(m_CurlMF[1][0]) << " , "
              << AvgInt(m_CurlMF[1][1]) << " , " << AvgInt(m_CurlMF[1][2])
              << " ) " << std::endl;
    std::cout << "*** Curl of MF3 = ( " << AvgInt(m_CurlMF[2][0]) << " , "
              << AvgInt(m_CurlMF[2][1]) << " , " << AvgInt(m_CurlMF[2][2])
              << " ) " << std::endl;
}

void MMFSystem::ComputeMFtrace()
{
    int MFdim = 3;

    int nq              = m_fields[0]->GetNpoints();
    int nTraceNumPoints = GetTraceTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> Fwdtmp(nq);
    Array<OneD, NekDouble> Bwdtmp(nq);

    m_MFtraceFwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(MFdim);
    m_MFtraceBwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(MFdim);

    for (int j = 0; j < MFdim; ++j)
    {
        m_MFtraceFwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        m_MFtraceBwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    }

    // m_MFtraceFwd[0] = e^1_{Fwd}, m_MFtraceFwd[1] = e^2_{Fwd}
    for (int j = 0; j < MFdim; ++j)
    {
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_MFtraceFwd[j][i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_MFtraceBwd[j][i] = Array<OneD, NekDouble>(nTraceNumPoints);

            Vmath::Vcopy(nq, &m_movingframes[j][i * nq], 1, &tmp[0], 1);

            m_fields[0]->GetFwdBwdTracePhys(tmp, Fwdtmp, Bwdtmp);

            CopyBoundaryTrace(Fwdtmp, Bwdtmp, eFwdEQBwd);

            Vmath::Vcopy(nTraceNumPoints, &Fwdtmp[0], 1, &m_MFtraceFwd[j][i][0],
                         1);
            Vmath::Vcopy(nTraceNumPoints, &Bwdtmp[0], 1, &m_MFtraceBwd[j][i][0],
                         1);
        }
    }

    std::cout << "*** MFtraceFwd = ( " << VectorAvgMagnitude(m_MFtraceFwd[0])
              << " , " << VectorAvgMagnitude(m_MFtraceFwd[1]) << " , "
              << VectorAvgMagnitude(m_MFtraceFwd[2]) << " ) " << std::endl;
    std::cout << "*** MFtraceBwd = ( " << VectorAvgMagnitude(m_MFtraceBwd[0])
              << " , " << VectorAvgMagnitude(m_MFtraceBwd[1]) << " , "
              << VectorAvgMagnitude(m_MFtraceBwd[2]) << " ) " << std::endl;
}

void MMFSystem::DeriveCrossProductMF(
    Array<OneD, Array<OneD, NekDouble>> &CrossProductMF)
{
    int MFdim = 3;
    int nq    = GetTotPoints();

    CrossProductMF = Array<OneD, Array<OneD, NekDouble>>(MFdim);

    Array<OneD, Array<OneD, NekDouble>> MF1tmp(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> MF2tmp(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> MF3tmp(m_spacedim);
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> MFtmpCurl(MFdim);
    for (int j = 0; j < MFdim; ++j)
    {
        CrossProductMF[j] = Array<OneD, NekDouble>(nq * m_spacedim);
        MFtmpCurl[j]      = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int k = 0; k < m_spacedim; ++k)
        {
            MFtmpCurl[j][k] = Array<OneD, NekDouble>(nq);
        }
    }

    for (int k = 0; k < m_spacedim; ++k)
    {
        MF1tmp[k] = Array<OneD, NekDouble>(nq);
        MF2tmp[k] = Array<OneD, NekDouble>(nq);
        MF3tmp[k] = Array<OneD, NekDouble>(nq);

        Vmath::Vcopy(nq, &m_movingframes[0][k * nq], 1, &MF1tmp[k][0], 1);
        Vmath::Vcopy(nq, &m_movingframes[1][k * nq], 1, &MF2tmp[k][0], 1);
        Vmath::Vcopy(nq, &m_movingframes[2][k * nq], 1, &MF3tmp[k][0], 1);
    }

    VectorCrossProd(MF3tmp, MF1tmp, MFtmpCurl[0]);
    VectorCrossProd(MF2tmp, MF3tmp, MFtmpCurl[1]);
    VectorCrossProd(MF1tmp, MF2tmp, MFtmpCurl[2]);

    for (int j = 0; j < MFdim; ++j)
    {
        for (int k = 0; k < m_spacedim; ++k)
        {
            Vmath::Vcopy(nq, &MFtmpCurl[j][k][0], 1, &CrossProductMF[j][k * nq],
                         1);
        }
    }
}

void MMFSystem::ComputeNtimesMF()
{
    int MFdim           = 3;
    int nTracePointsTot = GetTraceNpoints();

    m_ntimesMFFwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(MFdim);
    m_ntimesMFBwd = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(MFdim);
    m_ntimes_ntimesMFFwd =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(MFdim);
    m_ntimes_ntimesMFBwd =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(MFdim);
    for (int j = 0; j < MFdim; ++j)
    {
        m_ntimesMFFwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        m_ntimesMFBwd[j] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        m_ntimes_ntimesMFFwd[j] =
            Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        m_ntimes_ntimesMFBwd[j] =
            Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

        for (int k = 0; k < m_spacedim; ++k)
        {
            m_ntimesMFFwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);
            m_ntimesMFBwd[j][k] = Array<OneD, NekDouble>(nTracePointsTot);
            m_ntimes_ntimesMFFwd[j][k] =
                Array<OneD, NekDouble>(nTracePointsTot);
            m_ntimes_ntimesMFBwd[j][k] =
                Array<OneD, NekDouble>(nTracePointsTot);
        }

        VectorCrossProd(m_traceNormals, m_MFtraceFwd[j], m_ntimesMFFwd[j]);
        VectorCrossProd(m_traceNormals, m_MFtraceBwd[j], m_ntimesMFBwd[j]);
        VectorCrossProd(m_traceNormals, m_ntimesMFFwd[j],
                        m_ntimes_ntimesMFFwd[j]);
        VectorCrossProd(m_traceNormals, m_ntimesMFBwd[j],
                        m_ntimes_ntimesMFBwd[j]);
    }

    std::cout << "*** m_ntimesMFFwd = ( " << VectorAvgMagnitude(m_ntimesMFFwd[0])
              << " , " << VectorAvgMagnitude(m_ntimesMFFwd[1]) << " , "
              << VectorAvgMagnitude(m_ntimesMFFwd[2]) << " ) " << std::endl;
    std::cout << "*** m_ntimesMFBwd = ( " << VectorAvgMagnitude(m_ntimesMFBwd[0])
              << " , " << VectorAvgMagnitude(m_ntimesMFBwd[1]) << " , "
              << VectorAvgMagnitude(m_ntimesMFBwd[2]) << " ) " << std::endl;
    std::cout << "*** m_ntimes_ntimesMFFwd = ( "
              << VectorAvgMagnitude(m_ntimes_ntimesMFFwd[0]) << " , "
              << VectorAvgMagnitude(m_ntimes_ntimesMFFwd[1]) << " , "
              << VectorAvgMagnitude(m_ntimes_ntimesMFFwd[2]) << " ) "
              << std::endl;
    std::cout << "*** m_ntimes_ntimesMFBwd = ( "
              << VectorAvgMagnitude(m_ntimes_ntimesMFBwd[0]) << " , "
              << VectorAvgMagnitude(m_ntimes_ntimesMFBwd[1]) << " , "
              << VectorAvgMagnitude(m_ntimes_ntimesMFBwd[2]) << " ) "
              << std::endl;
}

void MMFSystem::VectorDotProd(
    const Array<OneD, const Array<OneD, NekDouble>> &v1,
    const Array<OneD, const Array<OneD, NekDouble>> &v2,
    Array<OneD, NekDouble> &v3)
{
    int coordim = v1.size();
    int nq      = v1[0].size();

    v3 = Array<OneD, NekDouble>(nq, 0.0);
    for (int i = 0; i < coordim; ++i)
    {
        Vmath::Vvtvp(nq, &v1[i][0], 1, &v2[i][0], 1, &v3[0], 1, &v3[0], 1);
    }
}

/**
 * Computes the vector cross-product in 3D of \a v1 and \a v2, storing
 * the result in \a v3.
 * @param   v1          First input vector.
 * @param   v2          Second input vector.
 * @param   v3          Output vector computed to be orthogonal to
 *                      both \a v1 and \a v2.
 */
void MMFSystem::VectorCrossProd(
    const Array<OneD, const Array<OneD, NekDouble>> &v1,
    const Array<OneD, const Array<OneD, NekDouble>> &v2,
    Array<OneD, Array<OneD, NekDouble>> &v3)
{
    ASSERTL0(v1.size() == 3, "Input 1 has dimension not equal to 3.");
    ASSERTL0(v2.size() == 3, "Input 2 has dimension not equal to 3.");
    ASSERTL0(v3.size() == 3,
             "Output vector has dimension not equal to 3.");

    int nq = v1[0].size();
    Array<OneD, NekDouble> temp(nq);

    Vmath::Vmul(nq, v1[2], 1, v2[1], 1, temp, 1);
    Vmath::Vvtvm(nq, v1[1], 1, v2[2], 1, temp, 1, v3[0], 1);

    Vmath::Vmul(nq, v1[0], 1, v2[2], 1, temp, 1);
    Vmath::Vvtvm(nq, v1[2], 1, v2[0], 1, temp, 1, v3[1], 1);

    Vmath::Vmul(nq, v1[1], 1, v2[0], 1, temp, 1);
    Vmath::Vvtvm(nq, v1[0], 1, v2[1], 1, temp, 1, v3[2], 1);
}

void MMFSystem::VectorCrossProd(const Array<OneD, NekDouble> &v1,
                                const Array<OneD, NekDouble> &v2,
                                Array<OneD, NekDouble> &v3)
{
    ASSERTL0(v1.size() == 3, "Input 1 has dimension not equal to 3.");
    ASSERTL0(v2.size() == 3, "Input 2 has dimension not equal to 3.");
    ASSERTL0(v3.size() == 3,
             "Output vector has dimension not equal to 3.");

    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void MMFSystem::ComputeCurl(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)

{
    int nq = inarray[0].size();

    Array<OneD, NekDouble> tmpx, tmpy, tmpz;
    Array<OneD, NekDouble> Dtmpzdx, Dtmpydx, Dtmpxdy, Dtmpzdy, Dtmpxdz, Dtmpydz;

    tmpx = Array<OneD, NekDouble>(nq);
    tmpy = Array<OneD, NekDouble>(nq);
    tmpz = Array<OneD, NekDouble>(nq);

    Dtmpzdx = Array<OneD, NekDouble>(nq);
    Dtmpydx = Array<OneD, NekDouble>(nq);
    Dtmpxdy = Array<OneD, NekDouble>(nq);
    Dtmpzdy = Array<OneD, NekDouble>(nq);
    Dtmpxdz = Array<OneD, NekDouble>(nq);
    Dtmpydz = Array<OneD, NekDouble>(nq);

    for (int k = 0; k < m_spacedim; ++k)
    {
        Vmath::Vcopy(nq, &inarray[0][0], 1, &tmpx[0], 1);
        Vmath::Vcopy(nq, &inarray[1][0], 1, &tmpy[0], 1);
        Vmath::Vcopy(nq, &inarray[2][0], 1, &tmpz[0], 1);

        m_fields[0]->PhysDeriv(0, tmpz, Dtmpzdx);
        m_fields[0]->PhysDeriv(0, tmpy, Dtmpydx);
        m_fields[0]->PhysDeriv(1, tmpx, Dtmpxdy);
        m_fields[0]->PhysDeriv(1, tmpz, Dtmpzdy);
        m_fields[0]->PhysDeriv(2, tmpx, Dtmpxdz);
        m_fields[0]->PhysDeriv(2, tmpy, Dtmpydz);

        Vmath::Vsub(nq, &Dtmpzdy[0], 1, &Dtmpydz[0], 1, &outarray[0][0], 1);
        Vmath::Vsub(nq, &Dtmpxdz[0], 1, &Dtmpzdx[0], 1, &outarray[1][0], 1);
        Vmath::Vsub(nq, &Dtmpydx[0], 1, &Dtmpxdy[0], 1, &outarray[2][0], 1);
    }
}

Array<OneD, NekDouble> MMFSystem::CartesianToMovingframes(
    const Array<OneD, const Array<OneD, NekDouble>> &uvec, unsigned int field)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> outarray(nq, 0.0);

    // new u0 = ( [u v] \cdot e^1 )/ |e^1|^2
    Vmath::Vmul(nq, &m_movingframes[field][0], 1, &uvec[0][0], 1, &outarray[0],
                1);
    Vmath::Vvtvp(nq, &m_movingframes[field][nq], 1, &uvec[1][0], 1,
                 &outarray[0], 1, &outarray[0], 1);
    Vmath::Vvtvp(nq, &m_movingframes[field][2 * nq], 1, &uvec[2][0], 1,
                 &outarray[0], 1, &outarray[0], 1);

    return outarray;
}

// x = r \cos \theta \cos \varphi
// y = r \cos \theta \sin \varphi
// z = r \sin \theta
void MMFSystem::CartesianToSpherical(const NekDouble x0j, const NekDouble x1j,
                                     const NekDouble x2j, NekDouble &sin_varphi,
                                     NekDouble &cos_varphi,
                                     NekDouble &sin_theta, NekDouble &cos_theta)
{
    NekDouble radius;
    NekDouble radxy;
    NekDouble Tol = 0.0000001;

    NekDouble m_Xscale = 1.0;
    NekDouble m_Yscale = 1.0;
    NekDouble m_Zscale = 1.0;

    radius = sqrt(x0j * x0j / (m_Xscale * m_Xscale) +
                  x1j * x1j / (m_Yscale * m_Yscale) +
                  x2j * x2j / (m_Zscale * m_Zscale));
    radxy = sqrt(x0j * x0j / (m_Xscale * m_Xscale) +
                 x1j * x1j / (m_Yscale * m_Yscale));

    if (radxy > Tol)
    {
        sin_varphi = x1j / (radxy * m_Yscale);
        cos_varphi = x0j / (radxy * m_Xscale);
    }

    else
    {
        sin_varphi = 0.0;
        if (x2j > 0)
        {
            cos_varphi = 1.0;
        }

        else
        {
            cos_varphi = -1.0;
        }
    }

    sin_theta = x2j / (radius * m_Zscale);
    cos_theta = radxy / radius;
}

void MMFSystem::CopyBoundaryTrace(
    const Array<OneD, const NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd,
    const BoundaryCopyType BDCopyType, const int var,
    const std::string BDtype)
{
    int id1, id2, npts, nptselem, cnt = 0, bdrycnt = 0;
    Array<OneD, NekDouble> Dirichlet, x0, x1, x2;

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[var]->GetBndConditions().size(); ++n)
    {
        nptselem = m_fields[var]->GetBndCondExpansions()[n]->GetNpoints();

        Dirichlet = Array<OneD, NekDouble>(nptselem);
        x0        = Array<OneD, NekDouble>(nptselem);
        x1        = Array<OneD, NekDouble>(nptselem);
        x2        = Array<OneD, NekDouble>(nptselem);

        if (BDCopyType == eDirichlet)
        {
            m_fields[var]->GetBndCondExpansions()[n]->GetCoords(x0, x1, x2);
            LibUtilities::EquationSharedPtr ifunc =
                m_session->GetFunction("BoundaryConditions", 0);
            ifunc->Evaluate(x0, x1, x2, 0.0, Dirichlet);
        }

        for (int e = 0;
             e < m_fields[var]->GetBndCondExpansions()[n]->GetExpSize(); ++e)
        {
            npts = m_fields[var]
                       ->GetBndCondExpansions()[n]
                       ->GetExp(e)
                       ->GetNumPoints(0);
            id1 = m_fields[var]->GetBndCondExpansions()[n]->GetPhys_Offset(e);
            id2 = m_fields[var]->GetTrace()->GetPhys_Offset(
                m_fields[var]->GetTraceMap()->GetBndCondIDToGlobalTraceID(
                    cnt + e));

            if (m_fields[var]->GetBndConditions()[n]->GetUserDefined() ==
                        BDtype || BDtype == "NoUserDefined")
            {
                switch (BDCopyType)
                {
                    case eDirichlet:
                    {
                        Vmath::Vcopy(npts, &Dirichlet[id1], 1, &Bwd[id2], 1);
                        bdrycnt++;
                    }
                    break;

                    case eFwdEQBwd:
                    {
                        Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                        bdrycnt++;
                    }
                    break;

                    case eFwdEQNegBwd:
                    {
                        Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                        Vmath::Neg(npts, &Bwd[id2], 1);
                        bdrycnt++;
                    }
                    break;

                    default:
                        break;
                }
            }
        }

        cnt += m_fields[var]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

// Compute  (e^[dir] \cdot ( n \times e^3)) * (imFwd * Fwd^3 + imBwd * Bwd^3 )
void MMFSystem::ComputeNtimesFz(const int dir,
                                const Array<OneD, Array<OneD, NekDouble>> &Fwd,
                                const Array<OneD, Array<OneD, NekDouble>> &Bwd,
                                const Array<OneD, const NekDouble> &imFwd,
                                const Array<OneD, const NekDouble> &imBwd,
                                Array<OneD, NekDouble> &outarrayFwd,
                                Array<OneD, NekDouble> &outarrayBwd)

{
    int nTraceNumPoints = GetTraceTotPoints();

    NekDouble tmpFwd, tmpBwd;
    NekDouble Aver, ntimesz;

    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        Aver = 0.5 * (imFwd[i] + imBwd[i]);

        tmpFwd = 0.0;
        tmpBwd = 0.0;
        for (int k = 0; k < m_spacedim; ++k)
        {
            ntimesz = 0.5 * (imFwd[i] * Fwd[2][i] + imBwd[i] * Bwd[2][i]);

            tmpFwd +=
                m_MFtraceFwd[dir][k][i] * m_ntimesMFFwd[2][k][i] * ntimesz;
            tmpBwd +=
                m_MFtraceBwd[dir][k][i] * m_ntimesMFBwd[2][k][i] * ntimesz;
        }

        outarrayFwd[i] = tmpFwd / Aver;
        outarrayBwd[i] = tmpBwd / Aver;
    }
}

// Compute e^3 \cdot ( n1e1 \times ( imFwd EFwd + imBwd EBwd ) ) / 2{{Y_i}}
void MMFSystem::ComputeNtimesF12(const Array<OneD, Array<OneD, NekDouble>> &Fwd,
                                 const Array<OneD, Array<OneD, NekDouble>> &Bwd,
                                 const Array<OneD, const NekDouble> &im1Fwd,
                                 const Array<OneD, const NekDouble> &im1Bwd,
                                 const Array<OneD, const NekDouble> &im2Fwd,
                                 const Array<OneD, const NekDouble> &im2Bwd,
                                 Array<OneD, NekDouble> &outarrayFwd,
                                 Array<OneD, NekDouble> &outarrayBwd)
{
    int nTraceNumPoints = GetTraceTotPoints();

    NekDouble tmpFwd, tmpBwd, Aver1, Aver2, HFwdk, HBwdk;

    Array<OneD, NekDouble> z1HAver(m_spacedim);
    Array<OneD, NekDouble> z2HAver(m_spacedim);
    Array<OneD, NekDouble> n1e1(m_spacedim);
    Array<OneD, NekDouble> n2e2(m_spacedim);

    Array<OneD, NekDouble> n1e1_times_z1HAver(m_spacedim);
    Array<OneD, NekDouble> n2e2_times_z2HAver(m_spacedim);

    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        Aver1 = 0.5 * (im1Fwd[i] + im1Bwd[i]);
        Aver2 = 0.5 * (im2Fwd[i] + im2Bwd[i]);

        for (int k = 0; k < m_spacedim; k++)
        {
            // Compute \vec{HFwd} and \vec{HBwd}
            HFwdk = Fwd[0][i] * m_MFtraceFwd[0][k][i] +
                    Fwd[1][i] * m_MFtraceFwd[1][k][i];
            HBwdk = Bwd[0][i] * m_MFtraceBwd[0][k][i] +
                    Bwd[1][i] * m_MFtraceBwd[1][k][i];

            // Compute z_i {{ \vec{H} }}
            z1HAver[k] = 0.5 * (im1Fwd[i] * HFwdk + im1Bwd[i] * HBwdk);
            z2HAver[k] = 0.5 * (im2Fwd[i] * HFwdk + im2Bwd[i] * HBwdk);

            // Choose e^i for the one in anisotropy region
            n1e1[k] = m_ncdotMFFwd[0][i] * m_MFtraceFwd[0][k][i];
            n2e2[k] = m_ncdotMFFwd[1][i] * m_MFtraceFwd[1][k][i];
        }

        // Compute n1e1 \times z1HAver and n2e2 \times z2HAver
        VectorCrossProd(n1e1, z1HAver, n1e1_times_z1HAver);
        VectorCrossProd(n2e2, z2HAver, n2e2_times_z2HAver);

        // e^3 \cdot ( n1e1 \times z1HAver + n2e2 \times z2HAver)
        tmpFwd = 0.0;
        tmpBwd = 0.0;
        for (int k = 0; k < m_spacedim; k++)
        {
            tmpFwd += m_MFtraceFwd[2][k][i] * (n1e1_times_z1HAver[k] / Aver1 +
                                               n2e2_times_z2HAver[k] / Aver2);
            tmpBwd += m_MFtraceBwd[2][k][i] * (n1e1_times_z1HAver[k] / Aver1 +
                                               n2e2_times_z2HAver[k] / Aver2);
        }

        outarrayFwd[i] = tmpFwd;
        outarrayBwd[i] = tmpBwd;
    }
}

void MMFSystem::ComputeNtimestimesdFz(
    const int dir, const Array<OneD, Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, Array<OneD, NekDouble>> &Bwd,
    const Array<OneD, const NekDouble> &imFwd,
    const Array<OneD, const NekDouble> &imBwd,
    Array<OneD, NekDouble> &outarrayFwd, Array<OneD, NekDouble> &outarrayBwd)
{
    int nTraceNumPoints = GetTraceTotPoints();

    Array<OneD, NekDouble> dH(m_spacedim);
    Array<OneD, NekDouble> nFwd(m_spacedim);

    Array<OneD, NekDouble> eiFwd(m_spacedim);
    Array<OneD, NekDouble> eiBwd(m_spacedim);

    Array<OneD, NekDouble> eitimesdHFwd(m_spacedim);
    Array<OneD, NekDouble> eitimesdHBwd(m_spacedim);

    Array<OneD, NekDouble> ntimeseitimesdHFwd(m_spacedim);
    Array<OneD, NekDouble> ntimeseitimesdHBwd(m_spacedim);

    NekDouble Aver, HFwdk, HBwdk, tmpFwd, tmpBwd;
    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        Aver = 0.5 * (imFwd[i] + imBwd[i]);

        // Get [H]
        for (int k = 0; k < m_spacedim; k++)
        {
            HFwdk = Fwd[0][i] * m_MFtraceFwd[0][k][i] +
                    Fwd[1][i] * m_MFtraceFwd[1][k][i];
            HBwdk = Bwd[0][i] * m_MFtraceBwd[0][k][i] +
                    Bwd[1][i] * m_MFtraceBwd[1][k][i];
            dH[k] = HFwdk - HBwdk;

            eiFwd[k] = m_MFtraceFwd[dir][k][i];
            eiBwd[k] = m_MFtraceBwd[dir][k][i];

            nFwd[k] = m_traceNormals[k][i];
        }

        // MFtraceFwd (MFtraceBwd) \times [H]
        // VectorCrossProd(eiFwd, dH, eitimesdHFwd);
        // VectorCrossProd(eiBwd, dH, eitimesdHBwd);
        VectorCrossProd(nFwd, dH, eitimesdHFwd);
        VectorCrossProd(nFwd, dH, eitimesdHBwd);

        // n times eitimesdH
        VectorCrossProd(nFwd, eitimesdHFwd, ntimeseitimesdHFwd);
        VectorCrossProd(nFwd, eitimesdHBwd, ntimeseitimesdHBwd);

        // MFtraceFwd \cdot ntimeseitimesdH
        tmpFwd = 0.0;
        tmpBwd = 0.0;
        for (int k = 0; k < m_spacedim; k++)
        {
            tmpFwd += eiFwd[k] * ntimeseitimesdHFwd[k];
            tmpBwd += eiBwd[k] * ntimeseitimesdHBwd[k];
        }

        outarrayFwd[i] = 0.5 * m_alpha * tmpFwd / Aver;
        outarrayBwd[i] = 0.5 * m_alpha * tmpBwd / Aver;
    }
}

// Compute - \alpha [E3] / ( {{im1}} + {{im2}} )
void MMFSystem::ComputeNtimestimesdF12(
    const Array<OneD, Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, Array<OneD, NekDouble>> &Bwd,
    const Array<OneD, const NekDouble> &im1Fwd,
    const Array<OneD, const NekDouble> &im1Bwd,
    const Array<OneD, const NekDouble> &im2Fwd,
    const Array<OneD, const NekDouble> &im2Bwd,
    Array<OneD, NekDouble> &outarrayFwd, Array<OneD, NekDouble> &outarrayBwd)
{
    int nTraceNumPoints = GetTraceTotPoints();

    Array<OneD, NekDouble> directFwd(nTraceNumPoints);
    Array<OneD, NekDouble> directBwd(nTraceNumPoints);

    NekDouble Aver1, Aver2;
    for (int i = 0; i < nTraceNumPoints; ++i)
    {
        Aver1 = im1Fwd[i] + im1Bwd[i];
        Aver2 = im2Fwd[i] + im2Bwd[i];

        outarrayFwd[i] =
            -m_alpha * (Fwd[2][i] - Bwd[2][i]) * (1.0 / Aver1 + 1.0 / Aver2);
        outarrayBwd[i] =
            -m_alpha * (Fwd[2][i] - Bwd[2][i]) * (1.0 / Aver1 + 1.0 / Aver2);
    }
}

// dim = dimension, pol = polarization, 0 = TM, 1 = TE.
void MMFSystem::ComputeZimYim(Array<OneD, Array<OneD, NekDouble>> &epsvec,
                              Array<OneD, Array<OneD, NekDouble>> &muvec)

{
    int nTraceNumPoints = GetTraceNpoints();

    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        case eScatField1D:
        {
            Array<OneD, NekDouble> Fwdeps(nTraceNumPoints, 1.0);
            Array<OneD, NekDouble> Bwdeps(nTraceNumPoints, 1.0);
            Array<OneD, NekDouble> Fwdmu(nTraceNumPoints, 1.0);
            Array<OneD, NekDouble> Bwdmu(nTraceNumPoints, 1.0);
            m_fields[0]->GetFwdBwdTracePhys(epsvec[0], Fwdeps, Bwdeps);
            m_fields[0]->GetFwdBwdTracePhys(muvec[0], Fwdeps, Bwdeps);

            CopyBoundaryTrace(Fwdeps, Bwdeps, eFwdEQBwd, 0);
            CopyBoundaryTrace(Fwdmu, Bwdmu, eFwdEQBwd, 1);

            m_ZimFwd = Array<OneD, Array<OneD, NekDouble>>(1);
            m_ZimBwd = Array<OneD, Array<OneD, NekDouble>>(1);
            m_YimFwd = Array<OneD, Array<OneD, NekDouble>>(1);
            m_YimBwd = Array<OneD, Array<OneD, NekDouble>>(1);

            m_ZimFwd[0] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
            m_ZimBwd[0] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
            m_YimFwd[0] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
            m_YimBwd[0] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);

            // ZimFwd = sqrt( muFwd / epsFwd),  ZimBwd = sqrt( muBwd / epsBwd)
            for (int i = 0; i < nTraceNumPoints; ++i)
            {
                m_ZimFwd[0][i] = sqrt(Fwdmu[i] / Fwdeps[i]);
                m_ZimBwd[0][i] = sqrt(Bwdmu[i] / Bwdeps[i]);

                m_YimFwd[0][i] = 1.0 / m_ZimFwd[0][i];
                m_YimBwd[0][i] = 1.0 / m_ZimBwd[0][i];
            }

            std::cout << "*** ZimFwd = " << RootMeanSquare(m_ZimFwd[0])
                      << ", ZimBwd = " << RootMeanSquare(m_ZimBwd[0])
                      << ", YimFwd = " << RootMeanSquare(m_YimFwd[0])
                      << ", YimBwd = " << RootMeanSquare(m_YimBwd[0])
                      << std::endl;
        }
        break;

        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        case eTestMaxwell2DPMC:
        case eScatField2D:
        case eTotField2D:
        case eMaxwellSphere:
        case eELF2DSurface:
        {
            m_ZimFwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
            m_ZimBwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
            m_YimFwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);
            m_YimBwd = Array<OneD, Array<OneD, NekDouble>>(m_shapedim);

            for (int j = 0; j < m_shapedim; ++j)
            {
                m_ZimFwd[j] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
                m_ZimBwd[j] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
                m_YimFwd[j] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
                m_YimBwd[j] = Array<OneD, NekDouble>(nTraceNumPoints, 1.0);
            }

            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    Array<OneD, NekDouble> Fwdmu1(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdmu1(nTraceNumPoints);
                    Array<OneD, NekDouble> Fwdmu2(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdmu2(nTraceNumPoints);
                    Array<OneD, NekDouble> Fwdeps3(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdeps3(nTraceNumPoints);

                    m_fields[0]->GetFwdBwdTracePhys(muvec[0], Fwdmu1, Bwdmu1);
                    m_fields[0]->GetFwdBwdTracePhys(muvec[1], Fwdmu2, Bwdmu2);
                    m_fields[0]->GetFwdBwdTracePhys(epsvec[2], Fwdeps3,
                                                    Bwdeps3);

                    CopyBoundaryTrace(Fwdmu1, Bwdmu1, eFwdEQBwd, 0);
                    CopyBoundaryTrace(Fwdmu2, Bwdmu2, eFwdEQBwd, 1);
                    CopyBoundaryTrace(Fwdeps3, Bwdeps3, eFwdEQBwd, 2);

                    // ZimFwd = sqrt( muFwd / epsFwd),  ZimBwd = sqrt( muBwd /
                    // epsBwd)
                    for (int i = 0; i < nTraceNumPoints; ++i)
                    {
                        m_ZimFwd[0][i] = sqrt(Fwdmu2[i] / Fwdeps3[i]);
                        m_ZimBwd[0][i] = sqrt(Bwdmu2[i] / Bwdeps3[i]);

                        m_YimFwd[0][i] = 1.0 / m_ZimFwd[0][i];
                        m_YimBwd[0][i] = 1.0 / m_ZimBwd[0][i];

                        m_ZimFwd[1][i] = sqrt(Fwdmu1[i] / Fwdeps3[i]);
                        m_ZimBwd[1][i] = sqrt(Bwdmu1[i] / Bwdeps3[i]);

                        m_YimFwd[1][i] = 1.0 / m_ZimFwd[1][i];
                        m_YimBwd[1][i] = 1.0 / m_ZimBwd[1][i];
                    }
                }
                break; // eTransMagnetic

                case eTransElectric:
                {
                    Array<OneD, NekDouble> Fwdeps1(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdeps1(nTraceNumPoints);
                    Array<OneD, NekDouble> Fwdeps2(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdeps2(nTraceNumPoints);
                    Array<OneD, NekDouble> Fwdmu3(nTraceNumPoints);
                    Array<OneD, NekDouble> Bwdmu3(nTraceNumPoints);

                    m_fields[0]->GetFwdBwdTracePhys(epsvec[0], Fwdeps1,
                                                    Bwdeps1);
                    m_fields[0]->GetFwdBwdTracePhys(epsvec[1], Fwdeps2,
                                                    Bwdeps2);
                    m_fields[0]->GetFwdBwdTracePhys(muvec[2], Fwdmu3, Bwdmu3);

                    CopyBoundaryTrace(Fwdeps1, Bwdeps1, eFwdEQBwd, 0);
                    CopyBoundaryTrace(Fwdeps2, Bwdeps2, eFwdEQBwd, 1);
                    CopyBoundaryTrace(Fwdmu3, Bwdmu3, eFwdEQBwd, 2);

                    for (int i = 0; i < nTraceNumPoints; ++i)
                    {
                        m_ZimFwd[0][i] = sqrt(Fwdmu3[i] / Fwdeps2[i]);
                        m_ZimBwd[0][i] = sqrt(Bwdmu3[i] / Bwdeps2[i]);

                        m_YimFwd[0][i] = 1.0 / m_ZimFwd[0][i];
                        m_YimBwd[0][i] = 1.0 / m_ZimBwd[0][i];

                        m_ZimFwd[1][i] = sqrt(Fwdmu3[i] / Fwdeps1[i]);
                        m_ZimBwd[1][i] = sqrt(Bwdmu3[i] / Bwdeps1[i]);

                        m_YimFwd[1][i] = 1.0 / m_ZimFwd[1][i];
                        m_YimBwd[1][i] = 1.0 / m_ZimBwd[1][i];
                    }
                }
                break; // eTransELectric

                default:
                    break;
            } // PolType

            std::cout << "*** ZimFwd0 = [ "
                      << Vmath::Vmin(nTraceNumPoints, m_ZimFwd[0], 1) << " , "
                      << Vmath::Vmax(nTraceNumPoints, m_ZimFwd[0], 1)
                      << " ], ZimBwd0 = [ "
                      << Vmath::Vmin(nTraceNumPoints, m_ZimBwd[0], 1) << " , "
                      << Vmath::Vmax(nTraceNumPoints, m_ZimBwd[0], 1) << " ] "
                      << std::endl;
            std::cout << "*** ZimFwd1 = [ "
                      << Vmath::Vmin(nTraceNumPoints, m_ZimFwd[1], 1) << " , "
                      << Vmath::Vmax(nTraceNumPoints, m_ZimFwd[1], 1)
                      << " ], ZimBwd1 = [ "
                      << Vmath::Vmin(nTraceNumPoints, m_ZimBwd[1], 1) << " , "
                      << Vmath::Vmax(nTraceNumPoints, m_ZimBwd[1], 1) << " ] "
                      << std::endl;
        }
        break; // eMaxwell2D

        default:
            break;
    } // TestMaxwellType
}

// m_dedxi_cdot_e[m][j][k][] = de^m / d \xi^j \cdot e^n
void MMFSystem::Computedemdxicdote()
{
    int MFdim = 3;
    int nq    = GetTotPoints();

    // Initialization
    m_dedxi_cdot_e =
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>(MFdim);
    for (int indm = 0; indm < MFdim; ++indm)
    {
        m_dedxi_cdot_e[indm] =
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(MFdim);
        for (int indj = 0; indj < MFdim; ++indj)
        {
            m_dedxi_cdot_e[indm][indj] =
                Array<OneD, Array<OneD, NekDouble>>(MFdim);
            for (int indn = 0; indn < MFdim; ++indn)
            {
                m_dedxi_cdot_e[indm][indj][indn] =
                    Array<OneD, NekDouble>(nq, 0.0);
            }
        }
    }

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpderiv(nq);
    Array<OneD, NekDouble> dedt;
    for (int indm = 0; indm < MFdim; ++indm)
    {
        for (int indj = 0; indj < MFdim; ++indj)
        {
            for (int indn = 0; indn < MFdim; ++indn)
            {
                dedt = Array<OneD, NekDouble>(nq, 0.0);
                for (int k = 0; k < m_spacedim; ++k)
                {
                    // Compute d e^m / d \xi_j cdot e^n
                    Vmath::Vcopy(nq, &m_movingframes[indm][k * nq], 1, &tmp[0],
                                 1);
                    m_fields[0]->PhysDirectionalDeriv(m_movingframes[indj], tmp,
                                                      tmpderiv);

                    Vmath::Vvtvp(nq, &tmpderiv[0], 1,
                                 &m_movingframes[indn][k * nq], 1, &dedt[0], 1,
                                 &dedt[0], 1);
                }

                Vmath::Vcopy(nq, &dedt[0], 1,
                             &m_dedxi_cdot_e[indm][indj][indn][0], 1);
            }
        }
    }

    int indx = 0;
    std::cout << "*** m_dedxi_cdot_e[0]/dxi1 = ( "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][0]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][1]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][2]) << " )_1, "
              << std::endl;
    std::cout << "*** m_dedxi_cdot_e[0]/dxi2 = ( "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][0]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][1]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][2]) << " )_2 "
              << std::endl;

    indx = 1;
    std::cout << "*** m_dedxi_cdot_e[1]/dxi1 = ( "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][0]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][1]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][2]) << " )_1, "
              << std::endl;
    std::cout << "*** m_dedxi_cdot_e[1]/dxi2 = ( "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][0]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][1]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][2]) << " )_2 "
              << std::endl;

    indx = 2;
    std::cout << "*** m_dedxi_cdot_e[2]/dxi1 = ( "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][0]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][1]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][0][2]) << " )_1, "
              << std::endl;
    std::cout << "*** m_dedxi_cdot_e[2]/dxi2 = ( "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][0]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][1]) << " , "
              << RootMeanSquare(m_dedxi_cdot_e[indx][1][2]) << " )_2 "
              << std::endl;
}

void MMFSystem::AdddedtMaxwell(
    const Array<OneD, const Array<OneD, NekDouble>> &physarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = GetTotPoints();

    // m_dedxi_cdot_e[m][j][n][] = de^m / d \xi^j \cdot e^n
    Array<OneD, NekDouble> dedtej(nq);
    Array<OneD, NekDouble> de1dtcdotej(nq);
    Array<OneD, NekDouble> de2dtcdotej(nq);

    Array<OneD, NekDouble> normH(nq);
    Array<OneD, NekDouble> NormedH1(nq, 0.0);
    Array<OneD, NekDouble> NormednegH2(nq, 0.0);

    Vmath::Vmul(nq, &physarray[0][0], 1, &physarray[0][0], 1, &normH[0], 1);
    Vmath::Vvtvp(nq, &physarray[1][0], 1, &physarray[1][0], 1, &normH[0], 1,
                 &normH[0], 1);
    Vmath::Vsqrt(nq, normH, 1, normH, 1);

    NekDouble Tol = 0.001;
    for (int i = 0; i < nq; ++i)
    {
        if (normH[i] > Tol)
        {
            NormedH1[i]    = physarray[0][i] / normH[i];
            NormednegH2[i] = -1.0 * physarray[1][i] / normH[i];
        }
    }

    for (int j = 0; j < m_shapedim; ++j)
    {
        // Compute de1 / dt \cdot ej = (-H2 de^1/d\xi1 \cdot e^j + H1 de^1/d\xi2
        // \cdot e^j) / sqrt{ H1^2 + H2^2 }
        Vmath::Vmul(nq, &NormednegH2[0], 1, &m_dedxi_cdot_e[0][0][j][0], 1,
                    &de1dtcdotej[0], 1);
        Vmath::Vvtvp(nq, &NormedH1[0], 1, &m_dedxi_cdot_e[0][1][j][0], 1,
                     &de1dtcdotej[0], 1, &de1dtcdotej[0], 1);

        // Compute de2 / dt \cdot ej = (-H2 de2/d\xi1 \cdot e^j + H1 de2/d\xi2
        // \cdot e^j) / sqrt{ H1^2 + H2^2 }
        Vmath::Vmul(nq, &NormednegH2[0], 1, &m_dedxi_cdot_e[1][0][j][0], 1,
                    &de2dtcdotej[0], 1);
        Vmath::Vvtvp(nq, &NormedH1[0], 1, &m_dedxi_cdot_e[1][1][j][0], 1,
                     &de2dtcdotej[0], 1, &de2dtcdotej[0], 1);

        // Add dedt component: (H1 (de1/dt) + H2 (de2/dt) ) \cdot ej
        Vmath::Vmul(nq, &physarray[0][0], 1, &de1dtcdotej[0], 1, &dedtej[0], 1);
        Vmath::Vvtvp(nq, &physarray[1][0], 1, &de2dtcdotej[0], 1, &dedtej[0], 1,
                     &dedtej[0], 1);

        Vmath::Neg(nq, dedtej, 1);

        switch (m_PolType)
        {
            case SolverUtils::eTransMagnetic:
            {
                if (j == 0)
                {
                    Vmath::Vmul(nq, m_muvec[0], 1, dedtej, 1, dedtej, 1);
                }

                else if (j == 1)
                {
                    Vmath::Vmul(nq, m_muvec[1], 1, dedtej, 1, dedtej, 1);
                }
            }
            break;

            case SolverUtils::eTransElectric:
            {
                if (j == 0)
                {
                    Vmath::Vmul(nq, m_epsvec[0], 1, dedtej, 1, dedtej, 1);
                }

                else if (j == 1)
                {
                    Vmath::Vmul(nq, m_epsvec[1], 1, dedtej, 1, dedtej, 1);
                }
            }
            break;

            default:
                break;
        }

        Vmath::Vadd(nq, &dedtej[0], 1, &outarray[j][0], 1, &outarray[j][0], 1);
    }
}

void MMFSystem::UpwindMaxwellFlux1D(
    Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd)
{
    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvar            = 2;

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvar);

    for (i = 0; i < nvar; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
    }

    // get the physical values at the trace from the dependent variables
    for (i = 0; i < nvar; ++i)
    {
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);
    }

    // E = 0 at the boundaries
    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
                      "PEC");

    // d H / d n = 0 at the boundaries
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
                      "PEC");

    Array<OneD, NekDouble> dE(nTraceNumPoints);
    Array<OneD, NekDouble> dH(nTraceNumPoints);

    Vmath::Vsub(nTraceNumPoints, &Fwd[0][0], 1, &Bwd[0][0], 1, &dE[0], 1);
    Vmath::Vsub(nTraceNumPoints, &Fwd[1][0], 1, &Bwd[1][0], 1, &dH[0], 1);

    NekDouble nx, AverZ, AverY, AverZH, AverYE;
    for (i = 0; i < nTraceNumPoints; ++i)
    {
        nx     = m_traceNormals[0][i];
        AverZ  = m_ZimFwd[0][i] + m_ZimBwd[0][i];
        AverY  = m_YimFwd[0][i] + m_YimBwd[0][i];
        AverZH = m_ZimFwd[0][i] * Fwd[1][i] + m_ZimBwd[0][i] * Bwd[1][i];
        AverYE = m_YimFwd[0][i] * Fwd[0][i] + m_YimBwd[0][i] * Bwd[0][i];

        numfluxFwd[0][i] = nx / AverZ * (AverZH - nx * dE[i]);
        numfluxFwd[1][i] = nx / AverY * (AverYE - nx * dH[i]);

        numfluxBwd[0][i] = nx / AverZ * (AverZH - nx * dE[i]);
        numfluxBwd[1][i] = nx / AverY * (AverYE - nx * dH[i]);
    }
}

void MMFSystem::LaxFriedrichMaxwellFlux1D(
    Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd)
{
    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvar            = 2;

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvar);

    for (i = 0; i < nvar; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
    }

    // get the physical values at the trace from the dependent variables
    for (i = 0; i < nvar; ++i)
    {
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);
    }

    // E = 0 at the boundaries
    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
                      "PEC");

    // d H / d n = 0 at the boundaries
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
                      "PEC");

    Array<OneD, NekDouble> dE(nTraceNumPoints);
    Array<OneD, NekDouble> dH(nTraceNumPoints);

    Vmath::Vsub(nTraceNumPoints, &Fwd[0][0], 1, &Bwd[0][0], 1, &dE[0], 1);
    Vmath::Vsub(nTraceNumPoints, &Fwd[1][0], 1, &Bwd[1][0], 1, &dH[0], 1);

    NekDouble nx;
    for (i = 0; i < nTraceNumPoints; ++i)
    {
        nx = m_traceNormals[0][i];
        numfluxFwd[0][i] =
            0.5 * nx * ((Fwd[1][i] + Bwd[1][i]) - nx * m_YimFwd[0][i] * dE[i]);
        numfluxFwd[1][i] =
            0.5 * nx * ((Fwd[0][i] + Bwd[0][i]) - nx * m_ZimFwd[0][i] * dH[i]);

        numfluxBwd[0][i] =
            0.5 * nx * ((Fwd[1][i] + Bwd[1][i]) - nx * m_YimFwd[0][i] * dE[i]);
        numfluxBwd[1][i] =
            0.5 * nx * ((Fwd[0][i] + Bwd[0][i]) - nx * m_ZimFwd[0][i] * dH[i]);
    }
}

void MMFSystem::AverageMaxwellFlux1D(
    Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd)
{
    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvar            = 2;

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvar);

    for (i = 0; i < nvar; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
    }

    // get the physical values at the trace from the dependent variables
    for (i = 0; i < nvar; ++i)
    {
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);
    }

    // E = 0 at the boundaries
    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
                      "PEC");

    // d H / d n = 0 at the boundaries
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
                      "PEC");

    for (i = 0; i < nTraceNumPoints; ++i)
    {
        numfluxFwd[0][i] = 0.5 * m_traceNormals[0][i] * (Fwd[1][i] + Bwd[1][i]);
        numfluxFwd[1][i] = 0.5 * m_traceNormals[0][i] * (Fwd[0][i] + Bwd[0][i]);

        numfluxBwd[0][i] = 0.5 * m_traceNormals[0][i] * (Fwd[1][i] + Bwd[1][i]);
        numfluxBwd[1][i] = 0.5 * m_traceNormals[0][i] * (Fwd[0][i] + Bwd[0][i]);
    }
}

void MMFSystem::GetMaxwellFluxVector(
    const int var, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        case eScatField1D:
        {
            GetMaxwellFlux1D(var, physfield, flux);
        }
        break;

        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        case eTestMaxwell2DPMC:
        case eScatField2D:
        case eTotField2D:
        case eMaxwellSphere:
        case eELF2DSurface:
        {
            GetMaxwellFlux2D(var, physfield, flux);
        }
        break;

        default:
            break;
    }
}

void MMFSystem::GetMaxwellFlux1D(
    const int var, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    int nq = m_fields[0]->GetTotPoints();

    switch (var)
    {
        case 0:
        {
            // H in flux 0
            Vmath::Vcopy(nq, physfield[1], 1, flux[0], 1);

            // E in flux 1
            Vmath::Zero(nq, flux[1], 1);
        }
        break;

        case 1:
        {
            // E in flux 0
            Vmath::Vcopy(nq, physfield[0], 1, flux[0], 1);

            // H in flux 1
            Vmath::Zero(nq, flux[1], 1);
        }
        break;
        //----------------------------------------------------

        default:
            break;
    }
}

void MMFSystem::GetMaxwellFlux2D(
    const int var, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    int nq = m_fields[0]->GetTotPoints();

    NekDouble sign = 1.0;
    switch (m_PolType)
    {
        // TransMagnetic
        case 0:
        {
            sign = -1.0;
        }
        break;

        // TransElectric
        case 1:
        {
            sign = 1.0;
        }
        break;

        default:
            break;
    }

    switch (var)
    {
        case 0:
        {
            // -Ez in flux 1
            Vmath::Smul(nq, sign, physfield[2], 1, flux[0], 1);
            Vmath::Zero(nq, flux[1], 1);
        }
        break;

        case 1:
        {
            // Ez in flux 0
            Vmath::Zero(nq, flux[0], 1);
            Vmath::Smul(nq, -sign, physfield[2], 1, flux[1], 1);
        }
        break;

        case 2:
        {
            Vmath::Smul(nq, sign, physfield[0], 1, flux[0], 1);
            Vmath::Smul(nq, -sign, physfield[1], 1, flux[1], 1);
        }
        break;

        default:
            ASSERTL0(false, "GetFluxVector2D: illegal vector index");
    }
}

void MMFSystem::NumericalMaxwellFlux(
    Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd, const NekDouble time)
{

    switch (m_TestMaxwellType)
    {
        case eMaxwell1D:
        case eScatField1D:
        {
            switch (m_upwindType)
            {
                case eAverage:
                {
                    AverageMaxwellFlux1D(physfield, numfluxFwd, numfluxBwd);
                }
                break;

                case eLaxFriedrich:
                {
                    LaxFriedrichMaxwellFlux1D(physfield, numfluxFwd,
                                              numfluxBwd);
                }
                break;

                case eUpwind:
                {
                    UpwindMaxwellFlux1D(physfield, numfluxFwd, numfluxBwd);
                }
                break;

                default:
                {
                    ASSERTL0(false,
                             "populate switch statement for upwind flux");
                }
                break; // upwindType
            }
        }
        break; // eMaxwell1D

        case eTestMaxwell2DPEC:
        case eTestMaxwell2DPECAVGFLUX:
        case eTestMaxwell2DPMC:
        case eScatField2D:
        case eTotField2D:
        case eMaxwellSphere:
        case eELF2DSurface:
        {
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    NumericalMaxwellFluxTM(physfield, numfluxFwd, numfluxBwd,
                                           time);
                }
                break;

                case eTransElectric:
                {
                    NumericalMaxwellFluxTE(physfield, numfluxFwd, numfluxBwd,
                                           time);
                }
                break;

                default:
                    break;
            }
        }
        break; // eMaxwell2D

        default:
            break;
    } // m_TestMaxwellType
}

void MMFSystem::NumericalMaxwellFluxTM(
    Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd, const NekDouble time)
{
    int nq              = m_fields[0]->GetNpoints();
    int nTraceNumPoints = GetTraceTotPoints();
    int nvar            = physfield.size();

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        // get the physical values at the trace
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);
    }

    // E^|| = 0 at the PEC boundaries vs. H^|| = 0 at PMC boundaries
    Array<OneD, NekDouble> IncField(nq, 0.0);
    Array<OneD, NekDouble> IncFieldBwd(nTraceNumPoints, 0.0);
    Array<OneD, Array<OneD, NekDouble>> IncFieldFwd(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        IncFieldFwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        IncField = GetIncidentField(i, time);
        m_fields[0]->GetFwdBwdTracePhys(IncField, IncFieldFwd[i], IncFieldBwd);

        Vmath::Svtvp(nTraceNumPoints, 2.0, &IncFieldFwd[i][0], 1, &Fwd[i][0], 1,
                     &IncFieldFwd[i][0], 1);
    }

    // Total Field Formulation
    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 0,
                      "PEC_Forces");
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
                      "PEC_Forces");
    CopyBoundaryTrace(IncFieldFwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
                      "PEC_Forces");

    CopyBoundaryTrace(IncFieldFwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
                      "PMC_Forces");
    CopyBoundaryTrace(IncFieldFwd[1], Bwd[1], SolverUtils::eFwdEQNegBwd, 1,
                      "PMC_Forces");
    CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd, 2,
                      "PMC_Forces");

    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 0,
                      "PEC");
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
                      "PEC");
    CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
                      "PEC");

    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
                      "PMC");
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQNegBwd, 1,
                      "PMC");
    CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd, 2,
                      "PMC");

    Array<OneD, NekDouble> e1Fwd_cdot_ncrossdH(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e1Bwd_cdot_ncrossdH(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e2Fwd_cdot_ncrossdH(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e2Bwd_cdot_ncrossdH(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e3Fwd_cdot_dEe3(nTraceNumPoints, 0.0);
    Array<OneD, NekDouble> e3Bwd_cdot_dEe3(nTraceNumPoints, 0.0);

    // Compute  numfluxFwd[dir] = (eFwd^[dir] \cdot n \times e^3) * (YimFwd *
    // EFwd^3 + YimBwd * EBwd^3 )
    ComputeNtimesFz(0, Fwd, Bwd, m_YimFwd[0], m_YimBwd[0], numfluxFwd[0],
                    numfluxBwd[0]);
    ComputeNtimesFz(1, Fwd, Bwd, m_YimFwd[1], m_YimBwd[1], numfluxFwd[1],
                    numfluxBwd[1]);

    // Compute numfluxFwd[2] = eFwd^3 \cdot ( n1e1 \times ( ZimFwd HFwd + ZimBwd
    // HBwd ) ) / 2 {{Z_i}}
    ComputeNtimesF12(Fwd, Bwd, m_ZimFwd[0], m_ZimBwd[0], m_ZimFwd[1],
                     m_ZimBwd[1], numfluxFwd[2], numfluxBwd[2]);

    // Compute e1Fwd_cdot_ncrossdE = eFwd[dir] \cdot \alpha n \times n \times
    // [H] / 2 {{YimFwd}}
    ComputeNtimestimesdFz(0, Fwd, Bwd, m_YimFwd[0], m_YimBwd[0],
                          e1Fwd_cdot_ncrossdH, e1Bwd_cdot_ncrossdH);
    ComputeNtimestimesdFz(1, Fwd, Bwd, m_YimFwd[1], m_YimBwd[1],
                          e2Fwd_cdot_ncrossdH, e2Bwd_cdot_ncrossdH);

    // Compute  \alpha [E3] * ( 1/2{{Zim1}} + 1/2{{Zim2}} )
    ComputeNtimestimesdF12(Fwd, Bwd, m_ZimFwd[0], m_ZimBwd[0], m_ZimFwd[1],
                           m_ZimBwd[1], e3Fwd_cdot_dEe3, e3Bwd_cdot_dEe3);

    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxFwd[0], 1, e1Fwd_cdot_ncrossdH,
                 1, numfluxFwd[0], 1);
    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxFwd[1], 1, e2Fwd_cdot_ncrossdH,
                 1, numfluxFwd[1], 1);
    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxFwd[2], 1, e3Fwd_cdot_dEe3, 1,
                 numfluxFwd[2], 1);

    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxBwd[0], 1, e1Bwd_cdot_ncrossdH,
                 1, numfluxBwd[0], 1);
    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxBwd[1], 1, e2Bwd_cdot_ncrossdH,
                 1, numfluxBwd[1], 1);
    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxBwd[2], 1, e3Bwd_cdot_dEe3, 1,
                 numfluxBwd[2], 1);
}

void MMFSystem::NumericalMaxwellFluxTE(
    Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
    Array<OneD, Array<OneD, NekDouble>> &numfluxBwd, const NekDouble time)

{
    int nq              = m_fields[0]->GetNpoints();
    int nTraceNumPoints = GetTraceTotPoints();
    int nvar            = physfield.size();

    // Get temporary arrays
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        // get the physical values at the trace
        m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd[i], Bwd[i]);
    }

    // E = 0 at the PEC boundaries:
    Array<OneD, NekDouble> IncField(nq, 0.0);
    Array<OneD, NekDouble> IncFieldBwd(nTraceNumPoints, 0.0);
    Array<OneD, Array<OneD, NekDouble>> IncFieldFwd(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        IncFieldFwd[i] = Array<OneD, NekDouble>(nTraceNumPoints, 0.0);

        IncField = GetIncidentField(i, time);
        m_fields[0]->GetFwdBwdTracePhys(IncField, IncFieldFwd[i], IncFieldBwd);

        Vmath::Svtvp(nTraceNumPoints, 2.0, &IncFieldFwd[i][0], 1, &Fwd[i][0], 1,
                     &IncFieldFwd[i][0], 1);
    }

    CopyBoundaryTrace(IncFieldFwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
                      "PEC_Forces");
    CopyBoundaryTrace(IncFieldFwd[1], Bwd[1], SolverUtils::eFwdEQNegBwd, 1,
                      "PEC_Forces");
    CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd, 2,
                      "PEC_Forces");

    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 0,
                      "PMC_Forces");
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
                      "PMC_Forces");
    CopyBoundaryTrace(IncFieldFwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
                      "PMC_Forces");

    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQNegBwd, 0,
                      "PEC");
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQNegBwd, 1,
                      "PEC");
    CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQBwd, 2,
                      "PEC");

    CopyBoundaryTrace(Fwd[0], Bwd[0], SolverUtils::eFwdEQBwd, 0,
                      "PMC");
    CopyBoundaryTrace(Fwd[1], Bwd[1], SolverUtils::eFwdEQBwd, 1,
                      "PMC");
    CopyBoundaryTrace(Fwd[2], Bwd[2], SolverUtils::eFwdEQNegBwd, 2,
                      "PMC");

    Array<OneD, NekDouble> e1Fwd_cdot_ncrossdE(nTraceNumPoints);
    Array<OneD, NekDouble> e1Bwd_cdot_ncrossdE(nTraceNumPoints);
    Array<OneD, NekDouble> e2Fwd_cdot_ncrossdE(nTraceNumPoints);
    Array<OneD, NekDouble> e2Bwd_cdot_ncrossdE(nTraceNumPoints);
    Array<OneD, NekDouble> e3Fwd_cdot_dHe3(nTraceNumPoints);
    Array<OneD, NekDouble> e3Bwd_cdot_dHe3(nTraceNumPoints);

    // Compute  numfluxFwd[dir] = (eFwd^[dir] \cdot n \times e^3) * (ZimFwd *
    // HFwd^3 + ZimBwd * HBwd^3 )
    // Compute  numfluxBwd[dir] = (eBwd^[dir] \cdot n \times e^3) * (ZimFwd *
    // HFwd^3 + ZimBwd * HBwd^3 )
    ComputeNtimesFz(0, Fwd, Bwd, m_ZimFwd[0], m_ZimBwd[0], numfluxFwd[0],
                    numfluxBwd[0]);
    ComputeNtimesFz(1, Fwd, Bwd, m_ZimFwd[1], m_ZimBwd[1], numfluxFwd[1],
                    numfluxBwd[1]);

    // Compute numfluxFwd[2] = eFwd^3 \cdot ( n1e1 \times ( imFwd EFwd + imBwd
    // EBwd ) ) / 2 {{Y_i}}
    ComputeNtimesF12(Fwd, Bwd, m_YimFwd[0], m_YimBwd[0], m_YimFwd[1],
                     m_YimBwd[1], numfluxFwd[2], numfluxBwd[2]);

    // Compute e1Fwd_cdot_ncrossdE = eFwd[dir] \cdot \alpha n \times n \times
    // [E] / 2 {{ZimFwd}}
    ComputeNtimestimesdFz(0, Fwd, Bwd, m_ZimFwd[0], m_ZimBwd[0],
                          e1Fwd_cdot_ncrossdE, e1Bwd_cdot_ncrossdE);
    ComputeNtimestimesdFz(1, Fwd, Bwd, m_ZimFwd[1], m_ZimBwd[1],
                          e2Fwd_cdot_ncrossdE, e2Bwd_cdot_ncrossdE);

    // Compute  - \alpha [H3] * ( 1/2{{Yim1}} + 1/2{{Yim2}} )
    ComputeNtimestimesdF12(Fwd, Bwd, m_YimFwd[0], m_YimBwd[0], m_YimFwd[1],
                           m_YimBwd[1], e3Fwd_cdot_dHe3, e3Bwd_cdot_dHe3);

    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxFwd[0], 1, e1Fwd_cdot_ncrossdE, 1,
                 numfluxFwd[0], 1);
    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxFwd[1], 1, e2Fwd_cdot_ncrossdE, 1,
                 numfluxFwd[1], 1);
    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxFwd[2], 1, e3Fwd_cdot_dHe3, 1,
                 numfluxFwd[2], 1);

    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxBwd[0], 1, e1Bwd_cdot_ncrossdE, 1,
                 numfluxBwd[0], 1);
    Vmath::Svtvp(nTraceNumPoints, 1.0, numfluxBwd[1], 1, e2Bwd_cdot_ncrossdE, 1,
                 numfluxBwd[1], 1);
    Vmath::Svtvp(nTraceNumPoints, -1.0, numfluxBwd[2], 1, e3Bwd_cdot_dHe3, 1,
                 numfluxBwd[2], 1);
}

Array<OneD, NekDouble> MMFSystem::GetIncidentField(const int var,
                                                   const NekDouble time)
{
    int nq = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0, x1, x2);

    // GetSmoothFactor such that wave propages from the left to the object.
    // a = 0.1, ta = 1, f = 1.0./(1.0 + exp( -0.5.*(time-ta)/a ));
    Array<OneD, NekDouble> SmoothFactor(nq, 1.0);
    switch (m_SmoothFactor)
    {
        case 0:
        {
            for (int i = 0; i < nq; i++)
            {
                SmoothFactor[i] = 1.0 / (1.0 + exp(-1.0 * (time - 1.0) / 0.1));
            }
        }
        break;

        case 1:
        {
            NekDouble xmin = Vmath::Vmin(nq, x0, 1);
            NekDouble xp;
            for (int i = 0; i < nq; i++)
            {
                xp = x0[i] - xmin - time - m_SFinit;
                if (xp > 0.0)
                {
                    SmoothFactor[i] =
                        2.0 / (1.0 + exp(0.5 * (sqrt(xp * xp) - 0.1)));
                }

                else
                {
                    SmoothFactor[i] = 1.0;
                }
            }
        }
        break;

        default:
            break;
    }

    // Generate a factor for smoothly increasing wave
    Array<OneD, NekDouble> F1(nq);
    Array<OneD, NekDouble> F2(nq);
    Array<OneD, NekDouble> F3(nq);

    Array<OneD, NekDouble> dF1dt(nq);
    Array<OneD, NekDouble> dF2dt(nq);
    Array<OneD, NekDouble> dF3dt(nq);

    Array<OneD, NekDouble> F1int(nq);
    Array<OneD, NekDouble> F2int(nq);
    Array<OneD, NekDouble> F3int(nq);

    Array<OneD, NekDouble> outarray(nq);

    switch (m_IncType)
    {
        case ePlaneWave:
        {
            NekDouble cs, sn;
            NekDouble e1y, e2y;
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    // H = 0 \hat{x} - e^{ikx} \hat{y}
                    // E = e^{ikx} \hat{z}
                    // outarray1 = Hr1inc
                    // outarray2 = Hr2inc
                    // outarray3 = Ezrinc
                    for (int i = 0; i < nq; i++)
                    {
                        e1y = m_movingframes[0][nq + i];
                        e2y = m_movingframes[1][nq + i];

                        cs = SmoothFactor[i] * cos(m_Incfreq * (x0[i] - time));
                        sn = SmoothFactor[i] * sin(m_Incfreq * (x0[i] - time));

                        F1[i] = -1.0 * cs * e1y;
                        F2[i] = -1.0 * cs * e2y;
                        F3[i] = cs;

                        dF1dt[i] = -1.0 * m_Incfreq * sn * e1y;
                        dF2dt[i] = -1.0 * m_Incfreq * sn * e2y;
                        dF3dt[i] = 1.0 * m_Incfreq * sn;

                        F1int[i] = (1.0 / m_Incfreq) * sn * e1y;
                        F2int[i] = (1.0 / m_Incfreq) * sn * e2y;
                        F3int[i] = (-1.0 / m_Incfreq) * sn;
                    }
                }
                break;

                case eTransElectric:
                {
                    // E = 0 \hat{x} + e^{ikx} \hat{y}
                    // H = e^{ikx} \hat{z}
                    // outarray1 = Er1inc
                    // outarray2 = Er2inc
                    // outarray3 = Hzrinc
                    // outarray1 = Ei1inc
                    // outarray2 = Ei2inc
                    // outarray3 = Hziinc
                    for (int i = 0; i < nq; i++)
                    {
                        e1y = m_movingframes[0][nq + i];
                        e2y = m_movingframes[1][nq + i];

                        cs = SmoothFactor[i] * cos(m_Incfreq * (x0[i] - time));
                        sn = SmoothFactor[i] * sin(m_Incfreq * (x0[i] - time));

                        F1[i] = cs * e1y;
                        F2[i] = cs * e2y;
                        F3[i] = cs;

                        dF1dt[i] = m_Incfreq * sn * e1y;
                        dF2dt[i] = m_Incfreq * sn * e2y;
                        dF3dt[i] = m_Incfreq * sn;

                        F1int[i] = (-1.0 / m_Incfreq) * sn * e1y;
                        F2int[i] = (-1.0 / m_Incfreq) * sn * e2y;
                        F3int[i] = (-1.0 / m_Incfreq) * sn;
                    }
                }
                break;

                default:
                    break;
            }
        }
        break;

        case ePlaneWaveImag:
        {
            NekDouble cs, sn;
            NekDouble e1y, e2y;
            switch (m_PolType)
            {
                case eTransMagnetic:
                {
                    // H = 0 \hat{x} - e^{ikx} \hat{y}
                    // E = e^{ikx} \hat{z}
                    // outarray1 = Hr1inc
                    // outarray2 = Hr2inc
                    // outarray3 = Ezrinc
                    for (int i = 0; i < nq; i++)
                    {
                        e1y = m_movingframes[0][nq + i];
                        e2y = m_movingframes[1][nq + i];

                        cs = SmoothFactor[i] * cos(m_Incfreq * (x0[i] - time));
                        sn = SmoothFactor[i] * sin(m_Incfreq * (x0[i] - time));

                        F1[i] = -1.0 * sn * e1y;
                        F2[i] = -1.0 * sn * e2y;
                        F3[i] = sn;

                        dF1dt[i] = m_Incfreq * cs * e1y;
                        dF2dt[i] = m_Incfreq * cs * e2y;
                        dF3dt[i] = -1.0 * m_Incfreq * cs;

                        F1int[i] = (-1.0 / m_Incfreq) * cs * e1y;
                        F2int[i] = (-1.0 / m_Incfreq) * cs * e2y;
                        F3int[i] = (1.0 / m_Incfreq) * cs;
                    }
                }
                break;

                case eTransElectric:
                {
                    // E = 0 \hat{x} + e^{ikx} \hat{y}
                    // H = e^{ikx} \hat{z}
                    // outarray1 = Er1inc
                    // outarray2 = Er2inc
                    // outarray3 = Hzrinc
                    // outarray1 = Ei1inc
                    // outarray2 = Ei2inc
                    // outarray3 = Hziinc
                    for (int i = 0; i < nq; i++)
                    {
                        e1y = m_movingframes[0][nq + i];
                        e2y = m_movingframes[1][nq + i];

                        cs = SmoothFactor[i] * cos(m_Incfreq * (x0[i] - time));
                        sn = SmoothFactor[i] * sin(m_Incfreq * (x0[i] - time));

                        F1[i] = sn * e1y;
                        F2[i] = sn * e2y;
                        F3[i] = sn;

                        dF1dt[i] = -1.0 * m_Incfreq * cs * e1y;
                        dF2dt[i] = -1.0 * m_Incfreq * cs * e2y;
                        dF3dt[i] = -1.0 * m_Incfreq * cs;

                        F1int[i] = (1.0 / m_Incfreq) * cs * e1y;
                        F2int[i] = (1.0 / m_Incfreq) * cs * e2y;
                        F3int[i] = (1.0 / m_Incfreq) * cs;
                    }
                }
                break;

                default:
                    break;
            }
        }
        break;

        default:
            break;
    }

    switch (var)
    {
        case 0:
        {
            outarray = F1;
        }
        break;

        case 1:
        {
            outarray = F2;
        }
        break;

        case 2:
        {
            outarray = F3;
        }
        break;

        case 10:
        {
            outarray = dF1dt;
        }
        break;

        case 11:
        {
            outarray = dF2dt;
        }
        break;

        case 12:
        {
            outarray = dF3dt;
        }
        break;

        case 20:
        {
            outarray = F1int;
        }
        break;

        case 21:
        {
            outarray = F2int;
        }
        break;

        case 22:
        {
            outarray = F3int;
        }
        break;

        default:
        {
            Vmath::Zero(nq, outarray, 1);
        }
        break;
    }

    return outarray;
}

NekDouble MMFSystem::AvgInt(const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();
    Array<OneD, NekDouble> Ones(nq, 1.0);

    if (inarray.size() != nq)
    {
        ASSERTL0(false, "AvgInt Error: Vector size is not correct");
    }

    NekDouble jac = m_fields[0]->PhysIntegral(Ones);

    return (m_fields[0]->PhysIntegral(inarray)) / jac;
}

NekDouble MMFSystem::AvgAbsInt(const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();
    Array<OneD, NekDouble> Ones(nq, 1.0);
    Array<OneD, NekDouble> tmp(nq);

    if (inarray.size() != nq)
    {
        ASSERTL0(false, "AvgAbsInt Error: Vector size is not correct");
    }

    NekDouble jac = m_fields[0]->PhysIntegral(Ones);

    Vmath::Vabs(nq, inarray, 1, tmp, 1);
    return (m_fields[0]->PhysIntegral(tmp)) / jac;
}

NekDouble MMFSystem::AbsIntegral(const Array<OneD, const NekDouble> &inarray)
{
    int nq = m_fields[0]->GetNpoints();
    Array<OneD, NekDouble> tmp(nq);

    if (inarray.size() != nq)
    {
        ASSERTL0(false, "AbsIntegral Error: Vector size is not correct");
    }

    Vmath::Vabs(nq, inarray, 1, tmp, 1);
    return m_fields[0]->PhysIntegral(tmp);
}

NekDouble MMFSystem::RootMeanSquare(const Array<OneD, const NekDouble> &inarray)
{
    int Ntot = inarray.size();

    NekDouble reval = 0.0;
    for (int i = 0; i < Ntot; ++i)
    {
        reval = reval + inarray[i] * inarray[i];
    }
    reval = sqrt(reval / Ntot);

    return reval;
}

NekDouble MMFSystem::VectorAvgMagnitude(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray)
{
    int nq = inarray[0].size();

    Array<OneD, NekDouble> tmp(nq, 0.0);
    for (int k = 0; k < m_spacedim; k++)
    {
        Vmath::Vvtvp(nq, &inarray[k][0], 1, &inarray[k][0], 1, &tmp[0], 1,
                     &tmp[0], 1);
    }
    Vmath::Vsqrt(nq, tmp, 1, tmp, 1);

    return RootMeanSquare(tmp);
}

void MMFSystem::BubbleSort(Array<OneD, NekDouble> &refarray,
                           Array<OneD, NekDouble> &sortarray)
{
    int nq = refarray.size();

    bool swapped = true;
    int j        = 0;
    NekDouble tmp;

    while (swapped)
    {
        swapped = false;
        j++;
        for (int i = 0; i < nq - j; i++)
        {
            if (refarray[i] > refarray[i + 1])
            {
                tmp             = refarray[i];
                refarray[i]     = refarray[i + 1];
                refarray[i + 1] = tmp;

                tmp              = sortarray[i];
                sortarray[i]     = sortarray[i + 1];
                sortarray[i + 1] = tmp;

                swapped = true;
            }
        }
    }
}

void MMFSystem::GramSchumitz(
    const Array<OneD, const Array<OneD, NekDouble>> &v1,
    const Array<OneD, const Array<OneD, NekDouble>> &v2,
    Array<OneD, Array<OneD, NekDouble>> &outarray, bool KeepTheMagnitude)
{

    int nq = v1[0].size();
    Array<OneD, NekDouble> tmp(nq, 0.0);
    Array<OneD, NekDouble> mag(nq, 0.0);

    for (int i = 0; i < m_spacedim; ++i)
    {
        // u2 = v2 - < u1 , v2 > ( u1 / < u1, u1 > )
        Vmath::Vvtvp(nq, &v1[i][0], 1, &v2[i][0], 1, &tmp[0], 1, &tmp[0], 1);
        Vmath::Vvtvp(nq, &v1[i][0], 1, &v1[i][0], 1, &mag[0], 1, &mag[0], 1);
    }
    Vmath::Vdiv(nq, &tmp[0], 1, &mag[0], 1, &tmp[0], 1);
    Vmath::Neg(nq, &tmp[0], 1);

    // outarray = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);

    // u2 = v2 - < u1 , v2 > ( u1 / < u1, u1 > )
    for (int i = 0; i < m_spacedim; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(nq, 0.0);
        Vmath::Vvtvp(nq, &tmp[0], 1, &v1[i][0], 1, &v2[i][0], 1,
                     &outarray[i][0], 1);
    }

    if (KeepTheMagnitude)
    {
        Array<OneD, NekDouble> magorig(nq, 0.0);
        Array<OneD, NekDouble> magnew(nq, 0.0);

        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vmul(nq, &v2[0][0], 1, &v2[0][0], 1, &magorig[0], 1);
            Vmath::Vvtvp(nq, &v2[1][0], 1, &v2[1][0], 1, &magorig[0], 1,
                         &magorig[0], 1);
            Vmath::Vvtvp(nq, &v2[2][0], 1, &v2[2][0], 1, &magorig[0], 1,
                         &magorig[0], 1);

            Vmath::Vmul(nq, &outarray[0][0], 1, &outarray[0][0], 1, &magnew[0],
                        1);
            Vmath::Vvtvp(nq, &outarray[1][0], 1, &outarray[1][0], 1, &magnew[0],
                         1, &magnew[0], 1);
            Vmath::Vvtvp(nq, &outarray[2][0], 1, &outarray[2][0], 1, &magnew[0],
                         1, &magnew[0], 1);
        }

        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = 0; j < nq; ++j)
            {
                if (fabs(magnew[j]) > 0.000000001)
                {
                    outarray[i][j] =
                        outarray[i][j] * sqrt(magorig[j] / magnew[j]);
                }
            }
        }
    }
}

void MMFSystem::v_GenerateSummary(SummaryList &s)
{
    int nq = m_fields[0]->GetNpoints();
    UnsteadySystem::v_GenerateSummary(s);

    AddSummaryItem(s, "Total grids", nq);
    AddSummaryItem(s, "Shape Dimension", m_shapedim);
    AddSummaryItem(s, "Surface", SurfaceTypeMap[m_surfaceType]);
    if (m_surfaceType == eSphere)
    {
        NekDouble MeshError;

        Array<OneD, NekDouble> x(nq);
        Array<OneD, NekDouble> y(nq);
        Array<OneD, NekDouble> z(nq);
        Array<OneD, NekDouble> rad(nq, 0.0);

        m_fields[0]->GetCoords(x, y, z);

        Vmath::Vvtvp(nq, x, 1, x, 1, rad, 1, rad, 1);
        Vmath::Vvtvp(nq, y, 1, y, 1, rad, 1, rad, 1);
        Vmath::Vvtvp(nq, z, 1, z, 1, rad, 1, rad, 1);
        Vmath::Vsqrt(nq, rad, 1, rad, 1);

        Vmath::Sadd(nq, -1.0, rad, 1, rad, 1);
        Vmath::Vabs(nq, rad, 1, rad, 1);

        MeshError = m_fields[0]->PhysIntegral(rad);
        SolverUtils::AddSummaryItem(s, "Mesh Error", MeshError);
    }

    AddSummaryItem(s, "MMFdir", SpatialDomains::GeomMMFMap[m_MMFdir]);
    AddSummaryItem(s, "SmoothFactor", m_SmoothFactor);
    AddSummaryItem(s, "SFinit", m_SFinit);

    if (fabs(m_alpha - 1.0) > 0.001)
    {
        AddSummaryItem(s, "Alpha", m_alpha);
    }

    if (m_Incfreq > 0.0)
    {
        AddSummaryItem(s, "Incfreq", m_Incfreq);
    }

    if (m_MMFdir == 4)
    {
        AddSummaryItem(s, "MMFCircAxisX", m_MMFfactors[0]);
        AddSummaryItem(s, "MMFCircAxisY", m_MMFfactors[1]);
        AddSummaryItem(s, "MMFCircCentreX", m_MMFfactors[2]);
        AddSummaryItem(s, "MMFCircCentreY", m_MMFfactors[3]);
    }
}
}
}
