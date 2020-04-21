///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFR.cpp
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
// Description: LFR diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <SolverUtils/Diffusion/DiffusionLFR.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <MultiRegions/DisContField1D.h>
#include <boost/core/ignore_unused.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <iomanip>

using namespace std;

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLFR::type[] = {
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRDG", DiffusionLFR::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRSD", DiffusionLFR::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRHU", DiffusionLFR::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRcmin", DiffusionLFR::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRcinf", DiffusionLFR::create)};

        /**
         * @brief DiffusionLFR uses the Flux Reconstruction (FR) approach to
         * compute the diffusion term. The implementation is only for segments,
         * quadrilaterals and hexahedra at the moment.
         *
         * \todo Extension to triangles, tetrahedra and other shapes.
         * (Long term objective)
         */
        DiffusionLFR::DiffusionLFR(std::string diffType):m_diffType(diffType)
        {
        }

        /**
         * @brief Initiliase DiffusionLFR objects and store them before starting
         * the time-stepping.
         *
         * This routine calls the virtual functions #v_SetupMetrics and
         * #v_SetupCFunctions to initialise the objects needed by DiffusionLFR.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void DiffusionLFR::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            m_session = pSession;
            v_SetupMetrics(pSession, pFields);
            v_SetupCFunctions(pSession, pFields);

            // Initialising arrays
            int i, j;
            int nConvectiveFields = pFields.size();
            int nDim         = pFields[0]->GetCoordim(0);
            int nSolutionPts = pFields[0]->GetTotPoints();
            int nTracePts    = pFields[0]->GetTrace()->GetTotPoints();

            m_IF1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_DU1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_DFC1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_BD1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_D1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_DD1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_IF2 = Array<OneD, Array<OneD, NekDouble> >   (nConvectiveFields);
            m_DFC2 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_divFD = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
            m_divFC = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);

            m_tmp1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_tmp2 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);



            for (i = 0; i < nConvectiveFields; ++i)
            {
                m_IF1[i]  = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_DU1[i]  = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_DFC1[i] = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_BD1[i]  = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_D1[i]   = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_DD1[i]  = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_IF2[i]  = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_DFC2[i] = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_divFD[i]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                m_divFC[i]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);

                m_tmp1[i] = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_tmp2[i] = Array<OneD, Array<OneD, NekDouble> >(nDim);

                for (j = 0; j < nDim; ++j)
                {
                    m_IF1[i][j]  = Array<OneD, NekDouble>(nTracePts, 0.0);
                    m_DU1[i][j]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_DFC1[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_BD1[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_D1[i][j]   = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_DD1[i][j]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_DFC2[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);

                    m_tmp1[i][j]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_tmp2[i][j]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                }
            }

        }

        /**
         * @brief Setup the metric terms to compute the contravariant
         * fluxes. (i.e. this special metric terms transform the fluxes
         * at the interfaces of each element from the physical space to
         * the standard space).
         *
         * This routine calls the function #GetEdgeQFactors to compute and
         * store the metric factors following an anticlockwise conventions
         * along the edges/faces of each element. Note: for 1D problem
         * the transformation is not needed.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         *
         * \todo Add the metric terms for 3D Hexahedra.
         */
        void DiffusionLFR::v_SetupMetrics(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            boost::ignore_unused(pSession);

            int i, n;
            int nquad0, nquad1;
            int phys_offset;
            int nLocalSolutionPts;
            int nElements    = pFields[0]->GetExpSize();
            int nDimensions  = pFields[0]->GetCoordim(0);
            int nSolutionPts = pFields[0]->GetTotPoints();
            int nTracePts    = pFields[0]->GetTrace()->GetTotPoints();

            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDimensions);
            for (i = 0; i < nDimensions; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            pFields[0]->GetTrace()->GetNormals(m_traceNormals);

            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;

            m_jac  = Array<OneD, NekDouble>(nSolutionPts);

            Array<OneD, NekDouble> auxArray1;
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            LibUtilities::PointsKeyVector ptsKeys;

            switch (nDimensions)
            {
                case 1:
                {
                    for (n = 0; n < nElements; ++n)
                    {
                        ptsKeys = pFields[0]->GetExp(n)->GetPointsKeys();
                        nLocalSolutionPts = pFields[0]->GetExp(n)->GetTotPoints();
                        phys_offset = pFields[0]->GetPhys_Offset(n);
                        jac = pFields[0]->GetExp(n)
                                ->as<LocalRegions::Expansion1D>()->GetGeom1D()
                                ->GetMetricInfo()->GetJac(ptsKeys);
                        for (i = 0; i < nLocalSolutionPts; ++i)
                        {
                            m_jac[i+phys_offset] = jac[0];
                        }
                    }
                    break;
                }
                case 2:
                {
                    m_gmat = Array<OneD, Array<OneD, NekDouble> >(4);
                    m_gmat[0] = Array<OneD, NekDouble>(nSolutionPts);
                    m_gmat[1] = Array<OneD, NekDouble>(nSolutionPts);
                    m_gmat[2] = Array<OneD, NekDouble>(nSolutionPts);
                    m_gmat[3] = Array<OneD, NekDouble>(nSolutionPts);

                    m_Q2D_e0 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_Q2D_e1 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_Q2D_e2 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_Q2D_e3 = Array<OneD, Array<OneD, NekDouble> >(nElements);

                    for (n = 0; n < nElements; ++n)
                    {
                        base        = pFields[0]->GetExp(n)->GetBase();
                        nquad0      = base[0]->GetNumPoints();
                        nquad1      = base[1]->GetNumPoints();

                        m_Q2D_e0[n] = Array<OneD, NekDouble>(nquad0);
                        m_Q2D_e1[n] = Array<OneD, NekDouble>(nquad1);
                        m_Q2D_e2[n] = Array<OneD, NekDouble>(nquad0);
                        m_Q2D_e3[n] = Array<OneD, NekDouble>(nquad1);

                        // Extract the Q factors at each edge point
                        pFields[0]->GetExp(n)->GetEdgeQFactors(
                            0, auxArray1 = m_Q2D_e0[n]);
                        pFields[0]->GetExp(n)->GetEdgeQFactors(
                            1, auxArray1 = m_Q2D_e1[n]);
                        pFields[0]->GetExp(n)->GetEdgeQFactors(
                            2, auxArray1 = m_Q2D_e2[n]);
                        pFields[0]->GetExp(n)->GetEdgeQFactors(
                            3, auxArray1 = m_Q2D_e3[n]);

                        ptsKeys = pFields[0]->GetExp(n)->GetPointsKeys();
                        nLocalSolutionPts = pFields[0]->GetExp(n)->GetTotPoints();
                        phys_offset = pFields[0]->GetPhys_Offset(n);

                        jac  = pFields[0]->GetExp(n)
                                ->as<LocalRegions::Expansion2D>()->GetGeom2D()
                                ->GetMetricInfo()->GetJac(ptsKeys);
                        gmat = pFields[0]->GetExp(n)
                                ->as<LocalRegions::Expansion2D>()->GetGeom2D()
                                ->GetMetricInfo()->GetDerivFactors(ptsKeys);

                        if (pFields[0]->GetExp(n)
                                ->as<LocalRegions::Expansion2D>()->GetGeom2D()
                                ->GetMetricInfo()->GetGtype()
                            == SpatialDomains::eDeformed)
                        {
                            for (i = 0; i < nLocalSolutionPts; ++i)
                            {
                                m_jac[i+phys_offset]     = jac[i];
                                m_gmat[0][i+phys_offset] = gmat[0][i];
                                m_gmat[1][i+phys_offset] = gmat[1][i];
                                m_gmat[2][i+phys_offset] = gmat[2][i];
                                m_gmat[3][i+phys_offset] = gmat[3][i];
                            }
                        }
                        else
                        {
                            for (i = 0; i < nLocalSolutionPts; ++i)
                            {
                                m_jac[i+phys_offset]     = jac[0];
                                m_gmat[0][i+phys_offset] = gmat[0][0];
                                m_gmat[1][i+phys_offset] = gmat[1][0];
                                m_gmat[2][i+phys_offset] = gmat[2][0];
                                m_gmat[3][i+phys_offset] = gmat[3][0];
                            }
                        }
                    }
                    break;
                }
                case 3:
                {
                    ASSERTL0(false,"3DFR Metric terms not implemented yet");
                    break;
                }
                default:
                {
                    ASSERTL0(false, "Expansion dimension not recognised");
                    break;
                }
            }
        }

        /**
         * @brief Setup the derivatives of the correction functions. For more
         * details see J Sci Comput (2011) 47: 50–72
         *
         * This routine calls 3 different bases:
         *      #eDG_DG_Left - #eDG_DG_Left which recovers a nodal DG scheme,
         *      #eDG_SD_Left - #eDG_SD_Left which recovers the SD scheme,
         *      #eDG_HU_Left - #eDG_HU_Left which recovers the Huynh scheme.
         * The values of the derivatives of the correction function are then
         * stored into global variables and reused into the virtual functions
         * #v_DivCFlux_1D, #v_DivCFlux_2D, #v_DivCFlux_3D to compute the
         * the divergence of the correction flux for 1D, 2D or 3D problems.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void DiffusionLFR::v_SetupCFunctions(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            boost::ignore_unused(pSession);

            int i, n;
            NekDouble c0 = 0.0;
            NekDouble c1 = 0.0;
            NekDouble c2 = 0.0;
            int nquad0, nquad1, nquad2;
            int nmodes0, nmodes1, nmodes2;
            Array<OneD, LibUtilities::BasisSharedPtr> base;

            int nElements = pFields[0]->GetExpSize();
            int nDim      = pFields[0]->GetCoordim(0);

            switch (nDim)
            {
                case 1:
                {
                    m_dGL_xi1 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi1 = Array<OneD, Array<OneD, NekDouble> >(nElements);

                    for (n = 0; n < nElements; ++n)
                    {
                        base    = pFields[0]->GetExp(n)->GetBase();
                        nquad0  = base[0]->GetNumPoints();
                        nmodes0 = base[0]->GetNumModes();
                        Array<OneD, const NekDouble> z0;
                        Array<OneD, const NekDouble> w0;

                        base[0]->GetZW(z0, w0);

                        m_dGL_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGR_xi1[n] = Array<OneD, NekDouble>(nquad0);

                        // Auxiliary vectors to build up the auxiliary
                        // derivatives of the Legendre polynomials
                        Array<OneD,NekDouble> dLp0   (nquad0, 0.0);
                        Array<OneD,NekDouble> dLpp0  (nquad0, 0.0);
                        Array<OneD,NekDouble> dLpm0  (nquad0, 0.0);

                        // Degree of the correction functions
                        int p0 = nmodes0 - 1;

                        // Function sign to compute  left correction function
                        NekDouble sign0 = pow(-1.0, p0);

                        // Factorial factor to build the scheme
                        NekDouble ap0 = boost::math::tgamma(2 * p0 + 1)
                        / (pow(2.0, p0)
                           * boost::math::tgamma(p0 + 1)
                           * boost::math::tgamma(p0 + 1));

                        // Scalar parameter which recovers the FR schemes
                        if (m_diffType == "LFRDG")
                        {
                            c0 = 0.0;
                        }
                        else if (m_diffType == "LFRSD")
                        {
                            c0 = 2.0 * p0 / ((2.0 * p0 + 1.0) * (p0 + 1.0)
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_diffType == "LFRHU")
                        {
                            c0 = 2.0 * (p0 + 1.0) / ((2.0 * p0 + 1.0) * p0
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_diffType == "LFRcmin")
                        {
                            c0 = -2.0 / ((2.0 * p0 + 1.0)
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_diffType == "LFRinf")
                        {
                            c0 = 10000000000000000.0;
                        }

                        NekDouble etap0 = 0.5 * c0 * (2.0 * p0 + 1.0)
                        * (ap0 * boost::math::tgamma(p0 + 1))
                        * (ap0 * boost::math::tgamma(p0 + 1));

                        NekDouble overeta0 = 1.0 / (1.0 + etap0);

                        // Derivative of the Legendre polynomials
                        // dLp  = derivative of the Legendre polynomial of order p
                        // dLpp = derivative of the Legendre polynomial of order p+1
                        // dLpm = derivative of the Legendre polynomial of order p-1
                        Polylib::jacobd(nquad0, z0.data(), &(dLp0[0]),  p0,   0.0, 0.0);
                        Polylib::jacobd(nquad0, z0.data(), &(dLpp0[0]), p0+1, 0.0, 0.0);
                        Polylib::jacobd(nquad0, z0.data(), &(dLpm0[0]), p0-1, 0.0, 0.0);

                        // Building the DG_c_Left
                        for(i = 0; i < nquad0; ++i)
                        {
                            m_dGL_xi1[n][i]  = etap0 * dLpm0[i];
                            m_dGL_xi1[n][i] += dLpp0[i];
                            m_dGL_xi1[n][i] *= overeta0;
                            m_dGL_xi1[n][i]  = dLp0[i] - m_dGL_xi1[n][i];
                            m_dGL_xi1[n][i]  = 0.5 * sign0 * m_dGL_xi1[n][i];
                        }

                        // Building the DG_c_Right
                        for(i = 0; i < nquad0; ++i)
                        {
                            m_dGR_xi1[n][i]  = etap0 * dLpm0[i];
                            m_dGR_xi1[n][i] += dLpp0[i];
                            m_dGR_xi1[n][i] *= overeta0;
                            m_dGR_xi1[n][i] += dLp0[i];
                            m_dGR_xi1[n][i]  = 0.5 * m_dGR_xi1[n][i];
                        }
                    }
                    break;
                }
                case 2:
                {
                    LibUtilities::BasisSharedPtr BasisFR_Left0;
                    LibUtilities::BasisSharedPtr BasisFR_Right0;
                    LibUtilities::BasisSharedPtr BasisFR_Left1;
                    LibUtilities::BasisSharedPtr BasisFR_Right1;

                    m_dGL_xi1 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi1 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGL_xi2 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi2 = Array<OneD, Array<OneD, NekDouble> >(nElements);

                    for (n = 0; n < nElements; ++n)
                    {
                        base      = pFields[0]->GetExp(n)->GetBase();
                        nquad0    = base[0]->GetNumPoints();
                        nquad1    = base[1]->GetNumPoints();
                        nmodes0   = base[0]->GetNumModes();
                        nmodes1   = base[1]->GetNumModes();

                        Array<OneD, const NekDouble> z0;
                        Array<OneD, const NekDouble> w0;
                        Array<OneD, const NekDouble> z1;
                        Array<OneD, const NekDouble> w1;

                        base[0]->GetZW(z0, w0);
                        base[1]->GetZW(z1, w1);

                        m_dGL_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGR_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGL_xi2[n] = Array<OneD, NekDouble>(nquad1);
                        m_dGR_xi2[n] = Array<OneD, NekDouble>(nquad1);

                        // Auxiliary vectors to build up the auxiliary
                        // derivatives of the Legendre polynomials
                        Array<OneD,NekDouble> dLp0   (nquad0, 0.0);
                        Array<OneD,NekDouble> dLpp0  (nquad0, 0.0);
                        Array<OneD,NekDouble> dLpm0  (nquad0, 0.0);
                        Array<OneD,NekDouble> dLp1   (nquad1, 0.0);
                        Array<OneD,NekDouble> dLpp1  (nquad1, 0.0);
                        Array<OneD,NekDouble> dLpm1  (nquad1, 0.0);

                        // Degree of the correction functions
                        int p0 = nmodes0 - 1;
                        int p1 = nmodes1 - 1;

                        // Function sign to compute  left correction function
                        NekDouble sign0 = pow(-1.0, p0);
                        NekDouble sign1 = pow(-1.0, p1);

                        // Factorial factor to build the scheme
                        NekDouble ap0 = boost::math::tgamma(2 * p0 + 1)
                        / (pow(2.0, p0)
                           * boost::math::tgamma(p0 + 1)
                           * boost::math::tgamma(p0 + 1));

                        NekDouble ap1 = boost::math::tgamma(2 * p1 + 1)
                        / (pow(2.0, p1)
                           * boost::math::tgamma(p1 + 1)
                           * boost::math::tgamma(p1 + 1));

                        // Scalar parameter which recovers the FR schemes
                        if (m_diffType == "LFRDG")
                        {
                            c0 = 0.0;
                            c1 = 0.0;
                        }
                        else if (m_diffType == "LFRSD")
                        {
                            c0 = 2.0 * p0 / ((2.0 * p0 + 1.0) * (p0 + 1.0)
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = 2.0 * p1 / ((2.0 * p1 + 1.0) * (p1 + 1.0)
                                        * (ap1 * boost::math::tgamma(p1 + 1))
                                        * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_diffType == "LFRHU")
                        {
                            c0 = 2.0 * (p0 + 1.0) / ((2.0 * p0 + 1.0) * p0
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = 2.0 * (p1 + 1.0) / ((2.0 * p1 + 1.0) * p1
                                        * (ap1 * boost::math::tgamma(p1 + 1))
                                        * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_diffType == "LFRcmin")
                        {
                            c0 = -2.0 / ((2.0 * p0 + 1.0)
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = -2.0 / ((2.0 * p1 + 1.0)
                                        * (ap1 * boost::math::tgamma(p1 + 1))
                                        * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_diffType == "LFRinf")
                        {
                            c0 = 10000000000000000.0;
                            c1 = 10000000000000000.0;
                        }

                        NekDouble etap0 = 0.5 * c0 * (2.0 * p0 + 1.0)
                        * (ap0 * boost::math::tgamma(p0 + 1))
                        * (ap0 * boost::math::tgamma(p0 + 1));

                        NekDouble etap1 = 0.5 * c1 * (2.0 * p1 + 1.0)
                        * (ap1 * boost::math::tgamma(p1 + 1))
                        * (ap1 * boost::math::tgamma(p1 + 1));

                        NekDouble overeta0 = 1.0 / (1.0 + etap0);
                        NekDouble overeta1 = 1.0 / (1.0 + etap1);

                        // Derivative of the Legendre polynomials
                        // dLp  = derivative of the Legendre polynomial of order p
                        // dLpp = derivative of the Legendre polynomial of order p+1
                        // dLpm = derivative of the Legendre polynomial of order p-1
                        Polylib::jacobd(nquad0, z0.data(), &(dLp0[0]),  p0,   0.0, 0.0);
                        Polylib::jacobd(nquad0, z0.data(), &(dLpp0[0]), p0+1, 0.0, 0.0);
                        Polylib::jacobd(nquad0, z0.data(), &(dLpm0[0]), p0-1, 0.0, 0.0);

                        Polylib::jacobd(nquad1, z1.data(), &(dLp1[0]),  p1,   0.0, 0.0);
                        Polylib::jacobd(nquad1, z1.data(), &(dLpp1[0]), p1+1, 0.0, 0.0);
                        Polylib::jacobd(nquad1, z1.data(), &(dLpm1[0]), p1-1, 0.0, 0.0);

                        // Building the DG_c_Left
                        for(i = 0; i < nquad0; ++i)
                        {
                            m_dGL_xi1[n][i]  = etap0 * dLpm0[i];
                            m_dGL_xi1[n][i] += dLpp0[i];
                            m_dGL_xi1[n][i] *= overeta0;
                            m_dGL_xi1[n][i]  = dLp0[i] - m_dGL_xi1[n][i];
                            m_dGL_xi1[n][i]  = 0.5 * sign0 * m_dGL_xi1[n][i];
                        }

                        // Building the DG_c_Left
                        for(i = 0; i < nquad1; ++i)
                        {
                            m_dGL_xi2[n][i]  = etap1 * dLpm1[i];
                            m_dGL_xi2[n][i] += dLpp1[i];
                            m_dGL_xi2[n][i] *= overeta1;
                            m_dGL_xi2[n][i]  = dLp1[i] - m_dGL_xi2[n][i];
                            m_dGL_xi2[n][i]  = 0.5 * sign1 * m_dGL_xi2[n][i];
                        }

                        // Building the DG_c_Right
                        for(i = 0; i < nquad0; ++i)
                        {
                            m_dGR_xi1[n][i]  = etap0 * dLpm0[i];
                            m_dGR_xi1[n][i] += dLpp0[i];
                            m_dGR_xi1[n][i] *= overeta0;
                            m_dGR_xi1[n][i] += dLp0[i];
                            m_dGR_xi1[n][i]  = 0.5 * m_dGR_xi1[n][i];
                        }

                        // Building the DG_c_Right
                        for(i = 0; i < nquad1; ++i)
                        {
                            m_dGR_xi2[n][i]  = etap1 * dLpm1[i];
                            m_dGR_xi2[n][i] += dLpp1[i];
                            m_dGR_xi2[n][i] *= overeta1;
                            m_dGR_xi2[n][i] += dLp1[i];
                            m_dGR_xi2[n][i]  = 0.5 * m_dGR_xi2[n][i];
                        }
                    }
                    break;
                }
                case 3:
                {
                    m_dGL_xi1 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi1 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGL_xi2 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi2 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGL_xi3 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi3 = Array<OneD, Array<OneD, NekDouble> >(nElements);

                    for (n = 0; n < nElements; ++n)
                    {
                        base      = pFields[0]->GetExp(n)->GetBase();
                        nquad0    = base[0]->GetNumPoints();
                        nquad1    = base[1]->GetNumPoints();
                        nquad2    = base[2]->GetNumPoints();
                        nmodes0   = base[0]->GetNumModes();
                        nmodes1   = base[1]->GetNumModes();
                        nmodes2   = base[2]->GetNumModes();

                        Array<OneD, const NekDouble> z0;
                        Array<OneD, const NekDouble> w0;
                        Array<OneD, const NekDouble> z1;
                        Array<OneD, const NekDouble> w1;
                        Array<OneD, const NekDouble> z2;
                        Array<OneD, const NekDouble> w2;

                        base[0]->GetZW(z0, w0);
                        base[1]->GetZW(z1, w1);
                        base[1]->GetZW(z2, w2);

                        m_dGL_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGR_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGL_xi2[n] = Array<OneD, NekDouble>(nquad1);
                        m_dGR_xi2[n] = Array<OneD, NekDouble>(nquad1);
                        m_dGL_xi3[n] = Array<OneD, NekDouble>(nquad2);
                        m_dGR_xi3[n] = Array<OneD, NekDouble>(nquad2);

                        // Auxiliary vectors to build up the auxiliary
                        // derivatives of the Legendre polynomials
                        Array<OneD,NekDouble> dLp0   (nquad0, 0.0);
                        Array<OneD,NekDouble> dLpp0  (nquad0, 0.0);
                        Array<OneD,NekDouble> dLpm0  (nquad0, 0.0);
                        Array<OneD,NekDouble> dLp1   (nquad1, 0.0);
                        Array<OneD,NekDouble> dLpp1  (nquad1, 0.0);
                        Array<OneD,NekDouble> dLpm1  (nquad1, 0.0);
                        Array<OneD,NekDouble> dLp2   (nquad2, 0.0);
                        Array<OneD,NekDouble> dLpp2  (nquad2, 0.0);
                        Array<OneD,NekDouble> dLpm2  (nquad2, 0.0);

                        // Degree of the correction functions
                        int p0 = nmodes0 - 1;
                        int p1 = nmodes1 - 1;
                        int p2 = nmodes2 - 1;

                        // Function sign to compute  left correction function
                        NekDouble sign0 = pow(-1.0, p0);
                        NekDouble sign1 = pow(-1.0, p1);
                        NekDouble sign2 = pow(-1.0, p2);

                        // Factorial factor to build the scheme
                        NekDouble ap0 = boost::math::tgamma(2 * p0 + 1)
                        / (pow(2.0, p0)
                           * boost::math::tgamma(p0 + 1)
                           * boost::math::tgamma(p0 + 1));

                        // Factorial factor to build the scheme
                        NekDouble ap1 = boost::math::tgamma(2 * p1 + 1)
                        / (pow(2.0, p1)
                           * boost::math::tgamma(p1 + 1)
                           * boost::math::tgamma(p1 + 1));

                        // Factorial factor to build the scheme
                        NekDouble ap2 = boost::math::tgamma(2 * p2 + 1)
                        / (pow(2.0, p2)
                           * boost::math::tgamma(p2 + 1)
                           * boost::math::tgamma(p2 + 1));

                        // Scalar parameter which recovers the FR schemes
                        if (m_diffType == "LFRDG")
                        {
                            c0 = 0.0;
                            c1 = 0.0;
                            c2 = 0.0;
                        }
                        else if (m_diffType == "LFRSD")
                        {
                            c0 = 2.0 * p0 / ((2.0 * p0 + 1.0) * (p0 + 1.0)
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = 2.0 * p1 / ((2.0 * p1 + 1.0) * (p1 + 1.0)
                                        * (ap1 * boost::math::tgamma(p1 + 1))
                                        * (ap1 * boost::math::tgamma(p1 + 1)));

                            c2 = 2.0 * p2 / ((2.0 * p2 + 1.0) * (p2 + 1.0)
                                        * (ap2 * boost::math::tgamma(p2 + 1))
                                        * (ap2 * boost::math::tgamma(p2 + 1)));
                        }
                        else if (m_diffType == "LFRHU")
                        {
                            c0 = 2.0 * (p0 + 1.0) / ((2.0 * p0 + 1.0) * p0
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = 2.0 * (p1 + 1.0) / ((2.0 * p1 + 1.0) * p1
                                        * (ap1 * boost::math::tgamma(p1 + 1))
                                        * (ap1 * boost::math::tgamma(p1 + 1)));

                            c2 = 2.0 * (p2 + 1.0) / ((2.0 * p2 + 1.0) * p2
                                        * (ap2 * boost::math::tgamma(p2 + 1))
                                        * (ap2 * boost::math::tgamma(p2 + 1)));
                        }
                        else if (m_diffType == "LFRcmin")
                        {
                            c0 = -2.0 / ((2.0 * p0 + 1.0)
                                        * (ap0 * boost::math::tgamma(p0 + 1))
                                        * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = -2.0 / ((2.0 * p1 + 1.0)
                                        * (ap1 * boost::math::tgamma(p1 + 1))
                                        * (ap1 * boost::math::tgamma(p1 + 1)));

                            c2 = -2.0 / ((2.0 * p2 + 1.0)
                                        * (ap2 * boost::math::tgamma(p2 + 1))
                                        * (ap2 * boost::math::tgamma(p2 + 1)));
                        }
                        else if (m_diffType == "LFRinf")
                        {
                            c0 = 10000000000000000.0;
                            c1 = 10000000000000000.0;
                            c2 = 10000000000000000.0;
                        }

                        NekDouble etap0 = 0.5 * c0 * (2.0 * p0 + 1.0)
                        * (ap0 * boost::math::tgamma(p0 + 1))
                        * (ap0 * boost::math::tgamma(p0 + 1));

                        NekDouble etap1 = 0.5 * c1 * (2.0 * p1 + 1.0)
                        * (ap1 * boost::math::tgamma(p1 + 1))
                        * (ap1 * boost::math::tgamma(p1 + 1));

                        NekDouble etap2 = 0.5 * c2 * (2.0 * p2 + 1.0)
                        * (ap2 * boost::math::tgamma(p2 + 1))
                        * (ap2 * boost::math::tgamma(p2 + 1));

                        NekDouble overeta0 = 1.0 / (1.0 + etap0);
                        NekDouble overeta1 = 1.0 / (1.0 + etap1);
                        NekDouble overeta2 = 1.0 / (1.0 + etap2);

                        // Derivative of the Legendre polynomials
                        // dLp  = derivative of the Legendre polynomial of order p
                        // dLpp = derivative of the Legendre polynomial of order p+1
                        // dLpm = derivative of the Legendre polynomial of order p-1
                        Polylib::jacobd(nquad0, z0.data(), &(dLp0[0]),  p0,   0.0, 0.0);
                        Polylib::jacobd(nquad0, z0.data(), &(dLpp0[0]), p0+1, 0.0, 0.0);
                        Polylib::jacobd(nquad0, z0.data(), &(dLpm0[0]), p0-1, 0.0, 0.0);

                        Polylib::jacobd(nquad1, z1.data(), &(dLp1[0]),  p1,   0.0, 0.0);
                        Polylib::jacobd(nquad1, z1.data(), &(dLpp1[0]), p1+1, 0.0, 0.0);
                        Polylib::jacobd(nquad1, z1.data(), &(dLpm1[0]), p1-1, 0.0, 0.0);

                        Polylib::jacobd(nquad2, z2.data(), &(dLp2[0]),  p2,   0.0, 0.0);
                        Polylib::jacobd(nquad2, z2.data(), &(dLpp2[0]), p2+1, 0.0, 0.0);
                        Polylib::jacobd(nquad2, z2.data(), &(dLpm2[0]), p2-1, 0.0, 0.0);

                        // Building the DG_c_Left
                        for(i = 0; i < nquad0; ++i)
                        {
                            m_dGL_xi1[n][i]  = etap0 * dLpm0[i];
                            m_dGL_xi1[n][i] += dLpp0[i];
                            m_dGL_xi1[n][i] *= overeta0;
                            m_dGL_xi1[n][i]  = dLp0[i] - m_dGL_xi1[n][i];
                            m_dGL_xi1[n][i]  = 0.5 * sign0 * m_dGL_xi1[n][i];
                        }

                        // Building the DG_c_Left
                        for(i = 0; i < nquad1; ++i)
                        {
                            m_dGL_xi2[n][i]  = etap1 * dLpm1[i];
                            m_dGL_xi2[n][i] += dLpp1[i];
                            m_dGL_xi2[n][i] *= overeta1;
                            m_dGL_xi2[n][i]  = dLp1[i] - m_dGL_xi2[n][i];
                            m_dGL_xi2[n][i]  = 0.5 * sign1 * m_dGL_xi2[n][i];
                        }

                        // Building the DG_c_Left
                        for(i = 0; i < nquad2; ++i)
                        {
                            m_dGL_xi3[n][i]  = etap2 * dLpm2[i];
                            m_dGL_xi3[n][i] += dLpp2[i];
                            m_dGL_xi3[n][i] *= overeta2;
                            m_dGL_xi3[n][i]  = dLp2[i] - m_dGL_xi3[n][i];
                            m_dGL_xi3[n][i]  = 0.5 * sign2 * m_dGL_xi3[n][i];
                        }

                        // Building the DG_c_Right
                        for(i = 0; i < nquad0; ++i)
                        {
                            m_dGR_xi1[n][i]  = etap0 * dLpm0[i];
                            m_dGR_xi1[n][i] += dLpp0[i];
                            m_dGR_xi1[n][i] *= overeta0;
                            m_dGR_xi1[n][i] += dLp0[i];
                            m_dGR_xi1[n][i]  = 0.5 * m_dGR_xi1[n][i];
                        }

                        // Building the DG_c_Right
                        for(i = 0; i < nquad1; ++i)
                        {
                            m_dGR_xi2[n][i]  = etap1 * dLpm1[i];
                            m_dGR_xi2[n][i] += dLpp1[i];
                            m_dGR_xi2[n][i] *= overeta1;
                            m_dGR_xi2[n][i] += dLp1[i];
                            m_dGR_xi2[n][i]  = 0.5 * m_dGR_xi2[n][i];
                        }

                        // Building the DG_c_Right
                        for(i = 0; i < nquad2; ++i)
                        {
                            m_dGR_xi3[n][i]  = etap2 * dLpm2[i];
                            m_dGR_xi3[n][i] += dLpp2[i];
                            m_dGR_xi3[n][i] *= overeta2;
                            m_dGR_xi3[n][i] += dLp2[i];
                            m_dGR_xi3[n][i]  = 0.5 * m_dGR_xi3[n][i];
                        }
                    }
                    break;
                }
                default:
                {
                    ASSERTL0(false,"Expansion dimension not recognised");
                    break;
                }
            }
        }

        /**
         * @brief Calculate FR Diffusion for the linear problems
         * using an LDG interface flux.
         *
         */
        void DiffusionLFR::v_Diffuse(
            const std::size_t                                 nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            boost::ignore_unused(pFwd, pBwd);

            int i, j, n;
            int phys_offset;
            //Array<TwoD, const NekDouble> gmat;
            //Array<OneD, const NekDouble> jac;
            Array<OneD,       NekDouble> auxArray1, auxArray2;

            Array<OneD, LibUtilities::BasisSharedPtr> Basis;
            Basis = fields[0]->GetExp(0)->GetBase();

            int nElements    = fields[0]->GetExpSize();
            int nDim         = fields[0]->GetCoordim(0);
            int nSolutionPts = fields[0]->GetTotPoints();
            int nCoeffs      = fields[0]->GetNcoeffs();

            Array<OneD, Array<OneD, NekDouble> > outarrayCoeff(nConvectiveFields);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                outarrayCoeff[i]  = Array<OneD, NekDouble>(nCoeffs);
            }

            // Compute interface numerical fluxes for inarray in physical space
            v_NumFluxforScalar(fields, inarray, m_IF1);

            switch(nDim)
            {
                // 1D problems
                case 1:
                {
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        // Computing the physical first-order discountinuous
                        // derivative
                        for (n = 0; n < nElements; n++)
                        {
                            phys_offset = fields[0]->GetPhys_Offset(n);

                            fields[i]->GetExp(n)->PhysDeriv(0,
                                auxArray1 = inarray[i] + phys_offset,
                                auxArray2 = m_DU1[i][0] + phys_offset);
                        }

                        // Computing the standard first-order correction
                        // derivative
                        v_DerCFlux_1D(nConvectiveFields, fields, inarray[i],
                                      m_IF1[i][0], m_DFC1[i][0]);

                        // Back to the physical space using global operations
                        Vmath::Vdiv(nSolutionPts, &m_DFC1[i][0][0], 1,
                                    &m_jac[0], 1, &m_DFC1[i][0][0], 1);
                        Vmath::Vdiv(nSolutionPts, &m_DFC1[i][0][0], 1,
                                    &m_jac[0], 1, &m_DFC1[i][0][0], 1);

                        // Computing total first order derivatives
                        Vmath::Vadd(nSolutionPts, &m_DFC1[i][0][0], 1,
                                    &m_DU1[i][0][0], 1, &m_D1[i][0][0], 1);

                        Vmath::Vcopy(nSolutionPts, &m_D1[i][0][0], 1,
                                     &m_tmp1[i][0][0], 1);
                    }

                    // Computing interface numerical fluxes for m_D1
                    // in physical space
                    v_NumFluxforVector(fields, inarray, m_tmp1, m_IF2);

                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        // Computing the physical second-order discountinuous
                        // derivative
                        for (n = 0; n < nElements; n++)
                        {
                            phys_offset = fields[0]->GetPhys_Offset(n);

                            fields[i]->GetExp(n)->PhysDeriv(0,
                                auxArray1 = m_D1[i][0] + phys_offset,
                                auxArray2 = m_DD1[i][0] + phys_offset);
                        }

                        // Computing the standard second-order correction
                        // derivative
                        v_DerCFlux_1D(nConvectiveFields, fields, m_D1[i][0],
                                      m_IF2[i], m_DFC2[i][0]);

                        // Back to the physical space using global operations
                        Vmath::Vdiv(nSolutionPts, &m_DFC2[i][0][0], 1,
                                    &m_jac[0], 1, &m_DFC2[i][0][0], 1);
                        Vmath::Vdiv(nSolutionPts, &m_DFC2[i][0][0], 1,
                                    &m_jac[0], 1, &m_DFC2[i][0][0], 1);

                        // Adding the total divergence to outarray (RHS)
                        Vmath::Vadd(nSolutionPts, &m_DFC2[i][0][0], 1,
                                    &m_DD1[i][0][0], 1, &outarray[i][0], 1);

                        // Primitive Dealiasing 1D
                        if (!(Basis[0]->Collocation()))
                        {
                            fields[i]->FwdTrans(outarray[i], outarrayCoeff[i]);
                            fields[i]->BwdTrans(outarrayCoeff[i], outarray[i]);
                        }
                    }
                    break;
                }
                // 2D problems
                case 2:
                {
                    for(i = 0; i < nConvectiveFields; ++i)
                    {
                        for (j = 0; j < nDim; ++j)
                        {
                            // Temporary vectors
                            Array<OneD, NekDouble> u1_hat(nSolutionPts, 0.0);
                            Array<OneD, NekDouble> u2_hat(nSolutionPts, 0.0);

                            if (j == 0)
                            {
                                Vmath::Vmul(nSolutionPts, &inarray[i][0], 1,
                                            &m_gmat[0][0], 1, &u1_hat[0], 1);

                                Vmath::Vmul(nSolutionPts, &u1_hat[0], 1,
                                            &m_jac[0], 1, &u1_hat[0], 1);

                                Vmath::Vmul(nSolutionPts, &inarray[i][0], 1,
                                            &m_gmat[1][0], 1, &u2_hat[0], 1);

                                Vmath::Vmul(nSolutionPts, &u2_hat[0], 1,
                                            &m_jac[0], 1, &u2_hat[0], 1);
                            }
                            else if (j == 1)
                            {
                                Vmath::Vmul(nSolutionPts, &inarray[i][0], 1,
                                            &m_gmat[2][0], 1, &u1_hat[0], 1);

                                Vmath::Vmul(nSolutionPts, &u1_hat[0], 1,
                                            &m_jac[0], 1, &u1_hat[0], 1);

                                Vmath::Vmul(nSolutionPts, &inarray[i][0], 1,
                                            &m_gmat[3][0], 1, &u2_hat[0], 1);

                                Vmath::Vmul(nSolutionPts, &u2_hat[0], 1,
                                            &m_jac[0], 1, &u2_hat[0], 1);
                            }

                            for (n = 0; n < nElements; n++)
                            {
                                phys_offset = fields[0]->GetPhys_Offset(n);

                                fields[i]->GetExp(n)->StdPhysDeriv(0,
                                    auxArray1 = u1_hat + phys_offset,
                                    auxArray2 = m_tmp1[i][j] + phys_offset);

                                fields[i]->GetExp(n)->StdPhysDeriv(1,
                                    auxArray1 = u2_hat + phys_offset,
                                    auxArray2 = m_tmp2[i][j] + phys_offset);
                            }

                            Vmath::Vadd(nSolutionPts, &m_tmp1[i][j][0], 1,
                                        &m_tmp2[i][j][0], 1,
                                        &m_DU1[i][j][0], 1);

                            // Divide by the metric jacobian
                            Vmath::Vdiv(nSolutionPts, &m_DU1[i][j][0], 1,
                                        &m_jac[0], 1, &m_DU1[i][j][0], 1);

                            // Computing the standard first-order correction
                            // derivatives
                            v_DerCFlux_2D(nConvectiveFields, j, fields,
                                          inarray[i], m_IF1[i][j],
                                          m_DFC1[i][j]);
                        }

                        // Multiplying derivatives by B matrix to get auxiliary
                        // variables
                        for (j = 0; j < nSolutionPts; ++j)
                        {
                            // std(q1)
                            m_BD1[i][0][j] =
                            (m_gmat[0][j]*m_gmat[0][j] +
                             m_gmat[2][j]*m_gmat[2][j]) *
                            m_DFC1[i][0][j] +
                            (m_gmat[1][j]*m_gmat[0][j] +
                             m_gmat[3][j]*m_gmat[2][j]) *
                            m_DFC1[i][1][j];

                            // std(q2)
                            m_BD1[i][1][j] =
                            (m_gmat[1][j]*m_gmat[0][j] +
                             m_gmat[3][j]*m_gmat[2][j]) *
                            m_DFC1[i][0][j] +
                            (m_gmat[1][j]*m_gmat[1][j] +
                             m_gmat[3][j]*m_gmat[3][j]) *
                            m_DFC1[i][1][j];
                        }

                        // Multiplying derivatives by A^(-1) to get back
                        // into the physical space
                        for (j = 0; j < nSolutionPts; j++)
                        {
                            // q1 = A11^(-1)*std(q1) + A12^(-1)*std(q2)
                            m_DFC1[i][0][j] =
                            m_gmat[3][j] * m_BD1[i][0][j] -
                            m_gmat[2][j] * m_BD1[i][1][j];

                            // q2 = A21^(-1)*std(q1) + A22^(-1)*std(q2)
                            m_DFC1[i][1][j] =
                            m_gmat[0][j] * m_BD1[i][1][j] -
                            m_gmat[1][j] * m_BD1[i][0][j];
                        }

                        // Computing the physical first-order derivatives
                        for (j = 0; j < nDim; ++j)
                        {
                            Vmath::Vadd(nSolutionPts, &m_DU1[i][j][0], 1,
                                        &m_DFC1[i][j][0], 1,
                                        &m_D1[i][j][0], 1);
                        }
                    }

                    // Computing interface numerical fluxes for m_D1
                    // in physical space
                    v_NumFluxforVector(fields, inarray, m_D1, m_IF2);

                    // Computing the standard second-order discontinuous
                    // derivatives
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        Array<OneD, NekDouble> f_hat(nSolutionPts, 0.0);
                        Array<OneD, NekDouble> g_hat(nSolutionPts, 0.0);

                        for (j = 0; j < nSolutionPts; j++)
                        {
                            f_hat[j] = (m_D1[i][0][j] * m_gmat[0][j] +
                            m_D1[i][1][j] * m_gmat[2][j])*m_jac[j];

                            g_hat[j] = (m_D1[i][0][j] * m_gmat[1][j] +
                            m_D1[i][1][j] * m_gmat[3][j])*m_jac[j];
                        }

                        for (n = 0; n < nElements; n++)
                        {
                            phys_offset = fields[0]->GetPhys_Offset(n);

                            fields[0]->GetExp(n)->StdPhysDeriv(0,
                            auxArray1 = f_hat + phys_offset,
                            auxArray2 = m_DD1[i][0] + phys_offset);

                            fields[0]->GetExp(n)->StdPhysDeriv(1,
                            auxArray1 = g_hat + phys_offset,
                            auxArray2 = m_DD1[i][1] + phys_offset);
                        }

                        // Divergence of the standard discontinuous flux
                        Vmath::Vadd(nSolutionPts,
                                    &m_DD1[i][0][0], 1,
                                    &m_DD1[i][1][0], 1,
                                    &m_divFD[i][0], 1);

                        // Divergence of the standard correction flux
                        if (Basis[0]->GetPointsType() ==
                            LibUtilities::eGaussGaussLegendre &&
                            Basis[1]->GetPointsType() ==
                            LibUtilities::eGaussGaussLegendre)
                        {
                            v_DivCFlux_2D_Gauss(nConvectiveFields, fields,
                                                f_hat, g_hat,
                                                m_IF2[i], m_divFC[i]);
                        }
                        else
                        {
                            v_DivCFlux_2D(nConvectiveFields, fields,
                                          m_D1[i][0], m_D1[i][1],
                                          m_IF2[i], m_divFC[i]);
                        }

                        // Divergence of the standard final flux
                        Vmath::Vadd(nSolutionPts, &m_divFD[i][0], 1,
                                    &m_divFC[i][0], 1, &outarray[i][0], 1);

                        // Dividing by the jacobian to get back into
                        // physical space
                        Vmath::Vdiv(nSolutionPts, &outarray[i][0], 1,
                                    &m_jac[0], 1, &outarray[i][0], 1);

                        // Primitive Dealiasing 2D
                        if (!(Basis[0]->Collocation()))
                        {
                            fields[i]->FwdTrans(outarray[i], outarrayCoeff[i]);
                            fields[i]->BwdTrans(outarrayCoeff[i], outarray[i]);
                        }
                    }//Close loop on nConvectiveFields
                    break;
                }
                // 3D problems
                case 3:
                {
                    ASSERTL0(false, "3D FRDG case not implemented yet");
                    break;
                }
            }
        }

        /**
         * @brief Builds the numerical flux for the 1st order derivatives
         *
         */
        void DiffusionLFR::v_NumFluxforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nvariables = fields.size();
            int nDim       = fields[0]->GetCoordim(0);

            Array<OneD, NekDouble > Fwd     (nTracePts);
            Array<OneD, NekDouble > Bwd     (nTracePts);
            Array<OneD, NekDouble > Vn      (nTracePts, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePts, 0.0);

            // Get the normal velocity Vn
            for (i = 0; i < nDim; ++i)
            {
                Vmath::Svtvp(nTracePts, 1.0, m_traceNormals[i], 1,
                             Vn, 1, Vn, 1);
            }

            // Get the sign of (v \cdot n), v = an arbitrary vector
            // Evaluate upwind flux:
            // uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nvariables ; ++i)
                {
                    // Compute Fwd and Bwd value of ufield of i direction
                    fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);

                    // if Vn >= 0, flux = uFwd, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd

                    // else if Vn < 0, flux = uBwd, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

                    fields[i]->GetTrace()->Upwind(Vn, Fwd, Bwd, fluxtemp);

                    // Imposing weak boundary condition with flux
                    // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd

                    // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd

                    if(fields[0]->GetBndCondExpansions().size())
                    {
                        v_WeakPenaltyforScalar(fields, i, ufield[i], fluxtemp);
                    }

                    // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n}),
                    // i.e,
                    // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n})
                    // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n})

                    // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n}),
                    // i.e,
                    // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n})
                    // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n})

                    Vmath::Vcopy(nTracePts, fluxtemp, 1, uflux[i][j], 1);
                }
            }
        }

        /**
         * @brief Imposes appropriate bcs for the 1st order derivatives
         *
         */
        void DiffusionLFR::v_WeakPenaltyforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const Array<OneD, const NekDouble>                &ufield,
                  Array<OneD,       NekDouble>                &penaltyflux)
        {
            // Variable initialisation
            int i, e, id1, id2;
            int nBndEdgePts, nBndEdges;
            int cnt         = 0;
            int nBndRegions = fields[var]->GetBndCondExpansions().size();
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();

            // Extract physical values of the fields at the trace space
            Array<OneD, NekDouble > uplus(nTracePts);
            fields[var]->ExtractTracePhys(ufield, uplus);

            // Impose boundary conditions
            for (i = 0; i < nBndRegions; ++i)
            {
                // Number of boundary edges related to bcRegion 'i'
                nBndEdges = fields[var]->
                GetBndCondExpansions()[i]->GetExpSize();

                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < nBndEdges ; ++e)
                {
                    // Number of points on the expansion
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();

                    // Offset of the boundary expansion
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);

                    // Offset of the trace space related to boundary expansion
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondIDToGlobalTraceID(cnt++));

                    // Dirichlet bcs ==> uflux = gD
                    if (fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vcopy(nBndEdgePts,
                                     &(fields[var]->
                                       GetBndCondExpansions()[i]->
                                       GetPhys())[id1], 1,
                                     &penaltyflux[id2], 1);
                    }
                    // Neumann bcs ==> uflux = u+
                    else if ((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(nBndEdgePts,
                                     &uplus[id2], 1,
                                     &penaltyflux[id2], 1);
                    }
                }
            }
        }

        /**
         * @brief Build the numerical flux for the 2nd order derivatives
         *
         */
        void DiffusionLFR::v_NumFluxforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            boost::ignore_unused(ufield);

            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nvariables = fields.size();
            int nDim       = fields[0]->GetCoordim(0);

            NekDouble C11 = 0.0;
            Array<OneD, NekDouble > Fwd(nTracePts);
            Array<OneD, NekDouble > Bwd(nTracePts);
            Array<OneD, NekDouble > Vn (nTracePts, 0.0);

            Array<OneD, NekDouble > qFwd     (nTracePts);
            Array<OneD, NekDouble > qBwd     (nTracePts);
            Array<OneD, NekDouble > qfluxtemp(nTracePts, 0.0);
            Array<OneD, NekDouble > uterm(nTracePts);

            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
            {
                Vmath::Svtvp(nTracePts, 1.0, m_traceNormals[i], 1,
                             Vn, 1, Vn, 1);
            }

            // Evaulate upwind flux:
            // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
            for (i = 0; i < nvariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
                for (j = 0; j < nDim; ++j)
                {
                    //  Compute Fwd and Bwd value of ufield of jth direction
                    fields[i]->GetFwdBwdTracePhys(qfield[i][j], qFwd, qBwd);

                    // if Vn >= 0, flux = uFwd, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick
                    // qflux = qBwd = q+
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick
                    // qflux = qBwd = q-

                    // else if Vn < 0, flux = uBwd, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick
                    // qflux = qFwd = q-
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick
                    // qflux = qFwd = q+

                    fields[i]->GetTrace()->Upwind(Vn, qBwd, qFwd, qfluxtemp);

                    Vmath::Vmul(nTracePts, m_traceNormals[j], 1, qfluxtemp, 1,
                                qfluxtemp, 1);

                    /*
                    // Generate Stability term = - C11 ( u- - u+ )
                    fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);

                    Vmath::Vsub(nTracePts,
                                Fwd, 1, Bwd, 1,
                                uterm, 1);

                    Vmath::Smul(nTracePts,
                                -1.0 * C11, uterm, 1,
                                uterm, 1);

                    // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
                    Vmath::Vadd(nTracePts,
                                uterm, 1,
                                qfluxtemp, 1,
                                qfluxtemp, 1);
                    */

                    // Imposing weak boundary condition with flux
                    if (fields[0]->GetBndCondExpansions().size())
                    {
                        v_WeakPenaltyforVector(fields, i, j, qfield[i][j],
                                               qfluxtemp, C11);
                    }

                    // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                    // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
                    // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
                    Vmath::Vadd(nTracePts, qfluxtemp, 1, qflux[i], 1,
                                qflux[i], 1);
                }
            }
        }



        /**
         * @brief Imposes appropriate bcs for the 2nd order derivatives
         *
         */
        void DiffusionLFR::v_WeakPenaltyforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const int                                          dir,
            const Array<OneD, const NekDouble>                &qfield,
                  Array<OneD,       NekDouble>                &penaltyflux,
            NekDouble                                          C11)
        {
            boost::ignore_unused(C11);

            int i, e, id1, id2;
            int nBndEdges, nBndEdgePts;
            int nBndRegions = fields[var]->GetBndCondExpansions().size();
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, NekDouble > uterm(nTracePts);
            Array<OneD, NekDouble > qtemp(nTracePts);
            int cnt = 0;

            fields[var]->ExtractTracePhys(qfield, qtemp);

            for (i = 0; i < nBndRegions; ++i)
            {
                nBndEdges = fields[var]->
                GetBndCondExpansions()[i]->GetExpSize();

                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < nBndEdges ; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();

                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);

                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondIDToGlobalTraceID(cnt++));
                    
                    // For Dirichlet boundary condition: 
                    //qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
                    if (fields[var]->GetBndConditions()[i]->
                       GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vmul(nBndEdgePts, &m_traceNormals[dir][id2], 1,
                                    &qtemp[id2], 1, &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: qflux = g_N
                    else if ((fields[var]->GetBndConditions()[i])->
                        GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vmul(nBndEdgePts, &m_traceNormals[dir][id2], 1,
                                    &(fields[var]->GetBndCondExpansions()[i]->
                                      GetPhys())[id1], 1, &penaltyflux[id2], 1);
                    }
                }
            }
        }

        /**
         * @brief Compute the derivative of the corrective flux for 1D problems.
         *
         * @param nConvectiveFields Number of fields.
         * @param fields            Pointer to fields.
         * @param flux              Volumetric flux in the physical space.
         * @param iFlux             Numerical interface flux in physical space.
         * @param derCFlux          Derivative of the corrective flux.
         *
         */
        void DiffusionLFR::v_DerCFlux_1D(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, const NekDouble>                &flux,
            const Array<OneD, const NekDouble>                &iFlux,
                  Array<OneD,       NekDouble>                &derCFlux)
        {
            boost::ignore_unused(nConvectiveFields);

            int n;
            int nLocalSolutionPts, phys_offset, t_offset;

            Array<OneD,       NekDouble> auxArray1, auxArray2;
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;

            LibUtilities::PointsKeyVector ptsKeys;
            LibUtilities::BasisSharedPtr Basis;
            Basis = fields[0]->GetExp(0)->GetBasis(0);

            int nElements    = fields[0]->GetExpSize();
            int nSolutionPts = fields[0]->GetTotPoints();

            vector<bool> negatedFluxNormal =
                std::static_pointer_cast<MultiRegions::DisContField1D>(
                    fields[0])->GetNegatedFluxNormal();

            // Arrays to store the derivatives of the correction flux
            Array<OneD, NekDouble> DCL(nSolutionPts/nElements, 0.0);
            Array<OneD, NekDouble> DCR(nSolutionPts/nElements, 0.0);

            // Arrays to store the intercell numerical flux jumps
            Array<OneD, NekDouble>  JumpL(nElements);
            Array<OneD, NekDouble>  JumpR(nElements);

            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                &elmtToTrace = fields[0]->GetTraceMap()->GetElmtToTrace();

            for (n = 0; n < nElements; ++n)
            {
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                phys_offset = fields[0]->GetPhys_Offset(n);

                Array<OneD, NekDouble> tmparrayX1(nLocalSolutionPts, 0.0);
                NekDouble tmpFluxVertex = 0;
                Vmath::Vcopy(nLocalSolutionPts,
                             &flux[phys_offset], 1,
                             &tmparrayX1[0], 1);

                fields[0]->GetExp(n)->GetVertexPhysVals(0, tmparrayX1,
                                                        tmpFluxVertex);

                t_offset = fields[0]->GetTrace()
                    ->GetPhys_Offset(elmtToTrace[n][0]->GetElmtId());

                if(negatedFluxNormal[2*n])
                {
                    JumpL[n] =  iFlux[t_offset] - tmpFluxVertex;
                }
                else
                {
                    JumpL[n] =  -iFlux[t_offset] - tmpFluxVertex;
                }

                t_offset = fields[0]->GetTrace()
                    ->GetPhys_Offset(elmtToTrace[n][1]->GetElmtId());

                fields[0]->GetExp(n)->GetVertexPhysVals(1, tmparrayX1,
                                                        tmpFluxVertex);
                if(negatedFluxNormal[2*n+1])
                {
                    JumpR[n] =  -iFlux[t_offset] - tmpFluxVertex;
                }
                else
                {
                    JumpR[n] =  iFlux[t_offset] - tmpFluxVertex;
                }
            }

            for (n = 0; n < nElements; ++n)
            {
                ptsKeys = fields[0]->GetExp(n)->GetPointsKeys();
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                phys_offset       = fields[0]->GetPhys_Offset(n);
                jac               = fields[0]->GetExp(n)
                                ->as<LocalRegions::Expansion1D>()->GetGeom1D()
                                ->GetMetricInfo()->GetJac(ptsKeys);

                JumpL[n] = JumpL[n] * jac[0];
                JumpR[n] = JumpR[n] * jac[0];

                // Left jump multiplied by left derivative of C function
                Vmath::Smul(nLocalSolutionPts, JumpL[n],
                            &m_dGL_xi1[n][0], 1, &DCL[0], 1);

                // Right jump multiplied by right derivative of C function
                Vmath::Smul(nLocalSolutionPts, JumpR[n],
                            &m_dGR_xi1[n][0], 1, &DCR[0], 1);

                // Assembling derivative of the correction flux
                Vmath::Vadd(nLocalSolutionPts, &DCL[0], 1, &DCR[0], 1,
                            &derCFlux[phys_offset], 1);
            }
        }

        /**
         * @brief Compute the derivative of the corrective flux wrt a given
         * coordinate for 2D problems.
         *
         * There could be a bug for deformed elements since first derivatives
         * are not exactly the same results of DiffusionLDG as expected
         *
         * @param nConvectiveFields Number of fields.
         * @param fields            Pointer to fields.
         * @param flux              Volumetric flux in the physical space.
         * @param numericalFlux     Numerical interface flux in physical space.
         * @param derCFlux          Derivative of the corrective flux.
         *
         * \todo: Switch on shapes eventually here.
         */
        void DiffusionLFR::v_DerCFlux_2D(
            const int                                         nConvectiveFields,
            const int                                         direction,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, const NekDouble>                &flux,
            const Array<OneD,       NekDouble>                &iFlux,
                  Array<OneD,       NekDouble>                &derCFlux)
        {
            boost::ignore_unused(nConvectiveFields);

            int n, e, i, j, cnt;

            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;

            int nElements = fields[0]->GetExpSize();

            int trace_offset, phys_offset;
            int nLocalSolutionPts;
            int nquad0, nquad1;
            int nEdgePts;

            LibUtilities::PointsKeyVector ptsKeys;
            Array<OneD, NekDouble> auxArray1, auxArray2;
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
            &elmtToTrace = fields[0]->GetTraceMap()->GetElmtToTrace();

            // Loop on the elements
            for (n = 0; n < nElements; ++n)
            {
                // Offset of the element on the global vector
                phys_offset = fields[0]->GetPhys_Offset(n);
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                ptsKeys = fields[0]->GetExp(n)->GetPointsKeys();

                jac  = fields[0]->GetExp(n)->as<LocalRegions::Expansion2D>()
                                ->GetGeom2D()->GetMetricInfo()->GetJac(ptsKeys);

                base = fields[0]->GetExp(n)->GetBase();
                nquad0 = base[0]->GetNumPoints();
                nquad1 = base[1]->GetNumPoints();

                Array<OneD, NekDouble> divCFluxE0(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE1(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE2(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE3(nLocalSolutionPts, 0.0);

                // Loop on the edges
                for (e = 0; e < fields[0]->GetExp(n)->GetNedges(); ++e)
                {
                    // Number of edge points of edge 'e'
                    nEdgePts = fields[0]->GetExp(n)->GetEdgeNumPoints(e);

                    // Array for storing volumetric fluxes on each edge
                    Array<OneD, NekDouble> tmparray(nEdgePts, 0.0);

                    // Offset of the trace space correspondent to edge 'e'
                    trace_offset = fields[0]->GetTrace()->GetPhys_Offset(
                                                elmtToTrace[n][e]->GetElmtId());

                    // Get the normals of edge 'e'
                    //const Array<OneD, const Array<OneD, NekDouble> > &normals =
                    //fields[0]->GetExp(n)->GetEdgeNormal(e);

                    // Extract the edge values of the volumetric fluxes
                    // on edge 'e' and order them accordingly to the order
                    // of the trace space
                    fields[0]->GetExp(n)->GetEdgePhysVals(e, elmtToTrace[n][e],
                                                          flux + phys_offset,
                                                          auxArray1 = tmparray);

                    // Compute the fluxJumps per each edge 'e' and each
                    // flux point
                    Array<OneD, NekDouble> fluxJumps(nEdgePts, 0.0);
                    Vmath::Vsub(nEdgePts, &iFlux[trace_offset], 1,
                                &tmparray[0], 1, &fluxJumps[0], 1);

                    // Check the ordering of the fluxJumps and reverse
                    // it in case of backward definition of edge 'e'
                    if (fields[0]->GetExp(n)->GetEorient(e) ==
                        StdRegions::eBackwards)
                    {
                        Vmath::Reverse(nEdgePts,
                                       &fluxJumps[0], 1,
                                       &fluxJumps[0], 1);
                    }

                    // Deformed elements
                    if (fields[0]->GetExp(n)->as<LocalRegions::Expansion2D>()
                                 ->GetGeom2D()->GetMetricInfo()->GetGtype()
                        == SpatialDomains::eDeformed)
                    {
                        // Extract the Jacobians along edge 'e'
                        Array<OneD, NekDouble> jacEdge(nEdgePts, 0.0);

                        fields[0]->GetExp(n)->GetEdgePhysVals(
                                                    e, elmtToTrace[n][e],
                                                    jac, auxArray1 = jacEdge);

                        // Check the ordering of the fluxJumps and reverse
                        // it in case of backward definition of edge 'e'
                        if (fields[0]->GetExp(n)->GetEorient(e) ==
                            StdRegions::eBackwards)
                        {
                            Vmath::Reverse(nEdgePts,
                                           &jacEdge[0], 1,
                                           &jacEdge[0], 1);
                        }

                        // Multiply the fluxJumps by the edge 'e' Jacobians
                        // to bring the fluxJumps into the standard space
                        for (j = 0; j < nEdgePts; j++)
                        {
                            fluxJumps[j] = fluxJumps[j] * jacEdge[j];
                        }
                    }
                    // Non-deformed elements
                    else
                    {
                        // Multiply the fluxJumps by the edge 'e' Jacobians
                        // to bring the fluxJumps into the standard space
                        Vmath::Smul(nEdgePts, jac[0], fluxJumps, 1,
                                    fluxJumps, 1);
                    }

                    // Multiply jumps by derivatives of the correction functions
                    // All the quntities at this level should be defined into
                    // the standard space
                    switch (e)
                    {
                        case 0:
                            for (i = 0; i < nquad0; ++i)
                            {
                                // Multiply fluxJumps by Q factors
                                //fluxJumps[i] = -(m_Q2D_e0[n][i]) * fluxJumps[i];

                                for (j = 0; j < nquad1; ++j)
                                {
                                    cnt = i + j*nquad0;
                                    divCFluxE0[cnt] = fluxJumps[i] * m_dGL_xi2[n][j];
                                }
                            }
                            break;
                        case 1:
                            for (i = 0; i < nquad1; ++i)
                            {
                                // Multiply fluxJumps by Q factors
                                //fluxJumps[i] = (m_Q2D_e1[n][i]) * fluxJumps[i];

                                for (j = 0; j < nquad0; ++j)
                                {
                                    cnt = (nquad0)*i + j;
                                    divCFluxE1[cnt] = fluxJumps[i] * m_dGR_xi1[n][j];
                                }
                            }
                            break;
                        case 2:
                            for (i = 0; i < nquad0; ++i)
                            {
                                // Multiply fluxJumps by Q factors
                                //fluxJumps[i] = (m_Q2D_e2[n][i]) * fluxJumps[i];

                                for (j = 0; j < nquad1; ++j)
                                {
                                    cnt = j*nquad0 + i;
                                    divCFluxE2[cnt] = fluxJumps[i] * m_dGR_xi2[n][j];
                                }
                            }
                            break;
                        case 3:
                            for (i = 0; i < nquad1; ++i)
                            {
                                // Multiply fluxJumps by Q factors
                                //fluxJumps[i] = -(m_Q2D_e3[n][i]) * fluxJumps[i];
                                for (j = 0; j < nquad0; ++j)
                                {
                                    cnt = j + i*nquad0;
                                    divCFluxE3[cnt] = fluxJumps[i] * m_dGL_xi1[n][j];
                                }
                            }
                            break;

                        default:
                            ASSERTL0(false, "edge value (< 3) is out of range");
                            break;
                    }
                }

                // Summing the component of the flux for calculating the
                // derivatives per each direction
                if (direction == 0)
                {
                    Vmath::Vadd(nLocalSolutionPts, &divCFluxE1[0], 1,
                                &divCFluxE3[0], 1, &derCFlux[phys_offset], 1);
                }
                else if (direction == 1)
                {
                    Vmath::Vadd(nLocalSolutionPts, &divCFluxE0[0], 1,
                                &divCFluxE2[0], 1, &derCFlux[phys_offset], 1);
                }
            }
        }

        /**
         * @brief Compute the divergence of the corrective flux for 2D problems.
         *
         * @param nConvectiveFields   Number of fields.
         * @param fields              Pointer to fields.
         * @param fluxX1              X1 - volumetric flux in physical space.
         * @param fluxX2              X2 - volumetric flux in physical space.
         * @param numericalFlux       Interface flux in physical space.
         * @param divCFlux            Divergence of the corrective flux for
         *                            2D problems.
         *
         * \todo: Switch on shapes eventually here.
         */
        void DiffusionLFR::v_DivCFlux_2D(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, const NekDouble>                &fluxX1,
            const Array<OneD, const NekDouble>                &fluxX2,
            const Array<OneD, const NekDouble>                &numericalFlux,
                  Array<OneD,       NekDouble>                &divCFlux)
        {
            boost::ignore_unused(nConvectiveFields);

            int n, e, i, j, cnt;

            int nElements = fields[0]->GetExpSize();

            int nLocalSolutionPts;
            int nEdgePts;
            int trace_offset;
            int phys_offset;
            int nquad0;
            int nquad1;

            Array<OneD, NekDouble> auxArray1, auxArray2;
            Array<OneD, LibUtilities::BasisSharedPtr> base;

            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
            &elmtToTrace = fields[0]->GetTraceMap()->GetElmtToTrace();

            // Loop on the elements
            for(n = 0; n < nElements; ++n)
            {
                // Offset of the element on the global vector
                phys_offset = fields[0]->GetPhys_Offset(n);
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();

                base   = fields[0]->GetExp(n)->GetBase();
                nquad0 = base[0]->GetNumPoints();
                nquad1 = base[1]->GetNumPoints();

                Array<OneD, NekDouble> divCFluxE0(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE1(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE2(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE3(nLocalSolutionPts, 0.0);

                // Loop on the edges
                for(e = 0; e < fields[0]->GetExp(n)->GetNedges(); ++e)
                {
                    // Number of edge points of edge e
                    nEdgePts = fields[0]->GetExp(n)->GetEdgeNumPoints(e);

                    Array<OneD, NekDouble> tmparrayX1(nEdgePts, 0.0);
                    Array<OneD, NekDouble> tmparrayX2(nEdgePts, 0.0);
                    Array<OneD, NekDouble> fluxN    (nEdgePts, 0.0);
                    Array<OneD, NekDouble> fluxT    (nEdgePts, 0.0);
                    Array<OneD, NekDouble> fluxJumps(nEdgePts, 0.0);

                    // Offset of the trace space correspondent to edge e
                    trace_offset  = fields[0]->GetTrace()->GetPhys_Offset(
                                                elmtToTrace[n][e]->GetElmtId());

                    // Get the normals of edge e
                    const Array<OneD, const Array<OneD, NekDouble> > &normals =
                    fields[0]->GetExp(n)->GetEdgeNormal(e);

                    // Extract the edge values of flux-x on edge e and order
                    // them accordingly to the order of the trace space
                    fields[0]->GetExp(n)->GetEdgePhysVals(
                                                    e, elmtToTrace[n][e],
                                                    fluxX1 + phys_offset,
                                                    auxArray1 = tmparrayX1);

                    // Extract the edge values of flux-y on edge e and order
                    // them accordingly to the order of the trace space
                    fields[0]->GetExp(n)->GetEdgePhysVals(
                                                    e, elmtToTrace[n][e],
                                                    fluxX2 + phys_offset,
                                                    auxArray1 = tmparrayX2);

                    // Multiply the edge components of the flux by the normal
                    for (i = 0; i < nEdgePts; ++i)
                    {
                        fluxN[i] =
                        tmparrayX1[i]*m_traceNormals[0][trace_offset+i] +
                        tmparrayX2[i]*m_traceNormals[1][trace_offset+i];
                    }

                    // Subtract to the Riemann flux the discontinuous flux
                    Vmath::Vsub(nEdgePts,
                                &numericalFlux[trace_offset], 1,
                                &fluxN[0], 1, &fluxJumps[0], 1);

                    // Check the ordering of the jump vectors
                    if (fields[0]->GetExp(n)->GetEorient(e) ==
                        StdRegions::eBackwards)
                    {
                        Vmath::Reverse(nEdgePts,
                                       auxArray1 = fluxJumps, 1,
                                       auxArray2 = fluxJumps, 1);
                    }

                    for (i = 0; i < nEdgePts; ++i)
                    {
                        if (m_traceNormals[0][trace_offset+i] != normals[0][i]
                        || m_traceNormals[1][trace_offset+i] != normals[1][i])
                        {
                            fluxJumps[i] = -fluxJumps[i];
                        }
                    }

                    // Multiply jumps by derivatives of the correction functions
                    switch (e)
                    {
                        case 0:
                            for (i = 0; i < nquad0; ++i)
                            {
                                // Multiply fluxJumps by Q factors
                                fluxJumps[i] = -(m_Q2D_e0[n][i]) * fluxJumps[i];

                                for (j = 0; j < nquad1; ++j)
                                {
                                    cnt = i + j*nquad0;
                                    divCFluxE0[cnt] = fluxJumps[i] * m_dGL_xi2[n][j];
                                }
                            }
                            break;
                        case 1:
                            for (i = 0; i < nquad1; ++i)
                            {
                                // Multiply fluxJumps by Q factors
                                fluxJumps[i] = (m_Q2D_e1[n][i]) * fluxJumps[i];

                                for (j = 0; j < nquad0; ++j)
                                {
                                    cnt = (nquad0)*i + j;
                                    divCFluxE1[cnt] = fluxJumps[i] * m_dGR_xi1[n][j];
                                }
                            }
                            break;
                        case 2:
                            for (i = 0; i < nquad0; ++i)
                            {
                                // Multiply fluxJumps by Q factors
                                fluxJumps[i] = (m_Q2D_e2[n][i]) * fluxJumps[i];

                                for (j = 0; j < nquad1; ++j)
                                {
                                    cnt = j*nquad0 + i;
                                    divCFluxE2[cnt] = fluxJumps[i] * m_dGR_xi2[n][j];
                                }
                            }
                            break;
                        case 3:
                            for (i = 0; i < nquad1; ++i)
                            {
                                // Multiply fluxJumps by Q factors
                                fluxJumps[i] = -(m_Q2D_e3[n][i]) * fluxJumps[i];
                                for (j = 0; j < nquad0; ++j)
                                {
                                    cnt = j + i*nquad0;
                                    divCFluxE3[cnt] = fluxJumps[i] * m_dGL_xi1[n][j];
                                }
                            }
                            break;

                        default:
                            ASSERTL0(false,"edge value (< 3) is out of range");
                            break;
                    }
                }

                // Sum each edge contribution
                Vmath::Vadd(nLocalSolutionPts, &divCFluxE0[0], 1,
                            &divCFluxE1[0], 1, &divCFlux[phys_offset], 1);

                Vmath::Vadd(nLocalSolutionPts, &divCFlux[phys_offset], 1,
                            &divCFluxE2[0], 1, &divCFlux[phys_offset], 1);

                Vmath::Vadd(nLocalSolutionPts, &divCFlux[phys_offset], 1,
                            &divCFluxE3[0], 1, &divCFlux[phys_offset], 1);
            }
        }

        /**
         * @brief Compute the divergence of the corrective flux for 2D problems
         *        where POINTSTYPE="GaussGaussLegendre"
         *
         * @param nConvectiveFields   Number of fields.
         * @param fields              Pointer to fields.
         * @param fluxX1              X1-volumetric flux in physical space.
         * @param fluxX2              X2-volumetric flux in physical space.
         * @param numericalFlux       Interface flux in physical space.
         * @param divCFlux            Divergence of the corrective flux.
         *
         * \todo: Switch on shapes eventually here.
         */

        void DiffusionLFR::v_DivCFlux_2D_Gauss(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, const NekDouble> &fluxX1,
            const Array<OneD, const NekDouble> &fluxX2,
            const Array<OneD, const NekDouble> &numericalFlux,
                  Array<OneD,       NekDouble> &divCFlux)
        {
            boost::ignore_unused(nConvectiveFields);

            int n, e, i, j, cnt;

            int nElements   = fields[0]->GetExpSize();
            int nLocalSolutionPts;
            int nEdgePts;
            int trace_offset;
            int phys_offset;
            int nquad0;
            int nquad1;

            Array<OneD, NekDouble> auxArray1, auxArray2;
            Array<OneD, LibUtilities::BasisSharedPtr> base;

            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
            &elmtToTrace = fields[0]->GetTraceMap()->GetElmtToTrace();

            // Loop on the elements
            for(n = 0; n < nElements; ++n)
            {
                // Offset of the element on the global vector
                phys_offset = fields[0]->GetPhys_Offset(n);
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();

                base = fields[0]->GetExp(n)->GetBase();
                nquad0 = base[0]->GetNumPoints();
                nquad1 = base[1]->GetNumPoints();

                Array<OneD, NekDouble> divCFluxE0(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE1(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE2(nLocalSolutionPts, 0.0);
                Array<OneD, NekDouble> divCFluxE3(nLocalSolutionPts, 0.0);

                // Loop on the edges
                for(e = 0; e < fields[0]->GetExp(n)->GetNedges(); ++e)
                {
                    // Number of edge points of edge e
                    nEdgePts = fields[0]->GetExp(n)->GetEdgeNumPoints(e);

                    Array<OneD, NekDouble> fluxN     (nEdgePts, 0.0);
                    Array<OneD, NekDouble> fluxT     (nEdgePts, 0.0);
                    Array<OneD, NekDouble> fluxN_R   (nEdgePts, 0.0);
                    Array<OneD, NekDouble> fluxN_D   (nEdgePts, 0.0);
                    Array<OneD, NekDouble> fluxJumps (nEdgePts, 0.0);

                    // Offset of the trace space correspondent to edge e
                    trace_offset = fields[0]->GetTrace()->GetPhys_Offset(
                                                elmtToTrace[n][e]->GetElmtId());

                    // Get the normals of edge e
                    const Array<OneD, const Array<OneD, NekDouble> > &normals =
                    fields[0]->GetExp(n)->GetEdgeNormal(e);

                    // Extract the trasformed normal flux at each edge
                    switch (e)
                    {
                        case 0:
                            // Extract the edge values of transformed flux-y on
                            // edge e and order them accordingly to the order of
                            // the trace space
                            fields[0]->GetExp(n)->GetEdgePhysVals(
                                                        e, elmtToTrace[n][e],
                                                        fluxX2 + phys_offset,
                                                        auxArray1 = fluxN_D);

                            Vmath::Neg   (nEdgePts, fluxN_D, 1);

                            // Extract the physical Rieamann flux at each edge
                            Vmath::Vcopy(nEdgePts,
                                         &numericalFlux[trace_offset], 1,
                                         &fluxN[0], 1);

                            // Check the ordering of vectors
                            if (fields[0]->GetExp(n)->GetEorient(e) ==
                                StdRegions::eBackwards)
                            {
                                Vmath::Reverse(nEdgePts,
                                               auxArray1 = fluxN, 1,
                                               auxArray2 = fluxN, 1);

                                Vmath::Reverse(nEdgePts,
                                               auxArray1 = fluxN_D, 1,
                                               auxArray2 = fluxN_D, 1);
                            }

                            // Transform Riemann Fluxes in the standard element
                            for (i = 0; i < nquad0; ++i)
                            {
                                // Multiply Riemann Flux by Q factors
                                fluxN_R[i] = (m_Q2D_e0[n][i]) * fluxN[i];
                            }

                            for (i = 0; i < nEdgePts; ++i)
                            {
                                if (m_traceNormals[0][trace_offset+i] !=
                                    normals[0][i] ||
                                    m_traceNormals[1][trace_offset+i] !=
                                    normals[1][i])
                                {
                                    fluxN_R[i] = -fluxN_R[i];
                                }
                            }

                            // Subtract to the Riemann flux the discontinuous
                            // flux
                            Vmath::Vsub(nEdgePts,
                                        &fluxN_R[0], 1,
                                        &fluxN_D[0], 1, &fluxJumps[0], 1);

                            // Multiplicate the flux jumps for the correction
                            // function
                            for (i = 0; i < nquad0; ++i)
                            {
                                for (j = 0; j < nquad1; ++j)
                                {
                                    cnt = i + j*nquad0;
                                    divCFluxE0[cnt] =
                                    -fluxJumps[i] * m_dGL_xi2[n][j];
                                }
                            }
                            break;
                        case 1:
                            // Extract the edge values of transformed flux-x on
                            // edge e and order them accordingly to the order of
                            // the trace space
                            fields[0]->GetExp(n)->GetEdgePhysVals(
                                                        e, elmtToTrace[n][e],
                                                        fluxX1 + phys_offset,
                                                        auxArray1 = fluxN_D);

                            // Extract the physical Rieamann flux at each edge
                            Vmath::Vcopy(nEdgePts,
                                         &numericalFlux[trace_offset], 1,
                                         &fluxN[0], 1);

                            // Check the ordering of vectors
                            if (fields[0]->GetExp(n)->GetEorient(e) ==
                                StdRegions::eBackwards)
                            {
                                Vmath::Reverse(nEdgePts,
                                               auxArray1 = fluxN, 1,
                                               auxArray2 = fluxN, 1);

                                Vmath::Reverse(nEdgePts,
                                               auxArray1 = fluxN_D, 1,
                                               auxArray2 = fluxN_D, 1);
                            }

                            // Transform Riemann Fluxes in the standard element
                            for (i = 0; i < nquad1; ++i)
                            {
                                // Multiply Riemann Flux by Q factors
                                fluxN_R[i] = (m_Q2D_e1[n][i]) * fluxN[i];
                            }

                            for (i = 0; i < nEdgePts; ++i)
                            {
                                if (m_traceNormals[0][trace_offset+i] !=
                                    normals[0][i] ||
                                    m_traceNormals[1][trace_offset+i] !=
                                    normals[1][i])
                                {
                                    fluxN_R[i] = -fluxN_R[i];
                                }
                            }

                            // Subtract to the Riemann flux the discontinuous
                            // flux
                            Vmath::Vsub(nEdgePts,
                                        &fluxN_R[0], 1,
                                        &fluxN_D[0], 1, &fluxJumps[0], 1);

                            // Multiplicate the flux jumps for the correction
                            // function
                            for (i = 0; i < nquad1; ++i)
                            {
                                for (j = 0; j < nquad0; ++j)
                                {
                                    cnt = (nquad0)*i + j;
                                    divCFluxE1[cnt] =
                                    fluxJumps[i] * m_dGR_xi1[n][j];
                                }
                            }
                            break;
                        case 2:

                            // Extract the edge values of transformed flux-y on
                            // edge e and order them accordingly to the order of
                            // the trace space

                            fields[0]->GetExp(n)->GetEdgePhysVals(
                                                        e, elmtToTrace[n][e],
                                                        fluxX2 + phys_offset,
                                                        auxArray1 = fluxN_D);

                            // Extract the physical Rieamann flux at each edge
                            Vmath::Vcopy(nEdgePts,
                                         &numericalFlux[trace_offset], 1,
                                         &fluxN[0], 1);

                            // Check the ordering of vectors
                            if (fields[0]->GetExp(n)->GetEorient(e) ==
                                StdRegions::eBackwards)
                            {
                                Vmath::Reverse(nEdgePts,
                                               auxArray1 = fluxN, 1,
                                               auxArray2 = fluxN, 1);

                                Vmath::Reverse(nEdgePts,
                                               auxArray1 = fluxN_D, 1,
                                               auxArray2 = fluxN_D, 1);
                            }

                            // Transform Riemann Fluxes in the standard element
                            for (i = 0; i < nquad0; ++i)
                            {
                                // Multiply Riemann Flux by Q factors
                                fluxN_R[i] = (m_Q2D_e2[n][i]) * fluxN[i];
                            }

                            for (i = 0; i < nEdgePts; ++i)
                            {
                                if (m_traceNormals[0][trace_offset+i] !=
                                    normals[0][i] ||
                                    m_traceNormals[1][trace_offset+i] !=
                                    normals[1][i])
                                {
                                    fluxN_R[i] = -fluxN_R[i];
                                }
                            }

                            // Subtract to the Riemann flux the discontinuous
                            // flux

                            Vmath::Vsub(nEdgePts,
                                        &fluxN_R[0], 1,
                                        &fluxN_D[0], 1, &fluxJumps[0], 1);

                            // Multiplicate the flux jumps for the correction
                            // function
                            for (i = 0; i < nquad0; ++i)
                            {
                                for (j = 0; j < nquad1; ++j)
                                {
                                    cnt = j*nquad0 + i;
                                    divCFluxE2[cnt] =
                                    fluxJumps[i] * m_dGR_xi2[n][j];
                                }
                            }
                            break;
                        case 3:
                            // Extract the edge values of transformed flux-x on
                            // edge e and order them accordingly to the order of
                            // the trace space

                            fields[0]->GetExp(n)->GetEdgePhysVals(
                                                        e, elmtToTrace[n][e],
                                                        fluxX1 + phys_offset,
                                                        auxArray1 = fluxN_D);
                            Vmath::Neg   (nEdgePts, fluxN_D, 1);

                            // Extract the physical Rieamann flux at each edge
                            Vmath::Vcopy(nEdgePts,
                                         &numericalFlux[trace_offset], 1,
                                         &fluxN[0], 1);

                            // Check the ordering of vectors
                            if (fields[0]->GetExp(n)->GetEorient(e) ==
                                StdRegions::eBackwards)
                            {
                                Vmath::Reverse(nEdgePts,
                                               auxArray1 = fluxN, 1,
                                               auxArray2 = fluxN, 1);

                                Vmath::Reverse(nEdgePts,
                                               auxArray1 = fluxN_D, 1,
                                               auxArray2 = fluxN_D, 1);
                            }

                            // Transform Riemann Fluxes in the standard element
                            for (i = 0; i < nquad1; ++i)
                            {
                                // Multiply Riemann Flux by Q factors
                                fluxN_R[i] = (m_Q2D_e3[n][i]) * fluxN[i];
                            }

                            for (i = 0; i < nEdgePts; ++i)
                            {
                                if (m_traceNormals[0][trace_offset+i] !=
                                    normals[0][i] ||
                                    m_traceNormals[1][trace_offset+i] !=
                                    normals[1][i])
                                {
                                    fluxN_R[i] = -fluxN_R[i];
                                }
                            }

                            // Subtract to the Riemann flux the discontinuous
                            // flux

                            Vmath::Vsub(nEdgePts,
                                        &fluxN_R[0], 1,
                                        &fluxN_D[0], 1, &fluxJumps[0], 1);

                            // Multiplicate the flux jumps for the correction
                            // function
                            for (i = 0; i < nquad1; ++i)
                            {
                                for (j = 0; j < nquad0; ++j)
                                {
                                    cnt = j + i*nquad0;
                                    divCFluxE3[cnt] =
                                    -fluxJumps[i] * m_dGL_xi1[n][j];
                                }
                            }
                            break;
                        default:
                            ASSERTL0(false,"edge value (< 3) is out of range");
                            break;
                    }
                }


                // Sum each edge contribution
                Vmath::Vadd(nLocalSolutionPts, &divCFluxE0[0], 1,
                            &divCFluxE1[0], 1, &divCFlux[phys_offset], 1);

                Vmath::Vadd(nLocalSolutionPts, &divCFlux[phys_offset], 1,
                            &divCFluxE2[0], 1, &divCFlux[phys_offset], 1);

                Vmath::Vadd(nLocalSolutionPts, &divCFlux[phys_offset], 1,
                            &divCFluxE3[0], 1, &divCFlux[phys_offset], 1);
            }
        }

    }// close namespace SolverUtils
}// close namespace nektar++
