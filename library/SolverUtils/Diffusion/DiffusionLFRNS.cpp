///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFRNS.cpp
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
// Description: LFRNS diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <SolverUtils/Diffusion/DiffusionLFRNS.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <boost/core/ignore_unused.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLFRNS::type[] = {
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRDGNS", DiffusionLFRNS::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRSDNS", DiffusionLFRNS::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRHUNS", DiffusionLFRNS::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRcminNS", DiffusionLFRNS::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                                        "LFRcinfNS", DiffusionLFRNS::create)};

        /**
         * @brief DiffusionLFRNS uses the Flux Reconstruction (FR) approach to
         * compute the diffusion term. The implementation is only for segments,
         * quadrilaterals and hexahedra at the moment.
         *
         * \todo Extension to triangles, tetrahedra and other shapes.
         * (Long term objective)
         */
        DiffusionLFRNS::DiffusionLFRNS(std::string diffType):m_diffType(diffType)
        {
        }

        /**
         * @brief Initiliase DiffusionLFRNS objects and store them before
         * starting the time-stepping.
         *
         * This routine calls the virtual functions #v_SetupMetrics and
         * #v_SetupCFunctions to initialise the objects needed
         * by DiffusionLFRNS.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void DiffusionLFRNS::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            m_session = pSession;
            m_session->LoadParameter ("Gamma",         m_gamma, 1.4);
            m_session->LoadParameter ("GasConstant",   m_gasConstant, 287.058);
            m_session->LoadParameter ("Twall",         m_Twall, 300.15);
            m_session->LoadSolverInfo("ViscosityType", m_ViscosityType,
                                      "Constant");
            m_session->LoadParameter ("mu",            m_mu, 1.78e-05);
            m_session->LoadParameter ("thermalConductivity",
                                      m_thermalConductivity, 0.0257);
            m_session->LoadParameter ("rhoInf",        m_rhoInf, 1.225);
            m_session->LoadParameter ("pInf",          m_pInf, 101325);
            v_SetupMetrics(pSession, pFields);
            v_SetupCFunctions(pSession, pFields);

            // Initialising arrays
            int i, j;
            int nConvectiveFields = pFields.size();
            int nScalars     = nConvectiveFields - 1;
            int nDim         = pFields[0]->GetCoordim(0);
            int nSolutionPts = pFields[0]->GetTotPoints();
            int nTracePts    = pFields[0]->GetTrace()->GetTotPoints();

            m_spaceDim = nDim;
            if (pSession->DefinesSolverInfo("HOMOGENEOUS"))
            {
                m_spaceDim = 3;
            }

            m_diffDim = m_spaceDim - nDim;

            m_traceVel = Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);

            for (i = 0; i < m_spaceDim; ++i)
            {
                m_traceVel[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
            }

            m_IF1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nScalars);
            m_DU1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nScalars);
            m_DFC1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nScalars);
            m_tmp1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nScalars);
            m_tmp2 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nScalars);
            m_BD1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nScalars);


            m_DFC2 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_DD1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                            nConvectiveFields);
            m_viscFlux = Array<OneD, Array<OneD, NekDouble> > (
                                                            nConvectiveFields);
            m_divFD = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
            m_divFC = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);

            m_D1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                                (m_spaceDim);
            m_viscTensor = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                                (m_spaceDim);
            for (i = 0; i < nScalars; ++i)
            {
                m_IF1[i]  = Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                m_DU1[i]  = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_DFC1[i] = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_tmp1[i] = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_tmp2[i] = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_BD1[i]  = Array<OneD, Array<OneD, NekDouble> >(nDim);

                for (j = 0; j < nDim; ++j)
                {
                    m_IF1[i][j]  = Array<OneD, NekDouble>(nTracePts, 0.0);
                    m_DU1[i][j]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_DFC1[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_tmp1[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_tmp2[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_BD1[i][j]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                }
            }

            for (i = 0; i < nConvectiveFields; ++i)
            {
                m_DFC2[i]     = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_DD1[i]      = Array<OneD, Array<OneD, NekDouble> >(nDim);
                m_viscFlux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_divFD[i]    = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                m_divFC[i]    = Array<OneD, NekDouble>(nSolutionPts, 0.0);

                for (j = 0; j < nDim; ++j)
                {
                    m_DFC2[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                    m_DD1[i][j]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                }
            }

            for (i = 0; i < m_spaceDim; ++i)
            {
                m_D1[i] = Array<OneD, Array<OneD, NekDouble> >(nScalars);

                for (j = 0; j < nScalars; ++j)
                {
                    m_D1[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                }
            }

            for (i = 0; i < m_spaceDim; ++i)
            {
                m_viscTensor[i] = Array<OneD, Array<OneD, NekDouble> >(nScalars+1);

                for (j = 0; j < nScalars+1; ++j)
                {
                    m_viscTensor[i][j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
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
        void DiffusionLFRNS::v_SetupMetrics(
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
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
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
         *
         * The values of the derivatives of the correction function are then
         * stored into global variables and reused to compute the corrective
         * fluxes.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void DiffusionLFRNS::v_SetupCFunctions(
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
                        base      = pFields[0]->GetExp(n)->GetBase();
                        nquad0    = base[0]->GetNumPoints();
                        nmodes0   = base[0]->GetNumModes();
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
                        if (m_diffType == "LFRDGNS")
                        {
                            c0 = 0.0;
                        }
                        else if (m_diffType == "LFRSDNS")
                        {
                            c0 = 2.0 * p0 / ((2.0 * p0 + 1.0) * (p0 + 1.0)
                                             * (ap0 * boost::math::tgamma(p0 + 1))
                                             * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_diffType == "LFRHUNS")
                        {
                            c0 = 2.0 * (p0 + 1.0) / ((2.0 * p0 + 1.0) * p0
                                                     * (ap0 * boost::math::tgamma(p0 + 1))
                                                     * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_diffType == "LFRcminNS")
                        {
                            c0 = -2.0 / ((2.0 * p0 + 1.0)
                                         * (ap0 * boost::math::tgamma(p0 + 1))
                                         * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_diffType == "LFRcinfNS")
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
                        if (m_diffType == "LFRDGNS")
                        {
                            c0 = 0.0;
                            c1 = 0.0;
                        }
                        else if (m_diffType == "LFRSDNS")
                        {
                            c0 = 2.0 * p0 / ((2.0 * p0 + 1.0) * (p0 + 1.0)
                                             * (ap0 * boost::math::tgamma(p0 + 1))
                                             * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = 2.0 * p1 / ((2.0 * p1 + 1.0) * (p1 + 1.0)
                                             * (ap1 * boost::math::tgamma(p1 + 1))
                                             * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_diffType == "LFRHUNS")
                        {
                            c0 = 2.0 * (p0 + 1.0) / ((2.0 * p0 + 1.0) * p0
                                                     * (ap0 * boost::math::tgamma(p0 + 1))
                                                     * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = 2.0 * (p1 + 1.0) / ((2.0 * p1 + 1.0) * p1
                                                     * (ap1 * boost::math::tgamma(p1 + 1))
                                                     * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_diffType == "LFRcminNS")
                        {
                            c0 = -2.0 / ((2.0 * p0 + 1.0)
                                         * (ap0 * boost::math::tgamma(p0 + 1))
                                         * (ap0 * boost::math::tgamma(p0 + 1)));

                            c1 = -2.0 / ((2.0 * p1 + 1.0)
                                         * (ap1 * boost::math::tgamma(p1 + 1))
                                         * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_diffType == "LFRcinfNS")
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
                        if (m_diffType == "LFRDGNS")
                        {
                            c0 = 0.0;
                            c1 = 0.0;
                            c2 = 0.0;
                        }
                        else if (m_diffType == "LFRSDNS")
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
                        else if (m_diffType == "LFRHUNS")
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
                        else if (m_diffType == "LFRcminNS")
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
                        else if (m_diffType == "LFRcinfNS")
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
         * @brief Calculate FR Diffusion for the Navier-Stokes (NS) equations
         * using an LDG interface flux.
         *
         * The equations that need a diffusion operator are those related
         * with the velocities and with the energy.
         *
         */
        void DiffusionLFRNS::v_Diffuse(
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

            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            Array<OneD,       NekDouble> auxArray1, auxArray2;

            Array<OneD, LibUtilities::BasisSharedPtr> Basis;
            Basis = fields[0]->GetExp(0)->GetBase();

            int nElements    = fields[0]->GetExpSize();
            int nDim         = fields[0]->GetCoordim(0);
            int nScalars     = inarray.size();
            int nSolutionPts = fields[0]->GetTotPoints();
            int nCoeffs      = fields[0]->GetNcoeffs();

            Array<OneD, Array<OneD, NekDouble> > outarrayCoeff(nConvectiveFields);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                outarrayCoeff[i]  = Array<OneD, NekDouble>(nCoeffs);
            }

            // Compute interface numerical fluxes for inarray in physical space
            v_NumericalFluxO1(fields, inarray, m_IF1);

            switch(nDim)
            {
                // 1D problems
                case 1:
                {
                    for (i = 0; i < nScalars; ++i)
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

                    // Computing the viscous tensor
                    m_fluxVectorNS(inarray, m_D1, m_viscTensor);

                    // Compute u from q_{\eta} and q_{\xi}
                    // Obtain numerical fluxes
                    v_NumericalFluxO2(fields, inarray, m_viscTensor,
                                      m_viscFlux);

                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        // Computing the physical second-order discountinuous
                        // derivative
                        for (n = 0; n < nElements; n++)
                        {
                            phys_offset = fields[0]->GetPhys_Offset(n);

                            fields[i]->GetExp(n)->PhysDeriv(0,
                                auxArray1 = m_viscTensor[0][i] + phys_offset,
                                auxArray2 = m_DD1[i][0] + phys_offset);
                        }

                        // Computing the standard second-order correction
                        // derivative
                        v_DerCFlux_1D(nConvectiveFields, fields,
                                      m_viscTensor[0][i], m_viscFlux[i],
                                      m_DFC2[i][0]);

                        // Back to the physical space using global operations
                        Vmath::Vdiv(nSolutionPts, &m_DFC2[i][0][0], 1,
                                    &m_jac[0], 1, &m_DFC2[i][0][0], 1);
                        Vmath::Vdiv(nSolutionPts, &m_DFC2[i][0][0], 1,
                                    &m_jac[0], 1, &m_DFC2[i][0][0], 1);

                        // Adding the total divergence to outarray (RHS)
                        Vmath::Vadd(nSolutionPts, &m_DFC2[i][0][0], 1,
                                    &m_DD1[i][0][0], 1, &outarray[i][0], 1);

                        // Primitive Dealiasing 1D
                        if(!(Basis[0]->Collocation()))
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
                    for(i = 0; i < nScalars; ++i)
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
                                        &m_D1[j][i][0], 1);
                        }
                    }// Close loop on nScalars

                    // For 3D Homogeneous 1D only take derivatives in 3rd direction
                    if (m_diffDim == 1)
                    {
                        for (i = 0; i < nScalars; ++i)
                        {
                            m_D1[2][i] = m_homoDerivs[i];
                        }

                    }

                    // Computing the viscous tensor
                    m_fluxVectorNS(inarray, m_D1, m_viscTensor);

                    // Compute u from q_{\eta} and q_{\xi}
                    // Obtain numerical fluxes
                    v_NumericalFluxO2(fields, inarray, m_viscTensor,
                                      m_viscFlux);

                    // Computing the standard second-order discontinuous
                    // derivatives
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        // Temporary vectors
                        Array<OneD, NekDouble> f_hat(nSolutionPts, 0.0);
                        Array<OneD, NekDouble> g_hat(nSolutionPts, 0.0);

                        for (j = 0; j < nSolutionPts; j++)
                        {
                            f_hat[j] = (m_viscTensor[0][i][j] * m_gmat[0][j] +
                            m_viscTensor[1][i][j] * m_gmat[2][j])*m_jac[j];

                            g_hat[j] = (m_viscTensor[0][i][j] * m_gmat[1][j] +
                            m_viscTensor[1][i][j] * m_gmat[3][j])*m_jac[j];
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
                        Vmath::Vadd(nSolutionPts, &m_DD1[i][0][0], 1,
                                    &m_DD1[i][1][0], 1, &m_divFD[i][0], 1);

                        // Divergence of the standard correction flux
                        if (Basis[0]->GetPointsType() ==
                            LibUtilities::eGaussGaussLegendre &&
                            Basis[1]->GetPointsType() ==
                            LibUtilities::eGaussGaussLegendre)
                        {
                            v_DivCFlux_2D_Gauss(nConvectiveFields, fields,
                                                f_hat, g_hat,
                                                m_viscFlux[i], m_divFC[i]);
                        }
                        else
                        {
                            v_DivCFlux_2D(nConvectiveFields, fields,
                                          m_viscTensor[0][i], m_viscTensor[1][i],
                                          m_viscFlux[i], m_divFC[i]);

                        }


                        // Divergence of the standard final flux
                        Vmath::Vadd(nSolutionPts, &m_divFD[i][0], 1,
                                    &m_divFC[i][0], 1, &outarray[i][0], 1);

                        // Dividing by the jacobian to get back into
                        // physical space
                        Vmath::Vdiv(nSolutionPts, &outarray[i][0], 1,
                                    &m_jac[0], 1, &outarray[i][0], 1);

                        // Primitive Dealiasing 2D
                        if(!(Basis[0]->Collocation()))
                        {
                            fields[i]->FwdTrans(outarray[i], outarrayCoeff[i]);
                            fields[i]->BwdTrans(outarrayCoeff[i], outarray[i]);
                        }
                    }
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
        void DiffusionLFRNS::v_NumericalFluxO1(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                            &numericalFluxO1)
        {
            int i, j;
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            int nScalars  = inarray.size();
            int nDim      = fields[0]->GetCoordim(0);

            Array<OneD, NekDouble > Vn      (nTracePts, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePts, 0.0);

            // Get the normal velocity Vn
            for (i = 0; i < nDim; ++i)
            {
                fields[0]->ExtractTracePhys(inarray[i], m_traceVel[i]);
                Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1,
                             m_traceVel[i], 1, Vn, 1, Vn, 1);
            }

            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nScalars);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nScalars);
            Array<OneD, Array<OneD, NekDouble> > numflux(nScalars);

            for (i = 0; i < nScalars; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts);
                numflux[i] = Array<OneD, NekDouble>(nTracePts);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                fields[0]->GetTrace()->Upwind(Vn, Fwd[i], Bwd[i], numflux[i]);
            }

            // Modify the values in case of boundary interfaces
            if (fields[0]->GetBndCondExpansions().size())
            {
                v_WeakPenaltyO1(fields, inarray, numflux);
            }

            // Splitting the numerical flux into the dimensions
            for (j = 0; j < nDim; ++j)
            {
                for (i = 0; i < nScalars; ++i)
                {
                    Vmath::Vcopy(nTracePts,numflux[i], 1,
                                 numericalFluxO1[i][j], 1);
                }
            }
        }

        /**
         * @brief Imposes appropriate bcs for the 1st order derivatives
         *
         */
        void DiffusionLFRNS::v_WeakPenaltyO1(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &penaltyfluxO1)
        {
            int cnt;
            int i, j, e;
            int id1, id2;

            int nBndEdgePts, nBndEdges, nBndRegions;

            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            int nScalars  = inarray.size();

            Array<OneD, NekDouble> tmp1(nTracePts, 0.0);
            Array<OneD, NekDouble> tmp2(nTracePts, 0.0);
            Array<OneD, NekDouble> Tw(nTracePts, m_Twall);

            Array< OneD, Array<OneD, NekDouble > > scalarVariables(nScalars);
            Array< OneD, Array<OneD, NekDouble > > uplus(nScalars);

            // Extract internal values of the scalar variables for Neumann bcs
            for (i = 0; i < nScalars; ++i)
            {
                scalarVariables[i] = Array<OneD, NekDouble>(nTracePts, 0.0);

                uplus[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                fields[i]->ExtractTracePhys(inarray[i], uplus[i]);
            }

            // Compute boundary conditions for velocities
            for (i = 0; i < nScalars-1; ++i)
            {
                // Note that cnt has to loop on nBndRegions and nBndEdges
                // and has to be reset to zero per each equation
                cnt = 0;
                nBndRegions = fields[i+1]->
                GetBndCondExpansions().size();
                for (j = 0; j < nBndRegions; ++j)
                {
                    if (fields[i+1]->GetBndConditions()[j]->
                        GetBoundaryConditionType() == SpatialDomains::ePeriodic)
                    {
                        continue;
                    }

                    nBndEdges = fields[i+1]->
                    GetBndCondExpansions()[j]->GetExpSize();
                    for (e = 0; e < nBndEdges; ++e)
                    {
                        nBndEdgePts = fields[i+1]->
                        GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                        id1 = fields[i+1]->
                        GetBndCondExpansions()[j]->GetPhys_Offset(e);

                        id2 = fields[0]->GetTrace()->
                        GetPhys_Offset(fields[0]->GetTraceMap()->
                                       GetBndCondIDToGlobalTraceID(cnt++));

                        // Reinforcing bcs for velocity in case of Wall bcs
                        if (boost::iequals(fields[i]->GetBndConditions()[j]->
                            GetUserDefined(),"WallViscous") ||
                            boost::iequals(fields[i]->GetBndConditions()[j]->
                            GetUserDefined(),"WallAdiabatic"))
                        {
                            Vmath::Zero(nBndEdgePts,
                                        &scalarVariables[i][id2], 1);
                        }

                        // Imposing velocity bcs if not Wall
                        else if (fields[i]->GetBndConditions()[j]->
                                 GetBoundaryConditionType() ==
                                 SpatialDomains::eDirichlet)
                        {
                            Vmath::Vdiv(nBndEdgePts,
                                        &(fields[i+1]->
                                          GetBndCondExpansions()[j]->
                                          UpdatePhys())[id1], 1,
                                        &(fields[0]->
                                          GetBndCondExpansions()[j]->
                                          UpdatePhys())[id1], 1,
                                        &scalarVariables[i][id2], 1);
                        }

                        // For Dirichlet boundary condition: uflux = u_bcs
                        if (fields[i]->GetBndConditions()[j]->
                            GetBoundaryConditionType() ==
                            SpatialDomains::eDirichlet)
                        {
                            Vmath::Vcopy(nBndEdgePts,
                                         &scalarVariables[i][id2], 1,
                                         &penaltyfluxO1[i][id2], 1);
                        }

                        // For Neumann boundary condition: uflux = u_+
                        else if ((fields[i]->GetBndConditions()[j])->
                                 GetBoundaryConditionType() ==
                                 SpatialDomains::eNeumann)
                        {
                            Vmath::Vcopy(nBndEdgePts,
                                         &uplus[i][id2], 1,
                                         &penaltyfluxO1[i][id2], 1);
                        }

                        // Building kinetic energy to be used for T bcs
                        Vmath::Vmul(nBndEdgePts,
                                    &scalarVariables[i][id2], 1,
                                    &scalarVariables[i][id2], 1,
                                    &tmp1[id2], 1);

                        Vmath::Smul(nBndEdgePts, 0.5,
                                    &tmp1[id2], 1,
                                    &tmp1[id2], 1);

                        Vmath::Vadd(nBndEdgePts,
                                    &tmp2[id2], 1,
                                    &tmp1[id2], 1,
                                    &tmp2[id2], 1);
                    }
                }
            }

            // Compute boundary conditions  for temperature
            cnt = 0;
            nBndRegions = fields[nScalars]->
            GetBndCondExpansions().size();
            for (j = 0; j < nBndRegions; ++j)
            {
                nBndEdges = fields[nScalars]->
                GetBndCondExpansions()[j]->GetExpSize();

                if (fields[nScalars]->GetBndConditions()[j]->
                    GetBoundaryConditionType() == SpatialDomains::ePeriodic)
                {
                    continue;
                }

                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = fields[nScalars]->
                    GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                    id1 = fields[nScalars]->
                    GetBndCondExpansions()[j]->GetPhys_Offset(e);

                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondIDToGlobalTraceID(cnt++));
                    
                    // Imposing Temperature Twall at the wall 
                    if (boost::iequals(fields[i]->GetBndConditions()[j]->
                                       GetUserDefined(),"WallViscous"))
                    {
                        Vmath::Vcopy(nBndEdgePts,
                                     &Tw[0], 1,
                                     &scalarVariables[nScalars-1][id2], 1);
                    }
                    // Imposing Temperature through condition on the Energy
                    // for no wall boundaries (e.g. farfield)
                    else if (fields[i]->GetBndConditions()[j]->
                             GetBoundaryConditionType() ==
                             SpatialDomains::eDirichlet)
                    {
                        // Divide E by rho
                        Vmath::Vdiv(nBndEdgePts,
                                    &(fields[nScalars]->
                                      GetBndCondExpansions()[j]->
                                      GetPhys())[id1], 1,
                                    &(fields[0]->
                                      GetBndCondExpansions()[j]->
                                      GetPhys())[id1], 1,
                                    &scalarVariables[nScalars-1][id2], 1);

                        // Subtract kinetic energy to E/rho
                        Vmath::Vsub(nBndEdgePts,
                                    &scalarVariables[nScalars-1][id2], 1,
                                    &tmp2[id2], 1,
                                    &scalarVariables[nScalars-1][id2], 1);

                        // Multiply by constant factor (gamma-1)/R
                        Vmath::Smul(nBndEdgePts, (m_gamma - 1)/m_gasConstant,
                                    &scalarVariables[nScalars-1][id2], 1,
                                    &scalarVariables[nScalars-1][id2], 1);
                    }

                    // For Dirichlet boundary condition: uflux = u_bcs
                    if (fields[nScalars]->GetBndConditions()[j]->
                        GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet &&
                        !boost::iequals(
                            fields[nScalars]->GetBndConditions()[j]
                            ->GetUserDefined(), "WallAdiabatic"))
                    {
                        Vmath::Vcopy(nBndEdgePts,
                                     &scalarVariables[nScalars-1][id2], 1,
                                     &penaltyfluxO1[nScalars-1][id2], 1);

                    }

                    // For Neumann boundary condition: uflux = u_+
                    else if (((fields[nScalars]->GetBndConditions()[j])->
                              GetBoundaryConditionType() ==
                              SpatialDomains::eNeumann) ||
                             boost::iequals(
                                 fields[nScalars]->GetBndConditions()[j]
                                 ->GetUserDefined(), "WallAdiabatic"))
                    {
                        Vmath::Vcopy(nBndEdgePts,
                                     &uplus[nScalars-1][id2], 1,
                                     &penaltyfluxO1[nScalars-1][id2], 1);

                    }
                }
            }
        }

        /**
         * @brief Build the numerical flux for the 2nd order derivatives
         *
         */
        void DiffusionLFRNS::v_NumericalFluxO2(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nVariables = fields.size();
            int nDim       = fields[0]->GetCoordim(0);

            Array<OneD, NekDouble > Fwd(nTracePts);
            Array<OneD, NekDouble > Bwd(nTracePts);
            Array<OneD, NekDouble > Vn (nTracePts, 0.0);

            Array<OneD, NekDouble > qFwd     (nTracePts);
            Array<OneD, NekDouble > qBwd     (nTracePts);
            Array<OneD, NekDouble > qfluxtemp(nTracePts, 0.0);

            // Get the normal velocity Vn
            for (i = 0; i < nDim; ++i)
            {
                fields[0]->ExtractTracePhys(ufield[i], m_traceVel[i]);
                Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1,
                             m_traceVel[i], 1, Vn, 1, Vn, 1);
            }

            // Evaulate Riemann flux
            // qflux = \hat{q} \cdot u = q \cdot n
            // Notice: i = 1 (first row of the viscous tensor is zero)
            for (i = 1; i < nVariables; ++i)
            {
                qflux[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
                for (j = 0; j < nDim; ++j)
                {
                    // Compute qFwd and qBwd value of qfield in position 'ji'
                    fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);

                    // Get Riemann flux of qflux --> LDG implies upwind
                    fields[i]->GetTrace()->Upwind(Vn, qBwd, qFwd, qfluxtemp);

                    // Multiply the Riemann flux by the trace normals
                    Vmath::Vmul(nTracePts, m_traceNormals[j], 1, qfluxtemp, 1,
                                qfluxtemp, 1);

                    // Impose weak boundary condition with flux
                    if (fields[0]->GetBndCondExpansions().size())
                    {
                        v_WeakPenaltyO2(fields, i, j, qfield[j][i], qfluxtemp);
                    }

                    // Store the final flux into qflux
                    Vmath::Vadd(nTracePts, qfluxtemp, 1, qflux[i], 1,
                                qflux[i], 1);
                }
            }
        }


        /**
         * @brief Imposes appropriate bcs for the 2nd order derivatives
         *
         */
        void DiffusionLFRNS::v_WeakPenaltyO2(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const int                                          dir,
            const Array<OneD, const NekDouble>                &qfield,
                  Array<OneD,       NekDouble>                &penaltyflux)
        {
            int cnt = 0;
            int nBndEdges, nBndEdgePts;
            int i, e;
            int id2;

            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();
            int nBndRegions = fields[var]->GetBndCondExpansions().size();

            Array<OneD, NekDouble > uterm(nTracePts);
            Array<OneD, NekDouble > qtemp(nTracePts);

            // Extract the physical values of the solution at the boundaries
            fields[var]->ExtractTracePhys(qfield, qtemp);

            // Loop on the boundary regions to apply appropriate bcs
            for (i = 0; i < nBndRegions; ++i)
            {
                // Number of boundary regions related to region 'i'
                nBndEdges = fields[var]->
                GetBndCondExpansions()[i]->GetExpSize();

                if (fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::ePeriodic)
                {
                    continue;
                }

                // Weakly impose bcs by modifying flux values
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();

                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondIDToGlobalTraceID(cnt++));


                    // In case of Dirichlet bcs:
                    // uflux = gD
                    if(fields[var]->GetBndConditions()[i]->
                       GetBoundaryConditionType() == SpatialDomains::eDirichlet
                       && !boost::iequals(fields[var]->GetBndConditions()[i]
                                          ->GetUserDefined(), "WallAdiabatic"))
                    {
                        Vmath::Vmul(nBndEdgePts,
                                    &m_traceNormals[dir][id2], 1,
                                    &qtemp[id2], 1,
                                    &penaltyflux[id2], 1);
                    }
                    // 3.4) In case of Neumann bcs:
                    // uflux = u+
                    else if((fields[var]->GetBndConditions()[i])->
                            GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        ASSERTL0(false,
                                 "Neumann bcs not implemented for LFRNS");
                    }
                    else if(boost::iequals(fields[var]->GetBndConditions()[i]
                                           ->GetUserDefined(), "WallAdiabatic"))
                    {
                        if ((var == m_spaceDim + 1))
                        {
                            Vmath::Zero(nBndEdgePts, &penaltyflux[id2], 1);
                        }
                        else
                        {

                            Vmath::Vmul(nBndEdgePts,
                                        &m_traceNormals[dir][id2], 1,
                                        &qtemp[id2], 1,
                                        &penaltyflux[id2], 1);
                        }

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
        void DiffusionLFRNS::v_DerCFlux_1D(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, const NekDouble>                &flux,
            const Array<OneD, const NekDouble>                &iFlux,
                  Array<OneD,       NekDouble>                &derCFlux)
        {
            boost::ignore_unused(nConvectiveFields);

            int n;
            int nLocalSolutionPts, phys_offset;

            Array<OneD,       NekDouble> auxArray1, auxArray2;
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;

            LibUtilities::PointsKeyVector ptsKeys;
            LibUtilities::BasisSharedPtr Basis;
            Basis = fields[0]->GetExp(0)->GetBasis(0);

            int nElements = fields[0]->GetExpSize();
            int nPts      = fields[0]->GetTotPoints();

            // Arrays to store the derivatives of the correction flux
            Array<OneD, NekDouble> DCL(nPts/nElements, 0.0);
            Array<OneD, NekDouble> DCR(nPts/nElements, 0.0);

            // Arrays to store the intercell numerical flux jumps
            Array<OneD, NekDouble>  JumpL(nElements);
            Array<OneD, NekDouble>  JumpR(nElements);

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
                JumpL[n] =  iFlux[n] - tmpFluxVertex;

                fields[0]->GetExp(n)->GetVertexPhysVals(1, tmparrayX1,
                                                        tmpFluxVertex);
                JumpR[n] =  iFlux[n+1] - tmpFluxVertex;
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
        void DiffusionLFRNS::v_DerCFlux_2D(
            const int                                         nConvectiveFields,
            const int                                         direction,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, const NekDouble>                &flux,
            const Array<OneD,       NekDouble>                &iFlux,
                  Array<OneD,       NekDouble>                &derCFlux)
        {
            boost::ignore_unused(nConvectiveFields);

            int n, e, i, j, cnt;

            Array<OneD, const NekDouble> jac;

            int nElements   = fields[0]->GetExpSize();
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
                    fields[0]->GetExp(n)->GetEdgePhysVals(
                                                    e, elmtToTrace[n][e],
                                                    flux + phys_offset,
                                                    auxArray1 = tmparray);

                    // Splitting inarray into the 'j' direction
                    /*Vmath::Vmul(nEdgePts,
                     &m_traceNormals[direction][trace_offset], 1,
                     &tmparray[0], 1, &tmparray[0], 1);*/

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
                            Vmath::Reverse(nEdgePts, &jacEdge[0], 1,
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


                // Sum all the edge contributions since I am passing the
                // component of the flux x and y already. So I should not
                // need to sum E0-E2 to get the derivative wrt xi2 and E1-E3
                // to get the derivative wrt xi1

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
        void DiffusionLFRNS::v_DivCFlux_2D(
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

        void DiffusionLFRNS::v_DivCFlux_2D_Gauss(
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
                                if (m_traceNormals[0][trace_offset+i]
                                    != normals[0][i] ||
                                    m_traceNormals[1][trace_offset+i]
                                    != normals[1][i])
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
                                if (m_traceNormals[0][trace_offset+i]
                                    != normals[0][i] ||
                                    m_traceNormals[1][trace_offset+i]
                                    != normals[1][i])
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
                                if (m_traceNormals[0][trace_offset+i]
                                    != normals[0][i] ||
                                    m_traceNormals[1][trace_offset+i]
                                    != normals[1][i])
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
                                if (m_traceNormals[0][trace_offset+i]
                                    != normals[0][i] ||
                                    m_traceNormals[1][trace_offset+i]
                                    != normals[1][i])
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

    }//end of namespace SolverUtils
}//end of namespace Nektar
