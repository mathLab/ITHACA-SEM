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
// License for the specific language governing rights and limitations under
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

#include <SolverUtils/Diffusion/DiffusionLFR.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <iomanip>

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
         * This routine calls the virtual functions #v_SetupMetrics, 
         * #v_SetupCFunctions and #v_SetupInterpolationMatrices to 
         * initialise the objects needed by DiffusionLFR. 
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
            v_SetupInterpolationMatrices(pSession, pFields);
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
            int n;
            int nquad0, nquad1, nquad2;
            int nElements = pFields[0]->GetExpSize();            
            int nDim      = pFields[0]->GetCoordim(0);
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            Array<OneD, NekDouble> auxArray1;
            
            switch (nDim)
            {
                case 1:
                {
                    // nothing to do for 1D problems
                    break;
                }
                case 2:
                {
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
            int i, n;
            NekDouble c0, c1, c2;
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
         * @brief Setup the interpolation matrices to compute the solution 
         * as well as the fluxes at the interfaces in case of Gauss points.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         *
         * \todo Complete the implementation in a more efficient way.
         */
        void DiffusionLFR::v_SetupInterpolationMatrices(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            LibUtilities::BasisSharedPtr Basis;
            Basis = pFields[0]->GetExp(0)->GetBasis(0);
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            
            if (Basis->GetPointsType() == LibUtilities::eGaussGaussLegendre)
            {
                int n;
                int nElements = pFields[0]->GetExpSize();            
                int nDim      = pFields[0]->GetCoordim(0);
                
                switch (nDim)
                {
                    case 1:
                    {                                                      
                        for (n = 0; n < nElements; ++n)
                        {
                            base    = pFields[0]->GetExp(n)->GetBase();
                            Array<OneD, NekDouble> coords_m(3, 0.0);
                            Array<OneD, NekDouble> coords_p(3, 0.0);
                            coords_m[0] = -1.0;
                            coords_p[0] =  1.0;
                            
                            m_Ixm = base[0]->GetI(coords_m);;
                            m_Ixp = base[0]->GetI(coords_p);;
                            
                        }
                        break;
                    }
                    case 2:
                    {
                        ASSERTL0(false,"2DFR Gauss points not implemented yet");
                        
                        break;
                    }
                    case 3:
                    {
                        ASSERTL0(false,"3DFR Gauss points not implemented yet");
                        
                        break;
                    }
                    default:
                    {
                        ASSERTL0(false,"Expansion dimension not recognised");
                        break;
                    }
                }
            }
        }


        
        void DiffusionLFR::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            //cout<<setprecision(16);
            int i, j, n, z;
            int nLocalSolutionPts, phys_offset;
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            Array<OneD,       NekDouble> auxArray1, auxArray2, auxArray3;
                        
            LibUtilities::BasisSharedPtr Basis;
            Basis = fields[0]->GetExp(0)->GetBasis(0);
            
            int nElements    = fields[0]->GetExpSize();            
            int nDim         = fields[0]->GetCoordim(0);                        
            int nSolutionPts = fields[0]->GetTotPoints();
            int nTracePts    = fields[0]->GetTrace()->GetTotPoints();
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for (i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
                        
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > iuFluxO1(
                                                            nConvectiveFields);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                iuFluxO1[i] = Array<OneD, Array<OneD, NekDouble> >(nDim);

                for (j = 0; j < nDim; ++j)
                {
                    iuFluxO1[i][j] = Array<OneD, NekDouble>(nTracePts, 0.0);
                }
            }

            // Compute interface numerical fluxes for inarray in physical space 
            v_NumFluxforScalar(fields, inarray, iuFluxO1);

            switch(nDim)
            {
                // 1D problems 
                case 1:
                {
                    // Variable initialisation
                    Array<OneD, Array<OneD, NekDouble> > DinarrayO1(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > DCorrFluxO1(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > stdDCorrFluxO1(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > derivativesO1(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> >iqFluxO2(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > stdDCorrFluxO2(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > DderivativesO1(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > DCorrFluxO2(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > derivativesO2(
                                                            nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >tmp(
                                                            nConvectiveFields);
                    
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        DinarrayO1[i]     = Array<OneD, NekDouble>(
                                                            nSolutionPts, 0.0);
                        DCorrFluxO1[i]    = Array<OneD, NekDouble>(
                                                            nSolutionPts, 0.0);
                        stdDCorrFluxO1[i] = Array<OneD, NekDouble>(
                                                            nSolutionPts, 0.0);
                        derivativesO1[i]  = Array<OneD, NekDouble>(
                                                            nSolutionPts, 0.0);
                        iqFluxO2[i]       = Array<OneD, NekDouble>(
                                                            nTracePts, 0.0);
                        stdDCorrFluxO2[i] = Array<OneD, NekDouble>(
                                                            nSolutionPts, 0.0);
                        DderivativesO1[i] = Array<OneD, NekDouble>(
                                                            nSolutionPts, 0.0);
                        DCorrFluxO2[i]    = Array<OneD, NekDouble>(
                                                            nSolutionPts, 0.0);
                        derivativesO2[i]  = Array<OneD, NekDouble>(
                                                            nSolutionPts, 0.0);
                        
                        tmp[i] = Array<OneD, Array<OneD, NekDouble> >(1);
                        tmp[i][0] = Array<OneD, NekDouble>(nSolutionPts, 0.0);

                        // Computing the physical first-order discountinuous 
                        // derivative
                        for (n = 0; n < nElements; n++)
                        {
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            
                            fields[i]->GetExp(n)->PhysDeriv(0, 
                                auxArray1 = inarray[i] + phys_offset, 
                                auxArray2 = DinarrayO1[i] + phys_offset);
                        }
                        
                        // Computing the standard first-order correction 
                        // derivative
                        v_DerCFlux_1D(nConvectiveFields, fields, inarray[i], 
                                      iuFluxO1[i][0], stdDCorrFluxO1[i]);
                        
                        // Computing the first-order derivative of the auxiliary
                        // equation(s) in the physical space
                        for (n = 0; n < nElements; n++) 
                        {
                            nLocalSolutionPts = fields[0]->
                                                    GetExp(n)->GetTotPoints();
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            jac = fields[0]->GetExp(n)->GetGeom1D()->GetJac();
                            
                            Vmath::Smul(nLocalSolutionPts, 
                                        (1/jac[0])*(1/jac[0]), 
                                        &stdDCorrFluxO1[i][phys_offset], 1, 
                                        &DCorrFluxO1[i][phys_offset], 1);
                            
                            Vmath::Vadd(nLocalSolutionPts, 
                                        &DCorrFluxO1[i][phys_offset], 1, 
                                        &DinarrayO1[i][phys_offset], 1, 
                                        &derivativesO1[i][phys_offset], 1);
                            
                            Vmath::Vcopy(nLocalSolutionPts,
                                         &derivativesO1[i][phys_offset], 1,
                                         &tmp[i][0]
                                         [phys_offset], 1);
                        }
                    }

                    // Computing interface numerical fluxes for derivativesO1 
                    // in physical space 
                    v_NumFluxforVector(fields, inarray, tmp, iqFluxO2);
                                         
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        // Computing the physical second-order discountinuous 
                        // derivative
                        for (n = 0; n < nElements; n++)
                        {
                            phys_offset = fields[0]->GetPhys_Offset(n);

                            fields[i]->GetExp(n)->PhysDeriv(0, 
                                auxArray1 = derivativesO1[i] + phys_offset, 
                                auxArray2 = DderivativesO1[i] + phys_offset);
                        }
                        
                        // Computing the standard second-order correction 
                        // derivative
                        v_DerCFlux_1D(nConvectiveFields,
                                      fields, 
                                      derivativesO1[i], 
                                      iqFluxO2[i],
                                      stdDCorrFluxO2[i]);
                        
                        // Computing the second-order derivative
                        for (n = 0; n < nElements; n++) 
                        {
                            nLocalSolutionPts = fields[0]->
                                                    GetExp(n)->GetTotPoints();
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            jac = fields[0]->GetExp(n)->GetGeom1D()->GetJac();
                            
                            Vmath::Smul(nLocalSolutionPts, 
                                        (1/jac[0])*(1/jac[0]),
                                        &stdDCorrFluxO2[i][phys_offset], 1, 
                                        &DCorrFluxO2[i][phys_offset], 1);
                                                        
                            Vmath::Vadd(nLocalSolutionPts, 
                                        &DCorrFluxO2[i][phys_offset], 1, 
                                        &DderivativesO1[i][phys_offset], 1, 
                                        &outarray[i][phys_offset], 1); 
                        }
                    }
                    break;
                }
                // 2D problems 
                case 2:
                { 
                    // Variable initialisation
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  
                                        stdDinarrayO1(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  
                                        DinarrayO1(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > 
                                        stdDCorrFluxO1(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                        DCorrFluxO1(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                        derivativesO1(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                        BderivativesO1(nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> >
                                        iqFluxO2(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                        stdDderivativesO1(nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > 
                                        divFD(nConvectiveFields);
                    Array<OneD, Array<OneD, NekDouble> > 
                                        divFC(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                        tmp(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                        tmp1(nConvectiveFields);
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                        tmp2(nConvectiveFields);
                    
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        stdDinarrayO1[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        DinarrayO1[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        stdDCorrFluxO1[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        DCorrFluxO1[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        derivativesO1[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        BderivativesO1[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        tmp[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        tmp1[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        tmp2[i] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nDim);
                        
                        stdDderivativesO1[i] = 
                        Array<OneD, Array<OneD, NekDouble> >(nDim); 
                        iqFluxO2[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                        divFD[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                        divFC[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);

                        for (j = 0; j < nDim; ++j)
                        {
                            stdDinarrayO1[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                            DinarrayO1[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                            stdDCorrFluxO1[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                            DCorrFluxO1[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0); 
                            derivativesO1[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                            BderivativesO1[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                            tmp[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                            tmp1[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                            tmp2[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                            stdDderivativesO1[i][j] = 
                            Array<OneD, NekDouble>(nSolutionPts, 0.0);
                        }
                    }
                    /*
                    for(i = 0; i < nConvectiveFields; ++i)
                    {               
                        for (j = 0; j < nDim; ++j)
                        {
                            // Computing the physical first-order discountinuous
                            // derivatives 
                            for (n = 0; n < nElements; n++)
                            {
                                phys_offset = fields[0]->GetPhys_Offset(n);

                                fields[i]->GetExp(n)->PhysDeriv(j, 
                                auxArray1 = inarray[i] + phys_offset, 
                                auxArray2 = DinarrayO1[i][j] + phys_offset); 
                            }
                                                    
                            // Computing the standard first-order correction 
                            // derivatives
                            v_DerCFlux_2D(nConvectiveFields, j, fields, 
                                          inarray[i], iuFluxO1[i][j], 
                                          stdDCorrFluxO1[i][j]);
                        }
                     */
                    
                    for(i = 0; i < nConvectiveFields; ++i)
                    {               
                        for (j = 0; j < nDim; ++j)
                        {
                            // Computing the physical first-order discountinuous
                            // derivatives 
                            for (n = 0; n < nElements; n++)
                            { 
                                // Discontinuous flux
                                nLocalSolutionPts = fields[0]->GetExp(n)->
                                GetTotPoints();
                                phys_offset = fields[0]->GetPhys_Offset(n);
                                
                                jac  = fields[0]->GetExp(n)->GetGeom2D()->
                                GetJac();
                                gmat = fields[0]->GetExp(n)->GetGeom2D()->
                                GetGmat();
                                
                                // Temporary vectors
                                Array<OneD, NekDouble> u1_hat(
                                                    nLocalSolutionPts, 0.0);
                                Array<OneD, NekDouble> u2_hat(
                                                    nLocalSolutionPts, 0.0);
                                
                                switch (j)
                                {
                                    case 0:
                                    {
                                        if (fields[0]->GetExp(n)->GetGeom2D()->
                                            GetGtype() == 
                                            SpatialDomains::eDeformed)
                                        {
                                            for (z = 0; z < nLocalSolutionPts; 
                                                 z++)
                                            {
                                                u1_hat[z] =
                                                (inarray[i][z+phys_offset]
                                                 * gmat[0][z]) * jac[z];
                                                
                                                u2_hat[z] =
                                                (inarray[i][z+phys_offset]
                                                 * gmat[1][z]) * jac[z];
                                            }
                                        }
                                        else
                                        {
                                            for (z = 0; z < nLocalSolutionPts; 
                                                 z++)
                                            {
                                                u1_hat[z] =
                                                (inarray[i][z+phys_offset]
                                                 * gmat[0][0]) * jac[0];
                                                
                                                u2_hat[z] =
                                                (inarray[i][z+phys_offset]
                                                 * gmat[1][0])*jac[0];
                                            }
                                        }
                                        break;
                                    }
                                    case 1:
                                    {
                                        if (fields[0]->GetExp(n)->GetGeom2D()->
                                            GetGtype() == 
                                            SpatialDomains::eDeformed)
                                        {
                                            for (z = 0; z < nLocalSolutionPts; 
                                                 z++)
                                            {
                                                u1_hat[z] =
                                                (inarray[i][z+phys_offset]
                                                 * gmat[2][z]) * jac[z];
                                                
                                                u2_hat[z] =
                                                (inarray[i][z+phys_offset]
                                                 * gmat[3][z]) * jac[z];
                                            }
                                        }
                                        else
                                        {
                                            for (z = 0; z < nLocalSolutionPts; 
                                                 z++)
                                            {
                                                u1_hat[z] =
                                                (inarray[i][z+phys_offset]
                                                 * gmat[2][0]) * jac[0];
                                                
                                                u2_hat[z] =
                                                (inarray[i][z+phys_offset]
                                                 * gmat[3][0])*jac[0];
                                            }
                                        }
                                        break;
                                    }
                                }
                                
                                fields[i]->GetExp(n)->StdPhysDeriv(0,
                                    auxArray1 = u1_hat,
                                    auxArray2 = tmp1[i][j] + phys_offset);
                                
                                fields[i]->GetExp(n)->StdPhysDeriv(1,
                                    auxArray1 = u2_hat,
                                    auxArray2 = tmp2[i][j] + phys_offset);
                                
                                
                                /*fields[i]->GetExp(n)->PhysDeriv(j,
                                 auxArray1 = inarray[i] + phys_offset,
                                 auxArray2 = DinarrayO1[i][j] + phys_offset);*/ 
                            }
                            
                            Vmath::Vadd(nSolutionPts,
                                        auxArray1 = tmp1[i][j], 1,
                                        auxArray2 = tmp2[i][j], 1,
                                        DinarrayO1[i][j], 1);
                            
                            // Multiply by the metric terms
                            for (n = 0; n < nElements; ++n)
                            {
                                nLocalSolutionPts = fields[0]->
                                GetExp(n)->GetTotPoints();
                                phys_offset = fields[0]->GetPhys_Offset(n);
                                
                                jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                                gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                                                                
                                if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype()
                                    == SpatialDomains::eDeformed)
                                {
                                    for (z = 0; z < nLocalSolutionPts; z++)
                                    {
                                        DinarrayO1[i][j][phys_offset + z] =
                                        DinarrayO1[i][j][phys_offset + z]/jac[z];
                                    }
                                }
                                else
                                {
                                    Vmath::Smul(nLocalSolutionPts, 1/jac[0],
                                                auxArray1 = DinarrayO1[i][j]
                                                + phys_offset, 1,
                                                auxArray2 = DinarrayO1[i][j]
                                                + phys_offset, 1);
                                }
                            }
                            
                            
                            // Computing the standard first-order correction 
                            // derivatives
                            v_DerCFlux_2D(nConvectiveFields, j, fields, 
                                          inarray[i], iuFluxO1[i][j], 
                                          stdDCorrFluxO1[i][j]);
                        }
                        
                        // Multiplying derivatives by B matrix to get auxiliary 
                        // variables
                        for (n = 0; n < nElements; n++)
                        {
                            nLocalSolutionPts = fields[0]->GetExp(n)->
                                                                GetTotPoints();
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                            jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                            
                            if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype() 
                                == SpatialDomains::eDeformed)
                                
                            {                                
                                for (j = 0; j < nLocalSolutionPts; ++j)
                                {   
                                    // std(q1)
                                    BderivativesO1[i][0][phys_offset+j] =  
                                    (gmat[0][j]*gmat[0][j] + 
                                     gmat[2][j]*gmat[2][j]) *
                                    stdDCorrFluxO1[i][0][phys_offset+j] +
                                    (gmat[1][j]*gmat[0][j] + 
                                     gmat[3][j]*gmat[2][j]) *
                                    stdDCorrFluxO1[i][1][phys_offset+j];
                                    
                                    // std(q2)
                                    BderivativesO1[i][1][phys_offset+j] =  
                                    (gmat[1][j]*gmat[0][j] + 
                                     gmat[3][j]*gmat[2][j]) *
                                    stdDCorrFluxO1[i][0][phys_offset+j] +
                                    (gmat[1][j]*gmat[1][j] + 
                                     gmat[3][j]*gmat[3][j]) *
                                    stdDCorrFluxO1[i][1][phys_offset+j];
                                }
                            }
                            else
                            {                                
                                // std(q1)
                                Vmath::Smul(nLocalSolutionPts, 
                                (gmat[0][0]*gmat[0][0] + gmat[2][0]*gmat[2][0]), 
                                &stdDCorrFluxO1[i][0][phys_offset], 1, 
                                &BderivativesO1[i][0][phys_offset], 1);
                                
                                Vmath::Smul(nLocalSolutionPts, 
                                (gmat[1][0]*gmat[0][0] + gmat[3][0]*gmat[2][0]), 
                                &stdDCorrFluxO1[i][1][phys_offset], 1, 
                                &tmp[i][0][phys_offset], 1);
                                
                                Vmath::Vadd(nLocalSolutionPts,
                                &BderivativesO1[i][0][phys_offset], 1,
                                &tmp[i][0][phys_offset], 1,
                                &BderivativesO1[i][0][phys_offset], 1);
                                
                                // std(q2)
                                Vmath::Smul(nLocalSolutionPts, 
                                (gmat[1][0]*gmat[0][0] + gmat[3][0]*gmat[2][0]), 
                                &stdDCorrFluxO1[i][0][phys_offset], 1, 
                                &BderivativesO1[i][1][phys_offset], 1);
                                
                                Vmath::Smul(nLocalSolutionPts, 
                                (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0]), 
                                &stdDCorrFluxO1[i][1][phys_offset], 1, 
                                &tmp[i][1][phys_offset], 1);
                                
                                Vmath::Vadd(nLocalSolutionPts,
                                &BderivativesO1[i][1][phys_offset], 1,
                                &tmp[i][1][phys_offset], 1,
                                &BderivativesO1[i][1][phys_offset], 1);
                            }
                        }
                        
                        // Multiplying derivatives by A^(-1) to get back 
                        // into the physical space 
                        for (n = 0; n < nElements; ++n)
                        {
                            nLocalSolutionPts = fields[0]->GetExp(n)->
                                                            GetTotPoints();
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                            
                            if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype() 
                                == SpatialDomains::eDeformed)
                            {
                                for (j = 0; j < nLocalSolutionPts; j++)
                                {
                                    // q1 = A11^(-1)*std(q1) + A12^(-1)*std(q2)
                                    stdDCorrFluxO1[i][0][phys_offset+j] = 
                                    gmat[3][j] * 
                                    BderivativesO1[i][0][phys_offset+j] - 
                                    gmat[2][j] * 
                                    BderivativesO1[i][1][phys_offset+j];
                                    
                                    // q2 = A21^(-1)*std(q1) + A22^(-1)*std(q2)
                                    stdDCorrFluxO1[i][1][phys_offset+j] = 
                                    gmat[0][j] * 
                                    BderivativesO1[i][1][phys_offset+j] - 
                                    gmat[1][j] * 
                                    BderivativesO1[i][0][phys_offset+j];
                                }
                            }
                            else
                            {
                                // q1
                                Vmath::Smul(nLocalSolutionPts, gmat[3][0], 
                                        &BderivativesO1[i][0][phys_offset], 1, 
                                        &stdDCorrFluxO1[i][0][phys_offset], 1);
                                
                                Vmath::Smul(nLocalSolutionPts, gmat[2][0], 
                                        &BderivativesO1[i][1][phys_offset], 1, 
                                        &tmp[i][0][phys_offset], 1);
                                
                                Vmath::Vsub(nLocalSolutionPts, 
                                        &stdDCorrFluxO1[i][0][phys_offset], 1,
                                        &tmp[i][0][phys_offset], 1,
                                        &stdDCorrFluxO1[i][0][phys_offset], 1);
                                
                                // q2
                                Vmath::Smul(nLocalSolutionPts, gmat[1][0], 
                                        &BderivativesO1[i][0][phys_offset], 1, 
                                        &stdDCorrFluxO1[i][1][phys_offset], 1);
                                
                                Vmath::Smul(nLocalSolutionPts, gmat[0][0], 
                                        &BderivativesO1[i][1][phys_offset], 1, 
                                        &tmp[i][1][phys_offset], 1);
                                
                                Vmath::Vsub(nLocalSolutionPts, 
                                        &tmp[i][1][phys_offset], 1,
                                        &stdDCorrFluxO1[i][1][phys_offset], 1,
                                        &stdDCorrFluxO1[i][1][phys_offset], 1);
                            }
                        }
                        
                        // Computing the physical first-order derivatives
                        for (j = 0; j < nDim; ++j)
                        {
                            Vmath::Vadd(nSolutionPts, 
                                        &DinarrayO1[i][j][0], 1,
                                        &stdDCorrFluxO1[i][j][0], 1, 
                                        &derivativesO1[i][j][0], 1);
                        }
                    }
                                           
                    // Computing interface numerical fluxes for derivativesO1 
                    // in physical space 
                    v_NumFluxforVector(fields, inarray, derivativesO1, 
                                       iqFluxO2);

                    // Computing the standard second-order discontinuous 
                    // derivatives 
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        for (n = 0; n < nElements; n++)
                        {
                            // Get number of solution points and phys_offset
                            nLocalSolutionPts = fields[0]->GetExp(n)->
                                                                GetTotPoints();
                            phys_offset = fields[0]->GetPhys_Offset(n);
                        
                            // Get jacobian and gmats
                            jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                            gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                        
                            // Temporary vectors
                            Array<OneD, NekDouble>f_hat(nLocalSolutionPts, 0.0);
                            Array<OneD, NekDouble>g_hat(nLocalSolutionPts, 0.0);
                        
                            if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype()
                                == SpatialDomains::eDeformed)
                            {
                                for (j = 0; j < nLocalSolutionPts; j++)
                                {
                                    f_hat[j] = 
                                    (derivativesO1[i][0][j+phys_offset] * 
                                     gmat[0][j] + 
                                     derivativesO1[i][1][j+phys_offset] * 
                                     gmat[2][j])*jac[j];
                                
                                    g_hat[j] = 
                                    (derivativesO1[i][0][j+phys_offset] * 
                                     gmat[1][j] + 
                                     derivativesO1[i][1][j+phys_offset] * 
                                     gmat[3][j])*jac[j];
                                }
                            }
                            else
                            {
                                for (j = 0; j < nLocalSolutionPts; j++)
                                {
                                    f_hat[j] = 
                                    (derivativesO1[i][0][j+phys_offset] * 
                                     gmat[0][0] + 
                                     derivativesO1[i][1][j+phys_offset] * 
                                     gmat[2][0])*jac[0];
                                
                                    g_hat[j] = 
                                    (derivativesO1[i][0][j+phys_offset] * 
                                     gmat[1][0] + 
                                     derivativesO1[i][1][j+phys_offset] * 
                                     gmat[3][0])*jac[0];
                                }
                            }
                        
                            fields[0]->GetExp(n)->StdPhysDeriv(0, 
                            auxArray1 = f_hat, 
                            auxArray2 = stdDderivativesO1[i][0] + phys_offset); 
                        
                            fields[0]->GetExp(n)->StdPhysDeriv(1, 
                            auxArray1 = g_hat, 
                            auxArray2 = stdDderivativesO1[i][1] + phys_offset);
                        }
                        
                        // Divergence of the standard discontinuous flux
                        Vmath::Vadd(nSolutionPts, 
                                    auxArray1 = stdDderivativesO1[i][0], 1,
                                    auxArray2 = stdDderivativesO1[i][1], 1, 
                                    auxArray3 = divFD[i], 1);
                        
                        // Divergence of the standard correction flux
                        v_DivCFlux_2D(nConvectiveFields,
                                      fields, 
                                      derivativesO1[i][0], 
                                      derivativesO1[i][1], 
                                      iqFluxO2[i], 
                                      divFC[i]);
                        
                        // Divergence of the standard final flux
                        Vmath::Vadd(nSolutionPts, 
                                    auxArray1 = divFD[i], 1,
                                    auxArray2 = divFC[i], 1, 
                                    auxArray3 = outarray[i], 1);
                                                
                        // Multiplying by the metric terms to get back into 
                        // physical space
                        for (n = 0; n < nElements; ++n)
                        {
                            nLocalSolutionPts = fields[0]->
                            GetExp(n)->GetTotPoints();
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            
                            jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                            gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                            
                            if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype() 
                                == SpatialDomains::eDeformed)
                            {
                                for (j = 0; j < nLocalSolutionPts; j++)
                                {
                                    outarray[i][phys_offset + j] = 
                                    outarray[i][phys_offset + j]/jac[j];
                                }
                            }
                            else
                            {
                                Vmath::Smul(nLocalSolutionPts, 1/jac[0], 
                                    auxArray1 = outarray[i] + phys_offset, 1, 
                                    auxArray2 = outarray[i] + phys_offset, 1);
                            }
                        }
                    }//Close loop on nConvectiveFields
                    break;
                }   
                // 3D-Problems 
                case 3:
                {
                    ASSERTL0(false, "3D FRDG case not implemented yet");
                    break;
                }
            }
        }
        
        void DiffusionLFR::v_NumFluxforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nvariables = fields.num_elements();
            int nDim       = fields[0]->GetCoordim(0);  
            NekDouble time = 0.0;
            
            Array<OneD, NekDouble > Fwd     (nTracePts);
            Array<OneD, NekDouble > Bwd     (nTracePts);
            Array<OneD, NekDouble > Vn      (nTracePts, 0.0);
            Array<OneD, NekDouble > fluxtemp(nTracePts, 0.0);
                        
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
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
                    
                    fields[i]->GetTrace()->Upwind(/*m_traceNormals[j]*/Vn, Fwd, Bwd, fluxtemp);
                    
                    // Imposing weak boundary condition with flux
                    // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd
                    
                    // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd
                    
                    if(fields[0]->GetBndCondExpansions().num_elements())
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
                    
                    /*
                    Vmath::Vmul(nTracePts, 
                                m_traceNormals[j], 1, 
                                fluxtemp, 1, 
                                uflux[i][j], 1);
                    */
                    Vmath::Vcopy(nTracePts, 
                                 fluxtemp, 1, 
                                 uflux[i][j], 1);
                }
            }
        }
        
        
        
        void DiffusionLFR::v_WeakPenaltyforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const Array<OneD, const NekDouble>                &ufield,
                  Array<OneD,       NekDouble>                &penaltyflux)
        {
            // Variable initialisation
            int i, j, e, id1, id2;
            int nBndEdgePts, nBndEdges;
            int cnt         = 0;
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();
            int nDim = fields[0]->GetCoordim(0);  
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
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    // Offset of the boundary expansion
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    // Offset of the trace space related to boundary expansion
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
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
        
        void DiffusionLFR::v_NumFluxforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nvariables = fields.num_elements();
            int nDim       = fields[0]->GetCoordim(0);  
            
            NekDouble C11 = 0.0;
            Array<OneD, NekDouble > Fwd(nTracePts);
            Array<OneD, NekDouble > Bwd(nTracePts);
            Array<OneD, NekDouble > Vn (nTracePts, 0.0);
            
            Array<OneD, NekDouble > qFwd     (nTracePts);
            Array<OneD, NekDouble > qBwd     (nTracePts);
            Array<OneD, NekDouble > qfluxtemp(nTracePts, 0.0);            
            Array<OneD, NekDouble > uterm(nTracePts);
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
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
                    
                    fields[i]->GetTrace()->Upwind(/*m_traceNormals[j]*/Vn, qBwd, qFwd, qfluxtemp);
                    
                    Vmath::Vmul(nTracePts, 
                                m_traceNormals[j], 1, 
                                qfluxtemp, 1, 
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
                    if (fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyforVector(fields, i, j, 
                                               qfield[i][j], 
                                               qfluxtemp, C11);
                    }
                    
                    // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                    // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
                    // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
                    Vmath::Vadd(nTracePts, 
                                qfluxtemp, 1, 
                                qflux[i], 1, 
                                qflux[i], 1);
                }
            }
        }
        
        
        
        /**
         * Diffusion: Imposing weak boundary condition for q with flux
         *  uflux = g_D  on Dirichlet boundary condition
         *  uflux = u_Fwd  on Neumann boundary condition
         */
        void DiffusionLFR::v_WeakPenaltyforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const int                                          dir,
            const Array<OneD, const NekDouble>                &qfield,
                  Array<OneD,       NekDouble>                &penaltyflux,
            NekDouble                                          C11)
        {
            int i, j, e, id1, id2;
            int nBndEdges, nBndEdgePts;
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();
            int nDim        = fields[0]->GetCoordim(0);
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();
            
            Array<OneD, NekDouble > uterm(nTracePts);
            Array<OneD, NekDouble > qtemp(nTracePts);
            int cnt = 0;
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            fields[var]->ExtractTracePhys(qfield, qtemp);
            
            for (i = 0; i < nBndRegions; ++i)
            {
                nBndEdges = fields[var]->
                GetBndCondExpansions()[i]->GetExpSize();
                
                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < nBndEdges ; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
                    
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // For Dirichlet boundary condition: 
                    //qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
                    if(fields[var]->GetBndConditions()[i]->
                       GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vmul(nBndEdgePts, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &qtemp[id2], 1, 
                                    &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: qflux = g_N
                    else if((fields[var]->GetBndConditions()[i])->
                            GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vmul(nBndEdgePts,
                                    &m_traceNormals[dir][id2], 1, 
                                    &(fields[var]->
                                      GetBndCondExpansions()[i]->
                                      GetPhys())[id1], 1, 
                                    &penaltyflux[id2], 1);
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
            int i, n;
            int nLocalSolutionPts, phys_offset;
            
            Array<OneD,       NekDouble> auxArray1, auxArray2, auxArray3;
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            
            LibUtilities::BasisSharedPtr Basis;
            Basis = fields[0]->GetExp(0)->GetBasis(0);
            
            int nElements    = fields[0]->GetExpSize();            
            int nSolutionPts = fields[0]->GetTotPoints();
            
            // Offsets for interface operations
            int offsetStart, offsetEnd;
            
            // Arrays to store the intercell numerical flux jumps
            Array<OneD, Array<OneD, NekDouble> > JumpL(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > JumpR(nConvectiveFields);
            
            // Arrays to store the derivatives of the correction flux
            Array<OneD, NekDouble> DCL(nSolutionPts/nElements, 0.0); 
            Array<OneD, NekDouble> DCR(nSolutionPts/nElements, 0.0);
            
            // The dimension of each column of the jump arrays
            for(i = 0; i < nConvectiveFields; ++i)
            {
                JumpL[i] = Array<OneD, NekDouble>(nElements);
                JumpR[i] = Array<OneD, NekDouble>(nElements);
            }
            
            // Interpolation routine for Gauss points and fluxJumps computation                   
            if (Basis->GetPointsType() == LibUtilities::eGaussGaussLegendre)
            {
                Array<OneD, NekDouble> interpolatedFlux_m(nElements);
                Array<OneD, NekDouble> interpolatedFlux_p(nElements);
                
                
                for (n = 0; n < nElements; ++n)
                {
                    nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                    
                    Array<OneD, NekDouble> physvals(nLocalSolutionPts, 0.0);
                    physvals = flux + n*nLocalSolutionPts;
                    
                    interpolatedFlux_m[n] = Blas::Ddot(nLocalSolutionPts,
                                                       m_Ixm->GetPtr(), 1,
                                                       physvals, 1);
                    
                    interpolatedFlux_p[n] = Blas::Ddot(nLocalSolutionPts, 
                                                       m_Ixp->GetPtr(), 1, 
                                                       physvals, 1);   
                    
                    JumpL[0][n] = iFlux[n]   - interpolatedFlux_m[n];
                    JumpR[0][n] = iFlux[n+1] - interpolatedFlux_p[n];
                }
            }
            
            // FluxJumps computation without interpolation                  
            else
            {
                for(n = 0; n < nElements; ++n)
                {
                    nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                    
                    offsetStart = fields[0]->GetPhys_Offset(n);
                    offsetEnd   = offsetStart + nLocalSolutionPts - 1;
                    
                    JumpL[0][n] = iFlux[n]   - flux[offsetStart];
                    JumpR[0][n] = iFlux[n+1] - flux[offsetEnd];
                }
            }
            
            for (n = 0; n < nElements; ++n)
            {
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                phys_offset       = fields[0]->GetPhys_Offset(n);
                jac               = fields[0]->GetExp(n)->GetGeom1D()->GetJac();
                
                JumpL[0][n] = JumpL[0][n] * jac[0];
                JumpR[0][n] = JumpR[0][n] * jac[0];
                
                // Left jump multiplied by left derivative of C function
                Vmath::Smul(nLocalSolutionPts, 
                            JumpL[0][n], 
                            auxArray1 = m_dGL_xi1[n], 1, 
                            DCL, 1);
                
                // Right jump multiplied by right derivative of C function
                Vmath::Smul(nLocalSolutionPts, 
                            JumpR[0][n], 
                            auxArray1 = m_dGR_xi1[n], 1, 
                            DCR, 1);
                
                // Assembling derivative of the correction flux
                Vmath::Vadd(nLocalSolutionPts, 
                            DCL, 1, 
                            DCR, 1, 
                            auxArray1 = derCFlux + phys_offset, 1);
            }
        }
        
        /**
         * @brief Compute the derivative of the corrective flux wrt a given 
         * coordinate for 2D problems.
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
            int n, e, i, j, cnt;
            
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            
            int nElements   = fields[0]->GetExpSize();
            int nDim = fields[0]->GetCoordim(0);  
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();
            
            int trace_offset, phys_offset;
            int nLocalSolutionPts;
            int nquad0, nquad1;
            int nEdgePts;  
            
            Array<OneD, NekDouble> auxArray1, auxArray2;
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
            &elmtToTrace = fields[0]->GetTraceMap()->GetElmtToTrace();
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for (i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
            // Loop on the elements
            for (n = 0; n < nElements; ++n)
            {
                // Offset of the element on the global vector
                phys_offset = fields[0]->GetPhys_Offset(n);
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                
                jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                
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
                    const Array<OneD, const Array<OneD, NekDouble> > &normals = 
                    fields[0]->GetExp(n)->GetEdgeNormal(e);
                    
                    // Extract the edge values of the volumetric fluxes 
                    // on edge 'e' and order them accordingly to the order 
                    // of the trace space 
                    fields[0]->GetExp(n)->GetEdgePhysVals(e, elmtToTrace[n][e],
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
                    if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype() 
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
                                    cnt = (nquad0 - 1) + j*nquad0 - i;
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
                                    cnt = (nquad0*nquad1 - nquad0) + j - i*nquad0;
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
                /*
                 derCFlux[phys_offset+i] = divCFluxE0[i] + divCFluxE1[i] + 
                 divCFluxE2[i] + divCFluxE3[i];
                 */
                
                if (direction == 0)
                {
                    for (i = 0; i < nLocalSolutionPts; ++i)
                    {
                        derCFlux[phys_offset+i] = divCFluxE3[i] + divCFluxE1[i];
                    }
                }
                else if (direction == 1)
                {
                    for (i = 0; i < nLocalSolutionPts; ++i)
                    {
                        derCFlux[phys_offset+i] = divCFluxE0[i] + divCFluxE2[i];
                    }
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
            int n, e, i, j, cnt;
            
            int nElements   = fields[0]->GetExpSize();
            int nDim = fields[0]->GetCoordim(0);  
            int nTracePts   = fields[0]->GetTrace()->GetTotPoints();
            
            int nLocalSolutionPts;
            int nEdgePts;  
            int trace_offset; 
            int phys_offset;
            int nquad0;
            int nquad1;
            
            Array<OneD, NekDouble> auxArray1, auxArray2;
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
            &elmtToTrace = fields[0]->GetTraceMap()->GetElmtToTrace();
            
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            
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
                                    cnt = (nquad0 - 1) + j*nquad0 - i;
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
                                    cnt = (nquad0*nquad1 - nquad0) + j - i*nquad0;
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
                for (i = 0; i < nLocalSolutionPts; ++i)
                {
                    divCFlux[phys_offset + i] = divCFluxE0[i] + 
                    divCFluxE1[i] +
                    divCFluxE2[i] + 
                    divCFluxE3[i];
                }
            }
        } 
        
    }// close namespace SolverUtils
}// close namespace nektar++
