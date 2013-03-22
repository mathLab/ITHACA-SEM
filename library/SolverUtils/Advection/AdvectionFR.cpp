///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionFR.cpp
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
// Description: FR advection class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/AdvectionFR.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <StdRegions/StdSegExp.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <boost/math/special_functions/gamma.hpp>


namespace Nektar
{
    namespace SolverUtils
    {        
        std::string AdvectionFR::type[] = {
            GetAdvectionFactory().RegisterCreatorFunction(
                                        "FRDG", AdvectionFR::create), 
            GetAdvectionFactory().RegisterCreatorFunction(
                                        "FRSD", AdvectionFR::create), 
            GetAdvectionFactory().RegisterCreatorFunction(
                                        "FRHU", AdvectionFR::create),
            GetAdvectionFactory().RegisterCreatorFunction(
                                        "FRcmin", AdvectionFR::create),
            GetAdvectionFactory().RegisterCreatorFunction(
                                        "FRcinf", AdvectionFR::create)};
        
        /**
         * @brief AdvectionFR uses the Flux Reconstruction (FR) approach to 
         * compute the advection term. The implementation is only for segments,
         * quadrilaterals and hexahedra at the moment.
         * 
         * \todo Extension to triangles, tetrahedra and other shapes. 
         * (Long term objective) 
         */
        AdvectionFR::AdvectionFR(std::string advType):m_advType(advType)
        {
        }
        
        /**
         * @brief Initiliase AdvectionFR objects and store them before starting 
         * the time-stepping.
         * 
         * This routine calls the virtual functions #v_SetupMetrics, 
         * #v_SetupCFunctions and #v_SetupInterpolationMatrices to 
         * initialise the objects needed by AdvectionFR. 
         * 
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void AdvectionFR::v_InitObject(
                LibUtilities::SessionReaderSharedPtr        pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
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
        void AdvectionFR::v_SetupMetrics(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            int n;
            int nquad0, nquad1;
            int nElements   = pFields[0]->GetExpSize();            
            int nDimensions = pFields[0]->GetCoordim(0);
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            Array<OneD, NekDouble> auxArray1;
                        
            switch (nDimensions)
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
        void AdvectionFR::v_SetupCFunctions(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {        
            int i, n;
            NekDouble c0, c1, c2;
            int nquad0, nquad1, nquad2;
            int nmodes0, nmodes1, nmodes2;
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            
            int nElements   = pFields[0]->GetExpSize();
            int nDimensions = pFields[0]->GetCoordim(0);
            
            switch (nDimensions)
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
                        if (m_advType == "FRDG")
                        {
                            c0 = 0.0;
                        }
                        else if (m_advType == "FRSD")
                        {
                            c0 = 2.0 * p0 / ((2.0 * p0 + 1.0) * (p0 + 1.0) 
                               * (ap0 * boost::math::tgamma(p0 + 1)) 
                               * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_advType == "FRHU")
                        {
                            c0 = 2.0 * (p0 + 1.0) / ((2.0 * p0 + 1.0) * p0 
                               * (ap0 * boost::math::tgamma(p0 + 1)) 
                               * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_advType == "FRcmin")
                        {
                            c0 = -2.0 / ((2.0 * p0 + 1.0) 
                               * (ap0 * boost::math::tgamma(p0 + 1))
                               * (ap0 * boost::math::tgamma(p0 + 1)));
                        }
                        else if (m_advType == "FRinf")
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
                        if (m_advType == "FRDG")
                        {
                            c0 = 0.0;
                            c1 = 0.0;
                        }
                        else if (m_advType == "FRSD")
                        {
                            c0 = 2.0 * p0 / ((2.0 * p0 + 1.0) * (p0 + 1.0) 
                               * (ap0 * boost::math::tgamma(p0 + 1)) 
                               * (ap0 * boost::math::tgamma(p0 + 1)));
                            
                            c1 = 2.0 * p1 / ((2.0 * p1 + 1.0) * (p1 + 1.0) 
                               * (ap1 * boost::math::tgamma(p1 + 1)) 
                               * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_advType == "FRHU")
                        {
                            c0 = 2.0 * (p0 + 1.0) / ((2.0 * p0 + 1.0) * p0 
                               * (ap0 * boost::math::tgamma(p0 + 1)) 
                               * (ap0 * boost::math::tgamma(p0 + 1)));
                            
                            c1 = 2.0 * (p1 + 1.0) / ((2.0 * p1 + 1.0) * p1 
                               * (ap1 * boost::math::tgamma(p1 + 1)) 
                               * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_advType == "FRcmin")
                        {
                            c0 = -2.0 / ((2.0 * p0 + 1.0) 
                               * (ap0 * boost::math::tgamma(p0 + 1))
                               * (ap0 * boost::math::tgamma(p0 + 1)));
                            
                            c1 = -2.0 / ((2.0 * p1 + 1.0) 
                               * (ap1 * boost::math::tgamma(p1 + 1))
                               * (ap1 * boost::math::tgamma(p1 + 1)));
                        }
                        else if (m_advType == "FRinf")
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
                        if (m_advType == "FRDG")
                        {
                            c0 = 0.0;
                            c1 = 0.0;
                            c2 = 0.0;
                        }
                        else if (m_advType == "FRSD")
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
                        else if (m_advType == "FRHU")
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
                        else if (m_advType == "FRcmin")
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
                        else if (m_advType == "FRinf")
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
                            m_dGL_xi3[n][i]  = 0.5 * sign1 * m_dGL_xi3[n][i];
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
        void AdvectionFR::v_SetupInterpolationMatrices(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            LibUtilities::BasisSharedPtr Basis;
            Basis = pFields[0]->GetExp(0)->GetBasis(0);
            Array<OneD, LibUtilities::BasisSharedPtr> base;

            if (Basis->GetPointsType() == LibUtilities::eGaussGaussLegendre)
            {
                int n;
                
                int nElements   = pFields[0]->GetExpSize();            
                int nDimensions = pFields[0]->GetCoordim(0);
                
                switch (nDimensions)
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
        
        /**
         * @brief Compute the advection term at each time-step using the Flux
         * Reconstruction approach (FR).
         *
         * @param nConvectiveFields   Number of fields (i.e. independent 
         *                            variables).
         * @param fields              Pointer to fields.
         * @param advVel              Advection velocities.
         * @param inarray             Solution at the previous time-step.
         * @param outarray            Advection term to be passed at the 
         *                            time integration class.
         *
         */
        void AdvectionFR::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int i, j, n;
            int nLocalSolutionPts, phys_offset;
            
            Array<OneD,       NekDouble> auxArray1, auxArray2, auxArray3;
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            
            LibUtilities::BasisSharedPtr Basis;
            Basis = fields[0]->GetExp(0)->GetBasis(0);
            
            int nElements       = fields[0]->GetExpSize();            
            int nDimensions     = fields[0]->GetCoordim(0);
            int nSolutionPts    = fields[0]->GetTotPoints();
            int nTracePts       = fields[0]->GetTrace()->GetTotPoints();
                                       
            // Store forwards/backwards space along trace space.
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);
            
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts);
                numflux[i] = Array<OneD, NekDouble>(nTracePts);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
            
            // Computing the Riemann flux at each trace point
            m_riemann->Solve(Fwd, Bwd, numflux);
            
            // Divergence of the flux (computing the RHS)  
            switch(nDimensions)
            {
                // 1D-Problems 
                case 1:
                {                    
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > 
                        fluxvector(nConvectiveFields);
                    Array<OneD, NekDouble> DfluxvectorX1(nSolutionPts);
                    Array<OneD, NekDouble> divFC        (nSolutionPts);

                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        fluxvector[i] = Array<OneD, Array<OneD, NekDouble> >(
                            nDimensions);
                        for (j = 0; j < nDimensions; ++j)
                        {
                            fluxvector[i][j] = Array<OneD, NekDouble>(
                                nSolutionPts);
                        }
                    }
                    
                    m_fluxVector(inarray, fluxvector);

                    // Get the discontinuous flux FD ("i" is used by inarray)
                    for(i = 0; i < nConvectiveFields; ++i)
                    {     
                        for (n = 0; n < nElements; n++)
                        {
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            
                            fields[i]->GetExp(n)->PhysDeriv(
                                0, fluxvector[i][0] + phys_offset, 
                                auxArray2 = DfluxvectorX1 + phys_offset);
                        }
                        
                        v_DivCFlux_1D(nConvectiveFields,
                                      fields, 
                                      fluxvector[i][0], 
                                      numflux[i], 
                                      divFC);
                        
                        // Computation of the advection term
                        for (n = 0; n < nElements; n++) 
                        {
                            nLocalSolutionPts = fields[0]->
                                            GetExp(n)->GetTotPoints();
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            jac = fields[0]->GetExp(n)->GetGeom1D()->GetJac();
                            
                            Vmath::Smul(
                                nLocalSolutionPts, 
                                1/jac[0], 
                                divFC + phys_offset, 1, 
                                auxArray2 = outarray[i] + phys_offset, 1);
                            
                            Vmath::Vadd(
                                nLocalSolutionPts, 
                                outarray[i] + phys_offset, 1, 
                                DfluxvectorX1 + phys_offset, 1, 
                                auxArray3 = outarray[i] + phys_offset, 1); 
                        }
                    }
                    break;
                }
                // 2D-Problems 
                case 2:
                {
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                        fluxvector(nConvectiveFields);
                    Array<OneD, NekDouble> DfluxvectorX1(nSolutionPts);
                    Array<OneD, NekDouble> DfluxvectorX2(nSolutionPts);
                    Array<OneD, NekDouble> divFD(nSolutionPts);
                    Array<OneD, NekDouble> divFC(nSolutionPts);

                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        fluxvector[i] =
                            Array<OneD, Array<OneD, NekDouble> >(nDimensions);
                        for (j = 0; j < nDimensions; ++j)
                        {
                            fluxvector[i][j] =
                                Array<OneD, NekDouble>(nSolutionPts);
                        }
                    }
                    
                    m_fluxVector(inarray, fluxvector);

                    // Get the discontinuous flux FD ("i" is used by inarray)
                    for(i = 0; i < nConvectiveFields; ++i)
                    {
                        for (n = 0; n < nElements; n++)
                        {
                            // Discontinuous flux
                            nLocalSolutionPts = fields[0]->GetExp(n)->
                            GetTotPoints();
                            phys_offset = fields[0]->GetPhys_Offset(n);
                            
                            jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                            gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                            
                            // Temporary vectors
                            Array<OneD, NekDouble> f_hat(nLocalSolutionPts);
                            Array<OneD, NekDouble> g_hat(nLocalSolutionPts);

                            if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype()
                                == SpatialDomains::eDeformed)
                            {
                                for (j = 0; j < nLocalSolutionPts; j++)
                                {
                                    f_hat[j] =
                                    (fluxvector[i][0][j+phys_offset]
                                     * gmat[0][j] +
                                     fluxvector[i][1][j+phys_offset]
                                     * gmat[2][j]) * jac[j];
                                    
                                    g_hat[j] =
                                    (fluxvector[i][0][j+phys_offset]
                                     * gmat[1][j] +
                                     fluxvector[i][1][j+phys_offset]
                                     * gmat[3][j]) * jac[j];
                                }
                            }
                            else
                            {
                                for (j = 0; j < nLocalSolutionPts; j++)
                                {
                                    f_hat[j] = 
                                    (fluxvector[i][0][j+phys_offset] 
                                     * gmat[0][0] + 
                                     fluxvector[i][1][j+phys_offset] 
                                     * gmat[2][0]) * jac[0];
                                    
                                    g_hat[j] = 
                                    (fluxvector[i][0][j+phys_offset] 
                                     * gmat[1][0] + 
                                     fluxvector[i][1][j+phys_offset] 
                                     * gmat[3][0])*jac[0];
                                }
                            }

                            fields[0]->GetExp(n)->StdPhysDeriv(0, f_hat,
                                auxArray2 = DfluxvectorX1 + phys_offset);
                            fields[0]->GetExp(n)->StdPhysDeriv(1, g_hat,
                                auxArray2 = DfluxvectorX2 + phys_offset);
                        }
                        
                        // Divergence of the discontinuous flux
                        Vmath::Vadd(nSolutionPts,
                                    DfluxvectorX1, 1,
                                    DfluxvectorX2, 1,
                                    divFD, 1);

                        // Divergence of the correction flux
                        v_DivCFlux_2D(nConvectiveFields,
                                      fields,
                                      fluxvector[i][0],
                                      fluxvector[i][1],
                                      numflux[i],
                                      divFC);
                        
                        // Divergence of the final flux
                        Vmath::Vadd(nSolutionPts,
                                    divFD, 1,
                                    divFC, 1,
                                    outarray[i], 1);

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
                                for (j = 0; j < nLocalSolutionPts; ++j)
                                {
                                    outarray[i][phys_offset+j] /= jac[j];
                                }
                            }
                            else
                            {
                                Vmath::Smul(
                                    nLocalSolutionPts, 
                                    1/jac[0], 
                                    outarray[i] + phys_offset, 1, 
                                    auxArray2 = outarray[i] + phys_offset, 1);
                            }
                        }
                        

                    } // close nConvectiveFields loop
                    break;
                }
                    
                    
                // 3D-Problems 
                case 3:
                {
                    ASSERTL0(false,"3D FRDG case not implemented yet");
                    break;
                }
            }
        }
        
        /**
         * @brief Compute the divergence of the corrective flux for 1D problems.
         *
         * @param nConvectiveFields   Number of fields (i.e. independent 
         *                            variables).
         * @param fields              Pointer to fields.
         * @param fluxX1              Volumetric flux in the physical space in 
         *                            direction X1.
         * @param numericalFlux       Riemann flux in the physical space.
         * @param divCFlux            Divergence of the corrective flux for 1D
         *                            Problems.
         *
         */
        void AdvectionFR::v_DivCFlux_1D(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1, 
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux)
        {
            int i, j, n;
            int nLocalSolutionPts, phys_offset;
            
            Array<OneD,       NekDouble> auxArray1, auxArray2, auxArray3;
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            
            LibUtilities::BasisSharedPtr Basis;
            Basis = fields[0]->GetExp(0)->GetBasis(0);
            
            int nElements       = fields[0]->GetExpSize();            
            int nDimensions     = fields[0]->GetCoordim(0);
            int nSolutionPts    = fields[0]->GetTotPoints();
            int nTracePts       = fields[0]->GetTrace()->GetTotPoints();
            
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
                    physvals = fluxX1 + n*nLocalSolutionPts;
                        
                    interpolatedFlux_m[n] = Blas::Ddot(
                                                    nLocalSolutionPts,
                                                    m_Ixm->GetPtr(), 1,
                                                    physvals, 1);
                        
                    interpolatedFlux_p[n] = Blas::Ddot(
                                                    nLocalSolutionPts, 
                                                    m_Ixp->GetPtr(), 1, 
                                                    physvals, 1);   
                        
                    JumpL[0][n] = numericalFlux[n]   - interpolatedFlux_m[n];
                    JumpR[0][n] = numericalFlux[n+1] - interpolatedFlux_p[n];
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
                    
                    JumpL[0][n] = numericalFlux[n]   - fluxX1[offsetStart];
                    JumpR[0][n] = numericalFlux[n+1] - fluxX1[offsetEnd];
                }
            }
            
            for (n = 0; n < nElements; ++n)
            {
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                phys_offset       = fields[0]->GetPhys_Offset(n);

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
            
                // Assembling divergence of the correction flux
                Vmath::Vadd(nLocalSolutionPts, 
                            DCL, 1, 
                            DCR, 1, 
                            auxArray1 = divCFlux + phys_offset, 1);
            }
        }
        
        /**
         * @brief Compute the divergence of the corrective flux for 2D problems.
         *
         * @param nConvectiveFields   Number of fields (i.e. independent 
         *                            variables).
         * @param fields              Pointer to fields.
         * @param fluxX1              Volumetric flux in the physical space in 
         *                            direction X1.
         * @param fluxX2              Volumetric flux in the physical space in 
         *                            direction X2.
         * @param numericalFlux       Riemann flux in the physical space.
         * @param divCFlux            Divergence of the corrective flux for 2D
         *                            problems.
         *
         * \todo: Switch on shapes eventually here.
         */
        void AdvectionFR::v_DivCFlux_2D(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1, 
                const Array<OneD, const NekDouble> &fluxX2, 
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux)
        {                   
            int n, e, i, j, cnt;
            
            int nElements   = fields[0]->GetExpSize();
            int nDimensions = fields[0]->GetCoordim(0);  
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
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDimensions);
            for(i = 0; i < nDimensions; ++i)
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
                    
                    NekDouble fac = fields[0]->GetExp(n)->EdgeNormalNegated(e) ?
                        -1.0 : 1.0;

                    for (i = 0; i < nEdgePts; ++i)
                    {
                        if (m_traceNormals[0][trace_offset+i] != fac*normals[0][i] 
                        || m_traceNormals[1][trace_offset+i] != fac*normals[1][i])
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
        
        /**
         * @brief Compute the divergence of the corrective flux for 3D problems.
         *
         * @param nConvectiveFields   Number of fields (i.e. independent 
         *                            variables).
         * @param fields              Pointer to fields.
         * @param fluxX1              Volumetric flux in the physical space in 
         *                            direction X1.
         * @param fluxX2              Volumetric flux in the physical space in 
         *                            direction X2.
         * @param fluxX3              Volumetric flux in the physical space in 
         *                            direction X3.
         * @param numericalFlux       Riemann flux in the physical space.
         * @param divCFlux            Divergence of the corrective flux for 3D
         *                            Problems.
         *
         * \todo: To be implemented. Switch on shapes eventually here.
         */
        void AdvectionFR::v_DivCFlux_3D(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, const NekDouble> &fluxX1, 
            const Array<OneD, const NekDouble> &fluxX2,
            const Array<OneD, const NekDouble> &fluxX3, 
            const Array<OneD, const NekDouble> &numericalFlux,
                  Array<OneD,       NekDouble> &divCFlux)
        {

        }
        
    }
}
