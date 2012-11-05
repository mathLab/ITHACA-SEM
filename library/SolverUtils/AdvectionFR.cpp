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

#include <SolverUtils/AdvectionFR.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <StdRegions/StdSegExp.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>



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
                                        "FRHU", AdvectionFR::create)};
        
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
            int nquad0, nquad1, nquad2;
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
            int i, j, n, p;
            int nquad0, nquad1, nquad2;
            int nmodes0, nmodes1, nmodes2;
            Array<OneD, NekDouble> auxArray1, auxArray2;
            LibUtilities::PointsKey FRPts0, FRPts1, FRPts2;
                        
            int nElements   = pFields[0]->GetExpSize();            
            int nDimensions = pFields[0]->GetCoordim(0);
            
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            
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
                        FRPts0    = base[0]->GetPointsKey();   
                        
                        LibUtilities::BasisSharedPtr BasisFR_Left0;
                        LibUtilities::BasisSharedPtr BasisFR_Right0;
                        
                        m_dGL_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGR_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        
                        if (m_advType == "FRDG")
                        {
                            // Derivatives of the correction functions (DG)
                            const LibUtilities::BasisKey FRBase_Left0(
                                LibUtilities::eDG_DG_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Right0(
                                LibUtilities::eDG_DG_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                        }
                        else if (m_advType == "FRSD")
                        {
                            // Derivatives of the correction functions (SD)
                            const LibUtilities::BasisKey FRBase_Left0(
                                LibUtilities::eDG_SD_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Right0(
                                LibUtilities::eDG_SD_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                        }
                        else if (m_advType == "FRHU")
                        {
                            // Derivatives of the correction functions (HU)
                            const LibUtilities::BasisKey  FRBase_Left0(
                                LibUtilities::eDG_HU_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey  FRBase_Right0(
                                LibUtilities::eDG_HU_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                        }
                        
                        // Storing the derivatives into two global variables 
                        m_dGL_xi1[n] = BasisFR_Left0 ->GetBdata();
                        m_dGR_xi1[n] = BasisFR_Right0->GetBdata();
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
                        FRPts0    = base[0]->GetPointsKey();   
                        FRPts1    = base[1]->GetPointsKey(); 
                        
                        m_dGL_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGR_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGL_xi2[n] = Array<OneD, NekDouble>(nquad1);
                        m_dGR_xi2[n] = Array<OneD, NekDouble>(nquad1);
                        
                        if (m_advType == "FRDG")
                        {
                            // Derivatives of the correction functions (DG)
                            const LibUtilities::BasisKey FRBase_Left0(
                                LibUtilities::eDG_DG_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Right0(
                                LibUtilities::eDG_DG_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Left1(
                                LibUtilities::eDG_DG_Left,  
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Right1(
                                LibUtilities::eDG_DG_Right, 
                                nmodes1, 
                                FRPts1);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = 
                                LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = 
                                LibUtilities::BasisManager()[FRBase_Right1];
                        }
                        else if (m_advType == "FRSD")
                        {
                            // Derivatives of the correction functions (SD)
                            const LibUtilities::BasisKey FRBase_Left0(
                                LibUtilities::eDG_SD_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Right0(
                                LibUtilities::eDG_SD_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Left1(
                                LibUtilities::eDG_SD_Left,  
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Right1(
                                LibUtilities::eDG_SD_Right, 
                                nmodes1, 
                                FRPts1);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = 
                                LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = 
                                LibUtilities::BasisManager()[FRBase_Right1];
                        }
                        else if (m_advType == "FRHU")
                        {
                            // Derivatives of the correction functions (HU)
                            const LibUtilities::BasisKey FRBase_Left0(
                                LibUtilities::eDG_HU_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Right0(
                                LibUtilities::eDG_HU_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Left1(
                                LibUtilities::eDG_HU_Left,  
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Right1(
                                LibUtilities::eDG_HU_Right, 
                                nmodes1, 
                                FRPts1);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = 
                                LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = 
                                LibUtilities::BasisManager()[FRBase_Right1];
                        }
                        
                        // Storing the derivatives into two global variables 
                        m_dGL_xi1[n] = BasisFR_Left0 ->GetBdata();
                        m_dGR_xi1[n] = BasisFR_Right0->GetBdata();
                        m_dGL_xi2[n] = BasisFR_Left1 ->GetBdata();
                        m_dGR_xi2[n] = BasisFR_Right1->GetBdata();
                    }
                    break;
                }
                case 3:
                {
                    LibUtilities::BasisSharedPtr BasisFR_Left0;
                    LibUtilities::BasisSharedPtr BasisFR_Right0;
                    LibUtilities::BasisSharedPtr BasisFR_Left1;
                    LibUtilities::BasisSharedPtr BasisFR_Right1;
                    LibUtilities::BasisSharedPtr BasisFR_Left2;
                    LibUtilities::BasisSharedPtr BasisFR_Right2;
                    
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
                        FRPts0    = base[0]->GetPointsKey();   
                        FRPts1    = base[1]->GetPointsKey(); 
                        FRPts2    = base[2]->GetPointsKey(); 
                        
                        m_dGL_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGR_xi1[n] = Array<OneD, NekDouble>(nquad0);
                        m_dGL_xi2[n] = Array<OneD, NekDouble>(nquad1);
                        m_dGR_xi2[n] = Array<OneD, NekDouble>(nquad1);
                        m_dGL_xi3[n] = Array<OneD, NekDouble>(nquad2);
                        m_dGR_xi3[n] = Array<OneD, NekDouble>(nquad2);
                        
                        if (m_advType == "FRDG")
                        {
                            // Derivatives of the correction functions (DG)
                            const LibUtilities::BasisKey FRBase_Left0(
                                LibUtilities::eDG_DG_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Right0(
                                LibUtilities::eDG_DG_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Left1(
                                LibUtilities::eDG_DG_Left,  
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Right1(
                                LibUtilities::eDG_DG_Right, 
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Left2(
                                LibUtilities::eDG_DG_Left,  
                                nmodes2, 
                                FRPts2);
                            
                            const LibUtilities::BasisKey FRBase_Right2(
                                LibUtilities::eDG_DG_Right, 
                                nmodes2, 
                                FRPts2);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = 
                                LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = 
                                LibUtilities::BasisManager()[FRBase_Right1];
                            BasisFR_Left2  = 
                                LibUtilities::BasisManager()[FRBase_Left2];
                            BasisFR_Right2 = 
                                LibUtilities::BasisManager()[FRBase_Right2];
                        }
                        else if (m_advType == "FRSD")
                        {
                            // Derivatives of the correction functions (SD)
                            const LibUtilities::BasisKey FRBase_Left0(
                                LibUtilities::eDG_SD_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Right0(
                                LibUtilities::eDG_SD_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Left1(
                                LibUtilities::eDG_SD_Left,  
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Right1(
                                LibUtilities::eDG_SD_Right, 
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Left2(
                                LibUtilities::eDG_SD_Left,  
                                nmodes2, 
                                FRPts2);
                            
                            const LibUtilities::BasisKey FRBase_Right2(
                                LibUtilities::eDG_SD_Right, 
                                nmodes2, 
                                FRPts2);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = 
                                LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = 
                                LibUtilities::BasisManager()[FRBase_Right1];
                            BasisFR_Left2  = 
                                LibUtilities::BasisManager()[FRBase_Left2];
                            BasisFR_Right2 = 
                                LibUtilities::BasisManager()[FRBase_Right2];
                        }
                        else if (m_advType == "FRHU")
                        {
                            // Derivatives of the correction functions (HU)
                            const LibUtilities::BasisKey FRBase_Left0(
                                LibUtilities::eDG_HU_Left,  
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Right0(
                                LibUtilities::eDG_HU_Right, 
                                nmodes0, 
                                FRPts0);
                            
                            const LibUtilities::BasisKey FRBase_Left1(
                                LibUtilities::eDG_HU_Left,  
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Right1(
                                LibUtilities::eDG_HU_Right, 
                                nmodes1, 
                                FRPts1);
                            
                            const LibUtilities::BasisKey FRBase_Left2(
                                LibUtilities::eDG_HU_Left,  
                                nmodes2, 
                                FRPts2);
                            
                            const LibUtilities::BasisKey FRBase_Right2(
                                LibUtilities::eDG_HU_Right, 
                                nmodes2, 
                                FRPts2);
                            
                            BasisFR_Left0  = 
                                LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = 
                                LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = 
                                LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = 
                                LibUtilities::BasisManager()[FRBase_Right1];
                            BasisFR_Left2  = 
                                LibUtilities::BasisManager()[FRBase_Left2];
                            BasisFR_Right2 = 
                                LibUtilities::BasisManager()[FRBase_Right2];
                        }
                        // Storing the derivatives into two global variables 
                        m_dGL_xi1[n] = BasisFR_Left0 ->GetBdata();
                        m_dGR_xi1[n] = BasisFR_Right0->GetBdata();
                        m_dGL_xi2[n] = BasisFR_Left1 ->GetBdata();
                        m_dGR_xi2[n] = BasisFR_Right1->GetBdata();
                        m_dGL_xi3[n] = BasisFR_Left2 ->GetBdata();
                        m_dGR_xi3[n] = BasisFR_Right2->GetBdata();
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
                int n, i, j, nLocalVertices;
                int nquad0, nquad1, nquad2;
                
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
            
            Array<OneD, Array<OneD, NekDouble> > fluxvector(nDimensions);            
            Array<OneD, Array<OneD, NekDouble> > Dfluxvector(nDimensions);            
            
            for(i = 0; i < nDimensions; ++i)
            {
                fluxvector[i]  = Array<OneD, NekDouble>(nSolutionPts);
                Dfluxvector[i] = Array<OneD, NekDouble>(nSolutionPts);
            }
            
            // Get the discontinuous flux FD ("i" is used by inarray)
            for(i = 0; i < nConvectiveFields; ++i)
            {                
                // Get the ith component of the flux vector in standard space
                m_fluxVector(i, inarray, fluxvector);
            }
            
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
                    // Divergence of the discontinuous flux
                    for (i = 0; i < nConvectiveFields; i++)
                    {
                        for (n = 0; n < nElements; n++)
                        {
                            phys_offset       = fields[0]->GetPhys_Offset(n);
                            
                            fields[i]->GetExp(n)->PhysDeriv(
                                0, auxArray1 = fluxvector[0] + phys_offset, 
                                auxArray2 = Dfluxvector[0] + phys_offset);
                        }
                    }
                    
                    // Divergence of the correction flux
                    Array<OneD, NekDouble> divFC(nSolutionPts, 0.0);
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        v_DivCFlux_1D(nConvectiveFields,
                                      fields, 
                                      fluxvector[0], 
                                      numflux[i], 
                                      divFC);
                    }

                    // Computation of the advection term
                    for (n = 0; n < nElements; n++) 
                    {
                        nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                        phys_offset = fields[0]->GetPhys_Offset(n);
                        jac = fields[0]->GetExp(n)->GetGeom1D()->GetJac();
                        
                        Vmath::Smul(nLocalSolutionPts, 
                                    1/jac[0], 
                                    auxArray1 = divFC + phys_offset, 1, 
                                    auxArray2 = outarray[0] + phys_offset, 1);
                        
                        Vmath::Vadd(nLocalSolutionPts, 
                                    auxArray1 = outarray[0] + phys_offset, 1, 
                                    auxArray2 = Dfluxvector[0] + phys_offset, 1, 
                                    auxArray3 = outarray[0] + phys_offset, 1); 
                    }
                    break;
                }
                // 2D-Problems 
                case 2:
                { 
                    // Divergence of the discontinuous flux
                    for (n = 0; n < nElements; ++n)
                    {
                        nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                        phys_offset = fields[0]->GetPhys_Offset(n);
                        
                        jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                        gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                        
                        Array<OneD, NekDouble> f_hat(nLocalSolutionPts, 0.0);
                        Array<OneD, NekDouble> g_hat(nLocalSolutionPts, 0.0);
                        
                        if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype()
                            == SpatialDomains::eDeformed)
                        {
                            for (j = 0; j < nLocalSolutionPts; j++)
                            {
                                f_hat[j] = 
                                (fluxvector[0][j+phys_offset]*gmat[0][j] + 
                                 fluxvector[1][j+phys_offset]*gmat[2][j])*jac[j];
                                
                                g_hat[j] = 
                                (fluxvector[0][j+phys_offset]*gmat[1][j] + 
                                 fluxvector[1][j+phys_offset]*gmat[3][j])*jac[j];
                            }
                        }
                        else
                        {
                            for (j = 0; j < nLocalSolutionPts; j++)
                            {
                                f_hat[j] = 
                                (fluxvector[0][j+phys_offset]*gmat[0][0] + 
                                 fluxvector[1][j+phys_offset]*gmat[2][0])*jac[0];
                                
                                g_hat[j] = 
                                (fluxvector[0][j+phys_offset]*gmat[1][0] + 
                                 fluxvector[1][j+phys_offset]*gmat[3][0])*jac[0];
                            }
                        }
                        
                        fields[0]->GetExp(n)->StdPhysDeriv(0, 
                                    auxArray1 = f_hat, 
                                    auxArray2 = Dfluxvector[0] + phys_offset); 
                        
                        fields[0]->GetExp(n)->StdPhysDeriv(1, 
                                    auxArray1 = g_hat, 
                                    auxArray2 = Dfluxvector[1] + phys_offset); 
                    }
                                                    
                    // Computation of the divergence of the discontinuous flux
                    Array<OneD, NekDouble> divFD(nSolutionPts, 0.0);
                    Vmath::Vadd(nSolutionPts, 
                                Dfluxvector[0], 1,
                                Dfluxvector[1], 1, 
                                auxArray1 = divFD, 1);
                    
                    // Computation of the divergence of the correction flux
                    Array<OneD, NekDouble> divFC(nSolutionPts, 0.0);
                    for (j = 0; j < nConvectiveFields; ++j)
                    {
                        v_DivCFlux_2D(nConvectiveFields,
                                      fields, 
                                      fluxvector[0], 
                                      fluxvector[1],
                                      numflux[j], 
                                      divFC);
                    }     
                        
                    // Computation of the divergence of the final flux
                    Vmath::Vadd(nSolutionPts, 
                                divFD, 1, 
                                divFC, 1, 
                                auxArray1 = outarray[0], 1);
                    
                    // Multiplication by the inverse of the jacobian
                    for (n = 0; n < nElements; ++n)
                    {
                        nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                        phys_offset = fields[0]->GetPhys_Offset(n);
                        
                        jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                        gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                        
                        if (fields[0]->GetExp(n)->GetGeom2D()->GetGtype() 
                            == SpatialDomains::eDeformed)
                        {
                            for (i = 0; i < nLocalSolutionPts; i++)
                            {
                                outarray[0][phys_offset + i] = 
                                outarray[0][phys_offset + i]/jac[i];
                            }
                        }
                        else
                        {
                            Vmath::Smul(
                                nLocalSolutionPts, 
                                1/jac[0], 
                                auxArray1 = outarray[0] + phys_offset, 1, 
                                auxArray2 = outarray[0] + phys_offset, 1);
                        }
                    }
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
         *                            Problems.
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
