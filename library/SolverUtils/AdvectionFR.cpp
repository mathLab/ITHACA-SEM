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
        
        AdvectionFR::AdvectionFR(std::string advType):m_advType(advType)
        {
        }
        
        void AdvectionFR::v_InitObject(
                LibUtilities::SessionReaderSharedPtr        pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            int i, j, n;
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
                    LibUtilities::BasisSharedPtr BasisFR_Left0;
                    LibUtilities::BasisSharedPtr BasisFR_Right0;
                    
                    m_dGL_xi1 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi1 = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    
                    for (n = 0; n < nElements; ++n)
                    {
                        base      = pFields[0]->GetExp(n)->GetBase();
                        nquad0    = base[0]->GetNumPoints();
                        nmodes0   = base[0]->GetNumModes();
                        FRPts0    = base[0]->GetPointsKey();
                        
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
                            
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
                        }
                        
                        
                        else if(m_advType == "FRSD")
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
                        
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
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
                            
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
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
                    
                    m_dGL_xi1     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi1     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGL_xi2     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi2     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    
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
                            
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = LibUtilities::BasisManager()[FRBase_Right1];
                        }
                        
                        
                        else if(m_advType == "FRSD")
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
                            
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = LibUtilities::BasisManager()[FRBase_Right1];
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
                            
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = LibUtilities::BasisManager()[FRBase_Right1];
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
                    
                    m_dGL_xi1     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi1     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGL_xi2     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi2     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGL_xi3     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    m_dGR_xi3     = Array<OneD, Array<OneD, NekDouble> >(nElements);
                    
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
                            
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = LibUtilities::BasisManager()[FRBase_Right1];
                            BasisFR_Left2  = LibUtilities::BasisManager()[FRBase_Left2];
                            BasisFR_Right2 = LibUtilities::BasisManager()[FRBase_Right2];
                        }
                        
                        
                        else if(m_advType == "FRSD")
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
                            
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = LibUtilities::BasisManager()[FRBase_Right1];
                            BasisFR_Left2  = LibUtilities::BasisManager()[FRBase_Left2];
                            BasisFR_Right2 = LibUtilities::BasisManager()[FRBase_Right2];
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
                            
                            BasisFR_Left0  = LibUtilities::BasisManager()[FRBase_Left0];
                            BasisFR_Right0 = LibUtilities::BasisManager()[FRBase_Right0];
                            BasisFR_Left1  = LibUtilities::BasisManager()[FRBase_Left1];
                            BasisFR_Right1 = LibUtilities::BasisManager()[FRBase_Right1];
                            BasisFR_Left2  = LibUtilities::BasisManager()[FRBase_Left2];
                            BasisFR_Right2 = LibUtilities::BasisManager()[FRBase_Right2];
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
            
            LibUtilities::BasisSharedPtr Basis;
            Basis = pFields[0]->GetExp(0)->GetBasis(0);
            
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
        
        void AdvectionFR::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int i, j, n, yepp;
            int nLocalSolutionPts, phys_offset;
            int offsetStart, offsetEnd;
            
            Array<OneD,       NekDouble> auxArray1, auxArray2, auxArray3;
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            
            LibUtilities::BasisSharedPtr Basis;
            Basis = fields[0]->GetExp(0)->GetBasis(0);
            
            int nElements       = fields[0]->GetExpSize();            
            int nDimensions     = fields[0]->GetCoordim(0);                        
            int nSolutionPts    = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
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
            
            // Divergence of the correction and final fluxes  
            switch(nDimensions)
            {
                case 1:
                {                    
                    // Divergence of the discontinuous flux at each solution point
                    for (i = 0; i < nConvectiveFields; i++)
                    {
                        for (n = 0; n < nElements; n++)
                        {
                            gmat = fields[0]->GetExp(n)->GetGeom1D()->GetGmat();
                            nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                            phys_offset       = fields[i]->GetPhys_Offset(n);
                            
                            fields[i]->GetExp(n)->PhysDeriv(
                                0, auxArray1 = fluxvector[0] + phys_offset, 
                                auxArray2 = Dfluxvector[0] + phys_offset);
                        }
                    }
                    
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
                    
                    if (Basis->GetPointsType() == LibUtilities::eGaussGaussLegendre)
                    {
                        // Interpolation routine for Gauss points
                        Array<OneD, NekDouble> interpolatedFlux_m(nElements);
                        Array<OneD, NekDouble> interpolatedFlux_p(nElements);
                        
                        for (j = 0; j < nConvectiveFields; ++j)
                        {
                            for (n = 0; n < nElements; ++n)
                            {
                                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                                
                                Array<OneD, NekDouble> physvals(nLocalSolutionPts, 0.0);
                                physvals = fluxvector[0] + n*nLocalSolutionPts;
                                
                                interpolatedFlux_m[n] = Blas::Ddot(
                                                            nLocalSolutionPts,
                                                            m_Ixm->GetPtr(), 1,
                                                            physvals, 1);
                                
                                interpolatedFlux_p[n] = Blas::Ddot(
                                                            nLocalSolutionPts, 
                                                            m_Ixp->GetPtr(), 1, 
                                                            physvals, 1);   
                                
                                JumpL[j][n] = numflux[j][n]   - interpolatedFlux_m[n];
                                JumpR[j][n] = numflux[j][n+1] - interpolatedFlux_p[n];
                            }
                        }
                    }
                    else
                    {
                        for(n = 0; n < nElements; ++n)
                        {
                            nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                            
                            offsetStart = fields[0]->GetPhys_Offset(n);
                            offsetEnd   = offsetStart + nLocalSolutionPts - 1;
                            
                            JumpL[0][n] = numflux[0][n]   - fluxvector[0][offsetStart];
                            JumpR[0][n] = numflux[0][n+1] - fluxvector[0][offsetEnd];
                        }
                    }
                    
                    for (n = 0; n < nElements; n++) 
                    {
                        nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                        phys_offset = fields[0]->GetPhys_Offset(n);
                        
                        gmat = fields[0]->GetExp(n)->GetGeom1D()->GetGmat();
                        jac = fields[0]->GetExp(n)->GetGeom1D()->GetJac();

                        Vmath::Smul(nLocalSolutionPts, 
                                    JumpL[0][n], 
                                    auxArray1 = m_dGL_xi1[n], 1, 
                                    DCL, 1);
                        
                        Vmath::Smul(nLocalSolutionPts, 
                                    JumpR[0][n], 
                                    auxArray1 = m_dGR_xi1[n], 1, 
                                    DCR, 1);
                        
                        Vmath::Vadd(nLocalSolutionPts, 
                                    DCL, 1, 
                                    DCR, 1, 
                                    auxArray1 = outarray[0] + phys_offset, 1);
                        
                        Vmath::Smul(nLocalSolutionPts, 
                                    1/jac[0], 
                                    auxArray1 = outarray[0] + phys_offset, 1, 
                                    auxArray2 = outarray[0] + phys_offset, 1);
                        
                        Vmath::Vadd(nLocalSolutionPts, 
                                    auxArray1 = outarray[0] + phys_offset, 1, 
                                    auxArray2 = Dfluxvector[0] + phys_offset, 1, 
                                    auxArray3 = outarray[0] + phys_offset, 1); 
                    }
                    break;
                }
                case 2:
                {                    
                    // Derivatives of the discontinuous flux at each solution point
                    for (i = 0; i < nConvectiveFields; ++i)
                    {
                        for (j = 0; j < nDimensions; ++j)
                        {
                            for (n = 0; n < nElements; ++n)
                            {
                                phys_offset = fields[i]->GetPhys_Offset(n);
                                
                                fields[i]->GetExp(n)->PhysDeriv(
                                    j, auxArray1 = fluxvector[j] + phys_offset, 
                                    auxArray2 = Dfluxvector[j] + phys_offset);
                            }
                        }
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
                        v_divCorrFlux(fields, 
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
                    break;
                }
                case 3:
                {
                    ASSERTL0(false,"3D FRDG case not implemented yet");
                    break;
                }
            }
        }
        
        
        
        void AdvectionFR::v_divCorrFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX, 
                const Array<OneD, const NekDouble> &fluxY, 
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux)
        {                   
            int n, e, i, j, cnt; 
            int nElements = fields[0]->GetExpSize();
            int nLocalSolutionPts;
            int nEdgePts;  
            int trace_offset; 
            int phys_offset;
            int nquad0;
            int nquad1;
            
            Array<OneD, NekDouble> auxArray1, auxArray2;
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
            &elmtToTrace = fields[0]->GetTraceMap()->GetElmtToTrace();
            
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            
            // Loop on the elements
            for(n = 0; n < nElements; ++n)
            {
                // Offset of the element on the global vector
                phys_offset = fields[0]->GetPhys_Offset(n);
                nLocalSolutionPts = fields[0]->GetExp(n)->GetTotPoints();
                
                base = fields[0]->GetExp(n)->GetBase();
                nquad0 = base[0]->GetNumPoints();
                nquad1 = base[1]->GetNumPoints();
                
                Array<OneD, NekDouble> divCFluxE0(nquad0 * nquad1, 0.0);
                Array<OneD, NekDouble> divCFluxE1(nquad0 * nquad1, 0.0);
                Array<OneD, NekDouble> divCFluxE2(nquad0 * nquad1, 0.0);
                Array<OneD, NekDouble> divCFluxE3(nquad0 * nquad1, 0.0);
                
                Array<OneD, NekDouble> derGE0(nquad1, 0.0);
                Array<OneD, NekDouble> derGE1(nquad0, 0.0);
                Array<OneD, NekDouble> derGE2(nquad1, 0.0);
                Array<OneD, NekDouble> derGE3(nquad0, 0.0);
                
                // Get the metric terms
                gmat = fields[0]->GetExp(n)->GetGeom2D()->GetGmat();
                jac  = fields[0]->GetExp(n)->GetGeom2D()->GetJac();
                
                // Loop on the edges
                for(e = 0; e < fields[0]->GetExp(n)->GetNedges(); ++e)
                {   
                    // Number of edge points of edge e
                    nEdgePts = fields[0]->GetExp(n)->GetEdgeNumPoints(e);
                    
                    Array<OneD, NekDouble> tmparrayX(nEdgePts, 0.0);
                    Array<OneD, NekDouble> tmparrayY(nEdgePts, 0.0);
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
                                                fluxX + phys_offset,
                                                auxArray1 = tmparrayX);
                    
                    // Extract the edge values of flux-y on edge e and order 
                    // them accordingly to the order of the trace space
                    fields[0]->GetExp(n)->GetEdgePhysVals(
                                                e, elmtToTrace[n][e],
                                                fluxY + phys_offset,
                                                auxArray1 = tmparrayY);
                    
                    // Multiply the edge components of the flux by the normal
                    for (i = 0; i < normals[0].num_elements(); ++i)
                    {
                        fluxN[i] = tmparrayX[i]*normals[0][i] + 
                        tmparrayY[i]*normals[1][i];
                        
                        fluxT[i] = -tmparrayX[i]*normals[1][i] + 
                        tmparrayY[i]*normals[0][i];
                    }
                    
                    // Subtract to the Riemann flux the discontinuous flux 
                    Vmath::Vsub(nEdgePts, 
                                &numericalFlux[trace_offset], 1, 
                                &fluxN[0], 1, &fluxJumps[0], 1);
                    
                    // Check the ordering of the jumps vector
                    if(fields[0]->GetExp(n)->GetEorient(e) == StdRegions::eBackwards)
                    {
                        Vmath::Reverse(nquad0, 
                                       auxArray1 = fluxJumps, 1,
                                       auxArray2 = fluxJumps, 1);
                    }
                    
                    // Multiply jumps by derivative of the correction functions
                    switch (e) 
                    {
                        case 0:
                            for (i = 0; i < nquad0; ++i)
                            {
                                for (j = 0; j < nquad1; ++j)
                                {
                                    cnt = i + j*nquad0;
                                    divCFluxE0[cnt] = -fluxJumps[i] * m_dGL_xi2[n][j];
                                }
                            }
                            break;
                        case 1:
                            for (i = 0; i < nquad1; ++i)
                            {
                                for (j = 0; j < nquad0; ++j)
                                {
                                    cnt = (nquad0)*i + nquad0-1 - j;
                                    divCFluxE1[cnt] = fluxJumps[i] * m_dGR_xi1[n][j];
                                }
                            }
                            break;
                        case 2:
                            for (i = 0; i < nquad0; ++i)
                            {
                                for (j = 0; j < nquad1; ++j)
                                {
                                    cnt = (nquad0*nquad1 - 1) - j*nquad0 - i;
                                    divCFluxE2[cnt] = fluxJumps[i] * m_dGR_xi2[n][j];
                                }
                            }
                            break;
                        case 3:  
                            for (i = 0; i < nquad1; ++i)
                            {
                                for (j = 0; j < nquad0; ++j)
                                {
                                    cnt = (nquad0*nquad1 - nquad0) + j - i*nquad0;
                                    divCFluxE3[cnt] = -fluxJumps[i] * m_dGL_xi1[n][j];
                                }
                            }
                            break;
                            
                        default:
                            ASSERTL0(false,"edge value (< 3) is out of range");
                            break;
                    }
                }
                
                // Multiply each edge contribution by the proper metrics
                Vmath::Smul(nLocalSolutionPts, 
                            1/jac[0], 
                            auxArray1 = divCFluxE0, 1, 
                            auxArray2 = divCFluxE0, 1);
                
                Vmath::Smul(nLocalSolutionPts, 
                            1/jac[0], 
                            auxArray1 = divCFluxE1, 1, 
                            auxArray2 = divCFluxE1, 1);
                
                Vmath::Smul(nLocalSolutionPts, 
                            1/jac[0], 
                            auxArray1 = divCFluxE2, 1, 
                            auxArray2 = divCFluxE2, 1);
                
                Vmath::Smul(nLocalSolutionPts, 
                            1/jac[0], 
                            auxArray1 = divCFluxE3, 1, 
                            auxArray2 = divCFluxE3, 1);

                // Sum each edge contribution
                for (i = 0; i < nquad0 * nquad1; ++i)
                {
                    divCFlux[phys_offset + i] = divCFluxE0[i] + divCFluxE1[i] +
                    divCFluxE2[i] + divCFluxE3[i];
                }
            }
        }
        
    }
}
