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
#include "LibUtilities/Foundations/ManagerAccess.h"  // for PointsManager, etc
#include <StdRegions/StdSegExp.h>

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
            LibUtilities::SessionReaderSharedPtr              pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
        {
            LibUtilities::BasisSharedPtr Basis;
            LibUtilities::BasisSharedPtr BasisFR_Left;
            LibUtilities::BasisSharedPtr BasisFR_Right;
            Basis = pFields[0]->GetExp(0)->GetBasis(0);
            
            // Number of modes
            int nModes  = Basis->GetNumModes();
            
            // Total number of quadrature points
            int nSolutionPts = Basis->GetNumPoints();
            
            // Type of points
            const LibUtilities::PointsKey FRpoints = Basis->GetPointsKey();
            
            if (m_advType == "FRDG")
            {
                // Computation of the derivatives of the correction functions (DG)
                const LibUtilities::BasisKey  FRBase_Left (LibUtilities::eDG_DG_Left,  
                                                           nSolutionPts, 
                                                           FRpoints);
                    
                const LibUtilities::BasisKey  FRBase_Right(LibUtilities::eDG_DG_Right, 
                                                           nSolutionPts, 
                                                           FRpoints);
                
                BasisFR_Left  = LibUtilities::BasisManager()[FRBase_Left];
                BasisFR_Right = LibUtilities::BasisManager()[FRBase_Right];
                
                // Storing the derivatives into two global variables 
                m_dGL = BasisFR_Left ->GetBdata();
                m_dGR = BasisFR_Right->GetBdata();
            }
        
                    
            else if(m_advType == "FRSD")
            {
                // Computation of the derivatives of the correction functions (SD)
                const LibUtilities::BasisKey  FRBase_Left (LibUtilities::eDG_SD_Left,  
                                                           nSolutionPts, 
                                                           FRpoints);
                    
                const LibUtilities::BasisKey  FRBase_Right(LibUtilities::eDG_SD_Right, 
                                                           nSolutionPts, 
                                                           FRpoints);
                
                BasisFR_Left  = LibUtilities::BasisManager()[FRBase_Left];
                BasisFR_Right = LibUtilities::BasisManager()[FRBase_Right];
                
                // Storing the derivatives into two global variables 
                m_dGL = BasisFR_Left ->GetBdata();
                m_dGR = BasisFR_Right->GetBdata();
                    
            }
                    
            else if (m_advType == "FRHU")
            {
                // Computation of the derivatives of the correction functions (HU)
                const LibUtilities::BasisKey  FRBase_Left (LibUtilities::eDG_HU_Left,  
                                                           nSolutionPts, 
                                                           FRpoints);
                    
                const LibUtilities::BasisKey  FRBase_Right(LibUtilities::eDG_HU_Right, 
                                                           nSolutionPts, 
                                                           FRpoints);
                
                BasisFR_Left  = LibUtilities::BasisManager()[FRBase_Left];
                BasisFR_Right = LibUtilities::BasisManager()[FRBase_Right];
                
                // Storing the derivatives into two global variables 
                m_dGL = BasisFR_Left ->GetBdata();
                m_dGR = BasisFR_Right->GetBdata();
            }
        }
        
        void AdvectionFR::v_Advect(
            const int                                          nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            // Counter variables
            int i, j;
            
            // Number of elements
            int nElements       = fields[0]->GetExpSize();
            
            // Number of spatial dimensions
            int nDimensions     = fields[0]->GetCoordim(0);
                        
            // Number of solution points
            int nSolutionPts    = fields[0]->GetTotPoints();
            
            // Number of coefficients
            int nCoeffs         = fields[0]->GetNcoeffs();
            
            // Number of trace points
            int nTracePts       = fields[0]->GetTrace()->GetTotPoints();
                        
            // Vector to store the discontinuos flux
            Array<OneD, Array<OneD, NekDouble> > fluxvector(nDimensions);
            
            // Vector to store the derivative of the discontinuous flux
            Array<OneD, Array<OneD, NekDouble> > divfluxvector(nDimensions);
            
            // Vector to store the solution in physical space
            Array<OneD, Array<OneD, NekDouble> > physfield (nConvectiveFields);
            
            // Vector to store the flux interpolated at the interfaces
            Array<OneD, Array<OneD, NekDouble> > interpolatedFlux (nConvectiveFields);
                        
            // Resize each column of the flux vector to the number of 
            // solution points
            for(i = 0; i < nDimensions; ++i)
            {
                fluxvector[i]       = Array<OneD, NekDouble>(nSolutionPts);
                divfluxvector[i]    = Array<OneD, NekDouble>(nSolutionPts);
            }
            
            for(i = 0; i < nConvectiveFields; ++i)
            {
                physfield[i]            = inarray[i];
                interpolatedFlux[i]     = Array<OneD, NekDouble>(2*nTracePts-2);
            }
            
            // Get the discontinuous flux FD
            for(i = 0; i < nConvectiveFields; ++i)
            {                
                // Get the ith component of the flux vector in physical space
                m_fluxVector(i, physfield, fluxvector);
            }
                        
            Array<OneD,NekDouble> tmpFD, tmpDFD;
            
            // Computation of the divergence of the discontinuous flux at each 
            // solution point
            LibUtilities::BasisSharedPtr Basis;
            Basis = fields[0]->GetExp(0)->GetBasis(0);
            
            for(j = 0; j < nConvectiveFields; j++)
            {
                for(i = 0; i < nElements; i++)
                {
                    int phys_offset = fields[j]->GetPhys_Offset(i);
                    fields[j]->GetExp(i)->StdPhysDeriv(j, tmpFD   = fluxvector[j] + phys_offset, 
                                                          tmpDFD  = divfluxvector[j] + phys_offset);
                }
            }
            
            
            // Interpolation routine for Gauss points
            if (Basis->GetPointsType() == LibUtilities::eGaussGaussLegendre)
            {
                for (j = 0; j < nConvectiveFields; j++)
                {
                    v_InterpToInterface(fields, fluxvector[j], interpolatedFlux[j]);
                }
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
            

            // Computing the Riemann flux at each flux (interface) point
            m_riemann->Solve(Fwd, Bwd, numflux);
                        
            // Arrays to store the intercell numerical flux jumps
            Array<OneD, Array<OneD, NekDouble> > numfluxjumpsLeft(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > numfluxjumpsRight(nConvectiveFields);
            
            int offsetStart, offsetEnd;
            
            // Dimension of each column of the jumps array has to be equal to 
            // the number of flux points
            switch(nDimensions)
            {
                case 1:
                {
                    // Temporay array for computing the jumps multiply by the 
                    // derivatives of the correction functions
                    Array<OneD,NekDouble> dercorrfluxLeft(nSolutionPts/nElements,  0.0); 
                    Array<OneD,NekDouble> dercorrfluxRight(nSolutionPts/nElements, 0.0);
                    
                    Array<OneD,NekDouble> tmp, tmparray;
                    
                    // The dimension of each column of the jumps arrays is equal
                    // to number of trace points minus one
                    for(i = 0; i < nConvectiveFields; ++i)
                    {
                        numfluxjumpsLeft[i]  = Array<OneD, NekDouble>(nTracePts - 1);
                        numfluxjumpsRight[i] = Array<OneD, NekDouble>(nTracePts - 1);
                    }
                    
                    // Loop to compute the left and the right jump of the flux
                    for(i = 0; i < nElements; i++)
                    {
                        offsetStart              = fields[0]->GetPhys_Offset(i);
                        offsetEnd                = offsetStart + nSolutionPts/nElements - 1;
                        numfluxjumpsLeft[0][i]   = numflux[0][i] - fluxvector[0][offsetStart];
                        numfluxjumpsRight[0][i]  = numflux[0][i+1] - fluxvector[0][offsetEnd];
                    }
                    
                
                    // Loop to compute the left and the right jumps of the flux 
                    // with the new interpolation order used for Gauss points
                    if (Basis->GetPointsType() == LibUtilities::eGaussGaussLegendre)
                    {
                        j = 0;
                        for (i = 0; i < nElements; i++)
                        {
                            numfluxjumpsLeft[0][i]   = numflux[0][i] - interpolatedFlux[0][j];
                            j++;
                     
                            numfluxjumpsRight[0][i]  = numflux[0][i+1] - interpolatedFlux[0][j];
                            j++;
                        }
                    }
                    
                    for (i = 0; i < nElements; i++) 
                    {
                        Vmath::Smul(nSolutionPts/nElements, 
                                    numfluxjumpsLeft[0][i], 
                                    tmp = m_dGL, 1, 
                                    dercorrfluxLeft, 1);
                        
                        Vmath::Smul(nSolutionPts/nElements, 
                                    numfluxjumpsRight[0][i], 
                                    tmp = m_dGR, 1, 
                                    dercorrfluxRight, 1);
                        
                        Vmath::Vadd(nSolutionPts/nElements, 
                                    dercorrfluxLeft, 1, 
                                    dercorrfluxRight, 1, 
                                    tmparray = outarray[0] + i*nSolutionPts/nElements, 1);
                        
                        Vmath::Vadd(nSolutionPts/nElements, 
                                    tmparray = outarray[0] + i*nSolutionPts/nElements, 1, 
                                    tmp = divfluxvector[0] + i*nSolutionPts/nElements, 1, 
                                    tmparray = outarray[0] + i*nSolutionPts/nElements, 1); 
                    }
                    break;
                }
                case 2:
                {
                    ASSERTL0(false,"2D FRDG case not implemented yet");
                    break;
                }
                case 3:
                {
                    ASSERTL0(false,"3D FRDG case not implemented yet");
                    break;
                }
            }
            
            // Array to store the Jacobian and its inverse
            Array<OneD, const NekDouble>jac(nElements);
            Array<OneD, NekDouble>      jacobian(nElements);
            Array<OneD, NekDouble>      tmparray;
            
            // Evaluation of the jacobian of each element
            for(i = 0; i < nElements; i++)
            {
                jac         = fields[0]->GetExp(i)->GetGeom1D()->GetJac();
                jacobian[i] = jac[0];
            }
            
            // Operations to compute the RHS
            for(i = 0; i < nConvectiveFields; ++i)
            {
                for(j = 0; j < nElements; j++)
                {
                    Vmath::Smul(nSolutionPts/nElements, 1/jacobian[j], 
                                tmparray = outarray[i] + j*nSolutionPts/nElements, 1.0, 
                                tmparray = outarray[i] + j*nSolutionPts/nElements, 1.0);
                }
            }
        }
        
        
        
        void AdvectionFR::v_InterpToInterface(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            Array<OneD, NekDouble>    &total,
            Array<OneD, NekDouble>    &InterfaceValue)
        {            
            // Counter variables
            int n, i, j, cnt, numero;
            
            // Number of elements
            int nElements = fields[0]->GetExpSize();
            
            // Initial index of each element
            int phys_offset;
            
            // Number of quadrature points of each element
            int nLocalSolutionPts;
            
            // Number of vertices of each element
            int nLocalVertices;
            
            // Index of each interface between adjacent elements
            int interface_offset;
            
            cnt = 0;
            
            for (n = 0; n < nElements; n++)
            {
                // Getting the basis
                LibUtilities::BasisSharedPtr Basis;
                Basis = fields[0]->GetExp(n)->GetBasis(0);
                StdRegions::StdSegExp StdSeg(Basis->GetBasisKey());
                
                // Set the offset of each element
                phys_offset = fields[0]->GetPhys_Offset(n);
                
                // Set the number of solution points of each element
                nLocalSolutionPts = fields[0]->GetExp(n)->GetNumPoints(0);
                
                // Set the number of vertices of each element
                nLocalVertices = fields[0]->GetExp(n)->GetNverts();
                
                Array<OneD, NekDouble> tmp(nLocalSolutionPts, 0.0);
                tmp = total + n*nLocalSolutionPts;
                
                for (i = 0; i < nLocalVertices; i++)
                {
                    Array<OneD, NekDouble> interface_coord(3, 0.0);
                    interface_offset = fields[0]->GetExp(n)->GetGeom1D()->GetVid(i);
                    
                    for (j = 0; j < (fields[0]->GetExp(n)->GetVertexNormal(i)).num_elements(); j++)
                    {
                        interface_coord[0] += (fields[0]->GetExp(n)->GetVertexNormal(i))[j][0];
                        
                        //InterfaceValue[cnt] = fields[0]->GetExp(n)->PhysEvaluate(interface_coord, total);
                        InterfaceValue[cnt] = StdSeg.PhysEvaluate(interface_coord, tmp);
                        cnt++;
                    }
                }
            }
        }// End v_InterpToInterface
    }
}
