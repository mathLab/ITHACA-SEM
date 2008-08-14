///////////////////////////////////////////////////////////////////////////////
//
// File Expansion.cpp
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
// Description: File for Expansion routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        Expansion::Expansion(void)
        {
        }
        
        // LocalRegions matrix generation which can be done at generic level
        DNekMatSharedPtr Expansion::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            
            DNekMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
                {
                    
                    ASSERTL1(v_IsBoundaryInteriorExpansion(),
                             "HybridDGHelmholtz matrix not set up "
                             "for non boundary-interior expansions");
                    
                    int i,j;
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    int       ncoeffs = v_GetNcoeffs();
                    
                    // Get basic Galerkin Helmholtz matrix 
                    DNekScalMat &Hmat = *v_GetLocMatrix(StdRegions::eHelmholtz,lambdaval);
                                        
                    int rows = Hmat.GetRows();
                    int cols = Hmat.GetColumns();
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                    DNekMat &Mat = *returnval;

                    //must be a better way of copying matrix
                    for(i = 0; i < cols; ++i)
                    {
                        for(j = 0; j < rows; ++j)
                        {
                            Mat(i,j) = Hmat(i,j);
                        }
                    }

                    //Mat = *Hmat.GetOwnedMatrix(); // might need to get hold of matrix and constant 
                    
                    Array<OneD,NekDouble> inarray(ncoeffs);
                    Array<OneD,NekDouble> outarray(ncoeffs);
                    
                    for(j = 0; j < ncoeffs; ++j)
                    {
                        Vmath::Zero(ncoeffs,&inarray[0],1);
                        Vmath::Zero(ncoeffs,&outarray[0],1);
                        inarray[j] = 1.0;
                        
                        v_AddHDGHelmholtzMatrixBoundaryTerms(tau,inarray,outarray);
                        
                        for(i = 0; i < ncoeffs; ++i)
                        {
                            Mat(i,j) += outarray[i];
                        }
                    }
                }
                break;
            case StdRegions::eHybridDGLamToU:
                {
                    int i,j,k;
                    int nbndry = v_NumDGBndryCoeffs();
                    int ncoeffs = v_GetNcoeffs();
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    
                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);                    
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    
                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Umat = *returnval;
                    
                    // Helmholtz matrix
                    DNekScalMat  &invHmat = *v_GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, lambdaval, tau);
                    
                    // for each degree of freedom of the lambda space
                    // calculate Umat entry 
                    // Generate Lambda to U_lambda matrix 
                    for(j = 0; j < nbndry; ++j)
                    {
                        Vmath::Zero(nbndry,&lambda[0],1);
                        Vmath::Zero(ncoeffs,&f[0],1);
                        lambda[j] = 1.0;
                        
                        v_SetTraceToGeomOrientation(lambda);
                        
                        v_AddHDGHelmholtzTraceTerms(tau,lambda,f);
                        
                        Ulam = invHmat*F; // generate Ulam from lambda
                        
                        // fill column of matrix
                        for(k = 0; k < ncoeffs; ++k)
                        {
                            Umat(k,j) = Ulam[k]; 
                        }
                    }
                }
                break;
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
                {
                    int i,j,k,dir;
                    int nbndry = v_NumDGBndryCoeffs();
                    int nquad  = v_GetNumPoints(0);
                    int ncoeffs = v_GetNcoeffs();

                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);                    
                    
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    
                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Qmat = *returnval;
                    
                    // Helmholtz matrix
                    DNekScalMat &invHmat = *v_GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, lambdaval,tau);
                    
                    // Lambda to U matrix
                    DNekScalMat &lamToU = *v_GetLocMatrix(StdRegions::eHybridDGLamToU, lambdaval, tau);
                    
                    // Inverse mass matrix 
                    DNekScalMat &invMass = *v_GetLocMatrix(StdRegions::eInvMass);
                    
                    //Weak Derivative matrix 
                    DNekScalMatSharedPtr Dmat;
                    switch(mkey.GetMatrixType())
                    {
                    case StdRegions::eHybridDGLamToQ0:
                        dir = 0;
                        Dmat = v_GetLocMatrix(StdRegions::eWeakDeriv0); 
                        break;
                    case StdRegions::eHybridDGLamToQ1:
                        dir = 1;
                        Dmat = v_GetLocMatrix(StdRegions::eWeakDeriv1); 
                        break;
                    case StdRegions::eHybridDGLamToQ2:
                        dir = 2;
                        Dmat = v_GetLocMatrix(StdRegions::eWeakDeriv2); 
                        break;
                    }
                
                    // for each degree of freedom of the lambda space
                    // calculate Qmat entry 
                    // Generate Lambda to Q_lambda matrix 
                    for(j = 0; j < nbndry; ++j)
                    {
                        Vmath::Zero(nbndry,&lambda[0],1);
                        lambda[j] = 1.0;
                        
                        // for lambda[j] = 1 this is the solution to ulam
                        for(k = 0; k < ncoeffs; ++k)
                        {
                            Ulam[k] = lamToU(k,j);
                        }
                        
                        // -D^T ulam
                        Vmath::Neg(ncoeffs,&ulam[0],1);
                        F = Transpose(*Dmat)*Ulam; 
                        
                        v_SetTraceToGeomOrientation(lambda);
                        
                        // + \tilde{G} \lambda
                        v_AddNormTraceInt(dir,lambda,f); 
                        
                        // multiply by inverse mass matrix
                        Ulam = invMass*F; 
                        
                        // fill column of matrix (Qmat is in column major format)
                        Vmath::Vcopy(ncoeffs,&ulam[0],1,&(Qmat.GetPtr())[0]+j*ncoeffs,1);
                    }
                }
                break;            
            default:
                ASSERTL0(false,"This matrix type cannot be generated from this class");
                break;
            }            
            return returnval;
        }            
    } //end of namespace
} //end of namespace

/** 
 *    $Log: Expansion.cpp,v $
 *
 **/
