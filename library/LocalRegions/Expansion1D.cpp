///////////////////////////////////////////////////////////////////////////////
//
// File Expansion1D.cpp
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
// Description: File for Expansion1D routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion1D.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        
        DNekMatSharedPtr Expansion1D::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            DNekMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
                {                    
                ASSERTL1(IsBoundaryInteriorExpansion(),
                         "HybridDGHelmholtz matrix not set up "
                         "for non boundary-interior expansions");
                    int       i;
                    NekDouble lambdaval = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    NekDouble tau       = mkey.GetConstFactor(StdRegions::eFactorTau);
                    int       ncoeffs   = GetNcoeffs();

                    int       coordim = GetCoordim();

                    DNekScalMat  &invMass = *GetLocMatrix(StdRegions::eInvMass);
                    StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                           StdRegions::eWeakDeriv1,
                                                           StdRegions::eWeakDeriv2};
                    DNekMat LocMat(ncoeffs,ncoeffs); 

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,ncoeffs);
                    DNekMat &Mat = *returnval;

                    Vmath::Zero(ncoeffs*ncoeffs,Mat.GetPtr(),1);

                    for(i=0;  i < coordim; ++i)
                    {
                        DNekScalMat &Dmat = *GetLocMatrix(DerivType[i]);

                        Mat = Mat + Dmat*invMass*Transpose(Dmat);
                    }

                    // Add end Mass Matrix Contribution
                    DNekScalMat  &Mass = *GetLocMatrix(StdRegions::eMass);
                    Mat = Mat + lambdaval*Mass;                    

                    Array<OneD,unsigned int> bmap;
                    GetBoundaryMap(bmap);

                    // Add tau*F_e using elemental mass matrices
                    for(i = 0; i < 2; ++i)
                    {
                        Mat(bmap[i],bmap[i]) = Mat(bmap[i],bmap[i]) +  tau; 
                    }
                }
                break;
            case StdRegions::eHybridDGLamToU:
                {
                    int j,k;
                    int nbndry = NumDGBndryCoeffs();
                    int ncoeffs = GetNcoeffs();
                    StdRegions::ConstFactorMap factors;
                    factors[StdRegions::eFactorLambda] = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    factors[StdRegions::eFactorTau] = mkey.GetConstFactor(StdRegions::eFactorTau);
                    
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
                    DNekScalMat  &invHmat = *GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, factors);
                    
                    // for each degree of freedom of the lambda space
                    // calculate Umat entry 
                    // Generate Lambda to U_lambda matrix 
                    for(j = 0; j < nbndry; ++j)
                    {
                        Vmath::Zero(nbndry,&lambda[0],1);
                        Vmath::Zero(ncoeffs,&f[0],1);
                        lambda[j] = 1.0;
                        
                        AddHDGHelmholtzTraceTerms(factors[StdRegions::eFactorTau],lambda,f);
                        
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
                    int j,k,dir;
                    int nbndry = NumDGBndryCoeffs();
                    int ncoeffs = GetNcoeffs();

                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);                    
                    
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    StdRegions::ConstFactorMap factors;
                    factors[StdRegions::eFactorLambda] = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    factors[StdRegions::eFactorTau] = mkey.GetConstFactor(StdRegions::eFactorTau);

                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Qmat = *returnval;
                    
                    // Lambda to U matrix
                    DNekScalMat &lamToU = *GetLocMatrix(StdRegions::eHybridDGLamToU, factors);
                    
                    // Inverse mass matrix 
                    DNekScalMat &invMass = *GetLocMatrix(StdRegions::eInvMass);
                    
                    //Weak Derivative matrix 
                    DNekScalMatSharedPtr Dmat;
                    switch(mkey.GetMatrixType())
                    {
                    case StdRegions::eHybridDGLamToQ0:
                        dir = 0;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv0);
                        break;
                    case StdRegions::eHybridDGLamToQ1:
                        dir = 1;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv1);
                        break;
                    case StdRegions::eHybridDGLamToQ2:
                        dir = 2;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv2);
                        break;
                    default:
                        ASSERTL0(false,"Direction not known");
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
                        
                        // + \tilde{G} \lambda
                        AddNormTraceInt(dir,lambda,f); 
                        
                        // multiply by inverse mass matrix
                        Ulam = invMass*F; 
                        
                        // fill column of matrix (Qmat is in column major format)
                        Vmath::Vcopy(ncoeffs,&ulam[0],1,&(Qmat.GetPtr())[0]+j*ncoeffs,1);
                    }
                }
                break;            
            case StdRegions::eHybridDGHelmBndLam:
                {
                    int j;
                    int nbndry = NumBndryCoeffs();

                    StdRegions::ConstFactorMap factors;
                    factors[StdRegions::eFactorLambda] = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    factors[StdRegions::eFactorTau] = mkey.GetConstFactor(StdRegions::eFactorTau);

                    Array<OneD,unsigned int> bmap;
                    Array<OneD, NekDouble>   lam(2);
                    GetBoundaryMap(bmap);
                    
                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,  nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    // Matrix to map Lambda to U
                    DNekScalMat &LamToU = *GetLocMatrix(StdRegions::eHybridDGLamToU, factors);
                
                    // Matrix to map Lambda to Q
                    DNekScalMat &LamToQ = *GetLocMatrix(StdRegions::eHybridDGLamToQ0, factors);

                    lam[0] = 1.0; lam[1] = 0.0;
                    for(j = 0; j < nbndry; ++j)
                    {
                        BndMat(0,j) = -LamToQ(bmap[0],j) - factors[StdRegions::eFactorTau]*(LamToU(bmap[0],j) - lam[j]);
                    }

                    lam[0] = 0.0; lam[1] = 1.0;
                    for(j = 0; j < nbndry; ++j)
                    {
                        BndMat(1,j) =  LamToQ(bmap[1],j) - factors[StdRegions::eFactorTau]*(LamToU(bmap[1],j) - lam[j]);
                    }
                }
                break;
            default:
                ASSERTL0(false,"This matrix type cannot be generated from this class");
                break;
            }
            
            return returnval;
        }
        
        void Expansion1D::AddNormTraceInt(const int dir, 
                                          Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray) 
        {
            
            int k;
            int nbndry = NumBndryCoeffs();
            int nquad  = GetNumPoints(0);
            const Array<OneD, const NekDouble> &Basis  = GetBasis(0)->GetBdata();
            Array<OneD, unsigned int> vmap;
            
            GetBoundaryMap(vmap);
            
            // add G \lambda term (can assume G is diagonal since one
            // of the basis is zero at boundary otherwise)
            for(k = 0; k < nbndry; ++k)
            {
                outarray[vmap[k]] += (Basis[(vmap[k]+1)*nquad-1]*Basis[(vmap[k]+1)*nquad-1] - Basis[vmap[k]*nquad]*Basis[vmap[k]*nquad])*inarray[vmap[k]];
            }
        }

        void Expansion1D::AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                                 const Array<OneD,const NekDouble> &inarray,  Array<OneD,NekDouble> &outarray)
        {
            int i,n;
            int nbndry  = NumBndryCoeffs();
            int nquad   = GetNumPoints(0);
            int ncoeffs = GetNcoeffs();
            int coordim = GetCoordim();
            Array<OneD, unsigned int> vmap;
            
            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");


            const Array<OneD, const NekDouble> &Basis  = GetBasis(0)->GetBdata();
            DNekScalMat  &invMass = *GetLocMatrix(StdRegions::eInvMass);

            GetBoundaryMap(vmap);

            // Add F = \tau <phi_i,phi_j> (note phi_i is zero if phi_j is non-zero)
            for(i = 0; i < nbndry; ++i)
            {
                outarray[vmap[i]] += tau*Basis[(vmap[i]+1)*nquad-1]*Basis[(vmap[i]+1)*nquad-1]*inarray[vmap[i]];
                outarray[vmap[i]] += tau*Basis[vmap[i]*nquad]*Basis[vmap[i]*nquad]*inarray[vmap[i]];
            }
            

            //===============================================================
            // Add -\sum_i D_i^T M^{-1} G_i + E_i M^{-1} G_i = 
            //                         \sum_i D_i M^{-1} G_i term

            StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                   StdRegions::eWeakDeriv1,
                                                   StdRegions::eWeakDeriv2};
            Array<OneD, NekDouble> tmpcoeff(ncoeffs,0.0);
            DNekVec                Coeffs  (ncoeffs,outarray,eWrapper);
            DNekVec                Tmpcoeff(ncoeffs,tmpcoeff,eWrapper);

            for(n = 0; n < coordim; ++n)
            {
                // evaluate M^{-1} G
                for(i = 0; i < ncoeffs; ++i)
                {
                    // lower boundary (negative normal) 
                    tmpcoeff[i] -= invMass(i,vmap[0])*Basis[vmap[0]*nquad]*Basis[vmap[0]*nquad]*inarray[vmap[0]];
                    
                    // upper boundary (positive normal) 
                    tmpcoeff[i] += invMass(i,vmap[1])*Basis[(vmap[1]+1)*nquad-1]*Basis[(vmap[1]+1)*nquad-1]*inarray[vmap[1]];
                    
                }
                
                DNekScalMat &Dmat = *GetLocMatrix(DerivType[n]);
                Coeffs = Coeffs  + Dmat*Tmpcoeff; 
            }
        }

        void Expansion1D::v_AddRobinMassMatrix(const int vert, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat)
        {
            ASSERTL0(IsBoundaryInteriorExpansion(),"Robin boundary conditions are only implemented for boundary-interior expanisons");
            ASSERTL1(inoutmat->GetRows() == inoutmat->GetColumns(),
                     "Assuming that input matrix was square");

            // Get local Element mapping for vertex point
            int map = GetVertexMap(vert);

            // Now need to identify a map which takes the local edge
            // mass matrix to the matrix stored in inoutmat;
            // This can currently be deduced from the size of the matrix       
            // - if inoutmat.m_rows() == v_NCoeffs() it is a full
            //   matrix system
            // - if inoutmat.m_rows() == v_NumBndCoeffs() it is a
            //  boundary CG system
             
             int rows = inoutmat->GetRows();
             
             if (rows == GetNcoeffs())
             {
                 // no need to do anything
             }
             else if(rows == NumBndryCoeffs())  // same as NumDGBndryCoeffs()
             {
                 int i;
                 Array<OneD,unsigned int> bmap;
                 GetBoundaryMap(bmap);
                 
                 for(i = 0; i < 2; ++i)
                 {
                     if(map == bmap[i])
                     {
                         map = i;
                         break;
                     }
                 }
                 ASSERTL1(i != 2,"Did not find number in map");
             }
             
             // assumes end points have unit magnitude
             (*inoutmat)(map,map) +=  primCoeffs[0];

        }

        /**
         * Given an edge and vector of element coefficients:
         * - maps those elemental coefficients corresponding to the edge into
         *   an edge-vector.
         * - resets the element coefficients
         * - multiplies the edge vector by the edge mass matrix
         * - maps the edge coefficients back onto the elemental coefficients
         */
        void Expansion1D::v_AddRobinEdgeContribution(const int vert, const Array<OneD, const NekDouble > &primCoeffs, Array<OneD, NekDouble> &coeffs)
        {
            ASSERTL1(IsBoundaryInteriorExpansion(),
                     "Not set up for non boundary-interior expansions");

            int map = GetVertexMap(vert);
            Vmath::Zero(GetNcoeffs(), coeffs, 1);
            coeffs[map] = primCoeffs[0];
        }


    } //end of namespace
} //end of namespace

