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
        
        DNekMatSharedPtr Expansion1D::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            
            DNekMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
                {                    
                ASSERTL1(v_IsBoundaryInteriorExpansion(),
                         "HybridDGHelmholtz matrix not set up "
                         "for non boundary-interior expansions");
                    int       i;
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    int       ncoeffs   = v_GetNcoeffs();

                    int       coordim = v_GetCoordim();
                    DNekScalMat  &invMass = *v_GetLocMatrix(StdRegions::eInvMass);
                    StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                           StdRegions::eWeakDeriv1,
                                                           StdRegions::eWeakDeriv2};
                    DNekMat LocMat(ncoeffs,ncoeffs); 

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,ncoeffs);
                    DNekMat &Mat = *returnval;

                    Vmath::Zero(ncoeffs*ncoeffs,Mat.GetPtr(),1);

                    for(i=0;  i < coordim; ++i)
                    {
                        DNekScalMat &Dmat = *v_GetLocMatrix(DerivType[i]);

                        Mat = Mat + Dmat*invMass*Transpose(Dmat);
                    }

                    // Add end Mass Matrix Contribution
                    DNekScalMat  &Mass = *v_GetLocMatrix(StdRegions::eMass);
                    Mat = Mat + lambdaval*Mass;                    

                    Array<OneD,unsigned int> bmap;
                    v_GetBoundaryMap(bmap);

                    // Add tau*F_e using elemental mass matrices
                    for(i = 0; i < 2; ++i)
                    {
                        Mat(bmap[i],bmap[i]) = Mat(bmap[i],bmap[i]) +  tau; 
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
                        
                        AddHDGHelmholtzTraceTerms(tau,lambda,f);
                        
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
                    int nbndry = v_NumBndryCoeffs();

                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);

                    Array<OneD,unsigned int> bmap;
                    Array<OneD, NekDouble>   lam(2);
                    v_GetBoundaryMap(bmap);
                    
                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,  nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    // Matrix to map Lambda to U
                    DNekScalMat &LamToU = *v_GetLocMatrix(StdRegions::eHybridDGLamToU, lambdaval,tau);
                
                    // Matrix to map Lambda to Q
                    DNekScalMat &LamToQ = *v_GetLocMatrix(StdRegions::eHybridDGLamToQ0, lambdaval,tau);

                    lam[0] = 1.0; lam[1] = 0.0;
                    for(j = 0; j < nbndry; ++j)
                    {
                        BndMat(0,j) = -LamToQ(bmap[0],j) - tau*(LamToU(bmap[0],j) - lam[j]);
                    }

                    lam[0] = 0.0; lam[1] = 1.0;
                    for(j = 0; j < nbndry; ++j)
                    {
                        BndMat(1,j) =  LamToQ(bmap[1],j) - tau*(LamToU(bmap[1],j) - lam[j]);
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
            int nbndry = v_NumBndryCoeffs();
            int nquad  = v_GetNumPoints(0);
            const Array<OneD, const NekDouble> &Basis  = v_GetBasis(0)->GetBdata();
            Array<OneD, unsigned int> vmap;
            
            v_GetBoundaryMap(vmap);
            
            // add G \lambda term (can assume G is diagonal since one
            // of the basis is zero at boundary otherwise)
            for(k = 0; k < nbndry; ++k)
            {
                outarray[vmap[k]] += (Basis[(vmap[k]+1)*nquad-1]*Basis[(vmap[k]+1)*nquad-1] - Basis[vmap[k]*nquad]*Basis[vmap[k]*nquad])*inarray[vmap[k]];
            }
        }

#if 0 
        void Expansion1D::AddHDGHelmholtzMatrixBoundaryTerms(const NekDouble tau, 
                                                             const Array<OneD, const NekDouble> &inarray,
                                                             Array<OneD,NekDouble> &outarray)
        {
            int i,j;
            int nbndry = v_NumBndryCoeffs();
            int nquad  = v_GetNumPoints(0);
            int ncoeffs = v_GetNcoeffs();
            Array<OneD, unsigned int> vmap;
            
            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");
            
            const LibUtilities::BasisSharedPtr base = v_GetBasis(0);
            const Array<OneD, const NekDouble> &Dbasis = base->GetDbdata();
            const Array<OneD, const NekDouble> &Basis  = base->GetBdata();
            
            SpatialDomains::GeomFactorsSharedPtr metricinfo = v_GetMetricInfo();
            Array<TwoD, const NekDouble>  gmat = metricinfo->GetGmat();
            NekDouble rx0,rx1;

            if(metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                rx0 = gmat[0][0];
                rx1 = gmat[0][nquad-1];
            }
            else
            {
                rx0 = rx1 = gmat[0][0];
            }

            v_GetBoundaryMap(vmap);

            // Add -D^T M^{-1}G operation =-<n phi_i, d\phi_j/dx>
            for(i = 0; i < nbndry; ++i)
            {
                for(j = 0; j < ncoeffs; ++j)
                {
                    outarray[vmap[i]] -= Basis[(vmap[i]+1)*nquad-1]*Dbasis[(j+1)*nquad-1]*rx1*inarray[j];
                    outarray[vmap[i]] += Basis[vmap[i]*nquad]*Dbasis[j*nquad]*rx0*inarray[j];
                }
            }

            AddHDGHelmholtzTraceTerms(tau,inarray,outarray);
        }
#endif        

        void Expansion1D::AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                                    const Array<OneD, const NekDouble> &inarray,
                                                    Array<OneD,NekDouble> &outarray)
        {
            int i,j,k,n;
            int nbndry  = v_NumBndryCoeffs();
            int nquad   = v_GetNumPoints(0);
            int ncoeffs = v_GetNcoeffs();
            int coordim = v_GetCoordim();
            NekDouble  val, val1;
            Array<OneD, unsigned int> vmap;
            
            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");


            const Array<OneD, const NekDouble> &Basis  = v_GetBasis(0)->GetBdata();
            DNekScalMat  &invMass = *v_GetLocMatrix(StdRegions::eInvMass); 

            v_GetBoundaryMap(vmap);

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
                
                DNekScalMat &Dmat = *v_GetLocMatrix(DerivType[n]);
                Coeffs = Coeffs  + Dmat*Tmpcoeff; 
            }
        }

    } //end of namespace
} //end of namespace

/** 
 *    $Log: Expansion1D.cpp,v $
 *    Revision 1.3  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.2  2008/08/18 08:30:36  sherwin
 *    Updates for HDG 1D work
 *
 *    Revision 1.1  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *
 **/
