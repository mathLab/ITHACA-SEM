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
            case StdRegions::eHybridDGHelmBndLam:
                {
                    int j;
                    int nbndry = v_NumBndryCoeffs();

                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);

                    Array<OneD,unsigned int> bmap;
                    v_GetBoundaryMap(bmap);
                    
                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,  nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    // Matrix to map Lambda to U
                    DNekScalMat &LamToU = *v_GetLocMatrix(StdRegions::eHybridDGLamToU, lambdaval,tau);
                
                    // Matrix to map Lambda to Q
                    DNekScalMat &LamToQ = *v_GetLocMatrix(StdRegions::eHybridDGLamToQ0, lambdaval,tau);
                    
                    for(j = 0; j < nbndry; ++j)
                    {
                        BndMat(0,j) = -LamToQ(bmap[0],j) + tau*LamToU(bmap[0],j);
                    }

                    for(j = 0; j < nbndry; ++j)
                    {
                        BndMat(1,j) =  LamToQ(bmap[1],j) - tau*LamToU(bmap[1],j);
                    }
                }
                break;
            default:
                returnval = Expansion::GenMatrix(mkey);
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
            
            // add G (\lambda - ulam) = G x F term (can
            // assume G is diagonal since one of the basis
            // is zero at boundary otherwise)
            for(k = 0; k < nbndry; ++k)
            {
                outarray[vmap[k]] += (Basis[(vmap[k]+1)*nquad-1]*Basis[(vmap[k]+1)*nquad-1] - Basis[vmap[k]*nquad]*Basis[vmap[k]*nquad])*inarray[vmap[k]];
            }
        }


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
        

        void Expansion1D::AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                                    const Array<OneD, const NekDouble> &inarray,
                                                    Array<OneD,NekDouble> &outarray)
        {
            int i,j,k,n;
            int nbndry = v_NumBndryCoeffs();
            int nquad  = v_GetNumPoints(0);
            int ncoeffs = v_GetNcoeffs();
            NekDouble  val, val1;
            Array<OneD, unsigned int> vmap;
            
            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");

            const Array<OneD, const NekDouble> &Dbasis = v_GetBasis(0)->GetDbdata();
            const Array<OneD, const NekDouble> &Basis  = v_GetBasis(0)->GetBdata();
            
            DNekScalMat  &invMass = *v_GetLocMatrix(StdRegions::eInvMass); 

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

            //Add -E^T M^{-1}D_i^e = -< d\phi_i/dx, n  phi_j>
            for(i = 0; i < ncoeffs; ++i)
            {
                for(j = 0; j < nbndry; ++j)
                {
                    outarray[i] -= Dbasis[(i+1)*nquad-1]*rx1*Basis[(vmap[j]+1)*nquad-1]*inarray[vmap[j]];
                    outarray[i] += Dbasis[i*nquad]*rx0*Basis[vmap[j]*nquad]*inarray[vmap[j]];
                }
            }
            
            // Add F = \tau <phi_i,phi_j> (note phi_i is zero if phi_j is non-zero)
            for(i = 0; i < nbndry; ++i)
            {
                outarray[vmap[i]] += tau*Basis[(vmap[i]+1)*nquad-1]*Basis[(vmap[i]+1)*nquad-1]*inarray[vmap[i]];
                outarray[vmap[i]] += tau*Basis[vmap[i]*nquad]*Basis[vmap[i]*nquad]*inarray[vmap[i]];
            }

            // Add E M^{-1} G term 
            for(i = 0; i < nbndry; ++i)
            {
                for(n = 0; n < nbndry; ++n)
                {
                    // evaluate M^{-1} G
                    val1 = 0.0;
                    for(k = 0; k < nbndry; ++k)
                    {
                        val = 0.0;
                        for(j = 0; j < nbndry; ++j)
                        {
                            val += (Basis[(vmap[k]+1)*nquad-1]*Basis[(vmap[j]+1)*nquad-1] - Basis[vmap[k]*nquad]*Basis[vmap[j]*nquad])*inarray[vmap[j]];
                        }
                        
                        val1 += invMass(vmap[n],vmap[k])*val; 
                    }

                    outarray[vmap[i]] += (Basis[(vmap[i]+1)*nquad-1]*Basis[(vmap[n]+1)*nquad-1] - Basis[vmap[i]*nquad]*Basis[vmap[n]*nquad])*val1; 
                }
            }
        }

    } //end of namespace
} //end of namespace

/** 
 *    $Log: Expansion1D.cpp,v $
 *
 **/
