///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion3D.cpp
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 3D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion3D.h>

#ifdef max
#undef max
#endif

namespace Nektar

{
    namespace StdRegions
    {
        StdExpansion3D::StdExpansion3D()
        {
        }
        
        StdExpansion3D::StdExpansion3D(int                           numcoeffs, 
                                       const LibUtilities::BasisKey &Ba,
                                       const LibUtilities::BasisKey &Bb, 
                                       const LibUtilities::BasisKey &Bc) :
            StdExpansion(numcoeffs,3,Ba,Bb,Bc)
        {
        }

        StdExpansion3D::StdExpansion3D(const StdExpansion3D &T):
            StdExpansion(T)
        {
        }
        
        StdExpansion3D::~StdExpansion3D()
        {
        }
        void StdExpansion3D::PhysTensorDeriv(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &out_dx,
                  Array<OneD,       NekDouble> &out_dy,
                  Array<OneD,       NekDouble> &out_dz)
        {
            const int nquad0 = m_base[0]->GetNumPoints();
            const int nquad1 = m_base[1]->GetNumPoints();
            const int nquad2 = m_base[2]->GetNumPoints();

            Array<OneD, NekDouble> wsp(nquad0*nquad1*nquad2);
            
            // copy inarray to wsp in case inarray is used as outarray
            Vmath::Vcopy(nquad0*nquad1*nquad2, &inarray[0], 1, &wsp[0], 1);
            
            if (out_dx.num_elements() > 0)
            {
                NekDouble  *D0 = &((m_base[0]->GetD())->GetPtr())[0];

                Blas::Dgemm('N','N', nquad0,nquad1*nquad2,nquad0,1.0,
                            D0,nquad0,&wsp[0],nquad0,0.0,&out_dx[0],nquad0);
            }

            if (out_dy.num_elements() > 0) 
            {
                NekDouble   *D1 = &((m_base[1]->GetD())->GetPtr())[0];
                for (int j = 0; j < nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', nquad0, nquad1,      nquad1,
                                1.0, &wsp[j*nquad0*nquad1],    nquad0,
                                D1,                            nquad1,
                                0.0, &out_dy[j*nquad0*nquad1], nquad0);
                }
            }

            if (out_dz.num_elements() > 0) 
            {
                NekDouble     *D2 = &((m_base[2]->GetD())->GetPtr())[0];

                Blas::Dgemm('N','T',nquad0*nquad1,nquad2,nquad2,1.0,
                            &wsp[0],nquad0*nquad1,D2,nquad2,0.0,&out_dz[0],
                            nquad0*nquad1);
            }
        }

        void StdExpansion3D::BwdTrans_SumFacKernel(
            const Array<OneD, const NekDouble>& base0,
            const Array<OneD, const NekDouble>& base1,
            const Array<OneD, const NekDouble>& base2,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray,
                  Array<OneD,       NekDouble>& wsp,
            bool                                doCheckCollDir0,
            bool                                doCheckCollDir1,
            bool                                doCheckCollDir2)
        {
            v_BwdTrans_SumFacKernel(base0, base1, base2, inarray, outarray, wsp, doCheckCollDir0, doCheckCollDir1, doCheckCollDir2);
        }

        void StdExpansion3D::IProductWRTBase_SumFacKernel(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& base2,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray,
                      Array<OneD, NekDouble> &wsp,
                bool doCheckCollDir0,
                bool doCheckCollDir1,
                bool doCheckCollDir2)
        {
            v_IProductWRTBase_SumFacKernel(base0, base1, base2, inarray, outarray, wsp, doCheckCollDir0, doCheckCollDir1, doCheckCollDir2);
        }

        NekDouble StdExpansion3D::v_PhysEvaluate(
            const Array<OneD, const NekDouble> &coords)
        {
            return PhysEvaluate(coords,m_phys);
        }
        
        NekDouble StdExpansion3D::v_PhysEvaluate(
            const Array<OneD, const NekDouble> &coords, 
            const Array<OneD, const NekDouble> &physvals)
        {
            NekDouble  value;
            ASSERTL2(coords[0] >= -1,"coord[0] < -1");
            ASSERTL2(coords[0] <=  1,"coord[0] >  1");
            ASSERTL2(coords[1] >= -1,"coord[1] < -1");
            ASSERTL2(coords[1] <=  1,"coord[1] >  1");
            ASSERTL2(coords[2] >= -1,"coord[2] < -1");
            ASSERTL2(coords[2] <=  1,"coord[2] >  1");
            
            int Qx = m_base[0]->GetNumPoints();
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();
            
            Array<OneD, NekDouble> sumFactorization_qr = Array<OneD, NekDouble>(Qy*Qz);
            Array<OneD, NekDouble> sumFactorization_r  = Array<OneD, NekDouble>(Qz);
            
            // Lagrangian interpolation matrix
            DNekMatSharedPtr I;
            NekDouble *interpolatingNodes = 0;
            
            // Interpolate first coordinate direction
            I = m_base[0]->GetI(coords);
            interpolatingNodes = &I->GetPtr()[0];
            for(int i = 0; i < Qy*Qz;++i)
            {
                sumFactorization_qr[i] =  Blas::Ddot(Qx, interpolatingNodes, 1, &physvals[ i*Qx ], 1);
            }
            
            // Interpolate in second coordinate direction
            I = m_base[1]->GetI(coords+1);
            interpolatingNodes = &I->GetPtr()[0];
            for(int j =0; j < Qz; ++j)
            {
                sumFactorization_r[j] = Blas::Ddot(Qy, interpolatingNodes, 1, &sumFactorization_qr[ j*Qy ], 1);
            }
            
            // Interpolate in third coordinate direction
            I = m_base[2]->GetI(coords+2);
            interpolatingNodes = &I->GetPtr()[0];
            value = Blas::Ddot(Qz, interpolatingNodes, 1, &sumFactorization_r[0], 1);
            
            return value;
        }
        

        /**
         * @param   inarray     Input coefficients.
         * @param   output      Output coefficients.
         * @param   mkey        Matrix key
         */
        void StdExpansion3D::v_LaplacianMatrixOp_MatFree(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            if(mkey.GetNVarCoeff() == 0)
            {
                // This implementation is only valid when there are no
                // coefficients associated to the Laplacian operator
                int nqtot = GetTotPoints();

                const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
                const Array<OneD, const NekDouble>& base2  = m_base[2]->GetBdata();

                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(7*nqtot);
                Array<OneD,NekDouble> wsp1(wsp0+nqtot);

                if(!(m_base[0]->Collocation() && m_base[1]->Collocation() &&
                     m_base[2]->Collocation()))
                {
                    // LAPLACIAN MATRIX OPERATION
                    // wsp0 = u       = B   * u_hat
                    // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                    // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                    BwdTrans_SumFacKernel(base0,base1,base2,inarray,wsp0,wsp1,true,true,true);
                    LaplacianMatrixOp_MatFree_Kernel(wsp0,outarray,wsp1);
                }
                else
                {
                    LaplacianMatrixOp_MatFree_Kernel(inarray,outarray,wsp1);
                }
            }
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }
        }


        void StdExpansion3D::v_HelmholtzMatrixOp_MatFree(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            if(mkey.GetNVarCoeff() == 0)
            {
                int nquad0  = m_base[0]->GetNumPoints();
                int nquad1  = m_base[1]->GetNumPoints();
                int nquad2  = m_base[2]->GetNumPoints();
                int nmodes0 = m_base[0]->GetNumModes();
                int nmodes1 = m_base[1]->GetNumModes();
                int nmodes2 = m_base[2]->GetNumModes();
                int wspsize = max(nquad0*nmodes2*(nmodes1+nquad1),
                                  nquad0*nquad1*(nquad2+nmodes0)+
                                  nmodes0*nmodes1*nquad2);

                NekDouble lambda  = mkey.GetConstFactor(StdRegions::eFactorLambda);

                const Array<OneD, const NekDouble>& base0 = m_base[0]->GetBdata ();
                const Array<OneD, const NekDouble>& base1 = m_base[1]->GetBdata ();
                const Array<OneD, const NekDouble>& base2 = m_base[2]->GetBdata ();
                Array<OneD,NekDouble> wsp0(8*wspsize);
                Array<OneD,NekDouble> wsp1(wsp0+1*wspsize);
                Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);

                if(!(m_base[0]->Collocation() && m_base[1]->Collocation() &&
                     m_base[2]->Collocation()))
                {
                    // MASS MATRIX OPERATION
                    // The following is being calculated:
                    // wsp0     = B   * u_hat = u
                    // wsp1     = W   * wsp0
                    // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat
                    BwdTrans_SumFacKernel           (base0,base1,base2,inarray,
                                                     wsp0,wsp2,true,true,true);
                    MultiplyByQuadratureMetric      (wsp0,wsp1);
                    IProductWRTBase_SumFacKernel    (base0,base1,base2,wsp1,
                                                     outarray,wsp2,true,true,true);
                    LaplacianMatrixOp_MatFree_Kernel(wsp0,wsp1,wsp2);
                }
                else
                {
                    // specialised implementation for the classical spectral
                    // element method
                    MultiplyByQuadratureMetric      (inarray,outarray);
                    LaplacianMatrixOp_MatFree_Kernel(inarray,wsp1,wsp2);
                }

                // outarray = lambda * outarray + wsp1
                //          = (lambda * M + L ) * u_hat
                Vmath::Svtvp(m_ncoeffs,lambda,&outarray[0],1,&wsp1[0],1,
                             &outarray[0],1);
            }
            else
            {
                StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }

        }

        const NormalVector & StdExpansion3D::v_GetSurfaceNormal() const
        {
            return m_surfaceNormal;
        }
        
        const NormalVector & StdExpansion3D::v_GetFaceNormal(const int face) const
        {
            std::map<int, NormalVector>::const_iterator x;
            x = m_faceNormals.find(face);
            ASSERTL0 (x != m_faceNormals.end(),
                      "face normal not computed.");
            return x->second;
        }

        void StdExpansion3D::v_NegateFaceNormal(const int face)
        {
            m_negatedNormals[face] = true;
            for (int i = 0; i < GetCoordim(); ++i)
            {
                Vmath::Neg(m_faceNormals[face][i].num_elements(), 
                           m_faceNormals[face][i], 1);
            }
        }
    }//end namespace
}//end namespace
