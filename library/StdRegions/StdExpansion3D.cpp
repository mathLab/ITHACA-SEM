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

#include <boost/core/ignore_unused.hpp>

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

            if (out_dx.size() > 0)
            {
                NekDouble  *D0 = &((m_base[0]->GetD())->GetPtr())[0];

                Blas::Dgemm('N','N', nquad0,nquad1*nquad2,nquad0,1.0,
                            D0,nquad0,&wsp[0],nquad0,0.0,&out_dx[0],nquad0);
            }

            if (out_dy.size() > 0)
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

            if (out_dz.size() > 0)
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

        void StdExpansion3D::v_GenStdMatBwdDeriv(
            const int dir,
                  DNekMatSharedPtr &mat)
        {
            ASSERTL1((dir==0)||(dir==1)||(dir==2),"Invalid direction.");

            const int nq0 = m_base[0]->GetNumPoints();
            const int nq1 = m_base[1]->GetNumPoints();
            const int nq2 = m_base[2]->GetNumPoints();
            const int nq  = nq0*nq1*nq2;
            const int nm0 = m_base[0]->GetNumModes();
            const int nm1 = m_base[1]->GetNumModes();
 
            Array<OneD, NekDouble> alloc(4*nq + m_ncoeffs + nm0*nq2*(nq1+nm1),0.0);
            Array<OneD, NekDouble> tmp1 (alloc);               // Quad metric
            Array<OneD, NekDouble> tmp2 (alloc +   nq);        // Dir1 metric
            Array<OneD, NekDouble> tmp3 (alloc + 2*nq);        // Dir2 metric
            Array<OneD, NekDouble> tmp4 (alloc + 3*nq);        // Dir3 metric
            Array<OneD, NekDouble> tmp5 (alloc + 4*nq);        // iprod tmp
            Array<OneD, NekDouble> wsp  (tmp5  +   m_ncoeffs); // Wsp
            switch(dir)
            {
            case 0:
                for(int i=0; i<nq;i++)
                {
                    tmp2[i] =   1.0;
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetBdata(),
                                         tmp2,tmp5,wsp,
                                         false,true,true);

                    tmp2[i] =   0.0;
                    
                    for(int j=0; j<m_ncoeffs;j++)
                    {
                        (*mat)(j,i) =   tmp5[j];
                    }
                }
                break;
            case 1:
                for(int i=0; i<nq;i++)
                {
                    tmp2[i] =   1.0;
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetDbdata(),
                                         m_base[2]->GetBdata(),
                                         tmp2,tmp5,wsp,
                                         true,false,true);

                    tmp2[i] =   0.0;
                    
                    for(int j=0; j<m_ncoeffs;j++)
                    {
                        (*mat)(j,i) =   tmp5[j];
                    }
                }
                break;
            case 2:
                for(int i=0; i<nq;i++)
                {
                    tmp2[i] =   1.0;
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetDbdata(),
                                         tmp2,tmp5,wsp,
                                         true,true,false);
                    tmp2[i] =   0.0;
                    
                    for(int j=0; j<m_ncoeffs;j++)
                    {
                        (*mat)(j,i) =   tmp5[j];
                    }
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Not a 2D expansion.");
                break;
            }
        }

        NekDouble StdExpansion3D::v_PhysEvaluate(
            const Array<OneD, const NekDouble> &coords,
            const Array<OneD, const NekDouble> &physvals)
        {
            Array<OneD, NekDouble> eta(3);

            WARNINGL2(coords[0] >= -1 - NekConstants::kNekZeroTol,
                      "coord[0] < -1");
            WARNINGL2(coords[0] <=  1 + NekConstants::kNekZeroTol,
                      "coord[0] >  1");
            WARNINGL2(coords[1] >= -1 - NekConstants::kNekZeroTol,
                      "coord[1] < -1");
            WARNINGL2(coords[1] <=  1 + NekConstants::kNekZeroTol,
                      "coord[1] >  1");
            WARNINGL2(coords[2] >= -1 - NekConstants::kNekZeroTol,
                      "coord[2] < -1");
            WARNINGL2(coords[2] <=  1 + NekConstants::kNekZeroTol,
                      "coord[2] >  1");

            // Obtain local collapsed corodinate from Cartesian coordinate.
            LocCoordToLocCollapsed(coords, eta);

            const int nq0 = m_base[0]->GetNumPoints();
            const int nq1 = m_base[1]->GetNumPoints();
            const int nq2 = m_base[2]->GetNumPoints();

            Array<OneD, NekDouble> wsp1(nq1 * nq2), wsp2(nq2);

            // Construct the 2D square...
            const NekDouble *ptr = &physvals[0];
            for (int i = 0; i < nq1 * nq2; ++i, ptr += nq0)
            {
                wsp1[i] = StdExpansion::BaryEvaluate<0>(eta[0], ptr);
            }

            for (int i = 0; i < nq2; ++i)
            {
                wsp2[i] = StdExpansion::BaryEvaluate<1>(eta[1], &wsp1[i * nq1]);
            }

            return StdExpansion::BaryEvaluate<2>(eta[2], &wsp2[0]);
        }

        NekDouble StdExpansion3D::v_PhysEvaluate(
            const Array<OneD, DNekMatSharedPtr > &I,
            const Array<OneD, const NekDouble> &physvals)
        {
            NekDouble  value;

            int Qx = m_base[0]->GetNumPoints();
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            Array<OneD, NekDouble> sumFactorization_qr = Array<OneD, NekDouble>(Qy*Qz);
            Array<OneD, NekDouble> sumFactorization_r  = Array<OneD, NekDouble>(Qz);

            // Lagrangian interpolation matrix
            NekDouble *interpolatingNodes = 0;

            // Interpolate first coordinate direction
            interpolatingNodes = &I[0]->GetPtr()[0];

            Blas::Dgemv('T',Qx,Qy*Qz,1.0,&physvals[0],Qx,&interpolatingNodes[0], 1, 0.0, &sumFactorization_qr[0], 1);

            // Interpolate in second coordinate direction
            interpolatingNodes = &I[1]->GetPtr()[0];

            Blas::Dgemv('T',Qy,Qz,1.0,&sumFactorization_qr[0],Qy,&interpolatingNodes[0],1,0.0,&sumFactorization_r[0], 1);

            // Interpolate in third coordinate direction
            interpolatingNodes = &I[2]->GetPtr()[0];
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
            if ( mkey.GetNVarCoeff() == 0 &&
                !mkey.ConstFactorExists(eFactorSVVCutoffRatio))
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
                using std::max;

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

        NekDouble StdExpansion3D::v_Integral(
            const Array<OneD, const NekDouble>& inarray)
        {
            const int nqtot = GetTotPoints();
            Array<OneD, NekDouble> tmp(GetTotPoints());
            v_MultiplyByStdQuadratureMetric(inarray, tmp);
            return Vmath::Vsum(nqtot, tmp, 1);
        }

        int StdExpansion3D::v_GetNedges(void) const
        {
            NEKERROR(ErrorUtil::efatal, "This function is not valid or not defined");
            return 0;
        }

        int StdExpansion3D::v_GetEdgeNcoeffs(const int i) const
        {
            boost::ignore_unused(i);
            NEKERROR(ErrorUtil::efatal, "This function is not valid or not defined");
            return 0;
        }

        void StdExpansion3D::v_GetEdgeInteriorToElementMap(
               const int                  tid,
               Array<OneD, unsigned int> &maparray,
               Array<OneD,          int> &signarray,
               Orientation                traceOrient)
        {
            boost::ignore_unused(tid,maparray,signarray,traceOrient);
            NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
        }

        LibUtilities::BasisKey EvaluateQuadFaceBasisKey(
            const int                     facedir,
            const LibUtilities::BasisType faceDirBasisType,
            const int                     numpoints,
            const int                     nummodes)
        {
            boost::ignore_unused(facedir);

            switch(faceDirBasisType)
            {
                case LibUtilities::eModified_A:
                {
                    const LibUtilities::PointsKey pkey(
                        numpoints, LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(
                        LibUtilities::eModified_A, nummodes, pkey);
                }
                case LibUtilities::eModified_B:
                case LibUtilities::eModified_C:
                {
                    const LibUtilities::PointsKey pkey(
                        numpoints+1, LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(
                        LibUtilities::eModified_A, nummodes, pkey);
                }
                case LibUtilities::eGLL_Lagrange:
                {
                    const LibUtilities::PointsKey pkey(
                        numpoints, LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(
                        LibUtilities::eGLL_Lagrange, nummodes, pkey);
                }
                case LibUtilities::eOrtho_A:
                {
                    const LibUtilities::PointsKey pkey(
                        numpoints, LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(
                        LibUtilities::eOrtho_A, nummodes, pkey);
                }
                case LibUtilities::eOrtho_B:
                case LibUtilities::eOrtho_C:
                {
                    const LibUtilities::PointsKey pkey(
                        numpoints+1, LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(
                        LibUtilities::eOrtho_A, nummodes, pkey);
                }
                default:
                {
                    NEKERROR(ErrorUtil::efatal, "expansion type unknown");
                    break;
                }
            }

            // Keep things happy by returning a value.
            return LibUtilities::NullBasisKey;
        }

        LibUtilities::BasisKey EvaluateTriFaceBasisKey(
            const int                     facedir,
            const LibUtilities::BasisType faceDirBasisType,
            const int                     numpoints,
            const int                     nummodes)
        {
            switch(faceDirBasisType)
            {
                case LibUtilities::eModified_A:
                {
                    const LibUtilities::PointsKey pkey(
                        numpoints, LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(
                        LibUtilities::eModified_A, nummodes, pkey);
                }
                case LibUtilities::eModified_B:
                case LibUtilities::eModified_C:
                case LibUtilities::eModifiedPyr_C:
                {
                    switch (facedir)
                    {
                        case 0:
                        {
                            const LibUtilities::PointsKey pkey(
                                numpoints+1,
                                LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(
                                LibUtilities::eModified_A, nummodes, pkey);
                        }
                        case 1:
                        {
//                            const LibUtilities::PointsKey pkey(
//                                numpoints+1,
//                                LibUtilities::eGaussLobattoLegendre);
                            const LibUtilities::PointsKey pkey(
	 			numpoints,
				LibUtilities::eGaussRadauMAlpha1Beta0);
                            return LibUtilities::BasisKey(
                                LibUtilities::eModified_B, nummodes, pkey);
                        }
                        default:
                        {

                            NEKERROR(ErrorUtil::efatal,"invalid value to flag");
                            break;
                        }
                    }
                    break;
                }

                case LibUtilities::eGLL_Lagrange:
                {
                    switch (facedir)
                    {
                        case 0:
                        {
                            const LibUtilities::PointsKey pkey(
                                numpoints,
                                LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(
                                LibUtilities::eOrtho_A, nummodes, pkey);
                        }
                        case 1:
                        {
                            const LibUtilities::PointsKey pkey(
                                numpoints,
                                LibUtilities::eGaussRadauMAlpha1Beta0);
                            return LibUtilities::BasisKey(
                                LibUtilities::eOrtho_B, nummodes, pkey);
                        }
                        default:
                        {
                            NEKERROR(ErrorUtil::efatal,"invalid value to flag");
                            break;
                        }
                    }
                    break;
                }

                case LibUtilities::eOrtho_A:
                case LibUtilities::eOrtho_B:
                case LibUtilities::eOrtho_C:
                case LibUtilities::eOrthoPyr_C:
                {
                    switch (facedir)
                    {
                        case 0:
                        {
                            const LibUtilities::PointsKey pkey(
                                numpoints,
                                LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(
                                LibUtilities::eOrtho_A, nummodes, pkey);
                        }
                        case 1:
                        {
                            const LibUtilities::PointsKey pkey(
                                numpoints,
                                LibUtilities::eGaussRadauMAlpha1Beta0);
                            return LibUtilities::BasisKey(
                                LibUtilities::eOrtho_B, nummodes, pkey);
                        }
                        default:
                        {
                            NEKERROR(ErrorUtil::efatal,"invalid value to flag");
                            break;
                        }
                    }
                    break;
                }
                default:
                {
                    NEKERROR(ErrorUtil::efatal,"expansion type unknown");
                    break;
                }
            }

            // Keep things happy by returning a value.
            return LibUtilities::NullBasisKey;
        }
    }//end namespace
}//end namespace
