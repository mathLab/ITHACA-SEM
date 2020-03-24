///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion2D.cpp
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
// which are common to 2D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion2D.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
    namespace StdRegions
    {

        StdExpansion2D::StdExpansion2D()
        {
        }

        StdExpansion2D::StdExpansion2D(int numcoeffs,
                       const LibUtilities::BasisKey &Ba,
                                       const LibUtilities::BasisKey &Bb):
        StdExpansion(numcoeffs,2, Ba, Bb)
        {
        }

        StdExpansion2D::StdExpansion2D(const StdExpansion2D &T):
                StdExpansion(T)
        {
        }

        StdExpansion2D::~StdExpansion2D()
        {
        }

        //----------------------------
        // Differentiation Methods
        //----------------------------
        void StdExpansion2D::PhysTensorDeriv(const Array<OneD, const NekDouble>& inarray,
                         Array<OneD, NekDouble> &outarray_d0,
                         Array<OneD, NekDouble> &outarray_d1)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            if (outarray_d0.size() > 0) // calculate du/dx_0
            {
                DNekMatSharedPtr D0 = m_base[0]->GetD();
                if(inarray.data() == outarray_d0.data())
                {
                    Array<OneD, NekDouble> wsp(nquad0 * nquad1);
                    Vmath::Vcopy(nquad0 * nquad1,inarray.get(),1,wsp.get(),1);
                    Blas::Dgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                                &(D0->GetPtr())[0], nquad0, &wsp[0], nquad0, 0.0,
                                &outarray_d0[0], nquad0);
                }
                else
                {
                    Blas::Dgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                                &(D0->GetPtr())[0], nquad0, &inarray[0], nquad0, 0.0,
                                &outarray_d0[0], nquad0);
                }
            }

            if (outarray_d1.size() > 0) // calculate du/dx_1
            {
                DNekMatSharedPtr D1 = m_base[1]->GetD();
                if(inarray.data() == outarray_d1.data())
                {
                    Array<OneD, NekDouble> wsp(nquad0 * nquad1);
                    Vmath::Vcopy(nquad0 * nquad1,inarray.get(),1,wsp.get(),1);
                    Blas:: Dgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &wsp[0], nquad0,
                                 &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
                }
                else
                {
                    Blas:: Dgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &inarray[0], nquad0,
                                 &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
                }
            }
        }

        NekDouble StdExpansion2D::v_PhysEvaluate(const Array<OneD, const NekDouble>& coords, const Array<OneD, const NekDouble> & physvals)
        {
            Array<OneD, NekDouble> coll(2);
            Array<OneD, DNekMatSharedPtr>  I(2);

            ASSERTL2(coords[0] > -1 - NekConstants::kNekZeroTol, "coord[0] < -1");
            ASSERTL2(coords[0] <  1 + NekConstants::kNekZeroTol, "coord[0] >  1");
            ASSERTL2(coords[1] > -1 - NekConstants::kNekZeroTol, "coord[1] < -1");
            ASSERTL2(coords[1] <  1 + NekConstants::kNekZeroTol, "coord[1] >  1");

            LocCoordToLocCollapsed(coords,coll);

            I[0] = m_base[0]->GetI(coll);
            I[1] = m_base[1]->GetI(coll+1);

            return v_PhysEvaluate(I,physvals);
        }

        NekDouble StdExpansion2D::v_PhysEvaluate(
            const Array<OneD, DNekMatSharedPtr > &I,
            const Array<OneD, const NekDouble> &physvals)
        {
            NekDouble val;
            int i;
            int nq0 = m_base[0]->GetNumPoints();
            int nq1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> wsp1(nq1);

            // interpolate first coordinate direction
            for (i = 0; i < nq1;++i)
            {
                wsp1[i] = Blas::Ddot(nq0, &(I[0]->GetPtr())[0], 1,
                                     &physvals[i * nq0], 1);
            }

            // interpolate in second coordinate direction
            val = Blas::Ddot(nq1, I[1]->GetPtr(), 1, wsp1, 1);

            return val;
        }

        //////////////////////////////
        // Integration Methods
        //////////////////////////////

        NekDouble StdExpansion2D::Integral(const Array<OneD, const NekDouble>& inarray,
                       const Array<OneD, const NekDouble>& w0,
                       const Array<OneD, const NekDouble>& w1)
        {
            int i;
            NekDouble Int = 0.0;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> tmp(nquad0 * nquad1);

            // multiply by integration constants
            for (i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0, &inarray[0] + i*nquad0, 1, w0.get(),
                            1, &tmp[0] + i*nquad0, 1);
            }

            for (i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1, &tmp[0]+ i, nquad0, w1.get(), 1,
                            &tmp[0] + i, nquad0);
            }
            Int = Vmath::Vsum(nquad0 * nquad1, tmp, 1);

            return Int;
        }

        void StdExpansion2D::BwdTrans_SumFacKernel(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray,
                Array<OneD, NekDouble> &wsp,
                bool doCheckCollDir0,
                bool doCheckCollDir1)
        {
            v_BwdTrans_SumFacKernel(base0, base1, inarray, outarray, wsp, doCheckCollDir0, doCheckCollDir1);
        }

        void StdExpansion2D::IProductWRTBase_SumFacKernel(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray,
                Array<OneD, NekDouble> &wsp,
                bool doCheckCollDir0,
                bool doCheckCollDir1)
        {
            v_IProductWRTBase_SumFacKernel(base0, base1, inarray, outarray, wsp, doCheckCollDir0, doCheckCollDir1);
        }

        void StdExpansion2D::v_LaplacianMatrixOp_MatFree(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            if (mkey.GetNVarCoeff() == 0
                &&!mkey.ConstFactorExists(StdRegions::eFactorSVVCutoffRatio))
            {
                using std::max;

                // This implementation is only valid when there are no
                // coefficients associated to the Laplacian operator
                int       nquad0  = m_base[0]->GetNumPoints();
                int       nquad1  = m_base[1]->GetNumPoints();
                int       nqtot   = nquad0*nquad1;
                int       nmodes0 = m_base[0]->GetNumModes();
                int       nmodes1 = m_base[1]->GetNumModes();
                int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);

                const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();

                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(4*wspsize);      // size wspsize
                Array<OneD,NekDouble> wsp1(wsp0+wspsize);   // size 3*wspsize

                if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                {
                    // LAPLACIAN MATRIX OPERATION
                    // wsp0 = u       = B   * u_hat
                    // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                    // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                    BwdTrans_SumFacKernel(base0,base1,inarray,wsp0,wsp1,true,true);
                    LaplacianMatrixOp_MatFree_Kernel(wsp0, outarray, wsp1);
                }
                else
                {
                    LaplacianMatrixOp_MatFree_Kernel(inarray, outarray, wsp1);
                }
            }
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(
                    inarray,outarray,mkey);
            }
        }



        void StdExpansion2D::v_HelmholtzMatrixOp_MatFree(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            if (mkey.GetNVarCoeff() == 0
                &&!mkey.ConstFactorExists(StdRegions::eFactorSVVCutoffRatio))
            {
                using std::max;

                int       nquad0  = m_base[0]->GetNumPoints();
                int       nquad1  = m_base[1]->GetNumPoints();
                int       nqtot   = nquad0*nquad1;
                int       nmodes0 = m_base[0]->GetNumModes();
                int       nmodes1 = m_base[1]->GetNumModes();
                int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),
                                        nquad0*nmodes1);
                NekDouble lambda  =
                    mkey.GetConstFactor(StdRegions::eFactorLambda);

                const Array<OneD, const NekDouble>& base0 = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1 = m_base[1]->GetBdata();

                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(5*wspsize);      // size wspsize
                Array<OneD,NekDouble> wsp1(wsp0 + wspsize);  // size wspsize
                Array<OneD,NekDouble> wsp2(wsp0 + 2*wspsize);// size 3*wspsize

                if (!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                {
                    // MASS MATRIX OPERATION
                    // The following is being calculated:
                    // wsp0     = B   * u_hat = u
                    // wsp1     = W   * wsp0
                    // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat
                    BwdTrans_SumFacKernel       (base0, base1, inarray,
                                                 wsp0, wsp2,true,true);
                    MultiplyByQuadratureMetric  (wsp0, wsp1);
                    IProductWRTBase_SumFacKernel(base0, base1, wsp1, outarray,
                                                 wsp2, true, true);

                    LaplacianMatrixOp_MatFree_Kernel(wsp0, wsp1, wsp2);
                }
                else
                {
                    MultiplyByQuadratureMetric(inarray,outarray);
                    LaplacianMatrixOp_MatFree_Kernel(inarray, wsp1, wsp2);
                }

                // outarray = lambda * outarray + wsp1
                //          = (lambda * M + L ) * u_hat
                Vmath::Svtvp(m_ncoeffs, lambda, &outarray[0], 1,
                              &wsp1[0], 1, &outarray[0], 1);
            }
            else
            {
                StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(
                    inarray,outarray,mkey);
            }
        }

    } //end namespace
} //end namespace
