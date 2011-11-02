///////////////////////////////////////////////////////////////////////////////
//
// File QuadExp.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////
#include <LocalRegions/LocalRegions.h>
#include <LocalRegions/LocalRegions.hpp>
#include <stdio.h>
#include <LocalRegions/QuadExp.h>


namespace Nektar
{
    namespace LocalRegions
    {
        QuadExp::QuadExp(const LibUtilities::BasisKey &Ba,
                         const LibUtilities::BasisKey &Bb,
                         const SpatialDomains::QuadGeomSharedPtr &geom):
             StdExpansion  (Ba.GetNumModes()*Bb.GetNumModes(),2,Ba,Bb),
             Expansion     (),
             StdExpansion2D(Ba.GetNumModes()*Bb.GetNumModes(),Ba,Bb),
             StdQuadExp(Ba,Bb),
            m_geom(geom),
            m_metricinfo(m_geom->GetGeomFactors(m_base)),
            m_matrixManager(std::string("QuadExpMatrix")),
            m_staticCondMatrixManager(std::string("QuadExpStaticCondMatrix"))
        {
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i, StdRegions::eNoExpansionType,*this), boost::bind(&QuadExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoExpansionType,*this),
                                                          boost::bind(&QuadExp::CreateStaticCondMatrix, this, _1));
            }
        }

        QuadExp::QuadExp(const QuadExp &T):
            StdExpansion(T),
            Expansion   (),
            StdExpansion2D(T),
            StdQuadExp(T),
            m_geom(T.m_geom),
            m_metricinfo(T.m_metricinfo),
            m_matrixManager(std::string("QuadExpMatrix")),
            m_staticCondMatrixManager(std::string("QuadExpStaticCondMatrix"))
        {
        }

        // by default the StdQuadExp destructor will be called
        QuadExp::~QuadExp()
        {
        }

        //----------------------------
        // Integration Methods
        //----------------------------

        /** \brief Integrate the physical point list \a inarray over region
            and return the value

            Inputs:\n

            - \a inarray: definition of function to be returned at quadrature point
            of expansion.

            Outputs:\n

            - returns \f$\int^1_{-1}\int^1_{-1} u(\xi_1, \xi_2) J[i,j] d
            \xi_1 d \xi_2 \f$ where \f$inarray[i,j] =
            u(\xi_{1i},\xi_{2j}) \f$ and \f$ J[i,j] \f$ is the
            Jacobian evaluated at the quadrature point.
        */
        NekDouble QuadExp::Integral(const Array<OneD, const NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            NekDouble ival;
            Array<OneD,NekDouble> tmp(nquad0*nquad1);

            // multiply inarray with Jacobian

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1, jac, 1, inarray, 1, tmp, 1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1, jac[0], inarray, 1, tmp, 1);
            }

            // call StdQuadExp version;
            ival = StdQuadExp::Integral(tmp);
            return  ival;
        }

        void QuadExp::IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble> &outarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();

            Array<OneD,NekDouble> tmp(nquad0*nquad1+nquad1*order0);
            Array<OneD,NekDouble> wsp(tmp+nquad0*nquad1);

            MultiplyByQuadratureMetric(inarray,tmp);
            StdQuadExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                     m_base[1]->GetBdata(),
                                                     tmp,outarray,wsp,true,true);
        }

        void QuadExp::IProductWRTBase_MatOp(const Array<OneD, const NekDouble>& inarray,
                                            Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            MatrixKey      iprodmatkey(StdRegions::eIProductWRTBase,DetExpansionType(),*this);
            DNekScalMatSharedPtr& iprodmat = m_matrixManager[iprodmatkey];

            Blas::Dgemv('N',m_ncoeffs,nq,iprodmat->Scale(),(iprodmat->GetOwnedMatrix())->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);

        }

        void QuadExp::IProductWRTDerivBase_SumFac(const int dir,
                                                  const Array<OneD, const NekDouble>& inarray,
                                                  Array<OneD, NekDouble> & outarray)
        {
            ASSERTL1((dir==0)||(dir==1)||(dir==2),"Invalid direction.");
            ASSERTL1((dir==2)?(m_geom->GetCoordim()==3):true,"Invalid direction.");

            int    nquad0  = m_base[0]->GetNumPoints();
            int    nquad1  = m_base[1]->GetNumPoints();
            int    nqtot   = nquad0*nquad1;
            int    nmodes0 = m_base[0]->GetNumModes();

            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD, NekDouble> tmp1(2*nqtot+m_ncoeffs+nmodes0*nquad1);
            Array<OneD, NekDouble> tmp2(tmp1 +   nqtot);
            Array<OneD, NekDouble> tmp3(tmp1 + 2*nqtot);
            Array<OneD, NekDouble> tmp4(tmp1 + 2*nqtot+m_ncoeffs);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nqtot,&gmat[2*dir][0],  1,inarray.get(),1,tmp1.get(),1);
                Vmath::Vmul(nqtot,&gmat[2*dir+1][0],1,inarray.get(),1,tmp2.get(),1);
            }
            else
            {
                Vmath::Smul(nqtot, gmat[2*dir][0],  inarray.get(),1,tmp1.get(), 1);
                Vmath::Smul(nqtot, gmat[2*dir+1][0],inarray.get(),1,tmp2.get(), 1);
            }

            MultiplyByQuadratureMetric(tmp1,tmp1);
            MultiplyByQuadratureMetric(tmp2,tmp2);

            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),m_base[1]->GetBdata(), tmp1,tmp3,    tmp4,false,true);
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata() ,m_base[1]->GetDbdata(),tmp2,outarray,tmp4,true,false);
            Vmath::Vadd(m_ncoeffs, tmp3, 1, outarray, 1, outarray, 1);
        }

        void QuadExp::IProductWRTDerivBase_MatOp(const int dir,
                                                 const Array<OneD, const NekDouble>& inarray,
                                                 Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdRegions::MatrixType mtype;

            switch(dir)
            {
            case 0:
                {
                    mtype = StdRegions::eIProductWRTDerivBase0;
                }
                break;
            case 1:
                {
                    mtype = StdRegions::eIProductWRTDerivBase1;
                }
                break;
            case 2:
                {
                    mtype = StdRegions::eIProductWRTDerivBase2;
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }

            MatrixKey      iprodmatkey(mtype,DetExpansionType(),*this);
            DNekScalMatSharedPtr& iprodmat = m_matrixManager[iprodmatkey];

            Blas::Dgemv('N',m_ncoeffs,nq,iprodmat->Scale(),(iprodmat->GetOwnedMatrix())->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void QuadExp::GeneralMatrixOp_MatOp(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD,NekDouble> &outarray,
                                            const StdRegions::StdMatrixKey &mkey)
        {
            int nConsts = mkey.GetNconstants();
            DNekScalMatSharedPtr   mat;

            switch(nConsts)
            {
            case 0:
                {
                    mat = GetLocMatrix(mkey.GetMatrixType());
                }
                break;
            case 1:
                {
                    mat = GetLocMatrix(mkey.GetMatrixType(),mkey.GetConstant(0));
                }
                break;
            case 2:
                {
                    mat = GetLocMatrix(mkey.GetMatrixType(),mkey.GetConstant(0),mkey.GetConstant(1));
                }
                break;

            default:
                {
                    NEKERROR(ErrorUtil::efatal, "Unknown number of constants");
                }
                break;
            }

            if(inarray.get() == outarray.get())
            {
                Array<OneD,NekDouble> tmp(m_ncoeffs);
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,tmp.get(),1);

                Blas::Dgemv('N',m_ncoeffs,m_ncoeffs,mat->Scale(),(mat->GetOwnedMatrix())->GetPtr().get(),
                            m_ncoeffs, tmp.get(), 1, 0.0, outarray.get(), 1);
            }
            else
            {
                Blas::Dgemv('N',m_ncoeffs,m_ncoeffs,mat->Scale(),(mat->GetOwnedMatrix())->GetPtr().get(),
                            m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
            }
        }

        void QuadExp::MultiplyByQuadratureMetric(const Array<OneD, const NekDouble>& inarray,
                                                 Array<OneD, NekDouble> &outarray)
        {
            if(m_metricinfo->IsUsingQuadMetrics())
            {
                int    nqtot = m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints();
                const Array<OneD, const NekDouble>& metric = m_metricinfo->GetQuadratureMetrics();

                Vmath::Vmul(nqtot, metric, 1, inarray, 1, outarray, 1);
            }
            else
            {
                int    i;
                int    nquad0 = m_base[0]->GetNumPoints();
                int    nquad1 = m_base[1]->GetNumPoints();
                int    nqtot  = nquad0*nquad1;

                const Array<OneD, const NekDouble>& jac = m_metricinfo->GetJac();
                const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
                const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();

                if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                {
                    Vmath::Vmul(nqtot, jac, 1, inarray, 1, outarray, 1);
                }
                else
                {
                    Vmath::Smul(nqtot, jac[0], inarray, 1, outarray, 1);
                }

                // multiply by integration constants
                for(i = 0; i < nquad1; ++i)
                {
                    Vmath::Vmul(nquad0, outarray.get()+i*nquad0,1,
                                w0.get(),1,outarray.get()+i*nquad0,1);
                }

                for(i = 0; i < nquad0; ++i)
                {
                    Vmath::Vmul(nquad1,outarray.get()+i,nquad0,w1.get(),1,
                                outarray.get()+i,nquad0);
                }
            }
        }

        /**
         * @param   inarray     Input coefficients.
         * @param   output      Output coefficients.
         * @param   mkey        Matrix key
         */
        void QuadExp::LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                Array<OneD,NekDouble> &outarray,
                                                const StdRegions::StdMatrixKey &mkey)
        {
            if(mkey.GetNvariableLaplacianCoefficients() == 0)
            {
                // This implementation is only valid when there are no
                // coefficients associated to the Laplacian operator
                if(m_metricinfo->IsUsingLaplMetrics())
                {
                    int       nquad0  = m_base[0]->GetNumPoints();
                    int       nquad1  = m_base[1]->GetNumPoints();
                    int       nqtot   = nquad0*nquad1;
                    int       nmodes0 = m_base[0]->GetNumModes();
                    int       nmodes1 = m_base[1]->GetNumModes();
                    int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);

                    const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                    const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
                    const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
                    const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
                    const Array<TwoD, const NekDouble>& metric = m_metricinfo->GetLaplacianMetrics();

                    // Allocate temporary storage
                    Array<OneD,NekDouble> wsp0(3*wspsize);
                    Array<OneD,NekDouble> wsp1(wsp0+wspsize);
                    Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);

                    if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                    {
                        // LAPLACIAN MATRIX OPERATION
                        // wsp0 = u       = B   * u_hat
                        // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                        // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                        BwdTrans_SumFacKernel(base0,base1,inarray,wsp0,wsp1,true,true);
                        StdExpansion2D::PhysTensorDeriv(wsp0,wsp1,wsp2);
                    }
                    else
                    {
                        StdExpansion2D::PhysTensorDeriv(inarray,wsp1,wsp2);
                    }

                    // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
                    // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
                    // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
                    // especially for this purpose
                    if(!m_metricinfo->LaplacianMetricIsZero(1))
                    {
                        Vmath::Vvtvvtp(nqtot,&metric[0][0],1,&wsp1[0],1,&metric[1][0],1,&wsp2[0],1,&wsp0[0],1);
                        Vmath::Vvtvvtp(nqtot,&metric[1][0],1,&wsp1[0],1,&metric[2][0],1,&wsp2[0],1,&wsp2[0],1);
                    }
                    else
                    {
                        // special implementation in case g1 = 0 (which should hold for undistorted quads)
                        // wsp0 = k = g0 * wsp1 = g0 * du_dxi1
                        // wsp2 = l = g2 * wsp2 = g2 * du_dxi2
                        Vmath::Vmul(nqtot,&metric[0][0],1,&wsp1[0],1,&wsp0[0],1);
                        Vmath::Vmul(nqtot,&metric[2][0],1,&wsp2[0],1,&wsp2[0],1);
                    }

                    // outarray = m = (D_xi1 * B)^T * k
                    // wsp1     = n = (D_xi2 * B)^T * l
                    IProductWRTBase_SumFacKernel(dbase0,base1,wsp0,outarray,wsp1,false,true);
                    IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp1,    wsp0,true,false);

                    // outarray = outarray + wsp1
                    //          = L * u_hat
                    Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);
                }
                else
                {
                    int       nquad0  = m_base[0]->GetNumPoints();
                    int       nquad1  = m_base[1]->GetNumPoints();
                    int       nqtot   = nquad0*nquad1;
                    int       nmodes0 = m_base[0]->GetNumModes();
                    int       nmodes1 = m_base[1]->GetNumModes();
                    int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);

                    const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                    const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
                    const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
                    const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();

                    // Allocate temporary storage
                    Array<OneD,NekDouble> wsp0(6*wspsize);
                    Array<OneD,NekDouble> wsp1(wsp0+wspsize);
                    Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);
                    Array<OneD,NekDouble> wsp3(wsp0+3*wspsize);
                    Array<OneD,NekDouble> wsp4(wsp0+4*wspsize);
                    Array<OneD,NekDouble> wsp5(wsp0+5*wspsize);

                    if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                    {
                        // LAPLACIAN MATRIX OPERATION
                        // wsp0 = u       = B   * u_hat
                        // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                        // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                        BwdTrans_SumFacKernel(base0,base1,inarray,wsp0,wsp1,true,true);
                        StdExpansion2D::PhysTensorDeriv(wsp0,wsp1,wsp2);
                    }
                    else
                    {
                        StdExpansion2D::PhysTensorDeriv(inarray,wsp1,wsp2);
                    }

                    // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
                    // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
                    // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
                    // especially for this purpose
                    const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();
                    int    dim = m_geom->GetCoordim();
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        // wsp3 = g0*g0 + g2*g2
                        Vmath::Vvtvvtp(nqtot,&gmat[0][0],1,&gmat[0][0],1,&gmat[2][0],1,&gmat[2][0],1,&wsp3[0],1);
                        // wsp4 = g0*g1 + g2*g3;
                        Vmath::Vvtvvtp(nqtot,&gmat[0][0],1,&gmat[1][0],1,&gmat[2][0],1,&gmat[3][0],1,&wsp4[0],1);
                        // wsp5 = g1*g1 + g3*g3;
                        Vmath::Vvtvvtp(nqtot,&gmat[1][0],1,&gmat[1][0],1,&gmat[3][0],1,&gmat[3][0],1,&wsp5[0],1);

                        // If 3D coordinates, tag on extra terms
                        if(dim == 3)
                        {
                            // wsp3 += g4*g4
                            Vmath::Vvtvp(nqtot,&gmat[4][0],1,&gmat[4][0],1,&wsp3[0],1,&wsp3[0],1);
                            // wsp4 += g4*g5
                            Vmath::Vvtvp(nqtot,&gmat[4][0],1,&gmat[5][0],1,&wsp4[0],1,&wsp4[0],1);
                            // wsp5 += g5*g5
                            Vmath::Vvtvp(nqtot,&gmat[5][0],1,&gmat[5][0],1,&wsp5[0],1,&wsp5[0],1);
                        }

                        Vmath::Vvtvvtp(nqtot,&wsp3[0],1,&wsp1[0],1,&wsp4[0],1,&wsp2[0],1,&wsp0[0],1);
                        Vmath::Vvtvvtp(nqtot,&wsp4[0],1,&wsp1[0],1,&wsp5[0],1,&wsp2[0],1,&wsp2[0],1);
                    }
                    else
                    {
                        NekDouble g0 = gmat[0][0]*gmat[0][0] + gmat[2][0]*gmat[2][0];
                        NekDouble g1 = gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0];
                        NekDouble g2 = gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0];

                        if(dim == 3)
                        {
                            g0 += gmat[4][0]*gmat[4][0];
                            g1 += gmat[4][0]*gmat[5][0];
                            g2 += gmat[5][0]*gmat[5][0];
                        }

                        if(fabs(g1) < NekConstants::kGeomFactorsTol)
                        {
                            Vmath::Smul(nqtot,g0,&wsp1[0],1,&wsp0[0],1);
                            Vmath::Smul(nqtot,g2,&wsp2[0],1,&wsp2[0],1);
                        }
                        else
                        {
                            Vmath::Svtsvtp(nqtot,g0,&wsp1[0],1,g1,&wsp2[0],1,&wsp0[0],1);
                            Vmath::Svtsvtp(nqtot,g1,&wsp1[0],1,g2,&wsp2[0],1,&wsp2[0],1);
                        }
                    }

                    MultiplyByQuadratureMetric(wsp0,wsp0);
                    MultiplyByQuadratureMetric(wsp2,wsp2);

                    // outarray = m = (D_xi1 * B)^T * k
                    // wsp1     = n = (D_xi2 * B)^T * l
                    IProductWRTBase_SumFacKernel(dbase0,base1,wsp0,outarray,wsp1,false,true);
                    IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp1,    wsp0,true,false);

                    // outarray = outarray + wsp1
                    //          = L * u_hat
                    Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);
                }
            }
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }
        }

        /**
         * @param   inarray     Input array @f$ \mathbf{u} @f$.
         * @param   outarray    Output array @f$ \boldsymbol{\nabla^2u}
         *                          + \lambda \boldsymbol{u} @f$.
         * @param   mkey
         */
        void QuadExp::HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                      Array<OneD,NekDouble> &outarray,
                                                      const StdRegions::StdMatrixKey &mkey)
        {
            if(m_metricinfo->IsUsingLaplMetrics())
            {
                int       nquad0  = m_base[0]->GetNumPoints();
                int       nquad1  = m_base[1]->GetNumPoints();
                int       nqtot   = nquad0*nquad1;
                int       nmodes0 = m_base[0]->GetNumModes();
                int       nmodes1 = m_base[1]->GetNumModes();
                int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);
                NekDouble lambda  = mkey.GetConstant(0);

                const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
                const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
                const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
                const Array<TwoD, const NekDouble>& metric = m_metricinfo->GetLaplacianMetrics();

                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(4*wspsize);
                Array<OneD,NekDouble> wsp1(wsp0+wspsize);
                Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);
                Array<OneD,NekDouble> wsp3(wsp0+3*wspsize);

                if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                {
                    // MASS MATRIX OPERATION
                    // The following is being calculated:
                    // wsp0     = B   * u_hat = u
                    // wsp1     = W   * wsp0
                    // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat
                    BwdTrans_SumFacKernel       (base0,base1,inarray,wsp0,    wsp1,true,true);
                    MultiplyByQuadratureMetric  (wsp0,wsp2);
                    IProductWRTBase_SumFacKernel(base0,base1,wsp2,   outarray,wsp1,true,true);

                    // LAPLACIAN MATRIX OPERATION
                    // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                    // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                    StdExpansion2D::PhysTensorDeriv(wsp0,wsp1,wsp2);
                }
                else
                {
                    // specialised implementation for the classical spectral element method
                    StdExpansion2D::PhysTensorDeriv(inarray,wsp1,wsp2);
                    MultiplyByQuadratureMetric(inarray,outarray);
                }

                // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
                // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g1 * du_dxi1 + g2 * du_dxi2
                // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
                // especially for this purpose
                if(!m_metricinfo->LaplacianMetricIsZero(1))
                {
                    Vmath::Vvtvvtp(nqtot,&metric[0][0],1,&wsp1[0],1,&metric[1][0],1,&wsp2[0],1,&wsp0[0],1);
                    Vmath::Vvtvvtp(nqtot,&metric[1][0],1,&wsp1[0],1,&metric[2][0],1,&wsp2[0],1,&wsp2[0],1);
                }
                else
                {
                    // special implementation in case g1 = 0 (which should hold for undistorted quads)
                    // wsp0 = k = g0 * wsp1 = g0 * du_dxi1
                    // wsp2 = l = g2 * wsp2 = g2 * du_dxi2
                    Vmath::Vmul(nqtot,&metric[0][0],1,&wsp1[0],1,&wsp0[0],1);
                    Vmath::Vmul(nqtot,&metric[2][0],1,&wsp2[0],1,&wsp2[0],1);
                }

                // wsp1 = m = (D_xi1 * B)^T * k
                // wsp0 = n = (D_xi2 * B)^T * l
                IProductWRTBase_SumFacKernel(dbase0,base1,wsp0,wsp1,wsp3,false,true);
                IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp0,wsp3,true,false);

                // outarray = lambda * outarray + (wsp0 + wsp1)
                //          = (lambda * M + L ) * u_hat
                Vmath::Vstvpp(m_ncoeffs,lambda,&outarray[0],1,&wsp1[0],1,&wsp0[0],1,&outarray[0],1);
            }
            else
            {
                int       nquad0  = m_base[0]->GetNumPoints();
                int       nquad1  = m_base[1]->GetNumPoints();
                int       nqtot   = nquad0*nquad1;
                int       nmodes0 = m_base[0]->GetNumModes();
                int       nmodes1 = m_base[1]->GetNumModes();
                int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);
                NekDouble lambda  = mkey.GetConstant(0);

                const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
                const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
                const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();

                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(6*wspsize);
                Array<OneD,NekDouble> wsp1(wsp0+wspsize);
                Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);
                Array<OneD,NekDouble> wsp3(wsp0+3*wspsize);
                Array<OneD,NekDouble> wsp4(wsp0+4*wspsize);
                Array<OneD,NekDouble> wsp5(wsp0+5*wspsize);

                if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                {
                    // MASS MATRIX OPERATION
                    // The following is being calculated:
                    // wsp0     = B   * u_hat = u
                    // wsp1     = W   * wsp0
                    // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat
                    BwdTrans_SumFacKernel       (base0,base1,inarray,wsp0,    wsp1,true,true);
                    MultiplyByQuadratureMetric  (wsp0,wsp2);
                    IProductWRTBase_SumFacKernel(base0,base1,wsp2,   outarray,wsp1,true,true);

                    // LAPLACIAN MATRIX OPERATION
                    // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                    // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                    StdExpansion2D::PhysTensorDeriv(wsp0,wsp1,wsp2);
                }
                else
                {
                    // specialised implementation for the classical spectral element method
                    StdExpansion2D::PhysTensorDeriv(inarray,wsp1,wsp2);
                    MultiplyByQuadratureMetric(inarray,outarray);
                }

                // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
                // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g1 * du_dxi1 + g2 * du_dxi2
                // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
                // especially for this purpose
                const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();
                int    dim = m_geom->GetCoordim();
                if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                {
                    // wsp3 = g0
                    Vmath::Vvtvvtp(nqtot,&gmat[0][0],1,&gmat[0][0],1,&gmat[2][0],1,&gmat[2][0],1,&wsp3[0],1);
                    // wsp4 = g1;
                    Vmath::Vvtvvtp(nqtot,&gmat[0][0],1,&gmat[1][0],1,&gmat[2][0],1,&gmat[3][0],1,&wsp4[0],1);
                    // wsp5 = g2;
                    Vmath::Vvtvvtp(nqtot,&gmat[1][0],1,&gmat[1][0],1,&gmat[3][0],1,&gmat[3][0],1,&wsp5[0],1);

                    if(dim == 3)
                    {
                        Vmath::Vvtvp(nqtot,&gmat[4][0],1,&gmat[4][0],1,&wsp3[0],1,&wsp3[0],1);
                        Vmath::Vvtvp(nqtot,&gmat[4][0],1,&gmat[5][0],1,&wsp4[0],1,&wsp4[0],1);
                        Vmath::Vvtvp(nqtot,&gmat[5][0],1,&gmat[5][0],1,&wsp5[0],1,&wsp5[0],1);
                    }

                    Vmath::Vvtvvtp(nqtot,&wsp3[0],1,&wsp1[0],1,&wsp4[0],1,&wsp2[0],1,&wsp0[0],1);
                    Vmath::Vvtvvtp(nqtot,&wsp4[0],1,&wsp1[0],1,&wsp5[0],1,&wsp2[0],1,&wsp2[0],1);
                }
                else
                {
                    NekDouble g0 = gmat[0][0]*gmat[0][0] + gmat[2][0]*gmat[2][0];
                    NekDouble g1 = gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0];
                    NekDouble g2 = gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0];

                    if(dim == 3)
                    {
                        g0 += gmat[4][0]*gmat[4][0];
                        g1 += gmat[4][0]*gmat[5][0];
                        g2 += gmat[5][0]*gmat[5][0];
                    }

                    if(fabs(g1) < NekConstants::kGeomFactorsTol)
                    {
                        Vmath::Smul(nqtot,g0,&wsp1[0],1,&wsp0[0],1);
                        Vmath::Smul(nqtot,g2,&wsp2[0],1,&wsp2[0],1);
                    }
                    else
                    {
                        Vmath::Svtsvtp(nqtot,g0,&wsp1[0],1,g1,&wsp2[0],1,&wsp0[0],1);
                        Vmath::Svtsvtp(nqtot,g1,&wsp1[0],1,g2,&wsp2[0],1,&wsp2[0],1);
                    }
                }

                MultiplyByQuadratureMetric(wsp0,wsp0);
                MultiplyByQuadratureMetric(wsp2,wsp2);

                // wsp1 = m = (D_xi1 * B)^T * k
                // wsp0 = n = (D_xi2 * B)^T * l
                IProductWRTBase_SumFacKernel(dbase0,base1,wsp0,wsp1,wsp3,false,true);
                IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp0,wsp3,true,false);

                // outarray = lambda * outarray + (wsp0 + wsp1)
                //          = (lambda * M + L ) * u_hat
                Vmath::Vstvpp(m_ncoeffs,lambda,&outarray[0],1,&wsp1[0],1,&wsp0[0],1,&outarray[0],1);
            }
        }

        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////

        /**
            \brief Calculate the derivative of the physical points

            For quadrilateral region can use the Tensor_Deriv function
            defined under StdExpansion.

        **/
        void QuadExp::PhysDeriv(const Array<OneD, const NekDouble> & inarray,
                                Array<OneD,NekDouble> &out_d0,
                                Array<OneD,NekDouble> &out_d1,
                                Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int     nqtot = nquad0*nquad1;
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> diff0(2*nqtot);
            Array<OneD,NekDouble> diff1(diff0+nqtot);

            StdQuadExp::PhysDeriv(inarray, diff0, diff1);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul  (nqtot,gmat[0],1,diff0,1, out_d0, 1);
                    Vmath::Vvtvp (nqtot,gmat[1],1,diff1,1, out_d0, 1,
                                  out_d0,1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Vmul  (nqtot,gmat[2],1,diff0,1, out_d1, 1);
                    Vmath::Vvtvp (nqtot,gmat[3],1,diff1,1, out_d1, 1,
                                  out_d1,1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Vmul  (nqtot,gmat[4],1,diff0,1, out_d2, 1);
                    Vmath::Vvtvp (nqtot,gmat[5],1,diff1,1, out_d2, 1,
                                  out_d2,1);
                }
            }
            else // regular geometry
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul (nqtot, gmat[0][0], diff0, 1, out_d0, 1);
                    Blas::Daxpy (nqtot, gmat[1][0], diff1, 1, out_d0, 1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Smul (nqtot, gmat[2][0], diff0, 1, out_d1, 1);
                    Blas::Daxpy (nqtot, gmat[3][0], diff1, 1, out_d1, 1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Smul (nqtot, gmat[4][0], diff0,1, out_d2, 1);
                    Blas::Daxpy (nqtot, gmat[5][0], diff1,1, out_d2, 1);
                }
            }
        }

        void QuadExp::PhysDeriv(const int dir,
                                const Array<OneD, const NekDouble>& inarray,
                                Array<OneD, NekDouble> &outarray)
        {
            switch(dir)
            {
            case 0:
                {
                    PhysDeriv(inarray, outarray, NullNekDouble1DArray, NullNekDouble1DArray);
                }
                break;
            case 1:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray, outarray, NullNekDouble1DArray);
                }
                break;
            case 2:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray, NullNekDouble1DArray, outarray);
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
        }

        // Physical Derivation along direction vector
        void QuadExp::PhysDirectionalDeriv(const Array<OneD, const NekDouble> & inarray,
                                           const Array<OneD, const NekDouble>& direction,
                                           Array<OneD,NekDouble> &out)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1;

            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD,NekDouble> diff0(2*nqtot);
            Array<OneD,NekDouble> diff1(diff0+nqtot);

            StdQuadExp::PhysDeriv(inarray, diff0, diff1);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Array<OneD, Array<OneD, NekDouble> > tangmat(2);

                // d/dx_v^s = v_x*ds/dx + v_y*ds/dy + v_z*dx/dz
                for (int i=0; i< 2; ++i)
                {
                    tangmat[i] = Array<OneD, NekDouble>(nqtot,0.0);
                    for (int k=0; k<(m_geom->GetCoordim()); ++k)
                    {
                        Vmath::Vvtvp(nqtot,&gmat[2*k+i][0],1,&direction[k*nqtot],1,&tangmat[i][0],1,&tangmat[i][0],1);
                    }
                }

                /// D_v = d/dx_v^s + d/dx_v^r
                if(out.num_elements())
                {
                    Vmath::Vmul  (nqtot,&tangmat[0][0],1,&diff0[0],1, &out[0], 1);
                    Vmath::Vvtvp (nqtot,&tangmat[1][0],1,&diff1[0],1, &out[0], 1, &out[0],1);
                }

            }
            else
            {
                ASSERTL1(m_metricinfo->GetGtype() == SpatialDomains::eDeformed,"Wrong route");
            }
        }

        /** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a (this)->m_coeffs

            Inputs:\n

            - \a inarray: array of physical quadrature points to be transformed

            Outputs:\n

            - (this)->_coeffs: updated array of expansion coefficients.

        */
        void QuadExp::FwdTrans(const Array<OneD, const NekDouble> & inarray,
                               Array<OneD,NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, m_coeffs, 1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass,
                                              DetExpansionType(),*this);
                DNekScalMatSharedPtr& matsys = m_matrixManager[masskey];

                // copy inarray in case inarray == outarray
                NekVector<const NekDouble> in(m_ncoeffs,outarray,eCopy);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;
            }
        }

        void QuadExp::FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray,
                                              Array<OneD, NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
            {
                int i,j;
                int npoints[2] = {m_base[0]->GetNumPoints(),
                                  m_base[1]->GetNumPoints()};
                int nmodes[2]  = {m_base[0]->GetNumModes(),
                                  m_base[1]->GetNumModes()};

                fill(outarray.get(), outarray.get()+m_ncoeffs, 0.0 );

                Array<OneD, NekDouble> physEdge[4];
                Array<OneD, NekDouble> coeffEdge[4];
                StdRegions::EdgeOrientation orient[4];
                for(i = 0; i < 4; i++)
                {
                    physEdge[i]  = Array<OneD, NekDouble>(npoints[i%2]);
                    coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i%2]);
                    orient[i]    = GetEorient(i);
                }

                for(i = 0; i < npoints[0]; i++)
                {
                    physEdge[0][i] = inarray[i];
                    physEdge[2][i] = inarray[npoints[0]*npoints[1]-1-i];
                }

                for(i = 0; i < npoints[1]; i++)
                {
                    physEdge[1][i] = inarray[npoints[0]-1+i*npoints[0]];
                    physEdge[3][i] = inarray[(npoints[1]-1)*npoints[0]-i*npoints[0]];
                }

                for(i = 0; i < 4; i++)
                {
                    if( orient[i] == StdRegions::eBackwards )
                    {
                        reverse( (physEdge[i]).get() , (physEdge[i]).get() + npoints[i%2] );
                    }
                }

                SegExpSharedPtr segexp[4];
                for(i = 0; i < 4; i++)
                {
                    segexp[i] = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(m_base[i%2]->GetBasisKey(),GetGeom2D()->GetEdge(i));
                }

                Array<OneD, unsigned int> mapArray;
                Array<OneD, int>          signArray;
                NekDouble sign;

                for(i = 0; i < 4; i++)
                {
                    segexp[i%2]->FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);

                    GetEdgeToElementMap(i,orient[i],mapArray,signArray);
                    for(j=0; j < nmodes[i%2]; j++)
                    {
                        sign = (NekDouble) signArray[j];
                        outarray[ mapArray[j] ] = sign * coeffEdge[i][j];
                    }
                }

                if (m_ncoeffs > 4) {
                    Array<OneD, NekDouble> tmp0(m_ncoeffs);
                    Array<OneD, NekDouble> tmp1(m_ncoeffs);

                    StdRegions::StdMatrixKey  stdmasskey(StdRegions::eMass,DetExpansionType(),*this);
                    MassMatrixOp(outarray,tmp0,stdmasskey);
                    IProductWRTBase(inarray,tmp1);

                    Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

                    // get Mass matrix inverse (only of interior DOF)
                    // use block (1,1) of the static condensed system
                    // note: this block alreay contains the inverse matrix
                    MatrixKey             masskey(StdRegions::eMass,DetExpansionType(),*this);
                    DNekScalMatSharedPtr  matsys = (m_staticCondMatrixManager[masskey])->GetBlock(1,1);

                    int nBoundaryDofs = NumBndryCoeffs();
                    int nInteriorDofs = m_ncoeffs - nBoundaryDofs;

                    Array<OneD, NekDouble> rhs(nInteriorDofs);
                    Array<OneD, NekDouble> result(nInteriorDofs);

                    GetInteriorMap(mapArray);

                    for(i = 0; i < nInteriorDofs; i++)
                    {
                        rhs[i] = tmp1[ mapArray[i] ];
                    }

                    Blas::Dgemv('N', nInteriorDofs, nInteriorDofs, matsys->Scale(), &((matsys->GetOwnedMatrix())->GetPtr())[0],
                                nInteriorDofs,rhs.get(),1,0.0,result.get(),1);

                    for(i = 0; i < nInteriorDofs; i++)
                    {
                        outarray[ mapArray[i] ] = result[i];
                    }
                }
            }

        }

//        void QuadExp::GetSurfaceNormal(Array<OneD,NekDouble> &SurfaceNormal,
//                                       const int k)
//	    {
//            int m_num = m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints();
//
//            Vmath::Vcopy(m_num, GetSurfaceNormal()[k], 1,
//                                SurfaceNormal, 1);
//      	}

        void QuadExp::GetCoords(Array<OneD,NekDouble> &coords_0,
                                Array<OneD,NekDouble> &coords_1,
                                Array<OneD,NekDouble> &coords_2)
        {
            LibUtilities::BasisSharedPtr CBasis0;
            LibUtilities::BasisSharedPtr CBasis1;
            Array<OneD,NekDouble>  x;

            ASSERTL0(m_geom, "m_geom not defined");

            // get physical points defined in Geom
            m_geom->FillGeom();

            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements() != 0,
                         "output coords_2 is not defined");

                CBasis0 = m_geom->GetBasis(2,0);
                CBasis1 = m_geom->GetBasis(2,1);

                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(2);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*
                                m_base[1]->GetNumPoints(),
                                x,1,coords_2,1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp2D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(),&(m_geom->UpdatePhys(2))[0], m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),&coords_2[0]);
                }
            case 2:
                ASSERTL0(coords_1.num_elements(),
                         "output coords_1 is not defined");

                CBasis0 = m_geom->GetBasis(1,0);
                CBasis1 = m_geom->GetBasis(1,1);

                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(1);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*
                                m_base[1]->GetNumPoints(),
                                x,1,coords_1,1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp2D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), &(m_geom->UpdatePhys(1))[0],
                                           m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),&coords_1[0]);
                }
            case 1:
                ASSERTL0(coords_0.num_elements(),
                         "output coords_0 is not defined");

                CBasis0 = m_geom->GetBasis(0,0);
                CBasis1 = m_geom->GetBasis(0,1);

                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(0);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*
                                m_base[1]->GetNumPoints(),
                                x,1,coords_0,1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp2D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), &(m_geom->UpdatePhys(0))[0],
                                           m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),&coords_0[0]);
                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions are greater than 2");
                break;
            }
        }

        // get the coordinates "coords" at the local coordinates "Lcoords"

        void QuadExp::GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                               Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[1] <= 1.0 &&
                     Lcoords[1] >= -1.0 && Lcoords[1]  <=1.0,
                     "Local coordinates are not in region [-1,1]");

            m_geom->FillGeom();
            for(i = 0; i < m_geom->GetCoordDim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }

        void QuadExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
        {
            if(format==eTecplot)
            {
                int i,j;
                int nquad0 = m_base[0]->GetNumPoints();
                int nquad1 = m_base[1]->GetNumPoints();
                Array<OneD,NekDouble> coords[3];

                ASSERTL0(m_geom,"m_geom not defined");

                int     coordim  = m_geom->GetCoordim();

                coords[0] = Array<OneD,NekDouble>(nquad0*nquad1);
                coords[1] = Array<OneD,NekDouble>(nquad0*nquad1);
                coords[2] = Array<OneD,NekDouble>(nquad0*nquad1);

                GetCoords(coords[0],coords[1],coords[2]);

                if(dumpVar)
                {
                    outfile << "Variables = x";

                    if(coordim == 2)
                    {
                        outfile << ", y";
                    }
                    else if (coordim == 3)
                    {
                        outfile << ", y, z";
                    }
                    outfile << ", "<< var << std::endl << std::endl;
                }

                outfile << "Zone, I=" << nquad0 << ", J=" <<
                    nquad1 <<", F=Point" << std::endl;

                for(i = 0; i < nquad0*nquad1; ++i)
                {
                    for(j = 0; j < coordim; ++j)
                    {
                        outfile << coords[j][i] << " ";
                    }
                    outfile << m_phys[i] << std::endl;
                }
            }
            else if(format==eGmsh)
            {
                if(dumpVar)
                {
                    outfile<<"View.MaxRecursionLevel = 4;"<<endl;
                    outfile<<"View.TargetError = 0.00;"<<endl;
                    outfile<<"View.AdaptVisualizationGrid = 1;"<<endl;
                    outfile<<"View \" \" {"<<endl;
                }
                outfile<<"SQ("<<endl;
                // write the coordinates of the vertices of the quadrilateral
                Array<OneD,NekDouble> coordVert1(2);
                Array<OneD,NekDouble> coordVert2(2);
                Array<OneD,NekDouble> coordVert3(2);
                Array<OneD,NekDouble> coordVert4(2);
                coordVert1[0]=-1.0;
                coordVert1[1]=-1.0;
                coordVert2[0]=1.0;
                coordVert2[1]=-1.0;
                coordVert3[0]=1.0;
                coordVert3[1]=1.0;
                coordVert4[0]=-1.0;
                coordVert4[1]=1.0;
                outfile<<m_geom->GetCoord(0,coordVert1)<<", ";
                outfile<<m_geom->GetCoord(1,coordVert1)<<", 0.0,"<<endl;
                outfile<<m_geom->GetCoord(0,coordVert2)<<", ";
                outfile<<m_geom->GetCoord(1,coordVert2)<<", 0.0,"<<endl;
                outfile<<m_geom->GetCoord(0,coordVert3)<<", ";
                outfile<<m_geom->GetCoord(1,coordVert3)<<", 0.0,"<<endl;
                outfile<<m_geom->GetCoord(0,coordVert4)<<", ";
                outfile<<m_geom->GetCoord(1,coordVert4)<<", 0.0"<<endl;
                outfile<<")"<<endl;

                // calculate the coefficients (monomial format)
                int i,j;

                int nModes0 = m_base[0]->GetNumModes();
                int nModes1 = m_base[1]->GetNumModes();

                const LibUtilities::PointsKey Pkey1Gmsh(nModes0,LibUtilities::eGaussLobattoLegendre);
                const LibUtilities::PointsKey Pkey2Gmsh(nModes1,LibUtilities::eGaussLobattoLegendre);
                const LibUtilities::BasisKey  Bkey1Gmsh(m_base[0]->GetBasisType(),nModes0,Pkey1Gmsh);
                const LibUtilities::BasisKey  Bkey2Gmsh(m_base[1]->GetBasisType(),nModes1,Pkey2Gmsh);

                StdRegions::StdQuadExpSharedPtr EGmsh;
                EGmsh = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(Bkey1Gmsh,Bkey2Gmsh);

                int nMonomialPolynomials = EGmsh->GetNcoeffs();

                Array<OneD,NekDouble> xi1(nMonomialPolynomials);
                Array<OneD,NekDouble> xi2(nMonomialPolynomials);

                Array<OneD,NekDouble> x(nMonomialPolynomials);
                Array<OneD,NekDouble> y(nMonomialPolynomials);

                EGmsh->GetCoords(xi1,xi2);

                for(i=0;i<nMonomialPolynomials;i++)
                {
                    x[i] = xi1[i];//0.5*(1.0+xi1[i]);
                    y[i] = xi2[i];//0.5*(1.0+xi2[i]);
                }

                int cnt  = 0;
                Array<TwoD, int> exponentMap(nMonomialPolynomials,3,0);
                for(i = 0; i < nModes1; i++)
                {
                    for(j = 0; j < nModes0; j++)
                    {
                        exponentMap[cnt][0] = j;
                        exponentMap[cnt++][1] = i;
                    }
                }

                NekMatrix<NekDouble> vdm(nMonomialPolynomials,nMonomialPolynomials);
                for(i = 0 ; i < nMonomialPolynomials; i++)
                {
                    for(j = 0 ; j < nMonomialPolynomials; j++)
                    {
                        vdm(i,j) = pow(x[i],exponentMap[j][0])*pow(y[i],exponentMap[j][1]);
                    }
                }

                vdm.Invert();

                Array<OneD,NekDouble> rhs(nMonomialPolynomials);
                EGmsh->BwdTrans(m_coeffs,rhs);

                NekVector<const NekDouble> in(nMonomialPolynomials,rhs,eWrapper);
                NekVector<NekDouble> out(nMonomialPolynomials);
                out = vdm*in;

                //write the coefficients
                outfile<<"{";
                for(i = 0; i < nMonomialPolynomials; i++)
                {
                    outfile<<out[i];
                    if(i < nMonomialPolynomials - 1)
                    {
                        outfile<<", ";
                    }
                }
                outfile<<"};"<<endl;

                if(dumpVar)
                {
                    outfile<<"INTERPOLATION_SCHEME"<<endl;
                    outfile<<"{"<<endl;
                    for(i=0; i < nMonomialPolynomials; i++)
                    {
                        outfile<<"{";
                        for(j = 0; j < nMonomialPolynomials; j++)
                        {
                            if(i==j)
                            {
                                outfile<<"1.00";
                            }
                            else
                            {
                                outfile<<"0.00";
                            }
                            if(j < nMonomialPolynomials - 1)
                            {
                                outfile<<", ";
                            }
                        }
                        if(i < nMonomialPolynomials - 1)
                        {
                            outfile<<"},"<<endl;
                        }
                        else
                        {
                            outfile<<"}"<<endl<<"}"<<endl;
                        }
                    }

                    outfile<<"{"<<endl;
                    for(i=0; i < nMonomialPolynomials; i++)
                    {
                        outfile<<"{";
                        for(j = 0; j < 3; j++)
                        {
                            outfile<<exponentMap[i][j];
                            if(j < 2)
                            {
                                outfile<<", ";
                            }
                        }
                        if(i < nMonomialPolynomials - 1)
                        {
                            outfile<<"},"<<endl;
                        }
                        else
                        {
                            outfile<<"}"<<endl<<"};"<<endl;
                        }
                    }
                    outfile<<"};"<<endl;
                }
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }

        DNekMatSharedPtr QuadExp::CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            StdRegions::StdQuadExpSharedPtr tmp = MemoryManager<StdQuadExp>::AllocateSharedPtr(bkey0,bkey1);
            return tmp->GetStdMatrix(mkey);
        }

        NekDouble QuadExp::PhysEvaluate(const Array<OneD, const NekDouble> &coord)
        {
            PhysEvaluate(coord,m_phys);
        }

        NekDouble QuadExp::PhysEvaluate(const Array<OneD, const NekDouble> &coord, const Array<OneD, const NekDouble> & physvals)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(2);
            
            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);

            return StdQuadExp::PhysEvaluate(Lcoord, physvals);
        }

        DNekMatSharedPtr QuadExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            DNekMatSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
            case StdRegions::eHybridDGLamToU:
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
            case StdRegions::eHybridDGHelmBndLam:
                returnval = Expansion2D::v_GenMatrix(mkey);
                break;
            default:
                returnval = StdQuadExp::v_GenMatrix(mkey);
            }

            return returnval;
        }


        DNekScalMatSharedPtr QuadExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if((m_metricinfo->GetGtype() == SpatialDomains::eDeformed)||
                       (mkey.GetNvariableCoefficients()))
                    {
                        NekDouble        one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble        jac = (m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,mat);
                    }
                }
                break;
            case StdRegions::eInvMass:
                {
                    if((m_metricinfo->GetGtype() == SpatialDomains::eDeformed)||
                       (mkey.GetNvariableCoefficients()))
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetExpansionType(),
                                                         *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble fac = 1.0/(m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                    }
                }
                break;
            case StdRegions::eWeakDeriv0:
            case StdRegions::eWeakDeriv1:
            case StdRegions::eWeakDeriv2:
                {
                    if((m_metricinfo->GetGtype() == SpatialDomains::eDeformed)||
                       (mkey.GetNvariableCoefficients()))
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
                        int dir;

                        switch(mkey.GetMatrixType())
                        {
                        case StdRegions::eWeakDeriv0:
                            dir = 0;
                            break;
                        case StdRegions::eWeakDeriv1:
                            dir = 1;
                            break;
                        case StdRegions::eWeakDeriv2:
                            dir = 2;
                            break;
                        }

                        MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                            mkey.GetExpansionType(), *this);
                        MatrixKey deriv1key(StdRegions::eWeakDeriv1,
                                            mkey.GetExpansionType(), *this);

                        DNekMat &deriv0 = *GetStdMatrix(deriv0key);
                        DNekMat &deriv1 = *GetStdMatrix(deriv1key);

                        int rows = deriv0.GetRows();
                        int cols = deriv1.GetColumns();

                        DNekMatSharedPtr WeakDeriv = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                        (*WeakDeriv) = gmat[2*dir][0]*deriv0 + gmat[2*dir+1][0]*deriv1;
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,WeakDeriv);
                    }
                }
                break;
          case StdRegions::eWeakDirectionalDeriv:
                {
                    int dim = m_geom->GetCoordim();
                    int nqtot   = (m_base[0]->GetNumPoints())*(m_base[1]->GetNumPoints());
                    int nvarcoeffs = mkey.GetNvariableCoefficients();

                    NekDouble jac = (m_metricinfo->GetJac())[0];
                    Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();

                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr WeakDirectionalDeriv = GenMatrix(mkey);
                        Array<OneD, Array<OneD, const NekDouble> > Weight(1+2*dim);
                        Array<OneD, NekDouble> tmp;

                        // Store tangential basis in Weighted[0-dim]
                        Weight[0] = mkey.GetVariableCoefficient(0);

                        // Store gmat info in Weight[dim+1]
                        for (int k=0; k < 2*dim; ++k)
                        {
                            tmp = Array<OneD, NekDouble>(nqtot);
                            Vmath::Vcopy(nqtot, &gmat[k][0], 1, &tmp[0], 1);
                            Weight[k+1] = tmp;
                        }

                        StdRegions::StdMatrixKey  stdmasskey(StdRegions::eMassLevelCurvature,DetExpansionType(),*this,Weight);
                        DNekMatSharedPtr MassLevelCurvaturemat = GenMatrix(stdmasskey);

                        (*WeakDirectionalDeriv) = (*WeakDirectionalDeriv) + (*MassLevelCurvaturemat);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,WeakDirectionalDeriv);
                    }
                    else
                    {
                        Array<OneD, Array<OneD, const NekDouble> > Cxi(1);
                        Array<OneD, Array<OneD, const NekDouble> > Ceta(1);

                        // Directional Forcing is applied
                        Array<OneD, NekDouble> Cxi_wk  = Array<OneD, NekDouble> (nqtot,0.0);
                        Array<OneD, NekDouble> Ceta_wk = Array<OneD, NekDouble> (nqtot,0.0);

                        // Cxi = tan_{xi_x} * d \xi/dx + tan_{xi_y} * d \xi/dy + tan_{xi_z} * d \xi/dz
                        // Ceta = tan_{eta_x} * d \eta/dx + tan_{xi_y} * d \xi/dy + tan_{xi_z} * d \xi/dz
                        for (int k=0; k<dim; ++k)
                        {
                            Vmath::Svtvp(nqtot,gmat[2*k][0],&(mkey.GetVariableCoefficient(0))[k*nqtot],1,
                                         &Cxi[0][0],1,&Cxi_wk[0],1);
                            Vmath::Svtvp(nqtot,gmat[2*k+1][0],&(mkey.GetVariableCoefficient(0))[k*nqtot],1,
                                         &Ceta[0][0],1,&Ceta_wk[0],1);
                        }

                        // Assign value to const NekDouble array for input to key
                        Cxi[0]  = Cxi_wk;
                        Ceta[0] = Ceta_wk;

                        // derivxi = Cxi * ( B * D_{\xi} *B^T )
                        // deriveta = Ceta * ( B * D_{\eta} * B^T )
                        MatrixKey derivxikey(StdRegions::eWeakDeriv0, mkey.GetExpansionType(), *this, Cxi);
                        MatrixKey derivetakey(StdRegions::eWeakDeriv1, mkey.GetExpansionType(), *this, Ceta);

                        DNekMat &derivxi = *GetStdMatrix(derivxikey);
                        DNekMat &deriveta = *GetStdMatrix(derivetakey);

                        int rows = derivxi.GetRows();
                        int cols = deriveta.GetColumns();

                        DNekMatSharedPtr WeakDirectionalDeriv = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                        // D_t = D_xi
                        (*WeakDirectionalDeriv) = derivxi + deriveta;

                        // Add Weighted Mass with (\grad \cdot u )
                        Array<OneD, Array<OneD, const NekDouble> > Weight(1+2*dim);
                        Array<OneD, NekDouble> Weight_wk;

                        // Store tangential basis in Weighted[0-dim]
                        Weight[0] = mkey.GetVariableCoefficient(0);

                        // Store gmat info in Weight[dim+1]
                        for (int k=0; k < 2*dim; ++k)
                        {
                            Weight_wk = Array<OneD, NekDouble>(gmat[k].num_elements());
                            // assign constant value to Weight_wk
                            Weight_wk[0] = gmat[k][0];
                            // Case array into const NekDouble Array.
                            Weight[k+1] = Weight_wk;
                        }

                        StdRegions::StdMatrixKey  stdmasskey(StdRegions::eMassLevelCurvature,DetExpansionType(),*this,Weight);
                        DNekMatSharedPtr MassLevelCurvaturemat = GetStdMatrix(stdmasskey);

                        (*WeakDirectionalDeriv) = (*WeakDirectionalDeriv) + (*MassLevelCurvaturemat);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,WeakDirectionalDeriv);
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if( (m_metricinfo->GetGtype() == SpatialDomains::eDeformed) ||
                        (mkey.GetNvariableLaplacianCoefficients() > 0) )
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetExpansionType(), *this);

                        DNekMat &lap00 = *GetStdMatrix(lap00key);
                        DNekMat &lap01 = *GetStdMatrix(lap01key);
                        DNekMat &lap11 = *GetStdMatrix(lap11key);

                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                        // Additional terms if Quad embedded in 3D coordinate system
                        if (m_geom->GetCoordDim() == 3)
                        {
                            (*lap) = (gmat[0][0]*gmat[0][0]+gmat[2][0]*gmat[2][0]+gmat[4][0]*gmat[4][0])*lap00 +
                                (gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0] + gmat[4][0]*gmat[5][0])*(lap01 + Transpose(lap01)) +
                                (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0] + gmat[5][0]*gmat[5][0])*lap11;
                        }
                        else {
                            (*lap) = (gmat[0][0]*gmat[0][0]+gmat[2][0]*gmat[2][0])*lap00 +
                                (gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0])*(lap01 + Transpose(lap01)) +
                                (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0])*lap11;
                        }

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,lap);
                    }
                }
                break;
            case StdRegions::eInvLaplacianWithUnityMean:
                {
                    NekDouble one = 1.0;
                    MatrixKey lapkey(StdRegions::eLaplacian,mkey.GetExpansionType(), *this);
                    DNekMatSharedPtr lmat = GenMatrix(lapkey);

                    // replace first column with inner product wrt 1                    
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmp(nq);
                    Array<OneD, NekDouble> outarray(m_ncoeffs);
                    Vmath::Fill(nq,one,tmp,1);
                    v_IProductWRTBase(tmp, outarray);

                    Vmath::Vcopy(m_ncoeffs,&outarray[0],1,
                                 &(lmat->GetPtr())[0],m_ncoeffs);
                    
                    lmat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,lmat); //Populate  matrix.
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetConstant(0);
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetExpansionType(), *this);
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetExpansionType(), *this);
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);
                }
                break;
            case StdRegions::eIProductWRTBase:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,mat);
                    }
                }
                break;
            case StdRegions::eIProductWRTDerivBase0:
            case StdRegions::eIProductWRTDerivBase1:
            case StdRegions::eIProductWRTDerivBase2:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();
                        int dir;

                        switch(mkey.GetMatrixType())
                        {
                        case StdRegions::eIProductWRTDerivBase0:
                            dir = 0;
                            break;
                        case StdRegions::eIProductWRTDerivBase1:
                            dir = 1;
                            break;
                        case StdRegions::eIProductWRTDerivBase2:
                            dir = 2;
                            break;
                        }

                        MatrixKey iProdDeriv0Key(StdRegions::eIProductWRTDerivBase0,
                                                 mkey.GetExpansionType(), *this);
                        MatrixKey iProdDeriv1Key(StdRegions::eIProductWRTDerivBase1,
                                                 mkey.GetExpansionType(), *this);

                        DNekMat &stdiprod0 = *GetStdMatrix(iProdDeriv0Key);
                        DNekMat &stdiprod1 = *GetStdMatrix(iProdDeriv0Key);

                        int rows = stdiprod0.GetRows();
                        int cols = stdiprod1.GetColumns();

                        DNekMatSharedPtr mat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                        (*mat) = gmat[2*dir][0]*stdiprod0 + gmat[2*dir+1][0]*stdiprod1;

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,mat);
                    }
                }
                break;
            case StdRegions::eInvHybridDGHelmholtz:
                {
                    NekDouble one = 1.0;

                    int nvarcoeffs = mkey.GetNvariableCoefficients();
                    Array<OneD, Array<OneD, const NekDouble> > varcoeffs(nvarcoeffs);

                    if(nvarcoeffs>0)
                    {
                        for(int j=0; j<nvarcoeffs; j++)
                        {
                            varcoeffs[j] = mkey.GetVariableCoefficient(j);
                        }
                    }

                    StdRegions::StdMatrixKey hkey(StdRegions::eHybridDGHelmholtz,
                                                  DetExpansionType(),*this,
                                                  mkey.GetConstant(0),
                                                  mkey.GetConstant(1),
                                                  varcoeffs);

                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            default:
                {
                    NekDouble        one = 1.0;
                    DNekMatSharedPtr mat = GenMatrix(mkey);

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            }

            return returnval;
        }


        // Get edge values from the 2D Phys space along an edge
        // following a counter clockwise edge convention for definition
        // of edgedir, Note that point distribution is given by QuadExp.
        void QuadExp::GetEdgePhysVals(const int edge,
                                      const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            Array<OneD,const NekDouble> e_tmp;

            StdRegions::EdgeOrientation edgedir = GetEorient(edge);
            switch(edge)
            {
            case 0:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad0,inarray,1,outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad0,e_tmp = inarray+(nquad0-1),-1,
                                 outarray,1);
                }

                break;
            case 1:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad1,e_tmp = inarray+(nquad0-1),nquad0,
                                 outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad1,e_tmp = inarray+(nquad0*nquad1-1),
                                 -nquad0, outarray,1);
                }
                break;
            case 2:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad0,e_tmp = inarray+(nquad0*nquad1-1),-1,
                                 outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad0,e_tmp = inarray+nquad0*(nquad1-1),1,
                                 outarray,1);
                }
                break;
            case 3:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad1,e_tmp = inarray + nquad0*(nquad1-1),
                                 -nquad0,outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad1,inarray,nquad0,outarray,1);
                }
                break;
            default:
                ASSERTL0(false,"edge value (< 3) is out of range");
                break;
            }
        }



        void QuadExp::GetEdgePhysVals(const int edge, const StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                      const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            Array<OneD,const NekDouble> e_tmp;
            Array<OneD,NekDouble>       outtmp(max(nquad0,nquad1));


            // get points in Cartesian orientation
            switch(edge)
            {
            case 0:
                Vmath::Vcopy(nquad0,inarray,1,outtmp,1);
                break;
            case 1:
                Vmath::Vcopy(nquad1,e_tmp = inarray+(nquad0-1),nquad0,outtmp,1);
                break;
            case 2:
                Vmath::Vcopy(nquad0,e_tmp = inarray+nquad0*(nquad1-1),1,
                             outtmp,1);
                break;
            case 3:
                Vmath::Vcopy(nquad1,inarray,nquad0,outtmp,1);
                break;
            default:
                ASSERTL0(false,"edge value (< 3) is out of range");
                break;
            }

            // Interpolate if required
            LibUtilities::Interp1D(m_base[edge%2]->GetPointsKey(),outtmp,
                                   EdgeExp->GetBasis(0)->GetPointsKey(),outarray);

            //Reverse data if necessary
            if(GetCartesianEorient(edge) == StdRegions::eBackwards)
            {
                Vmath::Reverse(EdgeExp->GetNumPoints(0),&outarray[0],1,
                               &outarray[0],1);
            }

        }


        DNekScalBlkMatSharedPtr QuadExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;
            unsigned int exp_size[] = {nbdry,nint};
            int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks,nblks,exp_size,exp_size); //Really need a constructor which takes Arrays
            NekDouble factor = 1.0;

            switch(mkey.GetMatrixType())
            {
                // this can only use stdregions statically condensed system for mass matrix 
            case StdRegions::eMass:
                if((m_metricinfo->GetGtype() == SpatialDomains::eDeformed)||(mkey.GetNvariableCoefficients()))
                {
                    factor = 1.0;
                    goto UseLocRegionsMatrix;
                }
                else
                {
                    factor = (m_metricinfo->GetJac())[0];
                    goto UseStdRegionsMatrix;
                }
                break;
            default: // use Deformed case for both regular and deformed geometries
                factor = 1.0;
                goto UseLocRegionsMatrix;
                break;
            UseStdRegionsMatrix:
                {
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekBlkMatSharedPtr& mat = GetStdStaticCondMatrix(mkey);
                    DNekScalMatSharedPtr Atmp;
                    DNekMatSharedPtr     Asubmat;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(0,0)));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Asubmat = mat->GetBlock(0,1)));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(1,0)));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,Asubmat = mat->GetBlock(1,1)));
                }
                break;
            UseLocRegionsMatrix:
                {
                    int i,j;
                    int cnt = 0;
                    int cnt2 = 0;
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekScalMat &mat = *GetLocMatrix(mkey);
                    DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
                    DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
                    DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
                    DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);

                    Array<OneD,unsigned int> bmap(nbdry);
                    Array<OneD,unsigned int> imap(nint);
                    GetBoundaryMap(bmap);
                    GetInteriorMap(imap);

                    for(i = 0; i < nbdry; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*A)(i,j) = mat(bmap[i],bmap[j]);
                        }

                        for(j = 0; j < nint; ++j)
                        {
                            (*B)(i,j) = mat(bmap[i],imap[j]);
                        }
                    }

                    for(i = 0; i < nint; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*C)(i,j) = mat(imap[i],bmap[j]);
                        }

                        for(j = 0; j < nint; ++j)
                        {
                            (*D)(i,j) = mat(imap[i],imap[j]);
                        }
                    }

                    // Calculate static condensed system
                    if(nint)
                    {
                        D->Invert();
                        (*B) = (*B)*(*D);
                        (*A) = (*A) - (*B)*(*C);
                    }

                    DNekScalMatSharedPtr     Atmp;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,A));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,C));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,D));

                }
            }
            return returnval;
        }

        // Unpack data from input file assuming it comes from the same expansion type
        void QuadExp::v_ExtractDataToCoeffs(const std::vector<NekDouble> &data, 
                                            const int offset, 
                                            const std::vector<unsigned int > &nummodes, 
                                            const int nmode_offset,
                                            Array<OneD, NekDouble> &coeffs)
        {
            int data_order0 = nummodes[nmode_offset];
            int fillorder0  = std::min(m_base[0]->GetNumModes(),data_order0);

            int data_order1 = nummodes[nmode_offset+1];
            int order1      = m_base[1]->GetNumModes();
            int fillorder1  = min(order1,data_order1);
            
            switch(m_base[0]->GetBasisType())
            { 
            case LibUtilities::eModified_A:
                {
                    int i;
                    int cnt = 0;
                    int cnt1 = 0;

                    ASSERTL1(m_base[1]->GetBasisType() == LibUtilities::eModified_A,
                             "Extraction routine not set up for this basis");

                    Vmath::Zero(m_ncoeffs,coeffs,1);
                    for(i = 0; i < fillorder0; ++i)
                    {
                        Vmath::Vcopy(fillorder1,&data[offset+cnt],1,&coeffs[cnt1],1);
                        cnt  += data_order1;
                        cnt1 += order1;
                    }
                }
                break;
            default:
                ASSERTL0(false,"basis is either not set up or not hierarchicial");
            }
        }

        void QuadExp::v_ComputeEdgeNormal(const int edge)
        {
            int i;
            const SpatialDomains::GeomFactorsSharedPtr & geomFactors = GetGeom()->GetMetricInfo();
            SpatialDomains::GeomType type = geomFactors->GetGtype();
            const Array<TwoD, const NekDouble> & gmat = geomFactors->GetGmat();
            const Array<OneD, const NekDouble> & jac  = geomFactors->GetJac();
            int nqe = m_base[0]->GetNumPoints();
            int vCoordDim = GetCoordim();

            m_edgeNormals[edge] = Array<OneD, Array<OneD, NekDouble> >(vCoordDim);
            Array<OneD, Array<OneD, NekDouble> > &normal = m_edgeNormals[edge];
            for (i = 0; i < vCoordDim; ++i)
            {
                normal[i] = Array<OneD, NekDouble>(nqe);
            }

            // Regular geometry case
            if((type == SpatialDomains::eRegular)||(type == SpatialDomains::eMovingRegular))
            {
                NekDouble fac;
                // Set up normals
                switch(edge)
                {
                case 0:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,-gmat[2*i+1][0],normal[i],1);
                    }
                    break;
                case 1:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,gmat[2*i][0],normal[i],1);
                    }
                    break;
                case 2:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,gmat[2*i+1][0],normal[i],1);
                    }
                    break;
                case 3:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,-gmat[2*i][0],normal[i],1);
                    }
                    break;
                default:
                    ASSERTL0(false,"edge is out of range (edge < 4)");
                }

                // normalise
                fac = 0.0;
                for(i =0 ; i < vCoordDim; ++i)
                {
                    fac += normal[i][0]*normal[i][0];
                }
                fac = 1.0/sqrt(fac);
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Smul(nqe,fac,normal[i],1,normal[i],1);
                }
            }
            else   // Set up deformed normals
            {
                int j;

                int nquad0 = geomFactors->GetPointsKey(0).GetNumPoints();
                int nquad1 = geomFactors->GetPointsKey(1).GetNumPoints();

                LibUtilities::PointsKey from_key;

                Array<OneD,NekDouble> normals(vCoordDim*max(nquad0,nquad1),0.0);
                Array<OneD,NekDouble> edgejac(vCoordDim*max(nquad0,nquad1),0.0);

                // Extract Jacobian along edges and recover local
                // derivates (dx/dr) for polynomial interpolation by
                // multiplying m_gmat by jacobian
                switch(edge)
                {
                case 0:
                    for(j = 0; j < nquad0; ++j)
                    {
                        edgejac[j] = jac[j];
                        for(i = 0; i < vCoordDim; ++i)
                        {
                           normals[i*nquad0+j] = -gmat[2*i+1][j]*edgejac[j];
                        }
                   }
                    from_key = geomFactors->GetPointsKey(0);
                    break;
                case 1:
                    for(j = 0; j < nquad1; ++j)
                    {
                        edgejac[j] = jac[nquad0*j+nquad0-1];
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            normals[i*nquad1+j]  = gmat[2*i][nquad0*j + nquad0-1]*edgejac[j];
                        }
                    }
                    from_key = geomFactors->GetPointsKey(1);
                    break;
                case 2:
                    for(j = 0; j < nquad0; ++j)
                    {
                        edgejac[j] = jac[nquad0*(nquad1-1)+j];
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            normals[i*nquad0+j] = (gmat[2*i+1][nquad0*(nquad1-1)+j])*edgejac[j];
                        }
                    }
                    from_key = geomFactors->GetPointsKey(0);
                    break;
                case 3:
                    for(j = 0; j < nquad1; ++j)
                    {
                        edgejac[j] = jac[nquad0*j];
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            normals[i*nquad1+j] = -gmat[2*i][nquad0*j]*edgejac[j];
                        }
                    }
                    from_key = geomFactors->GetPointsKey(1);
                    break;
                default:
                    ASSERTL0(false,"edge is out of range (edge < 3)");
                }

                int nq  = from_key.GetNumPoints();
                Array<OneD,NekDouble> work(nqe,0.0);

                // interpolate Jacobian and invert
                LibUtilities::Interp1D(from_key,jac,m_base[0]->GetPointsKey(),work);
                Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);

                // interpolate
                for(i = 0; i < GetCoordim(); ++i)
                {
                    LibUtilities::Interp1D(from_key,&normals[i*nq],m_base[0]->GetPointsKey(),&normal[i][0]);
                    Vmath::Vmul(nqe,work,1,normal[i],1,normal[i],1);
                }

                //normalise normal vectors
                Vmath::Zero(nqe,work,1);
                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vvtvp(nqe,normal[i],1, normal[i],1,work,1,work,1);
                }

                Vmath::Vsqrt(nqe,work,1,work,1);
                Vmath::Sdiv(nqe,1.0,work,1,work,1);

                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vmul(nqe,normal[i],1,work,1,normal[i],1);
                }

                // Reverse direction so that points are in
                // anticlockwise direction if edge >=2
                if(edge >= 2)
                {
                    for(i = 0; i < GetCoordim(); ++i)
                    {
                        Vmath::Reverse(nqe,normal[i],1, normal[i],1);
                    }
                }
            }
            if(GetGeom()->GetEorient(edge) == StdRegions::eBackwards)
            {
                for(i = 0; i < vCoordDim; ++i)
                {
                    if(geomFactors->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Reverse(nqe, normal[i], 1, normal[i],1);
                    }
                    Vmath::Neg(nqe,normal[i],1);
                }
            }
        }
        
      void QuadExp::v_NormVectorIProductWRTBase(
                const Array<OneD, const NekDouble> &Fx,
                const Array<OneD, const NekDouble> &Fy,
                const Array<OneD, const NekDouble> &Fz,
                Array< OneD, NekDouble> &outarray,
                bool NegateNormal)
        {
            int nq = m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints();
            Array<OneD, NekDouble > Fn(nq);

            const Array<OneD, const Array<OneD, NekDouble> > &normals = GetSurfaceNormal(); // m_metricinfo->GetNormals();
            //for(int i=0; i<nq; i++)
            //{
            //cout<<"nx= "<<normals[0][i]<<"  ny="<<normals[1][i]<<"  nz="<<normals[2][i]<<endl;
            //}
            Vmath::Vmul (nq,&Fx[0],1,&normals[0][0], 1,&Fn[0],1);
            Vmath::Vvtvp(nq,&Fy[0],1,&normals[0][1],1,&Fn[0],1,&Fn[0],1);
            Vmath::Vvtvp(nq,&Fz[0],1,&normals[0][3],1,&Fn[0],1,&Fn[0],1);
            //cout<<"NegateNormal =  "<<NegateNormal<<endl;
            if(NegateNormal == true)
            {
                Vmath::Neg(nq,Fn,1);
            }

            IProductWRTBase(Fn,outarray);
	    }

    }//end of namespace
}//end of namespace
