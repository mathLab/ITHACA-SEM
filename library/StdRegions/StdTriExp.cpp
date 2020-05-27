///////////////////////////////////////////////////////////////////////////////
//
// File StdTriExp.cpp
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
// Description: Triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdSegExp.h>       // for StdSegExp, etc

using namespace std;

namespace Nektar
{
    namespace StdRegions
    {

        StdTriExp::StdTriExp()
        {
        }


        StdTriExp::StdTriExp(
            const LibUtilities::BasisKey &Ba,
            const LibUtilities::BasisKey &Bb) :
            StdExpansion (LibUtilities::StdTriData::getNumberOfCoefficients(
                              Ba.GetNumModes(),
                              Bb.GetNumModes()),
                          2,Ba,Bb),
            StdExpansion2D(LibUtilities::StdTriData::getNumberOfCoefficients(
                               Ba.GetNumModes(),
                               Bb.GetNumModes()),
                           Ba,Bb)
        {
            ASSERTL0(Ba.GetNumModes() <= Bb.GetNumModes(),
                     "order in 'a' direction is higher than order "
                     "in 'b' direction");
        }

        StdTriExp::StdTriExp(const StdTriExp &T):
            StdExpansion(T),
            StdExpansion2D(T)
        {
        }

        StdTriExp::~StdTriExp()
        {
        }

        //-------------------------------
        // Integration Methods
        //-------------------------------
        NekDouble StdTriExp::v_Integral(
            const Array<OneD, const NekDouble>& inarray)
        {
            int    i;
            int nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> w1_tmp(nquad1);

            Array<OneD, const NekDouble> w0 = m_base[0]->GetW();
            Array<OneD, const NekDouble> w1 = m_base[1]->GetW();
            Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();

            switch(m_base[1]->GetPointsType())
            {
            case LibUtilities::eGaussRadauMAlpha1Beta0: // (0,1) Jacobi Inner product
                {
                    Vmath::Smul(nquad1, 0.5, w1, 1, w1_tmp,1);
                    break;
                }
            default:
                {
                    // include jacobian factor on whatever coordinates are defined.
                    for(i = 0; i < nquad1; ++i)
                    {
                        w1_tmp[i] = 0.5*(1-z1[i])*w1[i];
                    }
                    break;
                }
            }

            return StdExpansion2D::Integral(inarray,w0,w1_tmp);
        }

        //-----------------------------
        // Differentiation Methods
        //-----------------------------

        /**
         * \brief Calculate the derivative of the physical points.
         *
         * \f$ \frac{\partial u}{\partial  x_1} =  \left .
         * \frac{2.0}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
         * \right |_{\eta_2}\f$
         *
         * \f$ \frac{\partial u}{\partial  x_2} =  \left .
         * \frac{1+\eta_1}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
         * \right |_{\eta_2}  + \left . \frac{\partial u}{\partial d\eta_2}
         * \right |_{\eta_1}  \f$
         */
        void StdTriExp::v_PhysDeriv(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& out_d0,
                  Array<OneD,       NekDouble>& out_d1,
                  Array<OneD,       NekDouble>& out_d2)
        {
            boost::ignore_unused(out_d2);

            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> wsp(std::max(nquad0, nquad1));

            const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();

            // set up geometric factor: 2/(1-z1)
            Vmath::Sadd(nquad1, -1.0, z1, 1, wsp, 1);
            Vmath::Sdiv(nquad1, -2.0, wsp, 1, wsp, 1);

            if (out_d0.size() > 0)
            {
                PhysTensorDeriv(inarray, out_d0, out_d1);

                for (i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,wsp[i],&out_d0[0]+i*nquad0,1);
                }

                // if no d1 required do not need to calculate both deriv
                if (out_d1.size() > 0)
                {
                    // set up geometric factor: (1_z0)/(1-z1)
                    Vmath::Sadd(nquad0, 1.0, z0, 1, wsp, 1);
                    Vmath::Smul(nquad0, 0.5, wsp, 1, wsp, 1);

                    for (i = 0; i < nquad1; ++i)
                    {
                        Vmath::Vvtvp(nquad0,&wsp[0],1,&out_d0[0]+i*nquad0,
                                     1,&out_d1[0]+i*nquad0,
                                     1,&out_d1[0]+i*nquad0,1);
                    }
                }
            }
            else if (out_d1.size() > 0)
            {
                Array<OneD, NekDouble> diff0(nquad0*nquad1);
                PhysTensorDeriv(inarray, diff0, out_d1);

                for (i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,wsp[i],&diff0[0]+i*nquad0,1);
                }

                Vmath::Sadd(nquad0, 1.0, z0, 1, wsp, 1);
                Vmath::Smul(nquad0, 0.5, wsp, 1, wsp, 1);

                for (i = 0; i < nquad1; ++i)
                {
                    Vmath::Vvtvp(nquad0,&wsp[0],1,&diff0[0]+i*nquad0,
                                 1,&out_d1[0]+i*nquad0,
                                 1,&out_d1[0]+i*nquad0,1);
                }
            }
        }

        void StdTriExp::v_PhysDeriv(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            switch(dir)
            {
                case 0:
                {
                    v_PhysDeriv(inarray, outarray, NullNekDouble1DArray);
                    break;
                }
                case 1:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray, outarray);
                    break;
                }
                default:
                {
                    ASSERTL1(false,"input dir is out of range");
                    break;
                }
            }
        }

        void StdTriExp::v_StdPhysDeriv(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& out_d0,
                  Array<OneD,       NekDouble>& out_d1,
                  Array<OneD,       NekDouble>& out_d2)
        {
            boost::ignore_unused(out_d2);
            StdTriExp::v_PhysDeriv(inarray, out_d0, out_d1);
        }

        void StdTriExp::v_StdPhysDeriv(
            const int dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTriExp::v_PhysDeriv(dir,inarray,outarray);
        }


        //---------------------------------------
        // Transforms
        //---------------------------------------

        /**
         * \brief Backward tranform for triangular elements
         *
         * @note 'q' (base[1]) runs fastest in this element.
         */
        void StdTriExp::v_BwdTrans(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            v_BwdTrans_SumFac(inarray,outarray);
        }


        void StdTriExp::v_BwdTrans_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            Array<OneD, NekDouble> wsp(m_base[0]->GetNumPoints()*
                                       m_base[1]->GetNumModes());

            BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                  m_base[1]->GetBdata(),
                                  inarray,outarray,wsp);
        }

        void StdTriExp::v_BwdTrans_SumFacKernel(
            const Array<OneD, const NekDouble>& base0,
            const Array<OneD, const NekDouble>& base1,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray,
                  Array<OneD,       NekDouble>& wsp,
                  bool                          doCheckCollDir0,
                  bool                          doCheckCollDir1)
        {
            boost::ignore_unused(doCheckCollDir0, doCheckCollDir1);

            int  i;
            int  mode;
            int  nquad0  = m_base[0]->GetNumPoints();
            int  nquad1  = m_base[1]->GetNumPoints();
            int  nmodes0 = m_base[0]->GetNumModes();
            int  nmodes1 = m_base[1]->GetNumModes();

            ASSERTL1(wsp.size() >= nquad0*nmodes1,
                     "Workspace size is not sufficient");
            ASSERTL2((m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)||
                     (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                     "Basis[1] is not of general tensor type");

            for (i = mode = 0; i < nmodes0; ++i)
            {
                Blas::Dgemv('N', nquad1,nmodes1-i,1.0,base1.get()+mode*nquad1,
                            nquad1,&inarray[0]+mode,1,0.0,&wsp[0]+i*nquad1,1);
                mode += nmodes1-i;
            }

            // fix for modified basis by splitting top vertex mode
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                Blas::Daxpy(nquad1,inarray[1],base1.get()+nquad1,1,
                            &wsp[0]+nquad1,1);
            }

            Blas::Dgemm('N','T', nquad0,nquad1,nmodes0,1.0, base0.get(),nquad0,
                        &wsp[0], nquad1,0.0, &outarray[0], nquad0);
        }

        void StdTriExp::v_FwdTrans(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            NekVector<NekDouble> in(m_ncoeffs,outarray,eCopy);
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

            out = (*matsys)*in;
        }


        void StdTriExp::v_FwdTrans_BndConstrained(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int i,j;
            int npoints[2] = {m_base[0]->GetNumPoints(),
                              m_base[1]->GetNumPoints()};
            int nmodes[2]  = {m_base[0]->GetNumModes(),
                              m_base[1]->GetNumModes()};

            fill(outarray.get(), outarray.get()+m_ncoeffs, 0.0 );

            Array<OneD, NekDouble> physEdge[3];
            Array<OneD, NekDouble> coeffEdge[3];
            for(i = 0; i < 3; i++)
            {
                physEdge[i]  = Array<OneD, NekDouble>(npoints[i!=0]);
                coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i!=0]);
            }

            for(i = 0; i < npoints[0]; i++)
            {
                physEdge[0][i] = inarray[i];
            }

            for(i = 0; i < npoints[1]; i++)
            {
                physEdge[1][i] = inarray[npoints[0]-1+i*npoints[0]];
                physEdge[2][i] = inarray[(npoints[1]-1)*npoints[0]-i*npoints[0]];
            }

            StdSegExpSharedPtr segexp[2] = {
                MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(
                    m_base[0]->GetBasisKey()),
                MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(
                    m_base[1]->GetBasisKey())
            };

            Array<OneD, unsigned int> mapArray;
            Array<OneD, int>          signArray;
            NekDouble sign;

            for (i = 0; i < 3; i++)
            {
                //segexp[i!=0]->v_FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);
                segexp[i!=0]->FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);

                v_GetEdgeToElementMap(i,eForwards,mapArray,signArray);
                for (j = 0; j < nmodes[i != 0]; j++)
                {
                    sign = (NekDouble) signArray[j];
                    outarray[ mapArray[j] ] = sign * coeffEdge[i][j];
                }
            }

            Array<OneD, NekDouble> tmp0(m_ncoeffs);
            Array<OneD, NekDouble> tmp1(m_ncoeffs);

            StdMatrixKey      masskey(eMass,DetShapeType(),*this);
            MassMatrixOp(outarray,tmp0,masskey);
            v_IProductWRTBase(inarray,tmp1);

            Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

            // get Mass matrix inverse (only of interior DOF)
            // use block (1,1) of the static condensed system
            // note: this block alreay contains the inverse matrix
            DNekMatSharedPtr matsys =
                (m_stdStaticCondMatrixManager[masskey])->GetBlock(1,1);

            int nBoundaryDofs = v_NumBndryCoeffs();
            int nInteriorDofs = m_ncoeffs - nBoundaryDofs;

            Array<OneD, NekDouble> rhs   (nInteriorDofs);
            Array<OneD, NekDouble> result(nInteriorDofs);

            v_GetInteriorMap(mapArray);

            for (i = 0; i < nInteriorDofs; i++)
            {
                rhs[i] = tmp1[ mapArray[i] ];
            }

            Blas::Dgemv('N',nInteriorDofs,nInteriorDofs,
                        1.0,&(matsys->GetPtr())[0],nInteriorDofs,
                        rhs.get(),1,
                        0.0,result.get(),1);

            for (i = 0; i < nInteriorDofs; i++)
            {
                outarray[ mapArray[i] ] = result[i];
            }
        }

        //---------------------------------------
        // Inner product functions
        //---------------------------------------

        /**
         * \brief Calculate the inner product of inarray with respect to the
         * basis B=base0[p]*base1[pq] and put into outarray.
         *
         * \f$
         * \begin{array}{rcl}
         * I_{pq} = (\phi^A_q \phi^B_{pq}, u) &=&
         * \sum_{i=0}^{nq_0}\sum_{j=0}^{nq_1}
         * \phi^A_p(\eta_{0,i})\phi^B_{pq}(\eta_{1,j}) w^0_i w^1_j
         * u(\xi_{0,i} \xi_{1,j}) \\
         * & = & \sum_{i=0}^{nq_0} \phi^A_p(\eta_{0,i})
         * \sum_{j=0}^{nq_1} \phi^B_{pq}(\eta_{1,j}) \tilde{u}_{i,j}
         * \end{array}
         * \f$
         *
         * where
         *
         * \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
         *
         * which can be implemented as
         *
         * \f$  f_{pj} = \sum_{i=0}^{nq_0} \phi^A_p(\eta_{0,i})
         * \tilde{u}_{i,j}
         * \rightarrow {\bf B_1 U}  \f$
         * \f$  I_{pq} = \sum_{j=0}^{nq_1} \phi^B_{pq}(\eta_{1,j}) f_{pj}
         * \rightarrow {\bf B_2[p*skip] f[skip]}  \f$
         *
         * \b Recall: \f$ \eta_{1} = \frac{2(1+\xi_1)}{(1-\xi_2)}-1, \,
         * \eta_2 = \xi_2\f$
         *
         * \b Note: For the orthgonality of this expansion to be realised the
         * 'q' ordering must run fastest in contrast to the Quad and Hex
         * ordering where 'p' index runs fastest to be consistent with the
         * quadrature ordering.
         *
         * In the triangular space the i (i.e. \f$\eta_1\f$ direction)
         * ordering still runs fastest by convention.
         */
        void StdTriExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTriExp::v_IProductWRTBase_SumFac(inarray,outarray);
        }

        void StdTriExp::v_IProductWRTBase_MatOp(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void StdTriExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray,
            bool                                multiplybyweights)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();

            if(multiplybyweights)
            {
                Array<OneD,NekDouble> tmp(nquad0*nquad1+nquad1*order0);
                Array<OneD,NekDouble> wsp(tmp+nquad0*nquad1);

                // multiply by integration constants
                MultiplyByQuadratureMetric(inarray,tmp);
                IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                             m_base[1]->GetBdata(),
                                             tmp,outarray,wsp);
            }
            else
            {
                Array<OneD,NekDouble> wsp(nquad1*order0);
                IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                             m_base[1]->GetBdata(),
                                             inarray,outarray,wsp);
            }
        }

        void StdTriExp::v_IProductWRTBase_SumFacKernel(
            const Array<OneD, const NekDouble>& base0,
            const Array<OneD, const NekDouble>& base1,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray,
                  Array<OneD,       NekDouble>& wsp,
                  bool                          doCheckCollDir0,
                  bool                          doCheckCollDir1)
        {
            boost::ignore_unused(doCheckCollDir0, doCheckCollDir1);

            int i;
            int mode;
            int nquad0  = m_base[0]->GetNumPoints();
            int nquad1  = m_base[1]->GetNumPoints();
            int nmodes0 = m_base[0]->GetNumModes();
            int nmodes1 = m_base[1]->GetNumModes();

            ASSERTL1(wsp.size() >= nquad1*nmodes0,
                     "Workspace size is not sufficient");

            Blas::Dgemm('T','N',nquad1,nmodes0,nquad0,1.0,inarray.get(),nquad0,
                        base0.get(),nquad0,0.0,wsp.get(),nquad1);

            // Inner product with respect to 'b' direction
            for (mode=i=0; i < nmodes0; ++i)
            {
                Blas::Dgemv('T',nquad1,nmodes1-i,1.0, base1.get()+mode*nquad1,
                            nquad1,wsp.get()+i*nquad1,1, 0.0,
                            outarray.get() + mode,1);
                mode += nmodes1 - i;
            }

            // fix for modified basis by splitting top vertex mode
            if (m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                outarray[1] += Blas::Ddot(nquad1,base1.get()+nquad1,1,
                                          wsp.get()+nquad1,1);
            }
        }

        void StdTriExp::v_IProductWRTDerivBase(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTriExp::v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }

        void StdTriExp::v_IProductWRTDerivBase_MatOp(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int nq = GetTotPoints();
            MatrixType mtype = eIProductWRTDerivBase0;

            switch(dir)
            {
                case 0:
                {
                    mtype = eIProductWRTDerivBase0;
                    break;
                }
                case 1:
                {
                    mtype = eIProductWRTDerivBase1;
                    break;
                }
                default:
                {
                    ASSERTL1(false,"input dir is out of range");
                    break;
                }
            }

            StdMatrixKey      iprodmatkey(mtype,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void StdTriExp::v_IProductWRTDerivBase_SumFac(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int    i;
            int    nquad0  = m_base[0]->GetNumPoints();
            int    nquad1  = m_base[1]->GetNumPoints();
            int    nqtot   = nquad0*nquad1;
            int    nmodes0 = m_base[0]->GetNumModes();
            int    wspsize = max(max(nqtot,m_ncoeffs),nquad1*nmodes0);

            Array<OneD, NekDouble> gfac0(2*wspsize);
            Array<OneD, NekDouble> tmp0 (gfac0+wspsize);

            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();

            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac0[i] = 2.0/(1-z1[i]);
            }

            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Smul(nquad0,gfac0[i],&inarray[0]+i*nquad0,1,
                            &tmp0[0]+i*nquad0,1);
            }

            MultiplyByQuadratureMetric(tmp0,tmp0);

            switch(dir)
            {
                case 0:
                {
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 tmp0,outarray,gfac0);
                    break;
                }
                case 1:
                {
                    Array<OneD, NekDouble> tmp3(m_ncoeffs);
                    const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();

                    for (i = 0; i < nquad0; ++i)
                    {
                        gfac0[i] = 0.5*(1+z0[i]);
                    }

                    for (i = 0; i < nquad1; ++i)
                    {
                        Vmath::Vmul(nquad0,&gfac0[0],1,&tmp0[0]+i*nquad0,1,
                                    &tmp0[0]+i*nquad0,1);
                    }

                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 tmp0,tmp3,gfac0);

                    MultiplyByQuadratureMetric(inarray,tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetDbdata(),
                                                 tmp0,outarray,gfac0);
                    Vmath::Vadd(m_ncoeffs,&tmp3[0],1,&outarray[0],1,
                                &outarray[0],1);
                    break;
                }
                default:
                {
                    ASSERTL1(false, "input dir is out of range");
                    break;
                }
            }
        }

        //---------------------------------------
        // Evaluation functions
        //---------------------------------------

        void StdTriExp::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble>& xi,
                                                 Array<OneD, NekDouble>& eta)
        {

            // set up local coordinate system
            if (fabs(xi[1]-1.0) < NekConstants::kNekZeroTol)
            {
                eta[0] = -1.0;
                eta[1] =  1.0;
            }
            else
            {
                eta[0] = 2*(1+xi[0])/(1-xi[1])-1.0;
                eta[1] = xi[1];
            }
        }

        void StdTriExp::v_FillMode(
            const int mode, Array<OneD, NekDouble> &outarray)
        {
            int   i,m;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            int   order0 = m_base[0]->GetNumModes();
            int   order1 = m_base[1]->GetNumModes();
            int   mode0  = 0;
            Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();

            ASSERTL2(mode <= m_ncoeffs,
                     "calling argument mode is larger than "
                     "total expansion order");

            m = order1;
            for (i = 0; i < order0; ++i, m+=order1-i)
            {
                if (m > mode)
                {
                    mode0 = i;
                    break;
                }
            }

            // deal with top vertex mode in modified basis
            if (mode == 1 &&
                m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                Vmath::Fill(nquad0*nquad1 , 1.0, outarray, 1);
            }
            else
            {
                for (i = 0; i < nquad1; ++i)
                {
                    Vmath::Vcopy(nquad0,(NekDouble *)(base0.get()+mode0*nquad0),
                                 1,&outarray[0]+i*nquad0,1);
                }
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,(NekDouble *)(base1.get() + mode*nquad1),
                            1,&outarray[0]+i,nquad0,&outarray[0]+i,nquad0);
            }
        }


        int StdTriExp::v_GetNverts() const
        {
            return 3;
        }

        int StdTriExp::v_GetNedges() const
        {
            return 3;
        }

        LibUtilities::ShapeType StdTriExp::v_DetShapeType() const
        {
            return LibUtilities::eTriangle;
        }

        int StdTriExp::v_NumBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B,
                     "BasisType is not a boundary interior form");

            return 3 + (GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
        }

        int StdTriExp::v_NumDGBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B,
                     "BasisType is not a boundary interior form");

            return GetBasisNumModes(0) + 2*GetBasisNumModes(1);
        }

        int StdTriExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");

            if (i == 0)
            {
                return GetBasisNumModes(0);
            }
            else
            {
                return GetBasisNumModes(1);
            }
        }

        int StdTriExp::v_GetEdgeNumPoints(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 2),"edge id is out of range");

            if (i == 0)
            {
                return GetNumPoints(0);
            }
            else
            {
                return GetNumPoints(1);
            }
        }

        int StdTriExp::v_CalcNumberOfCoefficients(
            const std::vector<unsigned int> &nummodes,
            int                             &modes_offset)
        {
            int nmodes = LibUtilities::StdTriData::getNumberOfCoefficients(
                nummodes[modes_offset],
                nummodes[modes_offset+1]);
            modes_offset += 2;

            return nmodes;
        }

        LibUtilities::BasisType StdTriExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");

            if (i == 0)
            {
                return GetBasisType(0);
            }
            else
            {
                return GetBasisType(1);
            }
        }


        void StdTriExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                    Array<OneD, NekDouble> &coords_1,
                                    Array<OneD, NekDouble> &coords_2)
        {
            boost::ignore_unused(coords_2);

            Array<OneD, const NekDouble> z0 = m_base[0]->GetZ();
            Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();
            int nq0 = GetNumPoints(0);
            int nq1 = GetNumPoints(1);
            int i,j;

            for(i = 0; i < nq1; ++i)
            {
                for(j = 0; j < nq0; ++j)
                {
                    coords_0[i*nq0+j] = (1+z0[j])*(1-z1[i])/2.0 - 1.0;
                }
                Vmath::Fill(nq0,z1[i],&coords_1[0] + i*nq0,1);
            }
        }

        bool StdTriExp::v_IsBoundaryInteriorExpansion()
        {
            return m_base[0]->GetBasisType() == LibUtilities::eModified_A &&
                   m_base[1]->GetBasisType() == LibUtilities::eModified_B;
        }

        int StdTriExp::v_DetCartesianDirOfEdge(const int edge)
        {
            ASSERTL2(edge >= 0 && edge <= 2, "edge id is out of range");

            return edge == 0 ? 0 : 1;
        }

        const LibUtilities::BasisKey StdTriExp::v_DetEdgeBasisKey(
            const int i) const
        {
            ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");

            if (i == 0)
            {
                return GetBasis(0)->GetBasisKey();
            }
            else
            {
                switch(m_base[1]->GetBasisType())
                {
                case LibUtilities::eModified_B:
                {
                    switch(m_base[1]->GetPointsType())
                    {
                    case LibUtilities::eGaussRadauMAlpha1Beta0:
                    {
                        LibUtilities::PointsKey pkey(
                                m_base[1]->GetBasisKey().GetPointsKey().
                                GetNumPoints()+1,
                                LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                                LibUtilities::eModified_A,
                                m_base[1]->GetNumModes(),pkey);
                        break;
                    }

                    default:
                        NEKERROR(ErrorUtil::efatal,
                                 "unexpected points distribution");
                        break;
                    }
                    break;
                }
                default:
                    NEKERROR(ErrorUtil::efatal,
                             "Information not available to set edge key");
                    break;
                }
            }
            return LibUtilities::NullBasisKey;
        }



        //--------------------------
        // Mappings
        //--------------------------

        void StdTriExp::v_GetEdgeToElementMap(
            const int                  eid,
            const Orientation      edgeOrient,
            Array<OneD, unsigned int>& maparray,
            Array<OneD,          int>& signarray,
            int P)
        {
            ASSERTL1(GetEdgeBasisType(eid) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(eid) == LibUtilities::eModified_B,
                     "Mapping not defined for this type of basis");

            int i;
            int numModes=0;
            int order0 = m_base[0]->GetNumModes();
            int order1 = m_base[1]->GetNumModes();

            switch (eid)
            {
            case 0:
                numModes = order0;
                break;
            case 1:
            case 2:
                numModes = order1;
                break;
            default:
                ASSERTL0(false,"eid must be between 0 and 2");
            }

            bool checkForZeroedModes = false;
            if (P == -1)
            {
                P = numModes;
            }
            else if(P != numModes)
            {
                checkForZeroedModes = true;
            }


            if(maparray.size() != P)
            {
                maparray = Array<OneD, unsigned int>(P);
            }

            if(signarray.size() != P)
            {
                signarray = Array<OneD, int>(P,1);
            }
            else
            {
                fill(signarray.get() , signarray.get()+P, 1);
            }

            switch(eid)
            {
                case 0:
                {
                    int cnt = 0;
                    for(i = 0; i < P; cnt+=order1-i, ++i)
                    {
                        maparray[i] = cnt;
                    }

                    if(edgeOrient==eBackwards)
                    {
                        swap( maparray[0] , maparray[1] );

                        for(i = 3; i < P; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                }
                case 1:
                {
                    maparray[0] = order1;
                    maparray[1] = 1;
                    for(i = 2; i < P; i++)
                    {
                        maparray[i] = order1-1+i;
                    }

                    if(edgeOrient==eBackwards)
                    {
                        swap( maparray[0] , maparray[1] );

                        for(i = 3; i < P; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                }
                case 2:
                {
                    for(i = 0; i < P; i++)
                    {
                        maparray[i] = i;
                    }

                    if(edgeOrient==eBackwards)
                    {
                        swap( maparray[0] , maparray[1] );

                        for(i = 3; i < P; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                }
            default:
                ASSERTL0(false,"eid must be between 0 and 2");
                break;
            }


            if (checkForZeroedModes)
            {
                // Zero signmap and set maparray to zero if
                // elemental modes are not as large as face modes
                for (int j = numModes; j < P; j++)
                {
                    signarray[j] = 0.0;
                    maparray[j]  = maparray[0];
                }
            }
        }

        int StdTriExp::v_GetVertexMap(const int localVertexId,bool useCoeffPacking)
        {
            ASSERTL0(
                GetEdgeBasisType(localVertexId) == LibUtilities::eModified_A ||
                GetEdgeBasisType(localVertexId) == LibUtilities::eModified_B,
                "Mapping not defined for this type of basis");

            int localDOF = 0;
            if(useCoeffPacking == true)
            {
                switch(localVertexId)
                {
                case 0:
                    {
                        localDOF = 0;
                        break;
                    }
                case 1:
                    {
                        localDOF = 1;
                        break;
                    }
                case 2:
                    {
                        localDOF = m_base[1]->GetNumModes();
                        break;
                    }
                default:
                    {
                        ASSERTL0(false,"eid must be between 0 and 2");
                        break;
                    }
                }
            }
            else // follow book format for vertex indexing.
            {
                switch(localVertexId)
                {
                case 0:
                    {
                        localDOF = 0;
                        break;
                    }
                case 1:
                    {
                        localDOF = m_base[1]->GetNumModes();
                        break;
                    }
                case 2:
                    {
                        localDOF = 1;
                        break;
                    }
                default:
                    {
                        ASSERTL0(false,"eid must be between 0 and 2");
                        break;
                    }
                }
            }

            return localDOF;
        }

        void StdTriExp::v_GetEdgeInteriorMap(
            const int                  eid,
            const Orientation      edgeOrient,
            Array<OneD, unsigned int>& maparray,
            Array<OneD,          int>& signarray)
        {
            ASSERTL0(GetEdgeBasisType(eid) == LibUtilities::eModified_A||
                     GetEdgeBasisType(eid) == LibUtilities::eModified_B,
                     "Mapping not defined for this type of basis");
            int i;
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nEdgeIntCoeffs = GetEdgeNcoeffs(eid)-2;

            if(maparray.size() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }

            if(signarray.size() != nEdgeIntCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeIntCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nEdgeIntCoeffs, 1 );
            }

            switch(eid)
            {
                case 0:
                {
                    int cnt = 2*nummodes1 - 1;
                    for(i = 0; i < nEdgeIntCoeffs; cnt+=nummodes1-2-i, ++i)
                    {
                        maparray[i] = cnt;
                    }

                    if(edgeOrient==eBackwards)
                    {
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                }
                case 1:
                {
                    for(i = 0; i < nEdgeIntCoeffs; i++)
                    {
                        maparray[i] = nummodes1+1+i;
                    }

                    if(edgeOrient==eBackwards)
                    {
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                }
                case 2:
                {
                    for(i = 0; i < nEdgeIntCoeffs; i++)
                    {
                        maparray[i] = 2+i;
                    }

                    if(edgeOrient==eBackwards)
                    {
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                }
                default:
                {
                    ASSERTL0(false,"eid must be between 0 and 2");
                    break;
                }
            }
        }

        void StdTriExp::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                     GetBasisType(1) == LibUtilities::eModified_B,
                     "Expansion not of a proper type");

            int i,j;
            int cnt = 0;
            int nummodes0, nummodes1;
            int startvalue;
            if(outarray.size()!=GetNcoeffs()-NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(GetNcoeffs()-NumBndryCoeffs());
            }

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            startvalue = 2*nummodes1;

            for(i = 0; i < nummodes0-2; i++)
            {
                for(j = 0; j < nummodes1-3-i; j++)
                {
                    outarray[cnt++]=startvalue+j;
                }
                startvalue+=nummodes1-2-i;
            }
        }

        void StdTriExp::v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                     GetBasisType(1) == LibUtilities::eModified_B,
                     "Expansion not of expected type");
            int i;
            int cnt;
            int nummodes0, nummodes1;
            int value;

            if (outarray.size()!=NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
            }

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            value = 2*nummodes1-1;
            for(i = 0; i < value; i++)
            {
                outarray[i]=i;
            }
            cnt = value;

            for(i = 0; i < nummodes0-2; i++)
            {
                outarray[cnt++]=value;
                value += nummodes1-2-i;
            }
        }


        //---------------------------------------
        // Wrapper functions
        //---------------------------------------

        DNekMatSharedPtr StdTriExp::v_GenMatrix(const StdMatrixKey &mkey)
        {

            MatrixType mtype   = mkey.GetMatrixType();

            DNekMatSharedPtr Mat;

            switch(mtype)
            {
                case ePhysInterpToEquiSpaced:
                {
                    int nq0, nq1, nq;

                    nq0 = m_base[0]->GetNumPoints();
                    nq1 = m_base[1]->GetNumPoints();

                    // take definition from key
                    if(mkey.ConstFactorExists(eFactorConst))
                    {
                        nq = (int) mkey.GetConstFactor(eFactorConst);
                    }
                    else
                    {
                        nq = max(nq0,nq1);
                    }

                    int neq = LibUtilities::StdTriData::
                                                getNumberOfCoefficients(nq,nq);
                    Array<OneD, Array<OneD, NekDouble> > coords(neq);
                    Array<OneD, NekDouble>               coll  (2);
                    Array<OneD, DNekMatSharedPtr>        I     (2);
                    Array<OneD, NekDouble>               tmp   (nq0);


                    Mat = MemoryManager<DNekMat>::AllocateSharedPtr(neq,nq0*nq1);
                    int cnt = 0;

                    for(int i = 0; i < nq; ++i)
                    {
                        for(int j = 0; j < nq-i; ++j,++cnt)
                        {
                            coords[cnt]    = Array<OneD, NekDouble>(2);
                            coords[cnt][0] = -1.0 + 2*j/(NekDouble)(nq-1);
                            coords[cnt][1] = -1.0 + 2*i/(NekDouble)(nq-1);
                        }
                    }

                    for(int i = 0; i < neq; ++i)
                    {
                        LocCoordToLocCollapsed(coords[i],coll);

                        I[0] = m_base[0]->GetI(coll);
                        I[1] = m_base[1]->GetI(coll+1);

                        // interpolate first coordinate direction
                        for (int j  = 0; j < nq1; ++j)
                        {
                            NekDouble fac = (I[1]->GetPtr())[j];
                            Vmath::Smul(nq0, fac, I[0]->GetPtr(), 1, tmp, 1);

                            Vmath::Vcopy(nq0, &tmp[0], 1,
                                         Mat->GetRawPtr() + j*nq0*neq + i, neq);
                        }

                    }
                    break;
                }
                default:
                {
                    Mat = StdExpansion::CreateGeneralMatrix(mkey);
                    break;
                }
            }

            return Mat;
        }

        DNekMatSharedPtr StdTriExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return v_GenMatrix(mkey);
        }


        //---------------------------------------
        // Operator evaluation functions
        //---------------------------------------

        void StdTriExp::v_MassMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void StdTriExp::v_LaplacianMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdTriExp::v_LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void StdTriExp::v_LaplacianMatrixOp(
            const int                           k1,
            const int                           k2,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdExpansion::LaplacianMatrixOp_MatFree(
                k1,k2,inarray,outarray,mkey);
        }

        void StdTriExp::v_WeakDerivMatrixOp(
            const int                           i,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
        }

        void StdTriExp::v_HelmholtzMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdTriExp::v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
        }


        void StdTriExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                             const StdMatrixKey &mkey)
        {
            int qa = m_base[0]->GetNumPoints();
            int qb = m_base[1]->GetNumPoints();
            int nmodes_a = m_base[0]->GetNumModes();
            int nmodes_b = m_base[1]->GetNumModes();

            // Declare orthogonal basis.
            LibUtilities::PointsKey pa(qa,m_base[0]->GetPointsType());
            LibUtilities::PointsKey pb(qb,m_base[1]->GetPointsType());

            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A,nmodes_a,pa);
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B,nmodes_b,pb);
            StdTriExp OrthoExp(Ba,Bb);

            Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());

            // project onto physical space.
            OrthoExp.FwdTrans(array,orthocoeffs);

            if(mkey.ConstFactorExists(eFactorSVVPowerKerDiffCoeff)) // Rodrigo's power kern
            {
                NekDouble cutoff =  mkey.GetConstFactor(eFactorSVVCutoffRatio);
                NekDouble  SvvDiffCoeff  =
                    mkey.GetConstFactor(eFactorSVVPowerKerDiffCoeff)*
                    mkey.GetConstFactor(eFactorSVVDiffCoeff);

                int cnt = 0;
                for(int j = 0; j < nmodes_a; ++j)
                {
                    for(int k = 0; k < nmodes_b-j; ++k, ++cnt)
                    {
                        NekDouble fac = std::max(
                                    pow((1.0*j)/(nmodes_a-1),cutoff*nmodes_a),
                                    pow((1.0*k)/(nmodes_b-1),cutoff*nmodes_b));

                        orthocoeffs[cnt] *= (SvvDiffCoeff *fac);
                    }
                }
            }
            else if(mkey.ConstFactorExists(eFactorSVVDGKerDiffCoeff)) // Rodrigo/mansoor's DG kernel
            {
                NekDouble  SvvDiffCoeff  =
                    mkey.GetConstFactor(eFactorSVVDGKerDiffCoeff)*
                    mkey.GetConstFactor(eFactorSVVDiffCoeff);
                int max_ab = max(nmodes_a-kSVVDGFiltermodesmin,
                                 nmodes_b-kSVVDGFiltermodesmin);
                max_ab = max(max_ab,0);
                max_ab = min(max_ab,kSVVDGFiltermodesmax-kSVVDGFiltermodesmin);

                int cnt = 0;
                for(int j = 0; j < nmodes_a; ++j)
                {
                    for(int k = 0; k < nmodes_b-j; ++k, ++cnt)
                    {
                        int maxjk = max(j,k);
                        maxjk = min(maxjk,kSVVDGFiltermodesmax-1);

                        orthocoeffs[cnt] *= SvvDiffCoeff *
                            kSVVDGFilter[max_ab][maxjk];
                    }
                }
            }
            else
            {
                NekDouble  SvvDiffCoeff =
                    mkey.GetConstFactor(eFactorSVVDiffCoeff);

                int cutoff = (int) (mkey.GetConstFactor(eFactorSVVCutoffRatio)*
                                                        min(nmodes_a,nmodes_b));

                NekDouble epsilon = 1.0;
                int nmodes = min(nmodes_a,nmodes_b);

                int cnt = 0;

                // apply SVV filter (JEL)
                for(int j = 0; j < nmodes_a; ++j)
                {
                    for(int k = 0; k < nmodes_b-j; ++k)
                    {
                        if(j + k >= cutoff)
                        {
                            orthocoeffs[cnt] *= (SvvDiffCoeff
                                *exp(-(j+k-nmodes)*(j+k-nmodes)
                                    /((NekDouble)((j+k-cutoff+epsilon)
                                            *(j+k-cutoff+epsilon)))));
                        }
                        else
                        {
                            orthocoeffs[cnt] *= 0.0;
                        }
                        cnt++;
                    }
                }

            }

            // backward transform to physical space
            OrthoExp.BwdTrans(orthocoeffs,array);
        }

        void StdTriExp::v_ReduceOrderCoeffs(
            int                                 numMin,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            int n_coeffs = inarray.size();
            int nquad0   = m_base[0]->GetNumPoints();
            int nquad1   = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> coeff(n_coeffs);
            Array<OneD, NekDouble> coeff_tmp(n_coeffs,0.0);
            Array<OneD, NekDouble> tmp;
            Array<OneD, NekDouble> tmp2;
            int nqtot    = nquad0*nquad1;
            Array<OneD, NekDouble> phys_tmp(nqtot,0.0);

            int       nmodes0 = m_base[0]->GetNumModes();
            int       nmodes1 = m_base[1]->GetNumModes();
            int       numMin2 = nmodes0;
            int       i;

            const LibUtilities::PointsKey Pkey0(
                nmodes0, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey Pkey1(
                nmodes1, LibUtilities::eGaussLobattoLegendre);

            LibUtilities::BasisKey b0(m_base[0]->GetBasisType(),nmodes0,Pkey0);
            LibUtilities::BasisKey b1(m_base[1]->GetBasisType(),nmodes1,Pkey1);

            LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A,nmodes0,Pkey0);
            LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_B,nmodes1,Pkey1);

            StdRegions::StdTriExpSharedPtr m_OrthoTriExp;
            StdRegions::StdTriExpSharedPtr m_TriExp;

            m_TriExp      = MemoryManager<StdRegions::StdTriExp>
                ::AllocateSharedPtr(b0,      b1);
            m_OrthoTriExp = MemoryManager<StdRegions::StdTriExp>
                ::AllocateSharedPtr(bortho0, bortho1);

            m_TriExp     ->BwdTrans(inarray,phys_tmp);
            m_OrthoTriExp->FwdTrans(phys_tmp, coeff);

            for (i = 0; i < n_coeffs; i++)
            {
                if (i == numMin)
                {
                    coeff[i] = 0.0;
                    numMin += numMin2 - 1;
                    numMin2 -= 1.0;
                }
            }

            m_OrthoTriExp->BwdTrans(coeff,phys_tmp);
            m_TriExp     ->FwdTrans(phys_tmp, outarray);

        }

        void StdTriExp::v_GeneralMatrixOp_MatOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            DNekMatSharedPtr mat = m_stdMatrixManager[mkey];

            if(inarray.get() == outarray.get())
            {
                Array<OneD,NekDouble> tmp(m_ncoeffs);
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,tmp.get(),1);

                Blas::Dgemv('N', m_ncoeffs, m_ncoeffs, 1.0, mat->GetPtr().get(),
                            m_ncoeffs, tmp.get(), 1, 0.0, outarray.get(), 1);
            }
            else
            {
                Blas::Dgemv('N', m_ncoeffs, m_ncoeffs, 1.0, mat->GetPtr().get(),
                            m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
            }
        }

        //---------------------------------------
        // Private helper functions
        //---------------------------------------

        void StdTriExp::v_MultiplyByStdQuadratureMetric(
            const Array<OneD, const NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();

            const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();

            // multiply by integration constants
            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0,inarray.get()+i*nquad0,1,
                            w0.get(),1, outarray.get()+i*nquad0,1);
            }

            switch(m_base[1]->GetPointsType())
            {
                // Legendre inner product
                case LibUtilities::ePolyEvenlySpaced:
                case LibUtilities::eGaussLobattoLegendre:
                    for(i = 0; i < nquad1; ++i)
                    {
                        Blas::Dscal(nquad0,0.5*(1-z1[i])*w1[i],
                                    outarray.get()+i*nquad0,1);
                    }
                    break;

                // (1,0) Jacobi Inner product
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                    for(i = 0; i < nquad1; ++i)
                    {
                        Blas::Dscal(nquad0,0.5*w1[i],outarray.get()+i*nquad0,1);
                    }
                    break;

                default:
                    ASSERTL0(false, "Unsupported quadrature points type.");
                    break;
            }
        }

        void StdTriExp::v_GetSimplexEquiSpacedConnectivity(
            Array<OneD, int> &conn,
            bool              standard)
        {
            boost::ignore_unused(standard);

            int np1 = m_base[0]->GetNumPoints();
            int np2 = m_base[1]->GetNumPoints();
            int np = max(np1,np2);

            conn = Array<OneD, int>(3*(np-1)*(np-1));

            int row   = 0;
            int rowp1 = 0;
            int cnt = 0;
            for(int i = 0; i < np-1; ++i)
            {
                rowp1 += np-i;
                for(int j = 0; j < np-i-2; ++j)
                {
                    conn[cnt++] = row   +j;
                    conn[cnt++] = row   +j+1;
                    conn[cnt++] = rowp1 +j;

                    conn[cnt++] = rowp1 +j+1;
                    conn[cnt++] = rowp1 +j;
                    conn[cnt++] = row   +j+1;
                }

                conn[cnt++] = row  +np-i-2;
                conn[cnt++] = row  +np-i-1;
                conn[cnt++] = rowp1+np-i-2;

                row += np-i;
            }
        }


    }//end namespace
}//end namespace
