///////////////////////////////////////////////////////////////////////////////
//
// File StdQuadExp.cpp
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
// Description: Quadrilateral routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
    namespace StdRegions
    {


        StdQuadExp::StdQuadExp()
        {
        }

        /** \brief Constructor using BasisKey class for quadrature
         *  points and order definition
         */
        StdQuadExp::StdQuadExp(const LibUtilities::BasisKey &Ba,
                               const LibUtilities::BasisKey &Bb):
            StdExpansion  (Ba.GetNumModes()*Bb.GetNumModes(),2,Ba,Bb),
            StdExpansion2D(Ba.GetNumModes()*Bb.GetNumModes(),Ba,Bb)
        {
        }

        /** \brief Copy Constructor */
        StdQuadExp::StdQuadExp(const StdQuadExp &T):
            StdExpansion(T),
            StdExpansion2D(T)
        {
        }

        /** \brief Destructor */
        StdQuadExp::~StdQuadExp()
        {
        }

        /////////////////////////
        // Integration Methods //
        /////////////////////////

        NekDouble StdQuadExp::v_Integral(
                                  const Array<OneD, const NekDouble>& inarray)
        {
            Array<OneD, const NekDouble> w0 = m_base[0]->GetW();
            Array<OneD, const NekDouble> w1 = m_base[1]->GetW();

            return StdExpansion2D::Integral(inarray,w0,w1);
        }

        /////////////////////////////
        // Differentiation Methods //
        /////////////////////////////

        /** \brief Calculate the derivative of the physical points
         *
         *  For quadrilateral region can use the Tensor_Deriv function
         *  defined under StdExpansion.
         */

        void StdQuadExp::v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &out_d0,
                            Array<OneD, NekDouble> &out_d1,
                            Array<OneD, NekDouble> &out_d2)
        {
            boost::ignore_unused(out_d2);
            PhysTensorDeriv(inarray, out_d0, out_d1);
        }

        void StdQuadExp::v_PhysDeriv(const int dir,
                             const Array<OneD, const NekDouble>& inarray,
                             Array<OneD, NekDouble> &outarray)
        {
            switch(dir)
            {
            case 0:
                {
                    PhysTensorDeriv(inarray, outarray, NullNekDouble1DArray);
                }
                break;
            case 1:
                {
                    PhysTensorDeriv(inarray, NullNekDouble1DArray, outarray);
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
        }

        void StdQuadExp::v_StdPhysDeriv(const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &out_d0,
                            Array<OneD, NekDouble> &out_d1,
                            Array<OneD, NekDouble> &out_d2)
        {
            boost::ignore_unused(out_d2);
            //PhysTensorDeriv(inarray, out_d0, out_d1);
            StdQuadExp::v_PhysDeriv(inarray, out_d0, out_d1);
        }

        void StdQuadExp::v_StdPhysDeriv(const int dir,
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray)
        {
            //PhysTensorDeriv(inarray, outarray);
            StdQuadExp::v_PhysDeriv(dir,inarray,outarray);

        }


        ////////////////
        // Transforms //
        ////////////////

        void StdQuadExp::v_BwdTrans(const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray)
        {
            if(m_base[0]->Collocation() && m_base[1]->Collocation())
            {
                Vmath::Vcopy(m_base[0]->GetNumPoints() * m_base[1]->GetNumPoints(),
                               inarray, 1, outarray, 1);
            }
            else
            {
                StdQuadExp::v_BwdTrans_SumFac(inarray,outarray);
            }
        }

        void StdQuadExp::v_BwdTrans_SumFac(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> wsp(m_base[0]->GetNumPoints()*
                                       m_base[1]->GetNumModes());

            BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                  m_base[1]->GetBdata(),
                                  inarray,outarray,wsp,true,true);
        }

        // The arguments doCheckCollDir0 and doCheckCollDir1 allow you to specify whether
        // to check if the basis has collocation properties (i.e. for the classical spectral
        // element basis, In this case the 1D 'B' matrix is equal to the identity matrix
        // which can be exploited to speed up the calculations).
        // However, as this routine also allows to pass the matrix 'DB' (derivative of the basis),
        // the collocation property cannot always be used. Therefor follow this rule:
        // if base0 == m_base[0]->GetBdata() --> set doCheckCollDir0 == true;
        //    base1 == m_base[1]->GetBdata() --> set doCheckCollDir1 == true;
        //    base0 == m_base[0]->GetDbdata() --> set doCheckCollDir0 == false;
        //    base1 == m_base[1]->GetDbdata() --> set doCheckCollDir1 == false;
        void StdQuadExp::v_BwdTrans_SumFacKernel(
                            const Array<OneD, const NekDouble>& base0,
                            const Array<OneD, const NekDouble>& base1,
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray,
                            Array<OneD, NekDouble> &wsp,
                            bool doCheckCollDir0,
                            bool doCheckCollDir1)
        {
            int  nquad0  = m_base[0]->GetNumPoints();
            int  nquad1  = m_base[1]->GetNumPoints();
            int  nmodes0 = m_base[0]->GetNumModes();
            int  nmodes1 = m_base[1]->GetNumModes();

            bool colldir0 = doCheckCollDir0?(m_base[0]->Collocation()):false;
            bool colldir1 = doCheckCollDir1?(m_base[1]->Collocation()):false;

            if(colldir0 && colldir1)
            {
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,outarray.get(),1);
            }
            else if(colldir0)
            {
                Blas::Dgemm('N','T', nquad0, nquad1,nmodes1, 1.0, &inarray[0], nquad0,
                                base1.get(), nquad1, 0.0, &outarray[0], nquad0);
            }
            else if(colldir1)
            {
                Blas::Dgemm('N','N', nquad0,nmodes1,nmodes0,1.0, base0.get(),
                                nquad0, &inarray[0], nmodes0,0.0,&outarray[0], nquad0);
            }
            else
            {
                ASSERTL1(wsp.size()>=nquad0*nmodes1,"Workspace size is not sufficient");

                // Those two calls correpsond to the operation
                // out = B0*in*Transpose(B1);
                Blas::Dgemm('N','N', nquad0,nmodes1,nmodes0,1.0, base0.get(),
                                nquad0, &inarray[0], nmodes0,0.0,&wsp[0], nquad0);
                Blas::Dgemm('N','T', nquad0, nquad1,nmodes1, 1.0, &wsp[0], nquad0,
                                base1.get(), nquad1, 0.0, &outarray[0], nquad0);
            }
        }


        void StdQuadExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
            {
	        StdQuadExp::v_IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                StdMatrixKey     masskey(eInvMass,DetShapeType(),*this);
                DNekMatSharedPtr matsys = GetStdMatrix(masskey);

                // copy inarray in case inarray == outarray
                NekVector<NekDouble> in(m_ncoeffs,outarray,eCopy);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;
            }
        }

      void StdQuadExp::v_FwdTrans_BndConstrained(
                            const Array<OneD, const NekDouble>& inarray,
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
                for(i = 0; i < 4; i++)
                {
                    physEdge[i]  = Array<OneD, NekDouble>(npoints[i%2]);
                    coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i%2]);
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

                StdSegExpSharedPtr segexp[2] = {MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(m_base[0]->GetBasisKey()),
                                                MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(m_base[1]->GetBasisKey())};

                Array<OneD, unsigned int> mapArray;
                Array<OneD, int>          signArray;
                NekDouble sign;

                for(i = 0; i < 4; i++)
                {
                    segexp[i%2]->FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);

                    GetEdgeToElementMap(i,eForwards,mapArray,signArray);
                    for(j=0; j < nmodes[i%2]; j++)
                    {
                        sign = (NekDouble) signArray[j];
                        outarray[ mapArray[j] ] = sign * coeffEdge[i][j];
                    }
                }

                Array<OneD, NekDouble> tmp0(m_ncoeffs);
                Array<OneD, NekDouble> tmp1(m_ncoeffs);

                StdMatrixKey   masskey(eMass,DetShapeType(),*this);
                MassMatrixOp(outarray,tmp0,masskey);
                IProductWRTBase(inarray,tmp1);

                Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

                // get Mass matrix inverse (only of interior DOF)
                // use block (1,1) of the static condensed system
                // note: this block alreay contains the inverse matrix
                DNekMatSharedPtr  matsys = (m_stdStaticCondMatrixManager[masskey])->GetBlock(1,1);

                int nBoundaryDofs = NumBndryCoeffs();
                int nInteriorDofs = m_ncoeffs - nBoundaryDofs;


                Array<OneD, NekDouble> rhs(nInteriorDofs);
                Array<OneD, NekDouble> result(nInteriorDofs);

                GetInteriorMap(mapArray);

                for(i = 0; i < nInteriorDofs; i++)
                {
                    rhs[i] = tmp1[ mapArray[i] ];
                }

                Blas::Dgemv('N',nInteriorDofs,nInteriorDofs,1.0, &(matsys->GetPtr())[0],
                            nInteriorDofs,rhs.get(),1,0.0,result.get(),1);

                for(i = 0; i < nInteriorDofs; i++)
                {
                    outarray[ mapArray[i] ] = result[i];
                }
            }

        }

        /////////////////////////////
        // Inner Product Functions //
        /////////////////////////////

       /** \brief Calculate the inner product of inarray with respect to
        *  the basis B=base0*base1 and put into outarray
        *
        *  \f$
        *  \begin{array}{rcl}
        *  I_{pq} = (\phi_p \phi_q, u) & = & \sum_{i=0}^{nq_0}
        *  \sum_{j=0}^{nq_1}
        *  \phi_p(\xi_{0,i}) \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i}
        *  \xi_{1,j}) \\
        *  & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
        *  \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}
        *  \end{array}
        *  \f$
        *
        *  where
        *
        *  \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
        *
        *  which can be implemented as
        *
        *  \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j})
        *  \tilde{u}_{i,j} = {\bf B_1 U}  \f$
        *  \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} =
        *  {\bf B_0 F}  \f$
        */
        void StdQuadExp::v_IProductWRTBase(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &outarray)
        {
            if(m_base[0]->Collocation() && m_base[1]->Collocation())
            {
                MultiplyByQuadratureMetric(inarray,outarray);
            }
            else
            {
                StdQuadExp::v_IProductWRTBase_SumFac(inarray,outarray);
            }
        }

        void StdQuadExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble> &outarray,
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
                                             tmp,outarray,wsp,true,true);
            }
            else
            {
                Array<OneD,NekDouble> wsp(nquad1*order0);
                IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                             m_base[1]->GetBdata(),
                                             inarray,outarray,wsp,true,true);
            }
        }

        void StdQuadExp::v_IProductWRTBase_MatOp(
                             const Array<OneD, const NekDouble>& inarray,
                             Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void StdQuadExp::v_IProductWRTDerivBase(const int dir,
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> & outarray)
        {
            StdQuadExp::IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }

        void StdQuadExp::v_IProductWRTDerivBase_SumFac(const int dir,
                             const Array<OneD, const NekDouble>& inarray,
                             Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0((dir==0)||(dir==1),"input dir is out of range");

            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int     nqtot = nquad0*nquad1;
            int    order0 = m_base[0]->GetNumModes();

            Array<OneD,NekDouble> tmp(nqtot+nquad1*order0);
            Array<OneD,NekDouble> wsp(tmp+nqtot);

            // multiply by integration constants
            MultiplyByQuadratureMetric(inarray,tmp);

            if(dir) // dir == 1
            {
                IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                             m_base[1]->GetDbdata(),
                                             tmp,outarray,wsp,true,false);
            }
            else    // dir == 0
            {
                IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                             m_base[1]->GetBdata(),
                                             tmp,outarray,wsp,false,true);
            }
        }

        void StdQuadExp::v_IProductWRTDerivBase_MatOp(const int dir,
                             const Array<OneD, const NekDouble>& inarray,
                             Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0((dir==0)||(dir==1),"input dir is out of range");

            int nq = GetTotPoints();
            MatrixType mtype;

            if(dir) // dir == 1
            {
                mtype = eIProductWRTDerivBase1;
            }
            else    // dir == 0
            {
                mtype = eIProductWRTDerivBase0;
            }

            StdMatrixKey      iprodmatkey(mtype,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        // the arguments doCheckCollDir0 and doCheckCollDir1 allow you to specify whether
        // to check if the basis has collocation properties (i.e. for the classical spectral
        // element basis, In this case the 1D 'B' matrix is equal to the identity matrix
        // which can be exploited to speed up the calculations).
        // However, as this routine also allows to pass the matrix 'DB' (derivative of the basis),
        // the collocation property cannot always be used. Therefor follow this rule:
        // if base0 == m_base[0]->GetBdata() --> set doCheckCollDir0 == true;
        //    base1 == m_base[1]->GetBdata() --> set doCheckCollDir1 == true;
        //    base0 == m_base[0]->GetDbdata() --> set doCheckCollDir0 == false;
        //    base1 == m_base[1]->GetDbdata() --> set doCheckCollDir1 == false;
        void StdQuadExp::v_IProductWRTBase_SumFacKernel(
                             const Array<OneD, const NekDouble>& base0,
                             const Array<OneD, const NekDouble>& base1,
                             const Array<OneD, const NekDouble>& inarray,
                             Array<OneD, NekDouble> &outarray,
                             Array<OneD, NekDouble> &wsp,
                             bool doCheckCollDir0,
                             bool doCheckCollDir1)
            {
                int    nquad0 = m_base[0]->GetNumPoints();
                int    nquad1 = m_base[1]->GetNumPoints();
                int    nmodes0 = m_base[0]->GetNumModes();
                int    nmodes1 = m_base[1]->GetNumModes();

                bool colldir0 = doCheckCollDir0?(m_base[0]->Collocation()):false;
                bool colldir1 = doCheckCollDir1?(m_base[1]->Collocation()):false;

                if(colldir0 && colldir1)
                {
                    Vmath::Vcopy(m_ncoeffs,inarray.get(),1,outarray.get(),1);
                }
                else if(colldir0)
                {
                    Blas::Dgemm('N', 'N',nmodes0,nmodes1, nquad1,1.0, inarray.get(),
                                nmodes0, base1.get(), nquad1, 0.0,outarray.get(),nmodes0);
                }
                else if(colldir1)
                {
                    Blas::Dgemm('T','N',nmodes0,nquad1,nquad0,1.0,base0.get(),
                                nquad0,inarray.get(),nquad0,0.0,outarray.get(),nmodes0);
                }
                else
                {
                    ASSERTL1(wsp.size()>=nquad1*nmodes0,"Workspace size is not sufficient");

#if 1
                    Blas::Dgemm('T','N',nmodes0,nquad1,nquad0,1.0,base0.get(),
                                nquad0,inarray.get(),nquad0,0.0,wsp.get(),nmodes0);

#else
                    for(int i = 0; i < nmodes0; ++i)
                    {
                        for(int j = 0; j < nquad1; ++j)
                        {
                            wsp[j*nmodes0+i] = Blas::Ddot(nquad0,
                                                          base0.get()+i*nquad0,1,
                                                          inarray.get()+j*nquad0,1);
                        }
                    }
#endif
                    Blas::Dgemm('N', 'N',nmodes0,nmodes1, nquad1,1.0, wsp.get(),
                                nmodes0, base1.get(), nquad1, 0.0,outarray.get(),nmodes0);
                }
            }

        //////////////////////////
        // Evaluation functions //
        //////////////////////////


        void StdQuadExp::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble>& xi,
                                                 Array<OneD, NekDouble>& eta)
        {
            eta[0] = xi[0];
            eta[1] = xi[1];
        }

        /** \brief Fill outarray with mode \a mode of expansion
         *
         *  Note for quadrilateral expansions _base[0] (i.e. p)  modes run
         *  fastest
         */

        void StdQuadExp::v_FillMode(const int mode,
                            Array<OneD, NekDouble> &outarray)
        {
            int    i;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, const NekDouble> base0  = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> base1  = m_base[1]->GetBdata();
            int   btmp0 = m_base[0]->GetNumModes();
            int   mode0 = mode%btmp0;
            int   mode1 = mode/btmp0;


            ASSERTL2(mode1 == (int)floor((1.0*mode)/btmp0),
                     "Integer Truncation not Equiv to Floor");

            ASSERTL2(m_ncoeffs <= mode,
                     "calling argument mode is larger than total expansion order");

            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vcopy(nquad0,(NekDouble *)(base0.get() + mode0*nquad0),
                             1, &outarray[0]+i*nquad0,1);
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,(NekDouble *)(base1.get() + mode1*nquad1),1,
                            &outarray[0]+i,nquad0,&outarray[0]+i,nquad0);
            }
        }

        //////////////////////
        // Helper functions //
        //////////////////////

        int StdQuadExp::v_GetNverts() const
        {
            return 4;
        }

        int StdQuadExp::v_GetNedges() const
        {
            return 4;
        }

        int StdQuadExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 3),"edge id is out of range");

            if((i == 0)||(i == 2))
            {
                return  GetBasisNumModes(0);
            }
            else
            {
                return  GetBasisNumModes(1);
            }
        }

        int StdQuadExp::v_GetEdgeNumPoints(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 3),"edge id is out of range");

            if((i == 0)||(i == 2))
            {
                return  GetNumPoints(0);
            }
            else
            {
                return  GetNumPoints(1);
            }
        }

        LibUtilities::BasisType StdQuadExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 3),"edge id is out of range");

            if((i == 0)||(i == 2))
            {
                return  GetBasisType(0);
            }
            else
            {
                return  GetBasisType(1);
            }
        }

        int StdQuadExp::v_DetCartesianDirOfEdge(const int edge)
        {
            ASSERTL2((edge >= 0)&&(edge <= 3),"edge id is out of range");

            if((edge == 0)||(edge == 2))
            {
                return  0;
            }
            else
            {
                return  1;
            }
        }

        const LibUtilities::BasisKey StdQuadExp::v_DetEdgeBasisKey(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 3),"edge id is out of range");

            if((i == 0)||(i == 2))
            {
                return  GetBasis(0)->GetBasisKey();
            }
            else
            {
               return  GetBasis(1)->GetBasisKey();
            }

        }

        LibUtilities::ShapeType StdQuadExp::v_DetShapeType() const
        {
            return LibUtilities::eQuadrilateral;
        }


        int StdQuadExp::v_NumBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange||
                     GetBasisType(0) == LibUtilities::eGauss_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange||
                     GetBasisType(0) == LibUtilities::eGauss_Lagrange,
                      "BasisType is not a boundary interior form");

            return 4 + 2*(GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
        }

        int StdQuadExp::v_NumDGBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange ||
                     GetBasisType(0) == LibUtilities::eGauss_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange ||
                     GetBasisType(0) == LibUtilities::eGauss_Lagrange,
                     "BasisType is not a boundary interior form");

            return  2*GetBasisNumModes(0) + 2*GetBasisNumModes(1);
        }

        int StdQuadExp::v_CalcNumberOfCoefficients(
                           const std::vector<unsigned int> &nummodes,
                           int &modes_offset)
        {
            int nmodes = nummodes[modes_offset]*nummodes[modes_offset+1];
            modes_offset += 2;

            return nmodes;
        }

        bool StdQuadExp::v_IsBoundaryInteriorExpansion()
        {
            bool returnval = false;

            if((m_base[0]->GetBasisType() == LibUtilities::eModified_A)
               ||(m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange))
            {
                if((m_base[1]->GetBasisType() == LibUtilities::eModified_A)
                   ||(m_base[1]->GetBasisType() == LibUtilities::eGLL_Lagrange))
                {
                    returnval = true;
                }
            }

            return returnval;
        }

        void StdQuadExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
			    Array<OneD, NekDouble> &coords_1,
			    Array<OneD, NekDouble> &coords_2)
        {
            boost::ignore_unused(coords_2);
            Array<OneD, const NekDouble> z0 = m_base[0]->GetZ();
            Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();
            int nq0 = GetNumPoints(0);
            int nq1 = GetNumPoints(1);
            int i;

            for(i = 0; i < nq1; ++i)
            {
                Blas::Dcopy(nq0,z0.get(), 1,&coords_0[0] + i*nq0,1);
                Vmath::Fill(nq0,z1[i],&coords_1[0] + i*nq0,1);
            }
        }

        //////////////
        // Mappings //
        //////////////

        void StdQuadExp::v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            int i;
            int cnt=0;
            int nummodes0, nummodes1;
            int value1 = 0, value2 = 0;
            if(outarray.size()!=NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
            }

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            const LibUtilities::BasisType Btype0 = GetBasisType(0);
            const LibUtilities::BasisType Btype1 = GetBasisType(1);

            switch(Btype1)
            {
            case LibUtilities::eGLL_Lagrange:
            case LibUtilities::eGauss_Lagrange:
                value1 = nummodes0;
                break;
            case LibUtilities::eModified_A:
                value1 = 2*nummodes0;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < value1; i++)
            {
                outarray[i]=i;
            }
            cnt=value1;

            switch(Btype0)
            {
            case LibUtilities::eGLL_Lagrange:
            case LibUtilities::eGauss_Lagrange:
                value2 = value1+nummodes0-1;
                break;
            case LibUtilities::eModified_A:
                value2 = value1+1;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < nummodes1-2; i++)
            {
                outarray[cnt++]=value1+i*nummodes0;
                outarray[cnt++]=value2+i*nummodes0;
            }


            if(Btype1 == LibUtilities::eGLL_Lagrange || Btype1 == LibUtilities::eGauss_Lagrange )
            {
                for(i = nummodes0*(nummodes1-1);i < GetNcoeffs(); i++)
                {
                    outarray[cnt++] = i;
                }
            }
        }

        void StdQuadExp::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            int i,j;
            int cnt=0;
            int nummodes0, nummodes1;
            int startvalue = 0;
            if(outarray.size()!=GetNcoeffs()-NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(GetNcoeffs()-NumBndryCoeffs());
            }

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            const LibUtilities::BasisType Btype0 = GetBasisType(0);
            const LibUtilities::BasisType Btype1 = GetBasisType(1);

            switch(Btype1)
            {
            case LibUtilities::eGLL_Lagrange:
                startvalue = nummodes0;
                break;
            case LibUtilities::eModified_A:
                startvalue = 2*nummodes0;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            switch(Btype0)
            {
            case LibUtilities::eGLL_Lagrange:
                startvalue++;
                break;
            case LibUtilities::eModified_A:
                startvalue+=2;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < nummodes1-2; i++)
            {
                for(j = 0; j < nummodes0-2; j++)
                {
                    outarray[cnt++]=startvalue+j;
                }
                startvalue+=nummodes0;
            }
        }

        int StdQuadExp::v_GetVertexMap(int localVertexId, bool useCoeffPacking)
        {
            int localDOF = 0;

            if(useCoeffPacking == true)
            {
                switch(localVertexId)
                {
                case 0:
                    {
                        localDOF = 0;
                    }
                    break;
                case 1:
                    {
                        if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                        {
                            localDOF = m_base[0]->GetNumModes()-1;
                        }
                        else
                        {
                            localDOF = 1;
                        }
                    }
                    break;
                case 2:
                    {
                        if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                        {
                            localDOF = m_base[0]->GetNumModes() * (m_base[1]->GetNumModes()-1);
                        }
                        else
                        {
                            localDOF = m_base[0]->GetNumModes();
                        }
                    }
                    break;
                case 3:
                    {
                        if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                        {
                            localDOF = m_base[0]->GetNumModes()*m_base[1]->GetNumModes()-1;
                        }
                        else
                        {
                            localDOF = m_base[0]->GetNumModes()+1;
                        }
                    }
                break;
                default:
                    ASSERTL0(false,"eid must be between 0 and 3");
                    break;
                }

            }
            else
            {
                switch(localVertexId)
                {
                case 0:
                    {
                        localDOF = 0;
                    }
                    break;
                case 1:
                    {
                        if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                        {
                            localDOF = m_base[0]->GetNumModes()-1;
                        }
                        else
                        {
                            localDOF = 1;
                        }
                    }
                    break;
                case 2:
                    {
                        if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                        {
                            localDOF = m_base[0]->GetNumModes()*m_base[1]->GetNumModes()-1;
                        }
                        else
                        {
                            localDOF = m_base[0]->GetNumModes()+1;
                        }
                    }
                break;
                case 3:
                    {
                        if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                        {
                            localDOF = m_base[0]->GetNumModes() * (m_base[1]->GetNumModes()-1);
                        }
                        else
                        {
                            localDOF = m_base[0]->GetNumModes();
                        }
                    }
                    break;
                default:
                    ASSERTL0(false,"eid must be between 0 and 3");
                    break;
                }
            }
            return localDOF;
        }

        void StdQuadExp::v_GetEdgeInteriorMap(const int eid,
                             const Orientation edgeOrient,
                             Array<OneD, unsigned int> &maparray,
                             Array<OneD, int> &signarray)
        {
            int i;
            const int nummodes0 = m_base[0]->GetNumModes();
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nEdgeIntCoeffs = GetEdgeNcoeffs(eid)-2;
            const LibUtilities::BasisType bType = GetEdgeBasisType(eid);

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

            if(bType == LibUtilities::eModified_A)
            {
                switch(eid)
                {
                case 0:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = i+2;
                        }

                        if(edgeOrient==eBackwards)
                        {
                            for(i = 1; i < nEdgeIntCoeffs; i+=2)
                            {
                                signarray[i] = -1;
                            }
                        }
                    }
                    break;
                case 1:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = (i+2)*nummodes0 + 1;
                        }

                        if(edgeOrient==eBackwards)
                        {
                            for(i = 1; i < nEdgeIntCoeffs; i+=2)
                            {
                                signarray[i] = -1;
                            }
                        }
                    }
                    break;
                case 2:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = nummodes0+i+2;
                        }

                        if(edgeOrient==eBackwards)
                        {
                            for(i = 1; i < nEdgeIntCoeffs; i+=2)
                            {
                                signarray[i] = -1;
                            }
                        }
                    }
                    break;
                case 3:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = (i+2)*nummodes0;
                        }

                        if(edgeOrient==eBackwards)
                        {
                            for(i = 1; i < nEdgeIntCoeffs; i+=2)
                            {
                                signarray[i] = -1;
                            }
                        }
                    }
                    break;
                default:
                    ASSERTL0(false,"eid must be between 0 and 3");
                    break;
                }
            }
            else if(bType == LibUtilities::eGLL_Lagrange)
            {
                switch(eid)
                {
                case 0:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = i+1;
                        }
                    }
                    break;
                case 1:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = (i+2)*nummodes0 - 1;
                        }
                    }
                    break;
                case 2:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = nummodes0*(nummodes1-1) + i + 1;
                        }
                    }
                    break;
                case 3:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = nummodes0 * (i+1);
                        }
                    }
                    break;
                default:
                    ASSERTL0(false,"eid must be between 0 and 3");
                    break;
                }
                if(edgeOrient == eBackwards)
                {
                    reverse( maparray.get() , maparray.get()+nEdgeIntCoeffs );
                }
            }
            else
            {
                ASSERTL0(false,"Mapping not defined for this type of basis");
            }

        }

        void StdQuadExp::v_GetEdgeToElementMap(
             const int                  eid,
             const Orientation          edgeOrient,
             Array<OneD, unsigned int>& maparray,
             Array<OneD, int>&          signarray,
             int                        P)
        {
            int i;
            int numModes=0;
            int order0 = m_base[0]->GetNumModes();
            int order1 = m_base[1]->GetNumModes();

            switch (eid)
            {
            case 0:
            case 2:
                numModes = order0;
                break;
            case 1:
            case 3:
                numModes = order1;
                break;
            default:
                ASSERTL0(false,"eid must be between 0 and 3");
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
            const LibUtilities::BasisType bType = GetEdgeBasisType(eid);


            if (maparray.size() != P)
            {
                maparray = Array<OneD, unsigned int>(P);
            }

            if(signarray.size() != P)
            {
                signarray = Array<OneD, int>(P, 1);
            }
            else
            {
                fill(signarray.get(), signarray.get()+P, 1);
            }

            if (bType == LibUtilities::eModified_A)
            {
                switch (eid)
                {
                case 0:
                    {
                        for (i = 0; i < P; i++)
                        {
                            maparray[i] = i;
                        }

                        if (edgeOrient == eBackwards)
                        {
                            swap(maparray[0], maparray[1]);

                            for(i = 3; i < P; i+=2)
                            {
                                signarray[i] = -1;
                            }
                        }
                    }
                    break;
                case 1:
                    {
                        for (i = 0; i < P; i++)
                        {
                            maparray[i] = i*order0 + 1;
                        }

                        if (edgeOrient == eBackwards)
                        {
                            swap(maparray[0], maparray[1]);

                            for(i = 3; i < P; i+=2)
                            {
                                signarray[i] = -1;
                            }
                        }
                    }
                    break;
                case 2:
                    {
                        for (i = 0; i < P; i++)
                        {
                            maparray[i] = order0+i;
                        }

                        if (edgeOrient == eBackwards)
                        {
                            swap(maparray[0], maparray[1]);

                            for (i = 3; i < P; i+=2)
                            {
                                signarray[i] = -1;
                            }
                        }
                    }
                    break;
                case 3:
                    {
                        for (i = 0; i < P; i++)
                        {
                            maparray[i] = i*order0;
                        }

                        if (edgeOrient == eBackwards)
                        {
                            swap(maparray[0], maparray[1]);

                            for (i = 3; i < P; i+=2)
                            {
                                signarray[i] = -1;
                            }
                        }
                    }
                    break;
                default:
                    ASSERTL0(false, "eid must be between 0 and 3");
                    break;
                }
            }
            else if(bType == LibUtilities::eGLL_Lagrange ||
                    bType == LibUtilities::eGauss_Lagrange)
            {
                switch (eid)
                {
                case 0:
                    {
                        for (i = 0; i < P; i++)
                        {
                            maparray[i] = i;
                        }
                    }
                    break;
                case 1:
                    {
                        for (i = 0; i < P; i++)
                        {
                            maparray[i] = (i+1)*order0 - 1;
                        }
                    }
                    break;
                case 2:
                    {
                        for (i = 0; i < P; i++)
                        {
                            maparray[i] = order0*(order1-1) + i;
                        }
                    }
                    break;
                case 3:
                    {
                        for (i = 0; i < P; i++)
                        {
                            maparray[i] = order0*i;
                        }
                    }
                    break;
                default:
                    ASSERTL0(false, "eid must be between 0 and 3");
                    break;
                }
                if (edgeOrient == eBackwards)
                {
                    reverse(maparray.get(), maparray.get()+P);
                }
            }
            else
            {
                ASSERTL0(false, "Mapping not defined for this type of basis");
            }

            if (checkForZeroedModes)
            {
                if (bType == LibUtilities::eModified_A)
                {
                    // Zero signmap and set maparray to zero if
                    // elemental modes are not as large as face modes
                    for (int j = numModes; j < P; j++)
                    {
                        signarray[j] = 0.0;
                        maparray[j]  = maparray[0];
                    }
                }
                else
                {
                    ASSERTL0(false, "Different trace space edge dimension "
                                    "and element edge dimension not possible "
                                    "for GLL-Lagrange bases");
                }
            }
        }



        ///////////////////////
        // Wrapper Functions //
        ///////////////////////

        DNekMatSharedPtr StdQuadExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            int      i;
            int      order0    = GetBasisNumModes(0);
            int      order1    = GetBasisNumModes(1);
            MatrixType mtype   = mkey.GetMatrixType();

            DNekMatSharedPtr Mat;

            switch(mtype)
            {
                case ePhysInterpToEquiSpaced:
                {
                    int nq0 = m_base[0]->GetNumPoints();
                    int nq1 = m_base[1]->GetNumPoints();
                    int nq;

                    // take definition from key
                    if(mkey.ConstFactorExists(eFactorConst))
                    {
                        nq = (int) mkey.GetConstFactor(eFactorConst);
                    }
                    else
                    {
                        nq = max(nq0,nq1);
                    }

                    int neq = LibUtilities::StdQuadData::
                                                getNumberOfCoefficients(nq, nq);
                    Array<OneD, Array<OneD, NekDouble> > coords(neq);
                    Array<OneD, NekDouble>               coll  (2);
                    Array<OneD, DNekMatSharedPtr>        I     (2);
                    Array<OneD, NekDouble>               tmp   (nq0);

                    Mat = MemoryManager<DNekMat>::
                                        AllocateSharedPtr(neq, nq0 * nq1);
                    int cnt = 0;

                    for(int i = 0; i < nq; ++i)
                    {
                        for(int j = 0; j < nq; ++j,++cnt)
                        {
                            coords[cnt] = Array<OneD, NekDouble>(2);
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
                            Vmath::Smul(nq0,fac,I[0]->GetPtr(),1,tmp,1);

                            Vmath::Vcopy(nq0, &tmp[0], 1,
                                         Mat->GetRawPtr()+j*nq0*neq+i,neq);
                        }

                    }
                    break;
                }
                case eMass:
                {
                    Mat = StdExpansion::CreateGeneralMatrix(mkey);
                    // For Fourier basis set the imaginary component of mean mode
                    // to have a unit diagonal component in mass matrix
                    if(m_base[0]->GetBasisType() == LibUtilities::eFourier)
                    {
                        for(i = 0; i < order1; ++i)
                        {
                            (*Mat)(order0*i+1,i*order0+1) = 1.0;
                        }
                    }

                    if(m_base[1]->GetBasisType() == LibUtilities::eFourier)
                    {
                        for(i = 0; i < order0; ++i)
                        {
                            (*Mat)(order0+i ,order0+i) = 1.0;
                        }
                    }
                    break;
                }
                case eFwdTrans:
                {
                    Mat = MemoryManager<DNekMat>::AllocateSharedPtr(
                            m_ncoeffs,m_ncoeffs);
                    StdMatrixKey iprodkey(eIProductWRTBase,DetShapeType(),*this);
                    DNekMat &Iprod = *GetStdMatrix(iprodkey);
                    StdMatrixKey imasskey(eInvMass,DetShapeType(),*this);
                    DNekMat &Imass = *GetStdMatrix(imasskey);

                    (*Mat) = Imass*Iprod;
                    break;
                }
                case eGaussDG:
                {
                    ConstFactorMap factors = mkey.GetConstFactors();

                    int edge    = (int)factors[StdRegions::eFactorGaussEdge];
                    int dir     = (edge + 1) % 2;
                    int nCoeffs = m_base[dir]->GetNumModes();

                    const LibUtilities::PointsKey BS_p(
                        nCoeffs, LibUtilities::eGaussGaussLegendre);
                    const LibUtilities::BasisKey  BS_k(
                        LibUtilities::eGauss_Lagrange, nCoeffs, BS_p);

                    Array<OneD, NekDouble> coords(1, 0.0);
                    coords[0] = (edge == 0 || edge == 3) ? -1.0 : 1.0;

                    LibUtilities::BasisSharedPtr basis =
                        LibUtilities::BasisManager()[BS_k];
                    DNekMatSharedPtr             m_Ix  = basis->GetI(coords);

                    Mat = MemoryManager<DNekMat>::AllocateSharedPtr(
                        1.0, nCoeffs);
                    Vmath::Vcopy(nCoeffs, m_Ix->GetPtr(), 1, Mat->GetPtr(), 1);
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

        DNekMatSharedPtr StdQuadExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return GenMatrix(mkey);
        }

        ///////////////////////////////////
        // Operator evaluation functions //
        ///////////////////////////////////

        void StdQuadExp::v_GeneralMatrixOp_MatOp(
                             const Array<OneD, const NekDouble> &inarray,
                             Array<OneD,NekDouble> &outarray,
                             const StdMatrixKey &mkey)
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


        void StdQuadExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                              const StdMatrixKey &mkey)
        {
            // Generate an orthonogal expansion
            int qa = m_base[0]->GetNumPoints();
            int qb = m_base[1]->GetNumPoints();
            int nmodes_a = m_base[0]->GetNumModes();
            int nmodes_b = m_base[1]->GetNumModes();
            int nmodes = min(nmodes_a,nmodes_b);
            // Declare orthogonal basis.
            LibUtilities::PointsKey pa(qa,m_base[0]->GetPointsType());
            LibUtilities::PointsKey pb(qb,m_base[1]->GetPointsType());

            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A,nmodes_a,pa);
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A,nmodes_b,pb);
            StdQuadExp OrthoExp(Ba,Bb);

            //SVV parameters loaded from the .xml case file
            Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());

            // project onto modal  space.
            OrthoExp.FwdTrans(array,orthocoeffs);

            if(mkey.ConstFactorExists(eFactorSVVPowerKerDiffCoeff)) // Rodrigo's power kernel
            {
                NekDouble cutoff = mkey.GetConstFactor(eFactorSVVCutoffRatio);
                NekDouble  SvvDiffCoeff  =
                    mkey.GetConstFactor(eFactorSVVPowerKerDiffCoeff)*
                    mkey.GetConstFactor(eFactorSVVDiffCoeff);

                for(int j = 0; j < nmodes_a; ++j)
                {
                    for(int k = 0; k < nmodes_b; ++k)
                    {
                        // linear space but makes high modes very negative
                        NekDouble fac = std::max(
                                    pow((1.0*j)/(nmodes_a-1),cutoff*nmodes_a),
                                    pow((1.0*k)/(nmodes_b-1),cutoff*nmodes_b));

                        orthocoeffs[j*nmodes_b+k] *= SvvDiffCoeff * fac;
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

                for(int j = 0; j < nmodes_a; ++j)
                {
                    for(int k = 0; k < nmodes_b; ++k)
                    {
                        int maxjk = max(j,k);
                        maxjk = min(maxjk,kSVVDGFiltermodesmax-1);

                        orthocoeffs[j*nmodes_b+k] *= SvvDiffCoeff *
                            kSVVDGFilter[max_ab][maxjk];
                    }
                }
            }
            else
            {
                NekDouble  SvvDiffCoeff = mkey.GetConstFactor(eFactorSVVDiffCoeff);
                //Exponential Kernel implementation
                int cutoff = (int) (mkey.GetConstFactor(eFactorSVVCutoffRatio)*
                                                        min(nmodes_a,nmodes_b));

                //counters for scanning through orthocoeffs array
                int cnt = 0;

                //------"New" Version August 22nd '13--------------------
                for(int j = 0; j < nmodes_a; ++j)
                {
                    for(int k = 0; k < nmodes_b; ++k)
                    {
                        if(j + k >= cutoff)//to filter out only the "high-modes"
                        {
                             orthocoeffs[j*nmodes_b+k] *=
                                 (SvvDiffCoeff*exp(-(j+k-nmodes)*(j+k-nmodes)/
                                                    ((NekDouble)((j+k-cutoff+1)*
                                                         (j+k-cutoff+1)))));
                         }
                        else
                        {
                             orthocoeffs[j*nmodes_b+k] *= 0.0;
                        }
                        cnt++;
                    }
                }
            }

            // backward transform to physical space
            OrthoExp.BwdTrans(orthocoeffs,array);
        }

        void StdQuadExp::v_ExponentialFilter(
                                          Array<OneD, NekDouble> &array,
                                    const NekDouble        alpha,
                                    const NekDouble        exponent,
                                    const NekDouble        cutoff)
        {
            // Generate an orthogonal expansion
            int qa      = m_base[0]->GetNumPoints();
            int qb      = m_base[1]->GetNumPoints();
            int nmodesA = m_base[0]->GetNumModes();
            int nmodesB = m_base[1]->GetNumModes();
            int P  = nmodesA - 1;
            int Q  = nmodesB - 1;

            // Declare orthogonal basis.
            LibUtilities::PointsKey pa(qa,m_base[0]->GetPointsType());
            LibUtilities::PointsKey pb(qb,m_base[1]->GetPointsType());

            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, nmodesA, pa);
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, nmodesB, pb);
            StdQuadExp OrthoExp(Ba,Bb);

            // Cutoff
            int Pcut = cutoff*P;
            int Qcut = cutoff*Q;

            // Project onto orthogonal space.
            Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
            OrthoExp.FwdTrans(array,orthocoeffs);

            //
            NekDouble fac, fac1, fac2;
            for(int i = 0; i < nmodesA; ++i)
            {
                for(int j = 0; j < nmodesB; ++j)
                {
                    //to filter out only the "high-modes"
                    if(i > Pcut || j > Qcut)
                    {
                        fac1 = (NekDouble) (i - Pcut)/( (NekDouble)(P - Pcut) );
                        fac2 = (NekDouble) (j - Qcut)/( (NekDouble)(Q - Qcut) );
                        fac  = max(fac1, fac2);
                        fac  = pow(fac, exponent);
                        orthocoeffs[i*nmodesB+j] *= exp(-alpha*fac);
                    }
                }
            }

            // backward transform to physical space
            OrthoExp.BwdTrans(orthocoeffs,array);
        }

        void StdQuadExp::v_ReduceOrderCoeffs(
            int                                 numMin,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            int n_coeffs = inarray.size();


            Array<OneD, NekDouble> coeff(n_coeffs);
            Array<OneD, NekDouble> coeff_tmp(n_coeffs,0.0);
            Array<OneD, NekDouble> tmp;
            Array<OneD, NekDouble> tmp2;

            int       nmodes0 = m_base[0]->GetNumModes();
            int       nmodes1 = m_base[1]->GetNumModes();
            int       numMax  = nmodes0;

            Vmath::Vcopy(n_coeffs,inarray,1,coeff_tmp,1);

            const LibUtilities::PointsKey Pkey0(
                nmodes0, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey Pkey1(
                nmodes1, LibUtilities::eGaussLobattoLegendre);

            LibUtilities::BasisKey b0(m_base[0]->GetBasisType(),nmodes0,Pkey0);
            LibUtilities::BasisKey b1(m_base[1]->GetBasisType(),nmodes1,Pkey1);

            LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A,nmodes0,Pkey0);
            LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_A,nmodes1,Pkey1);

            LibUtilities::InterpCoeff2D(
                b0, b1, coeff_tmp, bortho0, bortho1, coeff);

            Vmath::Zero(n_coeffs,coeff_tmp,1);

            int cnt = 0;
            for (int i = 0; i < numMin+1; ++i)
            {
                Vmath::Vcopy(numMin,
                             tmp  = coeff+cnt,1,
                             tmp2 = coeff_tmp+cnt,1);

                cnt = i*numMax;
            }

            LibUtilities::InterpCoeff2D(
                bortho0, bortho1, coeff_tmp, b0, b1, outarray);
        }


        void StdQuadExp::v_MassMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void StdQuadExp::v_LaplacianMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdQuadExp::v_LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void StdQuadExp::v_LaplacianMatrixOp(const int k1, const int k2,
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdExpansion::LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,mkey);
        }

        void StdQuadExp::v_WeakDerivMatrixOp(const int i,
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
        }

        void StdQuadExp::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD,NekDouble> &outarray,
                           const StdMatrixKey &mkey)
        {
            StdQuadExp::v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
        }

        //up to here
        void StdQuadExp::v_MultiplyByStdQuadratureMetric(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            int i;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();

            // multiply by integration constants
            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0, inarray.get()+i*nquad0,1,
                            w0.get(),1,outarray.get()+i*nquad0,1);
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,outarray.get()+i,nquad0,w1.get(),1,
                            outarray.get()+i,nquad0);
            }
        }

        void StdQuadExp::v_GetSimplexEquiSpacedConnectivity(
            Array<OneD, int> &conn,
            bool              standard)
        {
            boost::ignore_unused(standard);

            int np1 = m_base[0]->GetNumPoints();
            int np2 = m_base[1]->GetNumPoints();
            int np = max(np1,np2);

            conn = Array<OneD, int>(6*(np-1)*(np-1));

            int row   = 0;
            int rowp1 = 0;
            int cnt = 0;
            for(int i = 0; i < np-1; ++i)
            {
                rowp1 += np;
                for(int j = 0; j < np-1; ++j)
                {
                    conn[cnt++] = row   +j;
                    conn[cnt++] = row   +j+1;
                    conn[cnt++] = rowp1 +j;

                    conn[cnt++] = rowp1 +j+1;
                    conn[cnt++] = rowp1 +j;
                    conn[cnt++] = row   +j+1;
                }
                row += np;
            }
        }

    } //end namespace
}//end namespace
