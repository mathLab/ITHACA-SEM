///////////////////////////////////////////////////////////////////////////////
//
// File StdHexExp.cpp
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
// Description: Heaxhedral methods
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdHexExp.h>

#ifdef max
#undef max
#endif

using namespace std;

namespace Nektar
{
    namespace StdRegions
    {

        StdHexExp::StdHexExp()
        {
        }


        StdHexExp::StdHexExp(const LibUtilities::BasisKey &Ba,
                        const LibUtilities::BasisKey &Bb,
                        const LibUtilities::BasisKey &Bc):
            StdExpansion(Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes(), 3,
                                                   Ba, Bb, Bc),
            StdExpansion3D(Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes(),
                           Ba, Bb, Bc)
        {
        }


        StdHexExp::StdHexExp(const StdHexExp &T):
            StdExpansion(T),
            StdExpansion3D(T)
        {
        }


        StdHexExp::~StdHexExp()
        {
        }

        bool StdHexExp::v_IsBoundaryInteriorExpansion()
        {
            return
                (m_base[0]->GetBasisType() == LibUtilities::eModified_A &&
                 m_base[1]->GetBasisType() == LibUtilities::eModified_A &&
                 m_base[2]->GetBasisType() == LibUtilities::eModified_A) ||
                (m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange &&
                 m_base[1]->GetBasisType() == LibUtilities::eGLL_Lagrange &&
                 m_base[1]->GetBasisType() == LibUtilities::eGLL_Lagrange);
        }


        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////
        /**
         * For Hexahedral region can use the PhysTensorDeriv function defined
         * under StdExpansion. Following tenserproduct:
         */
        void StdHexExp::v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                  Array<OneD, NekDouble> &out_d0,
                                  Array<OneD, NekDouble> &out_d1,
                                  Array<OneD, NekDouble> &out_d2)
        {
            StdExpansion3D::PhysTensorDeriv(inarray, out_d0, out_d1, out_d2);
        }


        /**
         * @param   dir         Direction in which to compute derivative.
         *                      Valid values are 0, 1, 2.
         * @param   inarray     Input array.
         * @param   outarray    Output array.
         */
        void StdHexExp::v_PhysDeriv(const int dir,
                               const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD,       NekDouble>& outarray)
        {
            switch(dir)
            {
            case 0:
                {
                    PhysDeriv(inarray, outarray, NullNekDouble1DArray,
                              NullNekDouble1DArray);
                }
                break;
            case 1:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray, outarray,
                              NullNekDouble1DArray);
                }
                break;
            case 2:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray,
                              NullNekDouble1DArray, outarray);
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
        }

        void StdHexExp::v_StdPhysDeriv(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &out_d0,
                  Array<OneD,       NekDouble> &out_d1,
                  Array<OneD,       NekDouble> &out_d2)
        {
            StdHexExp::v_PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }


        void StdHexExp::v_StdPhysDeriv(const int dir,
                               const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD,       NekDouble>& outarray)
        {
            StdHexExp::v_PhysDeriv(dir, inarray, outarray);
        }

        /**
         * Backward transformation is three dimensional tensorial expansion
         * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k})
         *  = \sum_{p=0}^{Q_x} \psi_p^a (\xi_{1i})
         *  \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
         *    \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{r}^a (\xi_{3k})
         *    \rbrace}
         *  \rbrace}. \f$
         * And sumfactorizing step of the form is as:\\
         * \f$ f_{r} (\xi_{3k})
         * = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{r}^a (\xi_{3k}),\\
         * g_{p} (\xi_{2j}, \xi_{3k})
         * = \sum_{r=0}^{Q_y} \psi_{p}^a (\xi_{2j}) f_{r} (\xi_{3k}),\\
         * u(\xi_{1i}, \xi_{2j}, \xi_{3k})
         * = \sum_{p=0}^{Q_x} \psi_{p}^a (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).
         * \f$
         *
         * @param   inarray     ?
         * @param   outarray    ?
         */
        void StdHexExp::v_BwdTrans(
                                const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1( (m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)  ||
                      (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                      "Basis[1] is not a general tensor type");

            ASSERTL1( (m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
                      (m_base[2]->GetBasisType() != LibUtilities::eModified_C),
                      "Basis[2] is not a general tensor type");

            if(m_base[0]->Collocation() && m_base[1]->Collocation()
                    && m_base[2]->Collocation())
            {
                Vmath::Vcopy(m_base[0]->GetNumPoints()
                                * m_base[1]->GetNumPoints()
                                * m_base[2]->GetNumPoints(),
                             inarray, 1, outarray, 1);
            }
            else
            {
                StdHexExp::BwdTrans_SumFac(inarray,outarray);
            }
        }


        /**
         *
         */
        void StdHexExp::v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> wsp(m_base[0]->GetNumPoints()*
                                       m_base[2]->GetNumModes()*
                                       (m_base[1]->GetNumModes() + m_base[1]->GetNumPoints())); // FIX THIS

            BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                    m_base[1]->GetBdata(),
                                    m_base[2]->GetBdata(),
                                    inarray,outarray,wsp,true,true,true);
        }


        /**
         * @param   base0       x-dirn basis matrix
         * @param   base1       y-dirn basis matrix
         * @param   base2       z-dirn basis matrix
         * @param   inarray     Input vector of modes.
         * @param   outarray    Output vector of physical space data.
         * @param   wsp         Workspace of size Q_x*P_z*(P_y+Q_y)
         * @param   doCheckCollDir0     Check for collocation of basis.
         * @param   doCheckCollDir1     Check for collocation of basis.
         * @param   doCheckCollDir2     Check for collocation of basis.
         * @todo    Account for some directions being collocated. See
         *          StdQuadExp as an example.
         */
        void StdHexExp::v_BwdTrans_SumFacKernel(
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
            int  nquad0  = m_base[0]->GetNumPoints();
            int  nquad1  = m_base[1]->GetNumPoints();
            int  nquad2  = m_base[2]->GetNumPoints();
            int  nmodes0 = m_base[0]->GetNumModes();
            int  nmodes1 = m_base[1]->GetNumModes();
            int  nmodes2 = m_base[2]->GetNumModes();

            // Check if using collocation, if requested.
            bool colldir0 = doCheckCollDir0?(m_base[0]->Collocation()):false;
            bool colldir1 = doCheckCollDir1?(m_base[1]->Collocation()):false;
            bool colldir2 = doCheckCollDir2?(m_base[2]->Collocation()):false;

            // If collocation in all directions, Physical values at quadrature
            // points is just a copy of the modes.
            if(colldir0 && colldir1 && colldir2)
            {
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,outarray.get(),1);
            }
            else
            {
                // Check sufficiently large workspace.
                ASSERTL1(wsp.size()>=nquad0*nmodes2*(nmodes1+nquad1),
                         "Workspace size is not sufficient");

                // Assign second half of workspace for 2nd DGEMM operation.
                Array<OneD, NekDouble> wsp2 = wsp + nquad0*nmodes1*nmodes2;

                // BwdTrans in each direction using DGEMM
                Blas::Dgemm('T','T', nmodes1*nmodes2, nquad0, nmodes0,
                            1.0, &inarray[0],   nmodes0,
                                 base0.get(),   nquad0,
                            0.0, &wsp[0],       nmodes1*nmodes2);
                Blas::Dgemm('T','T', nquad0*nmodes2,  nquad1, nmodes1,
                            1.0, &wsp[0],       nmodes1,
                                 base1.get(),   nquad1,
                            0.0, &wsp2[0],      nquad0*nmodes2);
                Blas::Dgemm('T','T', nquad0*nquad1,   nquad2, nmodes2,
                            1.0, &wsp2[0],      nmodes2,
                                 base2.get(),   nquad2,
                            0.0, &outarray[0],  nquad0*nquad1);
            }
        }


        /**
         * Solves the system
         * \f$ \mathbf{B^{\top}WB\hat{u}}=\mathbf{B^{\top}Wu^{\delta}} \f$
         *
         * @param   inarray     array of physical quadrature points to be
         *                      transformed, \f$ \mathbf{u^{\delta}} \f$.
         * @param   outarray    array of expansion coefficients,
         *                      \f$ \mathbf{\hat{u}} \f$.
         */
        void StdHexExp::v_FwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray)
        {
            // If using collocation expansion, coefficients match physical
            // data points so just do a direct copy.
            if( (m_base[0]->Collocation())
                    &&(m_base[1]->Collocation())
                    &&(m_base[2]->Collocation()) )
            {
                Vmath::Vcopy(GetNcoeffs(), &inarray[0], 1, &outarray[0], 1);
            }
            else
            {
                // Compute B^TWu
                IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
                DNekMatSharedPtr matsys = GetStdMatrix(masskey);

                // copy inarray in case inarray == outarray
                DNekVec in (m_ncoeffs,outarray);
                DNekVec out(m_ncoeffs,outarray,eWrapper);

                // Solve for coefficients.
                out = (*matsys)*in;

            }
        }

        /**
         * \f$
         * \begin{array}{rcl}
         * I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
         * \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
         * \psi_{p}^{a}(\xi_{1i}) \psi_{q}^{a}(\xi_{2j}) \psi_{r}^{a}(\xi_{3k})
         * w_i w_j w_k u(\xi_{1,i} \xi_{2,j} \xi_{3,k})
         *
         * J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\xi_{1,i})
         *                   \sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2,j})
         *                   \sum_{k=0}^{nq_2} \psi_{r}^a
         *                   u(\xi_{1i},\xi_{2j},\xi_{3k}) J_{i,j,k}
         * \end{array} \f$ \n
         * where
         * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3)
         *  = \psi_p^a( \xi_1) \psi_{q}^a(\xi_2) \psi_{r}^a(\xi_3) \f$ \n
         * which can be implemented as \n
         * \f$f_{r} (\xi_{3k})
         *  = \sum_{k=0}^{nq_3} \psi_{r}^a u(\xi_{1i},\xi_{2j}, \xi_{3k})
         * J_{i,j,k} = {\bf B_3 U}   \f$ \n
         * \f$ g_{q} (\xi_{3k})
         *  = \sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2j}) f_{r}(\xi_{3k})
         *  = {\bf B_2 F}  \f$ \n
         * \f$ (\phi_{pqr}, u)_{\delta}
         *  = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k})  g_{q} (\xi_{3k})
         *  = {\bf B_1 G} \f$
         *
         * @param   inarray     ?
         * @param   outarray    ?
         */
        void StdHexExp::v_IProductWRTBase(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray)
        {
            if(m_base[0]->Collocation() &&
               m_base[1]->Collocation() &&
               m_base[2]->Collocation())
            {
                MultiplyByQuadratureMetric(inarray,outarray);
            }
            else
            {
                StdHexExp::v_IProductWRTBase_SumFac(inarray,outarray);
            }
        }

        /**
         * Implementation of the local matrix inner product operation.
         */
        void StdHexExp::v_IProductWRTBase_MatOp(const Array<OneD, const NekDouble>& inarray,
                                               Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        /**
         * Implementation of the sum-factorization inner product operation.
         */
        void StdHexExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble> &outarray,
            bool                                multiplybyweights)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
            int    order1 = m_base[1]->GetNumModes();

            Array<OneD, NekDouble> wsp(nquad0*nquad1*(nquad2+order0) +
                                       order0*order1*nquad2);

            if(multiplybyweights)
            {
                Array<OneD, NekDouble> tmp(inarray.size());
                MultiplyByQuadratureMetric(inarray,tmp);

                StdHexExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                           m_base[1]->GetBdata(),
                                           m_base[2]->GetBdata(),
                                           tmp,outarray,wsp,true,true,true);
            }
            else
            {
                StdHexExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                           m_base[1]->GetBdata(),
                                           m_base[2]->GetBdata(),
                                           inarray,outarray,wsp,true,true,true);
            }
        }


        /**
         * Implementation of the sum-factorisation inner product operation.
         * @todo    Implement cases where only some directions are collocated.
         */
        void StdHexExp::v_IProductWRTBase_SumFacKernel(const Array<OneD, const NekDouble>& base0,
                                                     const Array<OneD, const NekDouble>& base1,
                                                     const Array<OneD, const NekDouble>& base2,
                                                     const Array<OneD, const NekDouble>& inarray,
                                                     Array<OneD, NekDouble> &outarray,
                                                     Array<OneD, NekDouble> &wsp,
                                                     bool doCheckCollDir0,
                                                     bool doCheckCollDir1,
                                                     bool doCheckCollDir2)
        {
            int  nquad0  = m_base[0]->GetNumPoints();
            int  nquad1  = m_base[1]->GetNumPoints();
            int  nquad2  = m_base[2]->GetNumPoints();
            int  nmodes0 = m_base[0]->GetNumModes();
            int  nmodes1 = m_base[1]->GetNumModes();
            int  nmodes2 = m_base[2]->GetNumModes();

            bool colldir0 = doCheckCollDir0?(m_base[0]->Collocation()):false;
            bool colldir1 = doCheckCollDir1?(m_base[1]->Collocation()):false;
            bool colldir2 = doCheckCollDir2?(m_base[2]->Collocation()):false;

            if(colldir0 && colldir1 && colldir2)
            {
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,outarray.get(),1);
            }
            else
            {
                ASSERTL1(wsp.size() >= nmodes0*nquad2*(nquad1+nmodes1),
                         "Insufficient workspace size");

                Array<OneD, NekDouble> tmp0 = wsp;
                Array<OneD, NekDouble> tmp1 = wsp + nmodes0*nquad1*nquad2;


               if(colldir0)
               {
                    // reshuffle data for next operation.
                    for(int n = 0; n < nmodes0; ++n)
                    {
                        Vmath::Vcopy(nquad1*nquad2,inarray.get()+n,nquad0,
                                     tmp0.get()+nquad1*nquad2*n,1);
                    }
                }
                else
                {
                    Blas::Dgemm('T', 'N', nquad1*nquad2, nmodes0, nquad0,
                                1.0, inarray.get(),  nquad0,
                                 base0.get(),    nquad0,
                                0.0, tmp0.get(),     nquad1*nquad2);
                }

                if(colldir1)
                {
                    // reshuffle data for next operation.
                    for(int n = 0; n < nmodes1; ++n)
                    {
                        Vmath::Vcopy(nquad2*nmodes0,tmp0.get()+n,nquad1,
                                     tmp1.get()+nquad2*nmodes0*n,1);
                    }
                }
                else
                {
                    Blas::Dgemm('T', 'N', nquad2*nmodes0, nmodes1, nquad1,
                                1.0, tmp0.get(),     nquad1,
                                base1.get(),    nquad1,
                                0.0, tmp1.get(),     nquad2*nmodes0);
                }

                if(colldir2)
                {
                    // reshuffle data for next operation.
                    for(int n = 0; n < nmodes2; ++n)
                    {
                        Vmath::Vcopy(nmodes0*nmodes1,tmp1.get()+n,nquad2,
                                     outarray.get()+nmodes0*nmodes1*n,1);
                    }
                }
                else
                {
                    Blas::Dgemm('T', 'N', nmodes0*nmodes1, nmodes2, nquad2,
                                1.0, tmp1.get(),     nquad2,
                                base2.get(),    nquad2,
                                0.0, outarray.get(), nmodes0*nmodes1);
                }
           }
        }

        void StdHexExp::v_IProductWRTDerivBase(const int dir,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> & outarray)
        {
            StdHexExp::IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }


        void StdHexExp::v_IProductWRTDerivBase_MatOp(const int dir,
                                                    const Array<OneD, const NekDouble>& inarray,
                                                    Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0((dir==0)||(dir==1)||(dir==2),"input dir is out of range");

            int nq = GetTotPoints();
            MatrixType mtype = eIProductWRTDerivBase0;

            switch (dir)
            {
                case 0:
                    mtype = eIProductWRTDerivBase0;
                    break;
                case 1:
                    mtype = eIProductWRTDerivBase1;
                    break;
                case 2:
                    mtype = eIProductWRTDerivBase2;
                    break;
            }

            StdMatrixKey      iprodmatkey(mtype,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }


        void StdHexExp::v_IProductWRTDerivBase_SumFac(const int dir,
                                                     const Array<OneD, const NekDouble>& inarray,
                                                     Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0((dir==0)||(dir==1)||(dir==2),"input dir is out of range");

            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
            int    order1 = m_base[1]->GetNumModes();

            // If outarray > inarray then no need for temporary storage.
            Array<OneD, NekDouble> tmp = outarray;
            if (outarray.size() < inarray.size())
            {
                tmp = Array<OneD, NekDouble>(inarray.size());
            }

            // Need workspace for sumfackernel though
            Array<OneD, NekDouble> wsp(order0*nquad2*(nquad1+order1));

            // multiply by integration constants
            MultiplyByQuadratureMetric(inarray,tmp);

            // perform sum-factorisation
            switch (dir)
            {
                case 0:
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp,outarray,wsp,
                                                 false,true,true);
                    break;
                case 1:
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetDbdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp,outarray,wsp,
                                                 true,false,true);
                    break;
                case 2:
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetDbdata(),
                                                 tmp,outarray,wsp,
                                                 true,true,false);
                    break;
            }
        }

        void StdHexExp::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble>& xi,
                                      Array<OneD, NekDouble>& eta)
        {
            eta[0] = xi[0];
            eta[1] = xi[1];
            eta[2] = xi[2];
        }

        /**
         * @note for hexahedral expansions _base[0] (i.e. p) modes run fastest.
         */
        void StdHexExp::v_FillMode(const int mode,
                                Array<OneD, NekDouble> &outarray)
        {
            int    i,j;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            int   nquad2 = m_base[2]->GetNumPoints();

            Array<OneD, const NekDouble> base0  = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> base1  = m_base[1]->GetBdata();
            Array<OneD, const NekDouble> base2  = m_base[2]->GetBdata();

            int   btmp0 = m_base[0]->GetNumModes();
            int   btmp1 = m_base[1]->GetNumModes();
            int   mode2 = mode/(btmp0*btmp1);
            int   mode1 = (mode-mode2*btmp0*btmp1)/btmp0;
            int   mode0 = (mode-mode2*btmp0*btmp1)%btmp0;

            ASSERTL2(mode2 == (int)floor((1.0*mode)/(btmp0*btmp1)),
                     "Integer Truncation not Equiv to Floor");
            ASSERTL2(mode1 == (int)floor((1.0*mode-mode2*btmp0*btmp1)
                                /(btmp0*btmp1)),
                     "Integer Truncation not Equiv to Floor");
            ASSERTL2(m_ncoeffs <= mode,
                     "calling argument mode is larger than total expansion "
                     "order");

            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vcopy(nquad0,(NekDouble *)(base0.get() + mode0*nquad0),1,
                             &outarray[0]+i*nquad0, 1);
            }

            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad0; ++i)
                {
                    Vmath::Vmul(nquad1,(NekDouble *)(base1.get() + mode1*nquad1),1,
                                &outarray[0]+i+j*nquad0*nquad1, nquad0,
                                &outarray[0]+i+j*nquad0*nquad1, nquad0);
                }
            }

            for(i = 0; i < nquad2; i++)
            {
                Blas::Dscal(nquad0*nquad1,base2[mode2*nquad2+i],
                            &outarray[0]+i*nquad0*nquad1,1);
            }
        }


        int StdHexExp::v_GetNverts() const
        {
            return 8;
        }


        int StdHexExp::v_GetNedges() const
        {
            return 12;
        }


        int StdHexExp::v_GetNfaces() const
        {
            return 6;
        }


        LibUtilities::ShapeType StdHexExp::v_DetShapeType() const
        {
            return LibUtilities::eHexahedron;
        }


        int StdHexExp::v_NumBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int nmodes0 = m_base[0]->GetNumModes();
            int nmodes1 = m_base[1]->GetNumModes();
            int nmodes2 = m_base[2]->GetNumModes();

            return ( 2*( nmodes0*nmodes1 + nmodes0*nmodes2
                        + nmodes1*nmodes2)
                     - 4*( nmodes0 + nmodes1 + nmodes2 ) + 8 );
        }

        int StdHexExp::v_NumDGBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int nmodes0 = m_base[0]->GetNumModes();
            int nmodes1 = m_base[1]->GetNumModes();
            int nmodes2 = m_base[2]->GetNumModes();

            return  2*( nmodes0*nmodes1 + nmodes0*nmodes2
                        + nmodes1*nmodes2 );
        }

        int StdHexExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 11),"edge id is out of range");

            if((i == 0)||(i == 2)||(i == 8)||(i == 10))
            {
                return  GetBasisNumModes(0);
            }
            else if((i == 1)||(i == 3)||(i == 9)||(i == 11))
            {
                return  GetBasisNumModes(1);
            }
            else
            {
                return GetBasisNumModes(2);
            }
        }

        int StdHexExp::v_GetTotalEdgeIntNcoeffs() const
        {
            return 4*(GetBasisNumModes(0)+GetBasisNumModes(1)+GetBasisNumModes(2));
        }


        int StdHexExp::v_GetFaceNcoeffs(const int i) const
        {
            ASSERTL2((i >= 0) && (i <= 5), "face id is out of range");
            if((i == 0) || (i == 5))
            {
                return GetBasisNumModes(0)*GetBasisNumModes(1);
            }
            else if((i == 1) || (i == 3))
            {
                return GetBasisNumModes(0)*GetBasisNumModes(2);
            }
            else
            {
                return GetBasisNumModes(1)*GetBasisNumModes(2);
            }
        }


        int StdHexExp::v_GetFaceIntNcoeffs(const int i) const
        {
            ASSERTL2((i >= 0) && (i <= 5), "face id is out of range");
            if((i == 0) || (i == 5))
            {
                return (GetBasisNumModes(0)-2)*(GetBasisNumModes(1)-2);
            }
            else if((i == 1) || (i == 3))
            {
                return (GetBasisNumModes(0)-2)*(GetBasisNumModes(2)-2);
            }
            else
            {
                return (GetBasisNumModes(1)-2)*(GetBasisNumModes(2)-2);
            }

        }

        int StdHexExp::v_GetTotalFaceIntNcoeffs() const
        {
            return 2*((GetBasisNumModes(0)-2)*(GetBasisNumModes(1)-2)+
                      (GetBasisNumModes(0)-2)*(GetBasisNumModes(2)-2)+
                (GetBasisNumModes(1)-2)*(GetBasisNumModes(2)-2));
        }

        int StdHexExp::v_GetFaceNumPoints(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 5, "face id is out of range");

            if (i == 0 || i == 5)
            {
                return m_base[0]->GetNumPoints()*
                       m_base[1]->GetNumPoints();
            }
            else if (i == 1 || i == 3)
            {
                return m_base[0]->GetNumPoints()*
                       m_base[2]->GetNumPoints();
            }
            else
            {
                return m_base[1]->GetNumPoints()*
                       m_base[2]->GetNumPoints();
            }
        }

        LibUtilities::PointsKey StdHexExp::v_GetFacePointsKey(
            const int i, const int j) const
        {
            ASSERTL2(i >= 0 && i <= 5, "face id is out of range");
            ASSERTL2(j == 0 || j == 1, "face direction is out of range");

            if (i == 0 || i == 5)
            {
                return m_base[j]->GetPointsKey();
            }
            else if (i == 1 || i == 3)
            {
                return m_base[2*j]->GetPointsKey();
            }
            else
            {
                return m_base[j+1]->GetPointsKey();
            }
        }

        int StdHexExp::v_CalcNumberOfCoefficients(const std::vector<unsigned int> &nummodes, int &modes_offset)
        {
            int nmodes = nummodes[modes_offset]*nummodes[modes_offset+1]*nummodes[modes_offset+2];
            modes_offset += 3;

            return nmodes;
        }


        const LibUtilities::BasisKey StdHexExp::v_DetFaceBasisKey(
            const int i, const int k) const
        {
            ASSERTL2(i >= 0 && i <= 5, "face id is out of range");
            ASSERTL2(k >= 0 && k <= 1, "basis key id is out of range");

            int dir = k;
            switch(i)
            {
                case 0:
                case 5:
                    dir = k;
                    break;
                case 1:
                case 3:
                    dir = 2*k;
                    break;
                case 2:
                case 4:
                    dir = k+1;
                    break;
            }

            return EvaluateQuadFaceBasisKey(k,
                                            m_base[dir]->GetBasisType(),
                                            m_base[dir]->GetNumPoints(),
                                            m_base[dir]->GetNumModes());
        }

        LibUtilities::BasisType StdHexExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 11),"edge id is out of range");

            if((i == 0)||(i == 2)||(i==8)||(i==10))
            {
                return  GetBasisType(0);
            }
            else if((i == 1)||(i == 3)||(i == 9)||(i == 11))
            {
                return  GetBasisType(1);
            }
            else
            {
                return GetBasisType(2);
            }
        }

        void StdHexExp::v_GetCoords( Array<OneD, NekDouble> & xi_x,
                                Array<OneD, NekDouble> & xi_y,
                                Array<OneD, NekDouble> & xi_z)
        {
            Array<OneD, const NekDouble> eta_x = m_base[0]->GetZ();
            Array<OneD, const NekDouble> eta_y = m_base[1]->GetZ();
            Array<OneD, const NekDouble> eta_z = m_base[2]->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates:
            // eta --> xi
            for( int k = 0; k < Qz; ++k ) {
                for( int j = 0; j < Qy; ++j ) {
                    for( int i = 0; i < Qx; ++i ) {
                        int s = i + Qx*(j + Qy*k);
                        xi_x[s] = eta_x[i];
                        xi_y[s] = eta_y[j];
                        xi_z[s] = eta_z[k];

                    }
                }
            }
        }

        void StdHexExp::v_GetFaceNumModes(
                    const int                  fid,
                    const Orientation          faceOrient,
                    int &numModes0,
                    int &numModes1)
        {
            int nummodes [3] = {m_base[0]->GetNumModes(),
                                m_base[1]->GetNumModes(),
                                m_base[2]->GetNumModes()};
            switch(fid)
            {
            case 0:
            case 5:
                {
                    numModes0 = nummodes[0];
                    numModes1 = nummodes[1];
                }
                break;
            case 1:
            case 3:
                {
                    numModes0 = nummodes[0];
                    numModes1 = nummodes[2];
                }
                break;
            case 2:
            case 4:
                {
                    numModes0 = nummodes[1];
                    numModes1 = nummodes[2];
                }
                break;
            default:
                {
                    ASSERTL0(false,"fid out of range");
                }
                break;
            }

            if ( faceOrient >= eDir1FwdDir2_Dir2FwdDir1 )
            {
                std::swap(numModes0, numModes1);
            }
        }

        /**
         * Only for basis type Modified_A or GLL_LAGRANGE in all directions.
         */
        void StdHexExp::v_GetFaceToElementMap(
            const int                  fid,
            const Orientation          faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray,
            int                        P,
            int                        Q)
        {
            int i,j;
            int nummodesA=0, nummodesB=0;

            ASSERTL1(GetEdgeBasisType(0) == GetEdgeBasisType(1) &&
                     GetEdgeBasisType(0) == GetEdgeBasisType(2),
                     "Method only implemented if BasisType is indentical in "
                     "all directions");
            ASSERTL1(GetEdgeBasisType(0) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "Method only implemented for Modified_A or GLL_Lagrange BasisType");

            const int nummodes0 = m_base[0]->GetNumModes();
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nummodes2 = m_base[2]->GetNumModes();

            switch(fid)
            {
            case 0:
            case 5:
                nummodesA = nummodes0;
                nummodesB = nummodes1;
                break;
            case 1:
            case 3:
                nummodesA = nummodes0;
                nummodesB = nummodes2;
                break;
            case 2:
            case 4:
                nummodesA = nummodes1;
                nummodesB = nummodes2;
                break;
            default:
                ASSERTL0(false,"fid must be between 0 and 5");
            }

            bool CheckForZeroedModes = false;

            if (P == -1)
            {
                P = nummodesA;
                Q = nummodesB;
            }

            if((P != nummodesA)||(Q != nummodesB))
            {
                CheckForZeroedModes = true;
            }

            bool modified = (GetEdgeBasisType(0) == LibUtilities::eModified_A);
            int nFaceCoeffs = P*Q;

            if(maparray.size() != nFaceCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceCoeffs);
            }

            if(signarray.size() != nFaceCoeffs)
            {
                signarray = Array<OneD, int>(nFaceCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nFaceCoeffs, 1 );
            }

            Array<OneD, int> arrayindx(nFaceCoeffs);

            for(i = 0; i < Q; i++)
            {
                for(j = 0; j < P; j++)
                {
                    if( faceOrient < eDir1FwdDir2_Dir2FwdDir1 )
                    {
                        arrayindx[i*P+j] = i*P+j;
                    }
                    else
                    {
                        arrayindx[i*P+j] = j*Q+i;
                    }
                }
            }

            int offset = 0;
            int jump1 = 1;
            int jump2 = 1;

            switch(fid)
            {
                case 5:
                {
                    if (modified)
                    {
                        offset = nummodes0*nummodes1;
                    }
                    else
                    {
                        offset = (nummodes2-1)*nummodes0*nummodes1;
                        jump1 = nummodes0;
                    }
                }
                /* Falls through. */
                case 0:
                {
                    jump1 = nummodes0;
                    break;
                }
                case 3:
                {
                    if (modified)
                    {
                        offset = nummodes0;
                    }
                    else
                    {
                        offset = nummodes0*(nummodes1-1);
                        jump1 = nummodes0*nummodes1;
                    }
                }
                /* Falls through. */
                case 1:
                {
                    jump1 = nummodes0*nummodes1;
                    break;
                }
                case 2:
                {
                    if (modified)
                    {
                        offset = 1;
                    }
                    else
                    {
                        offset = nummodes0-1;
                        jump1 = nummodes0*nummodes1;
                        jump2 = nummodes0;

                    }
                }
                /* Falls through. */
                case 4:
                {
                    jump1 = nummodes0*nummodes1;
                    jump2 = nummodes0;
                    break;
                }
                default:
                    ASSERTL0(false,"fid must be between 0 and 5");
            }

            for(i = 0; i < Q; i++)
            {
                for(j = 0; j < P; j++)
                {
                    maparray[ arrayindx[i*P+j] ]
                    = i*jump1 + j*jump2 + offset;
                }
            }


            if(CheckForZeroedModes)
            {
                if(modified)
                {
                    // zero signmap and set maparray to zero if elemental
                    // modes are not as large as face modesl
                    for(i = 0; i < nummodesB; i++)
                    {
                        for(j = nummodesA; j < P; j++)
                        {
                            signarray[arrayindx[i*P+j]] = 0.0;
                            maparray[arrayindx[i*P+j]]  = maparray[0];
                        }
                    }

                    for(i = nummodesB; i < Q; i++)
                    {
                        for(j = 0; j < P; j++)
                        {
                            signarray[arrayindx[i*P+j]] = 0.0;
                            maparray[arrayindx[i*P+j]]  = maparray[0];
                        }
                    }
                }
                else
                {
                    ASSERTL0(false,"Different trace space face dimention and element face dimention not possible for GLL-Lagrange bases");
                }
            }

            if( (faceOrient==eDir1FwdDir1_Dir2BwdDir2) ||
                (faceOrient==eDir1BwdDir1_Dir2BwdDir2) ||
                (faceOrient==eDir1BwdDir2_Dir2FwdDir1) ||
                (faceOrient==eDir1BwdDir2_Dir2BwdDir1) )
            {
                if(faceOrient<eDir1FwdDir2_Dir2FwdDir1)
                {
                    if (modified)
                    {
                        for(i = 3; i < Q; i+=2)
                        {
                            for(j = 0; j < P; j++)
                            {
                                signarray[ arrayindx[i*P+j] ] *= -1;
                            }
                        }

                        for(i = 0; i < P; i++)
                        {
                            swap( maparray[i] , maparray[i+P] );
                            swap( signarray[i] , signarray[i+P] );
                        }

                    }
                    else
                    {
                        for(i = 0; i < P; i++)
                        {
                            for(j = 0; j < Q/2; j++)
                            {
                                swap( maparray[i + j*P],
                                     maparray[i+P*Q
                                              -P -j*P] );
                                swap( signarray[i + j*P],
                                     signarray[i+P*Q
                                              -P -j*P]);
                            }
                        }
                    }
                }
                else
                {
                    if (modified)
                    {
                        for(i = 0; i < Q; i++)
                        {
                            for(j = 3; j < P; j+=2)
                            {
                                signarray[ arrayindx[i*P+j] ] *= -1;
                            }
                        }

                        for(i = 0; i < Q; i++)
                        {
                            swap( maparray[i] , maparray[i+Q] );
                            swap( signarray[i] , signarray[i+Q] );
                        }

                    }
                    else
                    {
                        for(i = 0; i < P; i++)
                        {
                            for(j = 0; j < Q/2; j++)
                            {
                                swap( maparray[i*Q + j],
                                     maparray[i*Q + Q -1 -j]);
                                swap( signarray[i*Q + j],
                                     signarray[i*Q + Q -1 -j]);
                            }
                        }
                    }
                }
            }

            if( (faceOrient==eDir1BwdDir1_Dir2FwdDir2) ||
                (faceOrient==eDir1BwdDir1_Dir2BwdDir2) ||
                (faceOrient==eDir1FwdDir2_Dir2BwdDir1) ||
                (faceOrient==eDir1BwdDir2_Dir2BwdDir1) )
            {
                if(faceOrient<eDir1FwdDir2_Dir2FwdDir1)
                {
                    if (modified)
                    {
                        for(i = 0; i < Q; i++)
                        {
                            for(j = 3; j < P; j+=2)
                            {
                                signarray[ arrayindx[i*P+j] ] *= -1;
                            }
                        }

                        for(i = 0; i < Q; i++)
                        {
                            swap( maparray[i*P],
                                 maparray[i*P+1]);
                            swap( signarray[i*P],
                                 signarray[i*P+1]);
                        }
                    }
                    else
                    {
                        for(i = 0; i < Q; i++)
                        {
                            for(j = 0; j < P/2; j++)
                            {
                                swap( maparray[i*P + j],
                                     maparray[i*P + P -1 -j]);
                                swap( signarray[i*P + j],
                                     signarray[i*P + P -1 -j]);
                            }
                        }
                    }



                }
                else
                {
                    if (modified)
                    {
                        for(i = 3; i < Q; i+=2)
                        {
                            for(j = 0; j < P; j++)
                            {
                                signarray[ arrayindx[i*P+j] ] *= -1;
                            }
                        }

                        for(i = 0; i < P; i++)
                        {
                            swap( maparray[i*Q],
                                 maparray[i*Q+1]);
                            swap( signarray[i*Q],
                                 signarray[i*Q+1]);
                        }
                    }
                    else
                    {
                        for(i = 0; i < Q; i++)
                        {
                            for(j = 0; j < P/2; j++)
                            {
                                swap( maparray[i + j*Q] ,
                                     maparray[i+P*Q - Q -j*Q] );
                                swap( signarray[i + j*Q] ,
                                     signarray[i+P*Q - Q -j*Q] );
                            }
                        }
                    }
                }
            }
        }



        /**
         * Expansions in each of the three dimensions must be of type
         * LibUtilities#eModified_A or LibUtilities#eGLL_Lagrange.
         *
         * @param   localVertexId   ID of vertex (0..7)
         * @returns Position of vertex in local numbering scheme.
         */
        int StdHexExp::v_GetVertexMap(const int localVertexId, bool useCoeffPacking)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            ASSERTL1((localVertexId>=0)&&(localVertexId<8),
                     "local vertex id must be between 0 and 7");

            int p = 0;
            int q = 0;
            int r = 0;

            // Retrieve the number of modes in each dimension.
            int nummodes [3] = {m_base[0]->GetNumModes(),
                                m_base[1]->GetNumModes(),
                                m_base[2]->GetNumModes()};

            if(useCoeffPacking == true) // follow packing of coefficients i.e q,r,p
            {
                if(localVertexId > 3)
                {
                    if( GetBasisType(2) == LibUtilities::eGLL_Lagrange)
                    {
                        r = nummodes[2]-1;
                    }
                    else
                    {
                        r = 1;
                    }
                }

                switch(localVertexId % 4)
                {
                case 0:
                    break;
                case 1:
                    {
                        if( GetBasisType(0) == LibUtilities::eGLL_Lagrange)
                        {
                            p = nummodes[0]-1;
                        }
                        else
                        {
                            p = 1;
                        }
                    }
                    break;
                case 2:
                    {
                        if( GetBasisType(1) == LibUtilities::eGLL_Lagrange)
                        {
                            q = nummodes[1]-1;
                        }
                        else
                        {
                            q = 1;
                        }
                    }
                    break;
                case 3:
                    {
                        if( GetBasisType(1) == LibUtilities::eGLL_Lagrange)
                        {
                            p = nummodes[0]-1;
                            q = nummodes[1]-1;
                        }
                        else
                        {
                            p = 1;
                            q = 1;
                        }
                    }
                    break;
                }
            }
            else
            {
                // Right face (vertices 1,2,5,6)
                if( (localVertexId % 4) % 3 > 0 )
                {
                    if( GetBasisType(0) == LibUtilities::eGLL_Lagrange)
                    {
                        p = nummodes[0]-1;
                    }
                    else
                    {
                        p = 1;
                    }
                }

                // Back face (vertices 2,3,6,7)
                if( localVertexId % 4 > 1 )
                {
                    if( GetBasisType(1) == LibUtilities::eGLL_Lagrange)
                    {
                        q = nummodes[1]-1;
                    }
                    else
                    {
                        q = 1;
                    }
                }

                // Top face (vertices 4,5,6,7)
                if( localVertexId > 3)
                {
                    if( GetBasisType(2) == LibUtilities::eGLL_Lagrange)
                    {
                        r = nummodes[2]-1;
                    }
                    else
                    {
                        r = 1;
                    }
                }
            }
            // Compute the local number.
            return r*nummodes[0]*nummodes[1] + q*nummodes[0] + p;
        }


        /**
         * @param   eid         The edge to compute the numbering for.
         * @param   edgeOrient  Orientation of the edge.
         * @param   maparray    Storage for computed mapping array.
         * @param   signarray   ?
         */
        void StdHexExp::v_GetEdgeInteriorMap(const int eid,
                                const Orientation edgeOrient,
                                Array<OneD, unsigned int> &maparray,
                                Array<OneD, int> &signarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            ASSERTL1((eid>=0)&&(eid<12),
                     "local edge id must be between 0 and 11");

            int nEdgeIntCoeffs = GetEdgeNcoeffs(eid)-2;

            if(maparray.size()!=nEdgeIntCoeffs)
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

            int nummodes [3] = {m_base[0]->GetNumModes(),
                                m_base[1]->GetNumModes(),
                                m_base[2]->GetNumModes()};

            const LibUtilities::BasisType bType [3] = {GetBasisType(0),
                                                       GetBasisType(1),
                                                       GetBasisType(2)};

            bool reverseOrdering = false;
            bool signChange = false;

            int IdxRange [3][2] = {{0,0},{0,0},{0,0}};

            switch(eid)
            {
            case 0:
            case 1:
            case 2:
            case 3:
                {
                    IdxRange[2][0] = 0;
                    IdxRange[2][1] = 1;
                }
                break;
            case 8:
            case 9:
            case 10:
            case 11:
                {
                    if( bType[2] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[2][0] = nummodes[2] - 1;
                        IdxRange[2][1] = nummodes[2];
                    }
                    else
                    {
                        IdxRange[2][0] = 1;
                        IdxRange[2][1] = 2;
                    }
                }
                break;
            case 4:
            case 5:
            case 6:
            case 7:
                {
                    if( bType[2] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[2][0] = 1;
                        IdxRange[2][1] = nummodes[2] - 1;

                        if(edgeOrient==eBackwards)
                        {
                            reverseOrdering = true;
                        }
                    }
                    else
                    {
                        IdxRange[2][0] = 2;
                        IdxRange[2][1] = nummodes[2];

                        if(edgeOrient==eBackwards)
                        {
                            signChange = true;
                        }
                    }
                }
                break;
            }

            switch(eid)
            {
            case 0:
            case 4:
            case 5:
            case 8:
                {
                    IdxRange[1][0] = 0;
                    IdxRange[1][1] = 1;
                }
                break;
            case 2:
            case 6:
            case 7:
            case 10:
                {
                    if( bType[1] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[1][0] = nummodes[1] - 1;
                        IdxRange[1][1] = nummodes[1];
                    }
                    else
                    {
                        IdxRange[1][0] = 1;
                        IdxRange[1][1] = 2;
                    }
                }
                break;
            case 1:
            case 9:
                {
                    if( bType[1] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[1][0] = 1;
                        IdxRange[1][1] = nummodes[1] - 1;

                        if(edgeOrient==eBackwards)
                        {
                            reverseOrdering = true;
                        }
                    }
                    else
                    {
                        IdxRange[1][0] = 2;
                        IdxRange[1][1] = nummodes[1];

                        if(edgeOrient==eBackwards)
                        {
                            signChange = true;
                        }
                    }
                }
                break;
            case 3:
            case 11:
                {
                    if( bType[1] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[1][0] = 1;
                        IdxRange[1][1] = nummodes[1] - 1;

                        if(edgeOrient==eForwards)
                        {
                            reverseOrdering = true;
                        }
                    }
                    else
                    {
                        IdxRange[1][0] = 2;
                        IdxRange[1][1] = nummodes[1];

                        if(edgeOrient==eForwards)
                        {
                            signChange = true;
                        }
                    }
                }
                break;
            }

            switch(eid)
            {
            case 3:
            case 4:
            case 7:
            case 11:
                {
                    IdxRange[0][0] = 0;
                    IdxRange[0][1] = 1;
                }
                break;
            case 1:
            case 5:
            case 6:
            case 9:
                {
                    if( bType[0] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[0][0] = nummodes[0] - 1;
                        IdxRange[0][1] = nummodes[0];
                    }
                    else
                    {
                        IdxRange[0][0] = 1;
                        IdxRange[0][1] = 2;
                    }
                }
                break;
            case 0:
            case 8:
                {
                    if( bType[0] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[0][0] = 1;
                        IdxRange[0][1] = nummodes[0] - 1;

                        if(edgeOrient==eBackwards)
                        {
                            reverseOrdering = true;
                        }
                    }
                    else
                    {
                        IdxRange[0][0] = 2;
                        IdxRange[0][1] = nummodes[0];

                        if(edgeOrient==eBackwards)
                        {
                            signChange = true;
                        }
                    }
                }
                break;
            case 2:
            case 10:
                {
                    if( bType[0] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[0][0] = 1;
                        IdxRange[0][1] = nummodes[0] - 1;

                        if(edgeOrient==eForwards)
                        {
                            reverseOrdering = true;
                        }
                    }
                    else
                    {
                        IdxRange[0][0] = 2;
                        IdxRange[0][1] = nummodes[0];

                        if(edgeOrient==eForwards)
                        {
                            signChange = true;
                        }
                    }
                }
                break;
            }

            int p,q,r;
            int cnt = 0;

            for(r = IdxRange[2][0]; r < IdxRange[2][1]; r++)
            {
                for(q = IdxRange[1][0]; q < IdxRange[1][1]; q++)
                {
                    for(p = IdxRange[0][0]; p < IdxRange[0][1]; p++)
                    {
                        maparray[cnt++]
                                = r*nummodes[0]*nummodes[1] + q*nummodes[0] + p;
                    }
                }
            }

            if( reverseOrdering )
            {
                reverse( maparray.get() , maparray.get()+nEdgeIntCoeffs );
            }

            if( signChange )
            {
                for(p = 1; p < nEdgeIntCoeffs; p+=2)
                {
                    signarray[p] = -1;
                }
            }
        }


        /**
         * Generate mapping describing which elemental modes lie on the
         * interior of a given face. Accounts for face orientation.
         */
        void StdHexExp::v_GetFaceInteriorMap(const int fid,
                                const Orientation faceOrient,
                                Array<OneD, unsigned int> &maparray,
                                Array<OneD, int>& signarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            ASSERTL1((fid>=0)&&(fid<6),
                     "local face id must be between 0 and 5");

            int nFaceIntCoeffs = GetFaceIntNcoeffs(fid);

            if(maparray.size()!=nFaceIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceIntCoeffs);
            }

            if(signarray.size() != nFaceIntCoeffs)
            {
                signarray = Array<OneD, int>(nFaceIntCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nFaceIntCoeffs, 1 );
            }

            int nummodes [3] = {m_base[0]->GetNumModes(),
                                m_base[1]->GetNumModes(),
                                m_base[2]->GetNumModes()};

            const LibUtilities::BasisType bType [3] = {GetBasisType(0),
                                                       GetBasisType(1),
                                                       GetBasisType(2)};

            int nummodesA = 0;
            int nummodesB = 0;

            // Determine the number of modes in face directions A & B based
            // on the face index given.
            switch(fid)
            {
            case 0:
            case 5:
                {
                    nummodesA = nummodes[0];
                    nummodesB = nummodes[1];
                }
                break;
            case 1:
            case 3:
                {
                    nummodesA = nummodes[0];
                    nummodesB = nummodes[2];
                }
                break;
            case 2:
            case 4:
                {
                    nummodesA = nummodes[1];
                    nummodesB = nummodes[2];
                }
            }

            int i,j;
            Array<OneD, int> arrayindx(nFaceIntCoeffs);

            // Create a mapping array to account for transposition of the
            // coordinates due to face orientation.
            for(i = 0; i < (nummodesB-2); i++)
            {
                for(j = 0; j < (nummodesA-2); j++)
                {
                    if( faceOrient < eDir1FwdDir2_Dir2FwdDir1 )
                    {
                        arrayindx[i*(nummodesA-2)+j] = i*(nummodesA-2)+j;
                    }
                    else
                    {
                        arrayindx[i*(nummodesA-2)+j] = j*(nummodesB-2)+i;
                    }
                }
            }

            int IdxRange [3][2];
            int Incr[3];

            Array<OneD, int> sign0(nummodes[0], 1);
            Array<OneD, int> sign1(nummodes[1], 1);
            Array<OneD, int> sign2(nummodes[2], 1);

            // Set the upper and lower bounds, and increment for the faces
            // involving the first coordinate direction.
            switch(fid)
            {
            case 0: // bottom face
                {
                    IdxRange[2][0] = 0;
                    IdxRange[2][1] = 1;
                    Incr[2] = 1;
                }
                break;
            case 5: // top face
                {
                    if( bType[2] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[2][0] = nummodes[2] - 1;
                        IdxRange[2][1] = nummodes[2];
                        Incr[2] = 1;
                    }
                    else
                    {
                        IdxRange[2][0] = 1;
                        IdxRange[2][1] = 2;
                        Incr[2] = 1;
                    }

                }
                break;
            default: // all other faces
                {
                    if( bType[2] == LibUtilities::eGLL_Lagrange)
                    {
                        if( ((int) (faceOrient-eDir1FwdDir1_Dir2FwdDir2)) % 2 )
                        {
                            IdxRange[2][0] = nummodes[2] - 2;
                            IdxRange[2][1] = 0;
                            Incr[2] = -1;

                        }
                        else
                        {
                            IdxRange[2][0] = 1;
                            IdxRange[2][1] = nummodes[2] - 1;
                            Incr[2] = 1;
                        }
                    }
                    else
                    {
                        IdxRange[2][0] = 2;
                        IdxRange[2][1] = nummodes[2];
                        Incr[2] = 1;

                        if( ((int) (faceOrient-eDir1FwdDir1_Dir2FwdDir2)) % 2 )
                        {
                            for(i = 3; i < nummodes[2]; i+=2)
                            {
                                sign2[i] = -1;
                            }
                        }
                    }
                }
            }

            // Set the upper and lower bounds, and increment for the faces
            // involving the second coordinate direction.
            switch(fid)
            {
            case 1:
                {
                    IdxRange[1][0] = 0;
                    IdxRange[1][1] = 1;
                    Incr[1] = 1;
                }
                break;
            case 3:
                {
                    if( bType[1] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[1][0] = nummodes[1] - 1;
                        IdxRange[1][1] = nummodes[1];
                        Incr[1] = 1;
                    }
                    else
                    {
                        IdxRange[1][0] = 1;
                        IdxRange[1][1] = 2;
                        Incr[1] = 1;
                    }
                }
                break;
            case 0:
            case 5:
                {
                    if( bType[1] == LibUtilities::eGLL_Lagrange)
                    {
                        if( ((int) (faceOrient-eDir1FwdDir1_Dir2FwdDir2)) % 2 )
                        {
                            IdxRange[1][0] = nummodes[1] - 2;
                            IdxRange[1][1] = 0;
                            Incr[1] = -1;

                        }
                        else
                        {
                            IdxRange[1][0] = 1;
                            IdxRange[1][1] = nummodes[1] - 1;
                            Incr[1] = 1;
                        }
                    }
                    else
                    {
                        IdxRange[1][0] = 2;
                        IdxRange[1][1] = nummodes[1];
                        Incr[1] = 1;

                        if( ((int) (faceOrient-eDir1FwdDir1_Dir2FwdDir2)) % 2 )
                        {
                            for(i = 3; i < nummodes[1]; i+=2)
                            {
                                sign1[i] = -1;
                            }
                        }
                    }
                }
                break;
            default: // case2: case4:
                {
                    if( bType[1] == LibUtilities::eGLL_Lagrange)
                    {
                        if( ((int) (faceOrient-eDir1FwdDir1_Dir2FwdDir2)) % 4 > 1 )
                        {
                            IdxRange[1][0] = nummodes[1] - 2;
                            IdxRange[1][1] = 0;
                            Incr[1] = -1;

                        }
                        else
                        {
                            IdxRange[1][0] = 1;
                            IdxRange[1][1] = nummodes[1] - 1;
                            Incr[1] = 1;
                        }
                    }
                    else
                    {
                        IdxRange[1][0] = 2;
                        IdxRange[1][1] = nummodes[1];
                        Incr[1] = 1;

                        if( ((int) (faceOrient-eDir1FwdDir1_Dir2FwdDir2)) % 4 > 1 )
                        {
                            for(i = 3; i < nummodes[1]; i+=2)
                            {
                                sign1[i] = -1;
                            }
                        }
                    }
                }
            }

            switch(fid)
            {
            case 4:
                {
                    IdxRange[0][0] = 0;
                    IdxRange[0][1] = 1;
                    Incr[0] = 1;
                }
                break;
            case 2:
                {
                    if( bType[0] == LibUtilities::eGLL_Lagrange)
                    {
                        IdxRange[0][0] = nummodes[0] - 1;
                        IdxRange[0][1] = nummodes[0];
                        Incr[0] = 1;
                    }
                    else
                    {
                        IdxRange[0][0] = 1;
                        IdxRange[0][1] = 2;
                        Incr[0] = 1;
                    }
                }
                break;
            default:
                {
                    if( bType[0] == LibUtilities::eGLL_Lagrange)
                    {
                        if( ((int) (faceOrient-eDir1FwdDir1_Dir2FwdDir2)) % 4 > 1 )
                        {
                            IdxRange[0][0] = nummodes[0] - 2;
                            IdxRange[0][1] = 0;
                            Incr[0] = -1;

                        }
                        else
                        {
                            IdxRange[0][0] = 1;
                            IdxRange[0][1] = nummodes[0] - 1;
                            Incr[0] = 1;
                        }
                    }
                    else
                    {
                        IdxRange[0][0] = 2;
                        IdxRange[0][1] = nummodes[0];
                        Incr[0] = 1;

                        if( ((int) (faceOrient-eDir1FwdDir1_Dir2FwdDir2)) % 4 > 1 )
                        {
                            for(i = 3; i < nummodes[0]; i+=2)
                            {
                                sign0[i] = -1;
                            }
                        }
                    }
                }
            }

            int p,q,r;
            int cnt = 0;

            for(r = IdxRange[2][0]; r != IdxRange[2][1]; r+=Incr[2])
            {
                for(q = IdxRange[1][0]; q != IdxRange[1][1]; q+=Incr[1])
                {
                    for(p = IdxRange[0][0]; p != IdxRange[0][1]; p+=Incr[0])
                    {
                        maparray [ arrayindx[cnt  ] ]
                                = r*nummodes[0]*nummodes[1] + q*nummodes[0] + p;
                        signarray[ arrayindx[cnt++] ]
                                = sign0[p] * sign1[q] * sign2[r];
                    }
                }
            }
        }


        /**
         * @param   outarray    Storage area for computed map.
         */
        void StdHexExp::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int i;
            int nummodes [3] = {m_base[0]->GetNumModes(),
                                m_base[1]->GetNumModes(),
                                m_base[2]->GetNumModes()};

            int nIntCoeffs = m_ncoeffs - NumBndryCoeffs();

            if(outarray.size() != nIntCoeffs)
            {
                outarray = Array<OneD, unsigned int>(nIntCoeffs);
            }

            const LibUtilities::BasisType Btype [3] = {GetBasisType(0),
                                                       GetBasisType(1),
                                                       GetBasisType(2)};

            int p,q,r;
            int cnt = 0;

            int IntIdx [3][2];

            for(i = 0; i < 3; i++)
            {
                if( Btype[i] == LibUtilities::eModified_A)
                {
                    IntIdx[i][0]  = 2;
                    IntIdx[i][1]  = nummodes[i];
                }
                else
                {
                    IntIdx[i][0]  = 1;
                    IntIdx[i][1]  = nummodes[i]-1;
                }
            }

            for(r = IntIdx[2][0]; r < IntIdx[2][1]; r++)
            {
                for( q = IntIdx[1][0]; q < IntIdx[1][1]; q++)
                {
                    for( p = IntIdx[0][0]; p < IntIdx[0][1]; p++)
                    {
                        outarray[cnt++] = r*nummodes[0]*nummodes[1] +
                            q*nummodes[0] + p;
                    }
                }
            }
        }


        /**
         * @param   outarray    Storage for computed map.
         */
        void StdHexExp::v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int i;
            int nummodes [3] = {m_base[0]->GetNumModes(),
                                m_base[1]->GetNumModes(),
                                m_base[2]->GetNumModes()};

            int nBndCoeffs = NumBndryCoeffs();

            if(outarray.size()!=nBndCoeffs)
            {
                outarray = Array<OneD, unsigned int>(nBndCoeffs);
            }

            const LibUtilities::BasisType Btype [3] = {GetBasisType(0),
                                                       GetBasisType(1),
                                                       GetBasisType(2)};

            int p,q,r;
            int cnt = 0;

            int BndIdx [3][2];
            int IntIdx [3][2];

            for(i = 0; i < 3; i++)
            {
                BndIdx[i][0] = 0;

                if( Btype[i] == LibUtilities::eModified_A)
                {
                    BndIdx[i][1] = 1;
                    IntIdx[i][0]  = 2;
                    IntIdx[i][1]  = nummodes[i];
                }
                else
                {
                    BndIdx[i][1] = nummodes[i]-1;
                    IntIdx[i][0]  = 1;
                    IntIdx[i][1]  = nummodes[i]-1;
                }
            }


            for(i = 0; i < 2; i++)
            {
                r = BndIdx[2][i];
                for( q = 0; q < nummodes[1]; q++)
                {
                    for( p = 0; p < nummodes[0]; p++)
                    {
                        outarray[cnt++] = r*nummodes[0]*nummodes[1]+q*nummodes[0] + p;
                    }
                }
            }

            for(r = IntIdx[2][0]; r < IntIdx[2][1]; r++)
            {
                for( i = 0; i < 2; i++)
                {
                    q = BndIdx[1][i];
                    for( p = 0; p < nummodes[0]; p++)
                    {
                        outarray[cnt++] = r*nummodes[0]*nummodes[1] +
                            q*nummodes[0] + p;
                    }
                }

                for( q = IntIdx[1][0]; q < IntIdx[1][1]; q++)
                {
                    for( i = 0; i < 2; i++)
                    {
                        p = BndIdx[0][i];
                        outarray[cnt++] = r*nummodes[0]*nummodes[1] +
                            q*nummodes[0] + p;
                    }
                }
            }

            sort(outarray.get(), outarray.get() + nBndCoeffs);
        }

        DNekMatSharedPtr StdHexExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            return StdExpansion::CreateGeneralMatrix(mkey);
        }


        DNekMatSharedPtr StdHexExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return StdExpansion::CreateGeneralMatrix(mkey);
        }


        void StdHexExp::v_MassMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
        }


        void StdHexExp::v_LaplacianMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdHexExp::v_LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
        }


        void StdHexExp::v_LaplacianMatrixOp(const int k1, const int k2,
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdExpansion::LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,
                                                        mkey);
        }


        void StdHexExp::v_WeakDerivMatrixOp(const int i,
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,
                                                        mkey);
        }

        void StdHexExp::v_HelmholtzMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            StdHexExp::v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
        }


        void StdHexExp::v_GeneralMatrixOp_MatOp(
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


        void StdHexExp::v_MultiplyByStdQuadratureMetric(const Array<OneD, const NekDouble>& inarray,
                                                 Array<OneD, NekDouble> &outarray)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            int nq01 = nquad0*nquad1;
            int nq12 = nquad1*nquad2;

            const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();
            const Array<OneD, const NekDouble>& w2 = m_base[2]->GetW();

            for(i = 0; i < nq12; ++i)
            {
                Vmath::Vmul(nquad0, inarray.get()+i*nquad0, 1,
                            w0.get(), 1, outarray.get()+i*nquad0,1);
            }

            for(i = 0; i < nq12; ++i)
            {
                Vmath::Smul(nquad0, w1[i%nquad1], outarray.get()+i*nquad0, 1,
                            outarray.get()+i*nquad0, 1);
            }

            for(i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nq01, w2[i], outarray.get()+i*nq01, 1,
                            outarray.get()+i*nq01, 1);
            }
        }

        void StdHexExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                              const StdMatrixKey &mkey)
        {
            // Generate an orthonogal expansion
            int qa = m_base[0]->GetNumPoints();
            int qb = m_base[1]->GetNumPoints();
            int qc = m_base[2]->GetNumPoints();
            int nmodes_a = m_base[0]->GetNumModes();
            int nmodes_b = m_base[1]->GetNumModes();
            int nmodes_c = m_base[2]->GetNumModes();
            // Declare orthogonal basis.
            LibUtilities::PointsKey pa(qa,m_base[0]->GetPointsType());
            LibUtilities::PointsKey pb(qb,m_base[1]->GetPointsType());
            LibUtilities::PointsKey pc(qc,m_base[2]->GetPointsType());

            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A,nmodes_a,pa);
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A,nmodes_b,pb);
            LibUtilities::BasisKey Bc(LibUtilities::eOrtho_A,nmodes_c,pc);
            StdHexExp OrthoExp(Ba,Bb,Bc);

            Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
            int i,j,k,cnt=0;

            // project onto modal  space.
            OrthoExp.FwdTrans(array,orthocoeffs);

            if(mkey.ConstFactorExists(eFactorSVVPowerKerDiffCoeff))
            {
                // Rodrigo's power kernel
                NekDouble cutoff = mkey.GetConstFactor(eFactorSVVCutoffRatio);
                NekDouble  SvvDiffCoeff  =
                    mkey.GetConstFactor(eFactorSVVPowerKerDiffCoeff)*
                    mkey.GetConstFactor(eFactorSVVDiffCoeff);

                for(int i = 0; i < nmodes_a; ++i)
                {
                    for(int j = 0; j < nmodes_b; ++j)
                    {
                        NekDouble fac1 = std::max(
                                   pow((1.0*i)/(nmodes_a-1),cutoff*nmodes_a),
                                   pow((1.0*j)/(nmodes_b-1),cutoff*nmodes_b));

                        for(int k = 0; k < nmodes_c; ++k)
                        {
                            NekDouble fac = std::max(fac1,
                                     pow((1.0*k)/(nmodes_c-1),cutoff*nmodes_c));

                            orthocoeffs[cnt]
                                *= SvvDiffCoeff * fac;
                            cnt++;
                        }
                    }
                }
            }
            else if(mkey.ConstFactorExists(eFactorSVVDGKerDiffCoeff))  // Rodrigo/Mansoor's DG Kernel
            {
                NekDouble  SvvDiffCoeff  =
                    mkey.GetConstFactor(eFactorSVVDGKerDiffCoeff)*
                    mkey.GetConstFactor(eFactorSVVDiffCoeff);

                int max_abc = max(nmodes_a-kSVVDGFiltermodesmin,
                                  nmodes_b-kSVVDGFiltermodesmin);
                max_abc = max(max_abc, nmodes_c-kSVVDGFiltermodesmin);
                // clamp max_abc
                max_abc = max(max_abc,0);
                max_abc = min(max_abc,kSVVDGFiltermodesmax-kSVVDGFiltermodesmin);

                for(int i = 0; i < nmodes_a; ++i)
                {
                    for(int j = 0; j < nmodes_b; ++j)
                    {
                        int maxij = max(i,j);

                        for(int k = 0; k < nmodes_c; ++k)
                        {
                            int maxijk = max(maxij,k);
                            maxijk = min(maxijk,kSVVDGFiltermodesmax-1);

                            orthocoeffs[cnt] *= SvvDiffCoeff *
                                kSVVDGFilter[max_abc][maxijk];
                            cnt++;
                        }
                    }
                }
            }
            else
            {

                int cutoff = (int) (mkey.GetConstFactor(eFactorSVVCutoffRatio)*min(nmodes_a,nmodes_b));
                NekDouble  SvvDiffCoeff  = mkey.GetConstFactor(eFactorSVVDiffCoeff);
                //  Filter just trilinear space
                int nmodes = max(nmodes_a,nmodes_b);
                nmodes = max(nmodes,nmodes_c);

                Array<OneD, NekDouble> fac(nmodes,1.0);
                for(j = cutoff; j < nmodes; ++j)
                {
                    fac[j] = fabs((j-nmodes)/((NekDouble) (j-cutoff+1.0)));
                    fac[j] *= fac[j]; //added this line to conform with equation
                }

                for(i = 0; i < nmodes_a; ++i)
                {
                    for(j = 0; j < nmodes_b; ++j)
                    {
                        for(k =  0; k < nmodes_c; ++k)
                        {
                            if((i >= cutoff)||(j >= cutoff)||(k >= cutoff))
                            {
                                orthocoeffs[i*nmodes_a*nmodes_b +
                                            j*nmodes_c + k] *=
                                    (SvvDiffCoeff*exp(-(fac[i]+fac[j]+fac[k])));
                            }
                            else
                            {
                                orthocoeffs[i*nmodes_a*nmodes_b + j*nmodes_c + k] *= 0.0;
                            }
                        }
                    }
                }
            }

            // backward transform to physical space
            OrthoExp.BwdTrans(orthocoeffs,array);
        }

        void StdHexExp::v_ExponentialFilter(
                                          Array<OneD, NekDouble> &array,
                                    const NekDouble        alpha,
                                    const NekDouble        exponent,
                                    const NekDouble        cutoff)
        {
            // Generate an orthogonal expansion
            int qa      = m_base[0]->GetNumPoints();
            int qb      = m_base[1]->GetNumPoints();
            int qc      = m_base[2]->GetNumPoints();
            int nmodesA = m_base[0]->GetNumModes();
            int nmodesB = m_base[1]->GetNumModes();
            int nmodesC = m_base[2]->GetNumModes();
            int P  = nmodesA - 1;
            int Q  = nmodesB - 1;
            int R  = nmodesC - 1;

            // Declare orthogonal basis.
            LibUtilities::PointsKey pa(qa,m_base[0]->GetPointsType());
            LibUtilities::PointsKey pb(qb,m_base[1]->GetPointsType());
            LibUtilities::PointsKey pc(qc,m_base[2]->GetPointsType());

            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, nmodesA, pa);
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, nmodesB, pb);
            LibUtilities::BasisKey Bc(LibUtilities::eOrtho_A, nmodesC, pc);
            StdHexExp OrthoExp(Ba,Bb,Bc);

            // Cutoff
            int Pcut = cutoff*P;
            int Qcut = cutoff*Q;
            int Rcut = cutoff*R;

            // Project onto orthogonal space.
            Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
            OrthoExp.FwdTrans(array,orthocoeffs);

            //
            NekDouble fac, fac1, fac2, fac3;
            int index = 0;
            for(int i = 0; i < nmodesA; ++i)
            {
                for(int j = 0; j < nmodesB; ++j)
                {
                    for(int k = 0; k < nmodesC; ++k, ++index)
                    {
                        //to filter out only the "high-modes"
                        if(i > Pcut || j > Qcut || k > Rcut)
                        {
                            fac1 = (NekDouble) (i - Pcut)/( (NekDouble)(P - Pcut) );
                            fac2 = (NekDouble) (j - Qcut)/( (NekDouble)(Q - Qcut) );
                            fac3 = (NekDouble) (k - Rcut)/( (NekDouble)(R - Rcut) );
                            fac  = max( max(fac1, fac2), fac3);
                            fac  = pow(fac, exponent);
                            orthocoeffs[index] *= exp(-alpha*fac);
                        }
                    }
                }
            }

            // backward transform to physical space
            OrthoExp.BwdTrans(orthocoeffs,array);
        }

    }
}
