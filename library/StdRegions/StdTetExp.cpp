///////////////////////////////////////////////////////////////////////////////
//
// File StdTetExp.h
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
// Description: Header field for tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <StdRegions/StdTetExp.h>

using namespace std;

namespace Nektar
{
    namespace StdRegions
    {
        StdTetExp::StdTetExp()
        {
        }


        StdTetExp::StdTetExp(const LibUtilities::BasisKey &Ba,
                             const LibUtilities::BasisKey &Bb,
                             const LibUtilities::BasisKey &Bc):
            StdExpansion(LibUtilities::StdTetData::getNumberOfCoefficients(
                             Ba.GetNumModes(),
                             Bb.GetNumModes(),
                             Bc.GetNumModes()),
                         3, Ba, Bb, Bc),
            StdExpansion3D(LibUtilities::StdTetData::getNumberOfCoefficients(
                               Ba.GetNumModes(),
                               Bb.GetNumModes(),
                               Bc.GetNumModes()),
                           Ba, Bb, Bc)
        {
            ASSERTL0(Ba.GetNumModes() <= Bb.GetNumModes(),
                     "order in 'a' direction is higher than order "
                     "in 'b' direction");
            ASSERTL0(Ba.GetNumModes() <= Bc.GetNumModes(),
                     "order in 'a' direction is higher than order "
                     "in 'c' direction");
            ASSERTL0(Bb.GetNumModes() <= Bc.GetNumModes(),
                     "order in 'b' direction is higher than order "
                     "in 'c' direction");
        }

        StdTetExp::StdTetExp(const StdTetExp &T):
            StdExpansion(T),
            StdExpansion3D(T)
        {
        }


        StdTetExp::~StdTetExp()
        {
        }

        //----------------------------
        // Differentiation Methods
        //----------------------------

        /**
         * \brief Calculate the derivative of the physical points
         *
         * The derivative is evaluated at the nodal physical points.
         * Derivatives with respect to the local Cartesian coordinates
         *
         * \f$\begin{Bmatrix} \frac {\partial} {\partial \xi_1} \\ \frac
         * {\partial} {\partial \xi_2} \\ \frac {\partial} {\partial \xi_3}
         * \end{Bmatrix} = \begin{Bmatrix} \frac 4 {(1-\eta_2)(1-\eta_3)}
         * \frac \partial {\partial \eta_1} \ \ \frac {2(1+\eta_1)}
         * {(1-\eta_2)(1-\eta_3)} \frac \partial {\partial \eta_1} + \frac 2
         * {1-\eta_3} \frac \partial {\partial \eta_3} \\ \frac {2(1 +
         * \eta_1)} {2(1 - \eta_2)(1-\eta_3)} \frac \partial {\partial \eta_1}
         * + \frac {1 + \eta_2} {1 - \eta_3} \frac \partial {\partial \eta_2}
         * + \frac \partial {\partial \eta_3} \end{Bmatrix}\f$
         **/
        void StdTetExp::v_PhysDeriv(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& out_dxi0,
                  Array<OneD,       NekDouble>& out_dxi1,
                  Array<OneD,       NekDouble>& out_dxi2 )
        {
            int    Q0 = m_base[0]->GetNumPoints();
            int    Q1 = m_base[1]->GetNumPoints();
            int    Q2 = m_base[2]->GetNumPoints();
            int    Qtot = Q0*Q1*Q2;

            // Compute the physical derivative
            Array<OneD, NekDouble> out_dEta0(3*Qtot,0.0);
            Array<OneD, NekDouble> out_dEta1 = out_dEta0 + Qtot;
            Array<OneD, NekDouble> out_dEta2 = out_dEta1 + Qtot;

            bool Do_2 = (out_dxi2.size() > 0)? true:false;
            bool Do_1 = (out_dxi1.size() > 0)? true:false;

            if(Do_2) // Need all local derivatives
            {
                PhysTensorDeriv(inarray, out_dEta0, out_dEta1, out_dEta2);
            }
            else if (Do_1) // Need 0 and 1 derivatives
            {
                PhysTensorDeriv(inarray, out_dEta0, out_dEta1, NullNekDouble1DArray);
            }
            else // Only need Eta0 derivaitve
            {
                PhysTensorDeriv(inarray, out_dEta0, NullNekDouble1DArray,
                                NullNekDouble1DArray);
            }

            Array<OneD, const NekDouble> eta_0, eta_1, eta_2;
            eta_0 = m_base[0]->GetZ();
            eta_1 = m_base[1]->GetZ();
            eta_2 = m_base[2]->GetZ();

            // calculate 2.0/((1-eta_1)(1-eta_2)) Out_dEta0

            NekDouble *dEta0 = &out_dEta0[0];
            NekDouble fac;
            for(int k=0; k< Q2; ++k)
            {
                for(int j=0; j<Q1; ++j,dEta0+=Q0)
                {
                    Vmath::Smul(Q0,2.0/(1.0-eta_1[j]),dEta0,1,dEta0,1);
                }
                fac = 1.0/(1.0-eta_2[k]);
                Vmath::Smul(Q0*Q1,fac,&out_dEta0[0]+k*Q0*Q1,1,&out_dEta0[0]+k*Q0*Q1,1);
            }

            if (out_dxi0.size() > 0)
            {
                // out_dxi0 = 4.0/((1-eta_1)(1-eta_2)) Out_dEta0
                Vmath::Smul(Qtot,2.0,out_dEta0,1,out_dxi0,1);
            }

            if (Do_1||Do_2)
            {
                Array<OneD, NekDouble> Fac0(Q0);
                Vmath::Sadd(Q0,1.0,eta_0,1,Fac0,1);


                // calculate 2.0*(1+eta_0)/((1-eta_1)(1-eta_2)) Out_dEta0
                for(int k = 0; k < Q1*Q2; ++k)
                {
                    Vmath::Vmul(Q0,&Fac0[0],1,&out_dEta0[0]+k*Q0,1,&out_dEta0[0]+k*Q0,1);
                }
                // calculate 2/(1.0-eta_2) out_dEta1
                for(int k = 0; k < Q2; ++k)
                {
                    Vmath::Smul(Q0*Q1,2.0/(1.0-eta_2[k]),&out_dEta1[0]+k*Q0*Q1,1,
                                &out_dEta1[0]+k*Q0*Q1,1);
                }

                if(Do_1)
                {
                    // calculate out_dxi1 = 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) Out_dEta0
                    // + 2/(1.0-eta_2) out_dEta1
                    Vmath::Vadd(Qtot,out_dEta0,1,out_dEta1,1,out_dxi1,1);
                }


                if(Do_2)
                {
                    // calculate (1 + eta_1)/(1 -eta_2)*out_dEta1
                    NekDouble *dEta1 = &out_dEta1[0];
                    for(int k=0; k< Q2; ++k)
                    {
                        for(int j=0; j<Q1; ++j,dEta1+=Q0)
                        {
                            Vmath::Smul(Q0,(1.0+eta_1[j])/2.0,dEta1,1,dEta1,1);
                        }
                    }

                    // calculate out_dxi2 =
                    // 2.0(1+eta_0)/((1-eta_1)(1-eta_2)) Out_dEta0 +
                    // (1 + eta_1)/(1 -eta_2)*out_dEta1 + out_dEta2
                    Vmath::Vadd(Qtot,out_dEta0,1,out_dEta1,1,out_dxi2,1);
                    Vmath::Vadd(Qtot,out_dEta2,1,out_dxi2 ,1,out_dxi2,1);

                }
            }
        }

        /**
         * @param   dir         Direction in which to compute derivative.
         *                      Valid values are 0, 1, 2.
         * @param   inarray     Input array.
         * @param   outarray    Output array.
         */
        void StdTetExp::v_PhysDeriv(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            switch(dir)
            {
                case 0:
                {
                    v_PhysDeriv(inarray, outarray, NullNekDouble1DArray,
                                NullNekDouble1DArray);
                    break;
                }
                case 1:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray, outarray,
                                NullNekDouble1DArray);
                    break;
                }
                case 2:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray,
                                NullNekDouble1DArray, outarray);
                    break;
                }
            default:
                {
                    ASSERTL1(false, "input dir is out of range");
                }
                break;
            }
        }

        void StdTetExp::v_StdPhysDeriv(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& out_d0,
                  Array<OneD,       NekDouble>& out_d1,
                  Array<OneD,       NekDouble>& out_d2)
        {
            StdTetExp::v_PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }

        void StdTetExp::v_StdPhysDeriv(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTetExp::v_PhysDeriv(dir, inarray, outarray);
        }

        //---------------------------------------
        // Transforms
        //---------------------------------------

        /**
         * @note 'r' (base[2]) runs fastest in this element
         *
         * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat
         *  u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
         *
         * Backward transformation is three dimensional tensorial expansion
         * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a
         * (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{pq}^b (\xi_{2j})
         * \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c (\xi_{3k})
         * \rbrace} \rbrace}. \f$ And sumfactorizing step of the form is as:\\
         *
         * \f$ f_{pq} (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c
         * (\xi_{3k}),\\
         *
         * g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y} \psi_{pq}^b
         * (\xi_{2j}) f_{pq} (\xi_{3k}),\\
         *
         * u(\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a
         * (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).  \f$
         */
        void StdTetExp::v_BwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,       NekDouble>& outarray)
        {
            ASSERTL1((m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)  ||
                     (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                     "Basis[1] is not a general tensor type");

            ASSERTL1((m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
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
                StdTetExp::v_BwdTrans_SumFac(inarray,outarray);
            }
        }


        /**
         * Sum-factorisation implementation of the BwdTrans operation.
         */
        void StdTetExp::v_BwdTrans_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();

            Array<OneD, NekDouble> wsp(nquad2*order0*(2*order1-order0+1)/2+
                                       nquad2*nquad1*order0);

            BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                  m_base[1]->GetBdata(),
                                  m_base[2]->GetBdata(),
                                  inarray,outarray,wsp,
                                  true,true,true);
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
        void StdTetExp::v_BwdTrans_SumFacKernel(
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
            boost::ignore_unused(doCheckCollDir0, doCheckCollDir1,
                                 doCheckCollDir2);

            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();

            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble > tmp  = wsp;
            Array<OneD, NekDouble > tmp1 = tmp + nquad2*order0*(2*order1-order0+1)/2;

            int i, j, mode, mode1, cnt;

            // Perform summation over '2' direction
            mode = mode1 = cnt = 0;
            for(i = 0; i < order0; ++i)
            {
                for(j = 0; j < order1-i; ++j, ++cnt)
                {
                    Blas::Dgemv('N', nquad2, order2-i-j,
                                1.0, base2.get()+mode*nquad2, nquad2,
                                     inarray.get()+mode1,     1,
                                0.0, tmp.get()+cnt*nquad2,    1);
                    mode  += order2-i-j;
                    mode1 += order2-i-j;
                }
                //increment mode in case order1!=order2
                for(j = order1-i; j < order2-i; ++j)
                {
                    mode += order2-i-j;
                }
            }

            // fix for modified basis by adding split of top singular
            // vertex mode - currently (1+c)/2 x (1-b)/2 x (1-a)/2
            // component is evaluated
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                // top singular vertex - (1+c)/2 x (1+b)/2 x (1-a)/2 component
                Blas::Daxpy(nquad2,inarray[1],base2.get()+nquad2,1,
                            &tmp[0]+nquad2,1);

                // top singular vertex - (1+c)/2 x (1-b)/2 x (1+a)/2 component
                Blas::Daxpy(nquad2,inarray[1],base2.get()+nquad2,1,
                            &tmp[0]+order1*nquad2,1);
            }

            // Perform summation over '1' direction
            mode = 0;
            for(i = 0; i < order0; ++i)
            {
                Blas::Dgemm('N', 'T', nquad1, nquad2, order1-i,
                            1.0, base1.get()+mode*nquad1,    nquad1,
                                 tmp.get()+mode*nquad2,      nquad2,
                            0.0, tmp1.get()+i*nquad1*nquad2, nquad1);
                mode  += order1-i;
            }

            // fix for modified basis by adding additional split of
            // top and base singular vertex modes as well as singular
            // edge
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                // use tmp to sort out singular vertices and
                // singular edge components with (1+b)/2 (1+a)/2 form
                for(i = 0; i < nquad2; ++i)
                {
                    Blas::Daxpy(nquad1,tmp[nquad2+i], base1.get()+nquad1,1,
                                &tmp1[nquad1*nquad2]+i*nquad1,1);
                }
            }

            // Perform summation over '0' direction
            Blas::Dgemm('N', 'T', nquad0, nquad1*nquad2, order0,
                        1.0, base0.get(),    nquad0,
                             tmp1.get(),     nquad1*nquad2,
                        0.0, outarray.get(), nquad0);
        }


        /**
         * @param   inarray     array of physical quadrature points to be
         *                      transformed.
         * @param   outarray    updated array of expansion coefficients.
         */
        void StdTetExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
        {            //int       numMax  = nmodes0;
            v_IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs, outarray);
            DNekVec out(m_ncoeffs, outarray, eWrapper);

            out = (*matsys)*in;
        }


        //---------------------------------------
        // Inner product functions
        //---------------------------------------

        /**
         * \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
         * \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2} \psi_{p}^{a}
         * (\eta_{1i}) \psi_{pq}^{b} (\eta_{2j}) \psi_{pqr}^{c} (\eta_{3k})
         * w_i w_j w_k u(\eta_{1,i} \eta_{2,j} \eta_{3,k}) J_{i,j,k}\\ & = &
         * \sum_{i=0}^{nq_0} \psi_p^a(\eta_{1,i}) \sum_{j=0}^{nq_1}
         * \psi_{pq}^b(\eta_{2,j}) \sum_{k=0}^{nq_2} \psi_{pqr}^c
         * u(\eta_{1i},\eta_{2j},\eta_{3k}) J_{i,j,k} \end{array} \f$ \n
         *
         * where
         *
         * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\eta_1)
         * \psi_{pq}^b (\eta_2) \psi_{pqr}^c (\eta_3) \f$
         *
         * which can be implemented as \n \f$f_{pqr} (\xi_{3k}) =
         * \sum_{k=0}^{nq_3} \psi_{pqr}^c u(\eta_{1i},\eta_{2j},\eta_{3k})
         *
         * J_{i,j,k} = {\bf B_3 U}   \f$ \n
         *
         * \f$ g_{pq} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{pq}^b (\xi_{2j})
         * f_{pqr} (\xi_{3k}) = {\bf B_2 F} \f$ \n
         *
         * \f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0} \psi_{p}^a
         * (\xi_{3k}) g_{pq} (\xi_{3k}) = {\bf B_1 G} \f$
         *
         * @param   inarray     Function evaluated at physical collocation
         *                      points.
         * @param   outarray    Inner product with respect to each basis
         *                      function over the element.
         */
        void StdTetExp::v_IProductWRTBase(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> & outarray)
        {
            ASSERTL1( (m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)  ||
                      (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                      "Basis[1] is not a general tensor type");

            ASSERTL1( (m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
                      (m_base[2]->GetBasisType() != LibUtilities::eModified_C),
                      "Basis[2] is not a general tensor type");

            if(m_base[0]->Collocation() && m_base[1]->Collocation())
            {
                MultiplyByQuadratureMetric(inarray,outarray);
            }
            else
            {
                StdTetExp::v_IProductWRTBase_SumFac(inarray,outarray);
            }
        }


        void StdTetExp::v_IProductWRTBase_MatOp(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }


        /**
         * @param   inarray     Function evaluated at physical collocation
         *                      points.
         * @param   outarray    Inner product with respect to each basis
         *                      function over the element.
         */
        void StdTetExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray,
            bool                                multiplybyweights)
        {
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();

            Array<OneD, NekDouble> wsp (nquad1*nquad2*order0 +
                                        nquad2*order0*(2*order1-order0+1)/2);

            if(multiplybyweights)
            {
                Array<OneD, NekDouble> tmp (nquad0*nquad1*nquad2);
                MultiplyByQuadratureMetric(inarray, tmp);

                StdTetExp::IProductWRTBase_SumFacKernel(
                              m_base[0]->GetBdata(),
                              m_base[1]->GetBdata(),
                              m_base[2]->GetBdata(),
                              tmp, outarray, wsp, true, true, true);
            }
            else
            {
                StdTetExp::IProductWRTBase_SumFacKernel(
                               m_base[0]->GetBdata(),
                               m_base[1]->GetBdata(),
                               m_base[2]->GetBdata(),
                               inarray, outarray, wsp, true, true, true);
            }
        }


        void StdTetExp::v_IProductWRTBase_SumFacKernel(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& base2,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,       NekDouble> &outarray,
                          Array<OneD,       NekDouble> &wsp,
                          bool                          doCheckCollDir0,
                          bool                          doCheckCollDir1,
                          bool                          doCheckCollDir2)
        {
            boost::ignore_unused(doCheckCollDir0, doCheckCollDir1,
                                 doCheckCollDir2);

            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();

            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble > tmp1 = wsp;
            Array<OneD, NekDouble > tmp2 = wsp + nquad1*nquad2*order0;

            int i,j, mode,mode1, cnt;

            // Inner product with respect to the '0' direction
            Blas::Dgemm('T', 'N', nquad1*nquad2, order0, nquad0,
                        1.0, inarray.get(), nquad0,
                             base0.get(),   nquad0,
                        0.0, tmp1.get(),    nquad1*nquad2);

            // Inner product with respect to the '1' direction
            for(mode=i=0; i < order0; ++i)
            {
                Blas::Dgemm('T', 'N', nquad2, order1-i, nquad1,
                            1.0, tmp1.get()+i*nquad1*nquad2, nquad1,
                                 base1.get()+mode*nquad1,    nquad1,
                            0.0, tmp2.get()+mode*nquad2,     nquad2);
                mode  += order1-i;
            }

            // fix for modified basis for base singular vertex
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                //base singular vertex and singular edge (1+b)/2
                //(1+a)/2 components (makes tmp[nquad2] entry into (1+b)/2)
                Blas::Dgemv('T', nquad1, nquad2,
                            1.0, tmp1.get()+nquad1*nquad2, nquad1,
                                 base1.get()+nquad1,       1,
                            1.0, tmp2.get()+nquad2,        1);
            }

            // Inner product with respect to the '2' direction
            mode = mode1 = cnt = 0;
            for(i = 0; i < order0; ++i)
            {
                for(j = 0; j < order1-i; ++j, ++cnt)
                {
                    Blas::Dgemv('T', nquad2, order2-i-j,
                                1.0, base2.get()+mode*nquad2, nquad2,
                                     tmp2.get()+cnt*nquad2,   1,
                                0.0, outarray.get()+mode1,    1);
                    mode  += order2-i-j;
                    mode1 += order2-i-j;
                }
                //increment mode in case order1!=order2
                for(j = order1-i; j < order2-i; ++j)
                {
                    mode += order2-i-j;
                }
            }

            // fix for modified basis for top singular vertex component
            // Already have evaluated (1+c)/2 (1-b)/2 (1-a)/2
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                // add in (1+c)/2 (1+b)/2   component
                outarray[1] += Blas::Ddot(nquad2,base2.get()+nquad2,1,
                                          &tmp2[nquad2],1);

                // add in (1+c)/2 (1-b)/2 (1+a)/2 component
                outarray[1] += Blas::Ddot(nquad2,base2.get()+nquad2,1,
                                          &tmp2[nquad2*order1],1);
            }
        }


        void StdTetExp::v_IProductWRTDerivBase(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTetExp::v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }


        void StdTetExp::v_IProductWRTDerivBase_MatOp(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
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

            StdMatrixKey     iprodmatkey(mtype,DetShapeType(),*this);
            DNekMatSharedPtr iprodmat = GetStdMatrix(iprodmatkey);

            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }


        /**
         * @param   inarray     Function evaluated at physical collocation
         *                      points.
         * @param   outarray    Inner product with respect to each basis
         *                      function over the element.
         */
        void StdTetExp::v_IProductWRTDerivBase_SumFac(
            const int dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int    i;
            int    nquad0  = m_base[0]->GetNumPoints();
            int    nquad1  = m_base[1]->GetNumPoints();
            int    nquad2  = m_base[2]->GetNumPoints();
            int    nqtot   = nquad0*nquad1*nquad2;
            int    nmodes0 = m_base[0]->GetNumModes();
            int    nmodes1 = m_base[1]->GetNumModes();
            int    wspsize = nquad0 + nquad1 + nquad2 + max(nqtot,m_ncoeffs)
                + nquad1*nquad2*nmodes0 + nquad2*nmodes0*(2*nmodes1-nmodes0+1)/2;

            Array<OneD, NekDouble> gfac0(wspsize);
            Array<OneD, NekDouble> gfac1(gfac0 + nquad0);
            Array<OneD, NekDouble> gfac2(gfac1 + nquad1);
            Array<OneD, NekDouble> tmp0 (gfac2 + nquad2);
            Array<OneD, NekDouble> wsp(tmp0 + max(nqtot,m_ncoeffs));

            const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();
            const Array<OneD, const NekDouble>& z2 = m_base[2]->GetZ();

            // set up geometric factor: (1+z0)/2
            for(i = 0; i < nquad0; ++i)
            {
                gfac0[i] = 0.5*(1+z0[i]);
            }

            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac1[i] = 2.0/(1-z1[i]);
            }

            // Set up geometric factor: 2/(1-z2)
            for(i = 0; i < nquad2; ++i)
            {
            	gfac2[i] = 2.0/(1-z2[i]);
            }

            // Derivative in first direction is always scaled as follows
            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Smul(nquad0,gfac1[i%nquad1],&inarray[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
            }
            for(i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nquad0*nquad1,gfac2[i],&tmp0[0]+i*nquad0*nquad1,1,&tmp0[0]+i*nquad0*nquad1,1);
            }

            MultiplyByQuadratureMetric(tmp0,tmp0);

            switch(dir)
            {
            case 0:
                {
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp0,outarray,wsp,
                                                 false, true, true);
                }
                break;
            case 1:
                {
                    Array<OneD, NekDouble> tmp3(m_ncoeffs);

                    for(i = 0; i < nquad1*nquad2; ++i)
                    {
                        Vmath::Vmul(nquad0,&gfac0[0],1,&tmp0[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
                    }

                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp0,tmp3,wsp,
                                                 false, true, true);

                    for(i = 0; i < nquad2; ++i)
                    {
                        Vmath::Smul(nquad0*nquad1,gfac2[i],&inarray[0]+i*nquad0*nquad1,1,&tmp0[0]+i*nquad0*nquad1,1);
                    }
                    MultiplyByQuadratureMetric(tmp0,tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetDbdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp0,outarray,wsp,
                                                 true, false, true);
                    Vmath::Vadd(m_ncoeffs,&tmp3[0],1,&outarray[0],1,&outarray[0],1);
                }
                break;
            case 2:
				{
                    Array<OneD, NekDouble> tmp3(m_ncoeffs);
                    Array<OneD, NekDouble> tmp4(m_ncoeffs);
                    for(i = 0; i < nquad1; ++i)
                    {
                        gfac1[i] = (1+z1[i])/2;
                    }

                    for(i = 0; i < nquad1*nquad2; ++i)
                    {
                        Vmath::Vmul(nquad0,&gfac0[0],1,&tmp0[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
                    }
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp0,tmp3,wsp,
                                                 false, true, true);

                    for(i = 0; i < nquad2; ++i)
                    {
                        Vmath::Smul(nquad0*nquad1,gfac2[i],&inarray[0]+i*nquad0*nquad1,1,&tmp0[0]+i*nquad0*nquad1,1);
                    }
                    for(i = 0; i < nquad1*nquad2; ++i)
                    {
                        Vmath::Smul(nquad0,gfac1[i%nquad1],&tmp0[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
                    }
                    MultiplyByQuadratureMetric(tmp0,tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetDbdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp0,tmp4,wsp,
                                                 true, false, true);

                    MultiplyByQuadratureMetric(inarray,tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetDbdata(),
                                                 tmp0,outarray,wsp,
                                                 true, true, false);

                    Vmath::Vadd(m_ncoeffs,&tmp3[0],1,&outarray[0],1,&outarray[0],1);
                    Vmath::Vadd(m_ncoeffs,&tmp4[0],1,&outarray[0],1,&outarray[0],1);
				}
                break;
            default:
                {
                    ASSERTL1(false, "input dir is out of range");
                }
                break;
            }
        }


        //---------------------------------------
        // Evaluation functions
        //---------------------------------------


        void StdTetExp::v_LocCoordToLocCollapsed(
                                        const Array<OneD, const NekDouble>& xi,
                                        Array<OneD, NekDouble>& eta)
        {
            if( fabs(xi[2]-1.0) < NekConstants::kNekZeroTol)
            {
                // Very top point of the tetrahedron
                eta[0] = -1.0;
                eta[1] = -1.0;
                eta[2] = xi[2];
            }
            else
            {
                if( fabs(xi[1]-1.0) <  NekConstants::kNekZeroTol )
                {
                    // Distant diagonal edge shared by all eta_x
                    // coordinate planes: the xi_y == -xi_z line
                    eta[0] = -1.0;
                }
                else if (fabs(xi[1] + xi[2]) < NekConstants::kNekZeroTol)
                {
                    eta[0] = -1.0;
                }
                else
                {
                    eta[0] = 2.0*(1.0+xi[0])/(-xi[1]-xi[2]) - 1.0;
                }
                eta[1] = 2.0*(1.0+xi[1])/(1.0-xi[2]) - 1.0;
                eta[2] = xi[2];
            }
        }

        void StdTetExp::v_FillMode(
            const int                     mode,
                  Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp(m_ncoeffs,0.0);
            tmp[mode] = 1.0;
            StdTetExp::v_BwdTrans(tmp, outarray);
        }

        void StdTetExp::v_GetFaceNumModes(
                    const int                  fid,
                    const Orientation          faceOrient,
                    int &numModes0,
                    int &numModes1)
        {
            boost::ignore_unused(faceOrient);

            int nummodes [3] = {m_base[0]->GetNumModes(),
                                m_base[1]->GetNumModes(),
                                m_base[2]->GetNumModes()};
            switch(fid)
            {
            case 0:
                {
                    numModes0 = nummodes[0];
                    numModes1 = nummodes[1];
                }
                break;
            case 1:
                {
                    numModes0 = nummodes[0];
                    numModes1 = nummodes[2];
                }
                break;
            case 2:
            case 3:
                {
                    numModes0 = nummodes[1];
                    numModes1 = nummodes[2];
                }
                break;
            }
        }

        //---------------------------
        // Helper functions
        //---------------------------

        int StdTetExp::v_GetNverts() const
        {
            return 4;
        }

        int StdTetExp::v_GetNedges() const
        {
            return 6;
        }

        int StdTetExp::v_GetNfaces() const
        {
            return 4;
        }

        LibUtilities::ShapeType StdTetExp::v_DetShapeType() const
        {
            return DetShapeType();
        }

        int StdTetExp::v_NumBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int P = m_base[0]->GetNumModes();
            int Q = m_base[1]->GetNumModes();
            int R = m_base[2]->GetNumModes();

            return LibUtilities::StdTetData::
                                        getNumberOfBndCoefficients(P, Q, R);
        }

        int StdTetExp::v_NumDGBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int P = m_base[0]->GetNumModes()-1;
            int Q = m_base[1]->GetNumModes()-1;
            int R = m_base[2]->GetNumModes()-1;


            return  (Q+1) + P*(1 + 2*Q - P)/2  // base face
                +   (R+1) + P*(1 + 2*R - P)/2  // front face
                + 2*(R+1) + Q*(1 + 2*R - Q);   // back two faces
        }

        int StdTetExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2((i >= 0) && (i <= 5), "edge id is out of range");
            int P = m_base[0]->GetNumModes();
            int Q = m_base[1]->GetNumModes();
            int R = m_base[2]->GetNumModes();

            if (i == 0)
            {
                return P;
            }
            else if (i == 1 || i == 2)
            {
                return Q;
            }
            else
            {
                return R;
            }
        }

        int StdTetExp::v_GetTotalEdgeIntNcoeffs() const
        {
            int P = m_base[0]->GetNumModes()-2;
            int Q = m_base[1]->GetNumModes()-2;
            int R = m_base[2]->GetNumModes()-2;

            return P+Q+4*R;
	}

        int StdTetExp::v_GetFaceNcoeffs(const int i) const
        {
            ASSERTL2((i >= 0) && (i <= 3), "face id is out of range");
            int nFaceCoeffs = 0;
            int nummodesA, nummodesB, P, Q;
            if (i == 0)
            {
                nummodesA = GetBasisNumModes(0);
                nummodesB = GetBasisNumModes(1);
            }
            else if ((i == 1) || (i == 2))
            {
                nummodesA = GetBasisNumModes(0);
                nummodesB = GetBasisNumModes(2);
            }
            else
            {
                nummodesA = GetBasisNumModes(1);
                nummodesB = GetBasisNumModes(2);
            }
            P = nummodesA - 1;
            Q = nummodesB - 1;
            nFaceCoeffs = Q+1 + (P*(1 + 2*Q - P))/2;
            return nFaceCoeffs;
        }

        int StdTetExp::v_GetFaceIntNcoeffs(const int i) const
        {
            ASSERTL2((i >= 0) && (i <= 3), "face id is out of range");
            int Pi = m_base[0]->GetNumModes() - 2;
            int Qi = m_base[1]->GetNumModes() - 2;
            int Ri = m_base[2]->GetNumModes() - 2;

            if((i == 0))
            {
                return Pi * (2*Qi - Pi - 1) / 2;
            }
            else if((i == 1))
            {
                return Pi * (2*Ri - Pi - 1) / 2;
            }
            else
            {
                return Qi * (2*Ri - Qi - 1) / 2;
            }
        }

        int StdTetExp::v_GetTotalFaceIntNcoeffs() const
        {
            int Pi = m_base[0]->GetNumModes() - 2;
            int Qi = m_base[1]->GetNumModes() - 2;
            int Ri = m_base[2]->GetNumModes() - 2;

            return Pi * (2*Qi - Pi - 1) / 2 +
	           Pi * (2*Ri - Pi - 1) / 2 +
	           Qi * (2*Ri - Qi - 1);
	}

        int StdTetExp::v_GetFaceNumPoints(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 3, "face id is out of range");

            if (i == 0)
            {
                return m_base[0]->GetNumPoints()*
                       m_base[1]->GetNumPoints();
            }
            else if (i == 1)
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

        LibUtilities::PointsKey StdTetExp::v_GetFacePointsKey(
            const int i, const int j) const
        {
            ASSERTL2(i >= 0 && i <= 3, "face id is out of range");
            ASSERTL2(j == 0 || j == 1, "face direction is out of range");

            if (i == 0)
            {
                return m_base[j]->GetPointsKey();
            }
            else if (i == 1)
            {
                return m_base[2*j]->GetPointsKey();
            }
            else
            {
                return m_base[j+1]->GetPointsKey();
            }
        }

        int StdTetExp::v_CalcNumberOfCoefficients(
            const std::vector<unsigned int>& nummodes,
                  int                      & modes_offset)
        {
            int nmodes = LibUtilities::StdTetData::getNumberOfCoefficients(
                nummodes[modes_offset],
                nummodes[modes_offset+1],
                nummodes[modes_offset+2]);
            modes_offset += 3;

            return nmodes;
        }

        const LibUtilities::BasisKey StdTetExp::v_DetFaceBasisKey(
            const int i, const int k) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            ASSERTL2(k == 0 || k == 1, "face direction out of range");

            int dir = k;
            switch(i)
            {
                case 0:
                    dir = k;
                    break;
                case 1:
                    dir = 2*k;
                    break;
                case 2:
                case 3:
                    dir = k+1;
                    break;
            }

            return EvaluateTriFaceBasisKey(k,
                                           m_base[dir]->GetBasisType(),
                                           m_base[dir]->GetNumPoints(),
                                           m_base[dir]->GetNumModes());

            // Should not get here.
            return LibUtilities::NullBasisKey;
        }

        LibUtilities::BasisType StdTetExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 5, "edge id is out of range");

            if (i == 0)
            {
                return GetBasisType(0);
            }
            else if (i == 1 || i == 2)
            {
                return GetBasisType(1);
            }
            else
            {
                return GetBasisType(2);
            }
        }

        void StdTetExp::v_GetCoords(
            Array<OneD, NekDouble> &xi_x,
            Array<OneD, NekDouble> &xi_y,
            Array<OneD, NekDouble> &xi_z)
        {
            Array<OneD, const NekDouble> eta_x = m_base[0]->GetZ();
            Array<OneD, const NekDouble> eta_y = m_base[1]->GetZ();
            Array<OneD, const NekDouble> eta_z = m_base[2]->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates: eta
            // --> xi
            for( int k = 0; k < Qz; ++k ) {
                for( int j = 0; j < Qy; ++j ) {
                    for( int i = 0; i < Qx; ++i ) {
                        int s = i + Qx*(j + Qy*k);
                        xi_x[s] = (eta_x[i] + 1.0) * (1.0 - eta_y[j]) * (1.0 - eta_z[k]) / 4  -  1.0;
                        xi_y[s] = (eta_y[j] + 1.0) * (1.0 - eta_z[k]) / 2  -  1.0;
                        xi_z[s] = eta_z[k];
                    }
                }
            }
        }

        bool StdTetExp::v_IsBoundaryInteriorExpansion()
        {
            return (m_base[0]->GetBasisType() == LibUtilities::eModified_A) &&
                   (m_base[1]->GetBasisType() == LibUtilities::eModified_B) &&
                   (m_base[2]->GetBasisType() == LibUtilities::eModified_C);
        }


        //--------------------------
        // Mappings
        //--------------------------

        void StdTetExp::v_GetEdgeToElementMap(
            const int                  eid,
            const Orientation          edgeOrient,
            Array<OneD, unsigned int>& maparray,
            Array<OneD,          int>& signarray,
            int                        P)
        {
            boost::ignore_unused(P);

            ASSERTL2(eid >= 0 && eid < 6, "Invalid edge");
            ASSERTL2(v_IsBoundaryInteriorExpansion(),
                     "Method only implemented for Modified_A BasisType (x "
                     "direction), Modified_B BasisType (y direction), and "
                     "Modified_C BasisType(z direction)");

            int edgeVerts[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
            int nEdgeModes = v_GetEdgeNcoeffs(eid);

            if (maparray.size() != nEdgeModes)
            {
                maparray  = Array<OneD, unsigned int>(nEdgeModes);
            }

            if (signarray.size() != nEdgeModes)
            {
                signarray = Array<OneD, int>(nEdgeModes, 1.0);
            }

            if (edgeOrient == StdRegions::eForwards)
            {
                maparray[0] = v_GetVertexMap(edgeVerts[eid][0]);
                maparray[1] = v_GetVertexMap(edgeVerts[eid][1]);
            }
            else if (edgeOrient == StdRegions::eBackwards)
            {
                maparray[0] = v_GetVertexMap(edgeVerts[eid][1]);
                maparray[1] = v_GetVertexMap(edgeVerts[eid][0]);
            }
            else
            {
                ASSERTL2(false, "Unsupported edge orientation");
            }

            Array<OneD, unsigned int> tmp1(nEdgeModes-2);
            Array<OneD,          int> tmp2(nEdgeModes-2);
            v_GetEdgeInteriorMap(eid, edgeOrient, tmp1, tmp2);

            for (int i = 0; i < nEdgeModes-2; ++i)
            {
                maparray[i+2] = tmp1[i];
                signarray[i+2] = tmp2[i];
            }
        }

       /**
         * Maps Expansion2D modes of a 2D face to the corresponding expansion
         * modes.
         */
        void StdTetExp::v_GetFaceToElementMap(
            const int                  fid,
            const Orientation          faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray,
            int                        P,
            int                        Q)
        {
            int nummodesA=0,nummodesB=0, i, j, k, idx;

            ASSERTL1(v_IsBoundaryInteriorExpansion(),
                     "Method only implemented for Modified_A BasisType (x "
                     "direction), Modified_B BasisType (y direction), and "
                     "Modified_C BasisType(z direction)");

            int nFaceCoeffs = 0;

            switch(fid)
            {
            case 0:
                nummodesA = m_base[0]->GetNumModes();
                nummodesB = m_base[1]->GetNumModes();
                break;
            case 1:
                nummodesA = m_base[0]->GetNumModes();
                nummodesB = m_base[2]->GetNumModes();
                        break;
            case 2:
            case 3:
                nummodesA = m_base[1]->GetNumModes();
                nummodesB = m_base[2]->GetNumModes();
                break;
            default:
                ASSERTL0(false,"fid must be between 0 and 3");
            }

            bool CheckForZeroedModes = false;
            if(P == -1)
            {
                P = nummodesA;
                Q = nummodesB;
            }
            else
            {
                CheckForZeroedModes = true;
            }

            nFaceCoeffs = P*(2*Q-P+1)/2;

            // Allocate the map array and sign array; set sign array to ones (+)
            if(maparray.size() != nFaceCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceCoeffs,1);
            }

            if(signarray.size() != nFaceCoeffs)
            {
                signarray = Array<OneD, int>(nFaceCoeffs,1);
            }
            else
            {
                fill(signarray.get(),signarray.get()+nFaceCoeffs, 1 );
            }

            switch (fid)
            {
            case 0:
                idx = 0;
                for (i = 0; i < P; ++i)
                {
                    for (j = 0; j < Q-i; ++j)
                    {
                        if ((int)faceOrient == 7 && i > 1)
                        {
                            signarray[idx] = (i%2 ? -1 : 1);
                        }
                        maparray[idx++] = GetMode(i,j,0);
                    }
                }
                break;
            case 1:
                idx = 0;
                for (i = 0; i < P; ++i)
                {
                    for (k = 0; k < Q-i; ++k)
                    {
                        if ((int)faceOrient == 7 && i > 1)
                        {
                            signarray[idx] = (i%2 ? -1: 1);
                        }
                        maparray[idx++] = GetMode(i,0,k);
                    }
                }
                break;
            case 2:
                idx = 0;
                for (j = 0; j < P-1; ++j)
                {
                    for (k = 0; k < Q-1-j; ++k)
                    {
                        if ((int)faceOrient == 7 && j > 1)
                        {
                            signarray[idx] = ((j+1)%2 ? -1: 1);
                        }
                        maparray[idx++] = GetMode(1,j,k);
                        // Incorporate modes from zeroth plane where needed.
                        if (j == 0 && k == 0)
                        {
                            maparray[idx++] = GetMode(0,0,1);
                        }
                        if (j == 0 && k == Q-2)
                        {
                            for (int r = 0; r < Q-1; ++r)
                            {
                                maparray[idx++] = GetMode(0,1,r);
                            }
                        }
                    }
                }
                break;
            case 3:
                idx = 0;
                for (j = 0; j < P; ++j)
                {
                    for (k = 0; k < Q-j; ++k)
                    {
                        if ((int)faceOrient == 7 && j > 1)
                        {
                            signarray[idx] = (j%2 ? -1: 1);
                        }
                        maparray[idx++] = GetMode(0,j,k);
                    }
                }
                break;
            default:
                ASSERTL0(false, "Element map not available.");
            }

            if ((int)faceOrient == 7)
            {
                swap(maparray[0], maparray[Q]);

                for (i = 1; i < Q-1; ++i)
                {
                    swap(maparray[i+1], maparray[Q+i]);
                }
            }

            if(CheckForZeroedModes)
            {
                // zero signmap and set maparray to zero if elemental
                // modes are not as large as face modesl
                idx = 0;
                for (j = 0; j < nummodesA; ++j)
                {
                    idx += nummodesB-j;
                    for (k = nummodesB-j; k < Q-j; ++k)
                    {
                        signarray[idx]  = 0.0;
                        maparray[idx++] = maparray[0];
                    }
                }

                for (j = nummodesA; j < P; ++j)
                {
                    for (k = 0; k < Q-j; ++k)
                    {
                        signarray[idx]  = 0.0;
                        maparray[idx++] = maparray[0];
                    }
                }
            }
        }

        int StdTetExp::v_GetVertexMap(const int localVertexId, bool useCoeffPacking)
        {
            ASSERTL0((GetEdgeBasisType(localVertexId)==LibUtilities::eModified_A)||
                     (GetEdgeBasisType(localVertexId)==LibUtilities::eModified_B)||
                     (GetEdgeBasisType(localVertexId)==LibUtilities::eModified_C),
                     "Mapping not defined for this type of basis");

            int localDOF = 0;
            if(useCoeffPacking == true) // follow packing of coefficients i.e q,r,p
            {
                switch(localVertexId)
                {
                case 0:
                    {
                        localDOF = GetMode(0,0,0);
                        break;
                    }
                case 1:
                    {
                        localDOF = GetMode(0,0,1);
                        break;
                    }
                case 2:
                    {
                        localDOF = GetMode(0,1,0);
                        break;
                    }
                case 3:
                    {
                        localDOF = GetMode(1,0,0);
                        break;
                    }
                default:
                    {
                        ASSERTL0(false,"Vertex ID must be between 0 and 3");
                        break;
                    }
                }
            }
            else
            {
                switch(localVertexId)
                {
                case 0:
                    {
                        localDOF = GetMode(0,0,0);
                        break;
                    }
                case 1:
                    {
                        localDOF = GetMode(1,0,0);
                        break;
                    }
                case 2:
                    {
                        localDOF = GetMode(0,1,0);
                        break;
                    }
                case 3:
                    {
                    localDOF = GetMode(0,0,1);
                    break;
                    }
                default:
                    {
                        ASSERTL0(false,"Vertex ID must be between 0 and 3");
                        break;
                    }
                }

            }

            return localDOF;
        }

        /**
         * Maps interior modes of an edge to the elemental modes.
         */
        void StdTetExp::v_GetEdgeInteriorMap(
            const int                  eid,
            const Orientation      edgeOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray)
        {
            int i;
            const int P = m_base[0]->GetNumModes();
            const int Q = m_base[1]->GetNumModes();
            const int R = m_base[2]->GetNumModes();

            const int nEdgeIntCoeffs = v_GetEdgeNcoeffs(eid)-2;

            if(maparray.size() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }
            else
            {
            	fill( maparray.get(), maparray.get() + nEdgeIntCoeffs, 0);
            }

            if(signarray.size() != nEdgeIntCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeIntCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nEdgeIntCoeffs, 1 );
            }

            switch (eid)
            {
                case 0:
                    for (i = 0; i < P-2; ++i)
                    {
                        maparray[i] = GetMode(i+2, 0, 0);
                    }
                    if(edgeOrient==eBackwards)
                    {
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                case 1:
                    for (i = 0; i < Q-2; ++i)
                    {
                        maparray[i] = GetMode(1, i+1, 0);
                    }
                    if(edgeOrient==eBackwards)
                    {
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                case 2:
                    for (i = 0; i < Q-2; ++i)
                    {
                        maparray[i] = GetMode(0, i+2, 0);
                    }
                    if(edgeOrient==eBackwards)
                    {
                    	for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                        	signarray[i] = -1;
                        }
                    }
                    break;
                case 3:
                    for (i = 0; i < R-2; ++i)
                    {
                    	maparray[i] = GetMode(0, 0, i+2);
                    }
                    if(edgeOrient==eBackwards)
                    {
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                case 4:
                    for (i = 0; i < R - 2; ++i)
                    {
                    	maparray[i] = GetMode(1, 0, i+1);
                    }
                    if(edgeOrient==eBackwards)
                    {
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                case 5:
                    for (i = 0; i < R-2; ++i)
                    {
                    	maparray[i] = GetMode(0, 1, i+1);
                    }
                    if(edgeOrient==eBackwards)
                    {
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }
                    break;
                default:
                    ASSERTL0(false, "Edge not defined.");
                    break;
            }
        }

        void StdTetExp::v_GetFaceInteriorMap(
            const int                  fid,
            const Orientation      faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray)
        {
            int i, j, idx, k;
            const int P = m_base[0]->GetNumModes();
            const int Q = m_base[1]->GetNumModes();
            const int R = m_base[2]->GetNumModes();

            const int nFaceIntCoeffs = v_GetFaceIntNcoeffs(fid);

            if(maparray.size() != nFaceIntCoeffs)
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

            switch (fid)
            {
                case 0:
                    idx = 0;
                    for (i = 2; i < P-1; ++i)
                    {
                        for (j = 1; j < Q-i; ++j)
                        {
                            if ((int)faceOrient == 7)
                            {
                                signarray[idx] = (i%2 ? -1 : 1);
                            }
                            maparray[idx++] = GetMode(i,j,0);
                        }
                    }
                    break;
                case 1:
                    idx = 0;
                    for (i = 2; i < P; ++i)
                    {
                        for (k = 1; k < R-i; ++k)
                        {
                            if ((int)faceOrient == 7)
                            {
                                signarray[idx] = (i%2 ? -1: 1);
                            }
                            maparray[idx++] = GetMode(i,0,k);
                        }
                    }
                    break;
                case 2:
                    idx = 0;
                    for (j = 1; j < Q-2; ++j)
                    {
                        for (k = 1; k < R-1-j; ++k)
                        {
                            if ((int)faceOrient == 7)
                            {
                                signarray[idx] = ((j+1)%2 ? -1: 1);
                            }
                            maparray[idx++] = GetMode(1,j,k);
                        }
                    }
                    break;
                case 3:
                    idx = 0;
                    for (j = 2; j < Q-1; ++j)
                    {
                        for (k = 1; k < R-j; ++k)
                        {
                            if ((int)faceOrient == 7)
                            {
                                signarray[idx] = (j%2 ? -1: 1);
                            }
                            maparray[idx++] = GetMode(0,j,k);
                        }
                    }
                    break;
                default:
                    ASSERTL0(false, "Face interior map not available.");
                    break;
            }
        }

        /**
         * List of all interior modes in the expansion.
         */
        void StdTetExp::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int P = m_base[0]->GetNumModes();
            int Q = m_base[1]->GetNumModes();
            int R = m_base[2]->GetNumModes();

            int nIntCoeffs = m_ncoeffs - NumBndryCoeffs();

            if(outarray.size() != nIntCoeffs)
            {
                outarray = Array<OneD, unsigned int>(nIntCoeffs);
            }

            int idx = 0;
            for (int i = 2; i < P-2; ++i)
            {
            	for (int j = 1; j < Q-i-1; ++j)
            	{
                    for (int k = 1; k < R-i-j; ++k)
                    {
                        outarray[idx++] = GetMode(i,j,k);
                    }
            	}
            }
        }

        /**
         * List of all boundary modes in the the expansion.
         */
        void StdTetExp::v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int P = m_base[0]->GetNumModes();
            int Q = m_base[1]->GetNumModes();
            int R = m_base[2]->GetNumModes();

            int i,j,k;
            int idx = 0;

            int nBnd = NumBndryCoeffs();

            if (outarray.size() != nBnd)
            {
                outarray = Array<OneD, unsigned int>(nBnd);
            }

            for (i = 0; i < P; ++i)
            {
            	// First two Q-R planes are entirely boundary modes
            	if (i < 2)
            	{
                    for (j = 0; j < Q-i; j++)
                    {
                        for (k = 0; k < R-i-j; ++k)
                        {
                            outarray[idx++] = GetMode(i,j,k);
                        }
                    }
            	}
            	// Remaining Q-R planes contain boundary modes on bottom and
            	// left edge.
            	else
            	{
                    for (k = 0; k < R-i; ++k)
                    {
                        outarray[idx++] = GetMode(i,0,k);
                    }
                    for (j = 1; j < Q-i; ++j)
                    {
                        outarray[idx++] = GetMode(i,j,0);
                    }
            	}
            }
        }


        //---------------------------------------
        // Wrapper functions
        //---------------------------------------
        DNekMatSharedPtr StdTetExp::v_GenMatrix(const StdMatrixKey &mkey)
        {

            MatrixType mtype   = mkey.GetMatrixType();

            DNekMatSharedPtr Mat;

            switch(mtype)
            {
            case ePhysInterpToEquiSpaced:
                {
                    int nq0 = m_base[0]->GetNumPoints();
                    int nq1 = m_base[1]->GetNumPoints();
                    int nq2 = m_base[2]->GetNumPoints();
                    int nq;

                    // take definition from key
                    if(mkey.ConstFactorExists(eFactorConst))
                    {
                        nq = (int) mkey.GetConstFactor(eFactorConst);
                    }
                    else
                    {
                        nq = max(nq0,max(nq1,nq2));
                    }

                    int neq = LibUtilities::StdTetData::
                                            getNumberOfCoefficients(nq,nq,nq);
                    Array<OneD, Array<OneD, NekDouble> > coords(neq);
                    Array<OneD, NekDouble>    coll(3);
                    Array<OneD, DNekMatSharedPtr> I(3);
                    Array<OneD, NekDouble> tmp(nq0);

                    Mat = MemoryManager<DNekMat>::
                                    AllocateSharedPtr(neq, nq0 * nq1 * nq2);
                    int cnt = 0;

                    for(int i = 0; i < nq; ++i)
                    {
                        for(int j = 0; j < nq-i; ++j)
                        {
                            for(int k = 0; k < nq-i-j; ++k,++cnt)
                            {
                                coords[cnt] = Array<OneD, NekDouble>(3);
                                coords[cnt][0] = -1.0 + 2*k/(NekDouble)(nq-1);
                                coords[cnt][1] = -1.0 + 2*j/(NekDouble)(nq-1);
                                coords[cnt][2] = -1.0 + 2*i/(NekDouble)(nq-1);
                            }
                        }
                    }

                    for(int i = 0; i < neq; ++i)
                    {
                        LocCoordToLocCollapsed(coords[i],coll);

                        I[0] = m_base[0]->GetI(coll);
                        I[1] = m_base[1]->GetI(coll+1);
                        I[2] = m_base[2]->GetI(coll+2);

                        // interpolate first coordinate direction
                        NekDouble fac;
                        for( int k = 0; k < nq2; ++k)
                        {
                            for (int j  = 0; j < nq1; ++j)
                            {

                                fac = (I[1]->GetPtr())[j]*(I[2]->GetPtr())[k];
                                Vmath::Smul(nq0,fac,I[0]->GetPtr(),1,tmp,1);

                                Vmath::Vcopy(nq0, &tmp[0], 1,
                                             Mat->GetRawPtr()+
                                             k*nq0*nq1*neq+
                                             j*nq0*neq+i,neq);
                            }
                        }
                    }
                }
                break;
            default:
                {
                    Mat = StdExpansion::CreateGeneralMatrix(mkey);
                }
                break;
            }

            return Mat;
        }

        DNekMatSharedPtr StdTetExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return v_GenMatrix(mkey);
        }


        //---------------------------------------
        // Private helper functions
        //---------------------------------------

        /**
         * @brief Compute the mode number in the expansion for a particular
         * tensorial combination.
         *
         * Modes are numbered with the r index travelling fastest, followed by
         * q and then p, and each q-r plane is of size
         * (Q+1)*(Q+2)/2+max(0,R-Q-p)*Q. For example, when P=2, Q=3 and R=4
         * the indexing inside each q-r plane (with r increasing upwards and q
         * to the right) is:
         *
         * p = 0:      p = 1:       p = 2:
         * ----------------------------------
         * 4
         * 3 8         17
         * 2 7 11      16 20        26
         * 1 6 10 13   15 19 22     25 28
         * 0 5 9  12   14 18 21 23  24 27 29
         *
         * Note that in this element, we must have that \f$ P \leq Q \leq
         * R\f$.
         */
        int StdTetExp::GetMode(const int I, const int J, const int K)
        {
            const int Q = m_base[1]->GetNumModes();
            const int R = m_base[2]->GetNumModes();

            int i,j,q_hat,k_hat;
            int cnt = 0;

            // Traverse to q-r plane number I
            for (i = 0; i < I; ++i)
            {
                // Size of triangle part
                q_hat = min(Q,R-i);
                // Size of rectangle part
                k_hat = max(R-Q-i,0);
                cnt  += q_hat*(q_hat+1)/2 + k_hat*Q;
            }

            // Traverse to q column J
            q_hat = R-I;
            for (j = 0; j < J; ++j)
            {
                cnt += q_hat;
                q_hat--;
            }

            // Traverse up stacks to K
            cnt += K;

            return cnt;
        }

        void StdTetExp::v_MultiplyByStdQuadratureMetric(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int i, j;

            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();

            const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();
            const Array<OneD, const NekDouble>& w2 = m_base[2]->GetW();

            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();
            const Array<OneD, const NekDouble>& z2 = m_base[2]->GetZ();

            // multiply by integration constants
            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vmul(nquad0,(NekDouble*)&inarray[0]+i*nquad0,1,
                            w0.get(),1, &outarray[0]+i*nquad0,1);
            }

            switch(m_base[1]->GetPointsType())
            {
                // (1,0) Jacobi Inner product.
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                    for(j = 0; j < nquad2; ++j)
                    {
                        for(i = 0; i < nquad1; ++i)
                        {
                            Blas::Dscal(nquad0,0.5*w1[i], &outarray[0]+i*nquad0+
                                        j*nquad0*nquad1,1);
                        }
                    }
                    break;

                default:
                    for(j = 0; j < nquad2; ++j)
                    {
                        for(i = 0; i < nquad1; ++i)
                        {
                            Blas::Dscal(nquad0,
                                        0.5*(1-z1[i])*w1[i],
                                        &outarray[0]+i*nquad0 + j*nquad0*nquad1,
                                        1 );
                        }
                    }
                    break;
            }

            switch(m_base[2]->GetPointsType())
            {
                // (2,0) Jacobi inner product.
            case LibUtilities::eGaussRadauMAlpha2Beta0:
                for(i = 0; i < nquad2; ++i)
                {
                    Blas::Dscal(nquad0*nquad1, 0.25*w2[i],
                                &outarray[0]+i*nquad0*nquad1, 1);
                }
                break;
                // (1,0) Jacobi inner product.
            case LibUtilities::eGaussRadauMAlpha1Beta0:
                for(i = 0; i < nquad2; ++i)
                {
                    Blas::Dscal(nquad0*nquad1, 0.25*(1-z2[i])*w2[i],
                                &outarray[0]+i*nquad0*nquad1, 1);
                }
                break;
            default:
                for(i = 0; i < nquad2; ++i)
                {
                    Blas::Dscal(nquad0*nquad1,0.25*(1-z2[i])*(1-z2[i])*w2[i],
                                &outarray[0]+i*nquad0*nquad1,1);
                }
                break;
            }
        }

        void StdTetExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                             const StdMatrixKey &mkey)
        {
            //To do : 1) add a test to ensure 0 \leq SvvCutoff \leq 1.
            //        2) check if the transfer function needs an analytical
            //           Fourier transform.
            //        3) if it doesn't : find a transfer function that renders
            //           the if( cutoff_a ...) useless to reduce computational
            //           cost.
            //        4) add SVVDiffCoef to both models!!

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
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B,nmodes_b,pb);
            LibUtilities::BasisKey Bc(LibUtilities::eOrtho_C,nmodes_c,pc);

            StdTetExp OrthoExp(Ba,Bb,Bc);


            Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
            int i,j,k,cnt = 0;

            // project onto physical space.
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
                    for(int j = 0; j < nmodes_b-j; ++j)
                    {
                        NekDouble fac1 = std::max(
                                   pow((1.0*i)/(nmodes_a-1),cutoff*nmodes_a),
                                   pow((1.0*j)/(nmodes_b-1),cutoff*nmodes_b));

                        for(int k = 0; k < nmodes_c-i-j; ++k)
                        {
                            NekDouble fac = std::max(fac1,
                                   pow((1.0*k)/(nmodes_c-1),cutoff*nmodes_c));

                            orthocoeffs[cnt] *= SvvDiffCoeff * fac;
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
                    for(int j = 0; j < nmodes_b-j; ++j)
                    {
                        int maxij = max(i,j);

                        for(int k = 0; k < nmodes_c-i-j; ++k)
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

                //SVV filter paramaters (how much added diffusion
                //relative to physical one and fraction of modes from
                //which you start applying this added diffusion)

                NekDouble  SvvDiffCoeff = mkey.GetConstFactor(StdRegions::eFactorSVVDiffCoeff);
                NekDouble  SVVCutOff = mkey.GetConstFactor(StdRegions::eFactorSVVCutoffRatio);

                //Defining the cut of mode
                int cutoff_a = (int) (SVVCutOff*nmodes_a);
                int cutoff_b = (int) (SVVCutOff*nmodes_b);
                int cutoff_c = (int) (SVVCutOff*nmodes_c);
                int nmodes = min(min(nmodes_a,nmodes_b),nmodes_c);
                NekDouble cutoff = min(min(cutoff_a,cutoff_b),cutoff_c);
                NekDouble epsilon = 1;


                //------"New" Version August 22nd '13--------------------
                for(i = 0; i < nmodes_a; ++i)
                {
                    for(j = 0; j < nmodes_b-i; ++j)
                    {
                        for(k = 0; k < nmodes_c-i-j; ++k)
                        {
                            if(i + j + k >= cutoff)
                            {
                                orthocoeffs[cnt] *= ((SvvDiffCoeff)*exp(-(i+j+k-nmodes)*(i+j+k-nmodes)/((NekDouble)((i+j+k-cutoff+epsilon)*(i+j+k-cutoff+epsilon)))));
                            }
                            else
                            {
                                orthocoeffs[cnt] *= 0.0;
                            }
                            cnt++;
                        }
                    }
                }
            }

            // backward transform to physical space
            OrthoExp.BwdTrans(orthocoeffs,array);
        }


        void StdTetExp::v_ReduceOrderCoeffs(
            int                                 numMin,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            int nquad0   = m_base[0]->GetNumPoints();
            int nquad1   = m_base[1]->GetNumPoints();
            int nquad2   = m_base[2]->GetNumPoints();
            int nqtot    = nquad0 * nquad1 * nquad2;
            int nmodes0  = m_base[0]->GetNumModes();
            int nmodes1  = m_base[1]->GetNumModes();
            int nmodes2  = m_base[2]->GetNumModes();
            int numMax   = nmodes0;

            Array<OneD, NekDouble> coeff     (m_ncoeffs);
            Array<OneD, NekDouble> coeff_tmp1(m_ncoeffs, 0.0);
            Array<OneD, NekDouble> coeff_tmp2(m_ncoeffs, 0.0);
            Array<OneD, NekDouble> phys_tmp  (nqtot,     0.0);
            Array<OneD, NekDouble> tmp, tmp2, tmp3, tmp4;

            Vmath::Vcopy(m_ncoeffs,inarray,1,coeff_tmp2,1);

            const LibUtilities::PointsKey Pkey0 = m_base[0]->GetPointsKey();
            const LibUtilities::PointsKey Pkey1 = m_base[1]->GetPointsKey();
            const LibUtilities::PointsKey Pkey2 = m_base[2]->GetPointsKey();

            LibUtilities::BasisKey bortho0(LibUtilities::eOrtho_A,
                                           nmodes0, Pkey0);
            LibUtilities::BasisKey bortho1(LibUtilities::eOrtho_B,
                                           nmodes1, Pkey1);
            LibUtilities::BasisKey bortho2(LibUtilities::eOrtho_C,
                                           nmodes2, Pkey2);

            Vmath::Zero(m_ncoeffs, coeff_tmp2, 1);

            StdRegions::StdTetExpSharedPtr OrthoTetExp;
            OrthoTetExp = MemoryManager<StdRegions::StdTetExp>
                ::AllocateSharedPtr(bortho0, bortho1, bortho2);

            BwdTrans(inarray,phys_tmp);
            OrthoTetExp->FwdTrans(phys_tmp, coeff);

            Vmath::Zero(m_ncoeffs,outarray,1);

            // filtering
            int cnt = 0;
            for (int u = 0; u < numMin; ++u)
            {
                for (int i = 0; i < numMin-u; ++i)
                {
                    Vmath::Vcopy(numMin - u - i, tmp  = coeff      + cnt, 1,
                                                 tmp2 = coeff_tmp1 + cnt, 1);
                    cnt += numMax - u - i;
                }
                for (int i = numMin; i < numMax-u; ++i)
                {
                    cnt += numMax - u - i;
                }
            }

            OrthoTetExp->BwdTrans(coeff_tmp1,phys_tmp);
            FwdTrans(phys_tmp, outarray);
        }


        void StdTetExp::v_GetSimplexEquiSpacedConnectivity(
            Array<OneD, int> &conn,
            bool              standard)
        {
            boost::ignore_unused(standard);

            int np0 = m_base[0]->GetNumPoints();
            int np1 = m_base[1]->GetNumPoints();
            int np2 = m_base[2]->GetNumPoints();
            int np = max(np0,max(np1,np2));


            conn = Array<OneD, int>(4*(np-1)*(np-1)*(np-1));

            int row   = 0;
            int rowp1 = 0;
            int plane = 0;
            int row1   = 0;
            int row1p1 = 0;
            int planep1= 0;
            int cnt = 0;
            for(int i = 0; i < np-1; ++i)
            {
                planep1 += (np-i)*(np-i+1)/2;
                row    = 0; // current plane row offset
                rowp1  = 0; // current plane row plus one offset
                row1   = 0; // next plane row offset
                row1p1 = 0; // nex plane row plus one offset
                for(int j = 0; j < np-i-1; ++j)
                {
                    rowp1 += np-i-j;
                    row1p1 += np-i-j-1;
                    for(int k = 0; k < np-i-j-2; ++k)
                    {
                        conn[cnt++] = plane   + row   +k+1;
                        conn[cnt++] = plane   + row   +k;
                        conn[cnt++] = plane   + rowp1 +k;
                        conn[cnt++] = planep1 + row1  +k;

                        conn[cnt++] = plane   + row   +k+1;
                        conn[cnt++] = plane   + rowp1 +k+1;
                        conn[cnt++] = planep1 + row1  +k+1;
                        conn[cnt++] = planep1 + row1  +k;

                        conn[cnt++] = plane   + rowp1 +k+1;
                        conn[cnt++] = plane   + row   +k+1;
                        conn[cnt++] = plane   + rowp1 +k;
                        conn[cnt++] = planep1 + row1  +k;

                        conn[cnt++] = planep1 + row1  +k;
                        conn[cnt++] = planep1 + row1p1+k;
                        conn[cnt++] = plane   + rowp1 +k;
                        conn[cnt++] = plane   + rowp1 +k+1;

                        conn[cnt++] = planep1 + row1  +k;
                        conn[cnt++] = planep1 + row1p1+k;
                        conn[cnt++] = planep1 + row1  +k+1;
                        conn[cnt++] = plane   + rowp1 +k+1;

                        if(k < np-i-j-3)
                        {
                            conn[cnt++] = plane   + rowp1  +k+1;
                            conn[cnt++] = planep1 + row1p1 +k+1;
                            conn[cnt++] = planep1 + row1   +k+1;
                            conn[cnt++] = planep1 + row1p1 +k;
                        }
                    }

                    conn[cnt++] = plane   + row   + np-i-j-1;
                    conn[cnt++] = plane   + row   + np-i-j-2;
                    conn[cnt++] = plane   + rowp1 + np-i-j-2;
                    conn[cnt++] = planep1 + row1  + np-i-j-2;

                    row  += np-i-j;
                    row1 += np-i-j-1;
                }
                plane += (np-i)*(np-i+1)/2;
            }
        }

    }//end namespace
}//end namespace
