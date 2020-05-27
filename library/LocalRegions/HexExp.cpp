///////////////////////////////////////////////////////////////////////////////
//
// File HexExp.cpp
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
// Description: Methods for Hex expansion in local regoins
//
///////////////////////////////////////////////////////////////////////////////


#include <LocalRegions/HexExp.h>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <SpatialDomains/HexGeom.h>

using namespace std;

namespace Nektar
{
    namespace LocalRegions
    {
        /**
         * @class HexExp
         * Defines a hexahedral local expansion.
         */

        /**
	 * \brief Constructor using BasisKey class for quadrature points and
	 * order definition
	 *
         * @param   Ba          Basis key for first coordinate.
         * @param   Bb          Basis key for second coordinate.
         * @param   Bc          Basis key for third coordinate.
         */
        HexExp::HexExp(const LibUtilities::BasisKey &Ba,
                       const LibUtilities::BasisKey &Bb,
                       const LibUtilities::BasisKey &Bc,
                       const SpatialDomains::HexGeomSharedPtr &geom):
            StdExpansion  (Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes(),3,Ba,Bb,Bc),
            StdExpansion3D(Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes(),Ba,Bb,Bc),
            StdHexExp(Ba,Bb,Bc),
            Expansion     (geom),
            Expansion3D   (geom),
            m_matrixManager(
                std::bind(&HexExp::CreateMatrix, this, std::placeholders::_1),
                std::string("HexExpMatrix")),
            m_staticCondMatrixManager(
                std::bind(&HexExp::CreateStaticCondMatrix, this, std::placeholders::_1),
                std::string("HexExpStaticCondMatrix"))
        {
        }


        /**
	 * \brief Copy Constructor
	 *
         * @param   T           HexExp to copy.
         */
        HexExp::HexExp(const HexExp &T):
            StdExpansion(T),
            StdExpansion3D(T),
            StdHexExp(T),
            Expansion(T),
            Expansion3D(T),
            m_matrixManager(T.m_matrixManager),
            m_staticCondMatrixManager(T.m_staticCondMatrixManager)
        {
        }

        /**
	 * \brief Destructor
	 */
        HexExp::~HexExp()
        {
        }


        //-----------------------------
        // Integration Methods
        //-----------------------------
        /**
	 * \brief Integrate the physical point list \a inarray over region
	 *
         * @param   inarray     definition of function to be returned at
         *                      quadrature points of expansion.
         * @returns \f$\int^1_{-1}\int^1_{-1} \int^1_{-1}
         *   u(\eta_1, \eta_2, \eta_3) J[i,j,k] d \eta_1 d \eta_2 d \eta_3 \f$
         * where \f$inarray[i,j,k] = u(\eta_{1i},\eta_{2j},\eta_{3k}) \f$
         * and \f$ J[i,j,k] \f$ is the Jacobian evaluated at the quadrature
         * point.
         */
        NekDouble HexExp::v_Integral(
                 const Array<OneD, const NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac(GetPointsKeys());
            NekDouble returnVal;
            Array<OneD,NekDouble> tmp(nquad0*nquad1*nquad2);

            // multiply inarray with Jacobian

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,
                            (NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,(NekDouble) jac[0],
                            (NekDouble*)&inarray[0],1,&tmp[0],1);
            }

            // call StdHexExp version;
            returnVal = StdHexExp::v_Integral(tmp);

            return  returnVal;
        }


        //-----------------------------
        // Differentiation Methods
        //-----------------------------
        /**
	 * \brief Calculate the derivative of the physical points
	 *
         * For Hexahedral region can use the Tensor_Deriv function defined
         * under StdExpansion.
         * @param   inarray     Input array
         * @param   out_d0      Derivative of \a inarray in first direction.
         * @param   out_d1      Derivative of \a inarray in second direction.
         * @param   out_d2      Derivative of \a inarray in third direction.
         */
        void HexExp::v_PhysDeriv(
                const Array<OneD, const NekDouble> & inarray,
                      Array<OneD,NekDouble> &out_d0,
                      Array<OneD,NekDouble> &out_d1,
                      Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            int    ntot   = nquad0 * nquad1 * nquad2;

            Array<TwoD, const NekDouble> df =
                                m_metricinfo->GetDerivFactors(GetPointsKeys());
            Array<OneD,NekDouble> Diff0 = Array<OneD,NekDouble>(ntot);
            Array<OneD,NekDouble> Diff1 = Array<OneD,NekDouble>(ntot);
            Array<OneD,NekDouble> Diff2 = Array<OneD,NekDouble>(ntot);

            StdHexExp::v_PhysDeriv(inarray, Diff0, Diff1, Diff2);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.size())
                {
                    Vmath::Vmul (ntot,&df[0][0],1,&Diff0[0],1, &out_d0[0], 1);
                    Vmath::Vvtvp(ntot,&df[1][0],1,&Diff1[0],1, &out_d0[0], 1,
                                                                 &out_d0[0],1);
                    Vmath::Vvtvp(ntot,&df[2][0],1,&Diff2[0],1, &out_d0[0], 1,
                                                                 &out_d0[0],1);
                }

                if(out_d1.size())
                {
                    Vmath::Vmul (ntot,&df[3][0],1,&Diff0[0],1, &out_d1[0], 1);
                    Vmath::Vvtvp(ntot,&df[4][0],1,&Diff1[0],1, &out_d1[0], 1,
                                                                 &out_d1[0],1);
                    Vmath::Vvtvp(ntot,&df[5][0],1,&Diff2[0],1, &out_d1[0], 1,
                                                                 &out_d1[0],1);
                }

                if(out_d2.size())
                {
                    Vmath::Vmul (ntot,&df[6][0],1,&Diff0[0],1, &out_d2[0], 1);
                    Vmath::Vvtvp(ntot,&df[7][0],1,&Diff1[0],1, &out_d2[0], 1,
                                                                 &out_d2[0],1);
                    Vmath::Vvtvp(ntot,&df[8][0],1,&Diff2[0],1, &out_d2[0], 1,
                                                                 &out_d2[0],1);
                }
            }
            else // regular geometry
            {
                if(out_d0.size())
                {
                    Vmath::Smul (ntot,df[0][0],&Diff0[0],1, &out_d0[0], 1);
                    Blas::Daxpy (ntot,df[1][0],&Diff1[0],1, &out_d0[0], 1);
                    Blas::Daxpy (ntot,df[2][0],&Diff2[0],1, &out_d0[0], 1);
                }

                if(out_d1.size())
                {
                    Vmath::Smul (ntot,df[3][0],&Diff0[0],1, &out_d1[0], 1);
                    Blas::Daxpy (ntot,df[4][0],&Diff1[0],1, &out_d1[0], 1);
                    Blas::Daxpy (ntot,df[5][0],&Diff2[0],1, &out_d1[0], 1);
                }

                if(out_d2.size())
                {
                    Vmath::Smul (ntot,df[6][0],&Diff0[0],1, &out_d2[0], 1);
                    Blas::Daxpy (ntot,df[7][0],&Diff1[0],1, &out_d2[0], 1);
                    Blas::Daxpy (ntot,df[8][0],&Diff2[0],1, &out_d2[0], 1);
                }
            }
        }


        /**
	 * \brief Calculate the derivative of the physical points in a single
         * direction.
	 *
         * @param   dir         Direction in which to compute derivative.
         *                      Valid values are 0, 1, 2.
         * @param   inarray     Input array.
         * @param   outarray    Output array.
         */
        void HexExp::v_PhysDeriv(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble>& outarray)
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


        void HexExp::v_PhysDirectionalDeriv(
            const Array<OneD, const NekDouble>& inarray,
            const Array<OneD, const NekDouble>& direction,
                  Array<OneD, NekDouble> & outarray)
        {

            int    shapedim = 3;
            int    nquad0   = m_base[0]->GetNumPoints();
            int    nquad1   = m_base[1]->GetNumPoints();
            int    nquad2   = m_base[2]->GetNumPoints();
            int    ntot     = nquad0 * nquad1 * nquad2;

            Array<TwoD, const NekDouble> df =
                    m_metricinfo->GetDerivFactors(GetPointsKeys());
            Array<OneD,NekDouble> Diff0 = Array<OneD,NekDouble>(ntot);
            Array<OneD,NekDouble> Diff1 = Array<OneD,NekDouble>(ntot);
            Array<OneD,NekDouble> Diff2 = Array<OneD,NekDouble>(ntot);

            StdHexExp::v_PhysDeriv(inarray, Diff0, Diff1, Diff2);

            Array<OneD, Array<OneD, NekDouble> > dfdir(shapedim);
            Expansion::ComputeGmatcdotMF(df,direction,dfdir);

            Vmath::Vmul (ntot, &dfdir[0][0], 1,
                               &Diff0[0],    1,
                               &outarray[0], 1 );
            Vmath::Vvtvp(ntot, &dfdir[1][0], 1,
                               &Diff1[0],    1,
                               &outarray[0], 1,
                               &outarray[0], 1 );
            Vmath::Vvtvp(ntot, &dfdir[2][0], 1,
                               &Diff2[0],    1,
                               &outarray[0], 1,
                               &outarray[0], 1 );
        }

        //-----------------------------
        // Transforms
        //-----------------------------

        /**
	 * \brief Forward transform from physical quadrature space stored in \a
         * inarray and evaluate the expansion coefficients and store in
         * \a (this)->_coeffs
	 *
         * @param   inarray     Input array
         * @param   outarray    Output array
         */
        void HexExp::v_FwdTrans(
                const Array<OneD, const NekDouble> & inarray,
                      Array<OneD,NekDouble> &outarray)
        {
            if( m_base[0]->Collocation() && m_base[1]->Collocation()
                    && m_base[2]->Collocation())
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&outarray[0],1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass,
                                              DetShapeType(),*this);
                DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];

                // copy inarray in case inarray == outarray
                DNekVec in (m_ncoeffs,outarray);
                DNekVec out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;
            }
        }


        //-----------------------------
        // Inner product functions
        //-----------------------------

        /**
	 * \brief Calculate the inner product of inarray with respect to the
	 * elements basis.
	 *
         * @param   inarray     Input array of physical space data.
         * @param   outarray    Output array of data.
         */
        void HexExp::v_IProductWRTBase(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray)
        {
            HexExp::v_IProductWRTBase_SumFac(inarray, outarray);
        }

        /**
	 * \brief Calculate the inner product of inarray with respect to the
	 * given basis B = base0 * base1 * base2.
	 *
         * \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta}
         * & = & \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
         *     \psi_{p}^{a} (\xi_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{r}^{a}
         *     (\xi_{3k}) w_i w_j w_k u(\xi_{1,i} \xi_{2,j} \xi_{3,k})
         * J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\xi_{1,i})
         *     \sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2}
         *     \psi_{r}^a u(\xi_{1i},\xi_{2j},\xi_{3k})
         * J_{i,j,k} \end{array} \f$ \n
         * where
         * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3)
         *    = \psi_p^a ( \xi_1) \psi_{q}^a (\xi_2) \psi_{r}^a (\xi_3) \f$ \n
         * which can be implemented as \n
         * \f$f_{r} (\xi_{3k})
         *    = \sum_{k=0}^{nq_3} \psi_{r}^a u(\xi_{1i},\xi_{2j},\xi_{3k})
         * J_{i,j,k} = {\bf B_3 U}   \f$ \n
         * \f$ g_{q} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j})
         *                          f_{r} (\xi_{3k})  = {\bf B_2 F}  \f$ \n
         * \f$ (\phi_{pqr}, u)_{\delta}
         *    = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k})
         *    = {\bf B_1 G} \f$
         *
         * @param   base0       Basis to integrate wrt in first dimension.
         * @param   base1       Basis to integrate wrt in second dimension.
         * @param   base2       Basis to integrate wrt in third dimension.
         * @param   inarray     Input array.
         * @param   outarray    Output array.
         * @param   coll_check  (not used)
         */
        void HexExp::v_IProductWRTBase_SumFac(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,       NekDouble> &outarray,
                bool multiplybyweights)
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

                MultiplyByQuadratureMetric(inarray, tmp);
               IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                             m_base[1]->GetBdata(),
                                             m_base[2]->GetBdata(),
                                             tmp,outarray,wsp,
                                             true,true,true);
            }
            else
            {
               IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                            m_base[1]->GetBdata(),
                                            m_base[2]->GetBdata(),
                                            inarray,outarray,wsp,
                                            true,true,true);

            }
        }

        void HexExp::v_IProductWRTDerivBase(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> & outarray)
        {
            HexExp::IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }


        /**
         * @brief Calculates the inner product \f$ I_{pqr} = (u,
         * \partial_{x_i} \phi_{pqr}) \f$.
         *
         * The derivative of the basis functions is performed using the chain
         * rule in order to incorporate the geometric factors. Assuming that
         * the basis functions are a tensor product
         * \f$\phi_{pqr}(\xi_1,\xi_2,\xi_3) =
         * \phi_1(\xi_1)\phi_2(\xi_2)\phi_3(\xi_3)\f$, in the hexahedral
         * element, this is straightforward and yields the result
         *
         * \f[
         * I_{pqr} = \sum_{k=1}^3 \left(u, \frac{\partial u}{\partial \xi_k}
         * \frac{\partial \xi_k}{\partial x_i}\right)
         * \f]
         *
         * @param dir       Direction in which to take the derivative.
         * @param inarray   The function \f$ u \f$.
         * @param outarray  Value of the inner product.
         */
        void HexExp::IProductWRTDerivBase_SumFac(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> & outarray)
        {
            ASSERTL1((dir==0)||(dir==1)||(dir==2),"Invalid direction.");

            const int nq0 = m_base[0]->GetNumPoints();
            const int nq1 = m_base[1]->GetNumPoints();
            const int nq2 = m_base[2]->GetNumPoints();
            const int nq  = nq0*nq1*nq2;
            const int nm0 = m_base[0]->GetNumModes();
            const int nm1 = m_base[1]->GetNumModes();

            const Array<TwoD, const NekDouble>& df =
                                m_metricinfo->GetDerivFactors(GetPointsKeys());

            Array<OneD, NekDouble> alloc(4*nq + m_ncoeffs + nm0*nq2*(nq1+nm1));
            Array<OneD, NekDouble> tmp1 (alloc);               // Quad metric
            Array<OneD, NekDouble> tmp2 (alloc +   nq);        // Dir1 metric
            Array<OneD, NekDouble> tmp3 (alloc + 2*nq);        // Dir2 metric
            Array<OneD, NekDouble> tmp4 (alloc + 3*nq);        // Dir3 metric
            Array<OneD, NekDouble> tmp5 (alloc + 4*nq);        // iprod tmp
            Array<OneD, NekDouble> wsp  (tmp5  +   m_ncoeffs); // Wsp

            MultiplyByQuadratureMetric(inarray, tmp1);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nq,&df[3*dir][0],  1,tmp1.get(),1,tmp2.get(),1);
                Vmath::Vmul(nq,&df[3*dir+1][0],1,tmp1.get(),1,tmp3.get(),1);
                Vmath::Vmul(nq,&df[3*dir+2][0],1,tmp1.get(),1,tmp4.get(),1);
            }
            else
            {
                Vmath::Smul(nq, df[3*dir][0],  tmp1.get(),1,tmp2.get(), 1);
                Vmath::Smul(nq, df[3*dir+1][0],tmp1.get(),1,tmp3.get(), 1);
                Vmath::Smul(nq, df[3*dir+2][0],tmp1.get(),1,tmp4.get(), 1);
            }

            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetBdata(),
                                         tmp2,outarray,wsp,
                                         false,true,true);

            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetDbdata(),
                                         m_base[2]->GetBdata(),
                                         tmp3,tmp5,wsp,
                                         true,false,true);
            Vmath::Vadd(m_ncoeffs, tmp5, 1, outarray, 1, outarray, 1);

            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetDbdata(),
                                         tmp4,tmp5,wsp,
                                         true,true,false);
            Vmath::Vadd(m_ncoeffs, tmp5, 1, outarray, 1, outarray, 1);
        }


        void HexExp::IProductWRTDerivBase_MatOp(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdRegions::MatrixType mtype = StdRegions::eIProductWRTDerivBase0;

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

            MatrixKey      iprodmatkey(mtype,DetShapeType(),*this);
            DNekScalMatSharedPtr iprodmat = m_matrixManager[iprodmatkey];

            Blas::Dgemv('N',m_ncoeffs,nq,iprodmat->Scale(),(iprodmat->GetOwnedMatrix())->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }


        /**
         *
         * @param dir       Vector direction in which to take the derivative.
         * @param inarray   The function \f$ u \f$.
         * @param outarray  Value of the inner product.
         */
        void HexExp::IProductWRTDirectionalDerivBase_SumFac(
                const Array<OneD, const NekDouble>& direction,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> & outarray)
        {
            int shapedim  = 3;
            const int nq0 = m_base[0]->GetNumPoints();
            const int nq1 = m_base[1]->GetNumPoints();
            const int nq2 = m_base[2]->GetNumPoints();
            const int nq  = nq0*nq1*nq2;
            const int nm0 = m_base[0]->GetNumModes();
            const int nm1 = m_base[1]->GetNumModes();

            const Array<TwoD, const NekDouble>& df =
            m_metricinfo->GetDerivFactors(GetPointsKeys());

            Array<OneD, NekDouble> alloc(4*nq + m_ncoeffs + nm0*nq2*(nq1+nm1));
            Array<OneD, NekDouble> tmp1 (alloc);               // Quad metric
            Array<OneD, NekDouble> tmp2 (alloc +   nq);        // Dir1 metric
            Array<OneD, NekDouble> tmp3 (alloc + 2*nq);        // Dir2 metric
            Array<OneD, NekDouble> tmp4 (alloc + 3*nq);        // Dir3 metric
            Array<OneD, NekDouble> tmp5 (alloc + 4*nq);        // iprod tmp
            Array<OneD, NekDouble> wsp  (tmp5  +   m_ncoeffs); // Wsp

            MultiplyByQuadratureMetric(inarray, tmp1);

            Array<OneD, Array<OneD, NekDouble> > dfdir(shapedim);
            Expansion::ComputeGmatcdotMF(df,direction,dfdir);

            Vmath::Vmul(nq,&dfdir[0][0],1,tmp1.get(),1,tmp2.get(),1);
            Vmath::Vmul(nq,&dfdir[1][0],1,tmp1.get(),1,tmp3.get(),1);
            Vmath::Vmul(nq,&dfdir[2][0],1,tmp1.get(),1,tmp4.get(),1);

            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetBdata(),
                                         tmp2,outarray,wsp,
                                         false,true,true);

            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetDbdata(),
                                         m_base[2]->GetBdata(),
                                         tmp3,tmp5,wsp,
                                         true,false,true);

            Vmath::Vadd(m_ncoeffs, tmp5, 1, outarray, 1, outarray, 1);

            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetDbdata(),
                                         tmp4,tmp5,wsp,
                                         true,true,false);

            Vmath::Vadd(m_ncoeffs, tmp5, 1, outarray, 1, outarray, 1);
        }


        //-----------------------------
        // Evaluation functions
        //-----------------------------


        /**
         * Given the local cartesian coordinate \a Lcoord evaluate the
         * value of physvals at this point by calling through to the
         * StdExpansion method
         */
        NekDouble HexExp::v_StdPhysEvaluate(
                const Array<OneD, const NekDouble> &Lcoord,
                const Array<OneD, const NekDouble> &physvals)
        {
            // Evaluate point in local coordinates.
            return StdHexExp::v_PhysEvaluate(Lcoord,physvals);
        }

        NekDouble HexExp::v_PhysEvaluate(
                const Array<OneD, const NekDouble> &coord,
                const Array<OneD, const NekDouble> & physvals)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(3);

            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);
            return StdHexExp::v_PhysEvaluate(Lcoord, physvals);
        }

        StdRegions::StdExpansionSharedPtr HexExp::v_GetStdExp(void) const
        {
            return MemoryManager<StdRegions::StdHexExp>
                ::AllocateSharedPtr(m_base[0]->GetBasisKey(),
                                    m_base[1]->GetBasisKey(),
                                    m_base[2]->GetBasisKey());
        }


        StdRegions::StdExpansionSharedPtr HexExp::v_GetLinStdExp(void) const
        {
            LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(),
                           2, m_base[0]->GetPointsKey());
            LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(),
                           2, m_base[1]->GetPointsKey());
            LibUtilities::BasisKey bkey2(m_base[2]->GetBasisType(),
                           2, m_base[2]->GetPointsKey());

            return MemoryManager<StdRegions::StdHexExp>
                ::AllocateSharedPtr( bkey0, bkey1, bkey2);
        }

        /**
	 * \brief Retrieves the physical coordinates of a given set of
         * reference coordinates.
	 *
         * @param   Lcoords     Local coordinates in reference space.
         * @param   coords      Corresponding coordinates in physical space.
         */
        void HexExp::v_GetCoord(
                const Array<OneD, const NekDouble> &Lcoords,
                      Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[0] <= 1.0 &&
                     Lcoords[1] >= -1.0 && Lcoords[1] <= 1.0 &&
                     Lcoords[2] >= -1.0 && Lcoords[2] <= 1.0,
                     "Local coordinates are not in region [-1,1]");

              m_geom->FillGeom();

            for(i = 0; i < m_geom->GetCoordim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }

        void HexExp::v_GetCoords(
            Array<OneD, NekDouble> &coords_0,
            Array<OneD, NekDouble> &coords_1,
            Array<OneD, NekDouble> &coords_2)
        {
            Expansion::v_GetCoords(coords_0, coords_1, coords_2);
        }

        //-----------------------------
        // Helper functions
        //-----------------------------

        /// Return the region shape using the enum-list of ShapeType
        LibUtilities::ShapeType HexExp::v_DetShapeType() const
        {
            return LibUtilities::eHexahedron;
        }


        void HexExp::v_ExtractDataToCoeffs(
                const NekDouble *data,
                const std::vector<unsigned int > &nummodes,
                const int mode_offset,
                NekDouble * coeffs,
                std::vector<LibUtilities::BasisType> &fromType)
        {
            int data_order0 = nummodes[mode_offset];
            int fillorder0  = min(m_base[0]->GetNumModes(),data_order0);
            int data_order1 = nummodes[mode_offset+1];
            int order1      = m_base[1]->GetNumModes();
            int fillorder1  = min(order1,data_order1);
            int data_order2 = nummodes[mode_offset+2];
            int order2      = m_base[2]->GetNumModes();
            int fillorder2  = min(order2,data_order2);

            // Check if same basis
            if (fromType[0] != m_base[0]->GetBasisType() ||
                fromType[1] != m_base[1]->GetBasisType() ||
                fromType[2] != m_base[2]->GetBasisType())
            {
                // Construct a hex with the appropriate basis type at our
                // quadrature points, and one more to do a forwards
                // transform. We can then copy the output to coeffs.
                StdRegions::StdHexExp tmpHex(
                    LibUtilities::BasisKey(
                        fromType[0], data_order0, m_base[0]->GetPointsKey()),
                    LibUtilities::BasisKey(
                        fromType[1], data_order1, m_base[1]->GetPointsKey()),
                    LibUtilities::BasisKey(
                        fromType[2], data_order2, m_base[2]->GetPointsKey()));
                StdRegions::StdHexExp tmpHex2(m_base[0]->GetBasisKey(),
                                              m_base[1]->GetBasisKey(),
                                              m_base[2]->GetBasisKey());

                Array<OneD, const NekDouble> tmpData(tmpHex.GetNcoeffs(), data);
                Array<OneD, NekDouble> tmpBwd(tmpHex2.GetTotPoints());
                Array<OneD, NekDouble> tmpOut(tmpHex2.GetNcoeffs());

                tmpHex.BwdTrans(tmpData, tmpBwd);
                tmpHex2.FwdTrans(tmpBwd, tmpOut);
                Vmath::Vcopy(tmpOut.size(), &tmpOut[0], 1, coeffs, 1);

                return;
            }

            switch(m_base[0]->GetBasisType())
            {
            case LibUtilities::eModified_A:
                {
                    int i,j;
                    int cnt  = 0;
                    int cnt1 = 0;

                    ASSERTL1(m_base[1]->GetBasisType() ==
                             LibUtilities::eModified_A,
                             "Extraction routine not set up for this basis");
                    ASSERTL1(m_base[2]->GetBasisType() ==
                             LibUtilities::eModified_A,
                             "Extraction routine not set up for this basis");

                    Vmath::Zero(m_ncoeffs,coeffs,1);
                    for(j = 0; j < fillorder0; ++j)
                    {
                        for(i = 0; i < fillorder1; ++i)
                        {
                            Vmath::Vcopy(fillorder2, &data[cnt],    1,
                                                     &coeffs[cnt1], 1);
                            cnt  += data_order2;
                            cnt1 += order2;
                        }

                        // count out data for j iteration
                        for(i = fillorder1; i < data_order1; ++i)
                        {
                            cnt += data_order2;
                        }

                        for(i = fillorder1; i < order1; ++i)
                        {
                            cnt1 += order2;
                        }
                    }
                    break;
                }
                case LibUtilities::eGLL_Lagrange:
                {
                    LibUtilities::PointsKey
                        p0(nummodes[0], LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey
                        p1(nummodes[1], LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey
                        p2(nummodes[2], LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey t0(
                        m_base[0]->GetNumModes(),
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey t1(
                        m_base[1]->GetNumModes(),
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey t2(
                        m_base[2]->GetNumModes(),
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::Interp3D(p0, p1, p2, data, t0, t1, t2, coeffs);
                }
                    break;
            default:
                ASSERTL0(false, "basis is either not set up or not "
                                "hierarchicial");
            }
        }

        bool HexExp::v_GetFaceDGForwards(const int i) const
        {
            StdRegions::Orientation fo = GetGeom3D()->GetForient(i);

            return fo == StdRegions::eDir1FwdDir1_Dir2FwdDir2 ||
                   fo == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                   fo == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                   fo == StdRegions::eDir1FwdDir2_Dir2BwdDir1;
        }

        void HexExp::v_GetFacePhysMap(const int               face,
                                      Array<OneD, int>        &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();

            int nq0 = 0;
            int nq1 = 0;

            switch(face)
            {
                case 0:
                    nq0 = nquad0;
                    nq1 = nquad1;

                    //Directions A and B positive
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    for (int i = 0; i < nquad0*nquad1; ++i)
                    {
                        outarray[i] = i;
                    }

                    break;
                case 1:
                    nq0 = nquad0;
                    nq1 = nquad2;
                    //Direction A and B positive
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    //Direction A and B positive
                    for (int k = 0; k < nquad2; k++)
                    {
                        for(int i = 0; i < nquad0; ++i)
                        {
                            outarray[k*nquad0 + i] = nquad0*nquad1*k + i;
                        }
                    }
                    break;
                case 2:
                    nq0 = nquad1;
                    nq1 = nquad2;

                    //Direction A and B positive
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    for (int i = 0; i < nquad1*nquad2; i++)
                    {
                        outarray[i] = nquad0-1 + i*nquad0;
                    }
                    break;
                case 3:
                    nq0 = nquad0;
                    nq1 = nquad2;

                    //Direction A and B positive
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    for (int k = 0; k < nquad2; k++)
                    {
                        for (int i = 0; i < nquad0; i++)
                        {
                            outarray[k*nquad0 + i] = (nquad0*(nquad1-1))+(k*nquad0*nquad1) + i;
                        }
                    }
                    break;
                case 4:
                    nq0 = nquad1;
                    nq1 = nquad2;

                    //Direction A and B positive
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    for (int i = 0; i < nquad1*nquad2; i++)
                    {
                        outarray[i] = i*nquad0;
                    }
                    break;
                case 5:
                    nq0 = nquad0;
                    nq1 = nquad1;
                    //Directions A and B positive
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    for (int i = 0; i < nquad0*nquad1; i++)
                    {
                        outarray[i] = nquad0*nquad1*(nquad2-1) + i;
                    }

                    break;
                default:
                    ASSERTL0(false,"face value (> 5) is out of range");
                    break;
            }

        }

        void HexExp::v_ComputeFaceNormal(const int face)
        {
            int i;
            const SpatialDomains::GeomFactorsSharedPtr & geomFactors =
                GetGeom()->GetMetricInfo();
            SpatialDomains::GeomType type = geomFactors->GetGtype();

            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
            for(i = 0; i < ptsKeys.size(); ++i)
            {
                // Need at least 2 points for computing normals
                if (ptsKeys[i].GetNumPoints() == 1)
                {
                    LibUtilities::PointsKey pKey(2, ptsKeys[i].GetPointsType());
                    ptsKeys[i] = pKey;
                }
            }

            const Array<TwoD, const NekDouble> & df   = geomFactors->GetDerivFactors(ptsKeys);
            const Array<OneD, const NekDouble> & jac  = geomFactors->GetJac(ptsKeys);

            LibUtilities::BasisKey tobasis0 = DetFaceBasisKey(face,0);
            LibUtilities::BasisKey tobasis1 = DetFaceBasisKey(face,1);

            // Number of quadrature points in face expansion.
            int nq_face = tobasis0.GetNumPoints()*tobasis1.GetNumPoints();

            int vCoordDim = GetCoordim();

            m_faceNormals[face] = Array<OneD, Array<OneD, NekDouble> >(vCoordDim);
            Array<OneD, Array<OneD, NekDouble> > &normal = m_faceNormals[face];
            for (i = 0; i < vCoordDim; ++i)
            {
                normal[i] = Array<OneD, NekDouble>(nq_face);
            }

            size_t nqb = nq_face;
            size_t nbnd= face;
            m_elmtBndNormDirElmtLen[nbnd] = Array<OneD, NekDouble> {nqb, 0.0};
            Array<OneD, NekDouble> &length = m_elmtBndNormDirElmtLen[nbnd];

            // Regular geometry case
            if((type == SpatialDomains::eRegular)||(type == SpatialDomains::eMovingRegular))
            {
                NekDouble fac;
                // Set up normals
                switch(face)
                {
                case 0:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        normal[i][0] = -df[3*i+2][0];
                    }
                    break;
                case 1:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        normal[i][0] = -df[3*i+1][0];
                    }
                    break;
                case 2:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        normal[i][0] = df[3*i][0];
                    }
                    break;
                case 3:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        normal[i][0] = df[3*i+1][0];
                    }
                    break;
                case 4:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        normal[i][0] = -df[3*i][0];
                    }
                    break;
                case 5:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        normal[i][0] = df[3*i+2][0];
                    }
                    break;
                default:
                    ASSERTL0(false,"face is out of range (edge < 5)");
                }

                // normalise
                fac = 0.0;
                for(i =0 ; i < vCoordDim; ++i)
                {
                    fac += normal[i][0]*normal[i][0];
                }
                fac = 1.0/sqrt(fac);

                Vmath::Fill(nqb, fac, length, 1);
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Fill(nq_face, fac*normal[i][0], normal[i], 1);
                }

            }
            else   // Set up deformed normals
            {
                int j, k;

                int nqe0 = ptsKeys[0].GetNumPoints();
                int nqe1 = ptsKeys[0].GetNumPoints();
                int nqe2 = ptsKeys[0].GetNumPoints();
                int nqe01 = nqe0*nqe1;
                int nqe02 = nqe0*nqe2;
                int nqe12 = nqe1*nqe2;

                int nqe;
                if (face == 0 || face == 5)
                {
                    nqe = nqe01;
                }
                else if (face == 1 || face == 3)
                {
                    nqe = nqe02;
                }
                else
                {
                    nqe = nqe12;
                }

                LibUtilities::PointsKey points0;
                LibUtilities::PointsKey points1;

                Array<OneD, NekDouble> faceJac(nqe);
                Array<OneD, NekDouble> normals(vCoordDim*nqe,0.0);

                // Extract Jacobian along face and recover local
                // derivates (dx/dr) for polynomial interpolation by
                // multiplying m_gmat by jacobian
                switch(face)
                {
                    case 0:
                        for(j = 0; j < nqe; ++j)
                        {
                            normals[j]       = -df[2][j]*jac[j];
                            normals[nqe+j]   = -df[5][j]*jac[j];
                            normals[2*nqe+j] = -df[8][j]*jac[j];
                            faceJac[j]       = jac[j];
                        }

                        points0 = ptsKeys[0];
                        points1 = ptsKeys[1];
                        break;
                    case 1:
                        for (j = 0; j < nqe0; ++j)
                        {
                            for(k = 0; k < nqe2; ++k)
                            {
                                int idx = j + nqe01*k;
                                normals[j+k*nqe0]       = -df[1][idx]*jac[idx];
                                normals[nqe+j+k*nqe0]   = -df[4][idx]*jac[idx];
                                normals[2*nqe+j+k*nqe0] = -df[7][idx]*jac[idx];
                                faceJac[j+k*nqe0]       = jac[idx];
                            }
                        }
                        points0 = ptsKeys[0];
                        points1 = ptsKeys[2];
                        break;
                    case 2:
                        for (j = 0; j < nqe1; ++j)
                        {
                            for(k = 0; k < nqe2; ++k)
                            {
                                int idx = nqe0-1+nqe0*j+nqe01*k;
                                normals[j+k*nqe1]       = df[0][idx]*jac[idx];
                                normals[nqe+j+k*nqe1]   = df[3][idx]*jac[idx];
                                normals[2*nqe+j+k*nqe1] = df[6][idx]*jac[idx];
                                faceJac[j+k*nqe1]       = jac[idx];
                            }
                        }
                        points0 = ptsKeys[1];
                        points1 = ptsKeys[2];
                        break;
                    case 3:
                        for (j = 0; j < nqe0; ++j)
                        {
                            for(k = 0; k < nqe2; ++k)
                            {
                                int idx = nqe0*(nqe1-1)+j+nqe01*k;
                                normals[j+k*nqe0]       = df[1][idx]*jac[idx];
                                normals[nqe+j+k*nqe0]   = df[4][idx]*jac[idx];
                                normals[2*nqe+j+k*nqe0] = df[7][idx]*jac[idx];
                                faceJac[j+k*nqe0]       = jac[idx];
                            }
                        }
                        points0 = ptsKeys[0];
                        points1 = ptsKeys[2];
                        break;
                    case 4:
                        for (j = 0; j < nqe1; ++j)
                        {
                            for(k = 0; k < nqe2; ++k)
                            {
                                int idx = j*nqe0+nqe01*k;
                                normals[j+k*nqe1]       = -df[0][idx]*jac[idx];
                                normals[nqe+j+k*nqe1]   = -df[3][idx]*jac[idx];
                                normals[2*nqe+j+k*nqe1] = -df[6][idx]*jac[idx];
                                faceJac[j+k*nqe1]       = jac[idx];
                            }
                        }
                        points0 = ptsKeys[1];
                        points1 = ptsKeys[2];
                        break;
                    case 5:
                        for (j = 0; j < nqe01; ++j)
                        {
                            int idx = j+nqe01*(nqe2-1);
                            normals[j]       = df[2][idx]*jac[idx];
                            normals[nqe+j]   = df[5][idx]*jac[idx];
                            normals[2*nqe+j] = df[8][idx]*jac[idx];
                            faceJac[j]       = jac[idx];
                        }
                        points0 = ptsKeys[0];
                        points1 = ptsKeys[1];
                        break;
                    default:
                    ASSERTL0(false,"face is out of range (face < 5)");
                }

                Array<OneD, NekDouble> work   (nq_face, 0.0);
                // Interpolate Jacobian and invert
                LibUtilities::Interp2D(points0, points1, faceJac,
                                       tobasis0.GetPointsKey(),
                                       tobasis1.GetPointsKey(),
                                       work);

                Vmath::Sdiv(nq_face,1.0,&work[0],1,&work[0],1);

                // interpolate
                for(i = 0; i < GetCoordim(); ++i)
                {
                    LibUtilities::Interp2D(points0, points1,
                                           &normals[i*nqe],
                                           tobasis0.GetPointsKey(),
                                           tobasis1.GetPointsKey(),
                                           &normal[i][0]);
                    Vmath::Vmul(nq_face,work,1,normal[i],1,normal[i],1);
                }

                //normalise normal vectors
                Vmath::Zero(nq_face,work,1);
                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vvtvp(nq_face,normal[i],1, normal[i],1,work,1,work,1);
                }

                Vmath::Vsqrt(nq_face,work,1,work,1);
                Vmath::Sdiv(nq_face,1.0,work,1,work,1);

                Vmath::Vcopy(nqb, work, 1, length, 1);

                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vmul(nq_face,normal[i],1,work,1,normal[i],1);
                }
            }
        }

        //-----------------------------
        // Operator creation functions
        //-----------------------------
        void HexExp::v_MassMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void HexExp::v_LaplacianMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            HexExp::v_LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void HexExp::v_LaplacianMatrixOp(
                const int k1,
                const int k2,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,
                                                        mkey);
        }

        void HexExp::v_WeakDerivMatrixOp(
                const int i,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
        }

        void HexExp::v_WeakDirectionalDerivMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::WeakDirectionalDerivMatrixOp_MatFree(inarray,
                                                                outarray,mkey);
        }

        void HexExp::v_MassLevelCurvatureMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::MassLevelCurvatureMatrixOp_MatFree(inarray,
                                                                outarray,mkey);
        }

        void HexExp::v_HelmholtzMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            HexExp::v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
        }


        void HexExp::v_GeneralMatrixOp_MatOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            //int nConsts = mkey.GetNconstants();
            DNekScalMatSharedPtr   mat = GetLocMatrix(mkey);

//            switch(nConsts)
//            {
//            case 0:
//                {
//                    mat = GetLocMatrix(mkey.GetMatrixType());
//                }
//                break;
//            case 1:
//                {
//                    mat = GetLocMatrix(mkey.GetMatrixType(),mkey.GetConstant(0));
//                }
//                break;
//            case 2:
//                {
//                    mat = GetLocMatrix(mkey.GetMatrixType(),mkey.GetConstant(0),mkey.GetConstant(1));
//                }
//                break;
//
//            default:
//                {
//                    NEKERROR(ErrorUtil::efatal, "Unknown number of constants");
//                }
//                break;
//            }

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

        /**
         * This function is used to compute exactly the advective numerical flux
         * on the interface of two elements with different expansions, hence an
         * appropriate number of Gauss points has to be used. The number of
         * Gauss points has to be equal to the number used by the highest
         * polynomial degree of the two adjacent elements
         *
         * @param   numMin     Is the reduced polynomial order
         * @param   inarray    Input array of coefficients
         * @param   dumpVar    Output array of reduced coefficients.
         */
        void HexExp::v_ReduceOrderCoeffs(
            int                                 numMin,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            int n_coeffs = inarray.size();
            int nmodes0  = m_base[0]->GetNumModes();
            int nmodes1  = m_base[1]->GetNumModes();
            int nmodes2  = m_base[2]->GetNumModes();
            int numMax   = nmodes0;

            Array<OneD, NekDouble> coeff     (n_coeffs);
            Array<OneD, NekDouble> coeff_tmp1(nmodes0*nmodes1, 0.0);
            Array<OneD, NekDouble> coeff_tmp2(n_coeffs,        0.0);
            Array<OneD, NekDouble> tmp, tmp2, tmp3, tmp4;

            Vmath::Vcopy(n_coeffs,inarray,1,coeff_tmp2,1);

            const LibUtilities::PointsKey Pkey0(
                nmodes0, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey Pkey1(
                nmodes1, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey Pkey2(
                nmodes2, LibUtilities::eGaussLobattoLegendre);

            LibUtilities::BasisKey b0(
                m_base[0]->GetBasisType(), nmodes0, Pkey0);
            LibUtilities::BasisKey b1(
                m_base[1]->GetBasisType(), nmodes1, Pkey1);
            LibUtilities::BasisKey b2(
                m_base[2]->GetBasisType(), nmodes2, Pkey2);
            LibUtilities::BasisKey bortho0(
                LibUtilities::eOrtho_A,    nmodes0, Pkey0);
            LibUtilities::BasisKey bortho1(
                LibUtilities::eOrtho_A,    nmodes1, Pkey1);
            LibUtilities::BasisKey bortho2(
                LibUtilities::eOrtho_A,    nmodes2, Pkey2);

            LibUtilities::InterpCoeff3D(
                b0,      b1,      b2,      coeff_tmp2,
                bortho0, bortho1, bortho2, coeff);

            Vmath::Zero(n_coeffs, coeff_tmp2, 1);

            int cnt = 0, cnt2 = 0;

            for (int u = 0; u < numMin+1; ++u)
            {
                for (int i = 0; i < numMin; ++i)
                {
                    Vmath::Vcopy(numMin,
                                 tmp  = coeff+cnt+cnt2,1,
                                 tmp2 = coeff_tmp1+cnt,1);

                    cnt = i*numMax;
                }

                Vmath::Vcopy(nmodes0*nmodes1,
                             tmp3 = coeff_tmp1,1,
                             tmp4 = coeff_tmp2+cnt2,1);

                cnt2 = u*nmodes0*nmodes1;
            }

            LibUtilities::InterpCoeff3D(
                bortho0, bortho1, bortho2, coeff_tmp2,
                b0,      b1,      b2,      outarray);
        }

        void HexExp::v_SVVLaplacianFilter(
                    Array<OneD, NekDouble> &array,
                    const StdRegions::StdMatrixKey &mkey)
        {
            int nq = GetTotPoints();

            // Calculate sqrt of the Jacobian
            Array<OneD, const NekDouble> jac =
                                    m_metricinfo->GetJac(GetPointsKeys());
            Array<OneD, NekDouble> sqrt_jac(nq);
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vsqrt(nq,jac,1,sqrt_jac,1);
            }
            else
            {
                Vmath::Fill(nq,sqrt(jac[0]),sqrt_jac,1);
            }

            // Multiply array by sqrt(Jac)
            Vmath::Vmul(nq,sqrt_jac,1,array,1,array,1);

            // Apply std region filter
            StdHexExp::v_SVVLaplacianFilter( array, mkey);

            // Divide by sqrt(Jac)
            Vmath::Vdiv(nq,array,1,sqrt_jac,1,array,1);
        }

        //-----------------------------
        // Matrix creation functions
        //-----------------------------
        DNekMatSharedPtr HexExp::v_GenMatrix(
               const StdRegions::StdMatrixKey &mkey)
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
            case StdRegions::eInvLaplacianWithUnityMean:
                returnval = Expansion3D::v_GenMatrix(mkey);
                break;
            default:
                returnval = StdHexExp::v_GenMatrix(mkey);
            }

            return returnval;
        }


        DNekMatSharedPtr HexExp::v_CreateStdMatrix(
                const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();

            StdRegions::StdHexExpSharedPtr tmp = MemoryManager<StdHexExp>
                                        ::AllocateSharedPtr(bkey0,bkey1,bkey2);

            return tmp->GetStdMatrix(mkey);
        }


        DNekScalMatSharedPtr HexExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,
                     "Geometric information is not set up");

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed ||
                            mkey.GetNVarCoeff())
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat
                                        = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(jac,mat);
                    }
                }
                break;
            case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,
                                                    DetShapeType(), *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble fac = 1.0/(m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat
                                        = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(fac,mat);
                    }
                }
                break;
            case StdRegions::eWeakDeriv0:
            case StdRegions::eWeakDeriv1:
            case StdRegions::eWeakDeriv2:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed ||
                            mkey.GetNVarCoeff())
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        Array<TwoD, const NekDouble> df
                                    = m_metricinfo->GetDerivFactors(ptsKeys);
                        int dir = 0;

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
                            default:
                                break;
                        }

                        MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                            mkey.GetShapeType(), *this);
                        MatrixKey deriv1key(StdRegions::eWeakDeriv1,
                                            mkey.GetShapeType(), *this);
                        MatrixKey deriv2key(StdRegions::eWeakDeriv2,
                                            mkey.GetShapeType(), *this);

                        DNekMat &deriv0 = *GetStdMatrix(deriv0key);
                        DNekMat &deriv1 = *GetStdMatrix(deriv1key);
                        DNekMat &deriv2 = *GetStdMatrix(deriv2key);

                        int rows = deriv0.GetRows();
                        int cols = deriv1.GetColumns();

                        DNekMatSharedPtr WeakDeriv = MemoryManager<DNekMat>
                                                ::AllocateSharedPtr(rows,cols);

                        (*WeakDeriv) = df[3*dir  ][0]*deriv0
                                     + df[3*dir+1][0]*deriv1
                                     + df[3*dir+2][0]*deriv2;

                        returnval = MemoryManager<DNekScalMat>
                                            ::AllocateSharedPtr(jac,WeakDeriv);
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed ||
                        mkey.GetNVarCoeff()||
                        mkey.ConstFactorExists(
                                StdRegions::eFactorSVVCutoffRatio))
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetShapeType(), *this);
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetShapeType(), *this);
                        MatrixKey lap02key(StdRegions::eLaplacian02,
                                           mkey.GetShapeType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetShapeType(), *this);
                        MatrixKey lap12key(StdRegions::eLaplacian12,
                                           mkey.GetShapeType(), *this);
                        MatrixKey lap22key(StdRegions::eLaplacian22,
                                           mkey.GetShapeType(), *this);

                        DNekMat &lap00 = *GetStdMatrix(lap00key);
                        DNekMat &lap01 = *GetStdMatrix(lap01key);
                        DNekMat &lap02 = *GetStdMatrix(lap02key);
                        DNekMat &lap11 = *GetStdMatrix(lap11key);
                        DNekMat &lap12 = *GetStdMatrix(lap12key);
                        DNekMat &lap22 = *GetStdMatrix(lap22key);

                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        Array<TwoD, const NekDouble> gmat
                                            = m_metricinfo->GetGmat(ptsKeys);

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>
                                                ::AllocateSharedPtr(rows,cols);

                        (*lap)  = gmat[0][0]*lap00
                                + gmat[4][0]*lap11
                                + gmat[8][0]*lap22
                                + gmat[3][0]*(lap01 + Transpose(lap01))
                                + gmat[6][0]*(lap02 + Transpose(lap02))
                                + gmat[7][0]*(lap12 + Transpose(lap12));

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(jac,lap);
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetShapeType(), *this);
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + lambda*MassMat;

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
                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,mat);
                    }
                }
                break;
            case StdRegions::eHybridDGHelmholtz:
            case StdRegions::eHybridDGLamToU:
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
            case StdRegions::eHybridDGHelmBndLam:
            case StdRegions::eInvLaplacianWithUnityMean:
                {
                    NekDouble one    = 1.0;

                    DNekMatSharedPtr mat = GenMatrix(mkey);
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            case StdRegions::eInvHybridDGHelmholtz:
                {
                    NekDouble one = 1.0;

//                    StdRegions::StdMatrixKey hkey(StdRegions::eHybridDGHelmholtz,
//                                                  DetShapeType(),*this,
//                                                  mkey.GetConstant(0),
//                                                  mkey.GetConstant(1));
                    MatrixKey hkey(StdRegions::eHybridDGHelmholtz, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            case StdRegions::ePreconLinearSpace:
                {
                    NekDouble one = 1.0;
                    MatrixKey helmkey(StdRegions::eHelmholtz, mkey.GetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalBlkMatSharedPtr helmStatCond = GetLocStaticCondMatrix(helmkey);
                    DNekScalMatSharedPtr A =helmStatCond->GetBlock(0,0);
                    DNekMatSharedPtr R=BuildVertexMatrix(A);

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,R);
                }
                break;
            case StdRegions::ePreconLinearSpaceMass:
                {
                    NekDouble one = 1.0;
                    MatrixKey masskey(StdRegions::eMass, mkey.GetShapeType(), *this);
                    DNekScalBlkMatSharedPtr massStatCond = GetLocStaticCondMatrix(masskey);
                    DNekScalMatSharedPtr A =massStatCond->GetBlock(0,0);
                    DNekMatSharedPtr R=BuildVertexMatrix(A);

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,R);
                }
                break;
            case StdRegions::ePreconR:
                {
                    NekDouble one = 1.0;
                    MatrixKey helmkey(StdRegions::eHelmholtz, mkey.GetShapeType(), *this,mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalBlkMatSharedPtr helmStatCond = GetLocStaticCondMatrix(helmkey);
                    DNekScalMatSharedPtr A =helmStatCond->GetBlock(0,0);

                    DNekScalMatSharedPtr Atmp;
                    DNekMatSharedPtr R=BuildTransformationMatrix(A,mkey.GetMatrixType());

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,R);
                }
                break;
            case StdRegions::ePreconRMass:
                {
                    NekDouble one = 1.0;
                    MatrixKey masskey(StdRegions::eMass, mkey.GetShapeType(), *this);
                    DNekScalBlkMatSharedPtr massStatCond = GetLocStaticCondMatrix(masskey);
                    DNekScalMatSharedPtr A =massStatCond->GetBlock(0,0);

                    DNekScalMatSharedPtr Atmp;
                    DNekMatSharedPtr R=BuildTransformationMatrix(A,mkey.GetMatrixType());

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,R);
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


        DNekScalBlkMatSharedPtr HexExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            unsigned int nbdry = NumBndryCoeffs();
            unsigned int nint = (unsigned int)(m_ncoeffs - nbdry);
            unsigned int exp_size[] = {nbdry,nint};
            unsigned int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks,nblks,exp_size,exp_size); //Really need a constructor which takes Arrays
            NekDouble factor = 1.0;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eLaplacian:
            case StdRegions::eHelmholtz: // special case since Helmholtz not defined in StdRegions

                // use Deformed case for both regular and deformed geometries
                factor = 1.0;
                goto UseLocRegionsMatrix;
                break;
            default:
                if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed ||
                        mkey.GetNVarCoeff())
                {
                    factor = 1.0;
                    goto UseLocRegionsMatrix;
                }
                else
                {
                    DNekScalMatSharedPtr mat = GetLocMatrix(mkey);
                    factor = mat->Scale();
                    goto UseStdRegionsMatrix;
                }
                break;
            UseStdRegionsMatrix:
                {
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekBlkMatSharedPtr  mat = GetStdStaticCondMatrix(mkey);
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


        DNekScalMatSharedPtr HexExp::v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }


        DNekScalBlkMatSharedPtr HexExp::v_GetLocStaticCondMatrix(
                const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
        }

        void HexExp::v_DropLocStaticCondMatrix(const MatrixKey &mkey)
        {
            m_staticCondMatrixManager.DeleteObject(mkey);
        }

        void HexExp::v_LaplacianMatrixOp_MatFree_Kernel(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                      Array<OneD,       NekDouble> &wsp)
        {
            // This implementation is only valid when there are no
            // coefficients associated to the Laplacian operator
            if (m_metrics.count(eMetricLaplacian00) == 0)
            {
                ComputeLaplacianMetric();
            }

            int       nquad0  = m_base[0]->GetNumPoints();
            int       nquad1  = m_base[1]->GetNumPoints();
            int       nquad2  = m_base[2]->GetNumPoints();
            int       nqtot   = nquad0*nquad1*nquad2;

            ASSERTL1(wsp.size() >= 6*nqtot,
                     "Insufficient workspace size.");

            const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble>& base2  = m_base[2]->GetBdata();
            const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase2 = m_base[2]->GetDbdata();
            const Array<OneD, const NekDouble>& metric00 = m_metrics[eMetricLaplacian00];
            const Array<OneD, const NekDouble>& metric01 = m_metrics[eMetricLaplacian01];
            const Array<OneD, const NekDouble>& metric02 = m_metrics[eMetricLaplacian02];
            const Array<OneD, const NekDouble>& metric11 = m_metrics[eMetricLaplacian11];
            const Array<OneD, const NekDouble>& metric12 = m_metrics[eMetricLaplacian12];
            const Array<OneD, const NekDouble>& metric22 = m_metrics[eMetricLaplacian22];

            // Allocate temporary storage
            Array<OneD,NekDouble> wsp0(wsp);
            Array<OneD,NekDouble> wsp1(wsp+1*nqtot);
            Array<OneD,NekDouble> wsp2(wsp+2*nqtot);
            Array<OneD,NekDouble> wsp3(wsp+3*nqtot);
            Array<OneD,NekDouble> wsp4(wsp+4*nqtot);
            Array<OneD,NekDouble> wsp5(wsp+5*nqtot);

            StdExpansion3D::PhysTensorDeriv(inarray,wsp0,wsp1,wsp2);

            // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
            // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
            // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
            // especially for this purpose
            Vmath::Vvtvvtp(nqtot,&metric00[0],1,&wsp0[0],1,&metric01[0],1,&wsp1[0],1,&wsp3[0],1);
            Vmath::Vvtvp  (nqtot,&metric02[0],1,&wsp2[0],1,&wsp3[0],1,&wsp3[0],1);
            Vmath::Vvtvvtp(nqtot,&metric01[0],1,&wsp0[0],1,&metric11[0],1,&wsp1[0],1,&wsp4[0],1);
            Vmath::Vvtvp  (nqtot,&metric12[0],1,&wsp2[0],1,&wsp4[0],1,&wsp4[0],1);
            Vmath::Vvtvvtp(nqtot,&metric02[0],1,&wsp0[0],1,&metric12[0],1,&wsp1[0],1,&wsp5[0],1);
            Vmath::Vvtvp  (nqtot,&metric22[0],1,&wsp2[0],1,&wsp5[0],1,&wsp5[0],1);

            // outarray = m = (D_xi1 * B)^T * k
            // wsp1     = n = (D_xi2 * B)^T * l
            IProductWRTBase_SumFacKernel(dbase0,base1,base2,wsp3,outarray,wsp0,false,true,true);
            IProductWRTBase_SumFacKernel(base0,dbase1,base2,wsp4,wsp2,    wsp0,true,false,true);
            Vmath::Vadd(m_ncoeffs,wsp2.get(),1,outarray.get(),1,outarray.get(),1);
            IProductWRTBase_SumFacKernel(base0,base1,dbase2,wsp5,wsp2,    wsp0,true,true,false);
            Vmath::Vadd(m_ncoeffs,wsp2.get(),1,outarray.get(),1,outarray.get(),1);
        }


        void HexExp::v_ComputeLaplacianMetric()
        {
            if (m_metrics.count(eMetricQuadrature) == 0)
            {
                ComputeQuadratureMetric();
            }

            const SpatialDomains::GeomType type = m_metricinfo->GetGtype();
            const unsigned int nqtot = GetTotPoints();
            const unsigned int dim = 3;
            const MetricType m[3][3] = { {eMetricLaplacian00, eMetricLaplacian01, eMetricLaplacian02},
                                       {eMetricLaplacian01, eMetricLaplacian11, eMetricLaplacian12},
                                       {eMetricLaplacian02, eMetricLaplacian12, eMetricLaplacian22}
            };

            for (unsigned int i = 0; i < dim; ++i)
            {
                for (unsigned int j = i; j < dim; ++j)
                {
                    m_metrics[m[i][j]] = Array<OneD, NekDouble>(nqtot);
                    const Array<TwoD, const NekDouble> &gmat =
                                 m_metricinfo->GetGmat(GetPointsKeys());
                    if (type == SpatialDomains::eDeformed)
                    {
                        Vmath::Vcopy(nqtot, &gmat[i*dim+j][0], 1,
                                            &m_metrics[m[i][j]][0], 1);
                    }
                    else
                    {
                        Vmath::Fill(nqtot, gmat[i*dim+j][0],
                                    &m_metrics[m[i][j]][0], 1);
                    }
                    MultiplyByQuadratureMetric(m_metrics[m[i][j]],
                                               m_metrics[m[i][j]]);

                }
            }
        }

    }//end of namespace
}//end of namespace
