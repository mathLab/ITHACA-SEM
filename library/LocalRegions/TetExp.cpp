///////////////////////////////////////////////////////////////////////////////
//
// File TetExp.cpp
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

#include <LocalRegions/TetExp.h>
#include <SpatialDomains/SegGeom.h>

#include <LibUtilities/Foundations/Interp.h>

namespace Nektar
{
    namespace LocalRegions
    {
        /**
         * @class TetExp
         * Defines a Tetrahedral local expansion.
         */

        /**
	 * \brief Constructor using BasisKey class for quadrature points and 
	 * order definition 
	 *
         * @param   Ba          Basis key for first coordinate.
         * @param   Bb          Basis key for second coordinate.
         * @param   Bc          Basis key for third coordinate.
         */
        TetExp::TetExp( const LibUtilities::BasisKey &Ba,
                        const LibUtilities::BasisKey &Bb,
                        const LibUtilities::BasisKey &Bc,
                        const SpatialDomains::TetGeomSharedPtr &geom
                        ):
            StdExpansion  (StdRegions::StdTetData::getNumberOfCoefficients(Ba.GetNumModes(),Bb.GetNumModes(),Bc.GetNumModes()),3,Ba,Bb,Bc),
            StdExpansion3D(StdRegions::StdTetData::getNumberOfCoefficients(Ba.GetNumModes(),Bb.GetNumModes(),Bc.GetNumModes()),Ba,Bb,Bc),
            StdRegions::StdTetExp(Ba,Bb,Bc),
            Expansion     (),
            Expansion3D   (),
            m_geom(geom),
            m_metricinfo(m_geom->GetGeomFactors(m_base)),
            m_matrixManager(
                    boost::bind(&TetExp::CreateMatrix, this, _1),
                    std::string("TetExpMatrix")),
            m_staticCondMatrixManager(
                    boost::bind(&TetExp::CreateStaticCondMatrix, this, _1),
                    std::string("TetExpStaticCondMatrix"))
        {
        }


        /**
	 * \brief Copy Constructor
	 */
        TetExp::TetExp(const TetExp &T):
            StdExpansion(T),
            StdExpansion3D(T),
            StdRegions::StdTetExp(T),
            Expansion(T),
            Expansion3D(T),
            m_geom(T.m_geom),
            m_metricinfo(T.m_metricinfo),
            m_matrixManager(T.m_matrixManager),
            m_staticCondMatrixManager(T.m_staticCondMatrixManager)
        {
        }


        /**
	 * \brief Destructor
	 */
        TetExp::~TetExp()
        {
        }


        //-----------------------------
        // Integration Methods
        //-----------------------------
        /**
         * \brief Integrate the physical point list \a inarray over region
         *
         * @param   inarray     Definition of function to be returned at
         *                      quadrature point of expansion.
         * @returns \f$\int^1_{-1}\int^1_{-1} \int^1_{-1}
         *   u(\eta_1, \eta_2, \eta_3) J[i,j,k] d \eta_1 d \eta_2 d \eta_3 \f$
         * where \f$inarray[i,j,k] = u(\eta_{1i},\eta_{2j},\eta_{3k})
         * \f$ and \f$ J[i,j,k] \f$ is the Jacobian evaluated at the quadrature
         * point.
         */
        NekDouble TetExp::v_Integral(
                  const Array<OneD, const NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            NekDouble retrunVal;
            Array<OneD,NekDouble> tmp(nquad0*nquad1*nquad2);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,
                            (NekDouble*)&inarray[0],1, &tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,(NekDouble) jac[0],
                            (NekDouble*)&inarray[0],1,&tmp[0],1);
            }

            // call StdTetExp version;
            retrunVal = StdTetExp::v_Integral(tmp);

            return retrunVal;
        }
        

        //-----------------------------
        // Differentiation Methods
        //-----------------------------
        /**
	 * \brief Differentiate \a inarray in the three coordinate directions.
	 * 
         * @param   inarray     Input array of values at quadrature points to
         *                      be differentiated.
         * @param   out_d0      Derivative in first coordinate direction.
         * @param   out_d1      Derivative in second coordinate direction.
         * @param   out_d2      Derivative in third coordinate direction.
         */
        void TetExp::v_PhysDeriv(
                 const Array<OneD, const NekDouble> & inarray,
                       Array<OneD,NekDouble> &out_d0,
                       Array<OneD,NekDouble> &out_d1,
                       Array<OneD,NekDouble> &out_d2)
        {
            int  TotPts = m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints()*
                m_base[2]->GetNumPoints();
            
            Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> Diff0 = Array<OneD,NekDouble>(3*TotPts);
            Array<OneD,NekDouble> Diff1 = Diff0 + TotPts;
            Array<OneD,NekDouble> Diff2 = Diff1 + TotPts;
            
            StdTetExp::v_PhysDeriv(inarray, Diff0, Diff1, Diff2);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul  (TotPts,&gmat[0][0],1,&Diff0[0],1, &out_d0[0], 1);
                    Vmath::Vvtvp (TotPts,&gmat[1][0],1,&Diff1[0],1, &out_d0[0], 1,&out_d0[0],1);
                    Vmath::Vvtvp (TotPts,&gmat[2][0],1,&Diff2[0],1, &out_d0[0], 1,&out_d0[0],1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Vmul  (TotPts,&gmat[3][0],1,&Diff0[0],1, &out_d1[0], 1);
                    Vmath::Vvtvp (TotPts,&gmat[4][0],1,&Diff1[0],1, &out_d1[0], 1,&out_d1[0],1);
                    Vmath::Vvtvp (TotPts,&gmat[5][0],1,&Diff2[0],1, &out_d1[0], 1,&out_d1[0],1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Vmul  (TotPts,&gmat[6][0],1,&Diff0[0],1, &out_d2[0], 1);
                    Vmath::Vvtvp (TotPts,&gmat[7][0],1,&Diff1[0],1, &out_d2[0], 1, &out_d2[0],1);
                    Vmath::Vvtvp (TotPts,&gmat[8][0],1,&Diff2[0],1, &out_d2[0], 1, &out_d2[0],1);
                }
            }
            else // regular geometry
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul (TotPts,gmat[0][0],&Diff0[0],1, &out_d0[0], 1);
                    Blas::Daxpy (TotPts,gmat[1][0],&Diff1[0],1, &out_d0[0], 1);
                    Blas::Daxpy (TotPts,gmat[2][0],&Diff2[0],1, &out_d0[0], 1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Smul (TotPts,gmat[3][0],&Diff0[0],1, &out_d1[0], 1);
                    Blas::Daxpy (TotPts,gmat[4][0],&Diff1[0],1, &out_d1[0], 1);
                    Blas::Daxpy (TotPts,gmat[5][0],&Diff2[0],1, &out_d1[0], 1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Smul (TotPts,gmat[6][0],&Diff0[0],1, &out_d2[0], 1);
                    Blas::Daxpy (TotPts,gmat[7][0],&Diff1[0],1, &out_d2[0], 1);
                    Blas::Daxpy (TotPts,gmat[8][0],&Diff2[0],1, &out_d2[0], 1);
                }
            }
        }


        //-----------------------------
        // Transforms
        //-----------------------------
        /**
	 * \brief Forward transform from physical quadrature space stored in
	 * \a inarray and evaluate the expansion coefficients and store
	 * in \a (this)->_coeffs
	 *
         * @param   inarray     Array of physical quadrature points to be
         *                      transformed.
         * @param   outarray    Array of coefficients to update.
         */
        void TetExp::v_FwdTrans( 
                 const Array<OneD, const NekDouble> & inarray,
                       Array<OneD,NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation())&&(m_base[2]->Collocation()))
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&m_coeffs[0],1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass,
                                              DetExpansionType(),*this);
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
	 * basis B=m_base0*m_base1*m_base2 and put into outarray:
	 *
         * \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta}
         *   & = & \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
         *     \psi_{p}^{a} (\eta_{1i}) \psi_{pq}^{b} (\eta_{2j}) \psi_{pqr}^{c}
         *     (\eta_{3k}) w_i w_j w_k u(\eta_{1,i} \eta_{2,j} \eta_{3,k})
         * J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\eta_{1,i})
         *   \sum_{j=0}^{nq_1} \psi_{pq}^b(\eta_{2,j}) \sum_{k=0}^{nq_2}
         *   \psi_{pqr}^c u(\eta_{1i},\eta_{2j},\eta_{3k}) J_{i,j,k}
         * \end{array} \f$ \n
         * where
         * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3)
         *   = \psi_p^a (\eta_1) \psi_{pq}^b (\eta_2) \psi_{pqr}^c (\eta_3) \f$
         * which can be implemented as \n
         * \f$f_{pqr} (\xi_{3k})
         *   = \sum_{k=0}^{nq_3} \psi_{pqr}^c u(\eta_{1i},\eta_{2j},\eta_{3k})
         * J_{i,j,k} = {\bf B_3 U}   \f$ \n
         * \f$ g_{pq} (\xi_{3k})
         *   = \sum_{j=0}^{nq_1} \psi_{pq}^b (\xi_{2j}) f_{pqr} (\xi_{3k})
         *   = {\bf B_2 F}  \f$ \n
         * \f$ (\phi_{pqr}, u)_{\delta}
         *   = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{pq} (\xi_{3k})
         *   = {\bf B_1 G} \f$
         */
        void TetExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            v_IProductWRTBase_SumFac(inarray, outarray);
        }

        void TetExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            const int nquad0 = m_base[0]->GetNumPoints();
            const int nquad1 = m_base[1]->GetNumPoints();
            const int nquad2 = m_base[2]->GetNumPoints();
            const int order0 = m_base[0]->GetNumModes();
            const int order1 = m_base[1]->GetNumModes();
            Array<OneD, NekDouble> wsp(nquad1*nquad2*order0 +
                                       nquad2*order0*(order1+1)/2);
            Array<OneD, NekDouble> tmp(nquad0*nquad1*nquad2);

            MultiplyByQuadratureMetric(inarray, tmp);
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetBdata(),
                                         tmp,outarray,wsp,
                                         true,true,true);
        }

        /**
         * @brief Calculates the inner product \f$ I_{pqr} = (u,
         * \partial_{x_i} \phi_{pqr}) \f$.
         * 
         * The derivative of the basis functions is performed using the chain
         * rule in order to incorporate the geometric factors. Assuming that
         * the basis functions are a tensor product
         * \f$\phi_{pqr}(\eta_1,\eta_2,\eta_3) =
         * \phi_1(\eta_1)\phi_2(\eta_2)\phi_3(\eta_3)\f$, this yields the
         * result
         * 
         * \f[
         * I_{pqr} = \sum_{j=1}^3 \left(u, \frac{\partial u}{\partial \eta_j}
         * \frac{\partial \eta_j}{\partial x_i}\right)
         * \f]
         * 
         * In the prismatic element, we must also incorporate a second set of
         * geometric factors which incorporate the collapsed co-ordinate
         * system, so that
         * 
         * \f[ \frac{\partial\eta_j}{\partial x_i} = \sum_{k=1}^3
         * \frac{\partial\eta_j}{\partial\xi_k}\frac{\partial\xi_k}{\partial
         * x_i} \f]
         * 
         * These derivatives can be found on p152 of Sherwin & Karniadakis.
         * 
         * @param dir       Direction in which to take the derivative.
         * @param inarray   The function \f$ u \f$.
         * @param outarray  Value of the inner product.
         */
        void TetExp::v_IProductWRTDerivBase(
            const int                           dir, 
            const Array<OneD, const NekDouble> &inarray, 
                  Array<OneD,       NekDouble> &outarray)
        {
            const int nquad0 = m_base[0]->GetNumPoints();
            const int nquad1 = m_base[1]->GetNumPoints();
            const int nquad2 = m_base[2]->GetNumPoints();
            const int order0 = m_base[0]->GetNumModes ();
            const int order1 = m_base[1]->GetNumModes ();
            const int nqtot  = nquad0*nquad1*nquad2;
            int i, j;

            const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble> &z1 = m_base[1]->GetZ();
            const Array<OneD, const NekDouble> &z2 = m_base[2]->GetZ();
            
            Array<OneD, NekDouble> h0   (nqtot);
            Array<OneD, NekDouble> h1   (nqtot);
            Array<OneD, NekDouble> h2   (nqtot);
            Array<OneD, NekDouble> h3   (nqtot);
            Array<OneD, NekDouble> tmp1 (nqtot);
            Array<OneD, NekDouble> tmp2 (nqtot);
            Array<OneD, NekDouble> tmp3 (nqtot);
            Array<OneD, NekDouble> tmp4 (nqtot);
            Array<OneD, NekDouble> tmp5 (nqtot);
            Array<OneD, NekDouble> tmp6 (m_ncoeffs);
            Array<OneD, NekDouble> wsp  (nquad1*nquad2*order0 +
                                         nquad2*order0*(order1+1)/2);
            
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            MultiplyByQuadratureMetric(inarray,tmp1);
            
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nqtot,&gmat[3*dir][0],  1,tmp1.get(),1,tmp2.get(),1);
                Vmath::Vmul(nqtot,&gmat[3*dir+1][0],1,tmp1.get(),1,tmp3.get(),1);
                Vmath::Vmul(nqtot,&gmat[3*dir+2][0],1,tmp1.get(),1,tmp4.get(),1);
            }
            else
            {
                Vmath::Smul(nqtot, gmat[3*dir  ][0],tmp1.get(),1,tmp2.get(), 1);
                Vmath::Smul(nqtot, gmat[3*dir+1][0],tmp1.get(),1,tmp3.get(), 1);
                Vmath::Smul(nqtot, gmat[3*dir+2][0],tmp1.get(),1,tmp4.get(), 1);
            }
            
            const int nq01 = nquad0*nquad1;
            const int nq12 = nquad1*nquad2;

            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Vmath::Fill(nquad0, 4.0/(1.0-z1[i])/(1.0-z2[j]),
                                &h0[0]+i*nquad0 + j*nq01,1);
                    Vmath::Fill(nquad0, 2.0/(1.0-z1[i])/(1.0-z2[j]),
                                &h1[0]+i*nquad0 + j*nq01,1);
                    Vmath::Fill(nquad0, 2.0/(1.0-z2[j]),
                                &h2[0]+i*nquad0 + j*nq01,1);
                    Vmath::Fill(nquad0, (1.0+z1[i])/(1.0-z2[j]),
                                &h3[0]+i*nquad0 + j*nq01,1);
                }
            }

            for(i = 0; i < nquad0; i++)
            {
                Blas::Dscal(nq12, 1+z0[i], &h1[0]+i, nquad0);
            }
            
            // Assemble terms for first IP.
            Vmath::Vvtvvtp(nqtot, &tmp2[0], 1, &h0[0], 1,
                                  &tmp3[0], 1, &h1[0], 1,
                                  &tmp5[0], 1);
            Vmath::Vvtvp  (nqtot, &tmp4[0], 1, &h1[0], 1,
                                  &tmp5[0], 1, &tmp5[0], 1);

            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                         m_base[1]->GetBdata (),
                                         m_base[2]->GetBdata (),
                                         tmp5,outarray,wsp,
                                         true,true,true);

            // Assemble terms for second IP.
            Vmath::Vvtvvtp(nqtot, &tmp3[0], 1, &h2[0], 1,
                                  &tmp4[0], 1, &h3[0], 1,
                                  &tmp5[0], 1);
            
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata (),
                                         m_base[1]->GetDbdata(),
                                         m_base[2]->GetBdata (),
                                         tmp5,tmp6,wsp,
                                         true,true,true);
            Vmath::Vadd(m_ncoeffs, tmp6, 1, outarray, 1, outarray, 1);

            // Do third IP.
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata (),
                                         m_base[1]->GetBdata (),
                                         m_base[2]->GetDbdata(),
                                         tmp4,tmp6,wsp,
                                         true,true,true);

            // Sum contributions.
            Vmath::Vadd(m_ncoeffs, tmp6, 1, outarray, 1, outarray, 1);
        }


        //-----------------------------
        // Evaluation functions
        //-----------------------------
        /**
         * @param   coord       Physical space coordinate
         * @returns Evaluation of expansion at given coordinate.
         */
        NekDouble TetExp::v_PhysEvaluate(
                  const Array<OneD, const NekDouble> &coord)
        {
            return PhysEvaluate(coord,m_phys);
        }


        /**
         * @param   coord       Physical space coordinate
         * @returns Evaluation of expansion at given coordinate.
         */
        NekDouble TetExp::v_PhysEvaluate(
                  const Array<OneD, const NekDouble> &coord,
                  const Array<OneD, const NekDouble> & physvals)
        {
            ASSERTL0(m_geom,"m_geom not defined");

            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(3);

            // Get the local (eta) coordinates of the point
            m_geom->GetLocCoords(coord,Lcoord);

            // Evaluate point in local (eta) coordinates.
            return StdExpansion3D::v_PhysEvaluate(Lcoord,physvals);
        }


        /** 
	 * \brief Get the x,y,z coordinates of each quadrature point.
	 */
        void TetExp::v_GetCoords(
                  Array<OneD,NekDouble> &coords_0,
                  Array<OneD,NekDouble> &coords_1,
                  Array<OneD,NekDouble> &coords_2)
        {
            LibUtilities::BasisSharedPtr CBasis0;
            LibUtilities::BasisSharedPtr CBasis1;
            LibUtilities::BasisSharedPtr CBasis2;
            Array<OneD,NekDouble>  x;

            ASSERTL0(m_geom, "m_geom not define");

            // get physical points defined in Geom
            m_geom->FillGeom();  //TODO: implement

            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements(), "output coords_2 is not defined");
                CBasis0 = m_geom->GetBasis(2,0);
                CBasis1 = m_geom->GetBasis(2,1);
                CBasis2 = m_geom->GetBasis(2,2);

                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey()))&&
                   (m_base[2]->GetBasisKey().SamePoints(CBasis2->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(2);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints()*m_base[2]->GetNumPoints(),
                                x, 1, coords_2, 1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp3D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), CBasis2->GetPointsKey(), &(m_geom->UpdatePhys(2))[0],
                             m_base[0]->GetPointsKey(), m_base[1]->GetPointsKey(), m_base[2]->GetPointsKey(), &coords_2[0]);
                }
            case 2:
                ASSERTL0(coords_1.num_elements(), "output coords_1 is not defined");

                CBasis0 = m_geom->GetBasis(1,0);
                CBasis1 = m_geom->GetBasis(1,1);
                CBasis2 = m_geom->GetBasis(1,2);

                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey()))&&
                   (m_base[2]->GetBasisKey().SamePoints(CBasis2->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(1);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints()*m_base[2]->GetNumPoints(),
                                x, 1, coords_1, 1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp3D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), CBasis2->GetPointsKey(), &(m_geom->UpdatePhys(1))[0],
                             m_base[0]->GetPointsKey(), m_base[1]->GetPointsKey(), m_base[2]->GetPointsKey(), &coords_1[0]);
                }
            case 1:
                ASSERTL0(coords_0.num_elements(), "output coords_0 is not defined");

                CBasis0 = m_geom->GetBasis(0,0);
                CBasis1 = m_geom->GetBasis(0,1);
                CBasis2 = m_geom->GetBasis(0,2);

                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey()))&&
                   (m_base[2]->GetBasisKey().SamePoints(CBasis2->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(0);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints()*m_base[2]->GetNumPoints(),
                                x, 1, coords_0, 1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp3D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), CBasis2->GetPointsKey(), &(m_geom->UpdatePhys(0))[0],
                             m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),m_base[2]->GetPointsKey(),&coords_0[0]);
                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions are greater than 3");
                break;
            }
        }


        /**
	 * \brief Get the coordinates "coords" at the local coordinates "Lcoords"
	 */
        void TetExp::v_GetCoord(
                  const Array<OneD, const NekDouble> &Lcoords, 
                        Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] <= -1.0 && Lcoords[0] >= 1.0 &&
                     Lcoords[1] <= -1.0 && Lcoords[1] >= 1.0 &&
                     Lcoords[2] <= -1.0 && Lcoords[2] >= 1.0,
                     "Local coordinates are not in region [-1,1]");

            // m_geom->FillGeom(); // TODO: implement FillGeom()

            for(i = 0; i < m_geom->GetCoordim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }


        //-----------------------------
        // Helper functions
        //-----------------------------
        void TetExp::v_WriteToFile(
                  std::ofstream &outfile, 
                  OutputFormat format, 
                  const bool dumpVar, 
                  std::string var)
        {
            int i,j;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            Array<OneD,NekDouble> coords[3];

            ASSERTL0(m_geom,"m_geom not defined");

            int     coordim  = m_geom->GetCoordim();

            coords[0] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            coords[1] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            coords[2] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);

            GetCoords(coords[0],coords[1],coords[2]);

            if(format==eTecplot)
            {
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
                    outfile << std::endl << std::endl;
                }

                outfile << "Zone, I=" << nquad0 << ", J=" << nquad1 << ", K=" << nquad2 << ", F=Point" << std::endl;

                for(i = 0; i < nquad0*nquad1*nquad2; ++i)
                {
                    for(j = 0; j < coordim; ++j)
                    {
                        outfile << coords[j][i] << " ";
                    }
                    outfile << std::endl;
                }
            }
            else if (format==eGnuplot)
            {
                for(int k = 0; k < nquad2; ++k)
                {
                    for(int j = 0; j < nquad1; ++j)
                    {
                        for(int i = 0; i < nquad0; ++i)
                        {
                            int n = (k*nquad1 + j)*nquad0 + i;
                            outfile <<  coords[0][n] <<  " " << coords[1][n] << " "
                                    << coords[2][n] << " "
                                    << m_phys[i + nquad0*(j + nquad1*k)] << endl;
                        }
                        outfile << endl;
                    }
                    outfile << endl;
                }
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }


       /**
        * \brief Return Shape of region, using  ShapeType enum list.
	*/
        StdRegions::ExpansionType TetExp::v_DetExpansionType() const
        {
            return StdRegions::eTetrahedron;
        }


        const SpatialDomains::GeomFactorsSharedPtr& TetExp::v_GetMetricInfo() const
        {
            return m_metricinfo;
        }


        const SpatialDomains::GeometrySharedPtr TetExp::v_GetGeom() const
        {
            return m_geom;
        }


        const SpatialDomains::Geometry3DSharedPtr& TetExp::v_GetGeom3D() const
        {
            return m_geom;
        }


        int TetExp::v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }


        StdRegions::Orientation TetExp::v_GetFaceOrient(int face)
        {
            return m_geom->GetFaceOrient(face);
        }

      
        /**
         * \brief Returns the physical values at the quadrature points of a face
         */
        void TetExp::v_GetFacePhysVals(
            const int                                face,
            const StdRegions::StdExpansionSharedPtr &FaceExp,
            const Array<OneD, const NekDouble>      &inarray,
                  Array<OneD,       NekDouble>      &outarray,
            StdRegions::Orientation                  orient)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();

            Array<OneD,NekDouble> o_tmp (GetFaceNumPoints(face));
            Array<OneD,NekDouble> o_tmp2(FaceExp->GetTotPoints());
            Array<OneD,NekDouble> o_tmp3;
            
            if (orient == StdRegions::eNoOrientation)
            {
                orient = GetFaceOrient(face);
            }
            
            switch(face)
            {
                case 0:
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad0*nquad1,inarray.get(),1,o_tmp.get(),1);
                    //interpolate
                    LibUtilities::Interp2D(m_base[0]->GetPointsKey(), m_base[1]->GetPointsKey(), o_tmp.get(),
                                           FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),o_tmp2.get());
                    break;
                }
                case 1:
                {
                    //Direction A and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,inarray.get()+(nquad0*nquad1*k),1,o_tmp.get()+(k*nquad0),1);
                    }
                    //interpolate
                    LibUtilities::Interp2D(m_base[0]->GetPointsKey(), m_base[2]->GetPointsKey(), o_tmp.get(),
                                           FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),o_tmp2.get());
                    break;
                }
                case 2:
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad1*nquad2,inarray.get()+(nquad0-1),nquad0,o_tmp.get(),1);
                    //interpolate
                    LibUtilities::Interp2D(m_base[1]->GetPointsKey(), m_base[2]->GetPointsKey(), o_tmp.get(),
                                           FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),o_tmp2.get());
                    break;
                }
                case 3:
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad1*nquad2,inarray.get(),nquad0,o_tmp.get(),1);
                    //interpolate
                    LibUtilities::Interp2D(m_base[1]->GetPointsKey(), m_base[2]->GetPointsKey(), o_tmp.get(),
                                           FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),o_tmp2.get());
                }
                break;
            default:
                ASSERTL0(false,"face value (> 3) is out of range");
                break;
            }

            int nq1 = FaceExp->GetNumPoints(0);
            int nq2 = FaceExp->GetNumPoints(1);
            
            if ((int)orient == 7)
            {
                for (int j = 0; j < nq2; ++j)
                {
                    Vmath::Vcopy(nq1, o_tmp2.get()+((j+1)*nq1-1), -1, outarray.get()+j*nq1, 1);
                }
            }
            else
            {
                Vmath::Vcopy(nq1*nq2, o_tmp2.get(), 1, outarray.get(), 1);
            }
        }


        /**
	 * \brief Compute the normal of a triangular face
	 */
        void TetExp::v_ComputeFaceNormal(const int face)
        {
            int i;
            const SpatialDomains::GeomFactorsSharedPtr &geomFactors = 
                GetGeom()->GetMetricInfo();
            SpatialDomains::GeomType            type = geomFactors->GetGtype();
            const Array<TwoD, const NekDouble> &gmat = geomFactors->GetGmat();
            const Array<OneD, const NekDouble> &jac  = geomFactors->GetJac();
            int nq = m_base[0]->GetNumPoints()*m_base[0]->GetNumPoints();
            int vCoordDim = GetCoordim();
            
            m_faceNormals[face] = Array<OneD, Array<OneD, NekDouble> >(vCoordDim);
            Array<OneD, Array<OneD, NekDouble> > &normal = m_faceNormals[face];
            for (i = 0; i < vCoordDim; ++i)
            {
                normal[i] = Array<OneD, NekDouble>(nq);
            }

            // Regular geometry case
            if (type == SpatialDomains::eRegular ||
                type == SpatialDomains::eMovingRegular)
            {
                NekDouble fac;
                
                // Set up normals
                switch (face)
                {
                    case 0:
                    {
                        for (i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,-gmat[3*i+2][0],normal[i],1);
                        }
                        
                        break;
                    }
                    case 1:
                    {
                        for (i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,-gmat[3*i+1][0],normal[i],1);
                        }
                        
                        break;
                    }
                    case 2:
                    {
                        for (i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,gmat[3*i][0]+gmat[3*i+1][0]+
                                        gmat[3*i+2][0],normal[i],1);
                        }
                        
                        break;
                    }
                    case 3:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,-gmat[3*i][0],normal[i],1);
                        }
                        break;
                    }
                    default:
                        ASSERTL0(false,"face is out of range (edge < 3)");
                }
                
                // normalise
                fac = 0.0;
                for (i = 0; i < vCoordDim; ++i)
                {
                    fac += normal[i][0]*normal[i][0];
                }
                fac = 1.0/sqrt(fac);
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Smul(nq,fac,normal[i],1,normal[i],1);
                }
	    }
            else
            {
                // Set up deformed normals
                int j, k;

                int nq0 = geomFactors->GetPointsKey(0).GetNumPoints();
                int nq1 = geomFactors->GetPointsKey(1).GetNumPoints();
                int nq2 = geomFactors->GetPointsKey(2).GetNumPoints();
                int nqtot;
                int nq01 =nq0*nq1;

                if (face == 0)
                {
                    nqtot = nq01;
                }
                else if (face == 1)
                {
                    nqtot = nq0*nq2;
                }
                else
                {
                    nqtot = nq1*nq2;
                }



                LibUtilities::PointsKey points0;
                LibUtilities::PointsKey points1;

                Array<OneD,NekDouble> work   (nq,              0.0);
                Array<OneD,NekDouble> normals(vCoordDim*nqtot, 0.0);
                
                // Extract Jacobian along face and recover local derivates
                // (dx/dr) for polynomial interpolation by multiplying m_gmat by
                // jacobian
                switch (face)
	        {
                case 0:
                {
                    for(j = 0; j < nq01; ++j)
                    {
                        normals[j]         = -gmat[2][j]*jac[j];
                        normals[nqtot+j]   = -gmat[5][j]*jac[j];
                        normals[2*nqtot+j] = -gmat[8][j]*jac[j];
                    }
                    
                        points0 = geomFactors->GetPointsKey(0);
                        points1 = geomFactors->GetPointsKey(1);
                        break;
                }
                
                case 1:
                {
                    for (j = 0; j < nq0; ++j)
                    {
                        for(k = 0; k < nq2; ++k)
                        {
                            normals[j+k*nq0]  = 
                                -gmat[1][j+nq01*k]*
                                jac[j+nq01*k];
                            normals[nqtot+j+k*nq0]  = 
                                -gmat[4][j+nq01*k]*
                                jac[j+nq01*k];
                            normals[2*nqtot+j+k*nq0]  = 
                                -gmat[7][j+nq01*k]*
                                jac[j+nq01*k];
                        } 
                    }
                    
                    points0 = geomFactors->GetPointsKey(0);
                    points1 = geomFactors->GetPointsKey(2);
                    break;
                }
                
                case 2:
                {
                    for (j = 0; j < nq1; ++j)
                    {
                        for(k = 0; k < nq2; ++k)
                        {
                            normals[j+k*nq1]  = 
                                (gmat[0][nq0-1+nq0*j+nq01*k]+
                                 gmat[1][nq0-1+nq0*j+nq01*k]+
                                 gmat[2][nq0-1+nq0*j+nq01*k])*
                                jac[nq0-1+nq0*j+nq01*k];
                            normals[nqtot+j+k*nq1]  = 
                                (gmat[3][nq0-1+nq0*j+nq01*k]+
                                 gmat[4][nq0-1+nq0*j+nq01*k]+
                                 gmat[5][nq0-1+nq0*j+nq01*k])*
                                jac[nq0-1+nq0*j+nq01*k];
                            normals[2*nqtot+j+k*nq1]  = 
                                (gmat[6][nq0-1+nq0*j+nq01*k]+
                                 gmat[7][nq0-1+nq0*j+nq01*k]+
                                 gmat[8][nq0-1+nq0*j+nq01*k])*
                                jac[nq0-1+nq0*j+nq01*k];
                        } 
                    }
                    
                    points0 = geomFactors->GetPointsKey(1);
                    points1 = geomFactors->GetPointsKey(2);
                    break;
                }
                
                case 3:
                {
                    for (j = 0; j < nq1; ++j)
                    {
                        for(k = 0; k < nq2; ++k)
                        {
                            normals[j+k*nq1]  = 
                                -gmat[0][j*nq0+nq01*k]*
                                jac[j*nq0+nq01*k];
                            normals[nqtot+j+k*nq1]  = 
                                -gmat[3][j*nq0+nq01*k]*
                                jac[j*nq0+nq01*k];
                            normals[2*nqtot+j+k*nq1]  = 
                                -gmat[6][j*nq0+nq01*k]*
                                jac[j*nq0+nq01*k];
                        } 
                    }
                        
                    points0 = geomFactors->GetPointsKey(1);
                    points1 = geomFactors->GetPointsKey(2);
                    break;
                }
                
                    default:
                        ASSERTL0(false,"face is out of range (face < 3)");
                }

                // Interpolate Jacobian and invert
                LibUtilities::Interp2D(points0, points1, jac, 
                                       m_base[0]->GetPointsKey(),
                                       m_base[0]->GetPointsKey(),
                                       work);
                Vmath::Sdiv(nq, 1.0, &work[0], 1, &work[0], 1);
                    
                // Interpolate normal and multiply by inverse Jacobian.
                for(i = 0; i < vCoordDim; ++i)
                {
                    LibUtilities::Interp2D(points0, points1,
                                           &normals[i*nqtot],
                                           m_base[0]->GetPointsKey(),
                                           m_base[0]->GetPointsKey(),
                                           &normal[i][0]);
                    Vmath::Vmul(nq,work,1,normal[i],1,normal[i],1);
                }

                // Normalise to obtain unit normals.
                Vmath::Zero(nq,work,1);
                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vvtvp(nq,normal[i],1,normal[i],1,work,1,work,1);
                }
                
                Vmath::Vsqrt(nq,work,1,work,1);
                Vmath::Sdiv (nq,1.0,work,1,work,1);

                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vmul(nq,normal[i],1,work,1,normal[i],1);
		}
            }
        }


        NekDouble TetExp::v_Linf(const Array<OneD, const NekDouble> &sol)
        {
            return Linf(sol);
        }


        NekDouble TetExp::v_Linf()
        {
            return Linf();
        }


        NekDouble TetExp::v_L2(const Array<OneD, const NekDouble> &sol)
        {
            return StdExpansion::L2(sol);
        }


        NekDouble TetExp::v_L2()
        {
            return StdExpansion::L2();
        }


        //-----------------------------
        // Operator creation functions
        //-----------------------------
        void TetExp::v_HelmholtzMatrixOp(
                  const Array<OneD, const NekDouble> &inarray,
                        Array<OneD,NekDouble> &outarray,
                  const StdRegions::StdMatrixKey &mkey)
        {
            TetExp::v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
        }


        /**
	 * \brief Evaluates the Helmholtz operator using a matrix-free approach.
	 *
         * To construct the Helmholtz operator in a physical tetrahedron
         * requires coordinate transforms from both the collapsed coordinate
         * system to the standard region and from the standard region to the
         * local region. This double application of the chain rule requires the
         * calculation of two sets of geometric factors:
         * @f[ h_{ij} = \frac{\partial \eta_i}{\partial \xi_j} @f]
         * and
         * @f[ g_{ij} = \frac{\partial \xi_i}{\partial x_j} @f]
         *
         * From the definition of the collapsed coordinates, the @f$h_{ij}@f$
         * terms are (Sherwin & Karniadakis, p152)
         * @f[
         *      \mathbf{H} = \left[\begin{array}{ccc}
         *          \frac{4}{(1-\eta_2)(1-\eta_3)} &
         *          \frac{2(1+\eta_1)}{(1-\eta_2)(1-\eta_3)} &
         *          \frac{2(1+\eta_1)}{(1-\eta_2)(1-\eta_3)} \\
         *          0 &
         *          \frac{2}{1-eta_3} &
         *          \frac{1+\eta_2}{1-\eta_3} \\
         *          0 &
         *          0 &
         *          1
         *      \end{array}\right]
         * @f]
         * This maps from the collapsed coordinate system to the standard
         * tetrahedral region. The mapping to the local region is then given
         * by the @f$g_{ij}@f$ computed in the GeomFactors3D class. The
         * cumulative factors for mapping the collapsed coordinate system to
         * the physical region are therefore given by
         * @f$\mathbf{F} = \mathbf{GH^{\top}}@f$, i.e.
         * @f[
         *      f_{ij} = \frac{\partial \eta_i}{\partial x_j}
         *              = \sum_k g_{ik} h_{kj}
         * @f]
         *
         * Finally, the evaluation of the Helmholtz matrix operator requires
         * the summation of these factors as follows. For the case of deformed
         * elements, these coefficients are vectors, whereas for regular
         * elements they are just scalars.
         * @f[
         *      \begin{array}{l}
         *      p_0 = \sum_k f_{1k}^2 \\
         *      p_1 = \sum_k f_{2k}^2 \\
         *      p_2 = \sum_k f_{3k}^2 \\
         *      p_3 = \sum_k f_{1k}f_{2k} \\
         *      p_4 = \sum_k f_{1k}f_{3k} \\
         *      p_5 = \sum_k f_{2k}f_{3k}
         *      \end{array}
         * @f]
         * to give the Helmholtz operator:
         * @f{align}
         *      \mathbf{L^e\hat{u}}
         *          = \mathbf{B^{\top}D_{\eta_1}^{\top}Wp_0D_{\eta_1}B\hat{u}}
         *          + \mathbf{B^{\top}D_{\eta_2}^{\top}Wp_1D_{\eta_2}B\hat{u}}
         *          + \mathbf{B^{\top}D_{\eta_3}^{\top}Wp_2D_{\eta_3}B\hat{u}}\\
         *          + \mathbf{B^{\top}D_{\eta_1}^{\top}Wp_3D_{\eta_2}B\hat{u}}
         *          + \mathbf{B^{\top}D_{\eta_1}^{\top}Wp_4D_{\eta_3}B\hat{u}}
         *          + \mathbf{B^{\top}D_{\eta_2}^{\top}Wp_5D_{\eta_3}B\hat{u}}
         * @f}
         * Therefore, we construct the operator as follows:
         * -# Apply the mass matrix for the @f$\lambda@f$ term
         *    @f$ \mathbf{B^{\top}WB\hat{u}} @f$.
         *    and compute the derivatives @f$ \mathbf{D_{\xi_i}B} @f$.
         * -# Compute the non-trivial @f$ \mathbf{H} @f$ matrix terms.
         * -# Compute the intermediate factors @f$ \mathbf{G} @f$ and
         *    @f$ f_{ij} @f$ and then compute the combined terms @f$ p_i @f$.
         * -# Apply quadrature weights and inner product with respect to the
         *    derivative bases.
         * -# Combine to produce the complete operator.
         */
        void TetExp::v_HelmholtzMatrixOp_MatFree(
                  const Array<OneD, const NekDouble> &inarray,
                        Array<OneD,NekDouble> &outarray,
                  const StdRegions::StdMatrixKey &mkey)
        {
            if(m_metricinfo->IsUsingLaplMetrics())
            {
                ASSERTL0(false,"Finish implementing TetExp Helmholtz for Lapl Metrics");
/*
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
*/            }
            else
            {
                int nquad0  = m_base[0]->GetNumPoints();
                int nquad1  = m_base[1]->GetNumPoints();
                int nquad2  = m_base[2]->GetNumPoints();
                int nqtot   = nquad0*nquad1*nquad2;
                int nmodes0 = m_base[0]->GetNumModes();
                int nmodes1 = m_base[1]->GetNumModes();
                int nmodes2 = m_base[2]->GetNumModes();
                int wspsize = max(nquad0*nmodes2*(nmodes1+nquad1),
                                  nquad2*nmodes0*nmodes1*(nmodes1+1)/2+
                                  nquad2*nquad1*nmodes0);
                
                NekDouble lambda  = mkey.GetConstFactor(StdRegions::eFactorLambda);

                const Array<OneD, const NekDouble>& base0 = m_base[0]->GetBdata ();
                const Array<OneD, const NekDouble>& base1 = m_base[1]->GetBdata ();
                const Array<OneD, const NekDouble>& base2 = m_base[2]->GetBdata ();
                Array<OneD,NekDouble> wsp (wspsize);
                Array<OneD,NekDouble> wsp0(nqtot);
                Array<OneD,NekDouble> wsp1(nqtot);
                
                if(!(m_base[0]->Collocation() && m_base[1]->Collocation() &&
                     m_base[2]->Collocation()))
                {
                    // MASS MATRIX OPERATION
                    // The following is being calculated:
                    // wsp0     = B   * u_hat = u
                    // wsp1     = W   * wsp0
                    // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat
                    BwdTrans_SumFacKernel           (base0,base1,base2,inarray,
                                                     wsp0,wsp,true,true,true);
                    MultiplyByQuadratureMetric      (wsp0,wsp1);
                    IProductWRTBase_SumFacKernel    (base0,base1,base2,wsp1,
                                                     outarray,wsp,true,true,true);
                    LaplacianMatrixOp_MatFree_Kernel(wsp0,wsp1,wsp);
                }
                else
                {
                    // specialised implementation for the classical spectral
                    // element method
                    MultiplyByQuadratureMetric      (inarray,outarray);
                    LaplacianMatrixOp_MatFree_Kernel(inarray,wsp1,wsp);
                }

                // outarray = lambda * outarray + wsp1
                //          = (lambda * M + L ) * u_hat
                Vmath::Svtvp(m_ncoeffs,lambda,&outarray[0],1,&wsp1[0],1,
                             &outarray[0],1); 
            }
        }


        void TetExp::v_LaplacianMatrixOp(
                  const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
                  const StdRegions::StdMatrixKey &mkey)
        {
            TetExp::v_LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void TetExp::v_LaplacianMatrixOp(
                  const int k1, 
                  const int k2,
                  const Array<OneD, const NekDouble> &inarray,
                        Array<OneD,NekDouble> &outarray,
                  const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,
                                                        mkey);
        }

        void TetExp::v_LaplacianMatrixOp_MatFree(
                  const Array<OneD, const NekDouble> &inarray,
                        Array<OneD,NekDouble> &outarray,
                  const StdRegions::StdMatrixKey &mkey)
        {
            if(mkey.GetNVarCoeff() == 0)
            {
                // This implementation is only valid when there are no coefficients
                // associated to the Laplacian operator
                if(m_metricinfo->IsUsingLaplMetrics())
                {
                    ASSERTL0(false,"Finish implementing HexExp for Lap metrics");
                    // Get this from HexExp
                }
                else
                {
                    int nquad0  = m_base[0]->GetNumPoints();
                    int nquad1  = m_base[1]->GetNumPoints();
                    int nquad2  = m_base[2]->GetNumPoints();
                    int nqtot   = nquad0*nquad1*nquad2;
                    int nmodes0 = m_base[0]->GetNumModes();
                    int nmodes1 = m_base[1]->GetNumModes();
                    int nmodes2 = m_base[2]->GetNumModes();
                    int wspsize = max(nquad0*nmodes2*(nmodes1+nquad1),
                                      nquad2*nmodes0*nmodes1*(nmodes1+1)/2+
                                      nquad2*nquad1*nmodes0);
                    
                    const Array<OneD, const NekDouble>& base0 = m_base[0]->GetBdata ();
                    const Array<OneD, const NekDouble>& base1 = m_base[1]->GetBdata ();
                    const Array<OneD, const NekDouble>& base2 = m_base[2]->GetBdata ();

                    Array<OneD,NekDouble> wsp (wspsize);
                    Array<OneD,NekDouble> wsp1(nqtot);
                    
                    // Backwards transform to obtain u = B * u_hat.
                    if(!(m_base[0]->Collocation() && m_base[1]->Collocation() &&
                         m_base[2]->Collocation()))
                    {
                        BwdTrans_SumFacKernel           (base0,base1,base2,inarray,
                                                         wsp1,wsp,true,true,true);
                        LaplacianMatrixOp_MatFree_Kernel(wsp1, outarray, wsp);
                    }
                    else
                    {
                        LaplacianMatrixOp_MatFree_Kernel(wsp1, outarray, wsp);
                    }
                }
            }
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }
        }


        //-----------------------------
        // Matrix creation functions
        //-----------------------------
        DNekMatSharedPtr TetExp::v_GenMatrix(
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
                returnval = Expansion3D::v_GenMatrix(mkey);
                break;
            default:
                returnval = StdTetExp::v_GenMatrix(mkey);
            }

            return returnval;
        }


        DNekScalMatSharedPtr TetExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed ||
                            mkey.GetNVarCoeff())
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
            case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
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
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat
                                                = m_metricinfo->GetGmat();
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
                            default:
                                break;
                        }

                        MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                            mkey.GetExpansionType(), *this);
                        MatrixKey deriv1key(StdRegions::eWeakDeriv1,
                                            mkey.GetExpansionType(), *this);
                        MatrixKey deriv2key(StdRegions::eWeakDeriv2,
                                            mkey.GetExpansionType(), *this);

                        DNekMat &deriv0 = *GetStdMatrix(deriv0key);
                        DNekMat &deriv1 = *GetStdMatrix(deriv1key);
                        DNekMat &deriv2 = *GetStdMatrix(deriv2key);

                        int rows = deriv0.GetRows();
                        int cols = deriv1.GetColumns();

                        DNekMatSharedPtr WeakDeriv = MemoryManager<DNekMat>
                                                ::AllocateSharedPtr(rows,cols);
                        (*WeakDeriv) = gmat[3*dir][0]*deriv0
                                                + gmat[3*dir+1][0]*deriv1
												+ gmat[3*dir+2][0]*deriv2;

                        returnval = MemoryManager<DNekScalMat>
                                            ::AllocateSharedPtr(jac,WeakDeriv);
                    }
                }
                break;
            case StdRegions::eLaplacian:
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
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap02key(StdRegions::eLaplacian02,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap12key(StdRegions::eLaplacian12,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap22key(StdRegions::eLaplacian22,
                                           mkey.GetExpansionType(), *this);

                        DNekMat &lap00 = *GetStdMatrix(lap00key);
                        DNekMat &lap01 = *GetStdMatrix(lap01key);
                        DNekMat &lap02 = *GetStdMatrix(lap02key);
                        DNekMat &lap11 = *GetStdMatrix(lap11key);
                        DNekMat &lap12 = *GetStdMatrix(lap12key);
                        DNekMat &lap22 = *GetStdMatrix(lap22key);

                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat
                                                    = m_metricinfo->GetGmat();

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>
                                                ::AllocateSharedPtr(rows,cols);

                        (*lap) = (gmat[0][0]*gmat[0][0] + gmat[3][0]*gmat[3][0]
                                        + gmat[6][0]*gmat[6][0])*lap00
                               + (gmat[1][0]*gmat[1][0] + gmat[4][0]*gmat[4][0]
                                        + gmat[7][0]*gmat[7][0])*lap11
                               + (gmat[2][0]*gmat[2][0] + gmat[5][0]*gmat[5][0]
                                        + gmat[8][0]*gmat[8][0])*lap22
                               + (gmat[0][0]*gmat[1][0] + gmat[3][0]*gmat[4][0]
                                        + gmat[6][0]*gmat[7][0])
                                 *(lap01 + Transpose(lap01))
                               + (gmat[0][0]*gmat[2][0] + gmat[3][0]*gmat[5][0]
                                        + gmat[6][0]*gmat[8][0])
                                 *(lap02 + Transpose(lap02))
                               + (gmat[1][0]*gmat[2][0] + gmat[4][0]*gmat[5][0]
                                        + gmat[7][0]*gmat[8][0])
                                 *(lap12 + Transpose(lap12));

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(jac,lap);
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    MatrixKey masskey(StdRegions::eMass, mkey.GetExpansionType(), *this);
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian, mkey.GetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, helm);
                }
                break;
             case StdRegions::ePreconditioner:
                {
                    LibUtilities::BasisKey TetBa = m_base[0]->GetBasisKey();
		    LibUtilities::BasisKey TetBb = m_base[1]->GetBasisKey();
		    LibUtilities::BasisKey TetBc = m_base[2]->GetBasisKey();

		    SpatialDomains::TetGeomSharedPtr EquilateralTetGeom=CreateEquilateralTetGeom();

		    //create TetExp with equilateral Tet geometry object
                    TetExp eqtet(TetBa,TetBb,TetBc,EquilateralTetGeom);
		
		    int nquad0 = m_base[0]->GetNumPoints();
		    int nquad1 = m_base[1]->GetNumPoints();
		    int nquad2 = m_base[2]->GetNumPoints();

		    int nq=nquad0*nquad1*nquad2;
		    Array<OneD,NekDouble> coords[3];

		    coords[0] = Array<OneD,NekDouble>(nq);
		    coords[1] = Array<OneD,NekDouble>(nq);
		    coords[2] = Array<OneD,NekDouble>(nq);
		    eqtet.GetCoords(coords[0],coords[1],coords[2]);

                    NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    MatrixKey masskey(StdRegions::eMass, mkey.GetExpansionType(), eqtet);
                    DNekScalMat &MassMat = *(eqtet.m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian, mkey.GetExpansionType(), eqtet, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LapMat = *(eqtet.m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, helm);
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
			case StdRegions::eHybridDGHelmholtz:
			case StdRegions::eHybridDGLamToU:
			case StdRegions::eHybridDGLamToQ0:
			case StdRegions::eHybridDGLamToQ1:
			case StdRegions::eHybridDGHelmBndLam:
				{
					NekDouble one    = 1.0;

					DNekMatSharedPtr mat = GenMatrix(mkey);
					returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
				}
				break;
			case StdRegions::eInvHybridDGHelmholtz:
				{
					NekDouble one = 1.0;

					MatrixKey hkey(StdRegions::eHybridDGHelmholtz, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
					DNekMatSharedPtr mat = GenMatrix(hkey);

					mat->Invert();
					returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
				}
				break;
			default:
				{
					//ASSERTL0(false, "Missing definition for " + (*StdRegions::MatrixTypeMap[mkey.GetMatrixType()]));
					NekDouble        one = 1.0;
					DNekMatSharedPtr mat = GenMatrix(mkey);

					returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
				}
				break;
			}

            return returnval;
        }


        DNekScalBlkMatSharedPtr TetExp::CreateStaticCondMatrix(
                 const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;

            unsigned int exp_size[] = {nbdry, nint};
            int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks, nblks, exp_size, exp_size); //Really need a constructor which takes Arrays
            NekDouble factor = 1.0;
            MatrixStorage AMatStorage = eFULL;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eLaplacian:
            case StdRegions::ePreconditioner:
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

                    //TODO: check below
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
                    DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry,AMatStorage);
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


        DNekMatSharedPtr TetExp::v_CreateStdMatrix(
                 const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
            StdRegions::StdTetExpSharedPtr tmp = MemoryManager<StdTetExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);

            return tmp->GetStdMatrix(mkey);
        }

        DNekScalMatSharedPtr TetExp::v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }

        DNekScalBlkMatSharedPtr TetExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
        }

        void TetExp::GeneralMatrixOp_MatOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey)
        {
            DNekScalMatSharedPtr   mat = GetLocMatrix(mkey);

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

        /// @todo add functionality for IsUsingQuadMetrics
        void TetExp::MultiplyByQuadratureMetric(
                  const Array<OneD, const NekDouble>& inarray,
                        Array<OneD, NekDouble> &outarray)
        {
            if(m_metricinfo->IsUsingQuadMetrics())
            {
                const Array<OneD, const NekDouble> &metric 
                    = m_metricinfo->GetQuadratureMetrics();

                Vmath::Vmul(metric.num_elements(), metric, 1, inarray, 1,
                            outarray, 1);
            }
            else
            {
                const int nqtot = m_base[0]->GetNumPoints() *
                                  m_base[1]->GetNumPoints() *
                                  m_base[2]->GetNumPoints();
                const Array<OneD, const NekDouble> &jac 
                    = m_metricinfo->GetJac();
                
                if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                {
                    Vmath::Vmul(nqtot, jac, 1, inarray, 1, outarray, 1);
                }
                else
                {
                    Vmath::Smul(nqtot, jac[0], inarray, 1, outarray, 1);
                }

                StdTetExp::MultiplyByQuadratureMetric(outarray, outarray);
            }
        }
        
        void TetExp::LaplacianMatrixOp_MatFree_Kernel(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
                  Array<OneD,       NekDouble> &wsp)
        {
            int nquad0  = m_base[0]->GetNumPoints();
            int nquad1  = m_base[1]->GetNumPoints();
            int nquad2  = m_base[2]->GetNumPoints();
            int nqtot   = nquad0*nquad1*nquad2;
            int i, j;

            // Allocate temporary storage
            Array<OneD,NekDouble> alloc(13*nqtot,0.0);
            Array<OneD,NekDouble> wsp1 (alloc         );// TensorDeriv 1
            Array<OneD,NekDouble> wsp2 (alloc+ 1*nqtot);// TensorDeriv 2
            Array<OneD,NekDouble> wsp3 (alloc+ 2*nqtot);// TensorDeriv 3
            Array<OneD,NekDouble> g0   (alloc+ 3*nqtot);// g0
            Array<OneD,NekDouble> g1   (alloc+ 4*nqtot);// g1
            Array<OneD,NekDouble> g2   (alloc+ 5*nqtot);// g2
            Array<OneD,NekDouble> g3   (alloc+ 6*nqtot);// g3
            Array<OneD,NekDouble> g4   (alloc+ 7*nqtot);// g4
            Array<OneD,NekDouble> g5   (alloc+ 8*nqtot);// g5
            Array<OneD,NekDouble> h0   (alloc+ 9*nqtot);// h0
            Array<OneD,NekDouble> h1   (alloc+10*nqtot);// h1
            Array<OneD,NekDouble> h2   (alloc+11*nqtot);// h2
            Array<OneD,NekDouble> h3   (alloc+12*nqtot);// h3
            // Reuse some of the storage as workspace
            Array<OneD,NekDouble> wsp4 (alloc+ 4*nqtot);// wsp4 == g1
            Array<OneD,NekDouble> wsp5 (alloc+ 5*nqtot);// wsp5 == g2
            Array<OneD,NekDouble> wsp6 (alloc+ 8*nqtot);// wsp6 == g5
            Array<OneD,NekDouble> wsp7 (alloc+ 9*nqtot);// wsp7 == h0
            Array<OneD,NekDouble> wsp8 (alloc+10*nqtot);// wsp8 == h1
            Array<OneD,NekDouble> wsp9 (alloc+11*nqtot);// wsp9 == h2
            
            const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble>& base2  = m_base[2]->GetBdata();
            const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase2 = m_base[2]->GetDbdata();
            
            // LAPLACIAN MATRIX OPERATION
            // wsp1 = du_dxi1 = D_xi1 * inarray = D_xi1 * u
            // wsp2 = du_dxi2 = D_xi2 * inarray = D_xi2 * u
            // wsp2 = du_dxi3 = D_xi3 * inarray = D_xi3 * u
            StdExpansion3D::PhysTensorDeriv(inarray,wsp1,wsp2,wsp3);
            
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();
            const Array<OneD, const NekDouble>& z0   = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z1   = m_base[1]->GetZ();
            const Array<OneD, const NekDouble>& z2   = m_base[2]->GetZ();
            
            // Step 2. Calculate the metric terms of the collapsed
            // coordinate transformation (Spencer's book P152)
            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Vmath::Fill(nquad0, 4.0/(1.0-z1[i])/(1.0-z2[j]), &h0[0]+i*nquad0 + j*nquad0*nquad1,1);
                    Vmath::Fill(nquad0, 2.0/(1.0-z1[i])/(1.0-z2[j]), &h1[0]+i*nquad0 + j*nquad0*nquad1,1);
                    Vmath::Fill(nquad0, 2.0/(1.0-z2[j]),             &h2[0]+i*nquad0 + j*nquad0*nquad1,1);
                    Vmath::Fill(nquad0, (1.0+z1[i])/(1.0-z2[j]),     &h3[0]+i*nquad0 + j*nquad0*nquad1,1);
                }
            }
            for(i = 0; i < nquad0; i++)
            {
                Blas::Dscal(nquad1*nquad2, 1+z0[i], &h1[0]+i, nquad0);
            }
            
            // Step 3. Construct combined metric terms for physical space to
            // collapsed coordinate system.
            // Order of construction optimised to minimise temporary storage
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // wsp4
                Vmath::Vadd(nqtot, &gmat[1][0], 1, &gmat[2][0], 1, &wsp4[0], 1);
                Vmath::Vvtvvtp(nqtot, &gmat[0][0], 1, &h0[0], 1, &wsp4[0], 1, &h1[0], 1, &wsp4[0], 1);
                // wsp5
                Vmath::Vadd(nqtot, &gmat[4][0], 1, &gmat[5][0], 1, &wsp5[0], 1);
                Vmath::Vvtvvtp(nqtot, &gmat[3][0], 1, &h0[0], 1, &wsp5[0], 1, &h1[0], 1, &wsp5[0], 1);
                // wsp6
                Vmath::Vadd(nqtot, &gmat[7][0], 1, &gmat[8][0], 1, &wsp6[0], 1);
                Vmath::Vvtvvtp(nqtot, &gmat[6][0], 1, &h0[0], 1, &wsp6[0], 1, &h1[0], 1, &wsp6[0], 1);
                
                // g0
                Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp4[0], 1, &wsp5[0], 1, &wsp5[0], 1, &g0[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp6[0], 1, &wsp6[0], 1, &g0[0],   1, &g0[0],   1);
                
                // g4
                Vmath::Vvtvvtp(nqtot, &gmat[2][0], 1, &wsp4[0], 1, &gmat[5][0], 1, &wsp5[0], 1, &g4[0], 1);
                Vmath::Vvtvp  (nqtot, &gmat[8][0], 1, &wsp6[0], 1, &g4[0], 1, &g4[0], 1);
                
                // overwrite h0, h1, h2
                // wsp7 (h2f1 + h3f2)
                Vmath::Vvtvvtp(nqtot, &gmat[1][0], 1, &h2[0], 1, &gmat[2][0], 1, &h3[0], 1, &wsp7[0], 1);
                // wsp8 (h2f4 + h3f5)
                Vmath::Vvtvvtp(nqtot, &gmat[4][0], 1, &h2[0], 1, &gmat[5][0], 1, &h3[0], 1, &wsp8[0], 1);
                // wsp9 (h2f7 + h3f8)
                Vmath::Vvtvvtp(nqtot, &gmat[7][0], 1, &h2[0], 1, &gmat[8][0], 1, &h3[0], 1, &wsp9[0], 1);
                
                // g3
                Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp7[0], 1, &wsp5[0], 1, &wsp8[0], 1, &g3[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp6[0], 1, &wsp9[0], 1, &g3[0],   1, &g3[0],   1);
                
                // overwrite wsp4, wsp5, wsp6
                // g1
                Vmath::Vvtvvtp(nqtot, &wsp7[0], 1, &wsp7[0], 1, &wsp8[0], 1, &wsp8[0], 1, &g1[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp9[0], 1, &wsp9[0], 1, &g1[0],   1, &g1[0],   1);
                
                // g5
                Vmath::Vvtvvtp(nqtot, &gmat[2][0], 1, &wsp7[0], 1, &gmat[5][0], 1, &wsp8[0], 1, &g5[0], 1);
                Vmath::Vvtvp  (nqtot, &gmat[8][0], 1, &wsp9[0], 1, &g5[0], 1, &g5[0], 1);
                
                // g2
                Vmath::Vvtvvtp(nqtot, &gmat[2][0], 1, &gmat[2][0], 1, &gmat[5][0], 1, &gmat[5][0], 1, &g2[0], 1);
                Vmath::Vvtvp  (nqtot, &gmat[8][0], 1, &gmat[8][0], 1, &g2[0],      1, &g2[0],      1);
            }
            else
            {
                // wsp4
                Vmath::Svtsvtp(nqtot, gmat[0][0], &h0[0], 1, gmat[1][0] + gmat[2][0], &h1[0], 1, &wsp4[0], 1);
                // wsp5
                Vmath::Svtsvtp(nqtot, gmat[3][0], &h0[0], 1, gmat[4][0] + gmat[5][0], &h1[0], 1, &wsp5[0], 1);
                // wsp6
                Vmath::Svtsvtp(nqtot, gmat[6][0], &h0[0], 1, gmat[7][0] + gmat[8][0], &h1[0], 1, &wsp6[0], 1);
                
                // g0
                Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp4[0], 1, &wsp5[0], 1, &wsp5[0], 1, &g0[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp6[0], 1, &wsp6[0], 1, &g0[0],   1, &g0[0],   1);
                
                // g4
                Vmath::Svtsvtp(nqtot, gmat[2][0], &wsp4[0], 1, gmat[5][0], &wsp5[0], 1, &g4[0], 1);
                Vmath::Svtvp  (nqtot, gmat[8][0], &wsp6[0], 1, &g4[0], 1, &g4[0], 1);
                
                // overwrite h0, h1, h2
                // wsp7 (h2f1 + h3f2)
                Vmath::Svtsvtp(nqtot, gmat[1][0], &h2[0], 1, gmat[2][0], &h3[0], 1, &wsp7[0], 1);
                // wsp8 (h2f4 + h3f5)
                Vmath::Svtsvtp(nqtot, gmat[4][0], &h2[0], 1, gmat[5][0], &h3[0], 1, &wsp8[0], 1);
                // wsp9 (h2f7 + h3f8)
                Vmath::Svtsvtp(nqtot, gmat[7][0], &h2[0], 1, gmat[8][0], &h3[0], 1, &wsp9[0], 1);
                
                // g3
                Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp7[0], 1, &wsp5[0], 1, &wsp8[0], 1, &g3[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp6[0], 1, &wsp9[0], 1, &g3[0],   1, &g3[0],   1);
                
                // overwrite wsp4, wsp5, wsp6
                // g1
                Vmath::Vvtvvtp(nqtot, &wsp7[0], 1, &wsp7[0], 1, &wsp8[0], 1, &wsp8[0], 1, &g1[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp9[0], 1, &wsp9[0], 1, &g1[0],   1, &g1[0],   1);
                
                // g5
                Vmath::Svtsvtp(nqtot, gmat[2][0], &wsp7[0], 1, gmat[5][0], &wsp8[0], 1, &g5[0], 1);
                Vmath::Svtvp  (nqtot, gmat[8][0], &wsp9[0], 1, &g5[0], 1, &g5[0], 1);
                
                // g2
                Vmath::Fill(nqtot, gmat[2][0]*gmat[2][0] + gmat[5][0]*gmat[5][0] + gmat[8][0]*gmat[8][0], &g2[0], 1);
            }
            
            // Compute component derivatives into wsp7, 8, 9
            Vmath::Vvtvvtp(nqtot,&g0[0],1,&wsp1[0],1,&g3[0],1,&wsp2[0],1,&wsp7[0],1);
            Vmath::Vvtvp  (nqtot,&g4[0],1,&wsp3[0],1,&wsp7[0],1,&wsp7[0],1);
            Vmath::Vvtvvtp(nqtot,&g1[0],1,&wsp2[0],1,&g3[0],1,&wsp1[0],1,&wsp8[0],1);
            Vmath::Vvtvp  (nqtot,&g5[0],1,&wsp3[0],1,&wsp8[0],1,&wsp8[0],1);
            Vmath::Vvtvvtp(nqtot,&g2[0],1,&wsp3[0],1,&g4[0],1,&wsp1[0],1,&wsp9[0],1);
            Vmath::Vvtvp  (nqtot,&g5[0],1,&wsp2[0],1,&wsp9[0],1,&wsp9[0],1);
            
            // Step 4.
            // Multiply by quadrature metric
            MultiplyByQuadratureMetric(wsp7,wsp7);
            MultiplyByQuadratureMetric(wsp8,wsp8);
            MultiplyByQuadratureMetric(wsp9,wsp9);
            
            IProductWRTBase_SumFacKernel(dbase0,base1,base2,wsp7,wsp1,    wsp,false,true,true);
            IProductWRTBase_SumFacKernel(base0,dbase1,base2,wsp8,wsp2,    wsp,true,false,true);
            IProductWRTBase_SumFacKernel(base0,base1,dbase2,wsp9,outarray,wsp,true,true,false);
            
            // Step 5.
            // Sum contributions from wsp1, wsp2 and outarray.
            Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);
            Vmath::Vadd(m_ncoeffs,wsp2.get(),1,outarray.get(),1,outarray.get(),1);
        }

        SpatialDomains::TetGeomSharedPtr TetExp::CreateEquilateralTetGeom()
        {
	    int i,j;
	    const int three=3;
            const int nVerts = 4;
            const NekDouble point[][3] = {
	      {-1,-1/sqrt(NekDouble(3)),-1/sqrt(NekDouble(6))},
	      {1,-1/sqrt(NekDouble(3)),-1/sqrt(NekDouble(6))},
	      {0,2/sqrt(NekDouble(3)),-1/sqrt(NekDouble(6))},
	      {0,0,3/sqrt(NekDouble(6))}};
        
            boost::shared_ptr<SpatialDomains::VertexComponent> verts[4];
	    for(i=0; i < nVerts; ++i)
	    {
	        verts[i] =  MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr( three, i, point[i][0], point[i][1], point[i][2] );
	    }

            //////////////////////////////
            // Set up Tetrahedron Edges //
            //////////////////////////////

           // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
           const int nEdges = 6;
           const int vertexConnectivity[][2] = {
             {0,1},{1,2},{0,2},{0,3},{1,3},{2,3}
           };

           // Populate the list of edges
	   SpatialDomains::SegGeomSharedPtr edges[nEdges];
           for(i=0; i < nEdges; ++i)
	   {
               boost::shared_ptr<SpatialDomains::VertexComponent> vertsArray[2];
	       for(j=0; j<2; ++j)
	       {
                   vertsArray[j] = verts[vertexConnectivity[i][j]];
	       }

               edges[i] = MemoryManager<SpatialDomains::SegGeom>
                                ::AllocateSharedPtr(i, three, vertsArray);
           }

           //////////////////////////////
           // Set up Tetrahedron faces //
           //////////////////////////////
 
	   const int nFaces = 4;
           const int edgeConnectivity[][3] = {
                 {0,1,2}, {0,4,3}, {1,5,4}, {2,5,3}
                 };
           const bool isEdgeFlipped[][3] = {
                 {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}
                 };

           // Populate the list of faces
	   SpatialDomains::TriGeomSharedPtr faces[nFaces];
           for(i=0; i < nFaces; ++i)
	   {
	       SpatialDomains::SegGeomSharedPtr edgeArray[3];
	       StdRegions::Orientation eorientArray[3];
               for(j=0; j < 3; ++j)
	       {
                   edgeArray[j] = edges[edgeConnectivity[i][j]];
                   eorientArray[j] = isEdgeFlipped[i][j] ? StdRegions::eBackwards : StdRegions::eForwards;
	       }
           

	       faces[i] = MemoryManager<SpatialDomains::TriGeom>
                                ::AllocateSharedPtr(i, edgeArray, eorientArray);
	   }

           SpatialDomains::TetGeomSharedPtr geom =
                         MemoryManager<SpatialDomains::TetGeom>::AllocateSharedPtr(faces);

	   geom->SetOwnData();

           return geom;
       }

    }//end of namespace
}//end of namespace
