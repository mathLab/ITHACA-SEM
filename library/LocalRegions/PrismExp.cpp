///////////////////////////////////////////////////////////////////////////////
//
// File PrismExp.cpp
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
// Description:  PrismExp routines
//
///////////////////////////////////////////////////////////////////////////////


#include <LocalRegions/PrismExp.h>
#include <LibUtilities/Foundations/Interp.h>

namespace Nektar
{
    namespace LocalRegions 
    {

        PrismExp::PrismExp(const LibUtilities::BasisKey &Ba,
                           const LibUtilities::BasisKey &Bb,
                           const LibUtilities::BasisKey &Bc,
                           const SpatialDomains::PrismGeomSharedPtr &geom):
            StdExpansion  (StdRegions::StdPrismData::getNumberOfCoefficients(
                               Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                           3, Ba, Bb, Bc),
            StdExpansion3D(StdRegions::StdPrismData::getNumberOfCoefficients(
                               Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()), 
                           Ba, Bb, Bc),
            StdPrismExp   (Ba, Bb, Bc),
            Expansion     (),
            Expansion3D   (),
            m_geom        (geom),
            m_metricinfo  (m_geom->GetGeomFactors(m_base)),
            m_matrixManager(
                    boost::bind(&PrismExp::CreateMatrix, this, _1),
                    std::string("PrismExpMatrix")),
            m_staticCondMatrixManager(
                    boost::bind(&PrismExp::CreateStaticCondMatrix, this, _1),
                    std::string("PrismExpStaticCondMatrix"))
        {
        }

        PrismExp::PrismExp(const PrismExp &T):
            StdExpansion(T),
            StdExpansion3D(T),
            StdRegions::StdPrismExp(T),
            Expansion(T),
            Expansion3D(T),
            m_geom(T.m_geom),
            m_metricinfo(T.m_metricinfo),
            m_matrixManager(T.m_matrixManager),
            m_staticCondMatrixManager(T.m_staticCondMatrixManager)
        {
        } 

        PrismExp::~PrismExp()
        {
        }


        //-------------------------------
        // Integration Methods
        //-------------------------------
        
        /** 
         * \brief Integrate the physical point list \a inarray over prismatic
         * region and return the value.
         *
         * Inputs:\n
         *
         * - \a inarray: definition of function to be returned at quadrature
         * point of expansion.
         *
         * Outputs:\n
         * 
         * - returns \f$\int^1_{-1}\int^1_{-1}\int^1_{-1} u(\bar \eta_1,
         *  \xi_2, \xi_3) J[i,j,k] d \bar \eta_1 d \xi_2 d \xi_3 \f$ \n \f$ =
         *  \sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1} \sum_{k=0}^{Q_3 - 1}
         *  u(\bar \eta_{1i}^{0,0}, \xi_{2j}^{0,0},\xi_{3k}^{1,0})w_{i}^{0,0}
         *  w_{j}^{0,0} \hat w_{k}^{1,0} \f$ \n where \f$ inarray[i,j, k] =
         *  u(\bar \eta_{1i}^{0,0}, \xi_{2j}^{0,0},\xi_{3k}^{1,0}) \f$, \n
         *  \f$\hat w_{i}^{1,0} = \frac {w_{j}^{1,0}} {2} \f$ \n and \f$
         *  J[i,j,k] \f$ is the Jacobian evaluated at the quadrature point.
        */
        NekDouble PrismExp::v_Integral(const Array<OneD, const NekDouble> &inarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,       NekDouble> tmp(nquad0*nquad1*nquad2);

            // Multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,(NekDouble)jac[0],(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            
            // Call StdPrismExp version.
            return StdPrismExp::v_Integral(tmp);
        }


        //----------------------------
        // Differentiation Methods
        //----------------------------
        void PrismExp::v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD,       NekDouble>& out_d0,
                                         Array<OneD,       NekDouble>& out_d1,
                                         Array<OneD,       NekDouble>& out_d2)
        {
            int nqtot = GetTotPoints();
            
            Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
            Array<OneD,       NekDouble> diff0(nqtot);
            Array<OneD,       NekDouble> diff1(nqtot);
            Array<OneD,       NekDouble> diff2(nqtot);

            StdPrismExp::v_PhysDeriv(inarray, diff0, diff1, diff2);
            
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul  (nqtot,&gmat[0][0],1,&diff0[0],1,&out_d0[0],1);
                    Vmath::Vvtvp (nqtot,&gmat[1][0],1,&diff1[0],1,&out_d0[0],1,&out_d0[0],1);
                    Vmath::Vvtvp (nqtot,&gmat[2][0],1,&diff2[0],1,&out_d0[0],1,&out_d0[0],1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Vmul  (nqtot,&gmat[3][0],1,&diff0[0],1,&out_d1[0],1);
                    Vmath::Vvtvp (nqtot,&gmat[4][0],1,&diff1[0],1,&out_d1[0],1,&out_d1[0],1);
                    Vmath::Vvtvp (nqtot,&gmat[5][0],1,&diff2[0],1,&out_d1[0],1,&out_d1[0],1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Vmul  (nqtot,&gmat[6][0],1,&diff0[0],1,&out_d2[0],1);
                    Vmath::Vvtvp (nqtot,&gmat[7][0],1,&diff1[0],1,&out_d2[0],1,&out_d2[0],1);
                    Vmath::Vvtvp (nqtot,&gmat[8][0],1,&diff2[0],1,&out_d2[0],1,&out_d2[0],1);
                }
            }
            else // regular geometry
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul (nqtot,gmat[0][0],&diff0[0],1,&out_d0[0],1);
                    Blas::Daxpy (nqtot,gmat[1][0],&diff1[0],1,&out_d0[0],1);
                    Blas::Daxpy (nqtot,gmat[2][0],&diff2[0],1,&out_d0[0],1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Smul (nqtot,gmat[3][0],&diff0[0],1,&out_d1[0],1);
                    Blas::Daxpy (nqtot,gmat[4][0],&diff1[0],1,&out_d1[0],1);
                    Blas::Daxpy (nqtot,gmat[5][0],&diff2[0],1,&out_d1[0],1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Smul (nqtot,gmat[6][0],&diff0[0],1,&out_d2[0],1);
                    Blas::Daxpy (nqtot,gmat[7][0],&diff1[0],1,&out_d2[0],1);
                    Blas::Daxpy (nqtot,gmat[8][0],&diff2[0],1,&out_d2[0],1);
                }
            }
        }

        //---------------------------------------
        // Transforms
        //---------------------------------------
        
        /** 
         * \brief Forward transform from physical quadrature space stored in
         * \a inarray and evaluate the expansion coefficients and store in \a
         * (this)->m_coeffs
         * 
         * Inputs:\n
         *
         * - \a inarray: array of physical quadrature points to be transformed
         *
         * Outputs:\n
         *
         * - (this)->_coeffs: updated array of expansion coefficients.
         */
        void PrismExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                        Array<OneD,       NekDouble>& outarray)
        {
            if(m_base[0]->Collocation() && 
               m_base[1]->Collocation() && 
               m_base[2]->Collocation())
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&m_coeffs[0],1);
            }
            else
            {
                v_IProductWRTBase(inarray, outarray);

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


        //---------------------------------------
        // Inner product functions
        //---------------------------------------

        /**
         * \brief Calculate the inner product of inarray with respect to the
         * basis B=base0*base1*base2 and put into outarray:
         *            
         * \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
         * \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2} \psi_{p}^{a}
         * (\bar \eta_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{pr}^{b} (\xi_{3k})
         * w_i w_j w_k u(\bar \eta_{1,i} \xi_{2,j} \xi_{3,k}) J_{i,j,k}\\ & =
         * & \sum_{i=0}^{nq_0} \psi_p^a(\bar \eta_{1,i}) \sum_{j=0}^{nq_1}
         * \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2} \psi_{pr}^b u(\bar
         * \eta_{1i},\xi_{2j},\xi_{3k}) J_{i,j,k} \end{array} \f$ \n
         *  
         * where
         *  
         * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\bar \eta_1)
         * \psi_{q}^a (\xi_2) \psi_{pr}^b (\xi_3) \f$ \n
         * 
         * which can be implemented as \n \f$f_{pr} (\xi_{3k}) =
         * \sum_{k=0}^{nq_3} \psi_{pr}^b u(\bar \eta_{1i},\xi_{2j},\xi_{3k})
         * J_{i,j,k} = {\bf B_3 U} \f$ \n \f$ g_{q} (\xi_{3k}) =
         * \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{pr} (\xi_{3k}) = {\bf
         * B_2 F} \f$ \n \f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0}
         * \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k}) = {\bf B_1 G} \f$
         */
        void PrismExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTBase_SumFac(inarray, outarray);
        }

        void PrismExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            const int nquad0 = m_base[0]->GetNumPoints();
            const int nquad1 = m_base[1]->GetNumPoints();
            const int nquad2 = m_base[2]->GetNumPoints();
            const int order0 = m_base[0]->GetNumModes();
            const int order1 = m_base[1]->GetNumModes();

            Array<OneD, NekDouble> tmp(nquad0*nquad1*nquad2);
            Array<OneD, NekDouble> wsp(order0*nquad2*(nquad1+order1));

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
         * In the tetrahedral element, we must also incorporate a second set
         * of geometric factors which incorporate the collapsed co-ordinate
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
        void PrismExp::v_IProductWRTDerivBase(
            const int                           dir, 
            const Array<OneD, const NekDouble> &inarray, 
                  Array<OneD,       NekDouble> &outarray)
        {
            v_IProductWRTDerivBase_SumFac(dir, inarray, outarray);
        }

        void PrismExp::v_IProductWRTDerivBase_SumFac(
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
            int i;
            
            const Array<OneD, const NekDouble> &z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble> &z2 = m_base[2]->GetZ();
            
            Array<OneD, NekDouble> gfac0(nquad0   );
            Array<OneD, NekDouble> gfac2(nquad2   );
            Array<OneD, NekDouble> tmp1 (nqtot    );
            Array<OneD, NekDouble> tmp2 (nqtot    );
            Array<OneD, NekDouble> tmp3 (nqtot    );
            Array<OneD, NekDouble> tmp4 (nqtot    );
            Array<OneD, NekDouble> tmp5 (nqtot    );
            Array<OneD, NekDouble> tmp6 (m_ncoeffs);
            Array<OneD, NekDouble> wsp  (order0*nquad2*(nquad1+order1));

            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            MultiplyByQuadratureMetric(inarray, tmp1);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nqtot,&gmat[3*dir][0],  1,tmp1.get(),1,tmp2.get(),1);
                Vmath::Vmul(nqtot,&gmat[3*dir+1][0],1,tmp1.get(),1,tmp3.get(),1);
                Vmath::Vmul(nqtot,&gmat[3*dir+2][0],1,tmp1.get(),1,tmp4.get(),1);
            }
            else
            {
                Vmath::Smul(nqtot, gmat[3*dir][0],  tmp1.get(),1,tmp2.get(), 1);
                Vmath::Smul(nqtot, gmat[3*dir+1][0],tmp1.get(),1,tmp3.get(), 1);
                Vmath::Smul(nqtot, gmat[3*dir+2][0],tmp1.get(),1,tmp4.get(), 1);
            }
            
            // set up geometric factor: (1+z0)/2
            for (i = 0; i < nquad0; ++i)
            {
                gfac0[i] = 0.5*(1+z0[i]);
            }
            
            // Set up geometric factor: 2/(1-z2)
            for (i = 0; i < nquad2; ++i)
            {
            	gfac2[i] = 2.0/(1-z2[i]);
            }
            
            const int nq01 = nquad0*nquad1;
            
            for (i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nq01,gfac2[i],&tmp2[0]+i*nq01,1,&tmp2[0]+i*nq01,1);
                Vmath::Smul(nq01,gfac2[i],&tmp4[0]+i*nq01,1,&tmp5[0]+i*nq01,1);
            }
            
            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vmul(nquad0,&gfac0[0],1,&tmp5[0]+i*nquad0,1,
                            &tmp5[0]+i*nquad0,1);
            }

            Vmath::Vadd(nqtot, &tmp2[0], 1, &tmp5[0], 1, &tmp2[0], 1);

            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                         m_base[1]->GetBdata (),
                                         m_base[2]->GetBdata (),
                                         tmp2,outarray,wsp,
                                         true,true,true);
            
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata (),
                                         m_base[1]->GetDbdata(),
                                         m_base[2]->GetBdata (),
                                         tmp3,tmp6,wsp,
                                         true,true,true);

            Vmath::Vadd(m_ncoeffs, tmp6, 1, outarray, 1, outarray, 1);
            
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata (),
                                         m_base[1]->GetBdata (),
                                         m_base[2]->GetDbdata(),
                                         tmp4,tmp6,wsp,
                                         true,true,true);

            Vmath::Vadd(m_ncoeffs, tmp6, 1, outarray, 1, outarray, 1);
        }
        
        //---------------------------------------
        // Evaluation functions
        //---------------------------------------
        
        void PrismExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                   Array<OneD, NekDouble> &coords_1,
                                   Array<OneD, NekDouble> &coords_2)
        {
            LibUtilities::BasisSharedPtr CBasis0;
            LibUtilities::BasisSharedPtr CBasis1;
            LibUtilities::BasisSharedPtr CBasis2;
            Array<OneD,NekDouble>  x;
            
            ASSERTL0(m_geom, "m_geom not defined");
            
            // get physical points defined in Geom
            m_geom->FillGeom();  

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
                        //Blas::Dcopy(GetTotPoints(), x, 1, coords_2, 1);
                        Blas::Dcopy(m_base[0]->GetNumPoints()*
                                    m_base[1]->GetNumPoints()*
                                    m_base[2]->GetNumPoints(),
                                    x, 1, coords_2, 1);
                    }
                    else // LibUtilities::Interpolate to Expansion point distribution
                    {
                        LibUtilities::Interp3D(CBasis0->GetPointsKey(),
                                               CBasis1->GetPointsKey(), 
                                               CBasis2->GetPointsKey(), 
                                               &(m_geom->UpdatePhys(2))[0],
                                               m_base[0]->GetPointsKey(), 
                                               m_base[1]->GetPointsKey(), 
                                               m_base[2]->GetPointsKey(), 
                                               &coords_2[0]);
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
                        Blas::Dcopy(GetTotPoints(), x, 1, coords_1, 1);
                    }
                    else // Interpolate to Expansion point distribution
                    {
                        LibUtilities::Interp3D(CBasis0->GetPointsKey(), 
                                               CBasis1->GetPointsKey(), 
                                               CBasis2->GetPointsKey(), 
                                               &(m_geom->UpdatePhys(1))[0],
                                               m_base[0]->GetPointsKey(), 
                                               m_base[1]->GetPointsKey(), 
                                               m_base[2]->GetPointsKey(), 
                                               &coords_1[0]);
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
                        Blas::Dcopy(GetTotPoints(), x, 1, coords_0, 1);
                    }
                    else // Interpolate to Expansion point distribution
                    {
                        LibUtilities::Interp3D(CBasis0->GetPointsKey(), 
                                               CBasis1->GetPointsKey(), 
                                               CBasis2->GetPointsKey(), 
                                               &(m_geom->UpdatePhys(0))[0],
                                               m_base[0]->GetPointsKey(),
                                               m_base[1]->GetPointsKey(),
                                               m_base[2]->GetPointsKey(),
                                               &coords_0[0]);
                    }
                    break;
                default:
                    ASSERTL0(false,"Number of dimensions are greater than 3");
                    break;
            }
        }

        /**
         * @brief Get the coordinates #coords at the local coordinates
         * #Lcoords.
         */
        void PrismExp::v_GetCoord(const Array<OneD, const NekDouble>& Lcoords, 
                                        Array<OneD,       NekDouble>& coords)
        {
            int i;
            
            ASSERTL1(Lcoords[0] <= -1.0 && Lcoords[0] >= 1.0 && 
                     Lcoords[1] <= -1.0 && Lcoords[1] >= 1.0 &&
                     Lcoords[2] <= -1.0 && Lcoords[2] >= 1.0,
                     "Local coordinates are not in region [-1,1]");

            m_geom->FillGeom();
            
            for(i = 0; i < m_geom->GetCoordim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }

        NekDouble PrismExp::v_PhysEvaluate(const Array<OneD, const NekDouble>& coord)
        {
            return PhysEvaluate(coord, m_phys);
        }

        NekDouble PrismExp::v_PhysEvaluate(const Array<OneD, const NekDouble>& coord, 
                                           const Array<OneD, const NekDouble>& physvals)
        {
            Array<OneD, NekDouble> Lcoord(3);

            ASSERTL0(m_geom,"m_geom not defined");
	
            m_geom->GetLocCoords(coord, Lcoord);

            return StdPrismExp::v_PhysEvaluate(Lcoord, physvals);
        }


        //---------------------------------------
        // Helper functions
        //---------------------------------------

        const SpatialDomains::GeomFactorsSharedPtr& PrismExp::v_GetMetricInfo() const
        {
            return m_metricinfo;
        }

        const SpatialDomains::GeometrySharedPtr PrismExp::v_GetGeom() const
        {
            return m_geom;
        }

        const SpatialDomains::Geometry3DSharedPtr& PrismExp::v_GetGeom3D() const
        {
            return m_geom;
        }
        
        int PrismExp::v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }

        StdRegions::Orientation PrismExp::v_GetFaceOrient(int face)
        {
            return m_geom->GetFaceOrient(face);
        }

        ///Returns the physical values at the quadrature points of a face
        void PrismExp::v_GetFacePhysVals(
            const int                                face,
            const StdRegions::StdExpansionSharedPtr &FaceExp,
            const Array<OneD, const NekDouble>      &inarray,
                  Array<OneD,       NekDouble>      &outarray,
            StdRegions::Orientation                  orient)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();

            Array<OneD,NekDouble>       o_tmp(nquad0*nquad1*nquad2);
            
            if (orient == StdRegions::eNoOrientation)
            {
                orient = GetFaceOrient(face);
            }

            switch(face)
            {
            case 0:
                if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad0*nquad1,&(inarray[0]),1,&(outarray[0]),1);
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
                {
                    //Direction A negative and B positive
                    for (int j=0; j<nquad1; j++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0-1)+j*nquad0,-1,&(outarray[0])+(j*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2)
                {
                    //Direction A positive and B negative
                    for (int j=0; j<nquad1; j++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+nquad0*(nquad1-1-j),1,&(outarray[0])+(j*nquad0),1);
                    }
                } 
                else if(orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2)
                {
                    //Direction A negative and B negative
                    for(int j=0; j<nquad1; j++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1-1-j*nquad0),-1,&(outarray[0])+(j*nquad0),1);
                    }
                }
		else if(orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A and B positive
		    for (int i=0; i<nquad0; i++)
                    {
                        Vmath::Vcopy(nquad1,&(inarray[0])+i,nquad0,&(outarray[0])+(i*nquad1),1);
                    }
		}
		else if(orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A positive and B negative
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad1,&(inarray[0])+(nquad0-1-i),nquad0,&(outarray[0])+(i*nquad1),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A negative and B positive
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad1,&(inarray[0])+i+nquad0*(nquad1-1),-nquad0,&(outarray[0])+(i*nquad1),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A and B negative
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad1,&(inarray[0])+(nquad0*nquad1-1-i),-nquad0,&(outarray[0])+(i*nquad1),1);
                    }
		} 
                o_tmp=outarray;
                //interpolate
                LibUtilities::Interp2D(m_base[0]->GetPointsKey(), m_base[1]->GetPointsKey(), o_tmp,
                             FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),outarray);
                break;
	    case 1:
                if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Direction A and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1*k),1,&(outarray[0])+(k*nquad0),1);
                    }
                    o_tmp=outarray;
                }
                else
                {
                    //Direction A negative and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0-1)+(nquad0*nquad1*k),-1,&(outarray[0])+(k*nquad0),1);
                    }
                    o_tmp=outarray;
                }

                //interpolate
                LibUtilities::Interp2D(m_base[0]->GetPointsKey(), m_base[2]->GetPointsKey(), o_tmp,
                             FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),outarray);
                break;
            case 2:
	        if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad0*nquad2,&(inarray[0])+(nquad0-1),
				 nquad0,&(outarray[0]),1);
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
                {
                    //Direction A negative and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1-1)+(k*nquad0*nquad1),
                                     -nquad0,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2)
                {
                    //Direction A positive and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0-1)+(nquad0*nquad1*(nquad2-1-k)),
                                     nquad0,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2)
                {
                    //Direction A negative and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1-1)+(nquad0*nquad1*(nquad2-1-k)),
                                     -nquad0,&(outarray[0])+(k*nquad0),1);
                    }
                }
		else if(orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A and B positive
		    for (int j=0; j<nquad1; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0-1)+(j*nquad0),
                                     nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		}
		else if(orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A positive and B negative
		    for (int j=0; j<nquad0; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0*nquad1-1-j*nquad0),
                                     nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A negative and B positive
		    for (int j=0; j<nquad0; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+nquad0*nquad1*(nquad2-1)+nquad0+j*nquad0,
                                     -nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A and B negative
		    for (int j=0; j<nquad0; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0*nquad1*nquad2-1-j*nquad0),
                                     -nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		}
                o_tmp=outarray;
                //interpolate
                LibUtilities::Interp2D(m_base[1]->GetPointsKey(), m_base[2]->GetPointsKey(), o_tmp,
                             FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),outarray);
                break;
		case 3:
	        if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Directions A and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*(nquad1-1))+(k*nquad0*nquad1),
                                     1,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else
                {
                    //Direction A negative and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1-1)+(k*nquad0*nquad1),
                                     -1,&(outarray[0])+(k*nquad0),1);
                    }
		}
                o_tmp=outarray;
                //interpolate
                LibUtilities::Interp2D(m_base[0]->GetPointsKey(), m_base[2]->GetPointsKey(), o_tmp,
                             FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),outarray);

                break;
            case 4:
                if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad1*nquad2,&(inarray[0]),nquad0,&(outarray[0]),1);
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
                {
                    //Direction A negative and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad1,&(inarray[0])+nquad0*(nquad1-1)+(k*nquad0*nquad1),
                                     -nquad0,&(outarray[0])+(k*nquad1),1);
                    }
                }
                else if(orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2)
                {
                    //Direction A positive and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad1,&(inarray[0])+(nquad0*nquad1*(nquad2-1-k)),
                                     nquad0,&(outarray[0])+(k*nquad1),1);
                    }
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2)
                {
                    //Direction A negative and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad1,&(inarray[0])+nquad0*(nquad1-1)+(nquad0*nquad1*(nquad2-1-k)),
                                     -nquad0,&(outarray[0])+(k*nquad1),1);
                    }
                }
		else if(orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A and B positive
		    for (int j=0; j<nquad1; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+j*nquad0,nquad0*nquad1,
                                     &(outarray[0])+(j*nquad2),1);
                    }
		}
		else if(orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A positive and B negative
		    for (int j=0; j<nquad1; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0*(nquad1-1)-j*nquad0),
                                     nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A negative and B positive
		    for (int j=0; j<nquad1; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+nquad0*nquad1*(nquad2-1)+j*nquad0,
                                     -nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A and B negative
		    for (int j=0; j<nquad1; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0*(nquad1*nquad2-1)-j*nquad0),
                                     -nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		} 
                o_tmp=outarray;
                //interpolate
                LibUtilities::Interp2D(m_base[1]->GetPointsKey(), m_base[2]->GetPointsKey(), o_tmp,
                             FaceExp->GetBasis(0)->GetPointsKey(),FaceExp->GetBasis(1)->GetPointsKey(),outarray);
                break;
            default:
                ASSERTL0(false,"face value (> 4) is out of range");
                break;
	    }
	}

        void PrismExp::v_ComputeFaceNormal(const int face)
        {
            const SpatialDomains::GeomFactorsSharedPtr &geomFactors = 
                GetGeom()->GetMetricInfo();
            SpatialDomains::GeomType type            = geomFactors->GetGtype();
            const Array<TwoD, const NekDouble> &gmat = geomFactors->GetGmat();
            const Array<OneD, const NekDouble> &jac  = geomFactors->GetJac();
            
            // Number of quadrature points in face expansion.
            int nq        = m_base[0]->GetNumPoints()*m_base[0]->GetNumPoints();
            int vCoordDim = GetCoordim();
            int i;
            
            m_faceNormals[face] = Array<OneD, Array<OneD, NekDouble> >(vCoordDim);
            Array<OneD, Array<OneD, NekDouble> > &normal = m_faceNormals[face];
            for (i = 0; i < vCoordDim; ++i)
            {
                normal[i] = Array<OneD, NekDouble>(nq);
            }
            
            // Regular geometry case
            if (type == SpatialDomains::eRegular      || 
                type == SpatialDomains::eMovingRegular)
            {
                NekDouble fac;
                // Set up normals
                switch(face)
                {
                    case 0:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,-gmat[3*i+2][0],normal[i],1);
                        }
                        break;
                    }
                    case 1:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,-gmat[3*i+1][0],normal[i],1);
                        }
                        break;
                    }
                    case 2:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,gmat[3*i][0]+gmat[3*i+2][0],normal[i],1);
                        }
                        break;
                    }
                    case 3:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,gmat[3*i+1][0],normal[i],1);
                        }
                        break;
                    }
                    case 4:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nq,-gmat[3*i][0],normal[i],1);
                        }
                        break;
                    }
                    default:
                        ASSERTL0(false,"face is out of range (face < 4)");
                }

                // Normalise resulting vector.
                fac = 0.0;
                for(i = 0; i < vCoordDim; ++i)
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
                // Set up deformed normals.
	        int j, k;

                int nquad0 = geomFactors->GetPointsKey(0).GetNumPoints();
                int nquad1 = geomFactors->GetPointsKey(1).GetNumPoints();
		int nquad2 = geomFactors->GetPointsKey(2).GetNumPoints();
		int nqtot;
                
                // Determine number of quadrature points on the face.
                if (face == 0)
                {
                    nqtot = nquad0*nquad1;
                }
                else if (face == 1 || face == 3)
                {
                    nqtot = nquad0*nquad2;
                }
                else
                {
                    nqtot = nquad1*nquad2;
                }
                
                LibUtilities::PointsKey points0;
		LibUtilities::PointsKey points1;
                
                Array<OneD, NekDouble> work   (nq,             0.0);
                Array<OneD, NekDouble> normals(vCoordDim*nqtot,0.0);
                
                // Extract Jacobian along face and recover local derivatives
                // (dx/dr) for polynomial interpolation by multiplying m_gmat by
                // jacobian
                switch(face)
	        {
                    case 0:
                    {
                        for(j = 0; j < nquad0*nquad1; ++j)
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
                        int tmp = nquad0*nquad1;
                        
                        for (j = 0; j < nquad0; ++j)
                        {
                            for(k = 0; k < nquad2; ++k)
                            {
                                normals[j+k*nquad0]          = 
                                    -gmat[1][j+tmp*k]*jac[j+tmp*k];
                                normals[nqtot+j+k*nquad0]    = 
                                    -gmat[4][j+tmp*k]*jac[j+tmp*k];
                                normals[2*nqtot+j+k*nquad0]  = 
                                    -gmat[7][j+tmp*k]*jac[j+tmp*k];
                            } 
                        }

                        points0 = geomFactors->GetPointsKey(0);
                        points1 = geomFactors->GetPointsKey(2);
                        break;
                    }
                    
                    case 2:
                    {
                        for (j = 0; j < nquad1; ++j)
                        {
                            for(k = 0; k < nquad2; ++k)
                            {
                                normals[j+k*nquad0]  = (gmat[0][nquad0-1+nquad0*j+nquad0*nquad1*k]
                                                        +gmat[2][nquad0-1+nquad0*j+nquad0*nquad1*k])*jac[nquad0-1+nquad0*j+nquad0*nquad1*k];
                                normals[nqtot+j+k*nquad0]  = (gmat[3][nquad0-1+nquad0*j+nquad0*nquad1*k]
                                                              +gmat[5][nquad0-1+nquad0*j+nquad0*nquad1*k])*jac[nquad0-1+nquad0*j+nquad0*nquad1*k];
                                normals[2*nqtot+j+k*nquad0]  = (gmat[6][nquad0-1+nquad0*j+nquad0*nquad1*k]
                                                                +gmat[8][nquad0-1+nquad0*j+nquad0*nquad1*k])*jac[nquad0-1+nquad0*j+nquad0*nquad1*k];
                            } 
                        }
                        
                        points0 = geomFactors->GetPointsKey(1);
                        points1 = geomFactors->GetPointsKey(2);
                        break;
                    }
                        
                    case 3:
                    {
                        for (j=0; j< nquad0; ++j)
                        {
                            for(k=0; k<nquad2; ++k)
                            {
                                normals[j+k*nquad0]  = gmat[1][nquad0*(nquad1-1)+j+nquad0*nquad1*k]*jac[nquad0*(nquad1-1)+j+nquad0*nquad1*k];
                                normals[nqtot+j+k*nquad0]  = gmat[4][nquad0*(nquad1-1)+j+nquad0*nquad1*k]*jac[nquad0*(nquad1-1)+j+nquad0*nquad1*k];
                                normals[2*nqtot+j+k*nquad0]  = gmat[7][nquad0*(nquad1-1)+j+nquad0*nquad1*k]*jac[nquad0*(nquad1-1)+j+nquad0*nquad1*k];
                            } 
                        }
                        
                        points0 = geomFactors->GetPointsKey(0);
                        points1 = geomFactors->GetPointsKey(2);
                        break;
                    }
                    
                    case 4:
                    {
                        for (j=0; j< nquad0; ++j)
                        {
                            for(k=0; k<nquad2; ++k)
                            {
                                normals[j+k*nquad0]  = -gmat[0][j*nquad0+nquad0*nquad1*k]*jac[j*nquad0+nquad0*nquad1*k];
                                normals[nqtot+j+k*nquad0]  = -gmat[3][j*nquad0+nquad0*nquad1*k]*jac[j*nquad0+nquad0*nquad1*k];
                                normals[2*nqtot+j+k*nquad0]  = -gmat[6][j*nquad0+nquad0*nquad1*k]*jac[j*nquad0+nquad0*nquad1*k];
                            } 
                        }
                        
                        points0 = geomFactors->GetPointsKey(1);
                        points1 = geomFactors->GetPointsKey(2);
                        break;
                    }
                        
                    default:
                        ASSERTL0(false,"face is out of range (face < 4)");
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


        void PrismExp::v_WriteToFile(std::ofstream &outfile, 
                                     OutputFormat   format, 
                                     const bool     dumpVar, 
                                     std::string    var)
        {
            if(format==eTecplot)
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
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }

        void PrismExp::v_MassMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
        {
            StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void PrismExp::v_LaplacianMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
        {
            PrismExp::LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void PrismExp::v_LaplacianMatrixOp(
            const int                           k1, 
            const int                           k2,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
        {
            StdExpansion::LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,mkey);
        }

        void PrismExp::v_LaplacianMatrixOp_MatFree(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
        {
            if(mkey.GetNVarCoeff() == 0)
            {
                // This implementation is only valid when there are no
                // coefficients associated to the Laplacian operator
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
                    int nmodes0 = m_base[0]->GetNumModes ();
                    int nmodes1 = m_base[1]->GetNumModes ();
                    int nqtot   = nquad0*nquad1*nquad2;
                    
                    const Array<OneD, const NekDouble>& base0 = m_base[0]->GetBdata ();
                    const Array<OneD, const NekDouble>& base1 = m_base[1]->GetBdata ();
                    const Array<OneD, const NekDouble>& base2 = m_base[2]->GetBdata ();

                    Array<OneD,NekDouble> wsp (nquad2*nmodes0*(nquad1+nmodes1));
                    Array<OneD,NekDouble> wsp1(nqtot);
                    
                    // Backwards transform to obtain u = B * u_hat.
                    BwdTrans_SumFacKernel   (base0,base1,base2,inarray,wsp1,wsp,true,true,true);
                    LaplacianMatrixOp_Kernel(wsp1, outarray, wsp);
                }
            }
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(
                    inarray,outarray,mkey);
            }
        }

        void PrismExp::v_HelmholtzMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
        {
            PrismExp::v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void PrismExp::v_GeneralMatrixOp_MatOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
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
        
        void PrismExp::v_HelmholtzMatrixOp_MatFree(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
        {
            if(m_metricinfo->IsUsingLaplMetrics())
            {
                ASSERTL0(false,"Finish implementing PrismExp Helmholtz for Lapl Metrics");
            }
            else
            {
                int nquad0  = m_base[0]->GetNumPoints();
                int nquad1  = m_base[1]->GetNumPoints();
                int nquad2  = m_base[2]->GetNumPoints();
                int nmodes0 = m_base[0]->GetNumModes ();
                int nmodes1 = m_base[1]->GetNumModes ();
                int nqtot   = nquad0*nquad1*nquad2;
                
                const Array<OneD, const NekDouble>& base0 = m_base[0]->GetBdata ();
                const Array<OneD, const NekDouble>& base1 = m_base[1]->GetBdata ();
                const Array<OneD, const NekDouble>& base2 = m_base[2]->GetBdata ();
                
                Array<OneD,NekDouble> wsp (nquad2*nmodes0*(nquad1+nmodes1));
                Array<OneD,NekDouble> wsp0(nqtot);
                Array<OneD,NekDouble> wsp1(nqtot);

                NekDouble lambda  = mkey.GetConstFactor(StdRegions::eFactorLambda);
                
                // MASS MATRIX OPERATION
                // The following is being calculated:
                // wsp0     = B   * u_hat = u
                // wsp1     = W   * wsp0
                // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat
                BwdTrans_SumFacKernel       (base0,base1,base2,inarray,
                                             wsp0,wsp,true,true,true);
                MultiplyByQuadratureMetric  (wsp0,wsp1);
                IProductWRTBase_SumFacKernel(base0,base1,base2,wsp1,
                                             outarray,wsp,true,true,true);
                LaplacianMatrixOp_Kernel    (wsp0,wsp1,wsp);
                
                // outarray = lambda * outarray + wsp1
                //          = (lambda * M + L ) * u_hat
                Vmath::Svtvp(m_ncoeffs,lambda,&outarray[0],1,&wsp1[0],1,
                             &outarray[0],1); 
           }
        }    
        
        //---------------------------------------
        // Matrix creation functions
        //---------------------------------------
        
        DNekMatSharedPtr PrismExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
                    returnval = StdPrismExp::v_GenMatrix(mkey);
                    break;
            }
            
            return returnval;            
        }
        
        DNekMatSharedPtr PrismExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
            StdRegions::StdPrismExpSharedPtr tmp = 
                MemoryManager<StdPrismExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);
            
            return tmp->GetStdMatrix(mkey); 
        }

        DNekScalMatSharedPtr PrismExp::v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }

        DNekScalBlkMatSharedPtr PrismExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
        }

        DNekScalMatSharedPtr PrismExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,
                     "Geometric information is not set up");
            
            switch(mkey.GetMatrixType())
            {
                case StdRegions::eMass:
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
                    break;
                }
                
                case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetExpansionType(),*this);
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
                    break;
                }
                
                case StdRegions::eWeakDeriv0:
                case StdRegions::eWeakDeriv1:
                case StdRegions::eWeakDeriv2:
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

                        (*WeakDeriv) = gmat[3*dir  ][0]*deriv0
                                     + gmat[3*dir+1][0]*deriv1
                                     + gmat[3*dir+2][0]*deriv2;

                        returnval = MemoryManager<DNekScalMat>
                            ::AllocateSharedPtr(jac,WeakDeriv);
                    }
                    break;
                }
                
                case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
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
                    break;
                }
                
                case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetExpansionType(), *this);    
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;
                    
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);            
                    break;
                }

                case StdRegions::eHybridDGHelmholtz:
                case StdRegions::eHybridDGLamToU:
                case StdRegions::eHybridDGLamToQ0:
                case StdRegions::eHybridDGLamToQ1:
                case StdRegions::eHybridDGHelmBndLam:
                {
                    NekDouble one    = 1.0;
                    
                    DNekMatSharedPtr mat = GenMatrix(mkey);
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    
                    break;
                }
                
                case StdRegions::eInvHybridDGHelmholtz:
                {
                    NekDouble one = 1.0;
                    
                    MatrixKey hkey(StdRegions::eHybridDGHelmholtz, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
//                    StdRegions::StdMatrixKey hkey(StdRegions::eHybridDGHelmholtz,
//                                                  DetExpansionType(),*this,
//                                                  mkey.GetConstant(0),
//                                                  mkey.GetConstant(1));
                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    break;
                }
                
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
                    break;
                }
                default:
                {
                    NekDouble        one = 1.0;
                    DNekMatSharedPtr mat = GenMatrix(mkey);

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
            }
                
            return returnval;
        }

        DNekScalBlkMatSharedPtr PrismExp::CreateStaticCondMatrix(const MatrixKey &mkey)
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

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eLaplacian:
            case StdRegions::eHelmholtz: // special case since Helmholtz not defined in StdRegions

                // use Deformed case for both regular and deformed geometries
                factor = 1.0;
                goto UseLocRegionsMatrix;
                break;
            default:
                if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
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
        
        void PrismExp::MultiplyByQuadratureMetric(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
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

                StdPrismExp::MultiplyByQuadratureMetric(outarray, outarray);
            }
        }
        
        /**
         * @brief Calculate the Laplacian multiplication in a matrix-free
         * manner.
         * 
         * This function is the kernel of the Laplacian matrix-free operator,
         * and is used in #v_HelmholtzMatrixOp_MatFree to determine the effect
         * of the Helmholtz operator in a similar fashion.
         * 
         * The majority of the calculation is precisely the same as in the
         * hexahedral expansion; however the collapsed co-ordinate system must
         * be taken into account when constructing the geometric factors. How
         * this is done is detailed more exactly in the tetrahedral expansion.
         * On entry to this function, the input #inarray must be in its
         * backwards-transformed state (i.e. \f$\mathbf{u} =
         * \mathbf{B}\hat{\mathbf{u}}\f$). The output is in coefficient space.
         * 
         * @see %TetExp::v_HelmholtzMatrixOp_MatFree
         */
        void PrismExp::LaplacianMatrixOp_Kernel(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
                  Array<OneD,       NekDouble> &wsp)
        {
            int nquad0  = m_base[0]->GetNumPoints();
            int nquad1  = m_base[1]->GetNumPoints();
            int nquad2  = m_base[2]->GetNumPoints();
            int nqtot   = nquad0*nquad1*nquad2;
            int i;
            
            // Set up temporary storage.
            Array<OneD,NekDouble> alloc(11*nqtot,0.0);
            Array<OneD,NekDouble> wsp1 (alloc        );  // TensorDeriv 1
            Array<OneD,NekDouble> wsp2 (alloc+ 1*nqtot); // TensorDeriv 2
            Array<OneD,NekDouble> wsp3 (alloc+ 2*nqtot); // TensorDeriv 3
            Array<OneD,NekDouble> g0   (alloc+ 3*nqtot); // g0
            Array<OneD,NekDouble> g1   (alloc+ 4*nqtot); // g1
            Array<OneD,NekDouble> g2   (alloc+ 5*nqtot); // g2
            Array<OneD,NekDouble> g3   (alloc+ 6*nqtot); // g3
            Array<OneD,NekDouble> g4   (alloc+ 7*nqtot); // g4
            Array<OneD,NekDouble> g5   (alloc+ 8*nqtot); // g5
            Array<OneD,NekDouble> h0   (alloc+ 3*nqtot); // h0   == g0
            Array<OneD,NekDouble> h1   (alloc+ 6*nqtot); // h1   == g3
            Array<OneD,NekDouble> wsp4 (alloc+ 4*nqtot); // wsp4 == g1
            Array<OneD,NekDouble> wsp5 (alloc+ 5*nqtot); // wsp5 == g2
            Array<OneD,NekDouble> wsp6 (alloc+ 8*nqtot); // wsp6 == g5
            Array<OneD,NekDouble> wsp7 (alloc+ 3*nqtot); // wsp7 == g0
            Array<OneD,NekDouble> wsp8 (alloc+ 9*nqtot); // wsp8
            Array<OneD,NekDouble> wsp9 (alloc+10*nqtot); // wsp9
                    
            const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble>& base2  = m_base[2]->GetBdata();
            const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase2 = m_base[2]->GetDbdata();
            
            // Step 1. LAPLACIAN MATRIX OPERATION
            // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
            // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
            // wsp3 = du_dxi3 = D_xi3 * wsp0 = D_xi3 * u
            StdExpansion3D::PhysTensorDeriv(inarray,wsp1,wsp2,wsp3);

            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();
            const Array<OneD, const NekDouble>& z0   = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z2   = m_base[2]->GetZ();
            
            // Step 2. Calculate the metric terms of the collapsed
            // coordinate transformation (Spencer's book P152)
            for (i = 0; i < nquad2; ++i)
            {
                Vmath::Fill(nquad0*nquad1, 2.0/(1.0-z2[i]), &h0[0]+i*nquad0*nquad1,1);
                Vmath::Fill(nquad0*nquad1, 2.0/(1.0-z2[i]), &h1[0]+i*nquad0*nquad1,1);
            }
            for (i = 0; i < nquad0; i++)
            {
                Blas::Dscal(nquad1*nquad2, 0.5*(1+z0[i]), &h1[0]+i, nquad0);
            }
            
            // Step 3. Construct combined metric terms for physical space to
            // collapsed coordinate system.  Order of construction optimised
            // to minimise temporary storage
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // wsp4 = d eta_1/d x_1
                Vmath::Vvtvvtp(nqtot, &gmat[0][0], 1, &h0[0], 1, &gmat[2][0], 1, &h1[0], 1, &wsp4[0], 1);
                // wsp5 = d eta_2/d x_1
                Vmath::Vvtvvtp(nqtot, &gmat[3][0], 1, &h0[0], 1, &gmat[5][0], 1, &h1[0], 1, &wsp5[0], 1);
                // wsp6 = d eta_3/d x_1
                Vmath::Vvtvvtp(nqtot, &gmat[6][0], 1, &h0[0], 1, &gmat[8][0], 1, &h1[0], 1, &wsp6[0], 1);
                
                // g0 (overwrites h0)
                Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp4[0], 1, &wsp5[0], 1, &wsp5[0], 1, &g0[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp6[0], 1, &wsp6[0], 1, &g0[0],   1, &g0[0],   1);
                
                // g3 (overwrites h1)
                Vmath::Vvtvvtp(nqtot, &gmat[1][0], 1, &wsp4[0], 1, &gmat[4][0], 1, &wsp5[0], 1, &g3[0], 1);
                Vmath::Vvtvp  (nqtot, &gmat[7][0], 1, &wsp6[0], 1, &g3[0], 1, &g3[0], 1);
                
                // g4
                Vmath::Vvtvvtp(nqtot, &gmat[2][0], 1, &wsp4[0], 1, &gmat[5][0], 1, &wsp5[0], 1, &g4[0], 1);
                Vmath::Vvtvp  (nqtot, &gmat[8][0], 1, &wsp6[0], 1, &g4[0], 1, &g4[0], 1);

                // Overwrite wsp4/5/6 with g1/2/5
                // g1
                Vmath::Vvtvvtp(nqtot, &gmat[1][0], 1, &gmat[1][0], 1, &gmat[4][0], 1, &gmat[4][0], 1, &g1[0], 1);
                Vmath::Vvtvp  (nqtot, &gmat[7][0], 1, &gmat[7][0], 1, &g1[0], 1, &g1[0], 1);
                
                // g2
                Vmath::Vvtvvtp(nqtot, &gmat[2][0], 1, &gmat[2][0], 1, &gmat[5][0], 1, &gmat[5][0], 1, &g2[0], 1);
                Vmath::Vvtvp  (nqtot, &gmat[8][0], 1, &gmat[8][0], 1, &g2[0], 1, &g2[0], 1);
                
                // g5
                Vmath::Vvtvvtp(nqtot, &gmat[1][0], 1, &gmat[2][0], 1, &gmat[4][0], 1, &gmat[5][0], 1, &g5[0], 1);
                Vmath::Vvtvp  (nqtot, &gmat[7][0], 1, &gmat[8][0], 1, &g5[0], 1, &g5[0], 1);
            }
            else
            {
                // wsp4 = d eta_1/d x_1
                Vmath::Svtsvtp(nqtot, gmat[0][0], &h0[0], 1, gmat[2][0], &h1[0], 1, &wsp4[0], 1);
                // wsp5 = d eta_2/d x_1
                Vmath::Svtsvtp(nqtot, gmat[3][0], &h0[0], 1, gmat[5][0], &h1[0], 1, &wsp5[0], 1);
                // wsp6 = d eta_3/d x_1
                Vmath::Svtsvtp(nqtot, gmat[6][0], &h0[0], 1, gmat[8][0], &h1[0], 1, &wsp6[0], 1);
                
                // g0 (overwrites h0)
                Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp4[0], 1, &wsp5[0], 1, &wsp5[0], 1, &g0[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp6[0], 1, &wsp6[0], 1, &g0[0],   1, &g0[0],   1);
                
                // g3 (overwrites h1)
                Vmath::Svtsvtp(nqtot, gmat[1][0], &wsp4[0], 1, gmat[4][0], &wsp5[0], 1, &g3[0], 1);
                Vmath::Svtvp  (nqtot, gmat[7][0], &wsp6[0], 1, &g3[0], 1, &g3[0], 1);
                
                // g4
                Vmath::Svtsvtp(nqtot, gmat[2][0], &wsp4[0], 1, gmat[5][0], &wsp5[0], 1, &g4[0], 1);
                Vmath::Svtvp  (nqtot, gmat[8][0], &wsp6[0], 1, &g4[0], 1, &g4[0], 1);
                
                // Overwrite wsp4/5/6 with g1/2/5
                // g1
                Vmath::Fill(nqtot, gmat[1][0]*gmat[1][0] + gmat[4][0]*gmat[4][0] + gmat[7][0]*gmat[7][0], &g1[0], 1);
                        
                // g2
                Vmath::Fill(nqtot, gmat[2][0]*gmat[2][0] + gmat[5][0]*gmat[5][0] + gmat[8][0]*gmat[8][0], &g2[0], 1);
                
                // g5
                Vmath::Fill(nqtot, gmat[1][0]*gmat[2][0] + gmat[4][0]*gmat[5][0] + gmat[7][0]*gmat[8][0], &g5[0], 1);
            }
            // Compute component derivatives into wsp7, 8, 9 (wsp7 overwrites
            // g0).
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
            
            // Perform inner product w.r.t derivative bases.
            IProductWRTBase_SumFacKernel(dbase0,base1,base2,wsp7,wsp1,    wsp,false,true,true);
            IProductWRTBase_SumFacKernel(base0,dbase1,base2,wsp8,wsp2,    wsp,true,false,true);
            IProductWRTBase_SumFacKernel(base0,base1,dbase2,wsp9,outarray,wsp,true,true,false);
            
            // Step 5.
            // Sum contributions from wsp1, wsp2 and outarray.
            Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);
            Vmath::Vadd(m_ncoeffs,wsp2.get(),1,outarray.get(),1,outarray.get(),1);
        }
    }//end of namespace
}//end of namespace
