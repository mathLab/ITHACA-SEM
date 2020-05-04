///////////////////////////////////////////////////////////////////////////////
//
// File PyrExp.cpp
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
// Description:  PyrExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/PyrExp.h>
#include <LibUtilities/Foundations/Interp.h>

using namespace std;

namespace Nektar
{
    namespace LocalRegions
    {

	PyrExp::PyrExp(const LibUtilities::BasisKey &Ba,
                       const LibUtilities::BasisKey &Bb,
                       const LibUtilities::BasisKey &Bc,
                       const SpatialDomains::PyrGeomSharedPtr &geom):
            StdExpansion  (LibUtilities::StdPyrData::getNumberOfCoefficients(
                               Ba.GetNumModes(),
                               Bb.GetNumModes(),
                               Bc.GetNumModes()),
                           3, Ba, Bb, Bc),
            StdExpansion3D(LibUtilities::StdPyrData::getNumberOfCoefficients(
                               Ba.GetNumModes(),
                               Bb.GetNumModes(),
                               Bc.GetNumModes()),
                           Ba, Bb, Bc),
            StdPyrExp     (Ba,Bb,Bc),
            Expansion     (geom),
            Expansion3D   (geom),
            m_matrixManager(
                    std::bind(&PyrExp::CreateMatrix, this, std::placeholders::_1),
                    std::string("PyrExpMatrix")),
            m_staticCondMatrixManager(
                    std::bind(&PyrExp::CreateStaticCondMatrix, this, std::placeholders::_1),
                    std::string("PyrExpStaticCondMatrix"))
        {
        }

        PyrExp::PyrExp(const PyrExp &T):
            StdExpansion  (T),
            StdExpansion3D(T),
            StdPyrExp     (T),
            Expansion     (T),
            Expansion3D   (T),
            m_matrixManager(T.m_matrixManager),
            m_staticCondMatrixManager(T.m_staticCondMatrixManager)
        {
        }

        PyrExp::~PyrExp()
        {
        }


        //----------------------------
        // Integration Methods
        //----------------------------

        /**
         * \brief Integrate the physical point list \a inarray over pyramidic
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
         * \eta_2, \eta_3) J[i,j,k] d \bar \eta_1 d \eta_2 d \eta_3\f$ \n \f$=
         * \sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1} \sum_{k=0}^{Q_3 - 1}
         * u(\bar \eta_{1i}^{0,0}, \eta_{2j}^{0,0},\eta_{3k}^{2,0})w_{i}^{0,0}
         * w_{j}^{0,0} \hat w_{k}^{2,0} \f$ \n where \f$inarray[i,j, k] =
         * u(\bar \eta_{1i},\eta_{2j}, \eta_{3k}) \f$, \n \f$\hat w_{k}^{2,0}
         * = \frac {w^{2,0}} {2} \f$ \n and \f$ J[i,j,k] \f$ is the Jacobian
         * evaluated at the quadrature point.
         */
        NekDouble PyrExp::v_Integral(const Array<OneD, const NekDouble> &inarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac(GetPointsKeys());
            Array<OneD,       NekDouble> tmp(nquad0*nquad1*nquad2);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,(NekDouble*)&inarray[0],1, &tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,(NekDouble) jac[0], (NekDouble*)&inarray[0],1,&tmp[0],1);
            }

            // call StdPyrExp version;
            return StdPyrExp::v_Integral(tmp);
        }


        //----------------------------
        // Differentiation Methods
        //----------------------------

        void PyrExp::v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                       Array<OneD,       NekDouble>& out_d0,
                                       Array<OneD,       NekDouble>& out_d1,
                                       Array<OneD,       NekDouble>& out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            Array<TwoD, const NekDouble> gmat =
                                m_metricinfo->GetDerivFactors(GetPointsKeys());
            Array<OneD,NekDouble> diff0(nquad0*nquad1*nquad2);
            Array<OneD,NekDouble> diff1(nquad0*nquad1*nquad2);
            Array<OneD,NekDouble> diff2(nquad0*nquad1*nquad2);

            StdPyrExp::v_PhysDeriv(inarray, diff0, diff1, diff2);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.size())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[0][0],1,&diff0[0],1, &out_d0[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[1][0],1,&diff1[0],1, &out_d0[0], 1,&out_d0[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[2][0],1,&diff2[0],1, &out_d0[0], 1,&out_d0[0],1);
                }

                if(out_d1.size())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[3][0],1,&diff0[0],1, &out_d1[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[4][0],1,&diff1[0],1, &out_d1[0], 1,&out_d1[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[5][0],1,&diff2[0],1, &out_d1[0], 1,&out_d1[0],1);
                }

                if(out_d2.size())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[6][0],1,&diff0[0],1, &out_d2[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[7][0],1,&diff1[0],1, &out_d2[0], 1, &out_d2[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[8][0],1,&diff2[0],1, &out_d2[0], 1, &out_d2[0],1);
                }
            }
            else // regular geometry
            {
                if(out_d0.size())
                {
                    Vmath::Smul  (nquad0*nquad1*nquad2,gmat[0][0],&diff0[0],1, &out_d0[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[1][0],&diff1[0],1, &out_d0[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[2][0],&diff2[0],1, &out_d0[0], 1);
                }

                if(out_d1.size())
                {
                    Vmath::Smul  (nquad0*nquad1*nquad2,gmat[3][0],&diff0[0],1, &out_d1[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[4][0],&diff1[0],1, &out_d1[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[5][0],&diff2[0],1, &out_d1[0], 1);
                }

                if(out_d2.size())
                {
                    Vmath::Smul  (nquad0*nquad1*nquad2,gmat[6][0],&diff0[0],1, &out_d2[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[7][0],&diff1[0],1, &out_d2[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[8][0],&diff2[0],1, &out_d2[0], 1);
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
        void PyrExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD,       NekDouble>& outarray)
        {
            if(m_base[0]->Collocation() &&
               m_base[1]->Collocation() &&
               m_base[2]->Collocation())
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&outarray[0],1);
            }
            else
            {
                v_IProductWRTBase(inarray,outarray);

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


        //---------------------------------------
        // Inner product functions
        //---------------------------------------

        /**
         * \brief Calculate the inner product of inarray with respect to the
         * basis B=base0*base1*base2 and put into outarray:
         *
         * \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
         * \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2} \psi_{p}^{a}
         * (\bar \eta_{1i}) \psi_{q}^{a} (\eta_{2j}) \psi_{pqr}^{c}
         * (\eta_{3k}) w_i w_j w_k u(\bar \eta_{1,i} \eta_{2,j} \eta_{3,k})
         * J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\bar \eta_{1,i})
         * \sum_{j=0}^{nq_1} \psi_{q}^a(\eta_{2,j}) \sum_{k=0}^{nq_2}
         * \psi_{pqr}^c u(\bar \eta_{1i},\eta_{2j},\eta_{3k}) J_{i,j,k}
         * \end{array} \f$ \n
         *
         * where
         *
         * \f$\phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\bar \eta_1)
         * \psi_{q}^a (\eta_2) \psi_{pqr}^c (\eta_3) \f$ \n
         *
         * which can be implemented as \n \f$f_{pqr} (\xi_{3k}) =
         * \sum_{k=0}^{nq_3} \psi_{pqr}^c u(\bar
         * \eta_{1i},\eta_{2j},\eta_{3k}) J_{i,j,k} = {\bf B_3 U} \f$ \n \f$
         * g_{pq} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{pqr}
         * (\xi_{3k}) = {\bf B_2 F} \f$ \n \f$ (\phi_{pqr}, u)_{\delta} =
         * \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{pq} (\xi_{3k}) = {\bf
         * B_1 G} \f$
         */

        void PyrExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTBase_SumFac(inarray, outarray);
        }

        void PyrExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
            Array<OneD,       NekDouble>& outarray,
            bool multiplybyweights)
        {
            const int nquad0 = m_base[0]->GetNumPoints();
            const int nquad1 = m_base[1]->GetNumPoints();
            const int nquad2 = m_base[2]->GetNumPoints();
            const int order0 = m_base[0]->GetNumModes();
            const int order1 = m_base[1]->GetNumModes();

            Array<OneD, NekDouble> wsp(order0*nquad2*(nquad1+order1));

            if(multiplybyweights)
            {
                Array<OneD, NekDouble> tmp(nquad0*nquad1*nquad2);

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
         * In the pyramid element, we must also incorporate a second set
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
        void PyrExp::v_IProductWRTDerivBase(
            const int                           dir,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            v_IProductWRTDerivBase_SumFac(dir, inarray, outarray);
        }

        void PyrExp::v_IProductWRTDerivBase_SumFac(
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
            const Array<OneD, const NekDouble> &z1 = m_base[1]->GetZ();
            const Array<OneD, const NekDouble> &z2 = m_base[2]->GetZ();

            Array<OneD, NekDouble> gfac0(nquad0   );
            Array<OneD, NekDouble> gfac1(nquad1   );
            Array<OneD, NekDouble> gfac2(nquad2   );
            Array<OneD, NekDouble> tmp1 (nqtot    );
            Array<OneD, NekDouble> tmp2 (nqtot    );
            Array<OneD, NekDouble> tmp3 (nqtot    );
            Array<OneD, NekDouble> tmp4 (nqtot    );
            Array<OneD, NekDouble> tmp5 (nqtot    );
            Array<OneD, NekDouble> tmp6 (m_ncoeffs);
            Array<OneD, NekDouble> wsp  (std::max(nqtot,order0*nquad2*(nquad1+order1)));

            const Array<TwoD, const NekDouble>& df =
                            m_metricinfo->GetDerivFactors(GetPointsKeys());

            MultiplyByQuadratureMetric(inarray, tmp1);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nqtot,&df[3*dir][0],  1,tmp1.get(),1,tmp2.get(),1);
                Vmath::Vmul(nqtot,&df[3*dir+1][0],1,tmp1.get(),1,tmp3.get(),1);
                Vmath::Vmul(nqtot,&df[3*dir+2][0],1,tmp1.get(),1,tmp4.get(),1);
            }
            else
            {
                Vmath::Smul(nqtot, df[3*dir][0],  tmp1.get(),1,tmp2.get(), 1);
                Vmath::Smul(nqtot, df[3*dir+1][0],tmp1.get(),1,tmp3.get(), 1);
                Vmath::Smul(nqtot, df[3*dir+2][0],tmp1.get(),1,tmp4.get(), 1);
            }

            // set up geometric factor: (1+z0)/2
            for (i = 0; i < nquad0; ++i)
            {
                gfac0[i] = 0.5*(1+z0[i]);
            }

            // set up geometric factor: (1+z1)/2
            for(i = 0; i < nquad1; ++i)
            {
                gfac1[i] = 0.5*(1+z1[i]);
            }

            // Set up geometric factor: 2/(1-z2)
            for (i = 0; i < nquad2; ++i)
            {
            	gfac2[i] = 2.0/(1-z2[i]);
            }

            const int nq01 = nquad0*nquad1;

            for (i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nq01,gfac2[i],&tmp2[0]+i*nq01,1,&tmp2[0]+i*nq01,1); // 2/(1-z2) for d/dxi_0
                Vmath::Smul(nq01,gfac2[i],&tmp3[0]+i*nq01,1,&tmp3[0]+i*nq01,1); // 2/(1-z2) for d/dxi_1
                Vmath::Smul(nq01,gfac2[i],&tmp4[0]+i*nq01,1,&tmp5[0]+i*nq01,1); // 2/(1-z2) for d/dxi_2
            }

            // (1+z0)/(1-z2) for d/d eta_0
            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vmul(nquad0,&gfac0[0],1,
                            &tmp5[0]+i*nquad0,1,
                             &wsp[0]+i*nquad0,1);
            }

            Vmath::Vadd(nqtot, &tmp2[0], 1, &wsp[0], 1, &tmp2[0], 1);

            // (1+z1)/(1-z2) for d/d eta_1
            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Smul(nquad0,gfac1[i%nquad1],
                            &tmp5[0]+i*nquad0,1,
                            &tmp5[0]+i*nquad0,1);
            }
            Vmath::Vadd(nqtot, &tmp3[0], 1, &tmp5[0], 1, &tmp3[0], 1);


            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                         m_base[1]->GetBdata (),
                                         m_base[2]->GetBdata (),
                                         tmp2,outarray,wsp,
                                         false,true,true);

            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata (),
                                         m_base[1]->GetDbdata(),
                                         m_base[2]->GetBdata (),
                                         tmp3,tmp6,wsp,
                                         true,false,true);

            Vmath::Vadd(m_ncoeffs, tmp6, 1, outarray, 1, outarray, 1);

            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata (),
                                         m_base[1]->GetBdata (),
                                         m_base[2]->GetDbdata(),
                                         tmp4,tmp6,wsp,
                                         true,true,false);

            Vmath::Vadd(m_ncoeffs, tmp6, 1, outarray, 1, outarray, 1);
        }


        //---------------------------------------
        // Evaluation functions
        //---------------------------------------

        StdRegions::StdExpansionSharedPtr PyrExp::v_GetStdExp(void) const
        {
            return MemoryManager<StdRegions::StdPyrExp>
                ::AllocateSharedPtr(m_base[0]->GetBasisKey(),
                                    m_base[1]->GetBasisKey(),
                                    m_base[2]->GetBasisKey());
        }

        StdRegions::StdExpansionSharedPtr PyrExp::v_GetLinStdExp(void) const
        {
            LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(),
                           2, m_base[0]->GetPointsKey());
            LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(),
                           2, m_base[1]->GetPointsKey());
            LibUtilities::BasisKey bkey2(m_base[2]->GetBasisType(),
                           2, m_base[2]->GetPointsKey());

            return MemoryManager<StdRegions::StdPyrExp>
                ::AllocateSharedPtr( bkey0, bkey1, bkey2);
        }

        /*
         * @brief Get the coordinates #coords at the local coordinates
         * #Lcoords
         */
        void PyrExp::v_GetCoord(const Array<OneD, const NekDouble>& Lcoords,
                                      Array<OneD,       NekDouble>& coords)
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

        void PyrExp::v_GetCoords(
            Array<OneD, NekDouble> &coords_1,
            Array<OneD, NekDouble> &coords_2,
            Array<OneD, NekDouble> &coords_3)
        {
            Expansion::v_GetCoords(coords_1, coords_2, coords_3);
        }


        void PyrExp::v_ExtractDataToCoeffs(
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

            // Check if not same order or basis and if not make temp
            // element to read in data
            if (fromType[0] != m_base[0]->GetBasisType() ||
                fromType[1] != m_base[1]->GetBasisType() ||
                fromType[2] != m_base[2]->GetBasisType() ||
                data_order0 != fillorder0 ||
                data_order1 != fillorder1 ||
                data_order2 != fillorder2)
            {
                // Construct a pyr with the appropriate basis type at our
                // quadrature points, and one more to do a forwards
                // transform. We can then copy the output to coeffs.
                StdRegions::StdPyrExp tmpPyr(
                    LibUtilities::BasisKey(
                        fromType[0], data_order0, m_base[0]->GetPointsKey()),
                    LibUtilities::BasisKey(
                        fromType[1], data_order1, m_base[1]->GetPointsKey()),
                    LibUtilities::BasisKey(
                        fromType[2], data_order2, m_base[2]->GetPointsKey()));

                StdRegions::StdPyrExp tmpPyr2(m_base[0]->GetBasisKey(),
                                              m_base[1]->GetBasisKey(),
                                              m_base[2]->GetBasisKey());

                Array<OneD, const NekDouble> tmpData(tmpPyr.GetNcoeffs(), data);
                Array<OneD, NekDouble> tmpBwd(tmpPyr2.GetTotPoints());
                Array<OneD, NekDouble> tmpOut(tmpPyr2.GetNcoeffs());

                tmpPyr.BwdTrans(tmpData, tmpBwd);
                tmpPyr2.FwdTrans(tmpBwd, tmpOut);
                Vmath::Vcopy(tmpOut.size(), &tmpOut[0], 1, coeffs, 1);

            }
            else
            {
                Vmath::Vcopy(m_ncoeffs,&data[0],1,coeffs,1);
            }
        }

        /**
         * Given the local cartesian coordinate \a Lcoord evaluate the
         * value of physvals at this point by calling through to the
         * StdExpansion method
         */
        NekDouble PyrExp::v_StdPhysEvaluate(
                const Array<OneD, const NekDouble> &Lcoord,
                const Array<OneD, const NekDouble> &physvals)
        {
            // Evaluate point in local coordinates.
            return StdPyrExp::v_PhysEvaluate(Lcoord,physvals);
        }

        NekDouble PyrExp::v_PhysEvaluate(const Array<OneD, const NekDouble>& coord,
                                         const Array<OneD, const NekDouble>& physvals)
        {
            Array<OneD,NekDouble> Lcoord(3);

            ASSERTL0(m_geom,"m_geom not defined");

            //TODO: check GetLocCoords()
            m_geom->GetLocCoords(coord, Lcoord);

            return StdPyrExp::v_PhysEvaluate(Lcoord, physvals);
        }


        //---------------------------------------
        // Helper functions
        //---------------------------------------

        int PyrExp::v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }

        void PyrExp::v_GetFacePhysMap(const int               face,
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
                if(outarray.size()!=nq0*nq1)
                {
                    outarray = Array<OneD, int>(nq0*nq1);
                }

                //Directions A and B positive
                for(int i = 0; i < nquad0*nquad1; ++i)
                {
                    outarray[i] = i;
                }

                break;
              case 1:
                    nq0 = nquad0;
                    nq1 = nquad2;
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    //Direction A and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        for(int i = 0; i < nquad0; ++i)
                        {
                            outarray[k*nquad0+i] = (nquad0*nquad1*k)+i;
                        }
                    }

                    break;
                case 2:
                    nq0 = nquad1;
                    nq1 = nquad2;
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    //Directions A and B positive
                    for(int j = 0; j < nquad1*nquad2; ++j)
                    {
                        outarray[j] = nquad0-1 + j*nquad0;

                    }
                    break;
                case 3:

                    nq0 = nquad0;
                    nq1 = nquad2;
                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    //Direction A and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        for(int i = 0; i < nquad0; ++i)
                        {
                            outarray[k*nquad0+i] = nquad0*(nquad1-1) + (nquad0*nquad1*k)+i;
                        }
                    }
                    break;
                case 4:
                    nq0 = nquad1;
                    nq1 = nquad2;

                    if(outarray.size()!=nq0*nq1)
                    {
                        outarray = Array<OneD, int>(nq0*nq1);
                    }

                    //Directions A and B positive
                    for(int j = 0; j < nquad1*nquad2; ++j)
                    {
                        outarray[j] = j*nquad0;

                    }
                    break;
                default:
                    ASSERTL0(false,"face value (> 4) is out of range");
                    break;
            }
        }

        void PyrExp::v_ComputeFaceNormal(const int face)
        {
            const SpatialDomains::GeomFactorsSharedPtr &geomFactors =
                GetGeom()->GetMetricInfo();

            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
            for(int i = 0; i < ptsKeys.size(); ++i)
            {
                // Need at least 2 points for computing normals
                if (ptsKeys[i].GetNumPoints() == 1)
                {
                    LibUtilities::PointsKey pKey(2, ptsKeys[i].GetPointsType());
                    ptsKeys[i] = pKey;
                }
            }

            SpatialDomains::GeomType type            = geomFactors->GetGtype();
            const Array<TwoD, const NekDouble> &df   = geomFactors->GetDerivFactors(ptsKeys);
            const Array<OneD, const NekDouble> &jac  = geomFactors->GetJac(ptsKeys);

            LibUtilities::BasisKey tobasis0 = DetFaceBasisKey(face,0);
            LibUtilities::BasisKey tobasis1 = DetFaceBasisKey(face,1);

            // Number of quadrature points in face expansion.
            int nq_face = tobasis0.GetNumPoints()*tobasis1.GetNumPoints();

            int vCoordDim = GetCoordim();
            int i;

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
                            normal[i][0] = -df[3*i+2][0];
                        }
                        break;
                    }
                    case 1:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            normal[i][0] = -df[3*i+1][0];
                        }
                        break;
                    }
                    case 2:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            normal[i][0] = df[3*i][0]+df[3*i+2][0];
                        }
                        break;
                    }
                    case 3:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            normal[i][0] = df[3*i+1][0]+df[3*i+2][0];
                        }
                        break;
                    }
                    case 4:
                    {
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            normal[i][0] = -df[3*i][0];
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

                Vmath::Fill(nqb, fac, length, 1);

                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Fill(nq_face,fac*normal[i][0],normal[i],1);
                }
            }
            else
            {
                // Set up deformed normals.
                int j, k;

                int nq0  = ptsKeys[0].GetNumPoints();
                int nq1  = ptsKeys[1].GetNumPoints();
                int nq2  = ptsKeys[2].GetNumPoints();
                int nq01 = nq0*nq1;
                int nqtot;

                // Determine number of quadrature points on the face.
                if (face == 0)
                {
                    nqtot = nq0*nq1;
                }
                else if (face == 1 || face == 3)
                {
                    nqtot = nq0*nq2;
                }
                else
                {
                    nqtot = nq1*nq2;
                }

                LibUtilities::PointsKey points0;
                LibUtilities::PointsKey points1;

                Array<OneD, NekDouble> faceJac(nqtot);
                Array<OneD, NekDouble> normals(vCoordDim*nqtot,0.0);

                // Extract Jacobian along face and recover local derivatives
                // (dx/dr) for polynomial interpolation by multiplying m_gmat by
                // jacobian
                switch(face)
                {
                    case 0:
                    {
                        for(j = 0; j < nq01; ++j)
                        {
                            normals[j]         = -df[2][j]*jac[j];
                            normals[nqtot+j]   = -df[5][j]*jac[j];
                            normals[2*nqtot+j] = -df[8][j]*jac[j];
                            faceJac[j]         = jac[j];
                        }

                        points0 = ptsKeys[0];
                        points1 = ptsKeys[1];
                        break;
                    }

                    case 1:
                    {
                        for (j = 0; j < nq0; ++j)
                        {
                            for(k = 0; k < nq2; ++k)
                            {
                                int tmp = j+nq01*k;
                                normals[j+k*nq0]          =
                                    -df[1][tmp]*jac[tmp];
                                normals[nqtot+j+k*nq0]    =
                                    -df[4][tmp]*jac[tmp];
                                normals[2*nqtot+j+k*nq0]  =
                                    -df[7][tmp]*jac[tmp];
                                faceJac[j+k*nq0] = jac[tmp];
                            }
                        }

                        points0 = ptsKeys[0];
                        points1 = ptsKeys[2];
                        break;
                    }

                    case 2:
                    {
                        for (j = 0; j < nq1; ++j)
                        {
                            for(k = 0; k < nq2; ++k)
                            {
                                int tmp = nq0-1+nq0*j+nq01*k;
                                normals[j+k*nq1]         =
                                    (df[0][tmp]+df[2][tmp])*jac[tmp];
                                normals[nqtot+j+k*nq1]   =
                                    (df[3][tmp]+df[5][tmp])*jac[tmp];
                                normals[2*nqtot+j+k*nq1] =
                                    (df[6][tmp]+df[8][tmp])*jac[tmp];
                                faceJac[j+k*nq1] = jac[tmp];
                            }
                        }

                        points0 = ptsKeys[1];
                        points1 = ptsKeys[2];
                        break;
                    }

                    case 3:
                    {
                        for (j = 0; j < nq0; ++j)
                        {
                            for(k = 0; k < nq2; ++k)
                            {
                                int tmp = nq0*(nq1-1) + j + nq01*k;
                                normals[j+k*nq0]         =
                                    (df[1][tmp]+df[2][tmp])*jac[tmp];
                                normals[nqtot+j+k*nq0]   =
                                    (df[4][tmp]+df[5][tmp])*jac[tmp];
                                normals[2*nqtot+j+k*nq0] =
                                    (df[7][tmp]+df[8][tmp])*jac[tmp];
                                faceJac[j+k*nq0] = jac[tmp];
                            }
                        }

                        points0 = ptsKeys[0];
                        points1 = ptsKeys[2];
                        break;
                    }

                    case 4:
                    {
                        for (j = 0; j < nq1; ++j)
                        {
                            for(k = 0; k < nq2; ++k)
                            {
                                int tmp = j*nq0+nq01*k;
                                normals[j+k*nq1]         =
                                    -df[0][tmp]*jac[tmp];
                                normals[nqtot+j+k*nq1]   =
                                    -df[3][tmp]*jac[tmp];
                                normals[2*nqtot+j+k*nq1] =
                                    -df[6][tmp]*jac[tmp];
                                faceJac[j+k*nq1] = jac[tmp];
                            }
                        }

                        points0 = ptsKeys[1];
                        points1 = ptsKeys[2];
                        break;
                    }

                    default:
                        ASSERTL0(false,"face is out of range (face < 4)");
                }

                Array<OneD, NekDouble> work   (nq_face, 0.0);
                // Interpolate Jacobian and invert
                LibUtilities::Interp2D(points0, points1, faceJac,
                                       tobasis0.GetPointsKey(),
                                       tobasis1.GetPointsKey(),
                                       work);
                Vmath::Sdiv(nq_face, 1.0, &work[0], 1, &work[0], 1);

                // Interpolate normal and multiply by inverse Jacobian.
                for(i = 0; i < vCoordDim; ++i)
                {
                    LibUtilities::Interp2D(points0, points1,
                                           &normals[i*nqtot],
                                           tobasis0.GetPointsKey(),
                                           tobasis1.GetPointsKey(),
                                           &normal[i][0]);
                    Vmath::Vmul(nq_face,work,1,normal[i],1,normal[i],1);
                }

                // Normalise to obtain unit normals.
                Vmath::Zero(nq_face,work,1);
                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vvtvp(nq_face,normal[i],1,normal[i],1,work,1,work,1);
                }

                Vmath::Vsqrt(nq_face,work,1,work,1);
                Vmath::Sdiv (nq_face,1.0,work,1,work,1);
                
                Vmath::Vcopy(nqb, work, 1, length, 1);

                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vmul(nq_face,normal[i],1,work,1,normal[i],1);
                }
            }
        }

        void PyrExp::v_SVVLaplacianFilter(
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
            StdPyrExp::v_SVVLaplacianFilter( array, mkey);

            // Divide by sqrt(Jac)
            Vmath::Vdiv(nq,array,1,sqrt_jac,1,array,1);
        }


        //---------------------------------------
        // Matrix creation functions
        //---------------------------------------

        DNekMatSharedPtr PyrExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
                returnval = StdPyrExp::v_GenMatrix(mkey);
            }

            return returnval;
        }

        DNekMatSharedPtr PyrExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
            StdRegions::StdPyrExpSharedPtr tmp =
                MemoryManager<StdPyrExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);

            return tmp->GetStdMatrix(mkey);
        }

        DNekScalMatSharedPtr PyrExp::v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }

        DNekScalBlkMatSharedPtr PyrExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
        }

        void PyrExp::v_DropLocStaticCondMatrix(const MatrixKey &mkey)
        {
            m_staticCondMatrixManager.DeleteObject(mkey);
        }

        DNekScalMatSharedPtr PyrExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

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
                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
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
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetShapeType(),
                                                         *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble fac = 1.0/(m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed ||
                        mkey.GetNVarCoeff() > 0 ||
                        mkey.ConstFactorExists(
                                StdRegions::eFactorSVVCutoffRatio))
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
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
                        Array<TwoD, const NekDouble> gmat =
                                            m_metricinfo->GetGmat(ptsKeys);

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
                                                ::AllocateSharedPtr(jac, lap);
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    MatrixKey masskey(StdRegions::eMass, mkey.GetShapeType(), *this);
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian, mkey.GetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols);

                    (*helm) = LapMat + factor*MassMat;

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, helm);
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
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }

        DNekScalBlkMatSharedPtr PyrExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            unsigned int nbdry = NumBndryCoeffs();
            unsigned int nint = (unsigned int)(m_ncoeffs - nbdry);
            unsigned int exp_size[] = {nbdry, nint};
            unsigned int nblks = 2;
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

        void PyrExp::v_ComputeLaplacianMetric()
        {
            if (m_metrics.count(eMetricQuadrature) == 0)
            {
                ComputeQuadratureMetric();
            }

            int i, j;
            const unsigned int nqtot = GetTotPoints();
            const unsigned int dim = 3;
            const MetricType m[3][3] = {
                { eMetricLaplacian00, eMetricLaplacian01, eMetricLaplacian02 },
                { eMetricLaplacian01, eMetricLaplacian11, eMetricLaplacian12 },
                { eMetricLaplacian02, eMetricLaplacian12, eMetricLaplacian22 }
            };

            for (unsigned int i = 0; i < dim; ++i)
            {
                for (unsigned int j = i; j < dim; ++j)
                {
                    m_metrics[m[i][j]] = Array<OneD, NekDouble>(nqtot);
                }
            }

            // Define shorthand synonyms for m_metrics storage
            Array<OneD,NekDouble> g0   (m_metrics[m[0][0]]);
            Array<OneD,NekDouble> g1   (m_metrics[m[1][1]]);
            Array<OneD,NekDouble> g2   (m_metrics[m[2][2]]);
            Array<OneD,NekDouble> g3   (m_metrics[m[0][1]]);
            Array<OneD,NekDouble> g4   (m_metrics[m[0][2]]);
            Array<OneD,NekDouble> g5   (m_metrics[m[1][2]]);

            // Allocate temporary storage
            Array<OneD,NekDouble> alloc(9*nqtot,0.0);
            Array<OneD,NekDouble> h0   (nqtot, alloc);
            Array<OneD,NekDouble> h1   (nqtot, alloc+ 1*nqtot);
            Array<OneD,NekDouble> h2   (nqtot, alloc+ 2*nqtot);
            Array<OneD,NekDouble> wsp1 (nqtot, alloc+ 3*nqtot);
            Array<OneD,NekDouble> wsp2 (nqtot, alloc+ 4*nqtot);
            Array<OneD,NekDouble> wsp3 (nqtot, alloc+ 5*nqtot);
            Array<OneD,NekDouble> wsp4 (nqtot, alloc+ 6*nqtot);
            Array<OneD,NekDouble> wsp5 (nqtot, alloc+ 7*nqtot);
            Array<OneD,NekDouble> wsp6 (nqtot, alloc+ 8*nqtot);

            const Array<TwoD, const NekDouble>& df =
                                m_metricinfo->GetDerivFactors(GetPointsKeys());
            const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();
            const Array<OneD, const NekDouble>& z2 = m_base[2]->GetZ();
            const unsigned int nquad0 = m_base[0]->GetNumPoints();
            const unsigned int nquad1 = m_base[1]->GetNumPoints();
            const unsigned int nquad2 = m_base[2]->GetNumPoints();

            // Populate collapsed coordinate arrays h0, h1 and h2.
            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Vmath::Fill(nquad0, 2.0/(1.0-z2[j]),         &h0[0]+i*nquad0 + j*nquad0*nquad1,1);
                    Vmath::Fill(nquad0, 1.0/(1.0-z2[j]),         &h1[0]+i*nquad0 + j*nquad0*nquad1,1);
                    Vmath::Fill(nquad0, (1.0+z1[i])/(1.0-z2[j]), &h2[0]+i*nquad0 + j*nquad0*nquad1,1);
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
                // f_{1k}
                Vmath::Vvtvvtp(nqtot, &df[0][0], 1, &h0[0], 1, &df[2][0], 1, &h1[0], 1, &wsp1[0], 1);
                Vmath::Vvtvvtp(nqtot, &df[3][0], 1, &h0[0], 1, &df[5][0], 1, &h1[0], 1, &wsp2[0], 1);
                Vmath::Vvtvvtp(nqtot, &df[6][0], 1, &h0[0], 1, &df[8][0], 1, &h1[0], 1, &wsp3[0], 1);

                // g0
                Vmath::Vvtvvtp(nqtot, &wsp1[0], 1, &wsp1[0], 1, &wsp2[0], 1, &wsp2[0], 1, &g0[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp3[0], 1, &wsp3[0], 1, &g0[0],   1, &g0[0],   1);

                // g4
                Vmath::Vvtvvtp(nqtot, &df[2][0], 1, &wsp1[0], 1, &df[5][0], 1, &wsp2[0], 1, &g4[0], 1);
                Vmath::Vvtvp  (nqtot, &df[8][0], 1, &wsp3[0], 1, &g4[0], 1, &g4[0], 1);

                // f_{2k}
                Vmath::Vvtvvtp(nqtot, &df[1][0], 1, &h0[0], 1, &df[2][0], 1, &h2[0], 1, &wsp4[0], 1);
                Vmath::Vvtvvtp(nqtot, &df[4][0], 1, &h0[0], 1, &df[5][0], 1, &h2[0], 1, &wsp5[0], 1);
                Vmath::Vvtvvtp(nqtot, &df[7][0], 1, &h0[0], 1, &df[8][0], 1, &h2[0], 1, &wsp6[0], 1);

                // g1
                Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp4[0], 1, &wsp5[0], 1, &wsp5[0], 1, &g1[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp6[0], 1, &wsp6[0], 1, &g1[0],   1, &g1[0],   1);

                // g3
                Vmath::Vvtvvtp(nqtot, &wsp1[0], 1, &wsp4[0], 1, &wsp2[0], 1, &wsp5[0], 1, &g3[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp3[0], 1, &wsp6[0], 1, &g3[0],   1, &g3[0],   1);

                // g5
                Vmath::Vvtvvtp(nqtot, &df[2][0], 1, &wsp4[0], 1, &df[5][0], 1, &wsp5[0], 1, &g5[0], 1);
                Vmath::Vvtvp  (nqtot, &df[8][0], 1, &wsp6[0], 1, &g5[0], 1, &g5[0], 1);

                // g2
                Vmath::Vvtvvtp(nqtot, &df[2][0], 1, &df[2][0], 1, &df[5][0], 1, &df[5][0], 1, &g2[0], 1);
                Vmath::Vvtvp  (nqtot, &df[8][0], 1, &df[8][0], 1, &g2[0], 1, &g2[0], 1);
            }
            else
            {
                // f_{1k}
                Vmath::Svtsvtp(nqtot, df[0][0], &h0[0], 1, df[2][0], &h1[0], 1, &wsp1[0], 1);
                Vmath::Svtsvtp(nqtot, df[3][0], &h0[0], 1, df[5][0], &h1[0], 1, &wsp2[0], 1);
                Vmath::Svtsvtp(nqtot, df[6][0], &h0[0], 1, df[8][0], &h1[0], 1, &wsp3[0], 1);

                // g0
                Vmath::Vvtvvtp(nqtot, &wsp1[0], 1, &wsp1[0], 1, &wsp2[0], 1, &wsp2[0], 1, &g0[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp3[0], 1, &wsp3[0], 1, &g0[0],   1, &g0[0],   1);

                // g4
                Vmath::Svtsvtp(nqtot, df[2][0], &wsp1[0], 1, df[5][0], &wsp2[0], 1, &g4[0], 1);
                Vmath::Svtvp  (nqtot, df[8][0], &wsp3[0], 1, &g4[0], 1, &g4[0], 1);

                // f_{2k}
                Vmath::Svtsvtp(nqtot, df[1][0], &h0[0], 1, df[2][0], &h2[0], 1, &wsp4[0], 1);
                Vmath::Svtsvtp(nqtot, df[4][0], &h0[0], 1, df[5][0], &h2[0], 1, &wsp5[0], 1);
                Vmath::Svtsvtp(nqtot, df[7][0], &h0[0], 1, df[8][0], &h2[0], 1, &wsp6[0], 1);

                // g1
                Vmath::Vvtvvtp(nqtot, &wsp4[0], 1, &wsp4[0], 1, &wsp5[0], 1, &wsp5[0], 1, &g1[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp6[0], 1, &wsp6[0], 1, &g1[0],   1, &g1[0],   1);

                // g3
                Vmath::Vvtvvtp(nqtot, &wsp1[0], 1, &wsp4[0], 1, &wsp2[0], 1, &wsp5[0], 1, &g3[0], 1);
                Vmath::Vvtvp  (nqtot, &wsp3[0], 1, &wsp6[0], 1, &g3[0],   1, &g3[0],   1);

                // g5
                Vmath::Svtsvtp(nqtot, df[2][0], &wsp4[0], 1, df[5][0], &wsp5[0], 1, &g5[0], 1);
                Vmath::Svtvp  (nqtot, df[8][0], &wsp6[0], 1, &g5[0], 1, &g5[0], 1);

                // g2
                Vmath::Fill(nqtot, df[2][0]*df[2][0] + df[5][0]*df[5][0] + df[8][0]*df[8][0], &g2[0], 1);
            }

            for (unsigned int i = 0; i < dim; ++i)
            {
                for (unsigned int j = i; j < dim; ++j)
                {
                    MultiplyByQuadratureMetric(m_metrics[m[i][j]],
                                               m_metrics[m[i][j]]);

                }
            }
        }

        void PyrExp::v_LaplacianMatrixOp_MatFree_Kernel(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
                  Array<OneD,       NekDouble> &wsp)
        {
            // This implementation is only valid when there are no coefficients
            // associated to the Laplacian operator
            if (m_metrics.count(eMetricLaplacian00) == 0)
            {
                ComputeLaplacianMetric();
            }

            int nquad0  = m_base[0]->GetNumPoints();
            int nquad1  = m_base[1]->GetNumPoints();
            int nq2  = m_base[2]->GetNumPoints();
            int nqtot   = nquad0*nquad1*nq2;

            ASSERTL1(wsp.size() >= 6*nqtot,
                     "Insufficient workspace size.");
            ASSERTL1(m_ncoeffs <= nqtot,
                     "Workspace not set up for ncoeffs > nqtot");

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
            Array<OneD,NekDouble> wsp0 (2*nqtot, wsp);
            Array<OneD,NekDouble> wsp1 (  nqtot, wsp+1*nqtot);
            Array<OneD,NekDouble> wsp2 (  nqtot, wsp+2*nqtot);
            Array<OneD,NekDouble> wsp3 (  nqtot, wsp+3*nqtot);
            Array<OneD,NekDouble> wsp4 (  nqtot, wsp+4*nqtot);
            Array<OneD,NekDouble> wsp5 (  nqtot, wsp+5*nqtot);

            // LAPLACIAN MATRIX OPERATION
            // wsp1 = du_dxi1 = D_xi1 * inarray = D_xi1 * u
            // wsp2 = du_dxi2 = D_xi2 * inarray = D_xi2 * u
            // wsp2 = du_dxi3 = D_xi3 * inarray = D_xi3 * u
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
    }//end of namespace
}//end of namespace
