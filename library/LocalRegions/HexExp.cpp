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
// Description: Methods for Hex expansion in local regoins
//
///////////////////////////////////////////////////////////////////////////////


#include <LocalRegions/HexExp.h>
#include <LibUtilities/Foundations/Interp.h>
#include <SpatialDomains/HexGeom.h>

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
            StdRegions::StdHexExp(Ba,Bb,Bc),
            Expansion     (),
            Expansion3D   (),
            m_geom(geom),
            m_metricinfo(m_geom->GetGeomFactors(m_base)),
            m_matrixManager(
                    boost::bind(&HexExp::CreateMatrix, this, _1),
                    std::string("HexExpMatrix")),
            m_staticCondMatrixManager(
                    boost::bind(&HexExp::CreateStaticCondMatrix, this, _1),
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
            StdRegions::StdHexExp(T),
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
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
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

            Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> Diff0 = Array<OneD,NekDouble>(ntot);
            Array<OneD,NekDouble> Diff1 = Array<OneD,NekDouble>(ntot);
            Array<OneD,NekDouble> Diff2 = Array<OneD,NekDouble>(ntot);

            StdHexExp::v_PhysDeriv(inarray, Diff0, Diff1, Diff2);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul (ntot,&gmat[0][0],1,&Diff0[0],1, &out_d0[0], 1);
                    Vmath::Vvtvp(ntot,&gmat[1][0],1,&Diff1[0],1, &out_d0[0], 1,
                                                                 &out_d0[0],1);
                    Vmath::Vvtvp(ntot,&gmat[2][0],1,&Diff2[0],1, &out_d0[0], 1,
                                                                 &out_d0[0],1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Vmul (ntot,&gmat[3][0],1,&Diff0[0],1, &out_d1[0], 1);
                    Vmath::Vvtvp(ntot,&gmat[4][0],1,&Diff1[0],1, &out_d1[0], 1,
                                                                 &out_d1[0],1);
                    Vmath::Vvtvp(ntot,&gmat[5][0],1,&Diff2[0],1, &out_d1[0], 1,
                                                                 &out_d1[0],1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Vmul (ntot,&gmat[6][0],1,&Diff0[0],1, &out_d2[0], 1);
                    Vmath::Vvtvp(ntot,&gmat[7][0],1,&Diff1[0],1, &out_d2[0], 1,
                                                                 &out_d2[0],1);
                    Vmath::Vvtvp(ntot,&gmat[8][0],1,&Diff2[0],1, &out_d2[0], 1,
                                                                 &out_d2[0],1);
                }
            }
            else // regular geometry
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul (ntot,gmat[0][0],&Diff0[0],1, &out_d0[0], 1);
                    Blas::Daxpy (ntot,gmat[1][0],&Diff1[0],1, &out_d0[0], 1);
                    Blas::Daxpy (ntot,gmat[2][0],&Diff2[0],1, &out_d0[0], 1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Smul (ntot,gmat[3][0],&Diff0[0],1, &out_d1[0], 1);
                    Blas::Daxpy (ntot,gmat[4][0],&Diff1[0],1, &out_d1[0], 1);
                    Blas::Daxpy (ntot,gmat[5][0],&Diff2[0],1, &out_d1[0], 1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Smul (ntot,gmat[6][0],&Diff0[0],1, &out_d2[0], 1);
                    Blas::Daxpy (ntot,gmat[7][0],&Diff1[0],1, &out_d2[0], 1);
                    Blas::Daxpy (ntot,gmat[8][0],&Diff2[0],1, &out_d2[0], 1);
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
                      Array<OneD,       NekDouble> &outarray)
        {
            const int nqtot = GetTotPoints();

            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,       NekDouble> tmp(nqtot);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nqtot, &jac[0], 1, (NekDouble*)&inarray[0], 1,
                                   &tmp[0], 1);
            }
            else
            {
                Vmath::Smul(nqtot, jac[0], (NekDouble*)&inarray[0], 1,
                            &tmp[0], 1);
            }

            StdHexExp::v_IProductWRTBase(tmp, outarray);
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

            int    nquad0  = m_base[0]->GetNumPoints();
            int    nquad1  = m_base[1]->GetNumPoints();
            int    nquad2  = m_base[2]->GetNumPoints();
            int    nqtot   = nquad0*nquad1*nquad2;
            int    nmodes0 = m_base[0]->GetNumModes();
            int    nmodes1 = m_base[1]->GetNumModes();
 
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD, NekDouble> alloc(3*nqtot + 2*m_ncoeffs +
                                         nmodes0*nquad2*(nquad1+nmodes1));
            Array<OneD, NekDouble> tmp1(alloc);               // Dir1 metric
            Array<OneD, NekDouble> tmp2(alloc +   nqtot);     // Dir2 metric
            Array<OneD, NekDouble> tmp3(alloc + 2*nqtot);     // Dir3 metric
            Array<OneD, NekDouble> tmp4(alloc + 3*nqtot);     // Dir1 iprod
            Array<OneD, NekDouble> tmp5(tmp4  +   m_ncoeffs); // Dir2 iprod
            Array<OneD, NekDouble> wsp (tmp4  + 2*m_ncoeffs); // Wsp

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nqtot,&gmat[3*dir][0],  1,inarray.get(),1,tmp1.get(),1);
                Vmath::Vmul(nqtot,&gmat[3*dir+1][0],1,inarray.get(),1,tmp2.get(),1);
                Vmath::Vmul(nqtot,&gmat[3*dir+2][0],1,inarray.get(),1,tmp3.get(),1);
            }
            else
            {
                Vmath::Smul(nqtot, gmat[3*dir][0],  inarray.get(),1,tmp1.get(), 1);
                Vmath::Smul(nqtot, gmat[3*dir+1][0],inarray.get(),1,tmp2.get(), 1);
                Vmath::Smul(nqtot, gmat[3*dir+2][0],inarray.get(),1,tmp3.get(), 1);
            }  

            MultiplyByQuadratureMetric(tmp1,tmp1);
            MultiplyByQuadratureMetric(tmp2,tmp2);
            MultiplyByQuadratureMetric(tmp3,tmp3);
            
            IProductWRTBase_SumFacKernel(   m_base[0]->GetDbdata(),
                                            m_base[1]->GetBdata(),
                                            m_base[2]->GetBdata(),
                                            tmp1,tmp4,wsp,
                                            false,true,true);
            IProductWRTBase_SumFacKernel(   m_base[0]->GetBdata(),
                                            m_base[1]->GetDbdata(),
                                            m_base[2]->GetBdata(),
                                            tmp2,tmp5,wsp,
                                            true,false,true);
            IProductWRTBase_SumFacKernel(   m_base[0]->GetBdata(),
                                            m_base[1]->GetBdata(),
                                            m_base[2]->GetDbdata(),
                                            tmp3,outarray,wsp,
                                            true,true,false);
                                            
            Vmath::Vadd(GetNcoeffs(), tmp4, 1, outarray, 1, outarray, 1);
            Vmath::Vadd(GetNcoeffs(), tmp5, 1, outarray, 1, outarray, 1);
        }


        void HexExp::IProductWRTDerivBase_MatOp(
                const int dir, 
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
            DNekScalMatSharedPtr iprodmat = m_matrixManager[iprodmatkey];
            
            Blas::Dgemv('N',m_ncoeffs,nq,iprodmat->Scale(),(iprodmat->GetOwnedMatrix())->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }


        //-----------------------------
        // Evaluation functions
        //-----------------------------
        /**
         * \brief Interpolate the solution at given coordinates
	 *
	 * Evaluate the expansion at an arbitrary physical coordinate.
         * @param   coord       An array with three elements containing the
         *                      x,y,z coordinates of the point at which to
         *                      evaluate the expansion.
         * @returns The value of the expansion at the given point.
         */
        NekDouble HexExp::v_PhysEvaluate(
                const Array<OneD, const NekDouble> &coord)
        {
            return PhysEvaluate(coord,m_phys);
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


        /**
	 * \brief Retrieve the local coordinates of each quadrature point.
         * The coordinates are put into the three arrays provided.
         * @param   coords_0    Coordinate component in first direction.
         * @param   coords_1    Coordinate component in second direction.
         * @param   coords_2    Coordinate component in third direction.
         */
        void HexExp::v_GetCoords(
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
            m_geom->FillGeom();
             
            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements(),
                         "output coords_2 is not defined");
                CBasis0 = m_geom->GetBasis(2,0);
                CBasis1 = m_geom->GetBasis(2,1);
                CBasis2 = m_geom->GetBasis(2,2);

                if( m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey())
                        && m_base[1]->GetBasisKey().SamePoints(
                                                        CBasis1->GetBasisKey())
                        && m_base[2]->GetBasisKey().SamePoints(
                                                        CBasis2->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(2);
                    Blas::Dcopy(GetTotPoints(), x, 1, coords_2, 1);
                }
                else // Interpolate to Expansion point distribution
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
                ASSERTL0(coords_1.num_elements(),
                         "output coords_1 is not defined");

                CBasis0 = m_geom->GetBasis(1,0);
                CBasis1 = m_geom->GetBasis(1,1);
                CBasis2 = m_geom->GetBasis(1,2);

                if( m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey())
                        && m_base[1]->GetBasisKey().SamePoints(
                                                        CBasis1->GetBasisKey())
                        && m_base[2]->GetBasisKey().SamePoints(
                                                        CBasis2->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(1);
                    Blas::Dcopy(GetTotPoints(),
                                x, 1, coords_1, 1);
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
                ASSERTL0(coords_0.num_elements(),
                         "output coords_0 is not defined");

                CBasis0 = m_geom->GetBasis(0,0);
                CBasis1 = m_geom->GetBasis(0,1);
                CBasis2 = m_geom->GetBasis(0,2);

                if( m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey())
                        && m_base[1]->GetBasisKey().SamePoints(
                                                        CBasis1->GetBasisKey())
                        && m_base[2]->GetBasisKey().SamePoints(
                                                        CBasis2->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(0);
                    Blas::Dcopy(GetTotPoints(),
                                x, 1, coords_0, 1);
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


        //-----------------------------
        // Helper functions
        //-----------------------------
        /**
         * Writes the expansion evaluated at the quadrature points to a text
         * file suitable for reading by a variety of plotting programs.
         * @param   outfile     Output stream to write data to.
         * @param   format      Chosen format of output data.
         * @param   dumpVar     If true, write out the variable names too.
         * @param   var         If dumpVar set, uses this variable name.
         */
        void HexExp::v_WriteToFile(
                std::ofstream &outfile,
                OutputFormat format,
                const bool dumpVar,
                std::string var)
        {
            if(format==eTecplot)
            {
                int i, j;
                int nquad0 = m_base[0]->GetNumPoints();
                int nquad1 = m_base[1]->GetNumPoints();
                int nquad2 = m_base[2]->GetNumPoints();
                Array<OneD,NekDouble> coords[3];

                ASSERTL0(m_geom,"m_geom not defined");

                coords[0] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
                coords[1] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
                coords[2] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);

                GetCoords(coords[0],coords[1],coords[2]);

                if(dumpVar)
                {
                    outfile << "Variables = x,  y,  z";
                    outfile << ", "<< var << std::endl << std::endl;
                }

                outfile << "Zone, I=" << nquad0 << ", J=" << nquad1
                        << ", K=" << nquad2 << ", F=Point" << std::endl;

                for(i = 0; i < nquad0*nquad1*nquad2; ++i)
                {
                    for(j = 0; j < 3; ++j)
                    {
                            outfile << coords[j][i] << " ";
                    }
                    outfile << m_phys[i] << std::endl;
                }
            }
            else if(format==eGnuplot)
            {
                int nquad0 = m_base[0]->GetNumPoints();
                int nquad1 = m_base[1]->GetNumPoints();
                int nquad2 = m_base[2]->GetNumPoints();
                Array<OneD,NekDouble> coords[3];

                ASSERTL0(m_geom,"m_geom not defined");

                coords[0] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
                coords[1] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
                coords[2] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);

                GetCoords(coords[0],coords[1],coords[2]);

                for(int k = 0; k < nquad2; ++k)
                {
                    for(int j = 0; j < nquad1; ++j)
                    {
                        for(int i = 0; i < nquad0; ++i)
                        {
                            int n = (k*nquad1 + j)*nquad0 + i;
                            outfile << coords[0][n] << " "
                                    << coords[1][n] << " "
                                    << coords[2][n] << " "
                                    << m_phys[i + nquad0*(j + nquad1*k)]
                                    << endl;
                        }
                        outfile << endl;
                    }
                    outfile << endl;
                }
            }

            else
            {
                ASSERTL0(false, "Output routine not implemented for requested "
                                "type of output");
            }
        }


        /// Return the region shape using the enum-list of ShapeType
        StdRegions::ExpansionType HexExp::v_DetExpansionType() const
        {
            return StdRegions::eHexahedron;
        }

        
        const SpatialDomains::GeomFactorsSharedPtr& HexExp::v_GetMetricInfo() const
        {
            return m_metricinfo;
        }

        
        /// Returns the HexGeom object associated with this expansion.
        const SpatialDomains::GeometrySharedPtr HexExp::v_GetGeom() const
        {
            return m_geom;
        }

        
        /// Returns the HexGeom object associated with this expansion.
        const SpatialDomains::Geometry3DSharedPtr& HexExp::v_GetGeom3D() const
        {
            return m_geom;
        }

        
        int HexExp::v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }

        
        StdRegions::Orientation HexExp::v_GetFaceOrient(int face)
        {
            return m_geom->GetFaceOrient(face);
        }

        
        bool HexExp::v_GetFaceDGForwards(const int i) const
        {
            StdRegions::Orientation fo = m_geom->GetFaceOrient(i);
            
            return fo == StdRegions::eDir1FwdDir1_Dir2FwdDir2 || 
                   fo == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                   fo == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                   fo == StdRegions::eDir1FwdDir2_Dir2BwdDir1;
        }

        ///Returns the physical values at the quadrature points of a face
        void HexExp::v_GetFacePhysVals(
            const int                                face,
            const StdRegions::StdExpansionSharedPtr &FaceExp,
            const Array<OneD, const NekDouble>      &inarray,
                  Array<OneD,       NekDouble>      &outarray,
            StdRegions::Orientation                  orient)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            
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
                break;
            case 1:
                if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Direction A and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1*k),
				     1,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
                {
                    //Direction A negative and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0-1)+(nquad0*nquad1*k),
                                     -1,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2)
                {
                    //Direction A positive and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1*(nquad2-1-k)),
                                     1,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2)
                {
                    //Direction A negative and B negative
                    for(int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0-1)+(nquad0*nquad1*(nquad2-1-k)),
                                     -1,&(outarray[0])+(k*nquad0),1);
                    }
                }
		else if(orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A and B positive
		    for (int i=0; i<nquad0; i++)
                    {
                        Vmath::Vcopy(nquad2,&(inarray[0])+i,nquad0*nquad1,
                                     &(outarray[0])+(i*nquad2),1);
                    }
		}
		else if(orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A positive and B negative
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0-1-i),nquad0*nquad1,
                                     &(outarray[0])+(i*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A negative and B positive
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+nquad0*nquad1*(nquad2-1)+i,
                                     -nquad0*nquad1,&(outarray[0])+(i*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A and B negative
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+nquad0*nquad1*nquad2+(nquad0-1-i),
                                     -nquad0*nquad1,&(outarray[0])+(i*nquad2),1);
                    }
		} 
                break;
            case 2:
	        if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad0*nquad1,&(inarray[0])+(nquad0-1),
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
                else if(orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
                {
                    //Direction A negative and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1-1)+(k*nquad0*nquad1),
                                     -1,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2)
                {
                    //Direction A positive and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*(nquad1-1))+(nquad0*nquad1*(nquad2-1-k)),
                                     1,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2)
                {
                    //Direction A negative and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1-1)+(nquad0*nquad1*(nquad2-1-k)),
                                     -1,&(outarray[0])+(k*nquad0),1);
                    }
                }
		else if(orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A and B positive
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+nquad0*(nquad1-1)+i,nquad0*nquad1,
                                     &(outarray[0])+(i*nquad2),1);
                    }
		}
		else if(orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A positive and B negative
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0*nquad1-1-i),nquad0*nquad1,
                                     &(outarray[0])+(i*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A negative and B positive
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+nquad0*(nquad1*nquad2-1)+i,-nquad0*nquad1,
                                     &(outarray[0])+(i*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A and B negative
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0*nquad1*nquad2-1-i),-nquad0*nquad1,
                                     &(outarray[0])+(i*nquad2),1);
                    }
		} 
                break;
            case 4:
                if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad0*nquad1,&(inarray[0]),nquad0,&(outarray[0]),1);
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
                {
                    //Direction A negative and B positive
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+nquad0*(nquad1-1)+(k*nquad0*nquad1),
                                     -nquad0,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2)
                {
                    //Direction A positive and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1*(nquad2-1-k)),
                                     nquad0,&(outarray[0])+(k*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2)
                {
                    //Direction A negative and B negative
                    for (int k=0; k<nquad2; k++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+nquad0*(nquad1-1)+(nquad0*nquad1*(nquad2-1-k)),
                                     -nquad0,&(outarray[0])+(k*nquad0),1);
                    }
                }
		else if(orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A and B positive
		    for (int j=0; j<nquad0; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+j*nquad0,nquad0*nquad1,
                                     &(outarray[0])+(j*nquad2),1);
                    }
		}
		else if(orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A positive and B negative
		    for (int j=0; j<nquad0; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0*(nquad1-1)-j*nquad0),
                                     nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A negative and B positive
		    for (int j=0; j<nquad0; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+nquad0*nquad1*(nquad2-1)+j*nquad0,
                                     -nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A and B negative
		    for (int j=0; j<nquad0; j++)
                    {
		        Vmath::Vcopy(nquad2,&(inarray[0])+(nquad0*(nquad1*nquad2-1)-j*nquad0),
                                     -nquad0*nquad1,&(outarray[0])+(j*nquad2),1);
                    }
		} 
                break;
            case 5:
                if(orient == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    //Directions A and B positive
                    Vmath::Vcopy(nquad0*nquad1,&(inarray[0])+nquad0*nquad1*(nquad2-1),1,&(outarray[0]),1);
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
                {
                    //Direction A negative and B positive
                    for (int j=0; j<nquad1; j++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+nquad0*nquad1*(nquad2-1)+(nquad0-1+j*nquad0),
                                     -1,&(outarray[0])+(j*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2)
                {
                    //Direction A positive and B negative
                    for (int j=0; j<nquad1; j++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+((nquad0*nquad1*nquad2-1)-(nquad0-1)-j*nquad0),
                                     1,&(outarray[0])+(j*nquad0),1);
                    }
                }
                else if(orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2)
                {
                    //Direction A negative and B negative
                    for (int j=0; j<nquad1; j++)
                    {
                        Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1*nquad2-1-j*nquad0),
                                     -1,&(outarray[0])+(j*nquad0),1);
                    }
                }
		else if(orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A and B positive
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad1,&(inarray[0])+nquad0*nquad1*(nquad2-1)+i,nquad0,
                                     &(outarray[0])+(i*nquad1),1);
                    }
		}
		else if(orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A positive and B negative
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad1,&(inarray[0])+nquad0*nquad1*(nquad2-1)+(nquad0-1-i),
                                     nquad0,&(outarray[0])+(i*nquad1),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1)
		{
		    //Transposed, Direction A negative and B positive
		    for (int i=0; i<nquad0; i++)
                    {
		        Vmath::Vcopy(nquad1,&(inarray[0])+nquad0*(nquad1*nquad2-1)+i,-nquad0,
                                     &(outarray[0])+(i*nquad1),1);
                    }
		} 
		else if(orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
		{
		    //Transposed, Direction A and B negative
		    for (int i=0; i<nquad0; i++)
                    {
		      Vmath::Vcopy(nquad1,&(inarray[0])+(nquad0*nquad1*nquad2-1-i),-nquad0,
                                   &(outarray[0])+(i*nquad1),1);
                    }
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
            const SpatialDomains::GeomFactorsSharedPtr & geomFactors = GetGeom()->GetMetricInfo();
            SpatialDomains::GeomType type = geomFactors->GetGtype();
            const Array<TwoD, const NekDouble> & gmat = geomFactors->GetGmat();
            const Array<OneD, const NekDouble> & jac  = geomFactors->GetJac();
            int nqe = m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints();
            int vCoordDim = GetCoordim();

            m_faceNormals[face] = Array<OneD, Array<OneD, NekDouble> >(vCoordDim);
            Array<OneD, Array<OneD, NekDouble> > &normal = m_faceNormals[face];
            for (i = 0; i < vCoordDim; ++i)
            {
                normal[i] = Array<OneD, NekDouble>(nqe);
            }
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
                        Vmath::Fill(nqe,-gmat[3*i+2][0],normal[i],1);
                    }
                    break;
                case 1:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,-gmat[3*i+1][0],normal[i],1);
                    }
                    break;
                case 2:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,gmat[3*i][0],normal[i],1);
                    }
                    break;
                case 3:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,gmat[3*i+1][0],normal[i],1);
                    }
                    break;
                case 4:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,-gmat[3*i][0],normal[i],1);
                    }
                    break;
                case 5:
                    for(i = 0; i < vCoordDim; ++i)
                    {
                        Vmath::Fill(nqe,gmat[3*i+2][0],normal[i],1);
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
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Smul(nqe,fac,normal[i],1,normal[i],1);
                }
		
	    }
            else   // Set up deformed normals
                {
	        int j, k;

                int nquad0 = geomFactors->GetPointsKey(0).GetNumPoints();
                int nquad1 = geomFactors->GetPointsKey(1).GetNumPoints();
		int nquad2 = geomFactors->GetPointsKey(2).GetNumPoints();
		int nqtot = nquad0*nquad1;

                LibUtilities::PointsKey points0;
		LibUtilities::PointsKey points1;

                int nq;
                Array<OneD,NekDouble> work(nqe,0.0);
                Array<OneD,NekDouble> normals(vCoordDim*nquad0*nquad1,0.0);

                // Extract Jacobian along face and recover local
                // derivates (dx/dr) for polynomial interpolation by
                // multiplying m_gmat by jacobian
                switch(face)
	        {
                case 0:
		    for(j = 0; j < nquad0*nquad1; ++j)
                    {
                        normals[j] = -gmat[2][j]*jac[j];
                        normals[nqtot+j] = -gmat[5][j]*jac[j];
                        normals[2*nqtot+j] = -gmat[8][j]*jac[j];
                    }
                    points0 = geomFactors->GetPointsKey(0);
                    points1 = geomFactors->GetPointsKey(1);
		    nq=points0.GetNumPoints()*points1.GetNumPoints();

                    // interpolate Jacobian and invert
		    LibUtilities::Interp2D(points0,points1,jac,m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),work);
                    Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);
                    
                    // interpolate
                    for(i = 0; i < GetCoordim(); ++i)
                    {
		        LibUtilities::Interp2D(points0,points1,&normals[i*nq],m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),&normal[i][0]); 
                        Vmath::Vmul(nqe,work,1,normal[i],1,normal[i],1);
                    }
                    break;
		case 1:
		    for (j=0; j< nquad0; ++j)
                    {
                        for(k=0; k<nquad2; ++k)
                        {
                            normals[j+k*nquad0]  = -gmat[1][j+nquad0*nquad1*k]*jac[j+nquad0*nquad1*k];
                            normals[nqtot+j+k*nquad0]  = -gmat[4][j+nquad0*nquad1*k]*jac[j+nquad0*nquad1*k];
                            normals[2*nqtot+j+k*nquad0]  = -gmat[7][j+nquad0*nquad1*k]*jac[j+nquad0*nquad1*k];
                        } 
                    }
                    points0 = geomFactors->GetPointsKey(0);
                    points1 = geomFactors->GetPointsKey(2);
		    nq=points0.GetNumPoints()*points1.GetNumPoints();

                    // interpolate Jacobian and invert
		    LibUtilities::Interp2D(points0,points1,jac,m_base[0]->GetPointsKey(),m_base[2]->GetPointsKey(),work);
                    Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);

                    // interpolate
                    for(i = 0; i < GetCoordim(); ++i)
                    {
		        LibUtilities::Interp2D(points0,points1,&normals[i*nq],m_base[0]->GetPointsKey(),m_base[2]->GetPointsKey(),&normal[i][0]); 
                        Vmath::Vmul(nqe,work,1,normal[i],1,normal[i],1);
                    }
                    break;
                case 2:
		    for (j=0; j< nquad1; ++j)
                    {
                        for(k=0; k<nquad2; ++k)
                        {
			  normals[j+k*nquad0]  = gmat[0][nquad0-1+nquad0*j+nquad0*nquad1*k]*jac[nquad0-1+nquad0*j+nquad0*nquad1*k];
			  normals[nqtot+j+k*nquad0]  = gmat[3][nquad0-1+nquad0*j+nquad0*nquad1*k]*jac[nquad0-1+nquad0*j+nquad0*nquad1*k];
			  normals[2*nqtot+j+k*nquad0]  = gmat[6][nquad0-1+nquad0*j+nquad0*nquad1*k]*jac[nquad0-1+nquad0*j+nquad0*nquad1*k];
                        } 
                    }
                    points0 = geomFactors->GetPointsKey(1);
                    points1 = geomFactors->GetPointsKey(2);
		    nq=points0.GetNumPoints()*points1.GetNumPoints();

                    // interpolate Jacobian and invert
		    LibUtilities::Interp2D(points0,points1,jac,m_base[1]->GetPointsKey(),m_base[2]->GetPointsKey(),work);
                    Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);

                    // interpolate
                    for(i = 0; i < GetCoordim(); ++i)
                    {
		        LibUtilities::Interp2D(points0,points1,&normals[i*nq],m_base[1]->GetPointsKey(),m_base[2]->GetPointsKey(),&normal[i][0]); 
                        Vmath::Vmul(nqe,work,1,normal[i],1,normal[i],1);
                    }
                    break;               
                case 3:
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
		    nq=points0.GetNumPoints()*points1.GetNumPoints();

                    // interpolate Jacobian and invert
		    LibUtilities::Interp2D(points0,points1,jac,m_base[0]->GetPointsKey(),m_base[2]->GetPointsKey(),work);
                    Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);

                    // interpolate
                    for(i = 0; i < GetCoordim(); ++i)
                    {
		        LibUtilities::Interp2D(points0,points1,&normals[i*nq],m_base[0]->GetPointsKey(),m_base[2]->GetPointsKey(),&normal[i][0]); 
                        Vmath::Vmul(nqe,work,1,normal[i],1,normal[i],1);
                    }
                    break;
                case 4:
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
		    nq=points0.GetNumPoints()*points1.GetNumPoints();

                    // interpolate Jacobian and invert
		    LibUtilities::Interp2D(points0,points1,jac,m_base[1]->GetPointsKey(),m_base[2]->GetPointsKey(),work);
                    Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);

                    // interpolate
                    for(i = 0; i < GetCoordim(); ++i)
                    {
		        LibUtilities::Interp2D(points0,points1,&normals[i*nq],m_base[1]->GetPointsKey(),m_base[2]->GetPointsKey(),&normal[i][0]); 
                        Vmath::Vmul(nqe,work,1,normal[i],1,normal[i],1);
                    }
                    break;
                case 5:
                    for (j=0; j< nquad0*nquad1; ++j)
                    {
                        normals[j]  = gmat[2][j+nquad0*nquad1*(nquad2-1)]*jac[j+nquad0*nquad1*(nquad2-1)];
                        normals[nqtot+j]  = gmat[5][j+nquad0*nquad1*(nquad2-1)]*jac[j+nquad0*nquad1*(nquad2-1)];
                        normals[2*nqtot+j]  = gmat[8][j+nquad0*nquad1*(nquad2-1)]*jac[j+nquad0*nquad1*(nquad2-1)];
                    }
                    points0 = geomFactors->GetPointsKey(0);
                    points1 = geomFactors->GetPointsKey(1);
                    nq=points0.GetNumPoints()*points1.GetNumPoints();
                    
                    // interpolate Jacobian and invert
		    LibUtilities::Interp2D(points0,points1,jac,m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),work);
                    Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);

                    // interpolate
                    for(i = 0; i < GetCoordim(); ++i)
                    {
		        LibUtilities::Interp2D(points0,points1,&normals[i*nq],m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),&normal[i][0]); 
                        Vmath::Vmul(nqe,work,1,normal[i],1,normal[i],1);
                    }
                    break;                    
		default:
                    ASSERTL0(false,"face is out of range (face < 5)");
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
            }
        }

        NekDouble HexExp::v_Linf()
        {
            return Linf();
        }

        NekDouble HexExp::v_L2(const Array<OneD, const NekDouble> &sol)
        {
            return StdExpansion::L2(sol);
        }

        NekDouble HexExp::v_L2()
        {
            return StdExpansion::L2();
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
         * @param   inarray     Input coefficients.
         * @param   output      Output coefficients.
         * @param   mkey        Matrix key 
         */
        void HexExp::v_LaplacianMatrixOp_MatFree(
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
/*                  int       nquad0  = m_base[0]->GetNumPoints();
                    int       nquad1  = m_base[1]->GetNumPoints();
                    int       nquad2  = m_base[2]->GetNumPoints();
                    int       nqtot   = nquad0*nquad1*nquad2; 
                    int       nmodes0 = m_base[0]->GetNumModes();
                    int       nmodes1 = m_base[1]->GetNumModes();
                    int       nmodes2 = m_base[2]->GetNumModes();
                    int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);
                    
                    const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                    const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
                    const Array<OneD, const NekDouble>& base2  = m_base[2]->GetBdata();
                    const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
                    const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
                    const Array<OneD, const NekDouble>& dbase2 = m_base[2]->GetDbdata();
                    const Array<TwoD, const NekDouble>& metric = m_metricinfo->GetLaplacianMetrics(); 
                    
                    // Allocate temporary storage
                    Array<OneD,NekDouble> wsp0(4*wspsize);
                    Array<OneD,NekDouble> wsp1(wsp0+wspsize);
                    Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);
                    Array<OneD,NekDouble> wsp3(wsp0+3*wspsize);
                    
                    if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                    {  
                        // LAPLACIAN MATRIX OPERATION
                        // wsp0 = u       = B   * u_hat 
                        // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                        // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                        StdHexExp::BwdTrans_SumFacKernel(base0,base1,base2,inarray,wsp0,wsp1,true,true,true);
                        StdExpansion3D::PhysTensorDeriv(wsp0,wsp1,wsp2,wsp3);
                    }
                    else
                    {
                        StdExpansion3D::PhysTensorDeriv(inarray,wsp1,wsp2,wsp3);                    
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
                    IProductWRTBase_SumFacKernel(dbase0,base1,base2,wsp0,outarray,wsp1,false,true,true);
                    IProductWRTBase_SumFacKernel(base0,dbase1,base2,wsp2,wsp1,    wsp0,true,false,true);
//                    IProductWRTBase_SumFacKernel(base0,base1,dbase2,wsp2,wsp1,    wsp0,true,true,false);
                    
                    // outarray = outarray + wsp1
                    //          = L * u_hat
                    Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);      
*/                }
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
                                      nquad0*nquad1*(nquad2+nmodes0)+
                                      nmodes0*nmodes1*nquad2);
                    
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


        void HexExp::v_HelmholtzMatrixOp_MatFree(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey)
        {
            if(m_metricinfo->IsUsingLaplMetrics())
            {
                ASSERTL0(false,"Finish implementing HexExp Helmholtz for Lapl Metrics");
/*              int       nquad0  = m_base[0]->GetNumPoints();
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
                                  nquad0*nquad1*(nquad2+nmodes0)+
                                  nmodes0*nmodes1*nquad2);
                
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
                        NekDouble jac = (m_metricinfo->GetJac())[0];
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
                                                    DetExpansionType(), *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble fac = 1.0/(m_metricinfo->GetJac())[0];
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

                        (*WeakDeriv) = gmat[3*dir  ][0]*deriv0
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
            case StdRegions::ePreconditioner:
            case StdRegions::eHelmholtz:
                {
                    NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);
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

//                    StdRegions::StdMatrixKey hkey(StdRegions::eHybridDGHelmholtz,
//                                                  DetExpansionType(),*this,
//                                                  mkey.GetConstant(0),
//                                                  mkey.GetConstant(1));
                    MatrixKey hkey(StdRegions::eHybridDGHelmholtz, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
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


        DNekScalBlkMatSharedPtr HexExp::CreateStaticCondMatrix(const MatrixKey &mkey)
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
        

        void HexExp::MultiplyByQuadratureMetric(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray)
        {        
            const int nqtot = m_base[0]->GetNumPoints() *
                              m_base[1]->GetNumPoints() *
                              m_base[2]->GetNumPoints();
               
            if(m_metricinfo->IsUsingQuadMetrics())
            {
                const Array<OneD, const NekDouble> &metric
                    = m_metricinfo->GetQuadratureMetrics();

                Vmath::Vmul(nqtot, metric, 1, inarray, 1, outarray, 1);
            }
            else
            {
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

                StdHexExp::MultiplyByQuadratureMetric(outarray, outarray);
            }
        }
        
        void HexExp::LaplacianMatrixOp_MatFree_Kernel(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                      Array<OneD,       NekDouble> &wsp)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            int nqtot  = nquad0*nquad1*nquad2; 
            
            const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble>& base2  = m_base[2]->GetBdata();
            const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase2 = m_base[2]->GetDbdata();

            // Allocate temporary storage
            Array<OneD,NekDouble> alloc(9*nqtot,0.0);
            Array<OneD,NekDouble> wsp1 (alloc        );// TensorDeriv 1
            Array<OneD,NekDouble> wsp2 (alloc+1*nqtot);// TensorDeriv 2
            Array<OneD,NekDouble> wsp3 (alloc+2*nqtot);// TensorDeriv 3
            Array<OneD,NekDouble> wsp4 (alloc+3*nqtot);// g0
            Array<OneD,NekDouble> wsp5 (alloc+4*nqtot);// g1
            Array<OneD,NekDouble> wsp6 (alloc+5*nqtot);// g2
            Array<OneD,NekDouble> wsp7 (alloc+6*nqtot);// g3
            Array<OneD,NekDouble> wsp8 (alloc+7*nqtot);// g4
            Array<OneD,NekDouble> wsp9 (alloc+8*nqtot);// g5
            
            StdExpansion3D::PhysTensorDeriv(inarray,wsp1,wsp2,wsp3);
            
            // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
            // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
            // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
            // especially for this purpose
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();
            
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // Compute geometric factor composites
                // wsp3 = g0*g0 + g3*g3 + g6*g6 
                Vmath::Vvtvvtp(nqtot,&gmat[0][0],1,&gmat[0][0],1,&gmat[3][0],1,&gmat[3][0],1,&wsp4[0],1);
                Vmath::Vvtvp  (nqtot,&gmat[6][0],1,&gmat[6][0],1,&wsp4[0],1,&wsp4[0],1);
                // wsp4 = g1*g1 + g4*g4 + g7*g7;
                Vmath::Vvtvvtp(nqtot,&gmat[1][0],1,&gmat[1][0],1,&gmat[4][0],1,&gmat[4][0],1,&wsp5[0],1);
                Vmath::Vvtvp  (nqtot,&gmat[7][0],1,&gmat[7][0],1,&wsp5[0],1,&wsp5[0],1);
                // wsp5 = g2*g2 + g5*g5 + g8*g8; 
                Vmath::Vvtvvtp(nqtot,&gmat[2][0],1,&gmat[2][0],1,&gmat[5][0],1,&gmat[5][0],1,&wsp6[0],1);
                Vmath::Vvtvp  (nqtot,&gmat[8][0],1,&gmat[8][0],1,&wsp6[0],1,&wsp6[0],1);
                // wsp6 = g0*g1 + g3*g4 + g6*g7
                Vmath::Vvtvvtp(nqtot,&gmat[0][0],1,&gmat[1][0],1,&gmat[3][0],1,&gmat[4][0],1,&wsp7[0],1);
                Vmath::Vvtvp  (nqtot,&gmat[6][0],1,&gmat[7][0],1,&wsp7[0],1,&wsp7[0],1);
                // wsp7 = g0*g2 + g3*g5 + g6*g8
                Vmath::Vvtvvtp(nqtot,&gmat[0][0],1,&gmat[2][0],1,&gmat[3][0],1,&gmat[5][0],1,&wsp8[0],1);
                Vmath::Vvtvp  (nqtot,&gmat[6][0],1,&gmat[8][0],1,&wsp8[0],1,&wsp8[0],1);
                // wsp8 = g1*g2 + g4*g5 + g7*g8
                Vmath::Vvtvvtp(nqtot,&gmat[1][0],1,&gmat[2][0],1,&gmat[4][0],1,&gmat[5][0],1,&wsp9[0],1);
                Vmath::Vvtvp  (nqtot,&gmat[7][0],1,&gmat[8][0],1,&wsp9[0],1,&wsp9[0],1);
                
                // Multiply wsp1,2,3 by the appropriate factor composites
                Vmath::Vvtvvtp(nqtot,&wsp4[0],1,&wsp1[0],1,&wsp7[0],1,&wsp2[0],1,&wsp4[0],1);
                Vmath::Vvtvp  (nqtot,&wsp8[0],1,&wsp3[0],1,&wsp4[0],1,&wsp4[0],1);
                Vmath::Vvtvvtp(nqtot,&wsp5[0],1,&wsp2[0],1,&wsp7[0],1,&wsp1[0],1,&wsp5[0],1);
                Vmath::Vvtvp  (nqtot,&wsp9[0],1,&wsp3[0],1,&wsp5[0],1,&wsp5[0],1);
                Vmath::Vvtvvtp(nqtot,&wsp6[0],1,&wsp3[0],1,&wsp8[0],1,&wsp1[0],1,&wsp6[0],1);
                Vmath::Vvtvp  (nqtot,&wsp9[0],1,&wsp2[0],1,&wsp6[0],1,&wsp6[0],1);
            }
            else
            {
                NekDouble g0 = gmat[0][0]*gmat[0][0] + gmat[3][0]*gmat[3][0] + gmat[6][0]*gmat[6][0];
                NekDouble g1 = gmat[1][0]*gmat[1][0] + gmat[4][0]*gmat[4][0] + gmat[7][0]*gmat[7][0];
                NekDouble g2 = gmat[2][0]*gmat[2][0] + gmat[5][0]*gmat[5][0] + gmat[8][0]*gmat[8][0];
                NekDouble g3 = gmat[0][0]*gmat[1][0] + gmat[3][0]*gmat[4][0] + gmat[6][0]*gmat[7][0];
                NekDouble g4 = gmat[0][0]*gmat[2][0] + gmat[3][0]*gmat[5][0] + gmat[6][0]*gmat[8][0];
                NekDouble g5 = gmat[1][0]*gmat[2][0] + gmat[4][0]*gmat[5][0] + gmat[7][0]*gmat[8][0];
                
                Vmath::Svtsvtp(nqtot,g0,&wsp1[0],1,g3,&wsp2[0],1,&wsp4[0],1);
                Vmath::Svtvp  (nqtot,g4,&wsp3[0],1,&wsp4[0],1,&wsp4[0],1);
                Vmath::Svtsvtp(nqtot,g1,&wsp2[0],1,g3,&wsp1[0],1,&wsp5[0],1);
                Vmath::Svtvp  (nqtot,g5,&wsp3[0],1,&wsp5[0],1,&wsp5[0],1);
                Vmath::Svtsvtp(nqtot,g2,&wsp3[0],1,g4,&wsp1[0],1,&wsp6[0],1);
                Vmath::Svtvp  (nqtot,g5,&wsp2[0],1,&wsp6[0],1,&wsp6[0],1);
            }
            
            MultiplyByQuadratureMetric(wsp4,wsp4);
            MultiplyByQuadratureMetric(wsp5,wsp5);
            MultiplyByQuadratureMetric(wsp6,wsp6);
            
            // outarray = m = (D_xi1 * B)^T * k 
            // wsp1     = n = (D_xi2 * B)^T * l 
            IProductWRTBase_SumFacKernel(dbase0,base1,base2,wsp4,outarray,wsp,false,true,true);
            IProductWRTBase_SumFacKernel(base0,dbase1,base2,wsp5,wsp1,    wsp,true,false,true);
            IProductWRTBase_SumFacKernel(base0,base1,dbase2,wsp6,wsp2,    wsp,true,true,false);
            
            // outarray = outarray + wsp1
            //          = L * u_hat
            Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);
            Vmath::Vadd(m_ncoeffs,wsp2.get(),1,outarray.get(),1,outarray.get(),1);     
        }
    }//end of namespace
}//end of namespace
