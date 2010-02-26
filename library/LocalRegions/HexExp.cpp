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

#include <LocalRegions/LocalRegions.h>
#include <LocalRegions/LocalRegions.hpp>
#include <stdio.h>
#include <LocalRegions/HexExp.h>

namespace Nektar
{
    namespace LocalRegions
    {
        /**
         * @class HexExp
         * Defines a hexahedral local expansion.
         */

        /**
         * @param   Ba          Basis key for first coordinate.
         * @param   Bb          Basis key for second coordinate.
         * @param   Bc          Basis key for third coordinate.
         */
        HexExp::HexExp( const LibUtilities::BasisKey &Ba,
                            const LibUtilities::BasisKey &Bb,
                            const LibUtilities::BasisKey &Bc,
                            const SpatialDomains::HexGeomSharedPtr &geom):
            StdRegions::StdHexExp(Ba,Bb,Bc),
            m_geom(geom),
            m_metricinfo(m_geom->GetGeomFactors(m_base)),
            m_matrixManager(std::string("HexExpMatrix")),
            m_staticCondMatrixManager(std::string("HexExpStaticCondMatrix"))
        {
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(
                                MatrixKey((StdRegions::MatrixType) i,
                                StdRegions::eNoExpansionType,*this),
                                boost::bind(&HexExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(
                                MatrixKey((StdRegions::MatrixType) i,
                                StdRegions::eNoExpansionType,*this),
                                boost::bind(&HexExp::CreateStaticCondMatrix,
                                this, _1));
            }
        }


        /**
         * @param   T           HexExp to copy.
         */
        HexExp::HexExp( const HexExp &T ):
            StdRegions::StdHexExp(T),
            m_geom(T.m_geom),
            m_metricinfo(T.m_metricinfo),
            m_matrixManager(std::string("HexExpMatrix")),
            m_staticCondMatrixManager(std::string("HexExpStaticCondMatrix"))
        {
        }


        HexExp::~HexExp()
        {
        }


        /**
         * @param   inarray     Input array of physical space data.
         * @param   outarray    Output array of data.
         */
        void HexExp::v_IProductWRTBase(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> & outarray)
        {
            HexExp::v_IProductWRTBase(m_base[0]->GetBdata(),
                            m_base[1]->GetBdata(),
                            m_base[2]->GetBdata(),
                            inarray,outarray,1);
        }


        /**
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
        void HexExp::v_IProductWRTBase(
                            const Array<OneD, const NekDouble>& base0,
                            const Array<OneD, const NekDouble>& base1,
                            const Array<OneD, const NekDouble>& base2,
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> & outarray,
                            int coll_check)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,NekDouble> tmp(nquad0*nquad1*nquad2);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,
                            (NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,jac[0],
                            (NekDouble*)&inarray[0],1,&tmp[0],1);
            }

            StdHexExp::v_IProductWRTBase(base0, base1, base2, tmp, outarray, 1);
        }


        /**
         * Forward transform from physical quadrature space stored in \a
         * inarray and evaluate the expansion coefficients and store in \a
         * (this)->m_coeffs.
         * @param   inarray     Input array
         * @param   outarray    Output array
         */
        void HexExp::v_FwdTrans( const Array<OneD, const NekDouble> & inarray,
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


        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////

        /**
         * For Hexahedral region can use the Tensor_Deriv function defined
         * under StdExpansion.
         * @param   inarray     Input array
         * @param   out_d0      Derivative of \a inarray in first direction.
         * @param   out_d1      Derivative of \a inarray in second direction.
         * @param   out_d2      Derivative of \a inarray in third direction.
         */
        void HexExp::v_PhysDeriv(const Array<OneD, const NekDouble> & inarray,
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

            StdHexExp::PhysDeriv(inarray, Diff0, Diff1, Diff2);

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
                if(out_d0.num_elements())                {
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
         * @param   dir         Direction in which to compute derivative.
         *                      Valid values are 0, 1, 2.
         * @param   inarray     Input array.
         * @param   outarray    Output array.
         */
        void HexExp::v_PhysDeriv(const int dir,
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


        /**
         * Evaluate the expansion at an arbitrary physical coordinate.
         * @param   coord       An array with three elements containing the
         *                      x,y,z coordinates of the point at which to
         *                      evaluate the expansion.
         * @returns The value of the expansion at the given point.
         */
        NekDouble HexExp::v_PhysEvaluate(
                            const Array<OneD, const NekDouble> &coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(3);

            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);
            return StdHexExp::v_PhysEvaluate(Lcoord);
        }


        /**
         * The coordinates are put into the three arrays provided.
         * @param   coords_0    Coordinate component in first direction.
         * @param   coords_1    Coordinate component in second direction.
         * @param   coords_2    Coordinate component in third direction.
         */
        void HexExp::v_GetCoords( Array<OneD,NekDouble> &coords_0,
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
         * @param   Lcoords     Local coordinates in reference space.
         * @param   coords      Corresponding coordinates in physical space.
         */
        void HexExp::v_GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                            Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[0] <= 1.0 &&
                     Lcoords[1] >= -1.0 && Lcoords[1] <= 1.0 &&
                     Lcoords[2] >= -1.0 && Lcoords[2] <= 1.0,
                     "Local coordinates are not in region [-1,1]");

              m_geom->FillGeom();

            for(i = 0; i < m_geom->GetCoordDim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }


        /**
         * Writes the expansion evaluated at the quadrature points to a text
         * file suitable for reading by a variety of plotting programs.
         * @param   outfile     Output stream to write data to.
         * @param   format      Chosen format of output data.
         * @param   dumpVar     If true, write out the variable names too.
         * @param   var         If dumpVar set, uses this variable name.
         */
        void HexExp::v_WriteToFile(std::ofstream &outfile,
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


        /**
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
            returnVal = StdHexExp::Integral(tmp);

            return  returnVal;
        }



        DNekMatSharedPtr HexExp::CreateStdMatrix(
                            const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();

            StdRegions::StdHexExpSharedPtr tmp = MemoryManager<StdHexExp>
                                        ::AllocateSharedPtr(bkey0,bkey1,bkey2);

            return tmp->GetStdMatrix(mkey);
        }

//         DNekMatSharedPtr HexExp::CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
//         {
//             // Need to check if matrix exists in stdMatrixManager.
//             // If not then make a local expansion with standard metric info
//             // and generate matrix. Otherwise direct call is OK.
//             if(!StdMatManagerAlreadyCreated(mkey))
//             {
//                 LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
//                 LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
//                 LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
//                 HexExpSharedPtr tmp = MemoryManager<HexExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);
//
//                 return tmp->StdHexExp::GetStdMatrix(mkey); //TODO: infinity recursive recursion -- this is not working
//             }
//             else
//             {
//                 return StdHexExp::GetStdMatrix(mkey);
//             }
//         }
/*
        DNekBlkMatSharedPtr HexExp::CreateStdStaticCondMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            // Need to check if matrix exists in stdMatrixManager.
            // If not then make a local expansion with standard metric info
            // and generate matrix. Otherwise direct call is OK.
            if(!StdMatManagerAlreadyCreated(mkey))
            {
                LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
                LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
                LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
                HexExpSharedPtr tmp = MemoryManager<HexExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);

                return tmp->StdHexExp::GetStdStaticCondMatrix(mkey);
            }
            else
            {
                return StdHexExp::GetStdStaticCondMatrix(mkey);
            }
        }*/


        DNekMatSharedPtr HexExp::GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
                returnval = Expansion3D::GenMatrix(mkey);
                break;
            default:
                returnval = StdHexExp::GenMatrix(mkey);
            }

            return returnval;
        }

        DNekScalMatSharedPtr HexExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;


            ASSERTL2(m_metricinfo->GetGtype() == SpatialDomains::eNoGeomType,
                     "Geometric information is not set up");

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat
                                        = GetStdMatrix(*mkey.GetStdMatKey());
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
                                        = GetStdMatrix(*mkey.GetStdMatKey());
                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(fac,mat);
                    }
                }
                break;
            case StdRegions::eWeakDeriv0:
            case StdRegions::eWeakDeriv1:
            case StdRegions::eWeakDeriv2:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());

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
                        }

                        MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                            mkey.GetExpansionType(), *this);
                        MatrixKey deriv1key(StdRegions::eWeakDeriv1,
                                            mkey.GetExpansionType(), *this);

                        DNekMat &deriv0
                                    = *GetStdMatrix(*deriv0key.GetStdMatKey());
                        DNekMat &deriv1
                                    = *GetStdMatrix(*deriv1key.GetStdMatKey());

                        int rows = deriv0.GetRows();
                        int cols = deriv1.GetColumns();

                        DNekMatSharedPtr WeakDeriv = MemoryManager<DNekMat>
                                                ::AllocateSharedPtr(rows,cols);
                        (*WeakDeriv) = gmat[2*dir][0]*deriv0
                                                + gmat[2*dir+1][0]*deriv1;

                        returnval = MemoryManager<DNekScalMat>
                                            ::AllocateSharedPtr(jac,WeakDeriv);
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        ASSERTL1(m_geom->GetCoordDim() == 2,
                                 "Standard Region Laplacian is only set up for "
                                 "Quads in two-dimensional");
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetExpansionType(), *this);

                        DNekMat &lap00
                                    = *GetStdMatrix(*lap00key.GetStdMatKey());
                        DNekMat &lap01
                                    = *GetStdMatrix(*lap01key.GetStdMatKey());
                        DNekMat &lap11
                                    = *GetStdMatrix(*lap11key.GetStdMatKey());

                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat
                                                    = m_metricinfo->GetGmat();

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>
                                                ::AllocateSharedPtr(rows,cols);

                        (*lap) = (gmat[0][0]*gmat[0][0]
                                + gmat[2][0]*gmat[2][0])*lap00
                                + (gmat[0][0]*gmat[1][0]
                                    + gmat[2][0]*gmat[3][0])
                                * (lap01 + Transpose(lap01))
                                + (gmat[1][0]*gmat[1][0]
                                    + gmat[3][0]*gmat[3][0])*lap11;

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(jac,lap);
                    }
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
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());
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

                    DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            case StdRegions::eInvHybridDGHelmholtz:
                {
                    NekDouble one = 1.0;

                    StdRegions::StdMatrixKey hkey(StdRegions::eHybridDGHelmholtz,
                                                  DetExpansionType(),*this,
                                                  mkey.GetConstant(0),
                                                  mkey.GetConstant(1));
                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            default:
                {
                    NekDouble        one = 1.0;
                    DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            }

            return returnval;
        }


        DNekScalBlkMatSharedPtr HexExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() == SpatialDomains::eNoGeomType,"Geometric information is not set up");

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
                    DNekScalMatSharedPtr& mat = GetLocMatrix(mkey);
                    factor = mat->Scale();
                    goto UseStdRegionsMatrix;
                }
                break;
            UseStdRegionsMatrix:
                {
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekBlkMatSharedPtr& mat = GetStdStaticCondMatrix(*(mkey.GetStdMatKey()));
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
/*
        void HexExp::v_IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble> &outarray)
        {
            cout << "Hex::IProduct::SumFac" << endl;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
            int    order1 = m_base[1]->GetNumModes();

            Array<OneD,NekDouble> tmp(nquad0*nquad1*nquad2);
            Array<OneD,NekDouble> wsp(nquad0*nquad1*(nquad2+order0)
                                                + order0*order1*nquad2);

            MultiplyByQuadratureMetric(inarray, tmp);

            StdHexExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                    m_base[1]->GetBdata(),
                                                    m_base[2]->GetBdata(),
                                                    tmp, outarray, wsp,
                                                    true, true, true);
        }
*/
    }//end of namespace
}//end of namespace

/**
 *    $Log: HexExp.cpp,v $
 *    Revision 1.27  2009/12/15 18:09:02  cantwell
 *    Split GeomFactors into 1D, 2D and 3D
 *    Added generation of tangential basis into GeomFactors
 *    Updated ADR2DManifold solver to use GeomFactors for tangents
 *    Added <GEOMINFO> XML session section support in MeshGraph
 *    Fixed const-correctness in VmathArray
 *    Cleaned up LocalRegions code to generate GeomFactors
 *    Removed GenSegExp
 *    Temporary fix to SubStructuredGraph
 *    Documentation for GlobalLinSys and GlobalMatrix classes
 *
 *    Revision 1.26  2009/10/30 14:00:06  pvos
 *    Multi-level static condensation updates
 *
 *    Revision 1.25  2009/05/01 13:23:21  pvos
 *    Fixed various bugs
 *
 *    Revision 1.24  2009/04/27 21:34:07  sherwin
 *    Updated WriteToField
 *
 *    Revision 1.23  2009/01/21 16:59:56  pvos
 *    Added additional geometric factors to improve efficiency
 *
 *    Revision 1.22  2008/09/23 18:20:25  pvos
 *    Updates for working ProjectContField3D demo
 *
 *    Revision 1.21  2008/09/20 11:34:52  ehan
 *    Fixed some errors
 *
 *    Revision 1.20  2008/09/17 17:29:58  ehan
 *    Fixed some errors to test the LocHexDemo.
 *
 *    Revision 1.19  2008/09/09 15:05:09  sherwin
 *    Updates related to cuved geometries. Normals have been removed from m_metricinfo and replaced with a direct evaluation call. Interp methods have been moved to LibUtilities
 *
 *    Revision 1.18  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *    Revision 1.17  2008/07/19 21:15:38  sherwin
 *    Removed MapTo function, made orientation anticlockwise, changed enum from BndSys to BndLam
 *
 *    Revision 1.16  2008/07/09 11:44:49  sherwin
 *    Replaced GetScaleFactor call with GetConstant(0)
 *
 *    Revision 1.15  2008/07/04 10:19:04  pvos
 *    Some updates
 *
 *    Revision 1.14  2008/06/14 01:20:21  ehan
 *    Clean up the codes
 *
 *    Revision 1.13  2008/06/06 23:23:39  ehan
 *    Added doxygen documentation
 *
 *    Revision 1.12  2008/06/05 20:17:11  ehan
 *    Fixed undefined function GetGtype() in the ASSERTL2().
 *
 *    Revision 1.11  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.10  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.9  2008/05/29 01:02:13  bnelson
 *    Added precompiled header support.
 *
 *    Revision 1.8  2008/04/06 05:59:04  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.7  2008/03/19 06:52:46  ehan
 *    Fixed recent changes of call by reference for the matrix shared pointer. Also fixed name of old functions from Get* to Create*.
 *
 *    Revision 1.5  2008/02/16 05:49:32  ehan
 *    Added PhysDeriv, GenMatrixInfo, standard matrix, and virtual functions.
 *
 *    Revision 1.4  2008/02/05 00:35:25  ehan
 *    Modified coordinate
 *
 *    Revision 1.3  2008/01/31 10:56:50  ehan
 *    Implemented IProductWRTBase, FwdTrans, GetCoord, GetStdMatrix, and GetStdStaticCondMatrix.
 *
 *    Revision 1.2  2007/07/20 00:45:50  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.1  2006/05/04 18:58:45  kirby
 *    *** empty log message ***
 *
 *    Revision 1.8  2006/03/12 07:43:31  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
