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
// Description: Heaxhedral methods
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdHexExp.h>

#ifdef max
#undef max
#endif

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


        StdHexExp::StdHexExp(const  LibUtilities::BasisKey &Ba,
                        const  LibUtilities::BasisKey &Bb,
                        const  LibUtilities::BasisKey &Bc,
                        NekDouble *coeffs,
                        NekDouble *phys)
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


        //////////////////////////////
        // Integration Methods
        //////////////////////////////
        /**
         * @param   fx          ?
         * @param   gy          ?
         * @param   hz          ?
         * @param   inarray     ?
         * @param   outarray    ?
         */
        void StdHexExp::TripleTensorProduct(
                                const Array<OneD, const NekDouble>& fx,
                                const Array<OneD, const NekDouble>& gy,
                                const Array<OneD, const NekDouble>& hz,
                                const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray )
        {

            // Using matrix operation, not sum-factorization.
            // Regarding the 3D array, inarray[k][j][i], x is changing the
            // fastest and z the slowest. Thus, the first x-vector of points
            // refers to the first row of the first stack. The first y-vector
            // refers to the first column of the first stack. The first z-
            // vector refers to the vector of stacks intersecting the first row
            // and first column. So in C++, i refers to column, j to row, and k
            // to stack. Contrasting this with the usual C++ matrix convention,
            // note that i does not refer to a C++ row, nor j to C++ column.

            int nx = fx.num_elements();
            int ny = gy.num_elements();
            int nz = hz.num_elements();

            // Multiply by integration constants...
            // Hadamard multiplication refers to elementwise multiplication of
            // two vectors.

            // Hadamard each row with the first vector (x-vector); the index i
            // is changing the fastest.
            // For each j and k, iterate over each row in all of the stacks at
            // once.
            for (int jk = 0; jk < ny*nz; ++jk)
            {
                Vmath::Vmul(
                    nx,                     // Size of first weight vector
                    &inarray[0] + jk*nx, 1, // Offset and stride of each row-
                                            //  vector (x is changing fastest)
                    fx.get(), 1,            // First weight vector (with stride
                                            //  of 1)
                    &outarray[0] + jk*nx, 1 // Output has same offset and
                                            //  stride as input
                    );
            }

            // Hadamard each column with the second vector (y-vector)
            // For each stack in the 3D-array,  do the following...
            for (int k = 0; k < nz; ++k)
            {
                // Iterate over each column in the current stack
                for (int i = 0; i < nx; ++i)
                {
                    Vmath::Vmul(
                        ny,                     // Size of second weight vector
                        &outarray[0] + i + nx*ny*k, nx,     // Offset and
                                                //  stride of each column-vector
                        gy.get(), 1,            // second weight vector (with
                                                //  stride of 1)
                        &outarray[0] + i + nx*ny*k, nx      // Output has same
                                                //  offset and stride as input
                        );
                }
            }

            // Hadamard each stack-vector with the third vector (z-vector)
            // Iterate over each element in the topmost stack
            for (int ij = 0; ij < nx*ny; ++ij)
            {
                Vmath::Vmul(
                    nz,                         // Size of third weight vector
                    &outarray[0] + ij, nx*ny,   // Offset and stride of each
                                                //  stack-vector
                    hz.get(), 1,                // Third weight vector (with
                                                //  stride of 1)
                    &outarray[0] + ij, nx*ny    // Output has same offset and
                                                //  stride as input
                    );
            }

        }


        /**
         * Inner-Product with respect to the weights: i.e., this is the triple
         * sum of the product of the four inputs over the Hexahedron
         * x-dimension is the row, it is the index that changes the fastest
         * y-dimension is the column
         * z-dimension is the stack, it is the index that changes the slowest
         */
        NekDouble StdHexExp::TripleInnerProduct(
                                     const Array<OneD, const NekDouble>& fxyz,
                                     const Array<OneD, const NekDouble>& wx,
                                     const Array<OneD, const NekDouble>& wy,
                                     const Array<OneD, const NekDouble>& wz
                                      )
        {
            int Qx = wx.num_elements();
            int Qy = wy.num_elements();
            int Qz = wz.num_elements();

            if( fxyz.num_elements() != Qx*Qy*Qz ) {
                cerr << "TripleTetrahedralInnerProduct expected "
                     << fxyz.num_elements()
                     << " quadrature points from the discretized input "
                        "function but got "
                     << Qx*Qy*Qz << " instead." << endl;
            }

            // Sum-factorizing over the stacks
            Array<OneD, NekDouble> A(Qx*Qy, 0.0);
            for( int i = 0; i < Qx; ++i ) {
                for( int j = 0; j < Qy; ++j ) {
                    for( int k = 0; k < Qz; ++k ) {
                        A[i + Qx*j] +=  fxyz[i + Qx*(j + Qy*k)] * wz[k];
                    }
                }
            }

            // Sum-factorizing over the columns
            Array<OneD, NekDouble> b(Qx, 0.0);
            for( int i = 0; i < Qx; ++i ) {
                for( int j = 0; j < Qy; ++j ) {
                    b[i] +=  A[i + Qx*j] * wy[j];
                }
            }

            // Sum-factorizing over the rows
            NekDouble c = 0;
            for( int i = 0; i < Qx; ++i ) {
                c +=  b[i] * wx[i];
            }

            return c;
        }


        /**
         * @param   inarray     ?
         * @param   wx          ?
         * @param   wy          ?
         * @param   wz          ?
         */
        NekDouble StdHexExp::Integral3D(
                                const Array<OneD, const NekDouble>& inarray,
                                const Array<OneD, const NekDouble>& wx,
                                const Array<OneD, const NekDouble>& wy,
                                const Array<OneD, const NekDouble>& wz)
        {
            return TripleInnerProduct( inarray, wx, wy, wz );

        }


        /**
         * @param   inarray     Definition of function to be returned at
         *                      quadrature point of expansion.
         * @returns
         *  \f$\int^1_{-1}\int^1_{-1}\int^1_{-1} u(\xi_1, \xi_2, \xi_3)
         *                      J[i,j,k] d  \xi_1 d \xi_2 d \xi_3 \f$ \n
         *  \f$ = \sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1}
         *          \sum_{k=0}^{Q_3 - 1} u(\xi_{1i}, \xi_{2j},\xi_{3k})
         *          w_{i} w_{j}  w_{k}   \f$ \n
         *  where \f$inarray[i,j, k] = u(\xi_{1i},\xi_{2j}, \xi_{3k}) \f$ \n
         *  and \f$ J[i,j,k] \f$ is the Jacobian evaluated at the quadrature
         *  point.
         */
        NekDouble StdHexExp::v_Integral(
                                const Array<OneD, const NekDouble>& inarray)
        {
            Array<OneD, const NekDouble> w0, w1, w2;

            w0 = m_base[0]->GetW();
            w1 = m_base[1]->GetW();
            w2 = m_base[2]->GetW();

            return Integral3D(inarray, w0, w1, w2);
        }


        //   I/O routine
        void StdHexExp::v_WriteCoeffsToFile(std::ofstream &outfile)
        {
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble> wsp
                            = Array<OneD, NekDouble>(order0*order1*order2, 0.0);

            NekDouble *mat = wsp.get();

            // put coeffs into matrix and reverse order so that r index is
            // fastest for Prism
            Vmath::Zero(order0*order1*order2, mat, 1);

            for(int i = 0, cnt=0; i < order0; ++i)
            {
                for(int j = 0; j < order1-i; ++j)
                {
                    for(int k = 0; k < order2-i-j; ++k, cnt++)
                    {
                        mat[i + order1*(j + order2*k)] = m_coeffs[cnt];
                    }
                }
            }

            outfile <<"Coeffs = [" << " ";

            for(int k = 0; k < order2; ++k)
            {
                for(int j = 0; j < order1; ++j)
                {
                    for(int i = 0; i < order0; ++i)
                    {
                        outfile << mat[i + order0*(j + order1*k)] <<" ";
                    }
                    outfile << std::endl;
                }
            }
            outfile << "]" ;
        }

        bool StdHexExp::v_IsBoundaryInteriorExpansion()
        {
            return (m_base[0]->GetBasisType() == LibUtilities::eModified_A) &&
                   (m_base[1]->GetBasisType() == LibUtilities::eModified_A) &&
                   (m_base[2]->GetBasisType() == LibUtilities::eModified_A);
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
            PhysTensorDeriv(inarray, out_d0, out_d1, out_d2);
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
        void StdHexExp::BwdTrans_SumFacKernel(
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
                ASSERTL1(wsp.num_elements()>=nquad0*nmodes2*(nmodes1+nquad1),
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
                  Array<OneD, NekDouble> &outarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
            int    order1 = m_base[1]->GetNumModes();

            Array<OneD, NekDouble> tmp(inarray.num_elements());
            Array<OneD, NekDouble> wsp(nquad0*nquad1*(nquad2+order0) + 
                                       order0*order1*nquad2);

            MultiplyByQuadratureMetric(inarray,tmp);

            StdHexExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetBdata(),
                                         m_base[2]->GetBdata(),
                                         tmp,outarray,wsp,true,true,true);
        }


        /**
         * Implementation of the sum-factorisation inner product operation.
         * @todo    Implement cases where only some directions are collocated.
         */
        void StdHexExp::IProductWRTBase_SumFacKernel(const Array<OneD, const NekDouble>& base0,
                                                     const Array<OneD, const NekDouble>& base1,
                                                     const Array<OneD, const NekDouble>& base2,
                                                     const Array<OneD, const NekDouble>& inarray,
                                                     Array<OneD, NekDouble> &outarray,
                                                     Array<OneD, NekDouble> &wsp,
                                                     bool doCheckCollDir0,
                                                     bool doCheckCollDir1,
                                                     bool doCheckCollDir2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            int    nmodes0 = m_base[0]->GetNumModes();
            int    nmodes1 = m_base[1]->GetNumModes();
            int    nmodes2 = m_base[2]->GetNumModes();

            bool colldir0 = doCheckCollDir0?(m_base[0]->Collocation()):false;
            bool colldir1 = doCheckCollDir1?(m_base[1]->Collocation()):false;
            bool colldir2 = doCheckCollDir2?(m_base[2]->Collocation()):false;

            if(colldir0 && colldir1 && colldir2)
            {
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,outarray.get(),1);
            }
            else
            {               
                ASSERTL1(wsp.num_elements() >= nmodes0*nquad2*(nquad1+nmodes1),
                         "Insufficient workspace size");

                Array<OneD, NekDouble> tmp0 = wsp;
                Array<OneD, NekDouble> tmp1 = wsp + nmodes0*nquad1*nquad2;

                Blas::Dgemm('T', 'N', nquad1*nquad2, nmodes0, nquad0,
                            1.0, inarray.get(),  nquad0,
                                 base0.get(),    nquad0,
                            0.0, tmp0.get(),     nquad1*nquad2);

                Blas::Dgemm('T', 'N', nquad2*nmodes0, nmodes1, nquad1,
                            1.0, tmp0.get(),     nquad1,
                                 base1.get(),    nquad1,
                            0.0, tmp1.get(),     nquad2*nmodes0);

                Blas::Dgemm('T', 'N', nmodes0*nmodes1, nmodes2, nquad2,
                            1.0, tmp1.get(),     nquad2,
                                 base2.get(),    nquad2,
                            0.0, outarray.get(), nmodes0*nmodes1);
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
            MatrixType mtype;

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
            if (outarray.num_elements() < inarray.num_elements())
            {
                tmp = Array<OneD, NekDouble>(inarray.num_elements());
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


        NekDouble StdHexExp::v_PhysEvaluate(
                const Array<OneD, const NekDouble>& Lcoords)
        {
            return StdExpansion3D::v_PhysEvaluate(Lcoords, m_phys);
        }


        NekDouble StdHexExp::v_PhysEvaluate(
                const Array<OneD, const NekDouble>& Lcoords,
                const Array<OneD, const NekDouble>& physvals)
        {
            return StdExpansion3D::v_PhysEvaluate(Lcoords, physvals);
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
        };


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
            ASSERTL2(i >= 0 && i <= 6, "face id is out of range");
            ASSERTL2(k >= 0 && k <= 1, "basis key id is out of range");
            
            //temporary solution, need to add conditions based on face id
            //also need to add check of the points type
            switch(i)
            {
                case 0:
                case 5:
                    switch(k)
                    {
                        case 0:
                            return GetBasis(0)->GetBasisKey();
                            break;
                        case 1:
                            return GetBasis(1)->GetBasisKey();
                            break;
                    }
                    break;
                case 1:
                case 3:
                    switch(k)
                    {
                        case 0:
                            return GetBasis(0)->GetBasisKey();
                            break;
                        case 1:
                            return GetBasis(2)->GetBasisKey();
                            break;
                    }
                    break;
                case 2:
                case 4:
                    switch(k)
                    {
                        case 0:
                            return GetBasis(1)->GetBasisKey();
                            break;
                        case 1:
                            return GetBasis(2)->GetBasisKey();
                            break;
                    }
                    break;
            }
            
            // Should never get here.
            return LibUtilities::NullBasisKey;
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


        void StdHexExp::v_WriteToFile(std::ofstream &outfile,
                                OutputFormat format,
                                const bool dumpVar,
                                std::string var)
        {
            if(format==eTecplot)
            {
                int  Qx = m_base[0]->GetNumPoints();
                int  Qy = m_base[1]->GetNumPoints();
                int  Qz = m_base[2]->GetNumPoints();

                Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
                eta_x = m_base[0]->GetZ();
                eta_y = m_base[1]->GetZ();
                eta_z = m_base[2]->GetZ();

                if(dumpVar)
                {
                    outfile << "Variables = z1,  z2,  z3";
                    outfile << ", "<< var << std::endl << std::endl;
                }
                outfile << "Zone, I=" << Qx <<", J=" << Qy <<", K=" << Qz
                        <<", F=Point" << endl;

                for(int k = 0; k < Qz; ++k)
                {
                    for(int j = 0; j < Qy; ++j)
                    {
                        for(int i = 0; i < Qx; ++i)
                        {
                            outfile <<  eta_x[i] <<  " " << eta_y[j] << " "
                                    << eta_z[k] << " "
                                    << m_phys[i + Qx*(j + Qy*k)] << endl;
                        }
                    }
                }
            }
            else if (format==eGnuplot)
            {
                int  Qx = m_base[0]->GetNumPoints();
                int  Qy = m_base[1]->GetNumPoints();
                int  Qz = m_base[2]->GetNumPoints();

                Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
                eta_x = m_base[0]->GetZ();
                eta_y = m_base[1]->GetZ();
                eta_z = m_base[2]->GetZ();

                for(int k = 0; k < Qz; ++k)
                {
                    for(int j = 0; j < Qy; ++j)
                    {
                        for(int i = 0; i < Qx; ++i)
                        {
                            outfile <<  eta_x[i] <<  " " << eta_y[j] << " "
                                    << eta_z[k] << " "
                                    << m_phys[i + Qx*(j + Qy*k)] << endl;
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


        /**
         * Only for basis type Modified_A in all directions.
         */
        void StdHexExp::v_GetFaceToElementMap(
            const int                  fid,
            const Orientation          faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray,
            int                        nummodesA,
            int                        nummodesB)
        {
            int i,j;
            const int nummodes0 = m_base[0]->GetNumModes();
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nummodes2 = m_base[2]->GetNumModes();

            ASSERTL1(GetEdgeBasisType(0) == GetEdgeBasisType(1) &&
                     GetEdgeBasisType(0) == GetEdgeBasisType(2),
                     "Method only implemented if BasisType is indentical in "
                     "all directions");
            ASSERTL1(GetEdgeBasisType(0) == LibUtilities::eModified_A,
                      "Method only implemented for Modified_A BasisType");

            if (nummodesA == -1)
            {
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
                }
            }

            int nFaceCoeffs = nummodesA*nummodesB;

            if(maparray.num_elements() != nFaceCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceCoeffs);
            }

            if(signarray.num_elements() != nFaceCoeffs)
            {
                signarray = Array<OneD, int>(nFaceCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nFaceCoeffs, 1 );
            }

            Array<OneD, int> arrayindx(nFaceCoeffs);

            for(i = 0; i < nummodesB; i++)
            {
                for(j = 0; j < nummodesA; j++)
                {
                    if( faceOrient < 9 )
                    {
                        arrayindx[i*nummodesA+j] = i*nummodesA+j;
                    }
                    else
                    {
                        arrayindx[i*nummodesA+j] = j*nummodesB+i;
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
                    offset = nummodes0*nummodes1;
                }
            case 0:
                {
                    jump1 = nummodes0;
                }
                break;
            case 3:
                {
                    offset = nummodes0;
                }
            case 1:
                {
                    jump1 = nummodes0*nummodes1;
                }
                break;
            case 2:
                {
                    offset = 1;
                }
            case 4:
                {
                    jump1 = nummodes0*nummodes1;
                    jump2 = nummodes0;
                }
                break;
            default:
                ASSERTL0(false,"fid must be between 0 and 5");
            }

            for(i = 0; i < nummodesB; i++)
            {
                for(j = 0; j < nummodesA; j++)
                {
                    maparray[ arrayindx[i*nummodesA+j] ]
                                = i*jump1 + j*jump2 + offset;
                }
            }

            if( (faceOrient==6) || (faceOrient==8) ||
                (faceOrient==11) || (faceOrient==12) )
            {

                if(faceOrient<9)
                {
                    for(i = 3; i < nummodesB; i+=2)
                    {
                        for(j = 0; j < nummodesA; j++)
                        {
                            signarray[ arrayindx[i*nummodesA+j] ] *= -1;
                        }
                    }

                    for(i = 0; i < nummodesA; i++)
                    {
                        swap( maparray[i] , maparray[i+nummodesA] );
                        swap( signarray[i] , signarray[i+nummodesA] );
                    }
                }
                else
                {
                    for(i = 0; i < nummodesB; i++)
                    {
                        for(j = 3; j < nummodesA; j+=2)
                        {
                            signarray[ arrayindx[i*nummodesA+j] ] *= -1;
                        }
                    }

                    for(i = 0; i < nummodesB; i++)
                    {
                        swap( maparray[i] , maparray[i+nummodesB] );
                        swap( signarray[i] , signarray[i+nummodesB] );
                    }
                }
            }

            if( (faceOrient==7) || (faceOrient==8) ||
                (faceOrient==10) || (faceOrient==12) )
            {
                if(faceOrient<9)
                {
                    for(i = 0; i < nummodesB; i++)
                    {
                        for(j = 3; j < nummodesA; j+=2)
                        {
                            signarray[ arrayindx[i*nummodesA+j] ] *= -1;
                        }
                    }

                    for(i = 0; i < nummodesB; i++)
                    {
                        swap( maparray[i*nummodesA], maparray[i*nummodesA+1]);
                        swap( signarray[i*nummodesA], signarray[i*nummodesA+1]);
                    }
                }
                else
                {
                    for(i = 3; i < nummodesB; i+=2)
                    {
                        for(j = 0; j < nummodesA; j++)
                        {
                            signarray[ arrayindx[i*nummodesA+j] ] *= -1;
                        }
                    }

                    for(i = 0; i < nummodesA; i++)
                    {
                        swap( maparray[i*nummodesB], maparray[i*nummodesB+1]);
                        swap( signarray[i*nummodesB], signarray[i*nummodesB+1]);
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
        int StdHexExp::v_GetVertexMap(const int localVertexId)
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

            if(maparray.num_elements()!=nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }

            if(signarray.num_elements() != nEdgeIntCoeffs)
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

            int IdxRange [3][2];

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

            if(maparray.num_elements()!=nFaceIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceIntCoeffs);
            }

            if(signarray.num_elements() != nFaceIntCoeffs)
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

            int nummodesA;
            int nummodesB;

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
                    if( faceOrient < 9 )
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
                        if( (((int) faceOrient)-5) % 2 )
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

                        if( (((int) faceOrient)-5) % 2 )
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
                        if( (((int) faceOrient)-5) % 2 )
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

                        if( (((int) faceOrient)-5) % 2 )
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
                        if( (((int) faceOrient)-5) % 4 > 1 )
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

                        if( (((int) faceOrient)-5) % 4 > 1 )
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
                        if( (((int) faceOrient)-5) % 4 > 1 )
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

                        if( (((int) faceOrient)-5) % 4 > 1 )
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

            if(outarray.num_elements() != nIntCoeffs)
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

            if(outarray.num_elements()!=nBndCoeffs)
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


        void StdHexExp::v_LaplacianMatrixOp_MatFree(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            ASSERTL0(false,"StdHexExp::v_LaplacianMatrixOp_MatFree needs fixing");
/*            if(mkey.GetNvariableLaplacianCoefficients() == 0)
            {
                // This implementation is only valid when there are no coefficients
                // associated to the Laplacian operator
                int       nquad0  = m_base[0]->GetNumPoints();
                int       nquad1  = m_base[1]->GetNumPoints();
                int       nqtot   = nquad0*nquad1;
                int       nmodes0 = m_base[0]->GetNumModes();
                int       nmodes1 = m_base[1]->GetNumModes();
                int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);

                const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
                const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
                const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();

                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(3*wspsize);
                Array<OneD,NekDouble> wsp1(wsp0+wspsize);
                Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);

                if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                {
                    // LAPLACIAN MATRIX OPERATION
                    // wsp0 = u       = B   * u_hat
                    // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                    // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                    BwdTrans_SumFacKernel(base0,base1,inarray,wsp0,wsp1,true,true);
                    StdExpansion2D::PhysTensorDeriv(wsp0,wsp1,wsp2);
                }
                else
                {
                    StdExpansion2D::PhysTensorDeriv(inarray,wsp1,wsp2);
                }

                // wsp1 = k = wsp1 * w0 * w1
                // wsp2 = l = wsp2 * w0 * w1
                MultiplyByQuadratureMetric(wsp1,wsp1);
                MultiplyByQuadratureMetric(wsp2,wsp2);

                // outarray = m = (D_xi1 * B)^T * k
                // wsp1     = n = (D_xi2 * B)^T * l
                IProductWRTBase_SumFacKernel(dbase0,base1,wsp1,outarray,wsp0,false,true);
                IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp1,    wsp0,true,false);

                // outarray = outarray + wsp1
                //          = L * u_hat
                Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);
            }
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }*/
        }

        void StdHexExp::v_HelmholtzMatrixOp_MatFree(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey)
        {
            ASSERTL0(false,"StdHexExp::v_HelmholtzMatrixOp_MatFree needs fixing");
/*            int       nquad0  = m_base[0]->GetNumPoints();
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

            // wsp1 = k = wsp1 * w0 * w1
            // wsp2 = l = wsp2 * w0 * w1
            MultiplyByQuadratureMetric(wsp1,wsp1);
            MultiplyByQuadratureMetric(wsp2,wsp2);

            // wsp1 = m = (D_xi1 * B)^T * k
            // wsp0 = n = (D_xi2 * B)^T * l
            IProductWRTBase_SumFacKernel(dbase0,base1,wsp1,wsp0,wsp3,false,true);
            IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp1,wsp3,true,false);

            // outarray = lambda * outarray + (wsp0 + wsp1)
            //          = (lambda * M + L ) * u_hat
            Vmath::Vstvpp(m_ncoeffs,lambda,&outarray[0],1,&wsp1[0],1,&wsp0[0],1,&outarray[0],1);
*/        }


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


        void StdHexExp::MultiplyByQuadratureMetric(const Array<OneD, const NekDouble>& inarray,
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
                Vmath::Smul(nquad0, w1[i%nquad2], outarray.get()+i*nquad0, 1,
                            outarray.get()+i*nquad0, 1);
            }

            for(i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nq01, w2[i], outarray.get()+i*nq01, 1,
                            outarray.get()+i*nq01, 1);
            }
        }

    }
}

