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
// Description: Header field for tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdTetExp.h>

namespace Nektar
{
    namespace StdRegions
    {
        /**
         *
         */
        StdTetExp::StdTetExp()
        {
        }


        /**
         * @param   Ba          ?
         * @param   Bb          ?
         * @param   Bc          ?
         */
        StdTetExp::StdTetExp(   const LibUtilities::BasisKey &Ba,
                                const LibUtilities::BasisKey &Bb,
                                const LibUtilities::BasisKey &Bc )
            : StdExpansion3D(StdTetData::getNumberOfCoefficients(
                        Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                    Ba, Bb, Bc)
        {
            if(Ba.GetNumModes() >  Bb.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher than order "
                                "in 'b' direction");
            }
            if(Ba.GetNumModes() >  Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher than order "
                                "in 'c' direction");
            }
            if(Bb.GetNumModes() >  Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'b' direction is higher than order "
                                "in 'c' direction");
            }
        }


        /**
         * @param   T           StdTetExp object to copy.
         */
        StdTetExp::StdTetExp(const StdTetExp &T):
            StdExpansion3D(T)
        {
        }


        /**
         *
         */
        StdTetExp::~StdTetExp()
        {
        }

        //////////////////////////////
        // Integration Methods
        //////////////////////////////
        void StdTetExp::TripleTensorProduct(
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
            // once
            for (int jk = 0; jk < ny*nz; ++jk)
            {
                Vmath::Vmul(
                        nx,                 // Size of first weight vector
                        &inarray[0] + jk*nx, 1, // Offset and stride of each
                                            //  row-vector (x is changing
                                            //  fastest)
                        fx.get(), 1,        // First weight vector (with stride
                                            //  of 1)
                        &outarray[0] + jk*nx, 1     // Output has same offset
                                            //  and stride as input
                        );
            }

            // Hadamard each column with the second vector (y-vector)
            // For each stack   in the 3D-array, do the following...
            for (int k = 0; k < nz; ++k)
            {
                // Iterate over each column in the current stack
                for (int i = 0; i < nx; ++i)
                {
                    Vmath::Vmul(
                        ny,                 // Size of second weight vector
                        &outarray[0] + i + nx*ny*k, nx, // Offset and stride of
                                            //  each column-vector
                        gy.get(), 1,        // second weight vector (with
                                            //  stride of 1)
                        &outarray[0] + i + nx*ny*k, nx  // Output has same
                                            //  offset and stride as input
                        );
                }
            }

            // Hadamard each stack-vector with the third vector (z-vector)
            // Iterate over each element in the topmost stack
            for (int ij = 0; ij < nx*ny; ++ij)
            {
                Vmath::Vmul(
                        nz,                 // Size of third weight vector
                        &outarray[0] + ij, nx*ny,   // Offset and stride of
                                            //  each stack-vector
                        hz.get(), 1,        // Third weight vector (with stride
                                            //  of 1)
                        &outarray[0] + ij, nx*ny    // Output has same offset
                                            //  and stride as input
                        );
            }

        }


        /**
         * Inner-Product with respect to the weights: i.e., this is the triple
         * sum of the product of the four inputs over the Hexahedron.
         * x-dimension is the row, it is the index that changes the fastest
         * y-dimension is the column
         * z-dimension is the stack, it is the index that changes the slowest
         */
        NekDouble StdTetExp::TripleInnerProduct(
                                const Array<OneD, const NekDouble>& fxyz,
                                const Array<OneD, const NekDouble>& wx,
                                const Array<OneD, const NekDouble>& wy,
                                const Array<OneD, const NekDouble>& wz)
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


        NekDouble StdTetExp::Integral3D(
                                const Array<OneD, const NekDouble>& inarray,
                                const Array<OneD, const NekDouble>& wx,
                                const Array<OneD, const NekDouble>& wy,
                                const Array<OneD, const NekDouble>& wz)
        {
            return TripleInnerProduct( inarray, wx, wy, wz );

        }



        /**
         * @param   inarray     definition of function to be returned at
         *                      quadrature point of expansion.
         * @returns \f$\int^1_{-1}\int^1_{-1}\int^1_{-1}
         * u(\eta_1, \eta_2, \eta_3) J[i,j,k] d \eta_1 d \eta_2 d \eta_3 \f$ \n
         * = \f$\sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1} \sum_{k=0}^{Q_3 - 1}
         *   u(\eta_{1i}^{0,0}, \eta_{2j}^{1,0},\eta_{3k}^{2,0})w_{i}^{1,0}
         *   \hat w_{j}^{1,0} \hat w_{k}^{2,0}    \f$ \n
         * where
         * \f$inarray[i,j, k]  = u(\eta_{1i},\eta_{2j}, \eta_{3k}) \f$, \n
         * \f$\hat w_{j}^{1,0} = \frac {w_{j}^{1,0}} {2}, \hat w_{k}^{2,0}
         *                     = \frac{w_{k}^{2,0}} {4} \f$ \n
         * and \f$ J[i,j,k] \f$ is the Jacobian evaluated at the quadrature
         * point.
         */
        NekDouble StdTetExp::Integral(
                                const Array<OneD, const NekDouble>& inarray)
        {
            // Using implementation from page 145 of Spencer Sherwin's book
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            // Get the point distributions:
            // x is assumed to be Gauss-Lobatto-Legendre (includes -1 and +1)
            Array<OneD, const NekDouble> y,z,wx,wy,wz;
            wx = m_base[0]->GetW();
            m_base[1]->GetZW(y,wy);
            m_base[2]->GetZW(z,wz);

            Array<OneD, NekDouble> wy_hat = Array<OneD, NekDouble>(Qy, 0.0);
            Array<OneD, NekDouble> wz_hat = Array<OneD, NekDouble>(Qz, 0.0);

            // Convert wy into wy_hat, which includes the 1/2 scale factor.
            // Nothing else need be done if the point distribution is
            // Jacobi (1,0) since (1-eta_y) is aready factored into the weights.
            switch(m_base[1]->GetPointsType())
            {
            case LibUtilities::eGaussLobattoLegendre:   // Legendre inner product (Falls-through to next case)
            case LibUtilities::eGaussRadauMLegendre:    // (0,0) Jacobi Inner product
                for(int j = 0; j < Qy; ++j)
                {
                    wy_hat[j] = 0.5*(1.0 - y[j]) * wy[j];
                }
                break;

            case LibUtilities::eGaussRadauMAlpha1Beta0: // (1,0) Jacobi Inner product
                Vmath::Smul( Qy, 0.5, (NekDouble *)wy.get(), 1, wy_hat.get(), 1 );
                break;
            }

            // Convert wz into wz_hat, which includes the 1/4 scale factor.
            // Nothing else need be done if the point distribution is
            // Jacobi (2,0) since (1-eta_z)^2 is aready factored into the weights.
            // Note by coincidence, xi_z = eta_z    (xi_z = z according to our notation)
            switch(m_base[2]->GetPointsType())
            {
            case LibUtilities::eGaussLobattoLegendre:   // Legendre inner product (Falls-through to next case)
            case LibUtilities::eGaussRadauMLegendre:    // (0,0) Jacobi Inner product
                for(int k = 0; k < Qz; ++k)
                {
                    wz_hat[k] = 0.25*(1.0 - z[k])*(1.0 - z[k]) * wz[k];
                }
                break;

            case LibUtilities::eGaussRadauMAlpha2Beta0: // (2,0) Jacobi Inner product
                Vmath::Smul( Qz, 0.25,(NekDouble *)wz.get(), 1, wz_hat.get(), 1 );
                break;
            }

            return Integral3D( inarray, wx, wy_hat, wz_hat );
        }


        /**
         * @note 'r' (base[2]) runs fastest in this element
         *
         * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) 
         *  =  \sum_{m(pqr)} \hat u_{pqr} 
         *          \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
         *
         * Backward transformation is three dimensional tensorial expansion
         * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x}
         * \psi_p^a (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{pq}^b
         * (\xi_{2j}) \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr}
         * \psi_{pqr}^c (\xi_{3k}) \rbrace} \rbrace}. \f$
         * And sumfactorizing step of the form is as:\\
         *
         * \f$ f_{pq} (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr}
         * \psi_{pqr}^c (\xi_{3k}),\\
         *
         * g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y} \psi_{pq}^b
         * (\xi_{2j}) f_{pq} (\xi_{3k}),\\
         *
         * u(\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x}
         * \psi_{p}^a (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).  \f$
         */
        void StdTetExp::v_BwdTrans(
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
                bool doMatOp = NekOptimize::ElementalOptimization<eStdTetExp, NekOptimize::eBwdTrans, 3>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes(),m_base[2]->GetNumModes());

                if(doMatOp)
                {
                    StdExpansion::BwdTrans_MatOp(inarray,outarray);
                }
                else
                {
                    StdTetExp::v_BwdTrans_SumFac(inarray,outarray);
                }  
            }
        }

        
        /**
         * Sum-factorisation implementation of the BwdTrans operation.
         */
        void StdTetExp::v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
        {
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();

            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();
            Array<OneD, const NekDouble> base2 = m_base[2]->GetBdata();

            Array<OneD, NekDouble > tmp(nquad2*order0*(order1+1)/2);
            Array<OneD, NekDouble > tmp1(nquad2*nquad1*order0);

            int i,j, mode,mode1, cnt;

            // Perform summation over '2' direction
            mode = mode1 = cnt = 0;
            for(i = 0; i < order0; ++i)
            {
                for(j = 0; j < order1-i ; ++j, ++cnt)
                {
                    Blas::Dgemv('N', nquad2,order2-i-j,1.0,
                                base2.get()+mode*nquad2,
                                nquad2,&inarray[0]+mode1,1,0.0,
                                &tmp[0]+cnt*nquad2,1);
                    mode  += order2-i-j;
                    mode1 += order2-i-j;
                }
                //increment mode1 in case order1!=order2
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
                Blas::Dgemm('N','T',nquad1,nquad2,order1-i,1.0,
                            base1.get()+mode*nquad1,nquad1,
                            &tmp[0]+mode*nquad2,nquad2,0.0,
                            &tmp1[0]+i*nquad1*nquad2,nquad1);
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
            Blas::Dgemm('N','T', nquad0,nquad1*nquad2,order0,1.0,
                        base0.get(),nquad0, &tmp1[0], nquad1*nquad2,
                        0.0, &outarray[0], nquad0);
        }


        /**
         * @param   inarray     array of physical quadrature points to be
         *                      transformed.
         * @param   outarray    updated array of expansion coefficients.
         */
        void StdTetExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
        {
            v_IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetExpansionType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs, outarray);
            DNekVec out(m_ncoeffs, outarray, eWrapper);

            out = (*matsys)*in;
        }


        /**
         * \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
         * \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
         * \psi_{p}^{a} (\eta_{1i}) \psi_{pq}^{b} (\eta_{2j}) \psi_{pqr}^{c} (\eta_{3k})
         * w_i w_j w_k u(\eta_{1,i} \eta_{2,j} \eta_{3,k})
         * J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\eta_{1,i})
         * \sum_{j=0}^{nq_1} \psi_{pq}^b(\eta_{2,j}) \sum_{k=0}^{nq_2} \psi_{pqr}^c u(\eta_{1i},\eta_{2j},\eta_{3k})
         * J_{i,j,k} \end{array} \f$ \n
         *
         * where
         *
         * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\eta_1) \psi_{pq}^b (\eta_2) \psi_{pqr}^c (\eta_3) \f$
         *
         * which can be implemented as \n
         * \f$f_{pqr} (\xi_{3k}) = \sum_{k=0}^{nq_3} \psi_{pqr}^c u(\eta_{1i},\eta_{2j},\eta_{3k})
         * J_{i,j,k} = {\bf B_3 U}   \f$ \n
         * \f$ g_{pq} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{pq}^b (\xi_{2j}) f_{pqr} (\xi_{3k})  = {\bf B_2 F}  \f$ \n
         * \f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{pq} (\xi_{3k})  = {\bf B_1 G} \f$
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
                bool doMatOp = NekOptimize::ElementalOptimization<eStdTetExp, NekOptimize::eIProductWRTBase, 3>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes(),m_base[2]->GetNumModes());
                
                if(doMatOp)
                {
                    StdTetExp::v_IProductWRTBase_MatOp(inarray,outarray);
                }
                else
                {
                    StdTetExp::v_IProductWRTBase_SumFac(inarray,outarray);
                }  
            }

        }


        void StdTetExp::v_IProductWRTBase_MatOp (
                                const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetExpansionType(),*this);
            DNekMatSharedPtr& iprodmat = GetStdMatrix(iprodmatkey);            
            
            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);

        }
        

        /**
         * @param   inarray     Function evaluated at physical collocation
         *                      points.
         * @param   outarray    Inner product with respect to each basis
         *                      function over the element.
         */
        void StdTetExp::v_IProductWRTBase_SumFac (
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> & outarray)
        {
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();

            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble > tmp (nquad0*nquad1*nquad2);
            Array<OneD, NekDouble > wsp (nquad1*nquad2*order0 
                                            + nquad2*order0*(order1+1)/2);
            
            MultiplyByQuadratureMetric(inarray, tmp);
            
            StdTetExp::IProductWRTBase_SumFacKernel(
                    m_base[0]->GetBdata(),
                    m_base[1]->GetBdata(),
                    m_base[2]->GetBdata(),
                    tmp, outarray, wsp, true, true, true);
        }
        
        
        void StdTetExp::IProductWRTBase_SumFacKernel (
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
            Blas::Dgemm('T','N', nquad1*nquad2, order0, nquad0, 1.0,
                        &inarray[0], nquad0, base0.get(), nquad0, 0.0,
                        &tmp1[0], nquad1*nquad2);


            // Inner product with respect to the '1' direction
            for(mode=i=0; i < order0; ++i)
            {
                Blas::Dgemm('T','N',nquad2,order1-i,nquad1,1.0,
                            &tmp1[0]+i*nquad1*nquad2, nquad1,
                            base1.get()+mode*nquad1, nquad1,
                            0.0, &tmp2[0]+mode*nquad2, nquad2);
                mode  += order1-i;
            }


            // fix for modified basis for base singular vertex
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                //base singular vertex and singular edge (1+b)/2
                //(1+a)/2 components (makes tmp[nquad2] entry into (1+b)/2)
                Blas::Dgemv('T', nquad1,nquad2, 1.0, &tmp1[0]+nquad1*nquad2,
                            nquad1, base1.get()+nquad1,1, 1.0, &tmp2[nquad2],1);
            }


            // Inner product with respect to the '2' direction
            mode = mode1 = cnt = 0;
            for(i = 0; i < order0; ++i)
            {
                for(j = 0; j < order1-i ; ++j, ++cnt)
                {
                    Blas::Dgemv('T', nquad2, order2-i-j,1.0,
                                base2.get()+mode*nquad2,
                                nquad2,&tmp2[0]+cnt*nquad2, 1,
                                0.0, &outarray[0]+mode1,1);
                    mode  += order2-i-j;
                    mode1 += order2-i-j;
                }

                //increment mode1 in case order1!=order2
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


        void StdTetExp::MultiplyByQuadratureMetric(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray)
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
            // Legendre inner product.
            case LibUtilities::eGaussLobattoLegendre:

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

            // (1,0) Jacobi Inner product.
            case LibUtilities::eGaussRadauMAlpha1Beta0:
                for(j = 0; j < nquad2; ++j)
                {
                    for(i = 0; i < nquad1; ++i)
                    {
                        Blas::Dscal(nquad0,0.5*w1[i], &outarray[0]+i*nquad0 +
                                    j*nquad0*nquad1,1);
                    }
                }
                break;
            }

            switch(m_base[2]->GetPointsType())
            {
            // Legendre inner product.
            case LibUtilities::eGaussLobattoLegendre:
                for(i = 0; i < nquad2; ++i)
                {
                    Blas::Dscal(nquad0*nquad1,0.25*(1-z2[i])*(1-z2[i])*w2[i],
                                &outarray[0]+i*nquad0*nquad1,1);
                }
                break;
            // (2,0) Jacobi inner product.
            case LibUtilities::eGaussRadauMAlpha2Beta0:
                for(i = 0; i < nquad2; ++i)
                {
                    Blas::Dscal(nquad0*nquad1, 0.25*w2[i],
                                &outarray[0]+i*nquad0*nquad1, 1);
                }
                break;
            }
        }
        
        
        void StdTetExp::FillMode(
                            const int mode, Array<OneD, NekDouble> &outarray)
        {
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                ASSERTL0(false,"This function will not work with modified "
                               "basis since we have not dealt with singular "
                               "vertces/edges");
            }

            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;


            // Index map from the rectangle to the triangle
            Array<OneD, int> mode_pq  = Array<OneD, int>( (P+1)*(Q+1), -1 );
            for( int p = 0, m = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q - p; ++q, ++m ) {
                    mode_pq[q + (Q+1)*p] = m;
                }
            }

            // Create an index map from the hexahedron to the tetrahedron.
            // The actual index is too difficult to compute explicitly.
            Array<OneD, int> mode_pqr
                                = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, m = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q - p; ++q ) {
                    for( int r = 0; r <= R - p - q; ++r, ++m ) {
                        mode_pqr[r + (R+1)*(q + (Q+1)*p)] = m;
                    }
                }
            }

            // Find the pqr matching the provided mode
            int mode_p=0, mode_q=0, mode_r=0;
            for( int p = 0, m = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q - p; ++q ) {
                    for( int r = 0; r <= R - p - q; ++r, ++m ) {
                        if( m == mode ) {
                            mode_p = p;
                            mode_q = q;
                            mode_r = r;
                        }
                    }
                }
            }

            const Array<OneD, const NekDouble>& bx = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble>& by = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble>& bz = m_base[2]->GetBdata();

            int p = mode_p, q = mode_q, r = mode_r;

            // Determine the index for specifying which mode to use in the basis
            int sigma_p   = Qx*p;
            int sigma_pq  = Qy*mode_pq[q + (Q+1)*p];
            int sigma_pqr = Qz*mode_pqr[r + (R+1)*(q + (Q+1)*p)];
            int sigma = Qx*Qy*Qz*mode;


            // Compute tensor product of inarray with the 3 basis functions
            Array<OneD, NekDouble> g_pqr
                                = Array<OneD, NekDouble>( Qx*Qy*Qz, 0.0 );
            for( int k = 0; k < Qz; ++k ) {
                for( int j = 0; j < Qy; ++j ) {
                    for( int i = 0; i < Qx; ++i ) {
                        int s = i + Qx*(j + Qy*(k + Qz*mode));
                        outarray[s] =
                            bx[i + sigma_p] *
                            by[j + sigma_pq] *
                            bz[k + sigma_pqr];
                    }
                }
            }

        }



        //-----------------------------
        /// Differentiation Methods
        //-----------------------------

    /**
            \brief Calculate the derivative of the physical points

        The derivative is evaluated at the nodal physical points.
            Derivatives with respect to the local Cartesian coordinates

        \f$\begin{Bmatrix} \frac {\partial} {\partial \xi_1}  \\ \frac {\partial} {\partial \xi_2} \\ \frac {\partial} {\partial \xi_3}  \end{Bmatrix}  = \begin{Bmatrix} \frac 4 {(1-\eta_2)(1-\eta_3)} \frac \partial {\partial \eta_1} \\
        \frac {2(1+\eta_1)} {(1-\eta_2)(1-\eta_3)} \frac \partial {\partial \eta_1} + \frac 2 {1-\eta_3} \frac \partial {\partial \eta_3} \\
        \frac {2(1 + \eta_1)} {2(1 - \eta_2)(1-\eta_3)} \frac \partial {\partial \eta_1} + \frac {1 + \eta_2} {1 - \eta_3} \frac \partial {\partial \eta_2} + \frac \partial {\partial \eta_3}
        \end{Bmatrix}\f$

        **/
        // PhysDerivative implementation based on Spen's book page 152.
        void StdTetExp::PhysDeriv(const Array<OneD, const NekDouble>& u_physical,
                                  Array<OneD, NekDouble> &out_dxi1,
                                  Array<OneD, NekDouble> &out_dxi2,
                                  Array<OneD, NekDouble> &out_dxi3 )
        {

            int    Qx = m_base[0]->GetNumPoints();
            int    Qy = m_base[1]->GetNumPoints();
            int    Qz = m_base[2]->GetNumPoints();

            // Compute the physical derivative
            Array<OneD, NekDouble> out_dEta1(Qx*Qy*Qz,0.0), out_dEta2(Qx*Qy*Qz,0.0), out_dEta3(Qx*Qy*Qz,0.0);
            PhysTensorDeriv(u_physical, out_dEta1, out_dEta2, out_dEta3);

            Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
            eta_x = m_base[0]->GetZ();
            eta_y = m_base[1]->GetZ();
            eta_z = m_base[2]->GetZ();


            for(int k=0, n=0; k<Qz; ++k)
                for(int j=0; j<Qy; ++j){
                    for(int i=0; i<Qx; ++i, ++n){
                        {
                            out_dxi1[n] = 4.0 / ((1.0 - eta_y[j])*(1.0 - eta_z[k]))*out_dEta1[n];
                            out_dxi2[n] = 2.0*(1.0 + eta_x[i]) / ((1.0-eta_y[j])*(1.0 - eta_z[k]))*out_dEta1[n]  +  2.0/(1.0 - eta_z[k])*out_dEta2[n];
                            out_dxi3[n] = 2.0*(1.0 + eta_x[i]) / ((1.0 - eta_y[j])*(1.0 - eta_z[k]))*out_dEta1[n] + (1.0 + eta_y[j]) / (1.0 - eta_z[k])*out_dEta2[n] + out_dEta3[n];

                            //cout << "eta_x["<<i<<"] = " <<  eta_x[i] << ",  eta_y["<<j<<"] = " << eta_y[j] << ", eta_z["<<k<<"] = " <<eta_z[k] << endl;
                            //cout << "out_dEta1["<<n<<"] = " << out_dEta1[n] << ",  out_dEta2["<<n<<"] = " << out_dEta2[n] << ", out_dEta3["<<n<<"] = " <<out_dEta3[n] << endl;
                            //cout << "out_dxi1["<<n<<"] = " << out_dxi1[n] << ",  out_dxi2["<<n<<"] = " << out_dxi2[n] << ", out_dxi3["<<n<<"] = " << out_dxi3[n] << endl;

                        }
                    }
                }
        }


        ///////////////////////////////
        /// Evaluation Methods
        ///////////////////////////////



        void StdTetExp::GetCoords( Array<OneD, NekDouble> & xi_x, Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z)
        {
            Array<OneD, const NekDouble> eta_x = m_base[0]->GetZ();
            Array<OneD, const NekDouble> eta_y = m_base[1]->GetZ();
            Array<OneD, const NekDouble> eta_z = m_base[2]->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates: eta --> xi
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


        void StdTetExp::GetFaceToElementMap(const int fid, const FaceOrientation faceOrient,
                                            Array<OneD, unsigned int> &maparray,
                                            Array<OneD, int>& signarray)
        {

            int i,j;
            const int nummodes0 = m_base[0]->GetNumModes();
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nummodes2 = m_base[2]->GetNumModes();
            int nummodesA, nummodesB, P, Q;

            const LibUtilities::BasisType bType0 = GetEdgeBasisType(0);
            const LibUtilities::BasisType bType1 = GetEdgeBasisType(1);
            const LibUtilities::BasisType bType2 = GetEdgeBasisType(3);

            ASSERTL1( (bType0==LibUtilities::eModified_A) && (bType1==LibUtilities::eModified_B)&& (bType2==LibUtilities::eModified_C),
                      "Method only implemented for Modified_A BasisType (x direction), Modified_B BasisType (y direction), and Modified_C BasisType(z direction)");

            int nFaceCoeffs = 0;
            if (fid == 0)
            {
                nummodesA = nummodes0;
                nummodesB = nummodes1;
            }
            else if (fid == 1)
            {
                nummodesA = nummodes0;
                nummodesB = nummodes2;
            }
            else if (fid == 2) {
                nummodesA = nummodes1;
                nummodesB = nummodes2;
            }
            else
            {
                nummodesA = nummodes1;
                nummodesB = nummodes2;
            }
            P = nummodesA - 1;
            Q = nummodesB - 1;
            nFaceCoeffs = Q+1 + (P*(1 + 2*Q - P))/2;

            // Allocate the map array and sign array; set sign array to ones (+)
            if(maparray.num_elements() != nFaceCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceCoeffs,1);
            }

            if(signarray.num_elements() != nFaceCoeffs)
            {
                signarray = Array<OneD, int>(nFaceCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nFaceCoeffs, 1 );
            }

            // A mapping array to account for reversal of modes.
            Array<OneD, int> arrayindex(nFaceCoeffs,-1);
            
            int limit = 0;      // Inner loop index limit
            int offsetA = 0;    // Element Offset in B direction
            int offsetB = 0;    // Element Offset in A direction
            int jumpA = 0;      // Jump in the A direction
            int jumpB = 0;      // Jump in the B direction
            int cnt = 0;        // Face index offset
            int n0 = 0;         // Computed depth in face orthogonal dirn in
                                //   slow face coordinate
            int n1 = 0;         // Computed depth in face orthogonal dirn in
                                //   fast face coordinate

            // Create array indexing the modes at a face level. This handles
            // face orientation.
            if (faceOrient < 4) // Not transposed
            {
                for(int b = 0; b < nummodesB; ++b)
                {
                    limit = (b > nummodesB - nummodesA) ? nummodesB - b : nummodesA;
                    jumpB = limit;
                    for(int a = 0; a < limit; ++a)
                    {
                        arrayindex[a + offsetB] = a + offsetB;
                    }
                    offsetB += jumpB;
                }
            }
            else // Transposed
            {
                for(int a = 0; a < nummodesA; ++a)
                {
                    limit = (a > nummodesA - nummodesB) ? nummodesA - a : nummodesB;
                    jumpA = limit;                    
                    for(int b = 0; b < nummodesB - a; ++b)
                    {
                        arrayindex[b + offsetA] = b + offsetA;
                    }
                    offsetA += jumpA;
                }
            }

            // Check the arrayindex was created properly and all entries filled.
            for (int i = 0; i < arrayindex.num_elements(); ++i)
            {
                ASSERTL0(arrayindex[i] != -1, "arrayindex is not set up properly.");
            }
            
            // Reset counters
            limit = 0;
            offsetA = 0;
            offsetB = 0;
            jumpB = 0;
            cnt = 0;

            // Now generate the indices of the modes at an element level.
            // Each face is handled separately. The elemental mode index is
            // computed for each face mode.
            switch(fid)
            {

            // Base triangle
            case 0:
                for(int b = 0; b < nummodesB; ++b) 
                {
                    limit    = (b > nummodesB - nummodesA) ? nummodesB - b : nummodesA;
                    jumpB    = limit;
                    for(int a = 0; a < limit; ++a) 
                    {
                        maparray[ arrayindex[a + cnt] ] = a + offsetB;
                    }
                    cnt     += limit;
                    offsetB += jumpB;
                }
                break;

            // Front triangle
            case 1:
                for(int b = 0; b < nummodesB; ++b) 
                {
                    n0       = (b > nummodesB - nummodes1) ? nummodesB - b : nummodes1;
                    limit    = (b > nummodesB - nummodesA) ? nummodesB - b : nummodesA;
                    
                    for(int a = 0; a < limit; ++a) 
                    {
                        maparray[ arrayindex[a + cnt] ] = a + offsetB;
                    }
                    cnt     += limit;
                    jumpB    = (limit * (1 + 2*n0 - limit)) / 2;
                    offsetB += jumpB;
                }
                break;

            // Slanted triangle
            case 2:
                for(int b = 0; b < nummodesB; ++b) 
                {
                    n0       = (b > nummodesB - nummodes0) ? nummodesB - b : nummodes0;
                    limit    = (b > nummodesB - nummodesA) ? nummodesB - b : nummodesA;
                    offsetA  = -1;
                    for(int a = 0; a < limit; ++a) 
                    {
                        n1       = (a > limit - n0) ? limit - a : n0;
                        offsetA += n1;
                        maparray[ arrayindex[a + cnt] ] = offsetB + offsetA;
                    }
                    cnt     += limit;
                    jumpB    = (n0 * (1 + 2*limit - n0)) / 2;
                    offsetB += jumpB;
                }
                break;

            // Rear triangle
            case 3:
                for(int b = 0; b < nummodesB; ++b) 
                {
                    n0       = (b > nummodesB - nummodes0) ? nummodesB - b : nummodes0;
                    limit    = (b > nummodesB - nummodesA) ? nummodesB - b : nummodesA;
                    offsetA  = 0;
                    for(int a = 0; a < limit; ++a) 
                    {
                        maparray[ arrayindex[a + cnt] ] = offsetB + offsetA;
                        n1       = (a > limit - n0) ? limit - a : n0;
                        offsetA += n1;
                    }
                    cnt     += limit;
                    jumpB    = (n0 * (1 + 2*limit - n0)) / 2;
                    offsetB += jumpB;
                }
                break;
            }

            // Finally set up sign arrays.
            /// @todo   Correctly set up sign arrays.
            if( (faceOrient==1) || (faceOrient==3) ||
                (faceOrient==6) || (faceOrient==7) )
            {

                if(faceOrient<4)
                {
                    for(i = 3; i < nummodesB; i+=2)
                    {
                        for(j = 0; j < nummodesA; j++)
                        {
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    }

                    for(i = 0; i < nummodesA; i++)
                    {
                        //swap( maparray[i] , maparray[i+nummodesA] );
                        //swap( signarray[i] , signarray[i+nummodesA] );
                    }
                }
                else
                {
                    for(i = 0; i < nummodesB; i++)
                    {
                        for(j = 3; j < nummodesA; j+=2)
                        {
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    }

                    for(i = 0; i < nummodesB; i++)
                    {
                        //swap( maparray[i] , maparray[i+nummodesB] );
                        //swap( signarray[i] , signarray[i+nummodesB] );
                    }
                }
            }

            if( (faceOrient==2) || (faceOrient==3) ||
                (faceOrient==5) || (faceOrient==7) )
            {
                if(faceOrient<4)
                {
                    for(i = 0; i < nummodesB; i++)
                    {
                        for(j = 3; j < nummodesA; j+=2)
                        {
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    }

                    for(i = 0; i < nummodesB; i++)
                    {
                        //swap( maparray[i*nummodesA] , maparray[i*nummodesA+1] );
                        //swap( signarray[i*nummodesA] , signarray[i*nummodesA+1] );
                    }
                }
                else
                {
                    for(i = 3; i < nummodesB; i+=2)
                    {
                        for(j = 0; j < nummodesA; j++)
                        {
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    }

                    for(i = 0; i < nummodesA; i++)
                    {
                        //swap( maparray[i*nummodesB] , maparray[i*nummodesB+1] );
                        //swap( signarray[i*nummodesB] , signarray[i*nummodesB+1] );
                    }
                }
            }

        }


        void StdTetExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
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
                outfile << "Zone, I=" << Qx <<", J=" << Qy <<", K=" << Qz <<", F=Point" << std::endl;

                for(int k = 0; k < Qz; ++k)
                {
                    for(int j = 0; j < Qy; ++j)
                    {
                        for(int i = 0; i < Qx; ++i)
                        {
                            //outfile << 0.5*(1+z0[i])*(1.0-z1[j])-1 <<  " " << z1[j] << " " << m_phys[j*nquad0+i] << std::endl;
                            outfile <<  (eta_x[i] + 1.0) * (1.0 - eta_y[j]) * (1.0 - eta_z[k]) / 4  -  1.0 <<  " " << eta_z[k] << " " << m_phys[i + Qx*(j + Qy*k)] << std::endl;
                        }
                    }
                }
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }


        //   I/O routine
        void StdTetExp::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble> wsp  = Array<OneD, NekDouble>(order0*order1*order2, 0.0);

            NekDouble *mat = wsp.get();

            // put coeffs into matrix and reverse order so that p index is fastest recall q is fastest for tri's
            Vmath::Zero(order0*order1*order2, mat, 1);

            for(int i = 0, cnt=0; i < order0; ++i)
            {
                for(int j = 0; j < order1-i; ++j)
                {
                    for(int k = 0; k < order2-i-j; ++k, cnt++)
                    {
                        //                         mat[i+j*order1] = m_coeffs[cnt];
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
                        //                         outfile << mat[j*order0+i] <<" ";
                        outfile << mat[i + order0*(j + order1*k)] <<" ";
                    }
                    outfile << std::endl;
                }
            }
            outfile << "]" ;
        }


        NekDouble StdTetExp::v_PhysEvaluate(const Array<OneD, const NekDouble>& xi)
        {
            ASSERTL0(xi[0] + xi[1] + xi[2] + 1 <= NekConstants::kNekZeroTol,
                     "Coordinate outside bounds of tetrahedron.");
            ASSERTL0((xi[0] + 1 + NekConstants::kNekZeroTol >= 0) 
                        || (xi[1] + 1 + NekConstants::kNekZeroTol >= 0) 
                        || (xi[2] + 1 + NekConstants::kNekZeroTol >= 0),
                     "Coordinate outside bounds of tetrahedron.")
                              
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);

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

            ASSERTL0((eta[0] + NekConstants::kNekZeroTol >= -1)
                        || (eta[1] + NekConstants::kNekZeroTol >= -1)
                        || (eta[2] + NekConstants::kNekZeroTol >= -1),
                     "Eta Coordinate outside bounds of tetrahedron.");
            ASSERTL0((eta[0] - NekConstants::kNekZeroTol <= 1)
                        || (eta[1] - NekConstants::kNekZeroTol <= 1)
                        || (eta[2] - NekConstants::kNekZeroTol <= 1),
                     "Eta Coordinate outside bounds of tetrahedron.");
            
            return  StdExpansion3D::v_PhysEvaluate(eta);
        }




    }//end namespace
}//end namespace

/**
 * $Log: StdTetExp.cpp,v $
 * Revision 1.27  2009/12/15 18:09:02  cantwell
 * Split GeomFactors into 1D, 2D and 3D
 * Added generation of tangential basis into GeomFactors
 * Updated ADR2DManifold solver to use GeomFactors for tangents
 * Added <GEOMINFO> XML session section support in MeshGraph
 * Fixed const-correctness in VmathArray
 * Cleaned up LocalRegions code to generate GeomFactors
 * Removed GenSegExp
 * Temporary fix to SubStructuredGraph
 * Documentation for GlobalLinSys and GlobalMatrix classes
 *
 * Revision 1.26  2009/11/02 19:15:43  cantwell
 * Moved ContField1D to inherit from DisContField1D.
 * Moved ContField3D to inherit from DisContField3D.
 * Incorporated GenExpList1D functionality into ExpList1D.
 * Tidied up and added documentation to various classes.
 * Moved Namespace documentation and introductions to separate files along with
 * doxygen configuration.
 * Added option to use system ZLIB library instead of libboost_zlib on UNIX.
 * Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
 * Updated Linux compiling instructions.
 * Updated regDemo to use Helmholtz2D-g when built as debug.
 *
 * Revision 1.25  2009/04/27 21:32:45  sherwin
 * Updated WriteToField method
 *
 * Revision 1.24  2009/04/20 16:11:47  sherwin
 * Mods to handle output and optimise DG work
 *
 * Revision 1.23  2009/04/03 14:57:34  sherwin
 * Linear Advection matrices added, corrected unsigned int intialisation
 *
 * Revision 1.22  2009/01/01 02:40:25  ehan
 * Added GetFaceToElementMap().
 *
 * Revision 1.21  2008/12/18 14:11:35  pvos
 * NekConstants Update
 *
 * Revision 1.20  2008/11/23 00:33:46  sherwin
 * Added blas based IProductWRTBase and BwdTrans routines
 *
 * Revision 1.19  2008/09/17 13:46:06  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.18  2008/07/19 21:12:54  sherwin
 * Removed MapTo function and made orientation convention anticlockwise in UDG routines
 *
 * Revision 1.17  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.16  2008/06/16 22:46:43  ehan
 * Populated the function GetFaceToElementMap(..)
 *
 * Revision 1.15  2008/06/06 23:22:47  ehan
 * Added doxygen documentation
 *
 * Revision 1.14  2008/06/05 15:06:06  pvos
 * Added documentation
 *
 * Revision 1.13  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.12  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.11  2008/05/15 22:42:15  ehan
 * Added WriteToFile() function and its virtual function
 *
 * Revision 1.10  2008/05/07 16:04:57  pvos
 * Mapping + Manager updates
 *
 * Revision 1.9  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.8  2008/03/25 08:41:00  ehan
 * Added MapTo() function.
 *
 * Revision 1.7  2008/02/01 20:05:49  ehan
 * Added doxygen comments.
 *
 * Revision 1.6  2007/11/08 14:33:05  ehan
 * Fixed PhysDerivative and PhysTensorDerivative3D,  and improved L1 error up to 1e-15
 * Reimplimented WriteToFile and WriteCoeffsToFile function.
 *
 * Revision 1.5  2007/10/29 20:35:07  ehan
 * Fixed floating point approximation up to 1-e15 for PhysEvaluate.
 *
 * Revision 1.4  2007/10/15 20:40:06  ehan
 * Completed Basis, Backward, and Forward transformation
 *
 * Revision 1.3  2007/07/20 02:16:55  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.2  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.1  2006/05/04 18:58:33  kirby
 * *** empty log message ***
 *
 * Revision 1.7  2006/03/06 17:12:46  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.6  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 **/





