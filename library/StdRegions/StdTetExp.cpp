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
        StdTetExp::StdTetExp()
        {
        }
        

        StdTetExp::StdTetExp(const LibUtilities::BasisKey &Ba,
                             const LibUtilities::BasisKey &Bb,
                             const LibUtilities::BasisKey &Bc):
            StdExpansion(StdTetData::getNumberOfCoefficients(
                             Ba.GetNumModes(),
                             Bb.GetNumModes(), 
                             Bc.GetNumModes()),
                         3, Ba, Bb, Bc),
            StdExpansion3D(StdTetData::getNumberOfCoefficients(
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

        //-------------------------------
        // Integration Methods
        //-------------------------------
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
            return TripleInnerProduct(inarray, wx, wy, wz);
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
        NekDouble StdTetExp::v_Integral(
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
            // Nothing else need be done if the point distribution is Jacobi
            // (1,0) since (1-eta_y) is aready factored into the weights.
            switch(m_base[1]->GetPointsType())
            {
                // Legendre inner product (Falls-through to next case)
                case LibUtilities::eGaussLobattoLegendre:
                // (0,0) Jacobi Inner product
                case LibUtilities::eGaussRadauMLegendre:
                for(int j = 0; j < Qy; ++j)
                {
                    wy_hat[j] = 0.5*(1.0 - y[j]) * wy[j];
                }
                break;
                
                // (1,0) Jacobi Inner product
                case LibUtilities::eGaussRadauMAlpha1Beta0: 
                Vmath::Smul( Qy, 0.5, (NekDouble *)wy.get(), 1, wy_hat.get(), 1 );
                break;
                
                default:
                    ASSERTL0(false, "Unsupported quadrature points type.");
                    break;
            }

            // Convert wz into wz_hat, which includes the 1/4 scale factor.
            // Nothing else need be done if the point distribution is Jacobi
            // (2,0) since (1-eta_z)^2 is aready factored into the weights.
            // Note by coincidence, xi_z = eta_z (xi_z = z according to our
            // notation)
            switch(m_base[2]->GetPointsType())
            {
                // Legendre inner product (Falls-through to next case)
                case LibUtilities::eGaussLobattoLegendre:
                // (0,0) Jacobi Inner product
                case LibUtilities::eGaussRadauMLegendre:
                    for(int k = 0; k < Qz; ++k)
                    {
                        wz_hat[k] = 0.25*(1.0 - z[k])*(1.0 - z[k]) * wz[k];
                    }
                    break;
                // (2,0) Jacobi Inner product
                case LibUtilities::eGaussRadauMAlpha2Beta0: 
                    Vmath::Smul(Qz, 0.25, (NekDouble *)wz.get(), 1, 
                                wz_hat.get(), 1 );
                    break;
                
                default:
                    ASSERTL0(false, "Unsupported quadrature points type.");
                    break;
            }

            return Integral3D(inarray, wx, wy_hat, wz_hat);
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

            bool Do_2 = (out_dxi2.num_elements() > 0)? true:false;
            bool Do_1 = (out_dxi1.num_elements() > 0)? true:false;
            
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
            
            if (out_dxi0.num_elements() > 0)
            {
                // out_dxi1 = 4.0/((1-eta_1)(1-eta_2)) Out_dEta0
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
                    
                    // calculate out_dxi1 =
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

            Array<OneD, NekDouble> wsp(nquad2*order0*order1*(order1+1)/2+
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
        void StdTetExp::BwdTrans_SumFacKernel(
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
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();

            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble > tmp  = wsp;
            Array<OneD, NekDouble > tmp1 = tmp + nquad2*order0*order1*(order1+1)/2;

            //Array<OneD, NekDouble > tmp(nquad2*order0*(order1+1)/2);
            //Array<OneD, NekDouble > tmp1(nquad2*nquad1*order0);

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
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetExpansionType(),*this);
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
                  Array<OneD,       NekDouble>& outarray)
        {
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();

            Array<OneD, NekDouble> tmp (nquad0*nquad1*nquad2);
            Array<OneD, NekDouble> wsp (nquad1*nquad2*order0 +
                                        nquad2*order0*(order1+1)/2);

            MultiplyByQuadratureMetric(inarray, tmp);

            StdTetExp::IProductWRTBase_SumFacKernel(
                    m_base[0]->GetBdata(),
                    m_base[1]->GetBdata(),
                    m_base[2]->GetBdata(),
                    tmp, outarray, wsp, true, true, true);
        }


        void StdTetExp::IProductWRTBase_SumFacKernel(
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

            StdMatrixKey      iprodmatkey(mtype,DetExpansionType(),*this);
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
                + nquad1*nquad2*nmodes0 + nquad2*nmodes0*(nmodes1+1)/2;

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

        NekDouble StdTetExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble>& xi)
        {
            return v_PhysEvaluate(xi, m_phys);
        }

        NekDouble StdTetExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble>& xi,
            const Array<OneD, const NekDouble>& physvals)
        {
            // Validation checks
            ASSERTL0(xi[0] + xi[1] + xi[2] <= -1 + NekConstants::kNekZeroTol,
                     "Coordinate outside bounds of tetrahedron.");
            ASSERTL0(xi[0] >= -1 && xi[1] >= -1 && xi[2] >= -1,
                     "Coordinate outside bounds of tetrahedron.");

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

            ASSERTL0((eta[0] + NekConstants::kNekZeroTol >= -1) ||
                     (eta[1] + NekConstants::kNekZeroTol >= -1) ||
                     (eta[2] + NekConstants::kNekZeroTol >= -1),
                     "Eta Coordinate outside bounds of tetrahedron.");
            ASSERTL0((eta[0] - NekConstants::kNekZeroTol <= 1) ||
                     (eta[1] - NekConstants::kNekZeroTol <= 1) ||
                     (eta[2] - NekConstants::kNekZeroTol <= 1),
                     "Eta Coordinate outside bounds of tetrahedron.");

            return StdExpansion3D::v_PhysEvaluate(eta, physvals);
        }
        
        void StdTetExp::v_FillMode(
            const int                     mode, 
                  Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp(m_ncoeffs,0.0);
            tmp[mode] = 1.0;
            StdTetExp::v_BwdTrans(tmp, outarray);
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

        ExpansionType StdTetExp::v_DetExpansionType() const
        {
            return DetExpansionType();
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

            int p_hat, k;
            // All modes in the first layer are boundary modes
            int tot = P*(P+1)/2 + (Q-P)*P;
            // Loop over each plane in the stack
            for (int i = 1; i < R - 1; ++i)
            {
                p_hat = min(P, R-i);
                k = min(Q-P, max(0, Q-i-1));
                // First two columns and bottom row are boundary modes
                tot += (p_hat + k) + (p_hat + k - 1) + p_hat - 2;
            }

            // Add on top vertex mode
            return tot + 1;
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
            int nmodes = StdRegions::StdTetData::getNumberOfCoefficients(
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
            int nummodes = GetBasis(0)->GetNumModes(); 
            //temporary solution, need to add conditions based on face id
            //also need to add check of the points type
            switch (k)
            {
                case 0:
                {
                    const LibUtilities::PointsKey pkey(nummodes+1,LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eModified_A,nummodes,pkey);
                }
                break;
                case 1:
                {
                    const LibUtilities::PointsKey pkey(nummodes,LibUtilities::eGaussRadauMAlpha1Beta0);
                    //const LibUtilities::PointsKey pkey(nummodes+1,LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eModified_B,nummodes,pkey);
                }
                break;
            }

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

        void StdTetExp::v_WriteToFile(
            std::ofstream &outfile,
            OutputFormat   format, 
            const bool     dumpVar, 
            std::string    var)
        {
            if (format == eTecplot)
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
                outfile << "Zone, I=" << Qx <<", J=" << Qy 
                        <<", K=" << Qz <<", F=Point" << std::endl;

                for(int k = 0; k < Qz; ++k)
                {
                    for(int j = 0; j < Qy; ++j)
                    {
                        for(int i = 0; i < Qx; ++i)
                        {
                            outfile << (eta_x[i] + 1.0) * (1.0 - eta_y[j]) * (1.0 - eta_z[k]) / 4  -  1.0 <<  " " << eta_z[k] << " " << m_phys[i + Qx*(j + Qy*k)] << std::endl;
                        }
                    }
                }
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented "
                         "for requested type of output");
            }
        }

        void StdTetExp::v_WriteCoeffsToFile(std::ofstream &outfile)
        {
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble> wsp(order0*order1*order2, 0.0);
            NekDouble *mat = wsp.get();

            // put coeffs into matrix and reverse order so that p index is
            // fastest recall q is fastest for tri's
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
        
        /**
         * Maps Expansion2D modes of a 2D face to the corresponding expansion
         * modes.
         */
        void StdTetExp::v_GetFaceToElementMap(
            const int                  fid, 
            const Orientation      faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray,
            int                        nummodesA,
            int                        nummodesB)
        {
            int P, Q, i, j, k, idx;

            ASSERTL1(v_IsBoundaryInteriorExpansion(),
                     "Method only implemented for Modified_A BasisType (x "
                     "direction), Modified_B BasisType (y direction), and "
                     "Modified_C BasisType(z direction)");

            int nFaceCoeffs = 0;
            
            if (nummodesA == -1)
            {
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
                }
            }
            
            P = nummodesA;
            Q = nummodesB;

            nFaceCoeffs = Q + ((P-1)*(1 + 2*(Q-1) - (P-1)))/2;
            
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
                            if (j == 0 && k == 0) {
                                maparray[idx++] = GetMode(0,0,1);
                            }
                            if (j == 0 && k == Q-2) {
                                for (int r = 0; r < Q-1; ++r) {
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
        }
        
        int StdTetExp::v_GetVertexMap(const int localVertexId)
        {
            ASSERTL0((GetEdgeBasisType(localVertexId)==LibUtilities::eModified_A)||
                     (GetEdgeBasisType(localVertexId)==LibUtilities::eModified_B)||
                     (GetEdgeBasisType(localVertexId)==LibUtilities::eModified_C),
                     "Mapping not defined for this type of basis");

            int localDOF;
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

            if(maparray.num_elements() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }
            else
            {
            	fill( maparray.get(), maparray.get() + nEdgeIntCoeffs, 0);
            }

            if(signarray.num_elements() != nEdgeIntCoeffs)
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

            if(maparray.num_elements() != nFaceIntCoeffs)
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

            if(outarray.num_elements() != nIntCoeffs)
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
			return StdExpansion::CreateGeneralMatrix(mkey);
        }

        DNekMatSharedPtr StdTetExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
			return StdExpansion::CreateGeneralMatrix(mkey);
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

        void StdTetExp::MultiplyByQuadratureMetric(
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
                            Blas::Dscal(nquad0,0.5*w1[i], &outarray[0]+i*nquad0+
                                        j*nquad0*nquad1,1);
                        }
                    }
                    break;
                
                default:
                    ASSERTL0(false, "Unsupported quadrature points type.");
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

                default:
                    ASSERTL0(false, "Unsupported quadrature points type.");
                    break;
            }
        }
    }//end namespace
}//end namespace
