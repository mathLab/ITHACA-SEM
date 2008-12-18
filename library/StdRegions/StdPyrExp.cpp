///////////////////////////////////////////////////////////////////////////////
//
// File StdPyrExp.cpp
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
// Description: pyramadic routines built upon StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdPyrExp.h>

namespace Nektar
{
    namespace StdRegions
    {
        namespace 
        {
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc ) 
            {
                int nCoef = 0;
                for( int a = 0; a < Na; ++a )
                {
                    for( int b = 0; b < Nb; ++b )
                    {
                        for( int c = 0; c < Nc - a - b; ++c )
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }
        
        StdPyrExp::StdPyrExp() // Deafult construct of standard expansion directly called. 
        {
        }
        
        StdPyrExp::StdPyrExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc) 
            : StdExpansion3D(getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()), Ba, Bb, Bc)
              //     : StdExpansion3D(Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes(), Ba, Bb, Bc)
        {

            if(Ba.GetNumModes() >  Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher than order in 'c' direction");
            }    
            if(Bb.GetNumModes() >  Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'b' direction is higher than order in 'c' direction");
            }    
        }

        StdPyrExp::StdPyrExp(const StdPyrExp &T)
            : StdExpansion3D(T)
        {
        }


        // Destructor
        StdPyrExp::~StdPyrExp()
        {   
        } 

        //////////////////////////////
        // Integration Methods
        //////////////////////////////

        namespace 
        {
            void TripleTensorProduct(   const Array<OneD, const NekDouble>& fx, 
                                        const Array<OneD, const NekDouble>& gy, 
                                        const Array<OneD, const NekDouble>& hz, 
                                        const Array<OneD, const NekDouble>& inarray, 
                                        Array<OneD, NekDouble> & outarray )
            {
            
                // Using matrix operation, not sum-factorization.
                // Regarding the 3D array, inarray[k][j][i], x is changing the fastest and z the slowest.
                // Thus, the first x-vector of points refers to the first row of the first stack. The first y-vector
                // refers to the first column of the first stack. The first z-vector refers to the vector of stacks
                // intersecting the first row and first column. So in C++, i refers to column, j to row, and k to stack.
                // Contrasting this with the usual C++ matrix convention, note that i does not refer to a C++ row, nor j to C++ column.

                int nx = fx.num_elements();
                int ny = gy.num_elements();
                int nz = hz.num_elements();


                // Multiply by integration constants...
                // Hadamard multiplication refers to elementwise multiplication of two vectors.

                // Hadamard each row with the first vector (x-vector); the index i is changing the fastest.
                for (int jk = 0; jk < ny*nz; ++jk)  // For each j and k, iterate over each row in all of the stacks at once
                {
                    Vmath::Vmul(
                                nx,                         // Size of first weight vector
                                &inarray[0] + jk*nx, 1,     // Offset and stride of each row-vector (x is changing fastest)
                                fx.get(), 1,                // First weight vector (with stride of 1)
                                &outarray[0] + jk*nx, 1     // Output has same offset and stride as input
                                );
                }

                // Hadamard each column with the second vector (y-vector)
                for (int k = 0; k < nz; ++k)                    // For each stack in the 3D-array, do the following...
                {
                    for (int i = 0; i < nx; ++i)                // Iterate over each column in the current stack
                    {
                        Vmath::Vmul(
                                    ny,                                 // Size of second weight vector
                                    &outarray[0] + i + nx*ny*k, nx,     // Offset and stride of each column-vector
                                    gy.get(), 1,                        // second weight vector (with stride of 1)
                                    &outarray[0] + i + nx*ny*k, nx      // Output has same offset and stride as input
                                    );
                    }
                }

                // Hadamard each stack-vector with the third vector (z-vector)
                for (int ij = 0; ij < nx*ny; ++ij)              // Iterate over each element in the topmost stack
                {
                    Vmath::Vmul(
                                nz,                                     // Size of third weight vector
                                &outarray[0] + ij, nx*ny,               // Offset and stride of each stack-vector
                                hz.get(), 1,                            // Third weight vector (with stride of 1)
                                &outarray[0] + ij, nx*ny                // Output has same offset and stride as input
                                );
                }

            }
            

            // Inner-Product with respect to the weights: i.e., this is the triple sum of the product 
            // of the four inputs over the Hexahedron
            // x-dimension is the row, it is the index that changes the fastest
            // y-dimension is the column
            // z-dimension is the stack, it is the index that changes the slowest
            NekDouble TripleInnerProduct( 
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
                    cerr << "TripleInnerProduct expected " << fxyz.num_elements() << 
                        " quadrature points from the discretized input function but got " << 
                        Qx*Qy*Qz << " instead." << endl;
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
        }


        NekDouble StdPyrExp::Integral3D(const Array<OneD, const NekDouble>& inarray, 
                                        const Array<OneD, const NekDouble>& wx,
                                        const Array<OneD, const NekDouble>& wy, 
                                        const Array<OneD, const NekDouble>& wz)
        {
            return TripleInnerProduct( inarray, wx, wy, wz );

        }
        
        /** \brief Integrate the physical point list \a inarray over pyramidic region and return the value
                
            Inputs:\n

            - \a inarray: definition of function to be returned at quadrature point of expansion.

            Outputs:\n

            - returns \f$\int^1_{-1}\int^1_{-1}\int^1_{-1} u(\bar \eta_1, \eta_2, \eta_3) J[i,j,k] d \bar \eta_1 d \eta_2 d \eta_3\f$ \n
            \f$= \sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1} \sum_{k=0}^{Q_3 - 1} u(\bar \eta_{1i}^{0,0}, \eta_{2j}^{0,0},\eta_{3k}^{2,0})w_{i}^{0,0} w_{j}^{0,0} \hat w_{k}^{2,0}    \f$ \n
            where \f$inarray[i,j, k] = u(\bar \eta_{1i},\eta_{2j}, \eta_{3k}) \f$, \n
            \f$\hat w_{k}^{2,0} = \frac {w^{2,0}} {2} \f$ \n
            and \f$ J[i,j,k] \f$ is the  Jacobian evaluated at the quadrature point.

        */
        NekDouble StdPyrExp::Integral(const Array<OneD, const NekDouble>& inarray)
        {
            // Using implementation from page 146 - 147 of Spencer Sherwin's book
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            // Get the point distributions:
            // x is assumed to be Gauss-Lobatto-Legendre (includes -1 and +1)
            // y is assumed to be Gauss-Lobatto-Legendre (includes -1 and +1)
            Array<OneD, const NekDouble> z, wx, wy, wz;
            wx = m_base[0]->GetW();
            wy = m_base[1]->GetW();
            m_base[2]->GetZW(z,wz);

            Array<OneD, NekDouble> wz_hat = Array<OneD, NekDouble>(Qz, 0.0);

            
            // Convert wz into wz_hat, which includes the 1/4 scale factor.
            // Nothing else need be done if the point distribution is
            // Jacobi (1,0) since (1 - xi_z) is aready factored into the weights.
            // Note by coincidence, xi_z = eta_z (xi_z = z according to our notation)
            switch(m_base[2]->GetPointsType())
            {
            
                // Common case
            case LibUtilities::eGaussRadauMAlpha2Beta0: // (2,0) Jacobi Inner product 
                Vmath::Smul( Qz, 0.25, (NekDouble *)wz.get(), 1, wz_hat.get(), 1 );
                break;
            
                // Corner cases
            case LibUtilities::eGaussLobattoLegendre:   // Legendre inner product (Falls-through to next case)
            case LibUtilities::eGaussRadauMLegendre:    // (0,0) Jacobi Inner product 
                for(int k = 0; k < Qz; ++k)
                {   
                    wz_hat[k] = 0.25*(1.0 - z[k])*(1.0 - z[k]) * wz[k];
                }
                break;
            }


            return Integral3D( inarray, wx, wy, wz_hat );
        }

        void StdPyrExp::IProductWRTBase(const Array<OneD, const NekDouble>& bx, 
                                        const Array<OneD, const NekDouble>& by, 
                                        const Array<OneD, const NekDouble>& bz, 
                                        const Array<OneD, const NekDouble>& inarray, 
                                        Array<OneD, NekDouble> & outarray)
        {
            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;


            // Create an index map from the hexahedron to the pyramid.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p - q ; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }

            // Compute innerproduct over each mode in the Pyramid domain
            for( int p = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p - q; ++r ) {

                        // Determine the index for specifying which mode to use in the basis                       
                        int mode_pqr   = pqr[r + (R+1)*(q + (Q+1)*p)];

                        // Compute tensor product of inarray with the 3 basis functions
                        Array<OneD, NekDouble> g_pqr = Array<OneD, NekDouble>( Qx*Qy*Qz, 0.0 );
                        for( int k = 0; k < Qz; ++k ) {
                            for( int j = 0; j < Qy; ++j ) {
                                for( int i = 0; i < Qx; ++i ) {
                                    int s = i + Qx*(j + Qy*k);
                                    g_pqr[s] += inarray[s] * 
                                        bx[i + Qx*p] * 
                                        by[j + Qy*q] * 
                                        bz[k + Qz*mode_pqr];
                                }
                            }
                        }

                        outarray[mode_pqr] = Integral( g_pqr );
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

            \f$\begin{matrix} \frac {\partial u^{\delta}} {\partial \xi_1} (\xi_{1i}, \xi_{2j}, \xi_{3k})  \\ 
            \frac {\partial u^{\delta}} {\partial \xi_2} (\xi_{1i}, \xi_{2j}, \xi_{3k})\\ 
            \frac {\partial u^{\delta}} {\partial \xi_3} (\xi_{1i}, \xi_{2j}, \xi_{3k})  \end{matrix} 
            = \begin{matrix} \sum_p^P u_{pjk} \frac {\partial h_p (\xi_1)}{\partial \xi_1} |_{\xi_{1i}},  \\ 
            \sum_p^Q u_{iqk} \frac {\partial h_q (\xi_2)}{\partial \xi_2} |_{\xi_{2j}},  \\
            \sum_p^R u_{ijr} \frac {\partial h_r (\xi_3)}{\partial \xi_3} |_{\xi_{3k}}
            \end{matrix}\f$

        **/
        // PhysDerivative implementation based on Spen's book page 152.    
        void StdPyrExp::PhysDeriv(const Array<OneD, const NekDouble>& u_physical, 
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
                    
                            out_dxi1[n] = 2.0/ (1.0 - eta_z[k])*out_dEta1[n];
                            out_dxi2[n] = 2.0/ (1.0 - eta_z[k])*out_dEta2[n];
                            out_dxi3[n] = (1.0 + eta_x[i]) / (1.0 - eta_z[k])*out_dEta1[n] + (1.0 + eta_y[j])/(1.0 - eta_z[k])*out_dEta2[n] + out_dEta3[n];
                                                                  
                            //cout << "eta_x["<<i<<"] = " <<  eta_x[i] << ",  eta_y["<<j<<"] = " << eta_y[j] << ", eta_z["<<k<<"] = " <<eta_z[k] << endl;
                            //cout << "out_dEta1["<<n<<"] = " << out_dEta1[n] << ",  out_dEta2["<<n<<"] = " << out_dEta2[n] << ", out_dEta3["<<n<<"] = " <<out_dEta3[n] << endl;
                            //cout << "out_dxi1["<<n<<"] = " << out_dxi1[n] << ",  out_dxi2["<<n<<"] = " << out_dxi2[n] << ", out_dxi3["<<n<<"] = " << out_dxi3[n] << endl;
                        
                        }
                    }
                }
                        
        }
        
        void StdPyrExp::FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;


            // Create an index map from the hexahedron to the prymaid.
            Array<OneD, int> mode_pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, m = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p - q ; ++r, ++m ) {
                        mode_pqr[r + (R+1)*(q + (Q+1)*p)] = m;
                    }
                }
            }

            // Find the pqr matching the provided mode
            int mode_p=0, mode_q=0, mode_r=0;
            for( int p = 0, m = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q ; ++q ) {
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
            int sigma_q   = Qy*q;  
            int sigma_pqr = Qz*mode_pqr[r + (R+1)*(q + (Q+1)*p)];       


            // Compute tensor product of inarray with the 3 basis functions
            for( int k = 0; k < Qz; ++k ) {
                for( int j = 0; j < Qy; ++j ) {
                    for( int i = 0; i < Qx; ++i ) {
                        int s = i + Qx*(j + Qy*(k + Qz*mode));
                        outarray[s] = 
                            bx[i + sigma_p] * 
                            by[j + sigma_q] * 
                            bz[k + sigma_pqr];
                    }
                }
            }

        }
        
        ///////////////////////////////
        /// Evaluation Methods
        ///////////////////////////////

	/** 
            \brief Backward transformation is evaluated at the quadrature points 
		
	    \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
	    
            Backward transformation is three dimensional tensorial expansion
		
	    \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
	    \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c (\xi_{3k}) \rbrace}
	    \rbrace}. \f$
	       
            And sumfactorizing step of the form is as:\\
            \f$ f_{pqr} (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c (\xi_{3k}),\\ 
            g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y} \psi_{p}^a (\xi_{2j}) f_{pqr} (\xi_{3k}),\\
            u(\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a (\xi_{1i}) g_{p} (\xi_{2j}, \xi_{3k}).
            \f$	
        **/
        void StdPyrExp::BwdTrans(const Array<OneD, const NekDouble>& inarray, Array<OneD, NekDouble> &outarray)
        {

            ASSERTL1( (m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)  ||
                      (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                      "Basis[1] is not a general tensor type");

            ASSERTL1( (m_base[2]->GetBasisType() != LibUtilities::eOrtho_C) ||
                      (m_base[2]->GetBasisType() != LibUtilities::eModified_C),
                      "Basis[2] is not a general tensor type");

            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;

            Array<OneD, const NekDouble> xBasis  = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> yBasis  = m_base[1]->GetBdata();
            Array<OneD, const NekDouble> zBasis  = m_base[2]->GetBdata();

            
            // Create an index map from the hexahedron to the pyramid.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p - q ; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }                
            cout << zBasis.num_elements() / Qz << "  " << m_ncoeffs << endl;
            // Sum-factorize the triple summation starting with the z-dimension
            for( int k = 0; k < Qz; ++k ) {

                // Create the matrix of coefficients summed over the z-modes
                Array<OneD, NekDouble> Ak((P+1)*(Q+1), 0.0);
                for( int p = 0; p <= P; ++p ) {
                    for( int q = 0; q <= Q; ++q ) {
                        for( int r = 0; r <= R - p - q; ++r ) {
                            int mode = pqr[r + (R+1)*(q + (Q+1)*p)];
                            cout << p << "   " << q << "   " << r << endl;

                            Ak[q + (Q+1)*p]   +=   inarray[mode]  *  zBasis[k + Qz*mode];                           
                        }
                    }
                }

                // Factorize the y-dimension
                for( int j = 0; j < Qy; ++j ) {

                    // Create the vector of coefficients summed over the y and z-modes
                    Array<OneD, NekDouble> bjk(P+1, 0.0);
                    for( int p = 0; p <= P; ++p ) {
                        for( int q = 0; q <= Q; ++q ) {
                            int mode = q;
                            bjk[p]   +=   Ak[q + (Q+1)*p]  *  yBasis[j + Qy*mode];
                        }
                    }

                    // Factorize the x-dimension
                    for( int i = 0; i < Qx; ++i ) {

                        NekDouble cijk = 0.0;
                        for( int p = 0; p <= P; ++p ) {
                            int mode = p;
                            cijk   +=   bjk[p]  *  xBasis[i + Qx*mode];
                        }
                        outarray[i + Qx*(j + Qy*k)] = cijk;
                    }
                }
            }
        }


	/** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a (this)->m_coeffs  
            
            Inputs:\n
            
            - \a inarray: array of physical quadrature points to be transformed
            
            Outputs:\n
            
            - (this)->_coeffs: updated array of expansion coefficients. 
            
        */    
        void StdPyrExp::FwdTrans( const Array<OneD, const NekDouble>& inarray,  Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetExpansionType(),*this);
            
         
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);
            
            //  cout << "matsys = \n" <<  *matsys << endl;

            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs, outarray);
            DNekVec out(m_ncoeffs, outarray, eWrapper);

            out = (*matsys)*in;
            
        }
        
        NekDouble StdPyrExp::PhysEvaluate(const Array<OneD, const NekDouble>& xi)
        {
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);

            if( fabs(xi[2]-1.0) < NekConstants::kNekZeroTol )    // NekConstants::kEvaluateTol = 1e-12
            {
                // Very top point of the pyramid
                eta[0] = -1.0;
                eta[1] = -1.0;
                eta[2] = xi[2];
            }
            else  //  Below the line-singularity -- Common case
            {
            
                eta[2] = xi[2]; // eta_z = xi_z
                eta[1] = 2.0*(1.0 + xi[1])/(1.0 - xi[2]) - 1.0; 
                eta[0] = 2.0*(1.0 + xi[0])/(1.0 - xi[2]) - 1.0;
           
            } 
            
            // cout << "xi_x = " << xi[0] << ", xi_y = " << xi[1] << ", xi_z = " << xi[2] << endl;
            // cout << "eta_x = " << eta[0] << ", eta_y = " << eta[1] << ", eta_z = " << eta[2] << endl;
            
            return  StdExpansion3D::PhysEvaluate(eta);  
        }
        
        void StdPyrExp::GetCoords( Array<OneD, NekDouble> & xi_x, Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z)
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

                        xi_z[s] = eta_z[k];                   
                        xi_y[s] = (1.0 + eta_y[j]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                        xi_x[s] = (1.0 + eta_x[i]) * (1.0 - eta_z[k]) / 2.0  -  1.0;


                        //NekDouble eta_x_bar = (1.0 + eta_x[i]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                        //    xi_x[s] = (1.0 + eta_x_bar) * (1.0 - eta_z[k]) / 2.0  -  1.0;                    
                       
                    }
                }
            }
        }

        void StdPyrExp::GetFaceToElementMap(const int fid, const FaceOrientation faceOrient,
                                            Array<OneD, unsigned int> &maparray,
                                            Array<OneD, int>& signarray)
        {
            //TODO implement 

        }

        void StdPyrExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar)
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
                    outfile << "Variables = z1,  z2,  z3,  Coeffs \n" << std::endl;      
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
        void StdPyrExp::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble> wsp  = Array<OneD, NekDouble>(order0*order1*order2, 0.0);

            NekDouble *mat = wsp.get(); 

            // put coeffs into matrix and reverse order so that r index is fastest for Prism 
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
               
    }//end namespace
}//end namespace

/** 
 * $Log: StdPyrExp.cpp,v $
 * Revision 1.16  2008/11/17 09:02:30  ehan
 * Added necessary mapping routines
 *
 * Revision 1.15  2008/09/23 18:19:26  pvos
 * Updates for working ProjectContField3D demo
 *
 * Revision 1.14  2008/09/17 13:46:06  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.13  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.12  2008/06/16 22:46:26  ehan
 * Populated the function GetFaceToElementMap(..)
 *
 * Revision 1.11  2008/06/06 23:22:24  ehan
 * Added doxygen documentation
 *
 * Revision 1.10  2008/06/05 15:06:06  pvos
 * Added documentation
 *
 * Revision 1.9  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.8  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.7  2008/05/15 22:41:46  ehan
 * Added WriteToFile() function and its virtual function
 *
 * Revision 1.6  2008/05/07 16:04:57  pvos
 * Mapping + Manager updates
 *
 * Revision 1.5  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.4  2008/02/01 20:05:35  ehan
 * Added doxygen comments
 *
 * Revision 1.3  2008/01/03 12:33:00  ehan
 * Fixed errors from StdMatrix to StdMatrixKey.
 *
 * Revision 1.2  2008/01/03 10:34:51  ehan
 * Added basis, derivative, physEval, interpolation, integration, backward transform, forward transform functions.
 *
 * Revision 1.1  2006/05/04 18:58:32  kirby
 * *** empty log message ***
 *
 * Revision 1.5  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 **/ 


  
