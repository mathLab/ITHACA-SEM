///////////////////////////////////////////////////////////////////////////////
//
// File StdPrismExp.cpp
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
// Description: Prismatic routines built upon StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdPrismExp.h>

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
                        for( int c = 0; c < Nc - a; ++c )
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }

        StdPrismExp::StdPrismExp() // Deafult construct of standard expansion directly called. 
        {
        }

        StdPrismExp::StdPrismExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc) 
            : StdExpansion3D(getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()), Ba, Bb, Bc)
            //  : StdExpansion3D(Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes(), Ba, Bb, Bc)
        {

            if(Ba.GetNumModes() >  Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher than order in 'c' direction");
            }    
        }

        StdPrismExp::StdPrismExp(const StdPrismExp &T)
            : StdExpansion3D(T)
        {
        }


        // Destructor
        StdPrismExp::~StdPrismExp()
        {   
        } 


        //////////////////////////////
        /// Integration Methods
        //////////////////////////////

        namespace 
        {
            void TripleTensorProduct(   const ConstArray<OneD, NekDouble>& fx, 
                                        const ConstArray<OneD, NekDouble>& gy, 
                                        const ConstArray<OneD, NekDouble>& hz, 
                                        const ConstArray<OneD, NekDouble>& inarray, 
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
                                        const ConstArray<OneD, NekDouble>& fxyz, 
                                        const ConstArray<OneD, NekDouble>& wx, 
                                        const ConstArray<OneD, NekDouble>& wy, 
                                        const ConstArray<OneD, NekDouble>& wz
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
  
        NekDouble StdPrismExp::Integral3D(const ConstArray<OneD, NekDouble>& inarray, 
                                        const ConstArray<OneD, NekDouble>& wx,
                                        const ConstArray<OneD, NekDouble>& wy, 
                                        const ConstArray<OneD, NekDouble>& wz)
        {
            return TripleInnerProduct( inarray, wx, wy, wz );

        }
        
        NekDouble StdPrismExp::Integral(const ConstArray<OneD, NekDouble>& inarray)
        {
            // Using implementation from page 146 of Spencer Sherwin's book
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            // Get the point distributions:
            // x is assumed to be Gauss-Lobatto-Legendre (includes -1 and +1)
            // y is assumed to be Gauss-Lobatto-Legendre (includes -1 and +1)
            ConstArray<OneD, NekDouble> z,wx,wy,wz;
            wx = ExpPointsProperties(0)->GetW();
            wy = ExpPointsProperties(1)->GetW();
            ExpPointsProperties(2)->GetZW(z,wz);

            Array<OneD, NekDouble> wz_hat = Array<OneD, NekDouble>(Qz, 0.0);

            
            // Convert wz into wz_hat, which includes the 1/2 scale factor.
            // Nothing else need be done if the point distribution is
            // Jacobi (1,0) since (1 - xi_z) is aready factored into the weights.
            // Note by coincidence, xi_y = eta_y, xi_z = eta_z (xi_z = z according to our notation)
            switch(m_base[2]->GetPointsType())
            {
            
            // Common case
            case LibUtilities::eGaussRadauMAlpha1Beta0: // (1,0) Jacobi Inner product 
                Vmath::Smul( Qz, 0.5, (NekDouble *)wz.get(), 1, wz_hat.get(), 1 );
                break;
            
            // Corner cases
            case LibUtilities::eGaussLobattoLegendre:   // Legendre inner product (Falls-through to next case)
            case LibUtilities::eGaussRadauMLegendre:    // (0,0) Jacobi Inner product 
                for(int k = 0; k < Qz; ++k)
                {   
                    wz_hat[k] = 0.5*(1.0 - z[k]) * wz[k];
                 }
                break;
            }


            return Integral3D( inarray, wx, wy, wz_hat );
        }
        
        void StdPrismExp::IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata(),inarray,outarray);
        }


        // Interior prism implementation based on Spen's book page 119. and 608.  
        void StdPrismExp::IProductWRTBase(  const ConstArray<OneD, NekDouble>& bx, 
                                            const ConstArray<OneD, NekDouble>& by, 
                                            const ConstArray<OneD, NekDouble>& bz, 
                                            const ConstArray<OneD, NekDouble>& inarray, 
                                            Array<OneD, NekDouble> & outarray )
        {
            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;


            // Create an index map from the hexahedron to the prsim.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p ; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }
            
            // Create an index map from the rectangle to the triangle. 
            Array<OneD, int> pr = Array<OneD, int>( (P+1)*(R+1), -1 );
                for( int p = 0, mode=0; p <= P; ++p ) {
                    for( int r = 0; r <= R - p ; ++r, ++mode ) {
                        int index = r + (R+1)*p;
                        pr[r + (R+1)*p] = mode;
                    }
                }

            // Compute innerproduct over each mode in the Prismatic domain
            for( int p = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p; ++r ) {

                        // Determine the index for specifying which mode to use in the basis                       
                        int mode_pqr   = pqr[r + (R+1)*(q + (Q+1)*p)];
                        int mode_pr    = pr[r + (R+1)*p];

                        // Compute tensor product of inarray with the 3 basis functions
                        Array<OneD, NekDouble> g_pqr = Array<OneD, NekDouble>( Qx*Qy*Qz, 0.0 );
                        for( int k = 0; k < Qz; ++k ) {
                            for( int j = 0; j < Qy; ++j ) {
                                for( int i = 0; i < Qx; ++i ) {
                                    int s = i + Qx*(j + Qy*k);
                                    g_pqr[s] += inarray[s] * 
                                            bx[i + Qx*p] * 
                                            by[j + Qy*q] * 
                                            bz[k + Qz*mode_pr];
                                }
                            }
                        }

                        outarray[mode_pqr] = Integral( g_pqr );
                    }
                }
            }
        }


        //-----------------------------
        // Differentiation Methods
        //-----------------------------
        
      void StdPrismExp::PhysDeriv( Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2)
        {
            PhysDeriv(this->m_phys, out_d0, out_d1, out_d2);
        }
               
      // PhysDerivative implementation based on Spen's book page 152.    
      void StdPrismExp::PhysDeriv(const ConstArray<OneD, NekDouble>& u_physical, 
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


            ConstArray<OneD, NekDouble> eta_x, eta_y, eta_z;
            eta_x = ExpPointsProperties(0)->GetZ();
            eta_y = ExpPointsProperties(1)->GetZ();
            eta_z = ExpPointsProperties(2)->GetZ();

            
            for(int k=0, n=0; k<Qz; ++k)
                for(int j=0; j<Qy; ++j){
                    for(int i=0; i<Qx; ++i, ++n){
                    {
                    
                    
                        out_dxi1[n] = 2.0/ (1.0 - eta_z[k])*out_dEta1[n];
                        out_dxi2[n] = out_dEta2[n];
                        out_dxi3[n] = (1.0 + eta_x[i]) / (1.0 - eta_z[k])*out_dEta1[n] + out_dEta3[n];
                                        
//                         out_dxi1[n] = out_dEta1[n];
//                         out_dxi2[n] = 2.0/ (1.0 - eta_z[k])*out_dEta2[n];                     
//                         out_dxi3[n] = (1.0 + eta_y[j]) / (1.0 - eta_z[k])*out_dEta2[n] + out_dEta3[n];

                                              
                        //cout << "eta_x["<<i<<"] = " <<  eta_x[i] << ",  eta_y["<<j<<"] = " << eta_y[j] << ", eta_z["<<k<<"] = " <<eta_z[k] << endl;
                        //cout << "out_dEta1["<<n<<"] = " << out_dEta1[n] << ",  out_dEta2["<<n<<"] = " << out_dEta2[n] << ", out_dEta3["<<n<<"] = " <<out_dEta3[n] << endl;
                        //cout << "out_dxi1["<<n<<"] = " << out_dxi1[n] << ",  out_dxi2["<<n<<"] = " << out_dxi2[n] << ", out_dxi3["<<n<<"] = " << out_dxi3[n] << endl;
                        
                    }
                }
            }
                        
        }

       void StdPrismExp::FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;


         // Create an index map from the hexahedron to the prsim.
            Array<OneD, int> mode_pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, m = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p ; ++r, ++m ) {
                       mode_pqr[r + (R+1)*(q + (Q+1)*p)] = m;
                    }
                }
            }
            
            Array<OneD, int> mode_pr = Array<OneD, int>( (P+1)*(R+1), -1 );
                for( int p = 0, m=0; p <= P; ++p ) {
                    for( int r = 0; r <= R - p ; ++r, ++m ) {
                        int index = r + (R+1)*p;
                        mode_pr[r + (R+1)*p] = m;
                    }
                }

            // Find the pqr matching the provided mode
            int mode_p=0, mode_q=0, mode_r=0;
            for( int p = 0, m = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q ; ++q ) {
                    for( int r = 0; r <= R - p; ++r, ++m ) {
                        if( m == mode ) {
                            mode_p = p;
                            mode_q = q;
                            mode_r = r;
                        }
                    }
                }
            }

            const ConstArray<OneD, NekDouble>& bx = m_base[0]->GetBdata();
            const ConstArray<OneD, NekDouble>& by = m_base[1]->GetBdata();
            const ConstArray<OneD, NekDouble>& bz = m_base[2]->GetBdata();

            int p = mode_p, q = mode_q, r = mode_r;
            
            // Determine the index for specifying which mode to use in the basis
            int sigma_p   = Qx*p;
            int sigma_q   = Qy*q;         
            int sigma_qr  = Qz*mode_pr[r + (R+1)*p];
            int sigma = Qx*Qy*Qz*mode;


            // Compute tensor product of inarray with the 3 basis functions
            for( int k = 0; k < Qz; ++k ) {
                for( int j = 0; j < Qy; ++j ) {
                    for( int i = 0; i < Qx; ++i ) {
                        int s = i + Qx*(j + Qy*(k + Qz*mode));
                        outarray[s] = 
                                bx[i + sigma_p] * 
                                by[j + sigma_q] * 
                                bz[k + sigma_qr];
                    }
                }
            }

        }


        ///////////////////////////////
        // Evaluation Methods
        ///////////////////////////////

        void StdPrismExp::BwdTrans(
            const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray)
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

            ConstArray<OneD, NekDouble> xBasis  = m_base[0]->GetBdata();
            ConstArray<OneD, NekDouble> yBasis  = m_base[1]->GetBdata();
            ConstArray<OneD, NekDouble> zBasis  = m_base[2]->GetBdata();

            
             // Create an index map from the hexahedron to the prsim.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p ; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }                
              
            Array<OneD, int> pr = Array<OneD, int>( (P+1)*(R+1), -1 );
                for( int p = 0, mode=0; p <= P; ++p ) {
                    for( int r = 0; r <= R - p ; ++r, ++mode ) {
                        int index = r + (R+1)*p;
                        pr[r + (R+1)*p] = mode;
                    }
            }
            
           // Sum-factorize the triple summation starting with the z-dimension
            for( int k = 0; k < Qz; ++k ) {

                // Create the matrix of coefficients summed over the z-modes
                Array<OneD, NekDouble> Ak((P+1)*(Q+1), 0.0);
                for( int p = 0; p <= P; ++p ) {
                    for( int q = 0; q <= Q; ++q ) {
                        for( int r = 0; r <= R - p; ++r ) {
                            int mode_pqr = pqr[r + (R+1)*(q + (Q+1)*p)];

                            int mode     = pr[r + (R+1)*p];
                            Ak[q + (Q+1)*p]   +=   inarray[mode_pqr]  *  zBasis[k + Qz*mode];                           
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

        void StdPrismExp::FwdTrans( const ConstArray<OneD, NekDouble>& inarray,  Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
            
         
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);
            
             //  cout << "matsys = \n" <<  *matsys << endl;

            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs, outarray);
            DNekVec out(m_ncoeffs, outarray, eWrapper);

            out = (*matsys)*in;
            
        }



        NekDouble StdPrismExp::PhysEvaluate(const ConstArray<OneD, NekDouble>& xi)
        {
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);

            if( fabs(xi[2]-1.0) < NekConstants::kEvaluateTol )    // NekConstants::kEvaluateTol = 1e-12
            {
                // Very top point of the prism
                eta[0] = -1.0;
                eta[1] = -1.0;
                eta[2] = xi[2];
            }
            else  //  Below the line-singularity -- Common case
            {
                 // Third basis function collapsed to "pr" direction instead of "qr" direction
                eta[2] = xi[2]; // eta_z = xi_z
                eta[1] = xi[1]; //eta_y = xi_y
                eta[0] = 2.0*(1.0 + xi[0])/(1.0 - xi[2]) - 1.0;
                
                // Third basis function collapsed to "qr" direction instead of "pr" direction
//                 eta[1] = 2.0*(1.0 + xi[1])/(1.0 - xi[2]) - 1.0;
//                 eta[0] = xi[0]; //eta_x = xi_x
                

            } 


            return  StdExpansion3D::PhysEvaluate3D(eta);  
        }
 
        void StdPrismExp::GetCoords( Array<OneD, NekDouble> & xi_x, Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z)
        {
            ConstArray<OneD, NekDouble> eta_x = ExpPointsProperties(0)->GetZ();
            ConstArray<OneD, NekDouble> eta_y = ExpPointsProperties(1)->GetZ();
            ConstArray<OneD, NekDouble> eta_z = ExpPointsProperties(2)->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates: eta --> xi
            for( int k = 0; k < Qz; ++k ) {
                for( int j = 0; j < Qy; ++j ) {
                    for( int i = 0; i < Qx; ++i ) {
                        int s = i + Qx*(j + Qy*k);

                        //NekDouble eta_x_bar = (1.0 + eta_x[i]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                        xi_z[s] = eta_z[k];                   
                        xi_y[s] = eta_y[j];
                        xi_x[s] = (1.0 + eta_x[i]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                        
                        // Third basis function collapsed to "qr" direction instead of "pr" direction
//                         xi_y[s] = (1.0 + eta_y[j]) * (1.0 - eta_z[k]) / 2.0  -  1.0;                    
//                         xi_x[s] = eta_x[i];

                    //    xi_x[s] = (1.0 + eta_x_bar) * (1.0 - eta_z[k]) / 2.0  -  1.0;                    
                       
                    }
                }
            }
        }



      

        void StdPrismExp::GenLapMatrix(double * outarray)
        {
            ASSERTL0(false, "Not implemented");
        }








//    StdMatrix StdPrismExp::s_elmtmats;
  
  }//end namespace
}//end namespace

/** 
 * $Log: StdPrismExp.cpp,v $
 * Revision 1.3  2008/01/03 10:54:43  ehan
 * Fixed standard prismatic domain.
 *
 * Revision 1.1  2006/05/04 18:58:32  kirby
 * *** empty log message ***
 *
 * Revision 1.5  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 **/ 



