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
        namespace
        {
            // Adds up the number of cells in a truncated Nc by Nc by Nc pyramid, 
            // where the longest Na rows and longest Nb columns are kept.
            // Example: (Na, Nb, Nc) = (3, 4, 5); The number of coefficients is the 
            // sum of the elements of the following matrix:
            //     |5  4  3  2  0|
            //     |4  3  2  0   |
            //     |3  2  0      |
            //     |0  0         |
            //     |0            |
            // Sum = 28 = number of tet coefficients
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc ) 
            {
                int nCoef = 0;
                for( int a = 0; a < Na; ++a )
                {
                    for( int b = 0; b < Nb - a; ++b )
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


        StdTetExp::StdTetExp() // default constructor of StdExpansion is directly called. 
        {
        } //default constructor

        StdTetExp::StdTetExp( const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc ) :
            StdExpansion3D( getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()), Ba, Bb, Bc)
        {

            if(Ba.GetNumModes() >  Bb.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher than order in 'b' direction");
            }
            if(Ba.GetNumModes() >  Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher than order in 'c' direction");
            }
            if(Bb.GetNumModes() >  Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'b' direction is higher than order in 'c' direction");
            }
        }



        StdTetExp::StdTetExp(const StdTetExp &T):
        StdExpansion3D(T)
        {
        }

        StdTetExp::~StdTetExp()
        {
        }

        //////////////////////////////
        // Integration Methods
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
            

            // Inner-Product with respect to the weights: i.e., this is the triple sum of the product of the four inputs over the Hexahedron
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
                    cerr << "TripleTetrahedralInnerProduct expected " << fxyz.num_elements() << 
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


        NekDouble StdTetExp::Integral3D(const ConstArray<OneD, NekDouble>& inarray, const ConstArray<OneD, NekDouble>& wx,
                                        const ConstArray<OneD, NekDouble>& wy, const ConstArray<OneD, NekDouble>& wz)
        {
            return TripleInnerProduct( inarray, wx, wy, wz );

         }


        NekDouble StdTetExp::Integral(const ConstArray<OneD, NekDouble>& inarray)
        {
            // Using implementation from page 145 of Spencer Sherwin's book
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            // Get the point distributions:
            // x is assumed to be Gauss-Lobatto-Legendre (includes -1 and +1)
            ConstArray<OneD, NekDouble> y,z,wx,wy,wz;
            wx = ExpPointsProperties(0)->GetW();
            ExpPointsProperties(1)->GetZW(y,wy);
            ExpPointsProperties(2)->GetZW(z,wz);

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
                    NekDouble eta = 2.0*(1.0 + y[j])/(1.0 - z[j]) - 1.0;
                    wy_hat[j] = 0.5*(1.0 - eta) * wy[j];
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

        void StdTetExp::IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata(),inarray,outarray);
        }

        void StdTetExp::IProductWRTBase(    const ConstArray<OneD, NekDouble>& bx, 
                                            const ConstArray<OneD, NekDouble>& by, 
                                            const ConstArray<OneD, NekDouble>& bz, 
                                            const ConstArray<OneD, NekDouble>& inarray, 
                                            Array<OneD, NekDouble> & outarray)
        {

            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;

            // Build an index map from the rectangle to the triangle
            Array<OneD, int> pq  = Array<OneD, int>( (P+1)*(Q+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q - p; ++q, ++mode ) {
                    pq[q + (Q+1)*p] = mode;
                }
            }

            // Create an index map from the hexahedron to the tetrahedron.
            // The actual index is too difficult to compute explicitly.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q - p; ++q ) {
                    for( int r = 0; r <= R - p - q; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }


            // Compute innerproduct over each mode in the Tetrahedral domain
            for( int p = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q - p; ++q ) {
                    for( int r = 0; r <= R - p - q; ++r ) {

                        // Determine the index for specifying which mode to use in the basis
                        int mode_p      = p;
                        int mode_pq     = pq[q + (Q+1)*p];
                        int mode_pqr    = pqr[r + (R+1)*(q + (Q+1)*p)];


                        // Compute tensor product of inarray with the 3 basis functions
                        Array<OneD, NekDouble> g_pqr = Array<OneD, NekDouble>( Qx*Qy*Qz, 0.0 );
                        for( int k = 0; k < Qz; ++k ) {
                            for( int j = 0; j < Qy; ++j ) {
                                for( int i = 0; i < Qx; ++i ) {
                                    int s = i + Qx*(j + Qy*k);
                                    g_pqr[s] += inarray[s] * 
                                            bx[i + Qx*mode_p] * 
                                            by[j + Qy*mode_pq] * 
                                            bz[k + Qz*mode_pqr];
                                }
                            }
                        }

                        outarray[mode_pqr] = Integral( g_pqr );
                    }
                }
            }


        }

       void StdTetExp::FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
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
            Array<OneD, int> mode_pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
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

            const ConstArray<OneD, NekDouble>& bx = m_base[0]->GetBdata();
            const ConstArray<OneD, NekDouble>& by = m_base[1]->GetBdata();
            const ConstArray<OneD, NekDouble>& bz = m_base[2]->GetBdata();



            int p = mode_p, q = mode_q, r = mode_r;

            // Determine the index for specifying which mode to use in the basis
            int sigma_p   = Qx*p;
            int sigma_pq  = Qy*mode_pq[q + (Q+1)*p];
            int sigma_pqr = Qz*mode_pqr[r + (R+1)*(q + (Q+1)*p)];
            int sigma = Qx*Qy*Qz*mode;


            // Compute tensor product of inarray with the 3 basis functions
            Array<OneD, NekDouble> g_pqr = Array<OneD, NekDouble>( Qx*Qy*Qz, 0.0 );
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
        // Differentiation Methods
        //-----------------------------

        void StdTetExp::PhysDeriv(const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &out_d0, 
            Array<OneD, NekDouble> &out_d1,
            Array<OneD, NekDouble> &out_d2)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();

            ConstArray<OneD, NekDouble> z0,z1,z2;
            Array<OneD, NekDouble> d0;
            Array<OneD, NekDouble> wsp1  = Array<OneD, NekDouble>(nquad0*nquad1, 0.0);
            NekDouble *gfac = wsp1.get();

            z0 = ExpPointsProperties(0)->GetZ();
            z1 = ExpPointsProperties(1)->GetZ();
            z2 = ExpPointsProperties(2)->GetZ();

            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac[i] = 2.0/(1-z1[i]);
            }

            if(out_d1.num_elements() > 0)// if no d1 required do not need to calculate both deriv
            {
                PhysTensorDeriv(inarray, out_d0, out_d1, out_d2);

                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,gfac[i],&out_d0[0]+i*nquad0,1);
                }
            }
            else
            {
                if(out_d0.num_elements() > 0)// need other local callopsed derivative for d1 
                {
                    d0 = Array<OneD, NekDouble>(nquad0*nquad1, 0.0);
                }

                PhysTensorDeriv(inarray, d0, out_d1, out_d2);

                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,gfac[i],&d0[0]+i*nquad0,1);
                }

                // set up geometric factor: (1_z0)/(1-z1)
                for(i = 0; i < nquad0; ++i)
                {
                    gfac[i] = 0.5*(1+z0[i]);
                }

                for(i = 0; i < nquad1; ++i) 
                {
                    Vmath::Vvtvp(nquad0,gfac,1,&d0[0]+i*nquad0,1,&out_d1[0]+i*nquad0,1,
                        &out_d1[0]+i*nquad0,1);
                }    
            }
        }


        ///////////////////////////////
        // Evaluation Methods
        ///////////////////////////////

        void StdTetExp::BwdTrans(
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



            // Build an index map from the rectangle to the triangle
            Array<OneD, int> pq  = Array<OneD, int>( (P+1)*(Q+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q - p; ++q, ++mode ) {
                    pq[q + (Q+1)*p] = mode;
                }
            }

            // Create an index map from the hexahedron to the tetrahedron.
            // Explicit computation of the actual index is error-prone and too difficult.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q - p; ++q ) {
                    for( int r = 0; r <= R - p - q; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }



            // Sum-factorize the triple summation starting with the z-dimension
            for( int k = 0; k < Qz; ++k ) {

                // Create the matrix of coefficients summed over the z-modes
                Array<OneD, NekDouble> Ak((P+1)*(Q+1), 0.0);
                for( int p = 0; p <= P; ++p ) {
                    for( int q = 0; q <= Q - p; ++q ) {
                        for( int r = 0; r <= R - p - q; ++r ) {
                            int mode = pqr[r + (R+1)*(q + (Q+1)*p)];
                            Ak[q + (Q+1)*p]   +=   inarray[mode]  *  zBasis[k + Qz*mode];
                        }
                    }
                }

                // Factorize the y-dimension
                for( int j = 0; j < Qy; ++j ) {

                    // Create the vector of coefficients summed over the y and z-modes
                    Array<OneD, NekDouble> bjk(P+1, 0.0);
                    for( int p = 0; p <= P; ++p ) {
                        for( int q = 0; q <= Q - p; ++q ) {
                            int mode = pq[q + (Q+1)*p];
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

        void StdTetExp::FwdTrans(
            const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs,outarray);
            DNekVec out(m_ncoeffs,outarray,eWrapper);

            out = (*matsys)*in;
        }

        NekDouble StdTetExp::PhysEvaluate(const ConstArray<OneD, NekDouble>& coords)
        {
            Array<OneD, NekDouble> coll = Array<OneD, NekDouble>(2);

            // set up local coordinate system 
            if((fabs(coords[0]+1.0) < NekConstants::kEvaluateTol)
                &&(fabs(coords[1]-1.0) < NekConstants::kEvaluateTol))
            {
                coll[0] = 0.0;
                coll[1] = 1.0;
            }
            else
            {
                coll[0] = 2*(1+coords[0])/(1-coords[1])-1.0; 
                coll[1] = coords[1]; 
            }

            return  StdExpansion3D::PhysEvaluate(coll); 
        }

//         NekDouble StdTetExp::PhysEvaluate3D(const ConstArray<OneD, NekDouble>& coords)
//         {
//             NekDouble val;
//             int i;
//             int nq0 = m_base[0]->GetNumPoints();
//             int nq1 = m_base[1]->GetNumPoints();
//             Array<OneD, NekDouble> wsp1 = Array<OneD, NekDouble>(nq1, 0.0);
// 
//             DNekMatSharedPtr I;
// 
//             ASSERTL2(coords[0] < -1, "coord[0] < -1");
//             ASSERTL2(coords[0] > 1, "coord[0] >  1");
//             ASSERTL2(coords[1] < -1, "coord[1] < -1");
//             ASSERTL2(coords[1] > 1, "coord[1] >  1");
// 
//             // interpolate first coordinate direction
//             I = ExpPointsProperties(0)->GetI(coords);
//             for (i = 0; i < nq1;++i)
//         {
//                 wsp1[i] = Blas::Ddot(nq0, &(I->GetPtr())[0], 1,&m_phys[i * nq0], 1);
//         }
// 
//             // interpolate in second coordinate direction
//             I = ExpPointsProperties(1)->GetI(coords+1);
// 
//             val = Blas::Ddot(nq1, &(I->GetPtr())[0], 1, &wsp1[0], 1);
// 
//             return val;
//         }


        void StdTetExp::WriteToFile(std::ofstream &outfile)
        {
            int  i,j;
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            ConstArray<OneD, NekDouble> z0,z1;

            z0 = ExpPointsProperties(0)->GetZ();
            z1 = ExpPointsProperties(1)->GetZ();

            outfile << "Variables = z1,  z2, Coeffs \n" << std::endl;      
            outfile << "Zone, I=" << nquad0 <<", J=" << nquad1 <<", F=Point" << std::endl;

            for(j = 0; j < nquad1; ++j)
            {
                for(i = 0; i < nquad0; ++i)
                {
                    outfile << 0.5*(1+z0[i])*(1.0-z1[j])-1 <<  " " << 
                        z1[j] << " " << m_phys[j*nquad0+i] << std::endl;
                }
            }

        }

        //   I/O routine
        void StdTetExp::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int  i,j;
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  cnt = 0;
            Array<OneD, NekDouble> wsp  = Array<OneD, NekDouble>(order0*order1, 0.0);

            NekDouble *mat = wsp.get(); 

            // put coeffs into matrix and reverse order so that p index is fastest
            // recall q is fastest for tri's

            Vmath::Zero(order0*order1,mat,1);

            for(i = 0; i < order0; ++i)
            {
                for(j = 0; j < order1-i; ++j,cnt++)
                {
                    mat[i+j*order1] = m_coeffs[cnt];
                }
            }

            outfile <<"Coeffs = [" << " "; 

            for(j = 0; j < order1; ++j)
            {
                for(i = 0; i < order0; ++i)
                {
                    outfile << mat[j*order0+i] <<" ";
                }
                outfile << std::endl; 
            }
            outfile << "]" ; 
        }

        void StdTetExp::GetCoords( Array<OneD, NekDouble> & xi_x, 
            Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z)
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
                        xi_x[s] = (eta_x[i] + 1.0) * (1.0 - eta_y[j]) * (1.0 - eta_z[k]) / 4  -  1.0;
                        xi_y[s] = (eta_y[j] + 1.0) * (1.0 - eta_z[k]) / 2  -  1.0;
                        xi_z[s] = eta_z[k];
                    }
                }
            }
        }


//            StdMatrix StdTetExp::s_elmtmats;
//     
//         void StdTetExp::SetInvInfo(StdMatContainer *mat, MatrixType Mform)
//         {
//             mat->SetLda(m_ncoeffs);
//             mat->SetMatForm(eSymmetric_Positive);
//             
//             if(GeoFacType() == eRegular)
//             {
//                 switch(Mform)
//                 {
//                 case eMassMatrix: // definitions need adding 
//                     mat->SetMatForm(eSymmetric);
//                     break;
//                 case eLapMatrix:
//                     mat->SetMatForm(eSymmetric);
//                     break;
//                 default:
//                     ASSERTL0(false, "MatrixType not known");
//                     break;
//                 }
//             }
//         }
      }//end namespace
}//end namespace

/** 
 * $Log: StdTetExp.cpp,v $
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





