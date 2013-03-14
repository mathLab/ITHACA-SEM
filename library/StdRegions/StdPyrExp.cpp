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
        
        StdPyrExp::StdPyrExp() // Deafult construct of standard expansion directly called. 
        {
        }
        
        StdPyrExp::StdPyrExp(const LibUtilities::BasisKey &Ba,
                             const LibUtilities::BasisKey &Bb,
                             const LibUtilities::BasisKey &Bc) 
            : StdExpansion  (LibUtilities::StdPyrData::getNumberOfCoefficients(Ba.GetNumModes(),
                                                                 Bb.GetNumModes(),
                                                                 Bc.GetNumModes()),
                             3, Ba, Bb, Bc),
              StdExpansion3D(LibUtilities::StdPyrData::getNumberOfCoefficients(Ba.GetNumModes(),
                                                                 Bb.GetNumModes(),
                                                                 Bc.GetNumModes()),
                             Ba, Bb, Bc)
        {

            if (Ba.GetNumModes() > Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher "
                         "than order in 'c' direction");
            }
            if (Bb.GetNumModes() > Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'b' direction is higher "
                         "than order in 'c' direction");
            }
        }

        StdPyrExp::StdPyrExp(const StdPyrExp &T)
            : StdExpansion  (T),
              StdExpansion3D(T)
        {
        }


        // Destructor
        StdPyrExp::~StdPyrExp()
        {   
        } 

        //---------------------------------------
        // Integration/public 3D methods
        //---------------------------------------
        void StdPyrExp::TripleTensorProduct(const Array<OneD, const NekDouble>& fx, 
                                              const Array<OneD, const NekDouble>& gy, 
                                              const Array<OneD, const NekDouble>& hz, 
                                              const Array<OneD, const NekDouble>& inarray, 
                                              Array<OneD, NekDouble> & outarray)
        {
            // Using matrix operation, not sum-factorization.  Regarding the
            // 3D array, inarray[k][j][i], x is changing the fastest and z the
            // slowest.  Thus, the first x-vector of points refers to the
            // first row of the first stack. The first y-vector refers to the
            // first column of the first stack. The first z-vector refers to
            // the vector of stacks intersecting the first row and first
            // column. So in C++, i refers to column, j to row, and k to
            // stack.  Contrasting this with the usual C++ matrix convention,
            // note that i does not refer to a C++ row, nor j to C++ column.

            int nx = fx.num_elements();
            int ny = gy.num_elements();
            int nz = hz.num_elements();
            
            // Multiply by integration constants...  Hadamard multiplication
            // refers to elementwise multiplication of two vectors.  Hadamard
            // each row with the first vector (x-vector); the index i is
            // changing the fastest.
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
        
        // Inner-Product with respect to the weights: i.e., this is the triple
        // sum of the product of the four inputs over the prism. x-dimension
        // is the row, it is the index that changes the fastest; y-dimension
        // is the column; z-dimension is the stack, it is the index that
        // changes the slowest.
        NekDouble StdPyrExp::TripleInnerProduct(const Array<OneD, const NekDouble>& fxyz, 
                                                  const Array<OneD, const NekDouble>& wx, 
                                                  const Array<OneD, const NekDouble>& wy, 
                                                  const Array<OneD, const NekDouble>& wz)
        {
            int Qx = wx.num_elements();
            int Qy = wy.num_elements();
            int Qz = wz.num_elements();

            if (fxyz.num_elements() != Qx*Qy*Qz) 
            {
                cerr << "TripleInnerProduct expected " << fxyz.num_elements() 
                     << " quadrature points from the discretized input function but got " 
                     << Qx*Qy*Qz << " instead." << endl;
            }
            
            // Sum-factorizing over the stacks
            Array<OneD, NekDouble> A(Qx*Qy, 0.0);
            for (int i = 0; i < Qx; ++i)
            {
                for (int j = 0; j < Qy; ++j)
                {
                    for (int k = 0; k < Qz; ++k)
                    {
                        A[i + Qx*j] += fxyz[i + Qx*(j + Qy*k)] * wz[k];
                    }
                }
            }
            
            // Sum-factorizing over the columns
            Array<OneD, NekDouble> b(Qx, 0.0);
            for (int i = 0; i < Qx; ++i)
            {
                for (int j = 0; j < Qy; ++j)
                {
                    b[i] += A[i + Qx*j] * wy[j];
                }
            }
            
            // Sum-factorizing over the rows
            NekDouble c = 0;
            for (int i = 0; i < Qx; ++i)
            {
                c += b[i] * wx[i];
            }
            
            return c;
        }

        NekDouble StdPyrExp::Integral3D(
            const Array<OneD, const NekDouble>& inarray, 
            const Array<OneD, const NekDouble>& wx,
            const Array<OneD, const NekDouble>& wy, 
            const Array<OneD, const NekDouble>& wz)
        {
            return TripleInnerProduct(inarray, wx, wy, wz);
        }
        
        void StdPyrExp::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int order0 = m_base[0]->GetNumModes();
            int order1 = m_base[1]->GetNumModes();
            int order2 = m_base[2]->GetNumModes();

            Array<OneD, NekDouble> wsp(order0*order1*order2, 0.0);

            NekDouble *mat = wsp.get(); 

            // put coeffs into matrix and reverse order so that r index is
            // fastest for Prism
            Vmath::Zero(order0*order1*order2, mat, 1);

            for(int i = 0, cnt=0; i < order0; ++i)
            {
                for(int j = 0; j < order1; ++j)
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
                        outfile << mat[i + order0*(j + order1*k)] << " ";
                    }
                    outfile << std::endl; 
                }
            }
            outfile << "]"; 
        }

        //---------------------------------------
        // Differentiation/integration Methods
        //---------------------------------------
        
        /**
         * \brief Calculate the derivative of the physical points 
         *  
         * The derivative is evaluated at the nodal physical points.
         * Derivatives with respect to the local Cartesian coordinates.
         *  
         * \f$\begin{Bmatrix} \frac {\partial} {\partial \xi_1} \\ \frac
         * {\partial} {\partial \xi_2} \\ \frac {\partial} {\partial \xi_3}
         * \end{Bmatrix} = \begin{Bmatrix} \frac 2 {(1-\eta_3)} \frac \partial
         * {\partial \bar \eta_1} \\ \frac {\partial} {\partial \xi_2} \ \
         * \frac {(1 + \bar \eta_1)} {(1 - \eta_3)} \frac \partial {\partial
         * \bar \eta_1} + \frac {\partial} {\partial \eta_3} \end{Bmatrix}\f$
         */
        void StdPyrExp::v_PhysDeriv(
            const Array<OneD, const NekDouble> &u_physical,
                  Array<OneD,       NekDouble> &out_dxi1,
                  Array<OneD,       NekDouble> &out_dxi2,
                  Array<OneD,       NekDouble> &out_dxi3)
        {
            // PhysDerivative implementation based on Spen's book page 152.
            int    Qx = m_base[0]->GetNumPoints();
            int    Qy = m_base[1]->GetNumPoints();
            int    Qz = m_base[2]->GetNumPoints();

            Array<OneD, NekDouble> dEta_bar1(Qx*Qy*Qz,0.0);
            Array<OneD, NekDouble> dXi2     (Qx*Qy*Qz,0.0);
            Array<OneD, NekDouble> dEta3    (Qx*Qy*Qz,0.0);
            PhysTensorDeriv(u_physical, dEta_bar1, dXi2, dEta3);

            Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
            eta_x = m_base[0]->GetZ();
            eta_y = m_base[1]->GetZ();
            eta_z = m_base[2]->GetZ();

            int i, j, k, n;

            for (k = 0, n = 0; k < Qz; ++k)
            {
                for (j = 0; j < Qy; ++j)
                {
                    for (i = 0; i < Qx; ++i, ++n)
                    {
                        if (out_dxi1.num_elements() > 0)
                            out_dxi1[n] = 2.0/(1.0 - eta_z[k]) * dEta_bar1[n];
                        if (out_dxi2.num_elements() > 0)
                            out_dxi2[n] = 2.0/(1.0 - eta_z[k]) * dXi2[n];
                        if (out_dxi3.num_elements() > 0)
                            out_dxi3[n] = (1.0+eta_x[i])/(1.0-eta_z[k])*dEta_bar1[n] +
                                (1.0+eta_y[j])/(1.0-eta_z[k])*dXi2[n] + dEta3[n];
                    } 
                }
            }
        }

        void StdPyrExp::v_PhysDeriv(const int dir,
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
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
        }


        void StdPyrExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                                             Array<OneD,       NekDouble> &out_d0,
                                             Array<OneD,       NekDouble> &out_d1,
                                             Array<OneD,       NekDouble> &out_d2)
        {
            StdPyrExp::v_PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }
        
        /** 
         * \brief Integrate the physical point list \a inarray over pyramidic
         * region and return the value.
         *
         * Inputs:\n
         *
         * - \a inarray: definition of function to be returned at quadrature
         *    point of expansion.
         *
         * Outputs:\n
         *
         * - returns \f$\int^1_{-1}\int^1_{-1}\int^1_{-1} u(\bar \eta_1,
         *   \eta_2, \eta_3) J[i,j,k] d \bar \eta_1 d \eta_2 d \eta_3\f$ \n
         *   \f$= \sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1} \sum_{k=0}^{Q_3 -
         *   1} u(\bar \eta_{1i}^{0,0},
         *   \eta_{2j}^{0,0},\eta_{3k}^{2,0})w_{i}^{0,0} w_{j}^{0,0} \hat
         *   w_{k}^{2,0} \f$ \n where \f$inarray[i,j, k] = u(\bar
         *   \eta_{1i},\eta_{2j}, \eta_{3k}) \f$, \n \f$\hat w_{k}^{2,0} =
         *   \frac {w^{2,0}} {2} \f$ \n and \f$ J[i,j,k] \f$ is the Jacobian
         *   evaluated at the quadrature point.
         */
        NekDouble StdPyrExp::v_Integral(
            const Array<OneD, const NekDouble>& inarray)
        {
            // Using implementation from page 146 of Spencer Sherwin's book.
            int Qz = m_base[2]->GetNumPoints();

            // Get the point distributions:
            // * x is assumed to be Gauss-Lobatto-Legendre (incl. -1 and 1)
            // * y is assumed to be Gauss-Lobatto-Legendre (incl. -1 and 1)
            Array<OneD, const NekDouble> z, wx, wy, wz;
            wx = m_base[0]->GetW();
            wy = m_base[1]->GetW();
            m_base[2]->GetZW(z,wz);

            Array<OneD, NekDouble> wz_hat = Array<OneD, NekDouble>(Qz, 0.0);

            // Convert wz into wz_hat, which includes the 1/2 scale factor.
            // Nothing else need be done if the point distribution is Jacobi
            // (2,0) since (1 - xi_z)^2 is already factored into the weights.
            // Note by coincidence, xi_y = eta_y, xi_z = eta_z (xi_z = z
            // according to our notation).
            switch(m_base[2]->GetPointsType())
            {
                // Common case
                case LibUtilities::eGaussRadauMAlpha2Beta0: // (2,0) Jacobi Inner product
                    Vmath::Smul(Qz, 0.25, (NekDouble *)wz.get(), 1, wz_hat.get(), 1);
                    break;
                
                // Corner cases
                case LibUtilities::eGaussLobattoLegendre:
                case LibUtilities::eGaussRadauMLegendre:
                    for (int k = 0; k < Qz; ++k)
                    {
                        wz_hat[k] = 0.25*(1.0-z[k])*(1.0-z[k]) * wz[k];
                    }
                    break;
                    
                default:
                    ASSERTL0(false, "Unsupported quadrature points type.");
                    break;
            }
            
            return Integral3D(inarray, wx, wy, wz_hat);
        }


        //---------------------------------------
        // Transforms
        //---------------------------------------
        
	/** 
         * \brief Backward transformation is evaluated at the quadrature
         * points.
         *
         * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat
         * u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
         * 
         * Backward transformation is three dimensional tensorial expansion
         *
         * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a
	 *  (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
	 *  \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c (\xi_{3k})
	 *  \rbrace} \rbrace}. \f$
	 *
         * And sumfactorizing step of the form is as:\ \ \f$ f_{pqr}
         * (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c
         * (\xi_{3k}),\\ g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y}
         * \psi_{p}^a (\xi_{2j}) f_{pqr} (\xi_{3k}),\\ u(\xi_{1i}, \xi_{2j},
         * \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a (\xi_{1i}) g_{p}
         * (\xi_{2j}, \xi_{3k}).  \f$
         **/
        void StdPyrExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray)
        {
            ASSERTL1(m_base[1]->GetBasisType() != LibUtilities::eOrtho_B   ||
                     m_base[1]->GetBasisType() != LibUtilities::eModified_B,
                     "Basis[1] is not a general tensor type");

            ASSERTL1(m_base[2]->GetBasisType() != LibUtilities::eOrtho_C   ||
                     m_base[2]->GetBasisType() != LibUtilities::eModified_C,
                     "Basis[2] is not a general tensor type");

            int     Qx = m_base[0]->GetNumPoints();
            int     Qy = m_base[1]->GetNumPoints();
            int     Qz = m_base[2]->GetNumPoints();

            int     P = m_base[0]->GetNumModes() - 1;
            int     Q = m_base[1]->GetNumModes() - 1;
            int     R = m_base[2]->GetNumModes() - 1;

            Array<OneD, const NekDouble> bx = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> by = m_base[1]->GetBdata();
            Array<OneD, const NekDouble> bz = m_base[2]->GetBdata();

            /*
            // Create an index map from the hexahedron to the pyramid.
            Array<OneD, int> pqr = Array<OneD, int>( (P+1)*(Q+1)*(R+1), -1 );
            for( int p = 0, mode = 0; p <= P; ++p ) {
                for( int q = 0; q <= Q; ++q ) {
                    for( int r = 0; r <= R - p - q ; ++r, ++mode ) {
                        pqr[r + (R+1)*(q + (Q+1)*p)] = mode;
                    }
                }
            }
            */

            for (int k = 0; k < Qz; ++k) {
                for (int j = 0; j < Qy; ++j) {
                    for (int i = 0; i < Qx; ++i) {
                        NekDouble sum = 0.0;
                        for (int r = 0; r <= R; ++r) {
                            for (int q = 0; q <= min(R-r,Q); ++q) {
                                for (int p = 0; p <= min(R-r,P); ++p) {
                                    int mode = GetMode(p,q,r);
                                    sum += inarray[mode]*
                                        bx[i + Qx*p]*
                                        by[j + Qy*q]*
                                        bz[k + Qz*GetTetMode(p,q,r)];
                                }
                            }
                        }
                        
                        // Add in contributions from singular vertices;
                        // i.e. (p,q,r) = (1,1,1),(0,1,1),(1,0,1)
                        int m = GetMode(0,0,1);
                        sum += inarray[m]*bz[k*Qz]*(bx[i+Qx]*by[j+Qy]+
                                                    bx[i   ]*by[j+Qy]+
                                                    bx[i+Qx]*by[j   ]);
                        outarray[i + Qx*(j + Qy*k)] = sum;
                    }
                }
            }
            
            /*
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
            */
        }


	/** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a (this)->m_coeffs  
            
            Inputs:\n
            
            - \a inarray: array of physical quadrature points to be transformed
            
            Outputs:\n
            
            - (this)->_coeffs: updated array of expansion coefficients. 
            
        */    
        void StdPyrExp::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray)
        {
            v_IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs, outarray);
            DNekVec out(m_ncoeffs, outarray, eWrapper);

            out = (*matsys)*in;
        }
        
        
        //---------------------------------------
        // Inner product functions
        //---------------------------------------

        /** \brief  Inner product of \a inarray over region with respect to the 
            expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata() and return in \a outarray 
            
            Wrapper call to StdPyrExp::IProductWRTBase
            
            Input:\n
            
            - \a inarray: array of function evaluated at the physical collocation points
            
            Output:\n
            
            - \a outarray: array of inner product with respect to each basis over region
            
        */
        void StdPyrExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble> &inarray, 
                  Array<OneD,       NekDouble> &outarray)
        {
            int P = m_base[0]->GetNumModes()-1;
            int Q = m_base[1]->GetNumModes()-1;
            int R = m_base[2]->GetNumModes()-1;

            int Qx = m_base[0]->GetNumPoints();
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            const Array<OneD, const NekDouble> &bx = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble> &by = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble> &bz = m_base[2]->GetBdata();

            for( int r = 0; r <= R; ++r ) {
                for( int q = 0; q <= min(R-r,Q); ++q ) {
                    for( int p = 0; p <= min(R-r,P); ++p ) {
                        // Compute tensor product of inarray with the 3 basis functions
                        Array<OneD, NekDouble> g_pqr = Array<OneD, NekDouble>( Qx*Qy*Qz, 0.0 );
                        for( int k = 0; k < Qz; ++k ) {
                            for( int j = 0; j < Qy; ++j ) {
                                for( int i = 0; i < Qx; ++i ) {
                                    int s = i + Qx*(j + Qy*k);
                                    cout << p << " " << q << " " << r << " " << GetTetMode(p,q,r) << endl;
                                    g_pqr[s] += inarray[s] * 
                                        bx[i + Qx*p] * 
                                        by[j + Qy*q] * 
                                        bz[k + Qz*GetTetMode(p,q,r)];
                                    
                                    if (p == 0 && q == 0 && r == 1)
                                    {
                                        g_pqr[s] += inarray[s] * bz[k+Qz]*(
                                            bx[i+Qx]*by[j+Qy] + 
                                            bx[i+Qx]*by[j   ] + 
                                            bx[i   ]*by[j+Qy]);
                                    }
                                }
                            }
                        }
                        
                        outarray[GetMode(p,q,r)] = Integral( g_pqr );
                   }
                }
            }
        }
        
        
        //---------------------------------------
        // Evaluation functions
        //---------------------------------------
        
        NekDouble StdPyrExp::v_PhysEvaluate(const Array<OneD, const NekDouble>& xi)
        {
            return PhysEvaluate(xi, m_phys);
        }

        NekDouble StdPyrExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble>& xi,
            const Array<OneD, const NekDouble>& physvals)
        {
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);

            if (fabs(xi[2]-1.0) < NekConstants::kNekZeroTol)
            {
                // Very top point of the pyramid
                eta[0] = -1.0;
                eta[1] = -1.0;
                eta[2] = xi[2];
            }
            else  
            {
                // Below the line-singularity -- Common case
                eta[2] = xi[2]; // eta_z = xi_z
                eta[1] = 2.0*(1.0 + xi[1])/(1.0 - xi[2]) - 1.0; 
                eta[0] = 2.0*(1.0 + xi[0])/(1.0 - xi[2]) - 1.0;
            } 
            
            return StdExpansion3D::v_PhysEvaluate(eta, physvals);
        }

        void StdPyrExp::v_GetCoords(Array<OneD, NekDouble> &xi_x, 
                                    Array<OneD, NekDouble> &xi_y,
                                    Array<OneD, NekDouble> &xi_z)
        {
            Array<OneD, const NekDouble> etaBar_x = m_base[0]->GetZ();
            Array<OneD, const NekDouble> eta_y    = m_base[1]->GetZ();
            Array<OneD, const NekDouble> eta_z    = m_base[2]->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates: eta --> xi
            for (int k = 0; k < Qz; ++k )
            {
                for (int j = 0; j < Qy; ++j) 
                {
                    for (int i = 0; i < Qx; ++i) 
                    {
                        int s = i + Qx*(j + Qy*k);

                        xi_z[s] = eta_z[k];
                        xi_y[s] = (1.0 + eta_y[j]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                        xi_x[s] = (1.0 + etaBar_x[i]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                    }
                }
            }
        }

        void StdPyrExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int Qx = m_base[0]->GetNumPoints();
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();
            int P  = m_base[0]->GetNumModes() - 1;
            int Q  = m_base[1]->GetNumModes() - 1;
            int R  = m_base[2]->GetNumModes() - 1;

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
        
        
        //---------------------------------------
        // Helper functions
        //---------------------------------------

        int StdPyrExp::v_GetNverts() const
        {
            return 5;
        }

        int StdPyrExp::v_GetNedges() const
        {
            return 8;
        }

        int StdPyrExp::v_GetNfaces() const
        {
            return 5;
        }
        
        LibUtilities::ShapeType StdPyrExp::v_DetShapeType() const
        {
            return LibUtilities::ePyramid;
        }

        int StdPyrExp::v_NumBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_B ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            
            int P = m_base[0]->GetNumModes() - 1;
            int Q = m_base[1]->GetNumModes() - 1;
            int R = m_base[2]->GetNumModes() - 1;
            
            return (P+1)*(Q+1)              // 1 rect. face in p-q plane
                + 2*(R+1) + P*(1+2*R-P)     // 2 tri. faces in p-r plane
                + 2*(R+1) + Q*(1+2*R-Q)     // 2 tri. faces in q-r plane
                - 2*(P+1)-2*(Q+1)-4*(R+1)   // subtract double counted edges
                + 5;                        // add vertices
        }

        int StdPyrExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 7, "edge id is out of range");
            
            if (i == 0 || i == 2)
            {
                return GetBasisNumModes(0);
            }
            else if (i == 1 || i == 3)
            {
                return GetBasisNumModes(1);
            }
            else
            {
                return GetBasisNumModes(2);
            }
        }

        int StdPyrExp::v_GetFaceNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            
            if (i == 0)
            {
                return GetBasisNumModes(0)*GetBasisNumModes(1);
            }
            else if (i == 1 || i == 3)
            {
                return GetBasisNumModes(0)*GetBasisNumModes(2);
            }
            else
            {
                return GetBasisNumModes(1)*GetBasisNumModes(2);
            }
        }
        
        int StdPyrExp::v_GetFaceIntNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");

            int P = m_base[0]->GetNumModes()-1;
            int Q = m_base[1]->GetNumModes()-1;
            int R = m_base[2]->GetNumModes()-1;

            if (i == 0)
            {
                return (P-1)*(Q-1);
            }
            else if (i == 1 || i == 3)
            {
                return (P-1) * (2*(R-1) - (P-1) - 1) / 2;
            }
            else
            {
                return (Q-1) * (2*(R-1) - (Q-1) - 1) / 2;
            }
        }

        int StdPyrExp::v_CalcNumberOfCoefficients(
            const std::vector<unsigned int> &nummodes, 
            int &modes_offset)
        {
            int nmodes = LibUtilities::StdPyrData::getNumberOfCoefficients(
                nummodes[modes_offset],
                nummodes[modes_offset+1],
                nummodes[modes_offset+2]);
            
            modes_offset += 3;
            return nmodes;
        }
        
        LibUtilities::BasisType StdPyrExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 7, "edge id is out of range");
            if (i == 0 || i == 2)
            {
                return GetBasisType(0);
            }
            else if (i == 1 || i == 3)
            {
                return GetBasisType(1);
            }
            else
            {
                return GetBasisType(2);
            }
        }


        //---------------------------------------
        // Mappings
        //---------------------------------------
        
        void StdPyrExp::v_GetFaceToElementMap(
            const int                  fid, 
            const Orientation      faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray,
            int                        nummodesA, 
            int                        nummodesB)
        {
            int i,j,P,Q;
            const int nummodes0 = m_base[0]->GetNumModes();
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nummodes2 = m_base[2]->GetNumModes();
            //int nummodesA, nummodesB, P, Q;

            const LibUtilities::BasisType bType0 = GetEdgeBasisType(0);
            const LibUtilities::BasisType bType1 = GetEdgeBasisType(1);
            const LibUtilities::BasisType bType2 = GetEdgeBasisType(4);
            
            ASSERTL1( (bType0==bType1),
                      "Method only implemented if BasisType is indentical in x and y directions");
            ASSERTL1( (bType0==LibUtilities::eModified_A) && (bType1==LibUtilities::eModified_A) && (bType2==LibUtilities::eModified_C),
                      "Method only implemented for Modified_A BasisType (x and y direction) and Modified_C BasisType (z direction)");

            bool isQuad = true;

            int nFaceCoeffs = 0;
            if( fid == 0 ) // Base quad 
            {
                nummodesA = nummodes0;
                nummodesB = nummodes1;
                P = nummodesA-1;
                Q = nummodesB-1;
                nFaceCoeffs = nummodesA*nummodesB;
            }
//             else if((fid == 2) || (fid == 4)) // front and back quad
//             {
//                 nummodesA = nummodes1;
//                 nummodesB = nummodes2;
//                 P = nummodesA-1;
//                 Q = nummodesB-1;
//                 nFaceCoeffs = nummodesA*nummodesB;
//             }
            else  // left and right triangles
            {
                nummodesA = nummodes0;
                nummodesB = nummodes2;
                P = nummodesA-1;
                Q = nummodesB-1;
                nFaceCoeffs = Q+1 + (P*(1 + 2*Q - P))/2;
                isQuad = false;
            }

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



            Array<OneD, int> arrayindex(nFaceCoeffs,-1);

            for(int a = 0; a < nummodesA; ++a)
            {
                for(int b = 0; isQuad ? (b <  nummodesB) : (b < nummodesB - a); ++b)
                {
                    if( faceOrient < 9 ) // Not transposed
                    {
                        arrayindex[b + nummodesB*a] = b + nummodesB*a;
                    }
                    else // Transposed
                    {
                        arrayindex[b + nummodesB*a] = a + nummodesA*b;
                    }
                }
            }


            int baseCoefficient = 0;
            
            switch(fid)
            {

            // Base quad
            case 0: 
                for(int a = 0; a < nummodesA; ++a) {
                    for(int b = 0; b < nummodesB; ++b) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = b + nummodesB*a;
                    }
                }
            break;
            
            // Rear triangle
            case 3:
                baseCoefficient = (nummodes1 - 1) * nummodes2;
                for(int a = 0; a < nummodesA; ++a) {
                    for(int b = 0; b <  nummodesB - a; ++b) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = baseCoefficient + b;
                    }
                    baseCoefficient += nummodes1*(nummodesB-1 - a)  +  1;
                }
            break;

            // Front triangle
            case 1: 
                for(int a = 0; a < nummodesA; ++a) {
                    for(int b = 0; b <  nummodesB - a; ++b) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = baseCoefficient + b;
                    }
                    baseCoefficient += nummodes1 * (nummodes2 - a);
                }
            break;


            // Vertical triangle
            case 4: 
                for(int a = 0, n = 0; a < nummodesA; ++a) {
                    for(int b = 0; b < nummodesB - a; ++b, ++n) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = n;   
                    }
                }
            break;
            

            // Slanted triangle
            case 2:
                for(int b = nummodesB-1; b >= 0; --b) {
                    for(int a = 0; a < nummodesA - b; ++a) {
                        ASSERTL0(arrayindex[b + nummodesB*a] != -1, "arrayindex is not set up properly.");
                        maparray[ arrayindex[b + nummodesB*a] ] = baseCoefficient + (a+1)*(b+1) - 1;
                    }
                    baseCoefficient += nummodesA*(b+1);
                }
            break;
            
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
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
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
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
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
                            if( arrayindex[i*nummodesA+j] >= 0 )
                                signarray[ arrayindex[i*nummodesA+j] ] *= -1;
                        }
                    }                 
                        
                    for(i = 0; i < nummodesB; i++)
                    {
                        swap( maparray[i*nummodesA] , maparray[i*nummodesA+1] );
                        swap( signarray[i*nummodesA] , signarray[i*nummodesA+1] );
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
                        swap( maparray[i*nummodesB] , maparray[i*nummodesB+1] );
                        swap( signarray[i*nummodesB] , signarray[i*nummodesB+1] );
                    }
                }
            }      

        }

        int StdPyrExp::v_GetVertexMap(int vId)
        {
            ASSERTL0(GetEdgeBasisType(vId) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(vId) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(vId) == LibUtilities::eModified_B,
                     "Mapping not defined for this type of basis");
            
            int l = 0;
            
            switch (vId)
            {
                case 0:
                    l = GetMode(0,0,0);
                    break;
                case 1:
                    l = GetMode(1,0,0);
                    break;
                case 2:
                    l = GetMode(1,1,0);
                    break;
                case 3:
                    l = GetMode(0,1,0);
                    break;
                case 4:
                    l = GetMode(0,0,1);
                    break;
                default:
                    ASSERTL0(false, "local vertex id must be between 0 and 4");
            }
            
            return l;
        }
        

        //---------------------------------------
        // Wrapper functions
        //---------------------------------------
        
        DNekMatSharedPtr StdPyrExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            return CreateGeneralMatrix(mkey);
        }
        
        DNekMatSharedPtr StdPyrExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return v_GenMatrix(mkey);
        }

        /*
        void StdPyrExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
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
        */

        /**
         * @brief Compute the local mode number in the expansion for a
         * particular tensorial combination.
         *
         * Modes are numbered with the r index travelling fastest, followed by
         * q and then p, with each plane of size (Q+1)*[(Q+2)/2+R-Q-p]. For
         * example, with P=3, Q=3, R=4, the indexing inside each q-r plane
         * (with r increasing upwards and q to the right) is:
         * 
         * 4            *            *            *
         * 3  8         17 21        *  *         *  *
         * 2  7  11     16 20 24     29 32 35     *  *  *
         * 1  6  10 13  15 19 23 26  28 31 34 37  39 41 43 45
         * 0  5   9 12  14 18 22 25  27 30 33 36  38 40 42 44
         * 
         * Note that in the pyramid, then number of modes needed to perform
         * the boundary-interior decomposition is larger than the space of
         * polynomials used to represent it.
         */
        int StdPyrExp::GetMode(int p, int q, int r)
        {
            int Q = m_base[1]->GetNumModes() - 1;
            int R = m_base[2]->GetNumModes() - 1;
            
            return p*(Q+1)*((Q+2)/2+(R-Q)) - (p-1)*p*(p+1)/6 + // skip to p-th plane
                r+q*(2*(R+1)+1-q)/2 -                          // normal tri indexing
                // account for offset (starred points)
                (q > 0 ? p*(p+1)/2 - (q-1 < p ? (p-q+1)*(p-q)/2 : 0) : 0);
        }

        /**
         * 3
         * 2  6          12
         * 1  5  8       11 14      17
         * 0  4  7  9    10 13 15   16  18   19
         */

        int StdPyrExp::GetTetMode(const int I, const int J, const int K)
        {
            const int P = m_base[0]->GetNumModes();
            const int Q = m_base[1]->GetNumModes();
            const int R = m_base[2]->GetNumModes();

            int i,j,q_hat,k_hat;
            int cnt = 0;
            
            // Skip along the stacks (K)
            for (i = 0; i < I; ++i)
            {
                q_hat = min(Q,P-i);
                k_hat = min(R-Q, max(0, R-i));
                cnt += q_hat*(q_hat+1)/2 - k_hat*Q;
            }
            
            // Skip across the columns (J)
            q_hat = min(Q,P-I);
            k_hat = min(R-Q, max(0, R-I));
            for (j = 0; j < J; ++j)
            {
                cnt += q_hat + k_hat - j;
            }
            
            // Skip up the columns (K)
            cnt += K;
            
            // Return the final mode number
            return cnt;
        }

        void StdPyrExp::MultiplyByQuadratureMetric(
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
            
            const Array<OneD, const NekDouble>& z2 = m_base[2]->GetZ();
            
            // Multiply by integration constants in x-direction
            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vmul(nquad0, inarray.get()+i*nquad0, 1,
                            w0.get(), 1, outarray.get()+i*nquad0,1);
            }
            
            // Multiply by integration constants in y-direction
            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,w1[i], &outarray[0]+i*nquad0 +
                                j*nquad0*nquad1,1);
                }
            }
            
            // Multiply by integration constants in z-direction; need to
            // incorporate factor [(1-eta_3)/2]^2 into weights, but only if
            // using GLL quadrature points.
            switch(m_base[2]->GetPointsType())
            {
                // Legendre inner product.
                case LibUtilities::eGaussLobattoLegendre:
                    for(i = 0; i < nquad2; ++i)
                    {
                        Blas::Dscal(nquad0*nquad1,0.125*(1-z2[i])*(1-z2[i])*w2[i],
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
                    ASSERTL0(false, "Quadrature point type not supported for this element.");
                    break;
            }
        }
    }
}
