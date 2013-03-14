///////////////////////////////////////////////////////////////////////////////
//
// File NodalUtil.cpp
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
// Description: 2D Nodal Triangle Fekete Utilities --
//              Basis function, Interpolation, Integral, Derivation, etc.                
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <limits>

#include <LibUtilities/Foundations/NodalUtil.h>

#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

namespace Nektar
{
    namespace LibUtilities
    {

        template< typename T > NekVector<T> GetColumn(const NekMatrix<T> & matA, int n)
        {
            NekVector<T> v(matA.GetRows());
            for( int i=0; i<matA.GetRows(); ++i )
            {
                v[i] = matA(i,n);
            }
            return v;
        }
    
        NekMatrix<NekDouble> & SetColumn(NekMatrix<NekDouble> & matA, int n, const NekVector<NekDouble> & x)
        {
            for( int i=0; i<int(matA.GetRows()); ++i )
            {
                matA(i,n) = x[i];
            }
            return matA;
        }
    
        // Standard basis vector in R^M
        NekVector<NekDouble> GetE(int rows, int n)
        {
            NekVector<NekDouble> e(rows, 0.0);
            e(n) = 1;
            return e;
        }
    
        NekMatrix<NekDouble> Invert(const NekMatrix<NekDouble> & matA)
        {
            int rows = matA.GetRows(), columns = matA.GetColumns();
            NekMatrix<NekDouble> matX(rows,columns);
    
            // The linear system solver
            LinearSystem matL( SharedNekMatrixPtr(new NekMatrix<NekDouble>(matA)) );
    
            // Solve each column for the identity matrix
            for( int j=0; j<columns; ++j )
            {
                SetColumn( matX, j, matL.Solve( GetE(rows,j) ) );
            }
            
            return matX;
        }
        
        NekMatrix<NekDouble> GetTranspose(const NekMatrix<NekDouble> & matA)
        {
            int rows = matA.GetRows(), columns = matA.GetColumns();
            NekMatrix<NekDouble> matX(rows,columns);
        
            for( int i=0; i<rows; ++i )
            {
                for( int j=0; j<columns; ++j )
                {
                    matX(j,i) = matA(i,j);
                }
            }
            return matX;
        }
    
        int GetSize(const Array<OneD, const NekDouble> & x)
        {
            return x.num_elements();
        }
        
        int GetSize(const NekVector<NekDouble> & x)
        {
            return x.GetRows();
        }
                
        NekVector<NekDouble> ToVector( const Array<OneD, const NekDouble> & x )
        {
            return NekVector<NekDouble>( GetSize(x), x.data() );
        }
        
        Array<OneD, NekDouble> ToArray( const NekVector<NekDouble> & x )
        {
            return Array<OneD, NekDouble>( GetSize(x), x.GetPtr() );
        }
    
        NekVector<NekDouble> Hadamard( const NekVector<NekDouble> & x, const NekVector<NekDouble> & y )
        {
            int size = GetSize(x);
            NekVector<NekDouble> z(size);
    
            for( int i=0; i<size; ++i )
            {
                z[i] = x[i] * y[i];
            }
            return z;
        }

        NekVector<NekDouble> VectorPower( const NekVector<NekDouble> & x, NekDouble p )
        {
            int size = GetSize(x);
            NekVector<NekDouble> z(size);
    
            for( int i=0; i<size; ++i )
            {
                z[i] = pow( x[i], p );
            }
            
            return z;
        }

        // Formatted version of matrix ostream output
        std::string MatrixToString( const NekMatrix<NekDouble> & A, int precision, NekDouble expSigFigs )
        {
            stringstream s;
            s << setprecision(precision);
            int M = int(A.GetRows()), N = int(A.GetColumns());
            
            for(int i=0; i<M; ++i )
            {
                for(int j=0; j<N; ++j)
                {
                    NekDouble a = MakeRound(expSigFigs * A(i, j)) / expSigFigs;
                    s << setw(7) << right << a;
                    if( j < N-1 )
                    {
                        s << ", ";
                    }
                }
                if( i < M-1 )
                {
                    s << "\n";
                }
            }
            return s.str();
        }

        // Formatted version of vector ostream output
        std::string VectorToString( const NekVector<NekDouble> & v, int precision, NekDouble expSigFigs )
        {
            stringstream s;
            s << setprecision(precision) << "[ ";
            int N = int(v.GetRows());
            for(int j=0; j<N; ++j )
            {
                NekDouble x = MakeRound(expSigFigs * v(j)) / expSigFigs;
                s << setw(7) << right << x;
                if( j < N-1 )
                {
                    s << ", ";
                }
            }
            s << " ]";
            return s.str();
        }
                    

        int GetTriNumPoints(int degree)
        {
            return (degree+1) * (degree+2) / 2;
        }
        
        int GetDegree(int nBasisFunctions)
        {
            int degree = int(MakeRound((-3.0 + sqrt(1.0 + 8.0*nBasisFunctions))/2.0));

            // TODO: Find out why ASSERTL0 and ASSERTL1 don't work
            ASSERTL1( GetTriNumPoints(degree) == nBasisFunctions, "The number of points defines an expansion of fractional degree, which is not supported." );
            return degree;
        }
        
        int GetTetNumPoints(int degree)
        {
            return (degree+1) * (degree+2) * (degree+3) / 6;
        }

        // Get Tetrahedral number, where Tn = (d+1)(d+2)(d+3)/6
        int GetTetDegree(int nBasisFunc)
        {
            NekDouble eq = pow( 81.0 * nBasisFunc + 3.0 * sqrt(-3.0 + 729.0 * nBasisFunc * nBasisFunc), 1.0/3.0);
            int degree = int(MakeRound(eq/3.0 + 1.0/eq - 1.0)) - 1;

            ASSERTL1( GetTetNumPoints(degree) == nBasisFunc, "The number of points defines an expansion of fractional degree, which is not supported." );
            return degree;
        }
        
        NekDouble MakeRound(NekDouble x)
        {
            return floor(x + 0.5);
        }
        
        NekVector<NekDouble> MakeDubinerQuadratureSystem(int nBasisFunctions)
        {
            // Make the vector of integrals: note that each modal basis function integrates to zero except for the 0th degree
            NekVector<NekDouble> g(nBasisFunctions, 0.0);
            g(0) = sqrt(2.0);
    
            return g;
        }

        NekVector<NekDouble> MakeTetQuadratureSystem(int nBasisFunctions)
        {
            NekVector<NekDouble> g(nBasisFunctions, 0.0);
            g(0) = 1.0;
            return g;
        }
    
    
        NekVector<NekDouble> JacobiPoly(int degree, const NekVector<NekDouble>& x, NekDouble alpha, NekDouble beta)
        {
            int size = GetSize(x);
            NekVector<NekDouble> y(size);
            
            if(degree == 0)
            {
                // Set y to ones
                y = NekVector<NekDouble>(size, 1.0);
    
            }
            else if (degree == 1)
            {           
                for(int el=0; el<size; ++el)
                {
                    y[el] = 0.5*(alpha - beta + (alpha + beta + 2.0) * x[el]);
                }
                
            }
            else if (degree > 1)
            {           
                NekDouble degm1 = degree - 1.0;
                NekDouble tmp = 2.0 * degm1 + alpha + beta;
                NekDouble a1 = 2.0 * (degm1 + 1.0) * (degm1 + alpha + beta + 1.0) * tmp;
                NekDouble a2 = (tmp + 1.0) * (alpha * alpha - beta * beta);
                NekDouble a3 = tmp * (tmp + 1.0) * (tmp + 2.0);
                NekDouble a4 = 2.0 * (degm1 + alpha) * (degm1 + beta) * (tmp + 2.0);
    
                // TODO: Make this efficient: Reuse the degree-2 for the degree-1 call.
                // This can be done in linear time, but it is currently implemented to run in exponential time.
                NekVector<NekDouble> z2 = JacobiPoly(degree-2, x, alpha, beta);
                NekVector<NekDouble> z1 = JacobiPoly(degree-1, x, alpha, beta);

                for (int el=0; el<size; ++el)
                {
                    y[el] = ((a2 + a3 * x[el]) * z1[el] - a4 * z2[el])/a1;
                }
                
            }
            else
            {
                cerr << "Bad degree" << endl;
            }
    
            return y;
        }
        
        NekDouble JacobiPoly(int degree, NekDouble x, NekDouble alpha, NekDouble beta)
        {
            NekDouble y = 0.0;
            
            if(degree == 0)
            {
                y = 1.0;
            }
            
            else if (degree == 1)
            {           
                y = 0.5*(alpha - beta + (alpha + beta + 2.0) * x);
            }
            
            else if (degree > 1)
            {           
                NekDouble degm1 = degree - 1.0;
                NekDouble tmp = 2.0 * degm1 + alpha + beta;
                NekDouble a1 = 2.0 * (degm1 + 1.0) * (degm1 + alpha + beta + 1.0) * tmp;
                NekDouble a2 = (tmp + 1.0) * (alpha * alpha - beta * beta);
                NekDouble a3 = tmp * (tmp + 1.0) * (tmp + 2.0);
                NekDouble a4 = 2.0 * (degm1 + alpha) * (degm1 + beta) * (tmp + 2.0);

                // TODO: Make this efficient: Reuse the degree-2 for the degree-1 call.
                // This can be done in linear time, but it is currently implemented to run in exponential time.
                NekDouble z2 = JacobiPoly(degree-2, x, alpha, beta);
                NekDouble z1 = JacobiPoly(degree-1, x, alpha, beta);

                y = ((a2 + a3 * x) * z1 - a4 * z2)/a1;
            }
            
            else
            {
                cerr << "Bad degree" << endl;
            }
    
            return y;
        }
    
    
        NekVector<NekDouble> LegendrePoly(int degree, const NekVector<NekDouble>& x)
        {
            int size = GetSize(x);
            NekVector<NekDouble> y(size);
            
            if(degree > 1)
            {
                NekVector<NekDouble> a0(size, 1.0);
                NekVector<NekDouble> a1 = x;
                NekVector<NekDouble> a2(size);
    
                for(int i=2; i<=degree; ++i)
                {
                    NekDouble b = NekDouble(2.0*i-1.0)/i;
                    NekDouble c = NekDouble(i-1.0)/i;
    
                    // multiply each elements in matrix
                    a2 = b * Hadamard(a1, x)  -  c * a0;
                    a0 = a1;
                    a1 = a2;
                }
                y = a2;    
            }
            else if( degree == 1 )
            {
                y = x;                
            }
            else
            {
                y = NekVector<NekDouble>(size, 1.0);
            }
            return y;
        }

        // Triangle orthonormal basis expansion
        NekVector<NekDouble> DubinerPoly(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            // Get the coordinate transform
            int size = GetSize(x);
            NekVector<NekDouble> eta_1(size);
    
            // Initialize the horizontal coordinate of the triangle (beta in Barycentric coordinates)
            for(int el=0; el<size; ++el)
            {
                if( y[el] < 1.0 - 1e-16 )
                {
                    eta_1[el] = 2.0*(1.0 + x[el]) / (1.0 - y[el]) - 1.0;
                }
                else
                {
                    eta_1[el] = -1.0; // When y is close to 1, then we have a removeable singularity
                }
            }
    
            // Initialize the vertical coordinate of the triangle (gamma in Barycentric coordinates)
            NekVector<NekDouble> eta_2 = y;
    
            // Orthogonal Dubiner polynomial
            int alpha = 2*p + 1; int beta = 0;
            NekVector<NekDouble> phi(size);
            NekVector<NekDouble> psi_a = LegendrePoly(p, eta_1);
            NekVector<NekDouble> upsilon = JacobiPoly(q, eta_2, alpha, beta);
            NekDouble w = sqrt((2.0 * p + 1.0) * (p + q + 1.0) / 2.0); // Normalizing Orthonormal weight
    
            for(int el=0; el<size; ++el)
            {
                NekDouble zeta   = pow((1.0 - eta_2[el])/2.0, p);
                NekDouble psi_b  = zeta * upsilon[el];
                phi[el]          = w * psi_a[el] * psi_b;
            }
            return phi;
        }


        // Tetrahedral orthogonal basis expansion
        NekVector<NekDouble> TetrahedralBasis(int p, int q, int r, const NekVector<NekDouble>& x,
                                              const NekVector<NekDouble>& y, const NekVector<NekDouble>& z)
        {
            // Get the coordinate transform
            int size = GetSize(x);
            NekVector<NekDouble> eta_1(size), eta_2(size);
  
            // Initialize the horizontal coordinate of the Tetrahedral 
            for(int el=0; el<size; ++el)
            {
                if( z[el] < -y[el] - numeric_limits<NekDouble>::epsilon() )
                {
                    eta_1[el] = 2.0*(1.0 + x[el])/(-y[el]-z[el]) - 1.0;
                }
                else
                {
                    eta_1[el] = -1.0;
                }
            }

             // Initialize the  coordinate of the Tetrahedral 
            for(int el=0; el<size; ++el)
            {
                if( z[el] < 1.0 - numeric_limits<NekDouble>::epsilon())
                {
                    eta_2[el] = 2.0*(1.0 + y[el]) / (1.0 - z[el]) - 1.0;
                }
                else
                {
                    eta_2[el] = -1.0;  // When z is close to 1, then we have a removeable singularity
                    eta_1[el] = -1.0;  // When z is close to 1, then we have a removeable singularity
                }
            }
            
            // Initialize the vertical coordinate of the Tetrahedral 
            NekVector<NekDouble> eta_3 = z;

            // Orthogonal basis polynomial
            int alpha = 2*p + 1; int beta = 0;
            int alpha_r = 2*p + 2*q + 2; int beta_r = 0;
            NekVector<NekDouble> phi(size);
            NekVector<NekDouble> psi_a    = LegendrePoly(p,eta_1);
            NekVector<NekDouble> psi_bpq  = JacobiPoly(q, eta_2, alpha, beta);
            NekVector<NekDouble> psi_cpqr = JacobiPoly(r, eta_3, alpha_r, beta_r);
            NekDouble w = 1.0;

            for(int el=0; el<size; ++el)
            {
                NekDouble zeta_1 = pow((1.0 - eta_2[el])/2.0, p);                
                NekDouble zeta_2 = pow((1.0 - eta_3[el])/2.0, p+q);
                
                NekDouble psi_b = zeta_1 * psi_bpq[el];
                NekDouble psi_c = zeta_2 * psi_cpqr[el];

                phi[el]         = w * psi_a[el] * psi_b * psi_c;

            }
               
            return phi;

        }

        NekMatrix<NekDouble> GetTetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z, int degree)
        {
           //cout << "Begin  GetTetVandermonde" << endl;
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) * (degree + 3)/6;
            
            NekMatrix<NekDouble> matV(rows, cols, 0.0);

            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p <= d; ++p)
                {
                    for(int q=0; q <= d - p; ++q, ++n)
                    {
                        int r = d - p - q;
                        NekVector<NekDouble> columnVector = TetrahedralBasis(p, q, r, x, y, z);

                       //  cout << "degree = " << degree << ", (d,p,q,r) = (" << d << ", " << p << ", " << q << ", " << r << ")" << endl;
                        // Set n-th column of V to the TetrahedralBasis vector
                        for(int i=0; i<rows; ++i)
                        {
                            matV(i,n) = columnVector[i];
                        }
                    }
                }
            }
            return matV;
        }


        NekMatrix<NekDouble> GetTetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z)
        {
            int rows = GetSize(x);
            int degree = GetTetDegree(rows);
            return GetTetVandermonde(x, y, z, degree);
        }

        
        NekMatrix<NekDouble> GetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) / 2;
            
            NekMatrix<NekDouble> matV(rows, cols, 0.0);
    
            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p<=d; ++p, ++n)
                {
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPoly(p, q, x, y);
    
                    // Set n-th column of V to the DubinerPoly vector
                    for(int i=0; i<rows; ++i)
                    {
                            matV(i,n) = columnVector[i];
                    }
                }
            }
    
            return matV;
        }
    
        NekMatrix<NekDouble> GetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetVandermonde( x, y, degree );
        }
    
        SharedNekMatrixPtr MakeVmatrixOfTet(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z)
        {
            // Construct the Vandermonde matrix of Tetrahedron
            SharedNekMatrixPtr vMat(new NekMatrix<NekDouble>( GetTetVandermonde(x, y, z)));
            return vMat;
        }
        
        SharedNekMatrixPtr MakeVmatrixOfDubinerPolynomial(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            // Construct the Vandermonde matrix of DubinerPolynomials
            SharedNekMatrixPtr vandermonde(new NekMatrix<NekDouble>( GetVandermonde(x, y) ));
    
            return vandermonde;
        }

        NekVector<NekDouble> MakeTetWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z)
        {

            NekVector<NekDouble> g = MakeTetQuadratureSystem(x.GetRows());
          
            SharedNekMatrixPtr matV = MakeVmatrixOfTet(x, y, z);

            // Initialize the linear system solver
            LinearSystem matL(matV);


            // Deduce the quadrature weights
            NekVector<NekDouble> w = matL.SolveTranspose(g);
            
            return w;
        }
    
        NekVector<NekDouble> MakeQuadratureWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
    
            NekVector<NekDouble> g = MakeDubinerQuadratureSystem(x.GetRows());
            
            SharedNekMatrixPtr   matV =  MakeVmatrixOfDubinerPolynomial( x, y);

            // Initialize the linear system solver
            LinearSystem matL(matV);

            // Deduce the quadrature weights
            return matL.SolveTranspose(g);
    
        }

         NekMatrix<NekDouble> GetTetInterpolationMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                       const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi)
        {
            int nNodes = GetSize(x);
            int degree = GetTetDegree(nNodes);

            NekMatrix<NekDouble> matS = GetTetVandermonde(x, y, z); // Square 'short' matrix
            NekMatrix<NekDouble> matT = GetTetVandermonde(xi, yi, zi, degree); // Coefficient interpolation matrix (tall)

            NekMatrix<NekDouble> invertMatrix = Invert(matS);

            // Get the interpolation matrix
            return matT*invertMatrix;
        }
            
        NekMatrix<NekDouble> GetInterpolationMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi)
        {
            int nNodes = GetSize(x);
            int degree = GetDegree(nNodes);
    
            NekMatrix<NekDouble> matS = GetVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> matT = GetVandermonde(xi, yi, degree); // Coefficient interpolation matrix (tall)
    
            NekMatrix<NekDouble> invertMatrix = Invert(matS);

            // Get the interpolation matrix
            return matT*invertMatrix;
        }
    
        NekVector<NekDouble> LegendrePolyDerivative(int degree, const NekVector<NekDouble>& x)
        {
            int size = GetSize(x);
            NekVector<NekDouble> y(size), b3(size), a2(size);
            
            if(degree >= 3)
            {
                NekVector<NekDouble> b0(size, 0.0), b1(size, 1.0), b2 = 3.0 * x;
                NekVector<NekDouble> a0(size, 1.0), a1 = x;
    
                for(int n=3; n<=degree; ++n)
                {
                    a2 = ((2.0*n - 3.0)/(n - 1.0)) * Hadamard(x, a1) - (n - 2.0)/(n - 1.0) * a0;
                    a0 = a1;
                    a1 = a2;
    
                    b3 = (2.0*n - 1.0)/n * (Hadamard(b2, x) + a2) - (n - 1.0)/n * b1;
                    b1 = b2;
                    b2 = b3;
                }
                y = b3;
                
            }
            else if(degree == 2)
            {
                y = 3.0 * x;                
            }
            else if(degree == 1)
            {
                y = NekVector<NekDouble>(size, 1.0);
                
            }
            else
            {
                y = NekVector<NekDouble>(size, 0.0);
            }
            return y;
        }
    
    
        NekVector<NekDouble> DubinerPolyXDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {            
            // Get the coordinate transform
            int size = GetSize(y);
            NekVector<NekDouble> eta_1(size);
            NekVector<NekDouble> psi_x(x.GetRows());
            if(p>0)
            {
                    // Initialize the horizontal coordinate of the triangle (beta in Barycentric coordinates)
                for(int el=0; el<size; ++el)
                {
                    if(y[el] < 1.0 - 1e-16)
                    {
                        eta_1[el] = 2.0*(1.0 + x[el]) / (1.0 - y[el]) - 1.0;
                    }
                    else
                    {
                        eta_1[el] = -1.0; // When y is close to 1, then we have a removeable singularity
                    }
                }
                // Initialize the vertical coordinate of the triangle (gamma in Barycentric coordinates)
                NekVector<NekDouble> eta_2 = y;
    
                // Orthogonal Dubiner polynomial x-derivative
                int alpha = 2*p + 1; int beta = 0;
                NekVector<NekDouble> psi_a_x = LegendrePolyDerivative(p, eta_1);
                NekVector<NekDouble> upsilon = JacobiPoly(q, eta_2, alpha, beta);
                NekDouble w = sqrt((2.0*p + 1.0) * (p + q + 1.0) / 2.0); // Normalizing Orthonormal weight
    
                for(int i=0; i<size; ++i)
                {
                    NekDouble zeta = pow((1.0 -  eta_2[i])/2.0, (p-1.0));
                    NekDouble psi_b = zeta*upsilon[i];
                    psi_x[i] = w * psi_a_x[i] * psi_b;
                }
            }
            else
            {
                for(int i=0; i<int(x.GetRows()); ++i)
                {
                    psi_x[i] = 0.0;
                }
            }
    
            return psi_x;
        }

        NekVector<NekDouble> TetXDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z)
        {
            // Get the coordinate transform
            int size = GetSize(y);
            NekVector<NekDouble> eta_1(size), eta_2(size);
            NekVector<NekDouble> psi_x(x.GetRows());
            if(p > 0){
                // Initialize the horizontal coordinate of the Tetrahedral (beta in Barycentric coordinate)
                for(int el=0; el<size; ++el)
                {
                    if( y[el] < -z[el] - numeric_limits<NekDouble>::epsilon())
                    {
                        eta_1[el] = 2.0*(1.0 + x[el])/(-y[el]-z[el]) - 1.0;
                    } else
                    {
                        eta_1[el] = -1.0;
                    }
                }
                
                    // Initialize the  coordinate of the Tetrahedral (gamma in Barycentric coordinate)
                for(int el=0; el<size; ++el)
                {
                    if( z[el] < 1.0 - numeric_limits<NekDouble>::epsilon())
                    {
                        eta_2[el] = 2.0*(1.0 + y[el]) / (1.0 - z[el]) - 1.0;
                    }
                    else
                    {
                        eta_2[el] = -1.0;  // When z is close to 1, then we have a removeable singularity
                    }
                }
                // Initialize the vertical coordinate of the Tetrahedral (delta in Barycentric coordinate)
                NekVector<NekDouble> eta_3 = z;
                            
                // Orthogonal basis polynomial x-derivative
                int alpha = 2*p + 1; int beta = 0;
                int alpha_r = 2*p + 2*q + 2; int beta_r = 0;
                NekVector<NekDouble> psi_a_x    = LegendrePolyDerivative(p,eta_1);
                NekVector<NekDouble> psi_bpq  = JacobiPoly(q, eta_2, alpha, beta);
                NekVector<NekDouble> psi_cpqr = JacobiPoly(r, eta_3, alpha_r, beta_r);

                for(int el=0; el<size; ++el)
                {
                    NekDouble jacobi_b =  pow((1.0-eta_2[el])/2.0, p-1.0);
                    NekDouble jacobi_c =  pow((1.0-eta_3[el])/2.0,p+q-1.0);
                    NekDouble psi_b = jacobi_b * psi_bpq[el];
                    NekDouble psi_c = jacobi_c * psi_cpqr[el];
                    psi_x[el] = psi_a_x[el] * psi_b * psi_c;
                }
            }
            else
            {
                for(int el=0; el<int(x.GetRows()); ++el)
                {
                    psi_x[el] = 0.0;
                }
            }
            return psi_x;
        }

        NekMatrix<NekDouble> GetVandermondeForTetXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                             const NekVector<NekDouble>& z, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) * (degree + 3)/6;
            NekMatrix<NekDouble> matVx(rows, cols, 0.0);
            
            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p <= d; ++p)
                {
                    for(int q=0; q <= d - p; ++q, ++n)
                    {
                        int r = d - p - q;
                        NekVector<NekDouble> columnVector = TetXDerivative(p, q, r, x, y, z);

                        // Set n-th column of V to the TetrahedralBasis vector
                        for(int i=0; i<rows; ++i)
                        {
                            matVx(i,n) = columnVector[i];
                        }
                    }
                }
            }
            return matVx;
        }

        NekMatrix<NekDouble> GetVandermondeForTetXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z)
        {
            int rows = GetSize(x);
            int degree = GetTetDegree(rows);
            return GetVandermondeForTetXDerivative(x, y, z, degree);
        }

        Points<NekDouble>::MatrixSharedPtrType GetTetXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
        const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi)
        {
            int nNodes = GetSize(x);
            int degree = GetTetDegree(nNodes);
            NekMatrix<NekDouble> matS  = GetTetVandermonde(x, y, z); // Square 'short' matrix
            NekMatrix<NekDouble> matTx = GetVandermondeForTetXDerivative(xi, yi, zi, degree); // Tall matrix
    
            NekMatrix<NekDouble> invertMatrix = Invert(matS);
    
            // Get the Derivation matrix    
            return Points<NekDouble>::MatrixSharedPtrType( new NekMatrix<NekDouble>(matTx*invertMatrix) );
        }
    
        
        NekMatrix<NekDouble> GetVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2)/2;
            NekMatrix<NekDouble> matVx(rows, cols, 0.0);
            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p<=d; ++p, ++n)
                {
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPolyXDerivative(p, q, x, y);
    
                    // Set n-th column of Vx to the DubinerPolyXDerivative vector
                    for(int i=0; i<rows; ++i)
                    {
                        matVx(i, n) = columnVector[i];
                    }
                }
            }
            
            return matVx;
        }
    
            
        NekMatrix<NekDouble> GetVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetVandermondeForXDerivative(x, y, degree);
        }
            
        Points<NekDouble>::MatrixSharedPtrType GetXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi)
        {
            int nNodes = GetSize(x);
            int degree = GetDegree(nNodes);
            NekMatrix<NekDouble> matS  = GetVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> matTx = GetVandermondeForXDerivative(xi, yi, degree); // Tall matrix
    
            NekMatrix<NekDouble> invertMatrix = Invert(matS);
    
            // Get the Derivation matrix    
            return Points<NekDouble>::MatrixSharedPtrType( new NekMatrix<NekDouble>(matTx*invertMatrix) );
        }
    
    
        NekVector<NekDouble> JacobiPolyDerivative(int degree, const NekVector<NekDouble>& x, int alpha, int beta)
        {
            int size = GetSize(x);
            NekVector<NekDouble> y(size);
    
            if(degree == 0)
            {
                y = NekVector<NekDouble>(size, 0.0);
    
            }
             else
             {
                y = 0.5 * (alpha + beta + degree + 1) * JacobiPoly(degree - 1, x, alpha + 1, beta + 1);
            }
            return y;
        }



         NekVector<NekDouble> TetYDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z)
         {
            // Get the coordinate transform
            int size = GetSize(y);
            NekVector<NekDouble> eta_1(size), eta_2(size);
            NekVector<NekDouble> eta_1_dy(size), eta_2_dy(size);
    
            // Initialize the collapsed horizontal coordinate of the Tetrahedral
            for(int el=0; el<size; ++el)
            {
                if( y[el] < -z[el] - numeric_limits<NekDouble>::epsilon())
                {
                    eta_1[el]       = 2.0*(1.0 + x[el])/(-y[el]-z[el]) - 1.0;
                    eta_1_dy[el]    = 2.0*(1.0 + x[el])/((y[el]+z[el])*(y[el]+z[el]));
                }
                else
                {
                    eta_1[el]       = -1.0;
                    eta_1_dy[el]    =  0.0; // No change in the squeeze direction
                }
            }
    
            // Initialize the collapsed depth coordinate of the Tetrahedral
            for(int el=0; el<size; ++el)
            {
                if( z[el] < 1.0 - numeric_limits<NekDouble>::epsilon())
                {
                    eta_2[el]       = 2.0*(1.0 + y[el]) / (1.0 - z[el]) - 1.0;
                    eta_2_dy[el]    = 2.0/(1.0 - z[el]);
                }
                else
                {
                        // When z is close to 1, then we have a removeable singularity
                    eta_2[el]       = -1.0;
                    eta_2_dy[el]    =  0.0; // No change in the squeeze direction
                }
            }
            // Initialize the collapsed vertical coordinate of the Tetrahedral
            NekVector<NekDouble> eta_3 = z;
    
            // Orthogonal basis expansion polynomials and their y-derivatives for the tetrahedron
            // phi(vec xi) = psi_a(eta_1(xi)) * psi_b(eta_2(xi)) * psi_c(eta_3(xi))
            int alpha_b = 2*p + 1; int beta_b = 0;
            int alpha_c = 2*p + 2*q + 2; int beta_c = 0;
            NekVector<NekDouble> ones( size, 1.0 );
    
            NekVector<NekDouble> jacobi_b   = VectorPower( (ones - eta_2)/2.0, p );
            NekVector<NekDouble> jacobi_c   = VectorPower( (ones - eta_3)/2.0, p+q );
    
            NekVector<NekDouble> J_q        = JacobiPoly(q, eta_2, alpha_b, beta_b);
            NekVector<NekDouble> J_r        = JacobiPoly(r, eta_3, alpha_c, beta_c);
    
            NekVector<NekDouble> psi_a      = LegendrePoly(p, eta_1);
            NekVector<NekDouble> LD         = LegendrePolyDerivative(p, eta_1);
            NekVector<NekDouble> JD         = JacobiPolyDerivative(q, eta_2, alpha_b, beta_b);
            NekVector<NekDouble> psi_b      = Hadamard( J_q, jacobi_b );
            NekVector<NekDouble> psi_c      = Hadamard( J_r, jacobi_c );
    
            // Compute the partials wrt y (psi_c_dy = 0)
            NekVector<NekDouble> secondComponentOfPsi_b_dy( size, 0.0 );

            if(p > 0)
            {
                for(int i=0; i<size; ++i)
                {
                    NekDouble jacobi_b_dy = -p / 2.0 * pow( (1.0 - eta_2[i])/2.0, p-1.0 );
                    secondComponentOfPsi_b_dy[i] = J_q[i] * jacobi_b_dy;
                }
            }
             else
            {
                for(int i=0; i<size; ++i)
                {
                    secondComponentOfPsi_b_dy[i] = 0.0;
                }
            }

            NekVector<NekDouble> psi_a_dy   = Hadamard( LegendrePolyDerivative(p, eta_1), eta_1_dy );
            NekVector<NekDouble> psi_b_dy   = Hadamard( Hadamard( JacobiPolyDerivative(q, eta_2, alpha_b, beta_b), jacobi_b)
                                              + secondComponentOfPsi_b_dy, eta_2_dy );

            NekVector<NekDouble> psi_dy(size);
            for(int k=0; k<size; ++k)
            {
                psi_dy[k] = (psi_a_dy[k]*psi_b[k]*psi_c[k])  +  (psi_a[k]*psi_b_dy[k]*psi_c[k]);
            }

            // Fix singularity at z=1
            for(int k=0; k<size; ++k)
            {
                if( z[k] >= 1.0 - numeric_limits<NekDouble>::epsilon() )
                {
                    if( p + q > 0 )
                    {
                        psi_dy[k] = ((2.0*p+q+2.0)*JacobiPoly(q, eta_2[k], alpha_b, 1) - (p+q+2.0)*J_q[k]) *
                        psi_a[k] * pow( (1.0-eta_3[k])/2.0, p+q-1 ) / 2.0 * ((p+q+r+3.0)*J_r[k] - (2.0*p+2.0*q+r+3.0)*JacobiPoly(r, eta_3[k], alpha_c, 1));
                    }
                    else
                    {
                        psi_dy[k] = 0.0;
                    }
                }
            }

//             cout << "(p,q,r) = (" << p << ", " << q << ", " << r << ")" << endl;
// 
//             cout << "psi_a    = " << VectorToString(psi_a) << endl;
//             cout << "psi_b    = " << VectorToString(psi_b) << endl;
//             cout << "psi_c    = " << VectorToString(psi_c) << endl;
//             cout << "psi_a_dy = " << VectorToString(psi_a_dy) << endl;
//             cout << "psi_b_dy = " << VectorToString(psi_b_dy) << endl;
//             cout << "(psi_a_dy*psi_b*psi_c)          = " << VectorToString(Hadamard(Hadamard(psi_a_dy,psi_b), psi_c)) << endl;
//             cout << "(psi_a*jacobi_b*psi_b_dy*psi_c) = " << VectorToString(Hadamard(Hadamard(psi_a,jacobi_b),Hadamard(psi_b_dy,psi_c))) << endl;
//             cout << "secondComponentOfPsi_b_dy = " << VectorToString(secondComponentOfPsi_b_dy) << endl;
//             
//             cout << "psi_dy   = " << VectorToString(psi_dy) << endl;
// // //            cout << "jacobi_c = " << VectorToString(jacobi_c) << endl;
// // //             cout << "eta_3 = " << VectorToString(eta_3) << endl;
// // //             cout << "ones - eta_3 = " << VectorToString(ones - eta_3) << endl;
// // //             cout << "(ones - eta_3)/2.0 = " << VectorToString((ones - eta_3)/2.0) << endl;


            return psi_dy;
        }

        
        NekMatrix<NekDouble> GetVandermondeForTetYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                              const NekVector<NekDouble>& z, int degree){
            int rows = GetSize(y);
            int cols = (degree + 1) * (degree + 2) * (degree + 3)/6;
            NekMatrix<NekDouble> matVy(rows, cols, 0.0);
            
            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p <= d; ++p)
                {
                    for(int q=0; q <= d - p; ++q, ++n)
                    {
                        int r = d - p - q;
                        NekVector<NekDouble> columnVector = TetYDerivative(p, q, r, x, y, z);

                        // Set n-th column of V to the TetrahedralBasis vector
                        for(int i=0; i<rows; ++i)
                        {
                            matVy(i,n) = columnVector[i];
                        }
                    }
                }
            }
            return matVy;
        }

        NekMatrix<NekDouble> GetVandermondeForTetYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                             const NekVector<NekDouble>& z)
        {
            int rows = GetSize(y);
            int degree = GetTetDegree(rows);
            return GetVandermondeForTetYDerivative(x, y, z, degree);
        }

        Points<NekDouble>::MatrixSharedPtrType GetTetYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
            const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi)
        {
            // int nNodes = GetSize(y);
            
            NekMatrix<NekDouble> matS  = GetTetVandermonde(x, y, z); // Square 'short' matrix

            NekMatrix<NekDouble> matTy = GetVandermondeForTetYDerivative(xi, yi, zi); // Tall matrix

            NekMatrix<NekDouble> invertMatrix = Invert(matS);
            
            NekMatrix<NekDouble> makeDerivativeMatrix = matTy*invertMatrix;

            return Points<NekDouble>::MatrixSharedPtrType(new NekMatrix<NekDouble>(makeDerivativeMatrix));

        }
    
         NekVector<NekDouble> TetZDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z)
         {
            int size = GetSize(z);
            NekVector<NekDouble> eta_1(size), eta_2(size);
            NekVector<NekDouble> eta_1_dz(size), eta_2_dz(size);
            
             // Initialize the collapsed horizontal coordinate of the Tetrahedral 
            for(int el=0; el<size; ++el)
            {
                if( y[el] < -z[el] - numeric_limits<NekDouble>::epsilon())
                {
                    eta_1[el]       = 2.0*(1.0 + x[el]) / (-y[el]-z[el]) - 1.0;
                    eta_1_dz[el]    = 2.0*(1.0 + x[el]) / ((y[el]+z[el])*(y[el]+z[el]));
                }
                else
                {
                    eta_1[el]       = -1.0;
                    eta_1_dz[el]    =  0.0; // No change in the squeeze direction
                }
            }

            // Initialize the collapsed depth coordinate of the Tetrahedral
            for(int el=0; el<size; ++el)
            {
                if( z[el] < 1.0 - numeric_limits<NekDouble>::epsilon())
                {
                    eta_2[el]       = 2.0*(1.0 + y[el]) / (1.0 - z[el])  -  1.0;
                    eta_2_dz[el]    = 2.0*(1.0 + y[el]) / ((1.0 - z[el])*(1.0 - z[el]));
                }
                else
                {    // When z is close to 1, then we have a removeable singularity
                    eta_2[el]       = -1.0;
                    eta_2_dz[el]    =  0.0; // No change in the squeeze direction
                }
            }


            // Initialize the collapsed vertical coordinate of the Tetrahedral 
            NekVector<NekDouble> eta_3 = z;
            NekVector<NekDouble> eta_3_dz(size, 1.0);

            
            // Orthogonal basis expansion polynomials and their z-derivatives for the tetrahedron
            // phi(vec xi) = psi_a(eta_1(xi)) * psi_b(eta_2(xi)) * psi_c(eta_3(xi))
            int alpha_b = 2*p + 1; int beta_b = 0;
            int alpha_c = 2*p + 2*q + 2; int beta_c = 0;
            NekVector<NekDouble> ones( size, 1.0 );

            NekVector<NekDouble> jacobi_b   = VectorPower( (ones - eta_2)/2.0, p );
            NekVector<NekDouble> jacobi_c   = VectorPower( (ones - eta_3)/2.0, p+q );
            
            NekVector<NekDouble> J_q        = JacobiPoly(q, eta_2, alpha_b, beta_b);
            NekVector<NekDouble> J_r        = JacobiPoly(r, eta_3, alpha_c, beta_c);
            
            NekVector<NekDouble> psi_a      = LegendrePoly(p, eta_1);
            NekVector<NekDouble> psi_b      = Hadamard( J_q, jacobi_b );
            NekVector<NekDouble> psi_c      = Hadamard( J_r, jacobi_c );

                        
            // Compute the partials wrt y and z (psi_c_dy = 0)
            NekVector<NekDouble> secondComponentOfPsi_b_dz( size, 0.0 );            
            NekVector<NekDouble> secondComponentOfPsi_c_dz( size, 0.0 );
            if(p > 0)
            {
                for(int i=0; i<size; ++i)
                {
                    NekDouble jacobi_b_dz = -p / 2.0 * pow( (1.0 - eta_2[i])/2.0, p-1.0 );
                    secondComponentOfPsi_b_dz[i] = J_q[i] * jacobi_b_dz;
                }
            }
            
            if(p + q > 0)
            {
                for(int i=0; i<size; ++i)
                {
                    NekDouble jacobi_c_dz = -(p+q)/2.0*pow( (1.0 - eta_3[i])/2.0, p+q-1.0);
                    secondComponentOfPsi_c_dz[i] = J_r[i] * jacobi_c_dz;
                }
            }
            
            NekVector<NekDouble> psi_a_dz   = Hadamard( LegendrePolyDerivative(p, eta_1), eta_1_dz );
            NekVector<NekDouble> psi_b_dz   = Hadamard( Hadamard(JacobiPolyDerivative(q, eta_2, alpha_b, beta_b), jacobi_b)
                                              + secondComponentOfPsi_b_dz, eta_2_dz );
            NekVector<NekDouble> psi_c_dz   = Hadamard( Hadamard(JacobiPolyDerivative(r, eta_3, alpha_c, beta_c), jacobi_c)
                                              + secondComponentOfPsi_c_dz, eta_3_dz );

            NekVector<NekDouble> psi_dz(size);
            
            for(int k=0; k<size; ++k)
            {
               psi_dz[k] = psi_a_dz[k] * psi_b[k] * psi_c[k]  +  psi_a[k] * psi_b_dz[k] * psi_c[k]  +  psi_a[k] * psi_b[k] * psi_c_dz[k];
            }
            
            return psi_dz;                                                                    
        }

        NekMatrix<NekDouble> GetVandermondeForTetZDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                              const NekVector<NekDouble>& z, int degree){
            int rows = GetSize(z);
            int cols = (degree + 1) * (degree + 2) * (degree + 3)/6;
            NekMatrix<NekDouble> matVz(rows, cols, 0.0);
            
            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p <= d; ++p)
                {
                    for(int q=0; q <= d - p; ++q, ++n)
                    {
                        int r = d - p - q;
                        NekVector<NekDouble> columnVector = TetZDerivative(p, q, r, x, y, z);

                        // Set n-th column of V to the TetrahedralBasis vector
                        for(int i=0; i<rows; ++i)
                        {
                            matVz(i,n) = columnVector[i];
                        }
                    }
                }
            }
            return matVz;
        }

        NekMatrix<NekDouble> GetVandermondeForTetZDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                             const NekVector<NekDouble>& z)
        {
            int rows = GetSize(z);
            int degree = GetTetDegree(rows);
            return GetVandermondeForTetZDerivative(x, y, z, degree);
        }

       Points<NekDouble>::MatrixSharedPtrType GetTetZDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
            const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi)
        {
            int nNodes = GetSize(z);
            int degree = GetTetDegree(nNodes);


            NekMatrix<NekDouble> matS  = GetTetVandermonde(x, y, z); // Square 'short' matrix

            NekMatrix<NekDouble> matTz = GetVandermondeForTetZDerivative(xi, yi, zi, degree); // Tall matrix

            NekMatrix<NekDouble> invertMatrix = Invert(matS);

            NekMatrix<NekDouble> makeDerivativeMatrix = matTz*invertMatrix;

            Points<NekDouble>::MatrixSharedPtrType TetZDerivative;
            TetZDerivative = Points<NekDouble>::MatrixSharedPtrType(new  NekMatrix<NekDouble>(makeDerivativeMatrix));

           return TetZDerivative;

        }
        
        NekVector<NekDouble> DubinerPolyYDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            // Get the coordinate transform
            int size = GetSize(y);
            NekVector<NekDouble> eta_1(size);
    
            // Initialize the horizontal coordinate of the triangle (beta in Barycentric coordinates)
            for(int el=0; el<size; ++el)
            {
                if(y[el] < 1.0 - 1e-16)
                {
                    eta_1[el] = 2.0*(1.0 + x[el]) / (1.0 - y[el]) - 1.0;
                }
                 else
                 {
                    eta_1[el] = -1.0; // When y is close to 1, then we have a removeable singularity
                }
            }
            // Initialize the vertical coordinate of the triangle (gamma in Barycentric coordinates)
            NekVector<NekDouble> eta_2 = y;
    
            // Orthogonal Dubiner y-polynomial
            int alpha = 2*p + 1; int beta = 0;
            
            NekVector<NekDouble> psi_a = LegendrePoly(p, eta_1);
            NekVector<NekDouble> psi_a_y = LegendrePolyDerivative(p, eta_1);
            NekVector<NekDouble> psi_b = JacobiPoly(q, eta_2, alpha, beta);
            NekVector<NekDouble> psi_b_y = JacobiPolyDerivative(q, eta_2, alpha, beta);
            NekVector<NekDouble> secondComponentOf_psi_b(size);
            NekVector<NekDouble> psi_y(size);
            
            NekVector<NekDouble> first_part_derivative(size);

            if(p > 0)
            {
                for(int i=0; i<size; ++i)
                {
                     
                     first_part_derivative[i] = (1.0 + eta_1[i])/2.0 * psi_a_y[i] * pow((1.0 - eta_2[i])/2.0, p-1.0) * psi_b[i];
                     secondComponentOf_psi_b[i] = -p/2.0 * pow(((1.0 - eta_2[i]) / 2.0), p - 1.0) * psi_b[i];

                }                
            }
            else
            {
                for(int i=0; i<size; ++i)
                {
                    secondComponentOf_psi_b[i] = 0.0;
                    first_part_derivative[i] = 0.0;
                }
            }
            for(int k=0; k<size; ++k)
            {
                     psi_y[k] = first_part_derivative[k] + psi_a[k] * ( pow((1.0 -  eta_2[k])/2.0, p) * psi_b_y[k] + secondComponentOf_psi_b[k] );

            }
            NekDouble w = sqrt((2.0*p + 1.0)*(p + q + 1.0)/2.0); // Normalizing Orthonormal weight
    
            return w * psi_y;
        }


        
        NekMatrix<NekDouble> GetVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree)
        {
            int rows = GetSize(y);
            int cols = (degree + 1) * (degree + 2)/2;
            NekMatrix<NekDouble> matVy(rows, cols, 0.0);
            
            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p<=d; ++p, ++n)
                {
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPolyYDerivative(p, q, x, y);
    
                    for(int i=0; i<rows; ++i)
                    {
                        matVy(i, n) = columnVector[i];
                    }
                }
            }
            return matVy;
        }
    
        
        Points<NekDouble>::MatrixSharedPtrType GetYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi)
        {
            int nNodes = GetSize(y);
            int degree = GetDegree(nNodes);
    
            NekMatrix<NekDouble> matS  = GetVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> matTy = GetVandermondeForYDerivative(xi, yi, degree); // Tall matrix
    
            NekMatrix<NekDouble> invertMatrix = Invert(matS);


            // Get the Derivation matrix
            return Points<NekDouble>::MatrixSharedPtrType( new NekMatrix<NekDouble>(matTy*invertMatrix) );
        }
        
    
        NekMatrix<NekDouble> GetVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            int rows = GetSize(y);
            int degree = GetDegree(rows);
            return GetVandermondeForYDerivative(x, y, degree);
        }
    
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1)*(degree + 2)/2;
            NekMatrix<NekDouble>matV(rows,cols, 0.0);

            for(int d=0, n=0; d <= degree; ++d)
            {
                for(int p=0; p <= d; ++p, ++n)
                {
                    int q = d - p;
                    
                    for(int i=0; i<rows; ++i)
                    {
                            matV(i, n) = pow(x[i], p) * pow(y[i], q);
                    }
                }
            }
            return matV;
        }
        
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& z, int degree)
        {            
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) * (degree + 3)/6;

            NekMatrix<NekDouble> matV(rows, cols, 0.0);

            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p <= d; ++p)
                {
                    for(int q=0; q <= d - p; ++q, ++n)
                    {
                        int r = d - p - q;

                        // Set n-th column of V to the monomial vector
                        for(int i=0; i<rows; ++i)
                        {
                            matV(i,n) = pow(x[i],p) * pow(y[i],q) * pow(z[i],r);
                        }
                    }
                }
            }
            return matV;
        }
        
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z)
        {
            int rows = GetSize(x);
            int degree = GetTetDegree(rows);
            return GetMonomialVandermonde(x, y, z, degree);
        }
        
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetMonomialVandermonde(x, y, degree);
        }
        
        NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) / 2;
            NekMatrix<NekDouble> matVx(rows, cols, 0.0);

            for(int d=0, n=0; d <= degree; ++d)
            {
                for(int p=0; p <= d; ++p, ++n)
                {
                    int q = d - p;

                    if(p > 0)
                    {
                        for(int i=0; i<rows; ++i)
                        {
                            matVx(i, n) = p * pow(x[i], p-1.0) * pow(y[i],q);
                        }
                    }
                    else
                    {
                        for(int j=0; j<rows; ++j)
                        {
                            matVx(j, n) = 0.0;
                        }
                    }
                }
            }
            return matVx;
        }
        
        NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetXDerivativeOfMonomialVandermonde(x, y, degree);
        }

        NekMatrix<NekDouble> GetTetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
            NekMatrix<NekDouble> matVx(rows, cols, 0.0);
            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p <= d; ++p)
                {
                    for(int q=0; q <= d - p; ++q, ++n)
                    {
                        int r = d - p - q;
                        if(p > 0)
                        {
                            // Set n-th column of V to the monomial vector
                            for(int i=0; i<rows; ++i)
                            {
                                matVx(i,n) = p * pow(x[i],p-1.0) * pow(y[i],q) * pow(z[i],r);
                            }
                        }
                        else{
                            for(int j=0; j<rows; ++j)
                            {
                                matVx(j, n) = 0.0;
                            }
                        }
                    }
                }
            }
            return matVx;
        }
        
        NekMatrix<NekDouble> GetTetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z)
        {
            int rows = GetSize(x);
            int degree = GetTetDegree(rows);
            return GetTetXDerivativeOfMonomialVandermonde(x, y, z, degree);
        }
        
        NekMatrix<NekDouble> GetTetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
            NekMatrix<NekDouble> matVy(rows, cols, 0.0);
            for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p <= d; ++p)
                {
                    for(int q=0; q <= d - p; ++q, ++n)
                    {
                        int r = d - p - q;
                        if(q > 0)
                        {
                            // Set n-th column of V to the monomial vector
                            for(int i=0; i<rows; ++i)
                            {
                                matVy(i,n) = q * pow(x[i],p) * pow(y[i],q-1.0) * pow(z[i],r);
                            }
                        }
                        else
                        {
                            for(int j=0; j<rows; ++j)
                            {
                                matVy(j, n) = 0.0;
                            }
                        }
                    }
                }
            }
            return matVy;
        }
        
        NekMatrix<NekDouble> GetTetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z)
        {
            int rows = GetSize(x);
            int degree = GetTetDegree(rows);
            return GetTetYDerivativeOfMonomialVandermonde(x, y, z, degree);
        }
        
        NekMatrix<NekDouble> GetTetZDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
            NekMatrix<NekDouble> matVz(rows, cols, 0.0);
             for(int d=0, n=0; d<=degree; ++d)
            {
                for(int p=0; p <= d; ++p)
                {
                    for(int q=0; q <= d - p; ++q, ++n)
                    {
                        int r = d - p - q;
                        if(r > 0)
                        {
                            // Set n-th column of V to the monomial vector
                            for(int i=0; i<rows; ++i)
                            {
                                matVz(i,n) = r * pow(x[i],p) * pow(y[i],q) * pow(z[i],r-1.0);
                            }
                        }
                        else{
                            for(int j=0; j<rows; ++j)
                            {
                                matVz(j, n) = 0.0;
                            }
                        }
                    }
                }
            }
            return matVz;
        }
        
        NekMatrix<NekDouble> GetTetZDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z)
        {
            int rows = GetSize(x);
            int degree = GetTetDegree(rows);
            return GetTetZDerivativeOfMonomialVandermonde(x, y, z, degree);
        }
        
        NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree)
        {
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) / 2;
            NekMatrix<NekDouble> matVy(rows, cols, 0.0);
            
            for(int d=0, n=0; d <= degree; ++d)
            {
                for(int p=0; p <= d; ++p, ++n)
                {
                    int q = d - p;
                    if(q > 0)
                    {
                        for(int i=0; i<rows; ++i)
                        {
                            matVy(i, n) = q * pow(x[i], p) * pow(y[i], q-1.0);
                        }
                     }
                     else
                     {
                        for(int j=0; j<rows; ++j)
                        {
                            matVy(j, n) = 0.0;
                        }
                    }
                }
            }
            return matVy;
        }
    
        NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y)
        {
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetYDerivativeOfMonomialVandermonde(x, y, degree);
        }
        
        NekVector<NekDouble> GetIntegralOfMonomialVandermonde(int degree)
        {
            int cols = (degree + 1) * (degree + 2) / 2;
            NekVector<NekDouble>integralMVandermonde(cols, 0.0);
            
            for(int d=0, n=0; d <= degree; ++d)
            {
                for(int p=0; p <= d; ++p, ++n)
                {
                    int q = d - p;       
                    int sp = 1 - 2 * (p % 2);
                    int sq = 1 - 2 * (q % 2);
        
                    integralMVandermonde(n) = NekDouble(sp * (p+1) + sq * (q+1) + sp * sq * (p+q+2)) / ((p+1) * (q+1) * (p+q+2));
                }
            }
            return integralMVandermonde;
        }        
    }
}

