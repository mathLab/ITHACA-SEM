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

#include <iostream>
#include <algorithm>
#include <limits>
#include <math.h>

#include <LibUtilities/Foundations/NodalUtil.h>

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar {
    namespace LibUtilities {

        template< typename T > NekVector<T> getColumn(const NekMatrix<T> & A, int n) {
            NekVector<T> v(A.GetRows());
            for( int i=0; i<A.GetRows(); ++i ) {
                v[i] = A(i,n);
            }
            return v;
        }
    
        NekMatrix<NekDouble> & setColumn(NekMatrix<NekDouble> & A, int n, const NekVector<NekDouble> & x) {
            for( int i=0; i<int(A.GetRows()); ++i ) {
                A(i,n) = x[i];
            }
            return A;
        }
    
        // Standard basis vector in R^M
            NekVector<NekDouble> getE(int M, int n) {
            NekVector<NekDouble> e(M, 0.0);
            e(n) = 1;
            return e;
        }
    
        NekMatrix<NekDouble> invert(const NekMatrix<NekDouble> & A) {
            int M = A.GetRows(), N = A.GetColumns();
            NekMatrix<NekDouble> X(M,N);
    
            // The linear system solver
            LinearSystem<NekMatrix<NekDouble> > L( SharedNekMatrixPtr(new NekMatrix<NekDouble>(A)) );
    
            // Solve each column for the identity matrix
            for( int j=0; j<N; ++j ) {
                setColumn( X, j, L.Solve( getE(M,j) ) );
            }
            
            return X;
        }
        
        NekMatrix<NekDouble> getTranspose(const NekMatrix<NekDouble> & A) {
            int M = A.GetRows(), N = A.GetColumns();
            NekMatrix<NekDouble> X(M,N);
        
            for( int i=0; i<M; ++i ) {
                for( int j=0; j<N; ++j ) {
                    X(j,i) = A(i,j);
                }
            }
            return X;
        }
    
        int getSize(const ConstArray<OneD, NekDouble> & x) {
            return x.num_elements();
        }
        
        int getSize(const NekVector<NekDouble> & x) {
            return x.GetRows();
        }
                
        NekVector<NekDouble> toVector( const ConstArray<OneD, NekDouble> & x ) {
            return NekVector<NekDouble>( getSize(x), x.data() );
        }
        
        Array<OneD, NekDouble> toArray( const NekVector<NekDouble> & x ) {
            return Array<OneD, NekDouble>( getSize(x), x.GetPtr() );
        }
    
        NekVector<NekDouble> hadamard( const NekVector<NekDouble> & x, const NekVector<NekDouble> & y ) {
            int N = getSize(x);
            NekVector<NekDouble> z(N);
    
            for( int i=0; i<N; ++i ) {
                z[i] = x[i] * y[i];
            }
            return z;
        }
    
        
        int getDegree(int nBasisFunctions){
            return (-3 + int(sqrt(1.0 + 8*nBasisFunctions)))/2;
        }

		NekDouble round(NekDouble x) {
			return floor(x + 0.5);
		}
        
        NekVector<NekDouble> makeDubinerQuadratureSystem(int nBasisFunctions){
            // Make the vector of integrals: note that each modal basis function integrates to zero except for the 0th degree
            NekVector<NekDouble> g(nBasisFunctions, 0.0);
            g(0) = sqrt(2.0);
    
            return g;
        }
    
    
        NekVector<NekDouble> JacobiPoly(int degree, const NekVector<NekDouble>& x, int alpha, int beta){
            int size = getSize(x);
            NekVector<NekDouble> y(size);
            
            if(degree == 0){
                // Set y to ones
                y = NekVector<NekDouble>(size, 1.0);
    
            } else if (degree == 1) {
            
                for(int el=0; el<size; ++el){
                    y[el] = 0.5*(alpha - beta + (alpha + beta + 2.0) * x[el]);
                }
                
            } else if (degree > 1) {
            
                NekDouble degm1 = degree - 1.0;
                NekDouble tmp = 2.0 * degm1 + alpha + beta;
                NekDouble a1 = 2.0 * (degm1 + 1.0) * (degm1 + alpha + beta + 1.0) * tmp;
                NekDouble a2 = (tmp + 1.0) * (alpha * alpha - beta * beta);
                NekDouble a3 = tmp * (tmp + 1.0) * (tmp + 2.0);
                NekDouble a4 = 2.0 * (degm1 + alpha) * (degm1 + beta) * (tmp + 2.0);
    
                NekVector<NekDouble> z1 = JacobiPoly(degree-1, x, alpha, beta);
                NekVector<NekDouble> z2 = JacobiPoly(degree-2, x, alpha, beta);
                for (int el=0; el<size; ++el) {
                    y[el] = ((a2 + a3 * x[el]) * z1[el] - a4 * z2[el])/a1;
                }
                
            } else {
                cerr << "Bad degree" << endl;
            }
    
            return y;
        }
    
    
        NekVector<NekDouble> LegendrePoly(int degree, const NekVector<NekDouble>& x){
            int size = getSize(x);
            NekVector<NekDouble> y(size);
            
            if(degree > 1){
                NekVector<NekDouble> a0(size, 1.0);
                NekVector<NekDouble> a1 = x;
                NekVector<NekDouble> a2(size);
    
                for(int i=2; i<=degree; ++i){
                    NekDouble b = NekDouble(2.0*i-1.0)/i;
                    NekDouble c = NekDouble(i-1.0)/i;
    
                    // multiply each elements in matrix
                    a2 = b * hadamard(a1, x)  -  c * a0;
                    a0 = a1;
                    a1 = a2;
                }
                y = a2;
    
            } else if( degree == 1 ) {
                y = x;
                
            } else {
                y = NekVector<NekDouble>(size, 1.0);
            }
            return y;
        }
        
        NekVector<NekDouble> DubinerPoly(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            // Get the coordinate transform
            int size = getSize(x);
            NekVector<NekDouble> eta_1(size);
    
            // Initialize the horizontal coordinate of the triangle (beta in Barycentric coordinates)
            for(int el=0; el<size; ++el){
                if( y[el] < 1.0 - 1e-16 ){
                    eta_1[el] = 2.0*(1.0 + x[el]) / (1.0 - y[el]) - 1.0;
                } else {
                    eta_1[el] = -1.0; // When y is close to 1, then we have a removeable singularity
                }
            }
    
            // Initialize the vertical coordinate of the triangle (gamma in Barycentric coordinates)
            NekVector<NekDouble> eta_2 = y;
    
            // Orthogonal Dubiner polynomial
            int alpha = 2 * p + 1; int beta = 0;
            NekVector<NekDouble> psi(size);
            NekVector<NekDouble> psi_a = LegendrePoly(p, eta_1);
            NekVector<NekDouble> upsilon = JacobiPoly(q, eta_2, alpha, beta);
            NekDouble w = sqrt((2.0 * p + 1.0) * (p + q + 1.0) / 2.0); // Normalizing Orthonormal weight
    
            for(int el=0; el<size; ++el){
                NekDouble zeta   = pow((1.0 - eta_2[el])/2.0, p);
                NekDouble psi_b  = zeta * upsilon[el];
                psi[el]          = w * psi_a[el] * psi_b;
            }
            return psi;
        }
    
        NekMatrix<NekDouble> getVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int M = getSize(x);
            int N = (degree + 1) * (degree + 2) / 2;
            
            NekMatrix<NekDouble> V(M, N, 0.0);
    
            for(int d=0, n=0; d<=degree; ++d){
                for(int p=0; p<=d; ++p, ++n){
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPoly(p, q, x, y);
    
                    // Set n-th column of V to the DubinerPoly vector
                    for(int i=0; i<M; ++i){
                            V(i,n) = columnVector[i];
                    }
                }
            }
    
            return V;
        }
    
        NekMatrix<NekDouble> getVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int M = getSize(x);
            int degree = getDegree(M);
            return getVandermonde( x, y, degree );
        }
    
    
        SharedNekMatrixPtr makeVmatrixOfDubinerPolynomial(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            
            // Construct the Vandermonde matrix of DubinerPolynomials
            SharedNekMatrixPtr V(new NekMatrix<NekDouble>( getVandermonde(x, y) ));
    
            return V;
        }
    
        NekVector<NekDouble> makeQuadratureWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
    
            NekVector<NekDouble> g = makeDubinerQuadratureSystem(x.GetRows());
            
            SharedNekMatrixPtr   V = makeVmatrixOfDubinerPolynomial(x, y);
    
            // Initialize the linear system solver
            LinearSystem<NekMatrix<NekDouble> > L(V);
    
            // Deduce the quadrature weights
            return L.SolveTranspose(g);
    
        }
    
        NekMatrix<NekDouble> getInterpolationMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi){
            int nNodes = getSize(x);
            int degree = getDegree(nNodes);
    
            NekMatrix<NekDouble> S = getVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> T = getVandermonde(xi, yi, degree); // Coefficient interpolation matrix (tall)
    
            NekMatrix<NekDouble> invertMatrix = invert(S);
            // Get the interpolation matrix
            return T*invertMatrix;
        }
    
        NekVector<NekDouble> LegendrePolyDerivative(int degree, const NekVector<NekDouble>& x){
            int size = getSize(x);
            NekVector<NekDouble> y(size), b3(size), a2(size);
            
            if(degree >= 3) {
                NekVector<NekDouble> b0(size, 0.0), b1(size, 1.0), b2 = 3.0 * x;
                NekVector<NekDouble> a0(size, 1.0), a1 = x;
    
                for(int n=3; n<=degree; ++n){
                    a2 = ((2.0*n - 3.0)/(n - 1.0)) * hadamard(x, a1) - (n - 2.0)/(n - 1.0) * a0;
                    a0 = a1;
                    a1 = a2;
    
                    b3 = (2.0*n - 1.0)/n * (hadamard(b2, x) + a2) - (n - 1.0)/n * b1;
                    b1 = b2;
                    b2 = b3;
                }
                y = b3;
                
            } else if(degree == 2){
                y = 3.0 * x;
                
            } else if(degree == 1){
                y = NekVector<NekDouble>(size, 1.0);
                
            } else {
                y = NekVector<NekDouble>(size, 0.0);
            }
            return y;
        }
    
    
        NekVector<NekDouble> DubinerPolyXDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            
            // Get the coordinate transform
            int size = getSize(y);
            NekVector<NekDouble> eta_1(size);
            NekVector<NekDouble> psi_x(x.GetRows());
            if(p>0){
                    // Initialize the horizontal coordinate of the triangle (beta in Barycentric coordinates)
                for(int el=0; el<size; ++el){
                    if(y[el] < 1.0 - 1e-16){
                        eta_1[el] = 2.0*(1.0 + x[el]) / (1.0 - y[el]) - 1.0;
                    } else {
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
    
                for(int i=0; i<size; ++i){
                    NekDouble zeta = pow((1.0 -  eta_2[i])/2.0, (p-1.0));
                    NekDouble psi_b = zeta*upsilon[i];
                    psi_x[i] = w * psi_a_x[i] * psi_b;
                }
            } else {
                for(int i=0; i<int(x.GetRows()); ++i){
                    psi_x[i] = 0.0;
                }
            }
    
            return psi_x;
        }
        
        NekMatrix<NekDouble> getVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int M = getSize(x);
            int N = (degree + 1) * (degree + 2)/2;
            NekMatrix<NekDouble> Vx(M, N, 0.0);
            for(int d=0, n=0; d<=degree; ++d){
                for(int p=0; p<=d; ++p, ++n){
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPolyXDerivative(p, q, x, y);
    
                    // Set n-th column of Vx to the DubinerPolyXDerivative vector
                    for(int i=0; i<M; ++i){
                        Vx(i, n) = columnVector[i];
                    }
                }
            }
            // cout << "Vx = \n"  << Vx << endl;
            return Vx;
        }
    
            
        NekMatrix<NekDouble> getVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int M = getSize(x);
            int degree = getDegree(M);
            return getVandermondeForXDerivative(x, y, degree);
        }
            
        Points<NekDouble>::MatrixSharedPtrType getXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi){
            int nNodes = getSize(x);
            int degree = getDegree(nNodes);
            NekMatrix<NekDouble> S  = getVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> Tx = getVandermondeForXDerivative(xi, yi, degree); // Tall matrix
    
            NekMatrix<NekDouble> invertMatrix = invert(S);
    
            // Get the Derivation matrix
    
            return Points<NekDouble>::MatrixSharedPtrType( new NekMatrix<NekDouble>(Tx*invertMatrix) );
        }
    
    
        NekVector<NekDouble> JacobiPolyDerivative(int degree, const NekVector<NekDouble>& x, int alpha, int beta){
            int size = getSize(x);
            NekVector<NekDouble> y(size);
    
            if(degree == 0){
                y = NekVector<NekDouble>(size, 0.0);
    
            } else {
                y = 0.5 * (alpha + beta + degree + 1) * JacobiPoly(degree - 1, x, alpha + 1, beta + 1);
            }
            return y;
        }
    
    
        NekVector<NekDouble> DubinerPolyYDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            // Get the coordinate transform
            int size = getSize(y);
            NekVector<NekDouble> eta_1(size), psi_y(size);
    
            // Initialize the horizontal coordinate of the triangle (beta in Barycentric coordinates)
            for(int el=0; el<size; ++el){
                if(y[el] < 1.0 - 1e-16){
                    eta_1[el] = 2.0*(1.0 + x[el]) / (1.0 - y[el]) - 1.0;
                } else {
                    eta_1[el] = -1.0; // When y is close to 1, then we have a removeable singularity
                }
            }
            // Initialize the vertical coordinate of the triangle (gamma in Barycentric coordinates)
            NekVector<NekDouble> eta_2 = y;
    
            // Orthogonal Dubiner y-polynomial
            int alpha = 2*p + 1; int beta = 0;
            
            NekVector<NekDouble> Pq = JacobiPoly(q, eta_2, alpha, beta);
            NekVector<NekDouble> upsilon = LegendrePolyDerivative(p, eta_1);
            NekVector<NekDouble> zeta = LegendrePoly(p, eta_1);
            NekVector<NekDouble> Jp = JacobiPolyDerivative(q, eta_2, alpha, beta);
            NekVector<NekDouble> Pq2(x.GetRows());
            NekVector<NekDouble> Ly(x.GetRows());
            
            if(p > 0){
                for(int i=0; i<size; ++i){
                    Pq2[i] = p/2.0 * pow(((1.0 - eta_2[i]) / 2.0), p - 1.0) * Pq[i];
                    Ly[i] = (1.0 + eta_1[i])/2.0 * upsilon[i] * pow((1.0 - eta_2[i])/2.0, p-1.0) * Pq[i];
                }
                
            } else {
                for(int i=0; i<size; ++i){
                    Pq2[i] = 0.0;
                    Ly[i] = 0.0;
                }
            }
            for(int k=0; k<size; ++k){
                psi_y[k] = Ly[k] + zeta[k] * (pow((1.0 -  eta_2[k])/2.0, p) * Jp[k] - Pq2[k]);
            }
            NekDouble w = sqrt((2.0*p + 1.0)*(p + q + 1.0)/2.0); // Normalizing Orthonormal weight
    
            return w * psi_y;
        }
        
        NekMatrix<NekDouble> getVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int M = getSize(y);
            int N = (degree + 1) * (degree + 2)/2;
            NekMatrix<NekDouble> Vy(M, N, 0.0);
            
            for(int d=0, n=0; d<=degree; ++d){
                for(int p=0; p<=d; ++p, ++n){
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPolyYDerivative(p, q, x, y);
    
                    for(int i=0; i<M; ++i){
                        Vy(i, n) = columnVector[i];
                    }
                }
            }
            return Vy;
        }
    
        
        Points<NekDouble>::MatrixSharedPtrType getYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi){
            int nNodes = getSize(y);
            int degree = getDegree(nNodes);
    
            NekMatrix<NekDouble> S  = getVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> Ty = getVandermondeForYDerivative(xi, yi, degree); // Tall matrix
    
            NekMatrix<NekDouble> invertMatrix = invert(S);
    
            // Get the Derivation matrix
            return Points<NekDouble>::MatrixSharedPtrType( new NekMatrix<NekDouble>(Ty*invertMatrix) );
        }
    
    
    
        NekMatrix<NekDouble> getVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int M = getSize(y);
            int degree = getDegree(M);
            return getVandermondeForYDerivative(x, y, degree);
        }
    
        NekMatrix<NekDouble> getMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int M = getSize(x);
            int N = (degree + 1)*(degree + 2)/2;
            NekMatrix<NekDouble>V(M,N, 0.0);
            for(int d=0, n=0; d <= degree; ++d){
                for(int p=0; p <= d; ++p, ++n){
                    int q = d - p;
                    
                    for(int i=0; i<M; ++i){
                            V(i, n) = pow(x[i], p) * pow(y[i], q);
                    }
                }
            }
            return V;
        }
        
        NekMatrix<NekDouble> getMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int M = getSize(x);
            int degree = getDegree(M);
            return getMonomialVandermonde(x, y, degree);
        }
        
        NekMatrix<NekDouble> getXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int M = getSize(x);
            int N = (degree + 1) * (degree + 2) / 2;
            NekMatrix<NekDouble> Vx(M, N, 0.0);
            for(int d=0, n=0; d <= degree; ++d){
                for(int p=0; p <= d; ++p, ++n){
                    int q = d - p;
                    if(p > 0){
                        for(int i=0; i<M; ++i){
                            Vx(i, n) = p * pow(x[i], p-1.0) * pow(y[i],q);
                        }
                    } else {
                        for(int j=0; j<M; ++j){
                            Vx(j, n) = 0.0;
                        }
                    }
                }
            }
            return Vx;
        }
        
        NekMatrix<NekDouble> getXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int M = getSize(x);
            int degree = getDegree(M);
            return getXDerivativeOfMonomialVandermonde(x, y, degree);
        }
    
        NekMatrix<NekDouble> getYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int M = getSize(x);
            int N = (degree + 1) * (degree + 2) / 2;
            NekMatrix<NekDouble> Vy(M, N, 0.0);
            for(int d=0, n=0; d <= degree; ++d){
                for(int p=0; p <= d; ++p, ++n){
                    int q = d - p;
                    if(q > 0){
                        for(int i=0; i<M; ++i){
                            Vy(i, n) = q * pow(x[i], p) * pow(y[i], q-1.0);
                        }
                        }else {
                        for(int j=0; j<M; ++j){
                            Vy(j, n) = 0.0;
                        }
                    }
                }
            }
            return Vy;
        }
    
        NekMatrix<NekDouble> getYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int M = getSize(x);
            int degree = getDegree(M);
            return getYDerivativeOfMonomialVandermonde(x, y, degree);
        }
        
        NekVector<NekDouble> getIntegralOfMonomialVandermonde(int degree){
            int N = (degree + 1) * (degree + 2) / 2;
            NekVector<NekDouble>VI(N, 0.0);
            for(int d=0, n=0; d <= degree; ++d){
                for(int p=0; p <= d; ++p, ++n){
                    int q = d - p;
        
                    int sp = 1 - 2 * (p % 2);
                    int sq = 1 - 2 * (q % 2);
        
                    VI(n) = NekDouble(sp * (p+1) + sq * (q+1) + sp * sq * (p+q+2)) / ((p+1) * (q+1) * (p+q+2));
                }
            }
            return VI;
        }        
    }
}

