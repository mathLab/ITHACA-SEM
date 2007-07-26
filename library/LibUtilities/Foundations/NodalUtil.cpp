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

        template< typename T > NekVector<T> GetColumn(const NekMatrix<T> & matA, int n) {
            NekVector<T> v(matA.GetRows());
            for( int i=0; i<matA.GetRows(); ++i ) {
                v[i] = matA(i,n);
            }
            return v;
        }
    
        NekMatrix<NekDouble> & SetColumn(NekMatrix<NekDouble> & matA, int n, const NekVector<NekDouble> & x) {
            for( int i=0; i<int(matA.GetRows()); ++i ) {
                matA(i,n) = x[i];
            }
            return matA;
        }
    
        // Standard basis vector in R^M
            NekVector<NekDouble> GetE(int rows, int n) {
            NekVector<NekDouble> e(rows, 0.0);
            e(n) = 1;
            return e;
        }
    
        NekMatrix<NekDouble> Invert(const NekMatrix<NekDouble> & matA) {
            int rows = matA.GetRows(), columns = matA.GetColumns();
            NekMatrix<NekDouble> matX(rows,columns);
    
            // The linear system solver
            LinearSystem<NekMatrix<NekDouble> > matL( SharedNekMatrixPtr(new NekMatrix<NekDouble>(matA)) );
    
            // Solve each column for the identity matrix
            for( int j=0; j<columns; ++j ) {
                SetColumn( matX, j, matL.Solve( GetE(rows,j) ) );
            }
            
            return matX;
        }
        
        NekMatrix<NekDouble> GetTranspose(const NekMatrix<NekDouble> & matA) {
            int rows = matA.GetRows(), columns = matA.GetColumns();
            NekMatrix<NekDouble> matX(rows,columns);
        
            for( int i=0; i<rows; ++i ) {
                for( int j=0; j<columns; ++j ) {
                    matX(j,i) = matA(i,j);
                }
            }
            return matX;
        }
    
        int GetSize(const ConstArray<OneD, NekDouble> & x) {
            return x.num_elements();
        }
        
        int GetSize(const NekVector<NekDouble> & x) {
            return x.GetRows();
        }
                
        NekVector<NekDouble> ToVector( const ConstArray<OneD, NekDouble> & x ) {
            return NekVector<NekDouble>( GetSize(x), x.data() );
        }
        
        Array<OneD, NekDouble> ToArray( const NekVector<NekDouble> & x ) {
            return Array<OneD, NekDouble>( GetSize(x), x.GetPtr() );
        }
    
        NekVector<NekDouble> Hadamard( const NekVector<NekDouble> & x, const NekVector<NekDouble> & y ) {
            int size = GetSize(x);
            NekVector<NekDouble> z(size);
    
            for( int i=0; i<size; ++i ) {
                z[i] = x[i] * y[i];
            }
            return z;
        }
    
        
        int GetDegree(int nBasisFunctions){
            return (-3 + int(sqrt(1.0 + 8*nBasisFunctions)))/2;
        }

		NekDouble MakeRound(NekDouble x) {
			return floor(x + 0.5);
		}
        
        NekVector<NekDouble> MakeDubinerQuadratureSystem(int nBasisFunctions){
            // Make the vector of integrals: note that each modal basis function integrates to zero except for the 0th degree
            NekVector<NekDouble> g(nBasisFunctions, 0.0);
            g(0) = sqrt(2.0);
    
            return g;
        }
    
    
        NekVector<NekDouble> JacobiPoly(int degree, const NekVector<NekDouble>& x, int alpha, int beta){
            int size = GetSize(x);
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
            int size = GetSize(x);
            NekVector<NekDouble> y(size);
            
            if(degree > 1){
                NekVector<NekDouble> a0(size, 1.0);
                NekVector<NekDouble> a1 = x;
                NekVector<NekDouble> a2(size);
    
                for(int i=2; i<=degree; ++i){
                    NekDouble b = NekDouble(2.0*i-1.0)/i;
                    NekDouble c = NekDouble(i-1.0)/i;
    
                    // multiply each elements in matrix
                    a2 = b * Hadamard(a1, x)  -  c * a0;
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
            int size = GetSize(x);
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
    
        NekMatrix<NekDouble> GetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) / 2;
            
            NekMatrix<NekDouble> matV(rows, cols, 0.0);
    
            for(int d=0, n=0; d<=degree; ++d){
                for(int p=0; p<=d; ++p, ++n){
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPoly(p, q, x, y);
    
                    // Set n-th column of V to the DubinerPoly vector
                    for(int i=0; i<rows; ++i){
                            matV(i,n) = columnVector[i];
                    }
                }
            }
    
            return matV;
        }
    
        NekMatrix<NekDouble> GetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetVandermonde( x, y, degree );
        }
    
    
        SharedNekMatrixPtr MakeVmatrixOfDubinerPolynomial(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            
            // Construct the Vandermonde matrix of DubinerPolynomials
            SharedNekMatrixPtr vandermonde(new NekMatrix<NekDouble>( GetVandermonde(x, y) ));
    
            return vandermonde;
        }
    
        NekVector<NekDouble> MakeQuadratureWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
    
            NekVector<NekDouble> g = MakeDubinerQuadratureSystem(x.GetRows());
            
            SharedNekMatrixPtr   matV = MakeVmatrixOfDubinerPolynomial(x, y);
    
            // Initialize the linear system solver
            LinearSystem<NekMatrix<NekDouble> > matL(matV);
    
            // Deduce the quadrature weights
            return matL.SolveTranspose(g);
    
        }
    
        NekMatrix<NekDouble> GetInterpolationMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi){
            int nNodes = GetSize(x);
            int degree = GetDegree(nNodes);
    
            NekMatrix<NekDouble> matS = GetVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> matT = GetVandermonde(xi, yi, degree); // Coefficient interpolation matrix (tall)
    
            NekMatrix<NekDouble> invertMatrix = Invert(matS);
            // Get the interpolation matrix
            return matT*invertMatrix;
        }
    
        NekVector<NekDouble> LegendrePolyDerivative(int degree, const NekVector<NekDouble>& x){
            int size = GetSize(x);
            NekVector<NekDouble> y(size), b3(size), a2(size);
            
            if(degree >= 3) {
                NekVector<NekDouble> b0(size, 0.0), b1(size, 1.0), b2 = 3.0 * x;
                NekVector<NekDouble> a0(size, 1.0), a1 = x;
    
                for(int n=3; n<=degree; ++n){
                    a2 = ((2.0*n - 3.0)/(n - 1.0)) * Hadamard(x, a1) - (n - 2.0)/(n - 1.0) * a0;
                    a0 = a1;
                    a1 = a2;
    
                    b3 = (2.0*n - 1.0)/n * (Hadamard(b2, x) + a2) - (n - 1.0)/n * b1;
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
            int size = GetSize(y);
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
        
        NekMatrix<NekDouble> GetVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2)/2;
            NekMatrix<NekDouble> matVx(rows, cols, 0.0);
            for(int d=0, n=0; d<=degree; ++d){
                for(int p=0; p<=d; ++p, ++n){
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPolyXDerivative(p, q, x, y);
    
                    // Set n-th column of Vx to the DubinerPolyXDerivative vector
                    for(int i=0; i<rows; ++i){
                        matVx(i, n) = columnVector[i];
                    }
                }
            }
            // cout << "Vx = \n"  << Vx << endl;
            return matVx;
        }
    
            
        NekMatrix<NekDouble> GetVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetVandermondeForXDerivative(x, y, degree);
        }
            
        Points<NekDouble>::MatrixSharedPtrType GetXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi){
            int nNodes = GetSize(x);
            int degree = GetDegree(nNodes);
            NekMatrix<NekDouble> matS  = GetVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> matTx = GetVandermondeForXDerivative(xi, yi, degree); // Tall matrix
    
            NekMatrix<NekDouble> invertMatrix = Invert(matS);
    
            // Get the Derivation matrix
    
            return Points<NekDouble>::MatrixSharedPtrType( new NekMatrix<NekDouble>(matTx*invertMatrix) );
        }
    
    
        NekVector<NekDouble> JacobiPolyDerivative(int degree, const NekVector<NekDouble>& x, int alpha, int beta){
            int size = GetSize(x);
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
            int size = GetSize(y);
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
            
            NekVector<NekDouble> vecPq = JacobiPoly(q, eta_2, alpha, beta);
            NekVector<NekDouble> upsilon = LegendrePolyDerivative(p, eta_1);
            NekVector<NekDouble> zeta = LegendrePoly(p, eta_1);
            NekVector<NekDouble> vecJp = JacobiPolyDerivative(q, eta_2, alpha, beta);
            NekVector<NekDouble> vecPq2(x.GetRows());
            NekVector<NekDouble> vecLy(x.GetRows());
            
            if(p > 0){
                for(int i=0; i<size; ++i){
                    vecPq2[i] = p/2.0 * pow(((1.0 - eta_2[i]) / 2.0), p - 1.0) * vecPq[i];
                    vecLy[i] = (1.0 + eta_1[i])/2.0 * upsilon[i] * pow((1.0 - eta_2[i])/2.0, p-1.0) * vecPq[i];
                }
                
            } else {
                for(int i=0; i<size; ++i){
                    vecPq2[i] = 0.0;
                    vecLy[i] = 0.0;
                }
            }
            for(int k=0; k<size; ++k){
                psi_y[k] = vecLy[k] + zeta[k] * (pow((1.0 -  eta_2[k])/2.0, p) * vecJp[k] - vecPq2[k]);
            }
            NekDouble w = sqrt((2.0*p + 1.0)*(p + q + 1.0)/2.0); // Normalizing Orthonormal weight
    
            return w * psi_y;
        }
        
        NekMatrix<NekDouble> GetVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int rows = GetSize(y);
            int cols = (degree + 1) * (degree + 2)/2;
            NekMatrix<NekDouble> matVy(rows, cols, 0.0);
            
            for(int d=0, n=0; d<=degree; ++d){
                for(int p=0; p<=d; ++p, ++n){
                    int q = d - p;
                    NekVector<NekDouble> columnVector = DubinerPolyYDerivative(p, q, x, y);
    
                    for(int i=0; i<rows; ++i){
                        matVy(i, n) = columnVector[i];
                    }
                }
            }
            return matVy;
        }
    
        
        Points<NekDouble>::MatrixSharedPtrType GetYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi){
            int nNodes = GetSize(y);
            int degree = GetDegree(nNodes);
    
            NekMatrix<NekDouble> matS  = GetVandermonde(x, y); // Square 'short' matrix
            NekMatrix<NekDouble> matTy = GetVandermondeForYDerivative(xi, yi, degree); // Tall matrix
    
            NekMatrix<NekDouble> invertMatrix = Invert(matS);
    
            // Get the Derivation matrix
            return Points<NekDouble>::MatrixSharedPtrType( new NekMatrix<NekDouble>(matTy*invertMatrix) );
        }
    
    
    
        NekMatrix<NekDouble> GetVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int rows = GetSize(y);
            int degree = GetDegree(rows);
            return GetVandermondeForYDerivative(x, y, degree);
        }
    
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int rows = GetSize(x);
            int cols = (degree + 1)*(degree + 2)/2;
            NekMatrix<NekDouble>matV(rows,cols, 0.0);
            for(int d=0, n=0; d <= degree; ++d){
                for(int p=0; p <= d; ++p, ++n){
                    int q = d - p;
                    
                    for(int i=0; i<rows; ++i){
                            matV(i, n) = pow(x[i], p) * pow(y[i], q);
                    }
                }
            }
            return matV;
        }
        
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetMonomialVandermonde(x, y, degree);
        }
        
        NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) / 2;
            NekMatrix<NekDouble> matVx(rows, cols, 0.0);
            for(int d=0, n=0; d <= degree; ++d){
                for(int p=0; p <= d; ++p, ++n){
                    int q = d - p;
                    if(p > 0){
                        for(int i=0; i<rows; ++i){
                            matVx(i, n) = p * pow(x[i], p-1.0) * pow(y[i],q);
                        }
                    } else {
                        for(int j=0; j<rows; ++j){
                            matVx(j, n) = 0.0;
                        }
                    }
                }
            }
            return matVx;
        }
        
        NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetXDerivativeOfMonomialVandermonde(x, y, degree);
        }
    
        NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree){
            int rows = GetSize(x);
            int cols = (degree + 1) * (degree + 2) / 2;
            NekMatrix<NekDouble> matVy(rows, cols, 0.0);
            for(int d=0, n=0; d <= degree; ++d){
                for(int p=0; p <= d; ++p, ++n){
                    int q = d - p;
                    if(q > 0){
                        for(int i=0; i<rows; ++i){
                            matVy(i, n) = q * pow(x[i], p) * pow(y[i], q-1.0);
                        }
                        }else {
                        for(int j=0; j<rows; ++j){
                            matVy(j, n) = 0.0;
                        }
                    }
                }
            }
            return matVy;
        }
    
        NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y){
            int rows = GetSize(x);
            int degree = GetDegree(rows);
            return GetYDerivativeOfMonomialVandermonde(x, y, degree);
        }
        
        NekVector<NekDouble> GetIntegralOfMonomialVandermonde(int degree){
            int cols = (degree + 1) * (degree + 2) / 2;
            NekVector<NekDouble>integralMVandermonde(cols, 0.0);
            for(int d=0, n=0; d <= degree; ++d){
                for(int p=0; p <= d; ++p, ++n){
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

