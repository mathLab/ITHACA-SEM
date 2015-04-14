///////////////////////////////////////////////////////////////////////////////
//
// File StdExpUtil.cpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <iosfwd>

#include <StdRegions/StdExpUtil.h>

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>


namespace Nektar
{
    namespace StdRegions
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
            LinearSystem<NekMatrix<NekDouble> > matL( SharedNekMatrixPtr(new NekMatrix<NekDouble>(matA)) );
    
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
    
        int GetSize(const ConstArray<OneD, NekDouble> & x)
        {
            return x.num_elements();
        }
        
        int GetSize(const NekVector<NekDouble> & x)
        {
            return x.GetRows();
        }
                
        NekVector<NekDouble> ToVector( const ConstArray<OneD, NekDouble> & x )
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
        
        
    
    } // end of namespace Stdregions
}  // end of namespace Nektar
