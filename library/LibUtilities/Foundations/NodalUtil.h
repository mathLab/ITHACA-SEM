///////////////////////////////////////////////////////////////////////////////
//
// File NodalUtil.h
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
// Description: 2D and 3D Nodal Triangle and Tetrahedron Utilities header file --
//              Basis function, Interpolation, Integral, Derivation, etc. 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALUTIL_H
#define NODALUTIL_H

#include <iosfwd>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.h>


namespace Nektar
{
    namespace LibUtilities
    {

        // /////////////////////////////////////
        // General matrix and vector stuff        
        template< typename T > NekVector<T> GetColumn(const NekMatrix<T> & matA, int n);
        NekMatrix<NekDouble> & SetColumn(NekMatrix<NekDouble> & matA, int n, const NekVector<NekDouble> & x);
        NekVector<NekDouble> GetE(int rows, int n);
        NekMatrix<NekDouble> Invert(const NekMatrix<NekDouble> & matA);
        NekMatrix<NekDouble> GetTranspose(const NekMatrix<NekDouble> & matA);
        int GetSize(const ConstArray<OneD, NekDouble> & x);
        int GetSize(const NekVector<NekDouble> & x);
        int GetDegree(int nBasisFunctions);
        NekDouble MakeRound(NekDouble);
        NekVector<NekDouble> ToVector( const ConstArray<OneD, NekDouble> & x );
        Array<OneD, NekDouble> ToArray( const NekVector<NekDouble> & x );
        NekVector<NekDouble> Hadamard( const NekVector<NekDouble> & x, const NekVector<NekDouble> & y );
        NekVector<NekDouble> VectorPower( const NekVector<NekDouble> & x, NekDouble p );
        std::string MatrixToString( const NekMatrix<NekDouble> & A, int precision = 2, double threshold = 1e12 );
        std::string VectorToString( const NekVector<NekDouble> & v, int precision = 2, double threshold = 1e12 );
        
        // ////////////////////////////////////////////////////////////////
        // Polynomials(Jacobi, Legendre, and Dubiner) and its derivations        
        NekVector<NekDouble> LegendrePoly(int degree, const NekVector<NekDouble>& x);
        NekVector<NekDouble> DubinerPoly(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> JacobiPoly(int degree, const NekVector<NekDouble>& x, NekDouble alpha, NekDouble beta);
        NekDouble JacobiPoly(int degree, NekDouble x, NekDouble alpha, NekDouble beta);
        NekVector<NekDouble> LegendrePolyDerivative(int degree, const NekVector<NekDouble>& x);
        NekVector<NekDouble> DubinerPolyXDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> DubinerPolyYDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> JacobiPolyDerivative(int degree, const NekVector<NekDouble>& x, int alpha, int beta);        

        // /////////////////////////////
        // Triangle stuff        
        NekVector<NekDouble> MakeDubinerQuadratureSystem(int nBasisFunctions);
        SharedNekMatrixPtr MakeVmatrixOfDubinerPolynomial(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> MakeQuadratureWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> GetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> GetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> GetVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> GetInterpolationMatrix(const NekVector<NekDouble>& x,  const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
        Points<NekDouble>::MatrixSharedPtrType GetYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
        Points<NekDouble>::MatrixSharedPtrType GetXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
        
        // /////////////////////////////
        // Tetrahedron stuff
        int GetTetDegree(int nBasisFunc);
        NekVector<NekDouble> TetrahedralBasis(int p, int q, int r, const NekVector<NekDouble>& x,
                                              const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetTetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z, int degree);
        NekMatrix<NekDouble> GetTetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        SharedNekMatrixPtr MakeVmatrixOfTet(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetTetInterpolationMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                       const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi);
        NekVector<NekDouble> MakeTetWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        NekVector<NekDouble> MakeTetQuadratureSystem(int nBasisFunctions);
        NekVector<NekDouble> TetXDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetVandermondeForTetXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetVandermondeForTetXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                             int degree);
        Points<NekDouble>::MatrixSharedPtrType GetTetXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                 const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi);

        NekVector<NekDouble> TetYDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetVandermondeForTetYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetVandermondeForTetYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                             int degree);
        Points<NekDouble>::MatrixSharedPtrType GetTetYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                 const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi);

        NekVector<NekDouble> TetZDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetVandermondeForTetZDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetVandermondeForTetZDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                             int degree);
        Points<NekDouble>::MatrixSharedPtrType GetTetZDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                 const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi);


        // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Monomial Vandermonde stuff for Triangle : Useful to test triangle(integration, interpolation, and derivation) 
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> GetIntegralOfMonomialVandermonde(int degree);
 
                 
        // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Monomial Vandermonde stuff for Tetrahedron : Useful to test tetrahedron(integration, interpolation, and derivation)
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& z, int degree);
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetTetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree);
        NekMatrix<NekDouble> GetTetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetTetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree);
        NekMatrix<NekDouble> GetTetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z);
        NekMatrix<NekDouble> GetTetZDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree);
        NekMatrix<NekDouble> GetTetZDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z);

    } // end of LibUtilities namespace
} // end of Nektar namespace

#endif //NODALUTIL_H
