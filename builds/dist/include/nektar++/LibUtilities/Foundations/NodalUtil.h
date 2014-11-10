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

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/Foundations/Points.h>


//#include <LibUtilities/BasicUtils/BasicUtilsFwd.hpp>  // for NekManager
#include <LibUtilities/BasicUtils/SharedArray.hpp>


namespace Nektar
{
    namespace LibUtilities
    {

        // /////////////////////////////////////
        // General matrix and vector stuff        
        template< typename T > NekVector<T> GetColumn(const NekMatrix<T> & matA, int n);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> & SetColumn(NekMatrix<NekDouble> & matA, int n, const NekVector<NekDouble> & x);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> GetE(int rows, int n);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> Invert(const NekMatrix<NekDouble> & matA);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTranspose(const NekMatrix<NekDouble> & matA);
        LIB_UTILITIES_EXPORT int GetSize(const Array<OneD, const NekDouble> & x);
        LIB_UTILITIES_EXPORT int GetSize(const NekVector<NekDouble> & x);
        LIB_UTILITIES_EXPORT int GetDegree(int nBasisFunctions);
        LIB_UTILITIES_EXPORT NekDouble MakeRound(NekDouble);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> ToVector( const Array<OneD, const NekDouble> & x );
        LIB_UTILITIES_EXPORT Array<OneD, NekDouble> ToArray( const NekVector<NekDouble> & x );
        LIB_UTILITIES_EXPORT NekVector<NekDouble> Hadamard( const NekVector<NekDouble> & x, const NekVector<NekDouble> & y );
        LIB_UTILITIES_EXPORT NekVector<NekDouble> VectorPower( const NekVector<NekDouble> & x, NekDouble p );
        LIB_UTILITIES_EXPORT std::string MatrixToString( const NekMatrix<NekDouble> & A, int precision = 2, NekDouble threshold = 1e12 );
        LIB_UTILITIES_EXPORT std::string VectorToString( const NekVector<NekDouble> & v, int precision = 2, NekDouble threshold = 1e12 );
        
        // ////////////////////////////////////////////////////////////////
        // Polynomials(Jacobi, Legendre, and Dubiner) and its derivations        
        LIB_UTILITIES_EXPORT NekVector<NekDouble> LegendrePoly(int degree, const NekVector<NekDouble>& x);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> DubinerPoly(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> JacobiPoly(int degree, const NekVector<NekDouble>& x, NekDouble alpha, NekDouble beta);
        LIB_UTILITIES_EXPORT NekDouble JacobiPoly(int degree, NekDouble x, NekDouble alpha, NekDouble beta);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> LegendrePolyDerivative(int degree, const NekVector<NekDouble>& x);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> DubinerPolyXDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> DubinerPolyYDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> JacobiPolyDerivative(int degree, const NekVector<NekDouble>& x, int alpha, int beta);        

        // /////////////////////////////
        // Triangle stuff        
        LIB_UTILITIES_EXPORT NekVector<NekDouble> MakeDubinerQuadratureSystem(int nBasisFunctions);
        LIB_UTILITIES_EXPORT SharedNekMatrixPtr MakeVmatrixOfDubinerPolynomial(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> MakeQuadratureWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetInterpolationMatrix(const NekVector<NekDouble>& x,  const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
        LIB_UTILITIES_EXPORT Points<NekDouble>::MatrixSharedPtrType GetYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
        LIB_UTILITIES_EXPORT Points<NekDouble>::MatrixSharedPtrType GetXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
        
        // /////////////////////////////
        // Tetrahedron stuff
        LIB_UTILITIES_EXPORT int GetTetDegree(int nBasisFunc);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> TetrahedralBasis(int p, int q, int r, const NekVector<NekDouble>& x,
                                              const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT SharedNekMatrixPtr MakeVmatrixOfTet(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetInterpolationMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                       const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> MakeTetWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> MakeTetQuadratureSystem(int nBasisFunctions);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> TetXDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForTetXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForTetXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                             int degree);
        LIB_UTILITIES_EXPORT Points<NekDouble>::MatrixSharedPtrType GetTetXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                 const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi);

        LIB_UTILITIES_EXPORT NekVector<NekDouble> TetYDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForTetYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForTetYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                             int degree);
        LIB_UTILITIES_EXPORT Points<NekDouble>::MatrixSharedPtrType GetTetYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                 const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi);

        LIB_UTILITIES_EXPORT NekVector<NekDouble> TetZDerivative(int p, int q, int r, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                            const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForTetZDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetVandermondeForTetZDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z,
                                                             int degree);
        LIB_UTILITIES_EXPORT Points<NekDouble>::MatrixSharedPtrType GetTetZDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                 const NekVector<NekDouble>& z, const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi, const NekVector<NekDouble>& zi);


        // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Monomial Vandermonde stuff for Triangle : Useful to test triangle(integration, interpolation, and derivation) 
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        LIB_UTILITIES_EXPORT NekVector<NekDouble> GetIntegralOfMonomialVandermonde(int degree);
 
                 
        // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Monomial Vandermonde stuff for Tetrahedron : Useful to test tetrahedron(integration, interpolation, and derivation)
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& z, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetZDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z, int degree);
        LIB_UTILITIES_EXPORT NekMatrix<NekDouble> GetTetZDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& z);

    } // end of LibUtilities namespace
} // end of Nektar namespace

#endif //NODALUTIL_H
