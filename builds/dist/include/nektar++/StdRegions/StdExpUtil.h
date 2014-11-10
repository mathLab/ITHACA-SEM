///////////////////////////////////////////////////////////////////////////////
//
// File StdExpUtil.h
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

#ifndef STDEXPUTIL_H
#define STDEXPUTIL_H

#include <iosfwd>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.h>

namespace Nektar
{
    namespace StdRegions
    {
    
        // /////////////////////////////////////
        // General matrix and vector stuff        
        template< typename T > NekVector<T> GetColumn(const NekMatrix<T> & matA, int n);
        NekMatrix<NekDouble> & SetColumn(NekMatrix<NekDouble> & matA, int n, const NekVector<NekDouble> & x);
        NekVector<NekDouble> GetE(int rows, int n);
        NekMatrix<NekDouble> Invert(const NekMatrix<NekDouble> & matA);
        NekMatrix<NekDouble> GetTranspose(const NekMatrix<NekDouble> & matA);
        int GetSize(const Array<OneD, const NekDouble> & x);
        int GetSize(const NekVector<NekDouble> & x);
        int GetDegree(int nBasisFunctions);
        NekDouble MakeRound(NekDouble);
        NekVector<NekDouble> ToVector( const Array<OneD, const NekDouble> & x );
        Array<OneD, NekDouble> ToArray( const NekVector<NekDouble> & x );
        NekVector<NekDouble> Hadamard( const NekVector<NekDouble> & x, const NekVector<NekDouble> & y );
        NekVector<NekDouble> VectorPower( const NekVector<NekDouble> & x, NekDouble p );
        std::string MatrixToString( const NekMatrix<NekDouble> & A, int precision = 2, NekDouble threshold = 1e12 );
        std::string VectorToString( const NekVector<NekDouble> & v, int precision = 2, NekDouble threshold = 1e12 );
        
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

    }
}

#endif
