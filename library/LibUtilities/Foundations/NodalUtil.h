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
// Description: 2D Nodal Triangle Fekete Utilities header file --
//              Basis function, Interpolation, Integral, Derivation, etc. 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALUTIL_H
#define NODALUTIL_H

#include <iostream>

#include <math.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar {
    namespace LibUtilities{

        template< typename T > NekVector<T> GetColumn(const NekMatrix<T> & matA, int n);
        NekMatrix<NekDouble> & SetColumn(NekMatrix<NekDouble> & matA, int n, const NekVector<NekDouble> & x);
        NekVector<NekDouble> GetE(int rows, int n);
        NekMatrix<NekDouble> Invert(const NekMatrix<NekDouble> & matA);
        NekMatrix<NekDouble> GetTranspose(const NekMatrix<NekDouble> & matA);
        int GetSize(const ConstArray<OneD, NekDouble> & x);
        int GetSize(const NekVector<NekDouble> & x);
        int GetDegree(int nBasisFunctions);
        NekVector<NekDouble> ToVector( const ConstArray<OneD, NekDouble> & x );
        Array<OneD, NekDouble> ToArray( const NekVector<NekDouble> & x );
        NekVector<NekDouble> Hadamard( const NekVector<NekDouble> & x, const NekVector<NekDouble> & y );

        NekVector<NekDouble> MakeDubinerQuadratureSystem(int nBasisFunctions);
        SharedNekMatrixPtr MakeVmatrixOfDubinerPolynomial(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> MakeQuadratureWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);

        NekVector<NekDouble> LegendrePoly(int degree, const NekVector<NekDouble>& x);
        NekVector<NekDouble> JacobiPoly(int degree, const NekVector<NekDouble>& x, int alpha, int beta);
        NekVector<NekDouble> DubinerPoly(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
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
                                                                                          
        NekVector<NekDouble> LegendrePolyDerivative(int degree, const NekVector<NekDouble>& x);
        NekVector<NekDouble> DubinerPolyXDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> DubinerPolyYDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> JacobiPolyDerivative(int degree, const NekVector<NekDouble>& x, int alpha, int beta);

        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> GetYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> GetIntegralOfMonomialVandermonde(int degree);
        NekDouble MakeRound(NekDouble);

    } // end of LibUtilities namespace
} // end of Nektar namespace

#endif //NODALUTIL_H
