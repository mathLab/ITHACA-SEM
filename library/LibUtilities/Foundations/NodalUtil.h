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

        template< typename T > NekVector<T> getColumn(const NekMatrix<T> & A, int n);
        NekMatrix<NekDouble> & setColumn(NekMatrix<NekDouble> & A, int n, const NekVector<NekDouble> & x);
        NekVector<NekDouble> getE(int M, int n);
        NekMatrix<NekDouble> invert(const NekMatrix<NekDouble> & A);
        NekMatrix<NekDouble> getTranspose(const NekMatrix<NekDouble> & A);
        int getSize(const ConstArray<OneD, NekDouble> & x);
        int getSize(const NekVector<NekDouble> & x);
        int getDegree(int nBasisFunctions);
        NekVector<NekDouble> toVector( const ConstArray<OneD, NekDouble> & x );
        Array<OneD, NekDouble> toArray( const NekVector<NekDouble> & x );
        NekVector<NekDouble> hadamard( const NekVector<NekDouble> & x, const NekVector<NekDouble> & y );

        NekVector<NekDouble> makeDubinerQuadratureSystem(int nBasisFunctions);
        SharedNekMatrixPtr makeVmatrixOfDubinerPolynomial(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> makeQuadratureWeights(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);

        NekVector<NekDouble> LegendrePoly(int degree, const NekVector<NekDouble>& x);
        NekVector<NekDouble> JacobiPoly(int degree, const NekVector<NekDouble>& x, int alpha, int beta);
        NekVector<NekDouble> DubinerPoly(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> getVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> getVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);      
        NekMatrix<NekDouble> getVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> getVandermondeForXDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);       
        NekMatrix<NekDouble> getVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> getVandermondeForYDerivative(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> getInterpolationMatrix(const NekVector<NekDouble>& x,  const NekVector<NekDouble>& y,
                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
        Points<NekDouble>::MatrixSharedPtrType getYDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
        Points<NekDouble>::MatrixSharedPtrType getXDerivativeMatrix(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y,
                                                                    const NekVector<NekDouble>& xi, const NekVector<NekDouble>& yi);
                                                                                          
        NekVector<NekDouble> LegendrePolyDerivative(int degree, const NekVector<NekDouble>& x);
        NekVector<NekDouble> DubinerPolyXDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> DubinerPolyYDerivative(int p, int q, const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> JacobiPolyDerivative(int degree, const NekVector<NekDouble>& x, int alpha, int beta);

        NekMatrix<NekDouble> getMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> getMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> getXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> getXDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekMatrix<NekDouble> getYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y, int degree);
        NekMatrix<NekDouble> getYDerivativeOfMonomialVandermonde(const NekVector<NekDouble>& x, const NekVector<NekDouble>& y);
        NekVector<NekDouble> getIntegralOfMonomialVandermonde(int degree);
        NekDouble round(NekDouble);

    } // end of LibUtilities namespace
} // end of Nektar namespace

#endif //NODALUTIL_H
