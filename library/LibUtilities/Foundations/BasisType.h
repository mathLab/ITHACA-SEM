///////////////////////////////////////////////////////////////////////////////
//
// File BasisType.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIS_TYPE_H
#define NEKTAR_LIB_UTILITIES_BASIS_TYPE_H


namespace Nektar
{    
    namespace LibUtilities
    {     
        enum BasisType
        {
            eNoBasisType,
            eOrtho_A,			//!< Principle Orthogonal Functions \f$\widetilde{\psi}^a_p(z_i)\f$
            eOrtho_B,			//!< Principle Orthogonal Functions \f$\widetilde{\psi}^b_{pq}(z_i)\f$
            eOrtho_C,			//!< Principle Orthogonal Functions \f$\widetilde{\psi}^c_{pqr}(z_i)\f$
            eModified_A,		//!< Principle Modified Functions \f$ \phi^a_p(z_i) \f$
            eModified_B,		//!< Principle Modified Functions \f$ \phi^b_{pq}(z_i) \f$
            eModified_C,		//!< Principle Modified Functions \f$ \phi^c_{pqr}(z_i) \f$
            eOrthoPyr_C,		//!< Principle Orthogonal Functions \f$\widetilde{\psi}^c_{pqr}(z_i) for Pyramids\f$
            eModifiedPyr_C,		//!< Principle Modified Functions \f$ \phi^c_{pqr}(z_i) for Pyramids\f$
            eFourier,			//!< Fourier Expansion \f$ \exp(i p\pi  z_i)\f$
            eGLL_Lagrange,		//!< Lagrange for SEM basis \f$ h_p(z_i) \f$
            eGauss_Lagrange,	//!< Lagrange Polynomials using the Gauss points \f$ h_p(z_i) \f$
            eLegendre,			//!< Legendre Polynomials \f$ L_p(z_i) = P^{0,0}_p(z_i)\f$. Same as Ortho_A
            eChebyshev,			//!< Chebyshev Polynomials \f$ T_p(z_i) = P^{-1/2,-1/2}_p(z_i)\f$
            eMonomial,			//!< Monomial polynomials \f$ L_p(z_i) = z_i^{p}\f$
            eFourierSingleMode, //!< Fourier ModifiedExpansion with just the first mode   \f$ \exp(i \pi  z_i)\f$
            eFourierHalfModeRe, //!< Fourier Modified expansions with just the real part of the first mode  \f$ Re[\exp(i \pi  z_i)]\f$    
            eFourierHalfModeIm, //!< Fourier Modified expansions with just the imaginary part of the first mode  \f$ Im[\exp(i \pi  z_i)]\f$    
            SIZE_BasisType		//!< Length of enum list
        };
    }
}

#endif
