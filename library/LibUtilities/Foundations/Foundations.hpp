///////////////////////////////////////////////////////////////////////////////
//
// File Foundations.hpp
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
// Description: Definition of enum lists and constants
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FOUNDATIONS_H
#define FOUNDATIONS_H

#include <string>

namespace Nektar
{
    namespace LibUtiltities
    {
        enum BasisType
        {
            eOrtho_A,     //!< Principle Orthogonal Functions \f$\widetilde{\psi}^a_p(z_i)\f$
            eOrtho_B,     //!< Principle Orthogonal Functions \f$\widetilde{\psi}^b_{pq}(z_i)\f$
            eOrtho_C,     //!< Principle Orthogonal Functions \f$\widetilde{\psi}^c_{pqr}(z_i)\f$
            eModified_A,  //!< Principle Modified Functions \f$ \phi^a_p(z_i) \f$
            eModified_B,  //!< Principle Modified Functions \f$ \phi^b_{pq}(z_i) \f$
            eModified_C,  //!< Principle Modified Functions \f$ \phi^c_{pqr}(z_i) \f$
            eFourier,     //!< Fourier Expansion \f$ \exp(i p\pi  z_i)\f$
            eGLL_Lagrange,//!< Lagrange for SEM basis \f$ h_p(z_i) \f$
            eLegendre,    //!< Legendre Polynomials \f$ L_p(z_i) = P^{0,0}_p(z_i)\f$. Same as Ortho_A
            eChebyshev,   //!< Chebyshev Polynomials \f$ T_p(z_i) = P^{-1/2,-1/2}_p(z_i)\f$
            SIZE_BasisType //!< Length of enum list
        };
    

        enum PointsType
        {
            eGauss,               //!< 1D Gauss quadrature points
            eLobatto,             //!< 1D Gauss Lobatto quadrature points
            eRadauM,              //!< 1D Gauss Radau quadrature points including (z=-1)
            eRadauP,              //!< 1D Gauss Radau quadrature points including (z=+1)
            ePolyEvenSp,          //!< 1D Evenly-spaced points using Lagrange polynomial
            eFourierEvenSp,       //!< 1D Evenly-spaced points using Fourier Fit
            eArbitrary,           //!< 1D Arbitrary point distribution
            eNodalTriElec,        //!< 2D Nodal Electrostatic Points on a Triangle
            eNodalTriFekete,      //!< 2D Nodal Fekete Points on a Triangle
            eNodalTetElec,        //!< 3D Nodal Electrostatic Points on a Tetrahedron
            SIZE_PointsType       //!< Length of enum list
        };


        enum PointsIdentifier
        {
            eWildcard,                  //!< Indentifier to use when PointType is sufficient to uniquely specify things
            eGaussChebyshevFirstKind,   //!< Gauss \f$ \alpha = -1/2, \beta = -1/2 \f$
            eGaussChebyshevSecondKind,  //!< Gauss \f$ \alpha =  1/2, \beta =  1/2 \f$
            eGaussLegendre,             //!< Gauss \f$ \alpha =    0, \beta =    0 \f$
            eGaussAlpha0Beta1,          //!< Gauss \f$ \alpha =    0, \beta =    1 \f$
            eGaussAlpha0Beta2,          //!< Gauss \f$ \alpha =    0, \beta =    2 \f$
            SIZE_PointsIdentifier       //!< Length of enum list
        };

    } // end of namespace
} // end of namespace

#endif //FOUNDATIONS_H

