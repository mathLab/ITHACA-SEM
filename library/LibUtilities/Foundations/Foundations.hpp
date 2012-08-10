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
            eFourier,			//!< Fourier Expansion \f$ \exp(i p\pi  z_i)\f$
            eGLL_Lagrange,		//!< Lagrange for SEM basis \f$ h_p(z_i) \f$
            eLegendre,			//!< Legendre Polynomials \f$ L_p(z_i) = P^{0,0}_p(z_i)\f$. Same as Ortho_A
            eChebyshev,			//!< Chebyshev Polynomials \f$ T_p(z_i) = P^{-1/2,-1/2}_p(z_i)\f$
            eMonomial,			//!< Monomial polynomials \f$ L_p(z_i) = z_i^{p}\f$
            eFourierSingleMode, //!< Fourier ModifiedExpansion with just the first mode   \f$ \exp(i \pi  z_i)\f$
			eFourierHalfModeRe, //!< Fourier Modified expansions with just the real part of the first mode  \f$ Re[\exp(i \pi  z_i)]\f$    
			eFourierHalfModeIm, //!< Fourier Modified expansions with just the imaginary part of the first mode  \f$ Im[\exp(i \pi  z_i)]\f$    
            eGauss_Lagrange,	//!< Lagrange Polynomials using the Gauss points \f$ h_p(z_i) \f$
            eDG_DG_Left,		//!< Derivative of the left correction function for DG FR  \f$ dGL_{p}(z_i) \f$
            eDG_DG_Right,		//!< Derivative of the Right correction function for DG FR \f$ dGR_{p}(z_i) \f$
            eDG_SD_Left,		//!< Derivative of the left correction function for SD FR  \f$ dGL_{p}(z_i) \f$
            eDG_SD_Right,		//!< Derivative of the Right correction function for SD FR \f$ dGR_{p}(z_i) \f$
            eDG_HU_Left,		//!< Derivative of the left correction function for HU FR  \f$ dGL_{p}(z_i) \f$
            eDG_HU_Right,		//!< Derivative of the Right correction function for HU FR \f$ dGR_{p}(z_i) \f$
            SIZE_BasisType		//!< Length of enum list
        };

        const char* const BasisTypeMap[] = 
        {
            "NoBasisType",
            "Ortho_A",
            "Ortho_B",
            "Ortho_C",
            "Modified_A",
            "Modified_B",
            "Modified_C",
            "Fourier",
            "GLL_Lagrange",
            "Legendre",
            "Chebyshev",
            "Monomial",
            "FourierSingleMode",
			"FourierHalfModeRe",
			"FourierHalfModeIm",
            "Gauss_Lagrange",
            "DG_DG_Left",
            "DG_DG_Right",
            "DG_SD_Left",
            "DG_SD_Right",
            "DG_HU_Left"
            "DG_HU_Right"
        };

        enum PointsType
        {
            eNoPointsType,
            eGaussGaussLegendre,            //!<  1D Gauss-Gauss-Legendre quadrature points
            eGaussRadauMLegendre,           //!<  1D Gauss-Radau-Legendre quadrature points, pinned at x=-1
            eGaussRadauPLegendre,           //!<  1D Gauss-Radau-Legendre quadrature points, pinned at x=1
            eGaussLobattoLegendre,          //!<  1D Gauss-Lobatto-Legendre quadrature points
            eGaussGaussChebyshev,           //!<  1D Gauss-Gauss-Chebyshev quadrature points
            eGaussRadauMChebyshev,          //!<  1D Gauss-Radau-Chebyshev quadrature points, pinned at x=-1
            eGaussRadauPChebyshev,          //!<  1D Gauss-Radau-Chebyshev quadrature points, pinned at x=1
            eGaussLobattoChebyshev,         //!<  1D Gauss-Lobatto-Legendre quadrature points
            eGaussRadauMAlpha0Beta1,        //!<  Gauss Radau pinned at x=-1, \f$ \alpha =    0, \beta =    1 \f$
            eGaussRadauMAlpha0Beta2,        //!<  Gauss Radau pinned at x=-1, \f$ \alpha =    0, \beta =    2 \f$
            eGaussRadauMAlpha1Beta0,        //!<  Gauss Radau pinned at x=-1, \f$ \alpha =    1, \beta =    0 \f$
            eGaussRadauMAlpha2Beta0,        //!<  Gauss Radau pinned at x=-1, \f$ \alpha =    2, \beta =    0 \f$
            eGaussKronrodLegendre,          //!<  1D Gauss-Kronrod-Legendre quadrature points
            eGaussRadauKronrodMLegendre,    //!<  1D Gauss-Radau-Kronrod-Legendre quadrature points, pinned at x=-1
            eGaussRadauKronrodMAlpha1Beta0, //!<  1D Gauss-Radau-Kronrod-Legendre pinned at x=-1, \f$ \alpha =    1, \beta =    0 \f$
            eGaussLobattoKronrodLegendre,   //!<  1D Lobatto Kronrod quadrature points
            ePolyEvenlySpaced,              //!<  1D Evenly-spaced points using Lagrange polynomial
            eFourierEvenlySpaced,           //!<  1D Evenly-spaced points using Fourier Fit
            eFourierSingleModeSpaced,       //!<  1D Non Evenly-spaced points for Single Mode analysis
            eBoundaryLayerPoints,           //!<  1D power law distribution for boundary layer points
            eNodalTriElec,                  //!<  2D Nodal Electrostatic Points on a Triangle
            eNodalTriFekete,                //!<  2D Nodal Fekete Points on a Triangle
            eNodalTriEvenlySpaced,          //!<  2D Evenly-spaced points on a Triangle
            eNodalTetEvenlySpaced,          //!<  3D Evenly-spaced points on a Tetrahedron
            eNodalTetElec,                  //!<  3D Nodal Electrostatic Points on a Tetrahedron
            eNodalPrismEvenlySpaced,        //!<  3D Evenly-spaced points on a Prism
            SIZE_PointsType                 //!<  Length of enum list
        };

        const std::string kPointsTypeStr[] = 
        {
            "NoPointsType",
            "GaussGaussLegendre",
            "GaussRadauMLegendre",
            "GaussRadauPLegendre",
            "GaussLobattoLegendre",
            "GaussGaussChebyshev",
            "GaussRadauMChebyshev",
            "GaussRadauPChebyshev",
            "GaussLobattoChebyshev",
            "GaussRadauMAlpha0Beta1",
            "GaussRadauMAlpha0Beta2",
            "GaussRadauMAlpha1Beta0",
            "GaussRadauMAlpha2Beta0",
            "GaussKronrodLegendre",
            "GaussRadauKronrodMLegendre",
            "GaussRadauKronrodMAlpha1Beta0",
            "GaussLobattoKronrodLegendre",
            "PolyEvenlySpaced",
            "FourierEvenlySpaced",
            "FourierSingleModeSpaced",
            "BoundaryLayerPoints",
            "NodalTriElec",
            "NodalTriFekete",
            "NodalTriEvenlySpaced",
            "NodalTetEvenlySpaced",
            "NodalTetElec",
            "NodalPrismEvenlySpaced"
        };

    } // end of namespace
} // end of namespace

#endif //FOUNDATIONS_H

