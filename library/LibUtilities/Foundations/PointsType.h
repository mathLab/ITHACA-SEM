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

#ifndef NEKTAR_LIB_UTILITIES_POINTS_TYPE_H
#define NEKTAR_LIB_UTILITIES_POINTS_TYPE_H

#include <vector>

namespace Nektar
{
    namespace LibUtilities
    {

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
            eBoundaryLayerPointsRev,        //!<  1D power law distribution for boundary layer points
            eNodalTriElec,                  //!<  2D Nodal Electrostatic Points on a Triangle
            eNodalTriFekete,                //!<  2D Nodal Fekete Points on a Triangle
            eNodalTriEvenlySpaced,          //!<  2D Evenly-spaced points on a Triangle
            eNodalTetEvenlySpaced,          //!<  3D Evenly-spaced points on a Tetrahedron
            eNodalTetElec,                  //!<  3D Nodal Electrostatic Points on a Tetrahedron
            eNodalPrismEvenlySpaced,        //!<  3D Evenly-spaced points on a Prism
            eNodalPrismElec,                //!<  3D electrostatically spaced points on a Prism
            eNodalTriSPI,                   //!<  2D Nodal Symmetric positive internal triangle (Whitherden, Vincent)
            eNodalTetSPI,                   //!<  3D Nodal Symmetric positive internal tet (Whitherden, Vincent)
            eNodalPrismSPI,                 //!<  3D prism SPI
            eNodalQuadElec,                 //!<  2D GLL for quad
            eNodalHexElec,                  //!<  3D GLL for hex
            SIZE_PointsType                 //!<  Length of enum list
        };

        static std::vector<LibUtilities::PointsType> NullPointsTypeVector;
    }
}

#endif
