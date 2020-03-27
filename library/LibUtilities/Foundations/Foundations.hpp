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

#ifndef FOUNDATIONS_H
#define FOUNDATIONS_H

#include <LibUtilities/Foundations/BasisType.h>
#include <LibUtilities/Foundations/PointsType.h>
#include <string>

namespace Nektar
{
    namespace LibUtilities
    {
        const char* const BasisTypeMap[] =
        {
            "NoBasisType",
            "Ortho_A",
            "Ortho_B",
            "Ortho_C",
            "Modified_A",
            "Modified_B",
            "Modified_C",
            "OrthoPyr_C",
            "ModifiedPyr_C",
            "Fourier",
            "GLL_Lagrange",
            "Gauss_Lagrange",
            "Legendre",
            "Chebyshev",
            "Monomial",
            "FourierSingleMode",
            "FourierHalfModeRe",
            "FourierHalfModeIm"
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
            "BoundaryLayerPointsRev",
            "NodalTriElec",
            "NodalTriFekete",
            "NodalTriEvenlySpaced",
            "NodalTetEvenlySpaced",
            "NodalTetElec",
            "NodalPrismEvenlySpaced",
            "NodalPrismElec",
            "NodalTriSPI",
            "NodalTetSPI",
            "NodalPrismSPI",
            "NodalQuadElec",
            "NodalHexElec"
        };
    } // end of namespace
} // end of namespace

#endif //FOUNDATIONS_H
