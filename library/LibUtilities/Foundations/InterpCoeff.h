///////////////////////////////////////////////////////////////////////////////
//
// File InterpCoeff.h
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
// Description: Definition of interpolation routines for Coefficients
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILTIES_FOUNDATIONS_INTERPCOEFF_H
#define NEKTAR_LIB_UTILTIES_FOUNDATIONS_INTERPCOEFF_H

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar { template <typename Dim, typename DataType> class Array; }

namespace Nektar
{    
    namespace LibUtilities
    {
        // Coefficient Space Interpolation methods

        // 1D Interpolation
        LIB_UTILITIES_EXPORT void InterpCoeff1D(
            const BasisKey                      &fbasis0,
            const Array<OneD, const NekDouble>  &from,
            const BasisKey                      &tbasis0,
                  Array<OneD,       NekDouble>  &to);

        // 2D Interpolation
        LIB_UTILITIES_EXPORT void InterpCoeff2D(
            const BasisKey                      &fbasis0,
            const BasisKey                      &fbasis1,
            const Array<OneD, const NekDouble>&  from,
            const BasisKey                      &tbasis0,
            const BasisKey                      &tbasis1,
                  Array<OneD,       NekDouble>  &to);
        
        LIB_UTILITIES_EXPORT void InterpCoeff2D(
            const BasisKey                      &fbasis0,
            const BasisKey                      &fbasis1,
            const NekDouble                     *from,
            const BasisKey                      &tbasis0,
            const BasisKey                      &tbasis1,
                  NekDouble                     *to);

        // 3D Interpolation
        LIB_UTILITIES_EXPORT void InterpCoeff3D(
            const BasisKey                      &fbasis0,
            const BasisKey                      &fbasis1,
            const BasisKey                      &fbasis2,
            const Array<OneD, const NekDouble>&  from,
            const BasisKey                      &tbasis0,
            const BasisKey                      &tbasis1,
            const BasisKey                      &tbasis2,
            Array<OneD, NekDouble>              &to);

        LIB_UTILITIES_EXPORT void InterpCoeff3D(
            const BasisKey                      &fbasis0,
            const BasisKey                      &fbasis1,
            const BasisKey                      &fbasis2,
            const NekDouble                     *from,
            const BasisKey                      &tbasis0,
            const BasisKey                      &tbasis1,
            const BasisKey                      &tbasis2,
                  NekDouble                     *to);
    }
}

#endif
