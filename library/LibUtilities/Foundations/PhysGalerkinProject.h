///////////////////////////////////////////////////////////////////////////////
//
// File PhysGalerkinProject.h
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
// Description: Definition of Physcal space (1D tensor based) Galerkin Projection
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILTIES_FOUNDATIONS_PHYSGALERKIN_H
#define NEKTAR_LIB_UTILTIES_FOUNDATIONS_PHYSGALERKIN_H

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar { template <typename Dim, typename DataType> class Array; }


namespace Nektar
{    
    namespace LibUtilities
    {     

        // Physical Space Interpolation methods
        
        LIB_UTILITIES_EXPORT void PhysGalerkinProject1D(const BasisKey &fbasis0, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0, 
                      Array<OneD, NekDouble> &to);

        LIB_UTILITIES_EXPORT void PhysGalerkinProject1D(const PointsKey &fpoints0, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0, 
                      Array<OneD, NekDouble> &to);

        LIB_UTILITIES_EXPORT void PhysGalerkinProject1D(const BasisKey &fbasis0, 
                      const NekDouble *from,  
                      const BasisKey &tbasis0, 
                      NekDouble *to);

        LIB_UTILITIES_EXPORT void PhysGalerkinProject1D(const PointsKey &fpoints0, 
                      const NekDouble *from,  
                      const PointsKey &tpoints0, 
                      NekDouble *to);

        // 2D PhysGalkerinProjectolation
        LIB_UTILITIES_EXPORT void PhysGalerkinProject2D(const BasisKey &fbasis0, 
                      const BasisKey &fbasis1, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0,
                      const BasisKey &tbasis1,
                      Array<OneD, NekDouble> &to);

        LIB_UTILITIES_EXPORT void PhysGalerkinProject2D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      Array<OneD, NekDouble> &to);

        LIB_UTILITIES_EXPORT void PhysGalerkinProject2D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const NekDouble *from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      NekDouble *to);


    } // end of namespace
} // end of namespace

#endif //FOUNDATIONS_H

