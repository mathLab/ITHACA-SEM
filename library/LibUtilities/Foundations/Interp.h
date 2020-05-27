///////////////////////////////////////////////////////////////////////////////
//
// File Interp.h
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
// Description: Definition of interpolation routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILTIES_FOUNDATIONS_INTERP_H
#define NEKTAR_LIB_UTILTIES_FOUNDATIONS_INTERP_H

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

//#include <LibUtilities/BasicUtils/BasicUtilsFwd.hpp>  // for SharedArray
namespace Nektar { template <typename Dim, typename DataType> class Array; }


namespace Nektar
{    
    namespace LibUtilities
    {     

        // Physical Space Interpolation methods
        


            /** \brief this function interpolates a 1D function \f$f\f$ evaluated
             *  at the quadrature points of the basis \a fbasis0 to the 
             *  function values at the quadrature points of the basis \a tbasis0
             *
             *  Given a function \f$ f\f$ evaluated at the \a Q quadrature points
             *  of the basis \a fbasis0, this routine calculates, using 
             *  \a (Q-1)th order polynomial interpolation, the function values
             *  at the \a Q2 quadrature points of the basis \a tbasis0.
             *
             *  \param fbasis0 the basis at which's quadrature points the 
             *  function is given
             *  \param from the array containg the function \f$ f\f$  evaluated
             *   at the quadrature points of \a fbasis0
             *  \param tbasis0 the basis to which's quadrature points the 
             *  function should be interpolated
             *  \param to the array containg the function \f$ f\f$  evaluated
             *   at the quadrature points of \a tbasis0 (output of the function)
             */
        // 1D Interpolation
        LIB_UTILITIES_EXPORT void Interp1D(const BasisKey &fbasis0, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0, 
                      Array<OneD, NekDouble> &to);

        LIB_UTILITIES_EXPORT void Interp1D(const PointsKey &fpoints0, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0, 
                      Array<OneD, NekDouble> &to);

        LIB_UTILITIES_EXPORT void Interp1D(const BasisKey &fbasis0, 
                      const NekDouble *from,  
                      const BasisKey &tbasis0, 
                      NekDouble *to);

        LIB_UTILITIES_EXPORT void Interp1D(const PointsKey &fpoints0, 
                      const NekDouble *from,  
                      const PointsKey &tpoints0, 
                      NekDouble *to);


            /** \brief this function interpolates a 2D function \f$f\f$ evaluated
             *  at the quadrature points of the 2D basis, constructed by 
             *  \a fbasis0 and \a fbasis1, to the function values at the 
             *  quadrature points of the 2D basis, constructed by \a tbasis0 and 
             *  \a tbasis1
             *
             *  Given a function \f$ f\f$ evaluated at the \a Q quadrature points
             *  of the first expansion basis, this routine calculates, using 
             *  \a (Q-1)th order polynomial interpolation, the function values
             *  at the \a Q2 quadrature points of the second basis.
             *
             *  \param fbasis0 the basis at which's quadrature points the 
             *  function is given
             *  \param from the array containg the function \f$ f\f$  evaluated
             *   at the quadrature points of \a fbasis0
             *  \param tbasis0 the basis to which's quadrature points the 
             *  function should be interpolated
             *  \param to the array containg the function \f$ f\f$  evaluated
             *   at the quadrature points of \a tbasis0 (output of the function)
             */
        // 2D Interpolation
        LIB_UTILITIES_EXPORT void Interp2D(const BasisKey &fbasis0, 
                      const BasisKey &fbasis1, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0,
                      const BasisKey &tbasis1,
                      Array<OneD, NekDouble> &to);

        LIB_UTILITIES_EXPORT void Interp2D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      Array<OneD, NekDouble> &to);

        LIB_UTILITIES_EXPORT void Interp2D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const NekDouble *from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      NekDouble *to);


        /** \brief this function interpolates a 3D function \f$f\f$ evaluated
         *  at the quadrature points of the 3D basis, constructed by 
         *  \a fbasis0, \a fbasis1, and \a fbasis2 to the function values at the 
         *  quadrature points of the 3D basis, constructed by \a tbasis0, 
         *  \a tbasis1, and \a tbasis2.
         *
         *  Given a function \f$ f\f$ evaluated at the \a Q quadrature points
         *  of the first expansion basis, this routine calculates, using 
         *  \a (Q-1)th order polynomial interpolation, the function values
         *  at the \a Q2 quadrature points of the second basis, and the function
         *  values at the \a Q3 quadrature points of the third basis.
         *
         *  \param fbasis0 the basis at which's quadrature points the 
         *  function is given
         *  \param from the array containg the function \f$ f\f$  evaluated
         *   at the quadrature points of \a fbasis0
         *  \param tbasis0 the basis to which's quadrature points the 
         *  function should be interpolated
         *  \param to the array containg the function \f$ f\f$  evaluated
         *   at the quadrature points of \a tbasis0 (output of the function)
         */
        // 3D interpolation
        LIB_UTILITIES_EXPORT void Interp3D(const BasisKey &fbasis0, 
                      const BasisKey &fbasis1, 
                      const BasisKey &fbasis2, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0,
                      const BasisKey &tbasis1,
                      const BasisKey &tbasis2,
                      Array<OneD, NekDouble> &to);


        LIB_UTILITIES_EXPORT void Interp3D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const PointsKey &fpoints2, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      const PointsKey &tpoints2,
                      Array<OneD, NekDouble> &to);

	    LIB_UTILITIES_EXPORT void Interp3D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1,
                      const PointsKey &fpoints2,  
                      const NekDouble *from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      const PointsKey &tpoints2,
                      NekDouble *to);



    } // end of namespace
} // end of namespace

#endif //FOUNDATIONS_H

