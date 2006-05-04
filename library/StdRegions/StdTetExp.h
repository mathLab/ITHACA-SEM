///////////////////////////////////////////////////////////////////////////////
//
// File StdTetExp.h
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
// Description: Header field for tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDTETEXP_H
#define STDTETEXP_H

#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdMatrix.h>

namespace Nektar
{
  namespace StdRegions
  {

    class StdTetExp: public StdExpansion3D
    {

    public:

      /** \brief Constructor using BasisKey class for quadrature
	  points and order definition */
      StdTetExp(const BasisKey &Ba, const BasisKey &Bb, const BasisKey &Bc);

      /** \brief Constructor using BasisKey class for quadrature
	  points and order definition where _coeffs and _phys are all
	  set. */
      StdTetExp(const BasisKey &Ba, const BasisKey &Bb, const BasisKey &Bc,
        double *coeffs, double *phys);

      /// Copy Constructor
      StdTetExp(const StdTetExp &T);

      /// Destructor
      ~StdTetExp();

      /// Return Shape of region, using  ShapeType enum list. i.e. Tetrahedron
      ShapeType DetShapeType()
      {
	return eTetrahedron;
      }


    protected:

      static StdMatrix s_elmtmats;

    private:

      virtual ShapeType v_DetShapeType()
      {
	return DetShapeType();
      }


  };

} //end of namespace
} //end of namespace

#endif //STDTETEXP_H

/**
 * $Log: StdTetExp.h,v $
 * Revision 1.24  2006/03/12 21:59:48  sherwin
 *
 * compiling version of LocalRegions
 *
 * Revision 1.23  2006/03/06 17:12:46  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.22  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.21  2006/03/01 08:25:05  sherwin
 *
 * First compiling version of StdRegions
 *
 **/

