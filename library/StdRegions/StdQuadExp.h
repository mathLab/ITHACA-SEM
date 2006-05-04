///////////////////////////////////////////////////////////////////////////////
//
// File StdQuadExp.h
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
// Description: Header field for Quadrilateral routines built upon
// StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDQUADEXP_H
#define STDQUADEXP_H

#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdMatrix.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdExpMap.h>

namespace Nektar
{
  namespace StdRegions
  {

    class StdQuadExp: public StdExpansion2D
    {

    public:

      StdQuadExp();

      /** \brief Constructor using BasisKey class for quadrature
    points and order definition */
      StdQuadExp(const BasisKey &Ba, const BasisKey &Bb);

      /** \brief Constructor using BasisKey class for quadrature
      points and order nition where _coeffs and _phys are all
      set. */
      StdQuadExp(const BasisKey &Ba, const BasisKey &Bb, double *coeffs,
         double *phys);

      /// Copy Constructor
      StdQuadExp(const StdQuadExp &T);

      /// Destructor
      ~StdQuadExp();

      // Return Shape of region, using ShapeType enum
      /// list. i.e. Quadrilateral
      ShapeType DetShapeType()
      {
    return eQuadrilateral;
      };

      void FillMode(int mode, double *array);

      //////////////////////////////
      // Integration Methods
      //////////////////////////////

      double Integral(const double *inarray);

      void IProductWRTBase(const double * inarray, double * outarray);

      void IProductWRTBase(const double *base0, const double *base1,
               const double *inarray, double *outarray,
               int coll_check);

      //----------------------------------
      // Local Matrix Routines
      //----------------------------------

      void GenMassMatrix(double * outarray);
      void GenLapMatrix(double * outarray);

      StdMatContainer * GetMassMatrix();
      StdMatContainer * GetLapMatrix();

      //----------------------------
      // Differentiation Methods
      //----------------------------

      void Deriv(double * outarray_d1, double *outarray_d2);
      void Deriv(const double *inarray, double * outarray_d1,
         double *outarray_d2);

      //----------------------------
      // Evaluations Methods
     //-----------------------------

      void BwdTrans(double * outarray);
      void FwdTrans(const double * inarray);

      double Evaluate(const double * coords);
      void MapTo(StdSegExp& edge, const int eid,
         const EdgeOrientation eorient, StdExpMap &Map);

    protected:

      static StdMatrix s_elmtmats;

    private:

      virtual ShapeType v_DetShapeType()
      {
	return DetShapeType();
      }

      virtual void v_FillMode(int mode, double *array)
      {
	FillMode(mode,array);
      }

      virtual double v_Integral(const double *inarray )
      {
	return Integral(inarray);
      }

      virtual void v_IProductWRTBase(const double * inarray, double * outarray)
      {
	IProductWRTBase(inarray,outarray);
      }

      virtual void v_GenMassMatrix(double * outarray)
      {
	GenMassMatrix(outarray);
      }

      virtual void v_GenLapMatrix(double * outarray)
      {
	GenLapMatrix(outarray);
      }

      virtual StdMatContainer *v_GetMassMatrix()
      {
	return GetMassMatrix();
      }
      
      virtual StdMatContainer * v_GetLapMatrix()
      {
	return GetLapMatrix();
      }

      virtual void v_Deriv(double * outarray_d0, double *outarray_d1)
      {
	Deriv(this->m_phys, outarray_d0, outarray_d1);
      }

      virtual void v_StdDeriv(double * outarray_d0, double *outarray_d1)
      {
	Deriv(this->m_phys, outarray_d0, outarray_d1);
      }

      virtual void v_Deriv(const double *inarray, double * outarray_d0,
			   double *outarray_d1)
      {
	Deriv(inarray, outarray_d0, outarray_d1);
      }

      virtual void v_StdDeriv(const double *inarray, double * outarray_d0,
			     double *outarray_d1)
      {
	Deriv(inarray, outarray_d0, outarray_d1);
      }

      virtual void v_BwdTrans(double * outarray)
      {
	BwdTrans(outarray);
      }
      
      virtual void v_FwdTrans(const double * inarray)
      {
	FwdTrans(inarray);
      }

      virtual double v_Evaluate(const double * coords)
      {
	return Evaluate(coords);
      }
      
    };
    
  } //end of namespace
} //end of namespace

#endif //STDQUADEXP_H

/**
 * $Log: StdQuadExp.h,v $
 * Revision 1.34  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.33  2006/03/05 23:17:53  sherwin
 *
 * Corrected to allow MMatrix1D and MMatrix2D to execute properly
 *
 * Revision 1.32  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.31  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.30  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/




