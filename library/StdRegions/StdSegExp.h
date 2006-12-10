///////////////////////////////////////////////////////////////////////////////
//
// File StdSegExp.h
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
// Description: Header file for Standard Segment Expansions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGIONS_STDSEGEXP_H
#define NEKTAR_LIBS_STDREGIONS_STDSEGEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion1D.h>
#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdExpMap.h>
#include <StdRegions/StdMatrix.h>

namespace Nektar
{
  namespace StdRegions
  {

    class StdSegExp: public StdExpansion1D
    {
    public:
      ///  Default constructor
      StdSegExp();

      /// \brief Constructor using BasisKey class for quadrature
      /// points and order definition
      StdSegExp(const BasisKey &Ba);

      /// \brief Constructor using BasisKey class for quadrature
      /// points and order definition where _coeffs and _phys are all
      /// set.
      StdSegExp(const BasisKey &Ba, double *coeffs, double *phys);

      ///Copy Constructor
      StdSegExp(const StdSegExp &T);

      ///Destructor
      ~StdSegExp();

      /// Return Shape of region, using  ShapeType enum list. i.e. Segment
      ShapeType DetShapeType()
      {
    return eSegment;
      };

      //----------------------------
      // Integration Methods
      //----------------------------

      /// \brief Integrate the physical point list \a inarray over region
      double Integral(const double *inarray);

      /// \brief Inner product of \a inarray over region with respect
      /// to the expansion basis (this)->_Base[0] and return in \a
      /// outarray
      void IProductWRTBase(const double * inarray, double * outarray);

      void FillMode(const int mode, double *outarray);

      //----------------------------------
      // Local Matrix Routines
      //----------------------------------

      /// \brief Generate local mass matrix \f$ {\bf M}[i][j] =
      /// \int^1_{-1} \phi_i(\xi_1) \phi_j(\xi_1) d\xi_1 \f$ in
      /// standard  region and store in \a outarray
      void GenMassMatrix(double * outarray);

      void GenLapMatrix(double * outarray);

      //// \brief Get the mass matrix attached to this expansion by
      //// using the StdMatrix manager _elmtmats and return the
      //// standard Matrix container
      StdMatContainer *GetMassMatrix();

      StdMatContainer *GetLapMatrix();

      //-----------------------------
      // Differentiation Methods
      //-----------------------------

      /// \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
      /// physical quadrature points in the expansion
      /// (i.e. (this)->_phys) and return in \a outarray.
      void Deriv(double * outarray);

      /// \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
      /// physical quadrature points given by \a inarray and return in \a
      /// outarray.
      void Deriv(const double *inarray, double * outarray);


      //----------------------------
      // Evaluations Methods
      //---------------------------

      ///  \brief Backward transform from coefficient space (stored in
      /// (this)->_coeffs) and evaluate at the physical quadrature points \a
      /// outarray
      void BwdTrans(double * outarray);

      /// \brief Forward transform from physical quadrature space
      /// stored in \a inarray and evaluate the expansion coefficients
      /// and store in \a (this)->_coeffs
      void FwdTrans(const double * inarray);

      /// \brief Single Point Evaluation: \f$ u(x) = \sum_p \phi_p(x)
      /// \hat{u}_p = \sum_p h_p(x) u(x_p)\f$
      double Evaluate(const double *Lcoords);

      void MapTo(EdgeOrientation dir, StdExpMap& Map);

      void SetInvInfo(StdMatContainer *mat, MatrixType Mform);

    protected:
      static StdMatrix s_elmtmats;

      /// \brief Inner product of \a inarray over region with respect
      /// to the expansion basis \a base and return in \a outarray
      void IProductWRTBase(const double *base, const double *inarray,
               double *outarray, int coll_check);

    private:

      virtual int v_GetNverts()
      {
	  return 2;
      }

      virtual ShapeType v_DetShapeType()
      {
	  return DetShapeType();
      };

      /// \brief Virtual call to integrate the physical point list \a
      /// inarray over region (see StdSegExp::Integral)
      virtual double v_Integral(const double *inarray )
      {
	  return Integral(inarray);
      }

      /// \brief Virtual call to StdSegExp::IProduct_WRT_B */
      virtual void v_IProductWRTBase(const double * inarray, double * outarray)
      {
	  IProductWRTBase(inarray,outarray);
      }

      virtual void v_FillMode(const int mode, double *outarray)
      {
	  FillMode(mode,outarray);
      }

      /// Virtual call to GenMassMatrix
      virtual void v_GenMassMatrix(double * outarray)
      {
	  GenMassMatrix(outarray);
      }
      
      virtual void v_GenLapMatrix(double * outarray)
      {
	  GenLapMatrix(outarray);
      }

      /// virtual call to GetMassMatrix
      virtual StdMatContainer *v_GetMassMatrix()
      {
	  return GetMassMatrix();
      }
      
      virtual StdMatContainer *v_GetLapMatrix()
      {
	  return GetLapMatrix();
      }
      
      /// Virtual call to StdSegExp::Deriv
      virtual void v_Deriv(double * outarray)
      {
	  Deriv(this->m_phys, outarray);
      }
      
      /// Virtual call to StdSegExp::Deriv
      virtual void v_StdDeriv(double * outarray)
      {
	  Deriv(this->m_phys, outarray);
      }
      
      /// Virtual call to StdSegExp::Deriv
      virtual void v_Deriv(const double *inarray, double * outarray)
      {
	  Deriv(inarray, outarray);
      }
      
      /// Virtual call to StdSegExp::Deriv
      virtual void v_StdDeriv(const double *inarray, double * outarray)
      {
	  Deriv(inarray, outarray);
      }
      
      /// Virtual call to StdSegExp::BwdTrans
      virtual void v_BwdTrans(double * outarray)
      {
	  BwdTrans(outarray);
      }
      
      /// Virtual call to StdSegExp::FwdTrans
      virtual void v_FwdTrans(const double * inarray)
      {
	  FwdTrans(inarray);
      }
      
      /// Virtual call to StdSegExp::Evaluate
      virtual double v_Evaluate(const double * Lcoords)
      {
	  return Evaluate(Lcoords);
      }
      
      virtual void v_MapTo(EdgeOrientation dir, StdExpMap &Map)
      {
	  MapTo(dir,Map);
      }
      

      virtual void v_SetInvInfo(StdMatContainer *mat, MatrixType Mform)
      {
	  SetInvInfo(mat,Mform);
      }
    };

  } //end of namespace
} //end of namespace

#endif //STDSEGEXP_H

/**
 * $Log: StdSegExp.h,v $
 * Revision 1.3  2006/07/02 17:16:19  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.2  2006/06/01 14:13:37  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:33  kirby
 * *** empty log message ***
 *
 * Revision 1.41  2006/03/13 18:29:35  sherwin
 *
 * Corrected error with definition of GetCoords
 *
 * Revision 1.40  2006/03/05 22:11:03  sherwin
 *
 * Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
 *
 * Revision 1.39  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.38  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.37  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/
