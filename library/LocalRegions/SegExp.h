///////////////////////////////////////////////////////////////////////////////
//
// File SegExp.h
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
// Description: Header file for SegExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SEGEXP_H
#define SEGEXP_H

#include <StdRegions/StdSegExp.h>
#include <SpatialDomains/SegGeom.h>

#include <LocalRegions/MetricRelatedInfo.h>
#include <boost/shared_ptr.hpp>

#include <fstream>

namespace Nektar
{
  namespace LocalRegions 
  {

    class SegExp: public StdRegions::StdSegExp
    {

    public:
      
      /// \brief Constructor using BasisKey class for quadrature
      /// points and order definition 
      SegExp(const StdRegions::BasisKey &Ba, SpatialDomains::SegGeomSharedPtr geom);
    
      /// \brief Constructor using BasisKey class for quadrature
      /// points and order definition where _coeffs and _phys are all set. 
      SegExp(const StdRegions::BasisKey &Ba, double *coeffs, double *phys, 
	     SpatialDomains::SegGeomSharedPtr geom);
    
      ///Copy Constructor
      SegExp(const SegExp &S);
    
      ///Destructor
      ~SegExp();
    
      /// Return Shape of region, using  ShapeType enum list. i.e. Segment  
      StdRegions::ShapeType DetShapeType() 
      { 
	return StdRegions::eSegment;
      }    

      MetricRelatedInfoSharedPtr GenGeoFac();    

      inline void SetGeoFac(MetricRelatedInfoSharedPtr minfo)
      {
	m_minfo = minfo;
      }
    
      void GetCoords(double **coords);
      
      void GetCoord(const double *Lcoords, double *coords);
    

      SpatialDomains::SegGeomSharedPtr GetGeom()
      {
	return m_geom;
      }

      void WriteToFile(FILE *outfile);
      void WriteToFile(std::ofstream &outfile, const int dumpVar);
    
    
      StdRegions::StdMatContainer * GetMassMatrix() 
      {
	return m_minfo->GetMassMatrix(this); 
      }

      //----------------------------
      // Integration Methods
      //----------------------------

      /// \brief Integrate the physical point list \a inarray over region
      double Integral(const double *inarray);

    
      /// \brief  Inner product of \a inarray over region with respect to the 
      /// expansion basis (this)->_Base[0] and return in \a outarray 
      void IProductWRTBase(const double * inarray, double * outarray);


      //-----------------------------
      // Differentiation Methods
      //-----------------------------
    
      /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
	  physical quadrature points in the expansion (i.e. (this)->_phys)
	  and return in \a outarray. */
      void Deriv(double * outarray);

      
      /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
	  physical quadrature points given by \a inarray and return in \a
	  outarray. */
      void Deriv(const double *inarray, double * outarray)
      {
	Deriv(1,inarray,&outarray);
      }


      void Deriv(const int n, double **outarray);
      void Deriv(const int n, const double *inarray, double ** outarray);

    //----------------------------
    // Evaluations Methods
    //---------------------------
    
    /** \brief Forward transform from physical quadrature space
	stored in \a inarray and evaluate the expansion coefficients and
	store in \a (this)->_coeffs  */
      void FwdTrans(const double * inarray);
    

      double Evaluate(const double *coord);

    protected:
      int m_id;
      int m_field;
    
      SpatialDomains::SegGeomSharedPtr m_geom;
      MetricRelatedInfoSharedPtr       m_minfo;
      
      /// \brief  Inner product of \a inarray over region with respect to
      /// the expansion basis \a base and return in \a outarray 
      inline void IProductWRTBase(const double *base, const double *inarray, 
				  double *outarray, int coll_check);
    
    private:
    
      virtual StdRegions::ShapeType v_DetShapeType() 
      {
		return DetShapeType();
      }

      virtual MetricRelatedInfoSharedPtr v_GenGeoFac()
      {
		return GenGeoFac();
      }

      virtual void  v_SetGeoFac(MetricRelatedInfoSharedPtr minfo)
      {
		SetGeoFac(minfo);
      } 

      virtual void v_GetCoord(const double *Lcoords, double *coords)
      {
		GetCoord(Lcoords, coords);
      }

      virtual void v_GetCoords(double **coords)
      {
		GetCoords(coords);
      }

      virtual int GetCoordim()
      {
		return m_geom->GetCoordim();
      }

      virtual void v_WriteToFile(FILE *outfile)
      {
		WriteToFile(outfile);
      }

      virtual void v_WriteToFile(std::ofstream &outfile, const int dumpVar)
      {
		WriteToFile(outfile,dumpVar);
      }

      virtual StdRegions::GeomType v_GeoFacType()
      {
		return m_minfo->GetGtype();
      }

      /// \brief Virtual call to integrate the physical point list \a inarray
      /// over region (see SegExp::Integral) 
      virtual double v_Integral(const double *inarray )
      {
		return Integral(inarray);
      }

      /** \brief Virtual call to SegExp::IProduct_WRT_B */
      virtual void v_IProductWRTBase(const double * inarray, double * outarray)
      {
		IProductWRTBase(inarray,outarray);
      }


      /// Virtual call to SegExp::Deriv
      virtual void v_Deriv(double *outarray)
      {
		Deriv(1, this->m_phys, &outarray);
      }

      /// Virtual call to SegExp::Deriv
      virtual void v_StdDeriv(double * outarray)
      {
		StdSegExp::Deriv(this->m_phys, outarray);
      }

      /// Virtual call to SegExp::Deriv
      virtual void v_Deriv(const double *inarray, double * outarray)
      {
		Deriv(inarray, outarray);
      }

      /// Virtual call to SegExp::Deriv
      virtual void v_StdDeriv(const double *inarray, double * outarray)
      {
		StdSegExp::Deriv(inarray, outarray);
      }

      virtual void v_Deriv(const int n,  double ** outarray)
      {
		Deriv(n, outarray);
      }

      virtual void v_Deriv(const int n, const double *inarray, 
			   double **outarray)
      {
		Deriv(n, inarray, outarray);
      }

      /// Virtual call to SegExp::FwdTrans
      virtual void v_FwdTrans(const double * inarray)
      {
		FwdTrans(inarray);
      }

      /// Virtual call to SegExp::Evaluate
      virtual double v_Evaluate(const double * coords)
      {
		return Evaluate(coords);
      }

      /** \brief Virtual function to evaluate the discrete \f$ L_\infty\f$
	  error \f$ |\epsilon|_\infty = \max |u - u_{exact}|\f$ where \f$
	  u_{exact}\f$ is given by the array \a sol. 
	  
	  The full function is defined in StdExpansion::Linf 
	  
	  Input: 
	  
	  - \a _phys: takes the physical value space array as
          approximate solution

	  - \a sol: array of solution function  at physical quadrature points
	  
	  output: 
	  
	  - returns the \f$ L_\infty \f$ error as a double. 
      */
      virtual double v_Linf(const double *sol)
      {
	return Linf(sol);
      }

      /** \brief Virtual function to evaluate the \f$ L_\infty \f$ norm of
	  the function defined at the physical points \a (this)->_phys. 
	  
	  The full function is defined in StdExpansion::Linf 
	  
	  Input: 
	  
	  - \a _phys: uses the physical value space array as discrete
          function to be evaulated.
	  
	  output: 
	  
	  - returns the \f$ L_\infty \f$  as a double. 
      */
      virtual double v_Linf()
      {
	return Linf();
      }
    
      /** \brief Virtual function to evaluate the \f$ L_2\f$, \f$ |
	  \epsilon |_{2} = \left [ \int^1_{-1} [u - u_{exact}]^2 dx
	  \right]^{1/2} d\xi_1 \f$ where \f$ u_{exact}\f$ is given by the
	  array sol.
	  
	  The full function is defined in StdExpansion::L2 
	  
	  Input: 
	  
	  - \a _phys: takes the physical value space array as
          approximate solution
	  - \a sol: array of solution function  at physical quadrature points
	  
	  output: 
	  
	  - returns the \f$ L_2 \f$ error as a double. 
      */
      virtual double v_L2(const double *sol)
      {
	return StdExpansion::L2(sol);
      }

      /** \brief Virtual function to evaluate the \f$ L_2\f$ norm of the
	  function defined at the physical points \a (this)->_phys.  
	  
	  The full function is defined in StdExpansion::L2 
	
	  Input: 
	  
	  - \a _phys: uses the physical value space array as discrete
          function to be evaulated.
	  
	  output: 
	
	  - returns the \f$ L_2 \f$  as a double. 
      */
      virtual double v_L2()
      {
	return StdExpansion::L2();
      }
    };

    // type defines for use of SegExp in a boost vector
    typedef boost::shared_ptr<SegExp> SegExpSharedPtr;
    typedef std::vector< SegExpSharedPtr > SegExpVector;
    typedef std::vector< SegExpSharedPtr >::iterator SegExpVectorIter;
    
  } //end of namespace
} //end of namespace

#endif // SEGEXP_H

//
// $Log: SegExp.h,v $
// Revision 1.3  2006/05/30 14:00:04  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.2  2006/05/29 17:05:49  sherwin
// Modified to put shared_ptr around geom definitions
//
// Revision 1.1  2006/05/04 18:58:46  kirby
// *** empty log message ***
//
// Revision 1.33  2006/03/13 18:20:33  sherwin
//
// Fixed error in ResetGmat
//
// Revision 1.32  2006/03/12 21:59:48  sherwin
//
// compiling version of LocalRegions
//
//
