///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/LocalRegions/TetExp.h,v $ 
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TETEXP_H
#define TETEXP_H

#include <StdRegions/StdTetExp.h>

#include <SpatialDomains/TetGeom.h>

namespace Nektar
{
  namespace LocalRegions 
  {

    class TetExp: public StdRegions::StdTetExp
    {

    public:
    
      /** \brief Constructor using BasisKey class for quadrature
	  points and order definition */
      TetExp(const StdRegions::BasisKey &Ba, const StdRegions::BasisKey &Bb, 
	     const StdRegions::BasisKey &Bc);

      /** \brief Constructor using BasisKey class for quadrature points
	  and order definition where _coeffs and _phys are all set. */
      TetExp(const StdRegions::BasisKey &Ba, const StdRegions::BasisKey &Bb, 
	     const StdRegions::BasisKey &Bc, double *coeffs, double *phys);
      
      /// Copy Constructor
      TetExp(TetExp &T);
    
      /// Destructor
      ~TetExp();
    
      /// Return Shape of region, using  ShapeType enum list. i.e. Tetrahedron
      StdRegions::ShapeType DetShapeType() 
      { 
		return StdRegions::eTetrahedron; 
      }
    
    protected:
      int m_id;
      int m_field;
      
      SpatialDomains::TetGeom * m_geom;

    private:

      virtual StdRegions::ShapeType v_DetShapeType() 
      {
		return DetShapeType();
      }  
    };
    
    // type defines for use of TetExp in a boost vector
    typedef boost::shared_ptr<TetExp> TetExpSharedPtr;
    typedef std::vector< TetExpSharedPtr > TetExpVector;
    typedef std::vector< TetExpSharedPtr >::iterator TetExpVectorIter;

  } //end of namespace
} //end of namespace

#endif // TETEXP_H

/** 
 *    $Log: TetExp.h,v $
 *    Revision 1.2  2006/05/30 14:00:04  sherwin
 *    Updates to make MultiRegions and its Demos work
 *
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.9  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.8  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
