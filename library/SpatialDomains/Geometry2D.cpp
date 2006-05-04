////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/Geometry2D.cpp,v $
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific 
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:  
//
//
////////////////////////////////////////////////////////////////////////////////
#include "pchSpatialDomains.h"

#include <SpatialDomains/Geometry2D.h>

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>

namespace Nektar
{
  namespace SpatialDomains
  {
        Geometry2D::Geometry2D()
	{
        };

        Geometry2D::Geometry2D(const int coordim)
	{
	  m_coordim = coordim;
        };


#if 0
        Geometry2D::Geometry2D(const int order0, const int nquad0, 
			       const int order1, const int nquad1){
	  ASSERTL1(shape == StdRegions::eTriangle || 
		   shape == StdRegions::eQuadrilateral,
		   "Geometry2D::Geometry2D","Shape should be a triangle "
		   "or quadrilateral");

	  if(Shape == StdRegions::eTriangle){
	    const StdRegions::BasisKey B0(StdRegions::eModified_A, order0,
					  StdRegions::eLobatto, nquad0,0,0);
	    const StdRegions::BasisKey B1(StdRegions::eModified_B, order1,
					  StdRegions::eRadauM, nquad1,1,0);
	    // this format is not allowed in ansi - follow Geometry1D
	    m_xmap[0] = new StdRegions::StdTriExp [dim](B0,B1);
	  }
	  else if (Shape == StdRegions::eQuadrilateral)
	  {
	    const StdRegions::BasisKey B0(StdRegions::eModified_A, order0,
					  StdRegions::eLobatto, nquad0,0,0);
	    const StdRegions::BasisKey B1(StdRegions::eModified_A, order1,
					  StdRegions::eLobatto, nquad1,0,0);
	    // this format is not allowed in ansi - follow Geometry1D
	    m_xmap[0] = new StdRegions::StdTriExp [dim](B0,B1);
	  }
	  else
	  {
	    ASSERTL0(0,"Geometry2D::Geometry2D","Shape should be a "
		     "triangle or quadrilateral");
	  }

	  for(int i = 1; i < dim; ++i){
	    m_xmap[i] = m_xmap[i-1]+1;
	  }

	  m_state = eNotFilled;
        };
#endif


        Geometry2D::~Geometry2D()
	{
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: Geometry2D.cpp,v $
// Revision 1.14  2006/04/09 02:08:34  jfrazier
// Added precompiled header.
//
// Revision 1.13  2006/03/13 18:04:07  sherwin
//
// Corrected silly error in calling new
//
// Revision 1.12  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.11  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.10  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
