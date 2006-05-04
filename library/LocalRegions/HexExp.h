///////////////////////////////////////////////////////////////////////////////
//
// File HexExp.h
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
// Description: Header routines for Hex expansion 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef HEXEXP_H
#define HEXEXP_H

#include <StdRegions/StdBasis.h>
#include <StdRegions/StdHexExp.h>
#include <SpatialDomains/HexGeom.h>

#include <StdRegions/StdRegions.hpp>

#include <LibUtilities/NekMemoryManager.hpp>

namespace Nektar
{
  namespace LocalRegions 
  {

 
    class HexExp: public StdRegions::StdHexExp
    {
    public:
    
      ///\brief Constructor using BasisKey class for quadrature
      /// points and order definition 
      HexExp(const StdRegions::BasisKey &Ba, const StdRegions::BasisKey &Bb, 
	     const StdRegions::BasisKey &Bc);
    
      /// \brief Constructor using BasisKey class for quadrature
      /// points and order definition where _coeffs and _phys are all set. 
      HexExp(const StdRegions::BasisKey &Ba, const StdRegions::BasisKey &Bb, 
	     const StdRegions::BasisKey &Bc, double *coeffs, double *phys);
    
      /// Copy Constructor
      HexExp(HexExp &T);

      /// Destructor
      ~HexExp();

      StdRegions::ShapeType DetShapeType() 
      { 
	return StdRegions::eHexahedron; 
      }
    
    protected:
      int m_id;
      int m_field;
      
      SpatialDomains::HexGeom * m_geom;

    private:
      /// Return Shape of region, using  ShapeType enum list. i.e. Hexahedron
      virtual StdRegions::ShapeType V_DetShapeType() 
      {
	DetShapeType();
      }
      

  };
  
  } //end of namespace
} //end of namespace

#endif //HEX_EXP_H

/** 
 *    $Log: HexExp.h,v $
 *    Revision 1.14  2006/05/02 21:21:11  sherwin
 *    Corrected libraries to compile new version of spatialdomains and demo Graph1D
 *
 *    Revision 1.13  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.12  2006/03/12 07:43:31  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
