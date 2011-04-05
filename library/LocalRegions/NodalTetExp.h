///////////////////////////////////////////////////////////////////////////////
//
// File NodalTetExp.h
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
// Description: Header for NodalTetExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTETEXP_H
#define NODALTETEXP_H

#include <LocalRegions/TetExp.h>
#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdTetExp.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/GeomFactors.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/LocalRegionsDeclspec.h>


namespace Nektar
{
  namespace LocalRegions 
  {
  
    class NodalTetExp: public TetExp
    {
//     public:
//       NodalTetExp(const StdRegions::Basis &Ba, const StdRegions::Basis &Bb, 
//           const StdRegions::Basis &Bc):
//     TetExp(Ba,Bb,Bc)
//       {
//       }
     	/** \brief Constructor using BasisKey class for quadrature
        points and order definition */
        LOCAL_REGIONS_EXPORT NodalTetExp(const LibUtilities::BasisKey &Ba,
                    const LibUtilities::BasisKey &Bb,
	            const LibUtilities::BasisKey &Bc,
                    const SpatialDomains::TetGeomSharedPtr &geom);

        LOCAL_REGIONS_EXPORT NodalTetExp(const LibUtilities::BasisKey &Ba,
	            const LibUtilities::BasisKey &Bb,
	            const LibUtilities::BasisKey &Bc);
	    
      /// Copy Constructor
          LOCAL_REGIONS_EXPORT NodalTetExp(const NodalTetExp &T);

      /// Destructor
          LOCAL_REGIONS_EXPORT ~NodalTetExp();
    
    protected:

    private:
            
    };

  } //end of namespace
} //end of namespace

#endif // NODALTETEXP_H


/** 
 *    $Log: NodalTetExp.h,v $
 *    Revision 1.2  2007/07/20 00:45:50  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.1  2006/05/04 18:58:45  kirby
 *    *** empty log message ***
 *
 *    Revision 1.7  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.6  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
