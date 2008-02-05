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


#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdHexExp.h>
#include <SpatialDomains/HexGeom.h>

#include <SpatialDomains/GeomFactors.h>

#include <LocalRegions/MatrixKey.h>


namespace Nektar
{
  namespace LocalRegions 
  {

 
    class HexExp: public StdRegions::StdHexExp
    {
    public:
    
      ///\brief Constructor using BasisKey class for quadrature
      /// points and order definition 
 
      HexExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb,
             const LibUtilities::BasisKey &Bc, const SpatialDomains::HexGeomSharedPtr &geom);

      HexExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, 
             const LibUtilities::BasisKey &Bc);

      /// Copy Constructor
      HexExp(const HexExp &T);

      /// Destructor
      ~HexExp();

      StdRegions::ShapeType DetShapeType() const
      { 
    	 return StdRegions::eHexahedron; 
      }

      //------------------------------
      //    Integration Method
      //------------------------------
      NekDouble Integral(const ConstArray<OneD,NekDouble> &inarray);

      void IProductWRTBase(const ConstArray<OneD,NekDouble> &inarray, Array<OneD,NekDouble> &outarray);

      void IProductWRTBase(const ConstArray<OneD,NekDouble> &base0, 
			   const ConstArray<OneD,NekDouble> &base1,
			   const ConstArray<OneD,NekDouble> &base2,  
			   const ConstArray<OneD,NekDouble> &inarray,
			   Array<OneD,NekDouble> &outarray );
      void FwdTrans(const ConstArray<OneD,NekDouble> & inarray, Array<OneD,NekDouble> &outarray);
 
      void GetCoord(const ConstArray<OneD,NekDouble> &Lcoords, Array<OneD,NekDouble> &coords);

      NekDouble PhysEvaluate(const ConstArray<OneD,NekDouble> &coord);




    protected:
       void GenMetricInfo();

        DNekMatSharedPtr GetStdMatrix(const StdRegions::StdMatrixKey &mkey);
        DNekScalMatSharedPtr  CreateMatrix(const MatrixKey &mkey);

        DNekBlkMatSharedPtr GetStdStaticCondMatrix(const StdRegions::StdMatrixKey &mkey);
        DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);

	SpatialDomains::HexGeomSharedPtr m_geom;
        SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

        LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
        LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;


    private:
      /// Return Shape of region, using  ShapeType enum list. i.e. Hexahedron
       HexExp();

        virtual StdRegions::ShapeType v_DetShapeType() const
        {
            return DetShapeType();
        }
    
        virtual SpatialDomains::GeomFactorsSharedPtr v_GetMetricInfo() const
        {
            return m_metricinfo;
        }

      
    };

    // type defines for use of HexExp in a boost vector
    typedef boost::shared_ptr<HexExp> HexExpSharedPtr;
    typedef std::vector< HexExpSharedPtr > HexExpVector;
    typedef std::vector< HexExpSharedPtr >::iterator HexExpVectorIter;

    
  } //end of namespace
} //end of namespace

#endif //HEX_EXP_H

/** 
 *    $Log: HexExp.h,v $
 *    Revision 1.8  2008/01/31 10:57:06  ehan
 *    Implemented IProductWRTBase, FwdTrans, GetCoord, GetStdMatrix, and GetStdStaticCondMatrix.
 *
 *    Revision 1.7  2007/07/22 23:04:17  bnelson
 *    Backed out Nektar::ptr.
 *
 *    Revision 1.6  2007/07/20 00:45:50  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.5  2007/01/15 21:12:25  sherwin
 *    First definition
 *
 *    Revision 1.4  2006/06/13 18:05:01  sherwin
 *    Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
 *    Revision 1.3  2006/06/01 13:42:07  kirby
 *    *** empty log message ***
 *
 *    Revision 1.2  2006/05/30 14:00:03  sherwin
 *    Updates to make MultiRegions and its Demos work
 *
 *    Revision 1.1  2006/05/04 18:58:45  kirby
 *    *** empty log message ***
 *
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
