///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/LocalRegions/TetExp.h,v $ 
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

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdTetExp.h>
#include <SpatialDomains/TetGeom.h>

#include <SpatialDomains/GeomFactors.h>
#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
  namespace LocalRegions 
  {

    class TetExp: public StdRegions::StdTetExp
    {

    public:

	/** \brief Constructor using BasisKey class for quadrature
        points and order definition */
        TetExp(const LibUtilities::BasisKey &Ba,
               const LibUtilities::BasisKey &Bb,
	       const LibUtilities::BasisKey &Bc,
               const SpatialDomains::TetGeomSharedPtr &geom);

        TetExp(const LibUtilities::BasisKey &Ba,
	       const LibUtilities::BasisKey &Bb,
	       const LibUtilities::BasisKey &Bc);

	    
      /// Copy Constructor
      TetExp(const TetExp &T);

      /// Destructor
      ~TetExp();

        void GetCoords(Array<OneD,NekDouble> &coords_0,
		       Array<OneD,NekDouble> &coords_1,
          	       Array<OneD,NekDouble> &coords_2);

	void GetCoord(const ConstArray<OneD,NekDouble> &Lcoords, Array<OneD,NekDouble> &coords);

        //----------------------------
        // Integration Methods
        //----------------------------

        /// \brief Integrate the physical point list \a inarray over region
       NekDouble Integral(const ConstArray<OneD,NekDouble> &inarray);

       void IProductWRTBase(const ConstArray<OneD,NekDouble> &inarray, Array<OneD,NekDouble> &outarray);
       void FwdTrans(const ConstArray<OneD,NekDouble> & inarray,Array<OneD,NekDouble> &outarray);

       NekDouble PhysEvaluate(const ConstArray<OneD,NekDouble> &coord);

      //-----------------------------
      // Differentiation Methods
      //-----------------------------

        void PhysDeriv(const ConstArray<OneD, NekDouble> &inarray, 
                       Array<OneD, NekDouble> &out_d0,
                       Array<OneD, NekDouble> &out_d1,
                       Array<OneD, NekDouble> &out_d2);


	/// Return Shape of region, using  ShapeType enum list. i.e. Tetrahedron
	StdRegions::ShapeType DetShapeType() 
	{ 
	return StdRegions::eTetrahedron; 
	}

        SpatialDomains::TetGeomSharedPtr GetGeom()
        {
            return m_geom;
        }

        void WriteToFile(FILE *outfile);
	void WriteToFile(std::ofstream &outfile, const int dumpVar);


    protected:

        void GenMetricInfo();

	void IProductWRTBase(const ConstArray<OneD,NekDouble> &base0, 
			     const ConstArray<OneD,NekDouble> &base1, 
			     const ConstArray<OneD,NekDouble> &base2, 
			     const ConstArray<OneD,NekDouble> &inarray,
			     Array<OneD,NekDouble> &outarray);


	DNekMatSharedPtr GetStdMatrix(const StdRegions::StdMatrixKey &mkey);
        DNekScalMatSharedPtr    CreateMatrix(const MatrixKey &mkey);

        DNekBlkMatSharedPtr GetStdStaticCondMatrix(const StdRegions::StdMatrixKey &mkey);
        DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);

	SpatialDomains::TetGeomSharedPtr m_geom;
        SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

	LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
        LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;


    private:
      TetExp();
	
      virtual StdRegions::ShapeType v_DetShapeType() const
      {
	 return DetShapeType();
      }

      virtual SpatialDomains::GeomFactorsSharedPtr v_GetMetricInfo() const
      {
         return m_metricinfo;
      }

//          virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
//                                  Array<OneD, NekDouble> &coords_1 = NullNekDouble1DArray,
//                                  Array<OneD, NekDouble> &coords_2 = NullNekDouble1DArray)
//          {
//              GetCoords(coords_0, coords_1, coords_2);
//          }

        virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                 Array<OneD, NekDouble> &coords_1,
                                 Array<OneD, NekDouble> &coords_2)
        {
            GetCoords(coords_0, coords_1, coords_2);
        }

        virtual void v_GetCoord(const ConstArray<OneD, NekDouble> &lcoord, 
                                Array<OneD, NekDouble> &coord)
        {
            GetCoord(lcoord, coord);
        }

        virtual int v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }

                 /// Virtual call to SegExp::PhysDeriv
        virtual void v_StdPhysDeriv(const ConstArray<OneD, NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &out_d0,
                                    Array<OneD, NekDouble> &out_d1,
                                    Array<OneD, NekDouble> &out_d2)
        {
            StdTetExp::PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }

        virtual void v_PhysDeriv(const ConstArray<OneD, NekDouble> &inarray,
                                 Array<OneD, NekDouble> &out_d0,
                                 Array<OneD, NekDouble> &out_d1,
                                 Array<OneD, NekDouble> &out_d2)
        {
            PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }

        virtual void v_WriteToFile(FILE *outfile)
        {
            WriteToFile(outfile);
        }

        virtual void v_WriteToFile(std::ofstream &outfile, const int dumpVar)
        {
            WriteToFile(outfile,dumpVar);
        }

        /** \brief Virtual call to integrate the physical point list \a inarray
        over region (see SegExp::Integral) */
        virtual NekDouble v_Integral(const ConstArray<OneD, NekDouble> &inarray )
        {
            return Integral(inarray);
        }

        /** \brief Virtual call to TriExp::IProduct_WRT_B */
        virtual void v_IProductWRTBase(const ConstArray<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);
        }

        virtual void v_IProductWRTDerivBase (const int dir,
                                             const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD, NekDouble> &outarray)
        {
            IProductWRTDerivBase(dir,inarray,outarray);
        }

        /// Virtual call to SegExp::FwdTrans
        virtual void v_FwdTrans(const ConstArray<OneD, NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray)
        {
            FwdTrans(inarray, outarray);
        }

        /// Virtual call to TetExp::Evaluate
        virtual NekDouble v_PhysEvaluate(const ConstArray<OneD, NekDouble> &coords)
        {
            return PhysEvaluate(coords);
        }

        virtual NekDouble v_Linf(const ConstArray<OneD, NekDouble> &sol)
        {
            return Linf(sol);
        }

        virtual NekDouble v_Linf()
        {
            return Linf();
        }

        virtual NekDouble v_L2(const ConstArray<OneD, NekDouble> &sol)
        {
            return StdExpansion::L2(sol);
        }


        virtual NekDouble v_L2()
        {
            return StdExpansion::L2();
        }

        virtual DNekScalMatSharedPtr v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }

        virtual DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
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
 *    Revision 1.7  2008/02/05 00:41:37  ehan
 *    Added initial tetrahedral expansion.
 *
 *    Revision 1.6  2007/07/22 23:04:19  bnelson
 *    Backed out Nektar::ptr.
 *
 *    Revision 1.5  2007/07/20 00:45:51  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.4  2007/01/15 21:12:26  sherwin
 *    First definition
 *
 *    Revision 1.3  2006/06/13 18:05:01  sherwin
 *    Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
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
