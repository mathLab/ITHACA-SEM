///////////////////////////////////////////////////////////////////////////////
//
// File TriExp.h
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

#ifndef TRIEXP_H
#define TRIEXP_H

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdTriExp.h>
#include <SpatialDomains/TriGeom.h>

#include <SpatialDomains/GeomFactors.h>

#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
    namespace LocalRegions
    {

    class TriExp: public StdRegions::StdTriExp
    {

    public:

        /** \brief Constructor using BasisKey class for quadrature
        points and order definition */
        TriExp(const LibUtilities::BasisKey &Ba,
                const LibUtilities::BasisKey &Bb,
                const SpatialDomains::TriGeomSharedPtr &geom);

        TriExp(const LibUtilities::BasisKey &Ba,
               const LibUtilities::BasisKey &Bb);


        /// Copy Constructor
        TriExp(const TriExp &T);

        /// Destructor
        ~TriExp();

        void GetCoords(Array<OneD,NekDouble> &coords_1,
                       Array<OneD,NekDouble> &coords_2, 
                       Array<OneD,NekDouble> &coords_3 = NullNekDouble1DArray);
        void GetCoord(const ConstArray<OneD,NekDouble>& Lcoords, 
                      Array<OneD,NekDouble> &coords);

        SpatialDomains::TriGeomSharedPtr GetGeom()
        {
            return m_geom;
        }

        void WriteToFile(FILE *outfile);
        void WriteToFile(std::ofstream &outfile, const int dumpVar);

        //----------------------------
        // Integration Methods
        //----------------------------

        /// \brief Integrate the physical point list \a inarray over region
        NekDouble Integral(const ConstArray<OneD, NekDouble> &inarray);

        /** \brief  Inner product of \a inarray over region with respect to the
        expansion basis (this)->_Base[0] and return in \a outarray */
        void IProductWRTBase(const ConstArray<OneD, NekDouble> &inarray, 
                             Array<OneD, NekDouble> &outarray);

        //-----------------------------
        // Differentiation Methods
        //-----------------------------

        void PhysDeriv(const ConstArray<OneD, NekDouble> &inarray, 
                       Array<OneD, NekDouble> &out_d0,
                       Array<OneD, NekDouble> &out_d1,
                       Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);  

        //----------------------------
        // Evaluations Methods
        //---------------------------

        /** \brief Forward transform from physical quadrature space
        stored in \a inarray and evaluate the expansion coefficients and
        store in \a (this)->_coeffs  */
        void FwdTrans(const ConstArray<OneD, NekDouble> &inarray, 
                      Array<OneD, NekDouble> &outarray);
        
        NekDouble PhysEvaluate(const ConstArray<OneD, NekDouble> &coord);
        
    protected:

        void GenMetricInfo();

        DNekMatSharedPtr GetStdMatrix(const StdRegions::StdMatrixKey &mkey);

        DNekScalMatSharedPtr    CreateMatrix(const MatrixKey &mkey);

        SpatialDomains::TriGeomSharedPtr m_geom;
        SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

        LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;

        /** \brief  Inner product of \a inarray over region with respect to
        the expansion basis \a base and return in \a outarray */
        inline void IProductWRTBase(const ConstArray<OneD, NekDouble> &base0, 
                                    const ConstArray<OneD, NekDouble> &base1, 
                                    const ConstArray<OneD, NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &outarray);

    private:

        virtual StdRegions::ShapeType v_DetShapeType()
        {
			return DetShapeType();
        }

        virtual SpatialDomains::GeomFactorsSharedPtr v_GetMetricInfo() const
        {
            return m_metricinfo;
        }

         virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                 Array<OneD, NekDouble> &coords_1 = NullNekDouble1DArray,
                                 Array<OneD, NekDouble> &coords_2 = NullNekDouble1DArray)
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
        virtual void v_IProductWRTBase(const ConstArray<OneD, NekDouble> &inarray,
                                       Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);
        }

         /// Virtual call to SegExp::PhysDeriv
        virtual void v_StdPhysDeriv(const ConstArray<OneD, NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &out_d0,
                                    Array<OneD, NekDouble> &out_d1)
        {
            StdTriExp::PhysDeriv(inarray, out_d0, out_d1);
        }
        
        virtual void v_PhysDeriv(const ConstArray<OneD, NekDouble> &inarray, 
                                 Array<OneD, NekDouble> &out_d0,
                                 Array<OneD, NekDouble> &out_d1,
                                 Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
	 {
             PhysDeriv(inarray, out_d0, out_d1);
         }

        /// Virtual call to SegExp::FwdTrans
        virtual void v_FwdTrans(const ConstArray<OneD, NekDouble> &inarray, 
                                Array<OneD, NekDouble> &outarray)
        {
            FwdTrans(inarray,outarray);
        }
        
        /// Virtual call to TriExp::Evaluate
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

        virtual DNekScalMatSharedPtr v_GetLocMatrix(StdRegions::MatrixType type)
        {
            MatrixKey masskey(type,DetShapeType(),*this);
            return m_matrixManager[masskey];
        }

    };

    // type defines for use of TriExp in a boost vector
    typedef boost::shared_ptr<TriExp> TriExpSharedPtr;
    typedef std::vector< TriExpSharedPtr > TriExpVector;
    typedef std::vector< TriExpSharedPtr >::iterator TriExpVectorIter;

    } //end of namespace
} //end of namespace

#endif // TRIEXP_H

/**
 *    $Log: TriExp.h,v $
 *    Revision 1.14  2007/07/10 17:17:26  sherwin
 *    Introduced Scaled Matrices into the MatrixManager
 *
 *    Revision 1.13  2007/06/07 15:54:19  pvos
 *    Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 *    Also made corrections to various ASSERTL2 calls
 *
 *    Revision 1.12  2007/06/01 17:08:07  pvos
 *    Modification to make LocalRegions/Project2D run correctly (PART1)
 *
 *    Revision 1.11  2007/05/31 19:13:12  pvos
 *    Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 *    Revision 1.10  2007/05/31 11:38:17  pvos
 *    Updated QuadExp and TriExp
 *
 *    Revision 1.9  2007/01/15 21:12:26  sherwin
 *    First definition
 *
 *    Revision 1.8  2006/12/10 18:59:47  sherwin
 *    Updates for Nodal points
 *
 *    Revision 1.7  2006/06/13 18:05:01  sherwin
 *    Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
 *    Revision 1.6  2006/06/05 00:08:48  bnelson
 *    Fixed a gcc 4.1.0 compilation problem.  TriExp::GenGeoFac not allowed in the class declaration, but GenGeoFac is.
 *
 *    Revision 1.5  2006/06/02 18:48:39  sherwin
 *    Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
 *
 *    Revision 1.4  2006/06/01 14:15:58  sherwin
 *    Added typdef of boost wrappers and made GeoFac a boost shared pointer.
 *
 *    Revision 1.3  2006/05/30 14:00:04  sherwin
 *    Updates to make MultiRegions and its Demos work
 *
 *    Revision 1.2  2006/05/29 17:05:49  sherwin
 *    Modified to put shared_ptr around geom definitions
 *
 *    Revision 1.1  2006/05/04 18:58:47  kirby
 *    *** empty log message ***
 *
 *    Revision 1.13  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.12  2006/03/12 07:43:33  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
