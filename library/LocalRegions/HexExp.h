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
#include <LocalRegions/QuadExp.h>

#include <SpatialDomains/GeomFactors.h>

#include <LocalRegions/MatrixKey.h>

#include <LocalRegions/Expansion3D.h>

namespace Nektar
{
    namespace LocalRegions 
    {

 
        class HexExp: public StdRegions::StdHexExp, public Expansion3D
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

            /// Return the region shape using the enum-list of ShapeType
            StdRegions::ExpansionType DetExpansionType() const
            { 
                return StdRegions::eHexahedron; 
            }

            const SpatialDomains::GeometrySharedPtr GetGeom() const
            {
                return m_geom;
            }

            const SpatialDomains::Geometry3DSharedPtr& GetGeom3D() const
            {
                return m_geom;
            }

            //------------------------------
            //    Integration Method
            //------------------------------
     
            /// \brief Integrate the physical point list \a inarray over region
            NekDouble Integral(const Array<OneD, const NekDouble> &inarray);
 
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(m_base[0]->GetBdata(), m_base[1]->GetBdata(), m_base[2]->GetBdata(), inarray,outarray,1);
            }
 

            void GetCoords(Array<OneD,NekDouble> &coords_1,
                           Array<OneD,NekDouble> &coords_2, 
                           Array<OneD,NekDouble> &coords_3);
            void GetCoord(const Array<OneD, const NekDouble> &Lcoords, Array<OneD,NekDouble> &coords);
      
            void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true);

            //-----------------------------
            // Differentiation Methods
            //-----------------------------

            void PhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                           Array<OneD, NekDouble> &out_d0,
                           Array<OneD, NekDouble> &out_d1,
                           Array<OneD, NekDouble> &out_d2);

            //----------------------------
            // Evaluations Methods
            //---------------------------
            /** \brief Forward transform from physical quadrature space
                stored in \a inarray and evaluate the expansion coefficients and
                store in \a (this)->_coeffs  */
            void FwdTrans(const Array<OneD, const NekDouble> & inarray, Array<OneD,NekDouble> &outarray);
            NekDouble PhysEvaluate(const Array<OneD, const NekDouble> &coord);

        protected:
            DNekMatSharedPtr GenMatrix(const StdRegions::StdMatrixKey &mkey);
            void GenMetricInfo();

            DNekMatSharedPtr CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekScalMatSharedPtr  CreateMatrix(const MatrixKey &mkey);
            DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);
//             DNekBlkMatSharedPtr CreateStdStaticCondMatrix(const StdRegions::StdMatrixKey &mkey);


            /** 
                \brief Calculate the inner product of inarray with respect to
                the basis B=base0*base1*base2 and put into outarray:
        
                \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
                \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
                \psi_{p}^{a} (\xi_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{r}^{a} (\xi_{3k})
                w_i w_j w_k u(\xi_{1,i} \xi_{2,j} \xi_{3,k})         
                J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\xi_{1,i})
                \sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2} \psi_{r}^a u(\xi_{1i},\xi_{2j},\xi_{3k})
                J_{i,j,k} \end{array} \f$ \n
        
                where
        
                \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a ( \xi_1) \psi_{q}^a (\xi_2) \psi_{r}^a (\xi_3) \f$ \n
        
                which can be implemented as \n
                \f$f_{r} (\xi_{3k}) = \sum_{k=0}^{nq_3} \psi_{r}^a u(\xi_{1i},\xi_{2j},\xi_{3k})
                J_{i,j,k} = {\bf B_3 U}   \f$ \n
                \f$ g_{q} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{r} (\xi_{3k})  = {\bf B_2 F}  \f$ \n
                \f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k})  = {\bf B_1 G} \f$

            **/
            void IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
                                 const Array<OneD, const NekDouble>& base1, 
                                 const Array<OneD, const NekDouble>& base2, 
                                 const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> & outarray, 
                                 int coll_check);

        private:
            SpatialDomains::Geometry3DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            HexExp();


            /// Return Shape of region, using  ShapeType enum list. i.e. Hexahedron
            virtual StdRegions::ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            }
    
            virtual const SpatialDomains::GeomFactorsSharedPtr& v_GetMetricInfo() const
            {
                return m_metricinfo;
            }

            virtual const SpatialDomains::GeometrySharedPtr v_GetGeom() const
            {
                return GetGeom();
            }

            virtual const SpatialDomains::Geometry3DSharedPtr& v_GetGeom3D() const
            {
                return GetGeom3D();
            }

            virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                     Array<OneD, NekDouble> &coords_1,
                                     Array<OneD, NekDouble> &coords_2)
            {
                GetCoords(coords_0, coords_1, coords_2);
            }

            virtual void v_GetCoord(const Array<OneD, const NekDouble> &lcoord, Array<OneD, NekDouble> &coord)
            {
                GetCoord(lcoord, coord);
            }

            virtual  int v_GetCoordim()
            {
                return m_geom->GetCoordim();
            }

            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true)
            {
                WriteToFile(outfile,format,dumpVar);
            }

            /** \brief Virtual call to integrate the physical point list \a inarray
                over region (see SegExp::Integral) */
            virtual NekDouble v_Integral(const Array<OneD, const NekDouble> &inarray )
            {
                return Integral(inarray);
            }

            /** \brief Virtual call to QuadExp::IProduct_WRT_B */
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray,outarray);
            }
        
            virtual void v_PhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1,
                                     Array<OneD, NekDouble> &out_d2)
            {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
            }

            /// Virtual call to SegExp::FwdTrans
            virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray,outarray);
            }
    
            /// Virtual call to QuadExp::Evaluate
            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble> &coords)
            {
                return PhysEvaluate(coords);
            }

            virtual NekDouble v_Linf()
            {
                return Linf();
            }
    
            virtual NekDouble v_L2(const Array<OneD, const NekDouble> &sol)
            {
                return StdExpansion::L2(sol);
            }

    
            virtual NekDouble v_L2()
            {
                return StdExpansion::L2();
            }

            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
            {
                return CreateStdMatrix(mkey);
            }

            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const MatrixKey &mkey)
            {
                return m_matrixManager[mkey];
            }
        
            virtual DNekScalBlkMatSharedPtr& v_GetLocStaticCondMatrix(const MatrixKey &mkey)
            {
                return m_staticCondMatrixManager[mkey];
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
 *    Revision 1.24  2008/09/23 18:20:25  pvos
 *    Updates for working ProjectContField3D demo
 *
 *    Revision 1.23  2008/09/17 17:30:28  ehan
 *    Fixed some errors to test the LocHexDemo.
 *
 *    Revision 1.22  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *    Revision 1.21  2008/07/29 22:25:34  sherwin
 *    general update for DG Advection including separation of GetGeom() into GetGeom1D,2D,3D()
 *
 *    Revision 1.20  2008/07/12 17:27:07  sherwin
 *    Update for AddBoundaryInt and moved various members to be private rather than protected
 *
 *    Revision 1.19  2008/07/04 10:19:04  pvos
 *    Some updates
 *
 *    Revision 1.18  2008/06/07 00:23:04  ehan
 *    Fixed warnings.
 *
 *    Revision 1.17  2008/06/06 23:24:03  ehan
 *    Added doxygen documentation
 *
 *    Revision 1.16  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.15  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.14  2008/04/06 05:59:04  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.13  2008/03/19 06:53:08  ehan
 *    Fixed recent changes of call by reference for the matrix shared pointer. Also fixed name of old functions from Get* to Create*.
 *
 *    Revision 1.12  2008/03/17 10:36:00  pvos
 *    Clean up of the code
 *
 *    Revision 1.11  2008/03/12 15:24:29  pvos
 *    Clean up of the code
 *
 *    Revision 1.10  2008/02/16 05:49:47  ehan
 *    Added PhysDeriv, GenMatrixInfo, standard matrix, and virtual functions.
 *
 *    Revision 1.9  2008/02/05 00:36:23  ehan
 *    Removed old comments
 *
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
