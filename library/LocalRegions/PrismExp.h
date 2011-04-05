///////////////////////////////////////////////////////////////////////////////
//
// File PrismExp.h
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
// Description: Header file for PrismExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef PRISMEXP_H
#define PRISMEXP_H

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdPrismExp.h>
#include <SpatialDomains/PrismGeom.h>

#include <SpatialDomains/GeomFactors3D.h>
#include <LocalRegions/MatrixKey.h>

#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

namespace Nektar
{
    namespace LocalRegions 
    {

        class PrismExp: public StdRegions::StdPrismExp, public Expansion3D
        {
        public:
    
            /// \brief Constructor using BasisKey class for quadrature
            /// points and order definition 
            LOCAL_REGIONS_EXPORT PrismExp(const LibUtilities::BasisKey &Ba,
                     const LibUtilities::BasisKey &Bb,
                     const LibUtilities::BasisKey &Bc,
                     const SpatialDomains::PrismGeomSharedPtr &geom);

	    
            /// Copy Constructor
            LOCAL_REGIONS_EXPORT PrismExp(const PrismExp &T);

            /// Destructor
            LOCAL_REGIONS_EXPORT ~PrismExp();


            LOCAL_REGIONS_EXPORT void GetCoords(Array<OneD,NekDouble> &coords_0,
                           Array<OneD,NekDouble> &coords_1,
                           Array<OneD,NekDouble> &coords_2);

            LOCAL_REGIONS_EXPORT void GetCoord(const Array<OneD, const NekDouble> &Lcoords, 
                          Array<OneD,NekDouble> &coords);

            //----------------------------
            // Integration Methods
            //----------------------------

            /// \brief Integrate the physical point list \a inarray over region
            LOCAL_REGIONS_EXPORT NekDouble Integral(const Array<OneD, const NekDouble> &inarray);            

            /** \brief  Inner product of \a inarray over region with respect to the 
                expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata() and return in \a outarray 
    
                Wrapper call to StdPrismExp::IProductWRTBase
    
                Input:\n
    
                - \a inarray: array of function evaluated at the physical collocation points
    
                Output:\n
    
                - \a outarray: array of inner product with respect to each basis over region

            */
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata(),inarray,outarray);
            }

            LOCAL_REGIONS_EXPORT void FwdTrans(const Array<OneD, const NekDouble> & inarray,Array<OneD,NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT NekDouble PhysEvaluate(const Array<OneD, const NekDouble> &coord);

            //-----------------------------
            // Differentiation Methods
            //-----------------------------

            LOCAL_REGIONS_EXPORT void PhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                           Array<OneD, NekDouble> &out_d0,
                           Array<OneD, NekDouble> &out_d1,
                           Array<OneD, NekDouble> &out_d2);

            /// Return Shape of region, using  ShapeType enum list. i.e. Prism
            StdRegions::ExpansionType DetExpansionType() const
            { 
                return StdRegions::ePrism; 
            }

            const SpatialDomains::GeometrySharedPtr GetGeom() const
            {
                return m_geom;
            }

            const SpatialDomains::Geometry3DSharedPtr& GetGeom3D() const
            {
                return m_geom;
            }

            LOCAL_REGIONS_EXPORT void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v");


        protected:
            DNekMatSharedPtr GenMatrix(const StdRegions::StdMatrixKey &mkey);

            /**
               \brief Calculate the inner product of inarray with respect to
               the basis B=base0*base1*base2 and put into outarray:
              
               \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
               \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
               \psi_{p}^{a} (\bar \eta_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{pr}^{b} (\xi_{3k})
               w_i w_j w_k u(\bar \eta_{1,i} \xi_{2,j} \xi_{3,k})      
               J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\bar \eta_{1,i})
               \sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2} \psi_{pr}^b u(\bar \eta_{1i},\xi_{2j},\xi_{3k})
               J_{i,j,k} \end{array} \f$ \n
            
               where
            
               \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\bar \eta_1) \psi_{q}^a (\xi_2) \psi_{pr}^b (\xi_3) \f$ \n
            
               which can be implemented as \n
               \f$f_{pr} (\xi_{3k}) = \sum_{k=0}^{nq_3} \psi_{pr}^b u(\bar \eta_{1i},\xi_{2j},\xi_{3k})
               J_{i,j,k} = {\bf B_3 U}   \f$ \n
               \f$ g_{q} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{pr} (\xi_{3k})  = {\bf B_2 F}  \f$ \n
               \f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k})  = {\bf B_1 G} \f$

            **/
            void IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
                                 const Array<OneD, const NekDouble>& base1, 
                                 const Array<OneD, const NekDouble>& base2, 
                                 const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> & outarray);

            DNekMatSharedPtr CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekScalMatSharedPtr  CreateMatrix(const MatrixKey &mkey);
            DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);

    
        private:
            SpatialDomains::Geometry3DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;


            PrismExp();

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
                                     Array<OneD, NekDouble> &coords_1 = NullNekDouble1DArray,
                                     Array<OneD, NekDouble> &coords_2 = NullNekDouble1DArray)
            {
                GetCoords(coords_0, coords_1, coords_2);
            }

            virtual void v_GetCoord(const Array<OneD, const NekDouble> &lcoord, 
                                    Array<OneD, NekDouble> &coord)
            {
                GetCoord(lcoord, coord);
            }

            virtual int v_GetCoordim()
            {
                return m_geom->GetCoordim();
            }

            virtual void v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1,
                                     Array<OneD, NekDouble> &out_d2)
            {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
            }

            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v")
            {
                WriteToFile(outfile,format,dumpVar,var);
            }

            /** \brief Virtual call to integrate the physical point list \a inarray
                over region (see SegExp::Integral) */
            virtual NekDouble v_Integral(const Array<OneD, const NekDouble> &inarray )
            {
                return Integral(inarray);
            }

            /** \brief Virtual call to TriExp::IProduct_WRT_B */
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            virtual void v_IProductWRTDerivBase (const int dir,
                                                 const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTDerivBase(dir,inarray,outarray);
            }

            /// Virtual call to SegExp::FwdTrans
            virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray, outarray);
            }
        
            /// Virtual call to PrismExp::Evaluate
            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble> &coords)
            {
                return PhysEvaluate(coords);
            }

            virtual NekDouble v_Linf(const Array<OneD, const NekDouble> &sol)
            {
                return Linf(sol);
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

        // type defines for use of PrismExp in a boost vector
        typedef boost::shared_ptr<PrismExp> PrismExpSharedPtr;
        typedef std::vector< PrismExpSharedPtr > PrismExpVector;
        typedef std::vector< PrismExpSharedPtr >::iterator PrismExpVectorIter;
    } //end of namespace
} //end of namespace

#define H_PRISMEXP
#endif

/** 
 *    $Log: PrismExp.h,v $
 *    Revision 1.20  2009/04/27 21:34:07  sherwin
 *    Updated WriteToField
 *
 *    Revision 1.19  2009/04/20 16:12:28  sherwin
 *    Updates related to output format and optimising DG solver
 *
 *    Revision 1.18  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *    Revision 1.17  2008/07/29 22:25:34  sherwin
 *    general update for DG Advection including separation of GetGeom() into GetGeom1D,2D,3D()
 *
 *    Revision 1.16  2008/07/12 17:27:07  sherwin
 *    Update for AddBoundaryInt and moved various members to be private rather than protected
 *
 *    Revision 1.15  2008/07/04 10:19:05  pvos
 *    Some updates
 *
 *    Revision 1.14  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.13  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.12  2008/04/06 05:59:05  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.11  2008/03/18 16:18:05  pvos
 *    Fixed some compiler complaints
 *
 *    Revision 1.10  2008/03/17 10:35:15  pvos
 *    Clean up of the code
 *
 *    Revision 1.9  2008/03/12 15:24:29  pvos
 *    Clean up of the code
 *
 *    Revision 1.8  2008/02/16 05:50:51  ehan
 *    Added PhysDeriv and virtual functions.
 *
 *    Revision 1.7  2008/02/05 00:39:04  ehan
 *    Added initial prismatic expansion.
 *
 *    Revision 1.6  2007/07/22 23:04:18  bnelson
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
 *    Revision 1.2  2006/05/30 14:00:03  sherwin
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
