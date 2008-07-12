///////////////////////////////////////////////////////////////////////////////
//
// File QuadExp.h
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

#ifndef QUADEXP_H
#define QUADEXP_H

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdQuadExp.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/SegGeom.h>
#include <LocalRegions/SegExp.h>

#include <SpatialDomains/GeomFactors.h>
#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp> 
#include <LibUtilities/BasicUtils/Vmath.hpp> 

namespace Nektar
{

    namespace LocalRegions
    {
        
    class QuadExp: public StdRegions::StdQuadExp
        {
        public:

            /** \brief Constructor using BasisKey class for quadrature
                points and order definition */
            QuadExp(const LibUtilities::BasisKey &Ba,
                    const LibUtilities::BasisKey &Bb,
                    const SpatialDomains::QuadGeomSharedPtr &geom);
        
            QuadExp(const LibUtilities::BasisKey &Ba,
                    const LibUtilities::BasisKey &Bb);

            /// \brief Constructor using BasisKey class for quadrature
            /// points and order definition where it has standard geometric factors 
            QuadExp(const LibUtilities::BasisKey &Ba);

            /// Copy Constructor
            QuadExp(const QuadExp &T);

            /// Destructor
            ~QuadExp();

            /// Return Shape of region, using ShapeType enum
            /// list. i.e. Quadrilateral
            StdRegions::ExpansionType DetExpansionType() const
            {
                return StdRegions::eQuadrilateral;
            }

            void GetCoords(Array<OneD,NekDouble> &coords_1,
                           Array<OneD,NekDouble> &coords_2, 
                           Array<OneD,NekDouble> &coords_3 = NullNekDouble1DArray);
            void GetCoord(const Array<OneD, const NekDouble>& Lcoords, 
                          Array<OneD,NekDouble> &coords);

            const SpatialDomains::Geometry2DSharedPtr& GetGeom()
            {
                return m_geom;
            }

            void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true);

            //----------------------------
            // Integration Methods
            //----------------------------

            /// \brief Integrate the physical point list \a inarray over region
            NekDouble Integral(const Array<OneD, const NekDouble> &inarray);

        
    
            /** \brief  Inner product of \a inarray over region with respect to the 
                expansion basis (this)->_Base[0] and return in \a outarray 
            
                Wrapper call to QuadExp::IProduct_WRT_B
        
                Input:\n
            
                - \a inarray: array of function evaluated at the physical
                collocation points
            
                Output:\n
            
                - \a outarray: array of inner product with respect to each
                basis over region
            
            */    
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                                inarray,outarray,1);
            }

            void IProductWRTDerivBase(const int dir,
                                      const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray);
        
            //-----------------------------
            // Differentiation Methods
            //-----------------------------

            void PhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                           Array<OneD, NekDouble> &out_d0,
                           Array<OneD, NekDouble> &out_d1,
                           Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);      
        
            void PhysDeriv(const int dir, 
                           const Array<OneD, const NekDouble>& inarray,
                           Array<OneD, NekDouble> &outarray);
        
            //----------------------------
            // Evaluations Methods
            //---------------------------
        
            /** \brief Forward transform from physical quadrature space
                stored in \a inarray and evaluate the expansion coefficients and
                store in \a (this)->_coeffs  */
            void FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                          Array<OneD, NekDouble> &outarray);
        
            NekDouble PhysEvaluate(const Array<OneD, const NekDouble> &coord);        
        
            SegExpSharedPtr GetEdgeExp(int edge);
            void GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray);
    
            void AddNormTraceInt(const int dir,
                                 Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,NekDouble> &outarray)
            {
                AddNormBoundaryInt(dir,inarray,outarray,true);
            }

            void AddNormBoundaryInt(const int dir,
                                    Array<OneD, const NekDouble> &inarray,
                                    Array<OneD,NekDouble> &outarray,
                                    bool InArrayIsTrace = false);

            void AddBoundaryInt(Array<OneD, const NekDouble> &inarray,
                                Array<OneD,NekDouble> &outarray,
                                bool InArrayIsTrace = true);

            void AddUDGHelmholtzBoundaryTerms(const NekDouble tau,
                                              const Array<OneD,const NekDouble> &inarray,
                                              Array<OneD,NekDouble> &outarray);
            
            void AddEdgeBoundaryInt(const int edge, SegExpSharedPtr &EdgeExp,
                                    Array <OneD,NekDouble > &outarray);          


            void AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray,
                                           bool InputDataIsCartesianOrient = true);

            void AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,SegExpSharedPtr> &EdgeExp,
                                           Array<OneD,NekDouble> &outarray);

            void AddUDGHelmholtzEdgeTerms(const NekDouble tau, 
                                          const int edge,
                                          SegExpSharedPtr &EdgeExp, 
                                          Array <OneD,NekDouble > &outarray);

            DNekMatSharedPtr GenMatrix(const StdRegions::StdMatrixKey &mkey);

            void LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray);

            void HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const double lambda);

        protected:

            void GenMetricInfo();

            DNekMatSharedPtr CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekScalMatSharedPtr  CreateMatrix(const MatrixKey &mkey);
            DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);

        
            /** 
                \brief Calculate the inner product of inarray with respect to
                the basis B=base0*base1 and put into outarray:
            
                \f$ \begin{array}{rcl} I_{pq} = (\phi_q \phi_q, u) & = &
                \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \phi_p(\xi_{0,i})
                \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} \xi_{1,j})
                J_{i,j}\\ & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
                \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}
                J_{i,j} \end{array} \f$
            
                where
            
                \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
            
                which can be implemented as
            
                \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j} = 
                {\bf B_1 U}  \f$
                \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} = 
                {\bf B_0 F}  \f$
            **/
            void IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
                                 const Array<OneD, const NekDouble>& base1,
                                 const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray,
                                 int coll_check);

        private:
            SpatialDomains::Geometry2DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            QuadExp();

            virtual StdRegions::ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            }

            virtual SpatialDomains::GeomFactorsSharedPtr v_GetMetricInfo() const
            {
                return m_metricinfo;
            }

            virtual const SpatialDomains::Geometry2DSharedPtr& v_GetGeom()
            {
                return GetGeom();
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

            virtual void v_IProductWRTDerivBase (const int dir,
                                                 const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTDerivBase(dir,inarray,outarray);
            }
        
            virtual void v_PhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1,
                                     Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                PhysDeriv(inarray, out_d0, out_d1);
            }

            virtual void v_PhysDeriv(const int dir, 
                                     const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &outarray)
            {
                PhysDeriv(dir,inarray,outarray);
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

            virtual DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
            {
                return GenMatrix(mkey);
            }
        
            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
            {
                return CreateStdMatrix(mkey);
            }

            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const MatrixKey &mkey)
            {
                return m_matrixManager[mkey];
            }
        
            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const StdRegions::MatrixType mtype, NekDouble lambdaval, NekDouble tau)
            {
                MatrixKey mkey(mtype,DetExpansionType(),*this,lambdaval,tau);
                return m_matrixManager[mkey];
            }

            virtual DNekScalBlkMatSharedPtr& v_GetLocStaticCondMatrix(const MatrixKey &mkey)
            {
                return m_staticCondMatrixManager[mkey];
            }


            virtual StdRegions::EdgeOrientation v_GetEorient(int edge)
            {
                return m_geom->GetEorient(edge);
            }


            virtual StdRegions::EdgeOrientation v_GetCartesianEorient(int edge)
            {
                return m_geom->GetCartesianEorient(edge);
            }

        
            virtual void v_AddNormTraceInt(const int dir,
                                           Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray)
            {
                AddNormBoundaryInt(dir,inarray,outarray,true);
            }

            virtual void v_AddNormBoundaryInt(const int dir,
                                              Array<OneD, const NekDouble> &inarray,
                                              Array<OneD,NekDouble> &outarray)
            {
                AddNormBoundaryInt(dir,inarray,outarray,false);
            }
        

            virtual void v_AddBoundaryInt(Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray)
            {
                AddBoundaryInt(inarray,outarray);
            }
        

            virtual void v_AddEdgeBoundaryInt(const int edge, 
                                              SegExpSharedPtr &EdgeExp,
                                              Array <OneD,NekDouble > &outarray)
            {
                AddEdgeBoundaryInt(edge,EdgeExp,outarray);
            }

            virtual void v_AddUDGHelmholtzBoundaryTerms(const NekDouble tau,
                                                        const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD,NekDouble> &outarray)
            {
                AddUDGHelmholtzBoundaryTerms(tau,inarray,outarray);
            }

            virtual void v_AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,NekDouble> &outarray)
            {
                AddUDGHelmholtzTraceTerms(tau,inarray,outarray);
            }

            virtual void v_AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,SegExpSharedPtr> &EdgeExp, 
                                                     Array <OneD,NekDouble > &outarray)
            {
                AddUDGHelmholtzTraceTerms(tau,inarray,EdgeExp,outarray);
            }
            
            virtual void v_GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                GetEdgePhysVals(edge,inarray,outarray);
            }

            virtual void v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
            {
                LaplacianMatrixOp(inarray,outarray);
            }
        
            virtual void v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const double lambda)
            {
                HelmholtzMatrixOp(inarray,outarray,lambda);
            }   
        
        };

        // type defines for use of QuadExp in a boost vector
        typedef boost::shared_ptr<QuadExp> QuadExpSharedPtr;
        typedef std::vector< QuadExpSharedPtr > QuadExpVector;
        typedef std::vector< QuadExpSharedPtr >::iterator QuadExpVectorIter;


    } //end of namespace
} //end of namespace

#endif

/**
 *    $Log: QuadExp.h,v $
 *    Revision 1.33  2008/07/12 17:27:07  sherwin
 *    Update for AddBoundaryInt and moved various members to be private rather than protected
 *
 *    Revision 1.32  2008/07/04 10:19:05  pvos
 *    Some updates
 *
 *    Revision 1.31  2008/07/02 14:09:18  pvos
 *    Implementation of HelmholtzMatOp and LapMatOp on shape level
 *
 *    Revision 1.30  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.29  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.28  2008/05/10 18:27:33  sherwin
 *    Modifications necessary for QuadExp Unified DG Solver
 *
 *    Revision 1.27  2008/04/06 05:59:05  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.26  2008/04/02 22:19:26  pvos
 *    Update for 2D local to global mapping
 *
 *    Revision 1.25  2008/03/12 15:24:29  pvos
 *    Clean up of the code
 *
 *    Revision 1.24  2008/03/06 23:29:23  ehan
 *    Included file <Vmath.hpp> and <VmathArray.hpp>.
 *
 *    Revision 1.23  2008/02/28 10:04:11  sherwin
 *    Modes for UDG codes
 *
 *    Revision 1.22  2008/01/25 16:47:04  sherwin
 *    Added UDG work
 *
 *    Revision 1.21  2007/11/08 16:54:27  pvos
 *    Updates towards 2D helmholtz solver
 *
 *    Revision 1.20  2007/07/28 05:09:33  sherwin
 *    Fixed version with updated MemoryManager
 *
 *    Revision 1.19  2007/07/22 23:04:18  bnelson
 *    Backed out Nektar::ptr.
 *
 *    Revision 1.18  2007/07/20 00:45:51  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.17  2007/07/13 09:02:22  sherwin
 *    Mods for Helmholtz solver
 *
 *    Revision 1.16  2007/07/11 19:26:04  sherwin
 *    update for new Manager structure
 *
 *    Revision 1.15  2007/07/10 17:17:26  sherwin
 *    Introduced Scaled Matrices into the MatrixManager
 *
 *    Revision 1.14  2007/06/07 15:54:19  pvos
 *    Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 *    Also made corrections to various ASSERTL2 calls
 *
 *    Revision 1.13  2007/06/01 17:08:07  pvos
 *    Modification to make LocalRegions/Project2D run correctly (PART1)
 *
 *    Revision 1.12  2007/05/31 19:13:12  pvos
 *    Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 *    Revision 1.11  2007/05/31 11:38:17  pvos
 *    Updated QuadExp and TriExp
 *
 *    Revision 1.10  2007/05/27 16:10:28  bnelson
 *    Update to new Array type.
 *
 *    Revision 1.9  2007/04/26 15:00:16  sherwin
 *    SJS compiling working version using SHaredArrays
 *
 *    Revision 1.8  2007/01/15 21:12:26  sherwin
 *    First definition
 *
 *    Revision 1.7  2006/06/13 18:05:01  sherwin
 *    Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
 *    Revision 1.6  2006/06/05 00:10:01  bnelson
 *    Fixed a gcc 4.1.0 compilation error (ClassName::method not allowed in class definition).
 *
 *    Revision 1.5  2006/06/02 18:48:39  sherwin
 *    Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
 *
 *    Revision 1.4  2006/06/01 14:15:57  sherwin
 *    Added typdef of boost wrappers and made GeoFac a boost shared pointer.
 *
 *    Revision 1.3  2006/05/30 14:00:03  sherwin
 *    Updates to make MultiRegions and its Demos work
 *
 *    Revision 1.2  2006/05/29 17:05:49  sherwin
 *    Modified to put shared_ptr around geom definitions
 *
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.16  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.15  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
