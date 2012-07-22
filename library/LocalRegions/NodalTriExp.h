///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriExp.h
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
// Description: Header for NodalTriExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTRIEXP_H
#define NODALTRIEXP_H

#include <LocalRegions/LocalRegions.hpp>
#include <StdRegions/StdNodalTriExp.h>
#include <SpatialDomains/TriGeom.h>

#include <SpatialDomains/GeomFactors2D.h>

#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
    namespace LocalRegions 
    {  
        
    class NodalTriExp: virtual public StdRegions::StdNodalTriExp, virtual public Expansion2D
        {
        public:
            /** \brief Constructor using BasisKey class for quadrature
                points and order definition */
            LOCAL_REGIONS_EXPORT NodalTriExp(const LibUtilities::BasisKey &Ba,
                        const LibUtilities::BasisKey &Bb,
                        const LibUtilities::PointsType Ntype,
                        const SpatialDomains::TriGeomSharedPtr &geom);
            
            /// Copy Constructor
            LOCAL_REGIONS_EXPORT NodalTriExp(const NodalTriExp &T); 
            
            /// Destructor
            LOCAL_REGIONS_EXPORT ~NodalTriExp();
        
            
            LOCAL_REGIONS_EXPORT void GetCoords(Array<OneD,NekDouble> &coords_1,
                           Array<OneD,NekDouble> &coords_2, 
                           Array<OneD,NekDouble> &coords_3 = NullNekDouble1DArray);
            LOCAL_REGIONS_EXPORT void GetCoord(const Array<OneD, const NekDouble>& Lcoords, 
                          Array<OneD,NekDouble> &coords);

            const SpatialDomains::GeometrySharedPtr GetGeom() const
            {
                return m_geom;
            }
            
            const SpatialDomains::Geometry2DSharedPtr& GetGeom2D() const
            {
                return m_geom;
            }

            LOCAL_REGIONS_EXPORT void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v");
            
            //----------------------------
            // Integration Methods
            //----------------------------
            
            /// \brief Integrate the physical point list \a inarray over region
            LOCAL_REGIONS_EXPORT NekDouble Integral(const Array<OneD, const NekDouble> &inarray);
            
            /** \brief  Inner product of \a inarray over region with respect to the
                expansion basis (this)->_Base[0] and return in \a outarray */
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                NodalTriExp::IProductWRTBase_SumFac(inarray,outarray);
            }

            void IProductWRTDerivBase(const int dir,
                                      const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray)
            {
                NodalTriExp::IProductWRTDerivBase_SumFac(dir,inarray,outarray);
            }
            
            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            
            LOCAL_REGIONS_EXPORT void PhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                           Array<OneD, NekDouble> &out_d0,
                           Array<OneD, NekDouble> &out_d1,
                           Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);  
            
            //----------------------------
            // Evaluations Methods
            //---------------------------
            
            /** \brief Forward transform from physical quadrature space
                stored in \a inarray and evaluate the expansion coefficients and
                store in \a (this)->_coeffs  */
            LOCAL_REGIONS_EXPORT void FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                          Array<OneD, NekDouble> &outarray);
            
            LOCAL_REGIONS_EXPORT NekDouble PhysEvaluate(const Array<OneD, const NekDouble> &coord);

            void MassMatrixOp(const Array<OneD, const NekDouble> &inarray, 
                              Array<OneD,NekDouble> &outarray,
                              const StdRegions::StdMatrixKey &mkey)
            {              
                StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
            }

            void LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {           
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }

            void LaplacianMatrixOp(const int k1, const int k2, 
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {           
                StdExpansion::LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,mkey);
            }

            void WeakDerivMatrixOp(const int i,
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {
                StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
            }
            
            void HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {
                StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }  
            
        protected:
            
            DNekMatSharedPtr CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekScalMatSharedPtr    CreateMatrix(const MatrixKey &mkey);
            DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);

            void MultiplyByQuadratureMetric(const Array<OneD, const NekDouble>& inarray,
                                            Array<OneD, NekDouble> &outarray);      

            void IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray, 
                                        Array<OneD, NekDouble> &outarray);
            void IProductWRTBase_MatOp(const Array<OneD, const NekDouble>& inarray, 
                                       Array<OneD, NekDouble> &outarray);

            void IProductWRTDerivBase_SumFac(const int dir,
                                             const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble> & outarray);
            void IProductWRTDerivBase_MatOp(const int dir,
                                             const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble> & outarray);


            void GeneralMatrixOp_MatOp(const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,NekDouble> &outarray,
                                       const StdRegions::StdMatrixKey &mkey); 
            
            virtual DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey);

        private:           
            SpatialDomains::Geometry2DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;
            

            //virtual StdRegions::ExpansionType v_DetExpansionType() const
            //{
            //   return DetExpansionType();
            //}
            
            virtual DNekMatSharedPtr v_GenNBasisTransMatrix()
            {
                return StdNodalTriExp::GenNBasisTransMatrix();
            }
            
            virtual const SpatialDomains::GeomFactorsSharedPtr& v_GetMetricInfo() const
            {
                return m_metricinfo;
            }

            virtual const SpatialDomains::GeometrySharedPtr v_GetGeom() const
            {
                return GetGeom();
            }

            virtual const SpatialDomains::Geometry2DSharedPtr& v_GetGeom2D() const
            {
                return GetGeom2D();
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
            
            virtual void v_GetNodalPoints(Array<OneD, const NekDouble> &x, 
                                          Array<OneD, const NekDouble> &y)
            {
                return StdNodalTriExp::GetNodalPoints(x,y);
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
            
            virtual void v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                                        Array<OneD, NekDouble> &out_d0,
                                        Array<OneD, NekDouble> &out_d1)
            {
                StdTriExp::v_PhysDeriv(inarray, out_d0, out_d1);
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
                Array<OneD,NekDouble> tmp;
                switch(dir)
                {
                case 0:
                    {
                        PhysDeriv(inarray, outarray, tmp);   
                    }
                    break;
                case 1:
                    {
                        PhysDeriv(inarray, tmp, outarray);   
                    }
                    break;
                default:
                    {
                        ASSERTL1(dir >= 0 &&dir < 2,"input dir is out of range");
                    }
                    break;
                }             
            }
            
            /// Virtual call to SegExp::FwdTrans
            virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray,outarray);
            }
            
            /// Virtual call to TriExp::Evaluate
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
            
            virtual DNekScalMatSharedPtr v_GetLocMatrix(const MatrixKey &mkey)
            {
                return m_matrixManager[mkey];
            }
            
//            virtual DNekScalMatSharedPtr v_GetLocMatrix(const StdRegions::MatrixType mtype, NekDouble lambdaval, NekDouble tau)
//            {
//                MatrixKey mkey(mtype,DetExpansionType(),*this,lambdaval,tau);
//                return m_matrixManager[mkey];
//            }
            
            virtual DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(const MatrixKey &mkey)
            {
                return m_staticCondMatrixManager[mkey];
            }

            virtual void v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                StdNodalTriExp::v_BwdTrans_SumFac(inarray,outarray);
            }    
            
            virtual void v_IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray, 
                                                  Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase_SumFac(inarray,outarray);
            }            
            
            virtual void v_IProductWRTDerivBase_SumFac(const int dir,
                                                       const Array<OneD, const NekDouble>& inarray, 
                                                       Array<OneD, NekDouble> &outarray)
            {
                IProductWRTDerivBase_SumFac(dir,inarray,outarray);
            }    

            virtual void v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray, 
                                        Array<OneD,NekDouble> &outarray,
                                        const StdRegions::StdMatrixKey &mkey)
            {
                MassMatrixOp(inarray,outarray,mkey);
            }  
            
            virtual void v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdRegions::StdMatrixKey &mkey)
            {
                LaplacianMatrixOp(inarray,outarray,mkey);
            }

            virtual void v_LaplacianMatrixOp(const int k1, const int k2, 
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdRegions::StdMatrixKey &mkey)
            {
                LaplacianMatrixOp(k1,k2,inarray,outarray,mkey);
            }

            virtual void v_WeakDerivMatrixOp(const int i,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdRegions::StdMatrixKey &mkey)
            {
                WeakDerivMatrixOp(i,inarray,outarray,mkey);
            }
            
            virtual void v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdRegions::StdMatrixKey &mkey)
            {
                HelmholtzMatrixOp(inarray,outarray,mkey);
            }  
            
            void v_ComputeEdgeNormal(const int edge);
        };
    
        // type defines for use of TriExp in a boost vector
        typedef boost::shared_ptr<NodalTriExp> NodalTriExpSharedPtr;
        typedef std::vector< NodalTriExpSharedPtr > NodalTriExpVector;
        typedef std::vector< NodalTriExpSharedPtr >::iterator NodalTriExpVectorIter;
    
    } //end of namespace
} //end of namespace

#endif // NODALTRIEXP_H

/** 
 *    $Log: NodalTriExp.h,v $
 *    Revision 1.27  2009/04/29 10:53:48  pvos
 *    Made demos working with nodal expansions
 *
 *    Revision 1.26  2009/04/27 21:34:07  sherwin
 *    Updated WriteToField
 *
 *    Revision 1.25  2009/01/21 16:59:57  pvos
 *    Added additional geometric factors to improve efficiency
 *
 *    Revision 1.24  2008/11/24 10:31:14  pvos
 *    Changed name from _PartitionedOp to _MatFree
 *
 *    Revision 1.23  2008/11/19 16:01:41  pvos
 *    Added functionality for variable Laplacian coeffcients
 *
 *    Revision 1.22  2008/11/05 16:08:15  pvos
 *    Added elemental optimisation functionality
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
 *    Revision 1.18  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.17  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.16  2008/04/06 05:59:05  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.15  2008/04/02 22:19:26  pvos
 *    Update for 2D local to global mapping
 *
 *    Revision 1.14  2008/03/18 14:12:53  pvos
 *    Update for nodal triangular helmholtz solver
 *
 *    Revision 1.13  2007/07/28 05:09:32  sherwin
 *    Fixed version with updated MemoryManager
 *
 *    Revision 1.12  2007/07/22 23:04:17  bnelson
 *    Backed out Nektar::ptr.
 *
 *    Revision 1.11  2007/07/20 00:45:50  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.10  2007/07/13 09:02:21  sherwin
 *    Mods for Helmholtz solver
 *
 *    Revision 1.9  2007/07/11 19:26:03  sherwin
 *    update for new Manager structure
 *
 *    Revision 1.8  2007/07/10 17:17:24  sherwin
 *    Introduced Scaled Matrices into the MatrixManager
 *
 *    Revision 1.7  2007/06/07 15:54:18  pvos
 *    Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 *    Also made corrections to various ASSERTL2 calls
 *
 *    Revision 1.6  2007/06/01 17:08:07  pvos
 *    Modification to make LocalRegions/Project2D run correctly (PART1)
 *
 *    Revision 1.5  2007/05/31 19:13:12  pvos
 *    Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 *    Revision 1.4  2006/12/10 18:59:46  sherwin
 *    Updates for Nodal points
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
 *    Revision 1.10  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.9  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
