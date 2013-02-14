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

#include <StdRegions/StdNodalTriExp.h>
#include <SpatialDomains/TriGeom.h>

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
            
            virtual void v_StdPhysDeriv(
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,       NekDouble> &out_d0,
                      Array<OneD,       NekDouble> &out_d1,
                      Array<OneD,       NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                StdTriExp::v_PhysDeriv(inarray, out_d0, out_d1, out_d2);
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

