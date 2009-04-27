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

#include <SpatialDomains/GeomFactors.h>

#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/GenSegExp.h>

#include <LibUtilities/BasicUtils/VmathArray.hpp> 
#include <LibUtilities/BasicUtils/Vmath.hpp> 

#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
    namespace LocalRegions
    {
        
        class QuadExp: public StdRegions::StdQuadExp, public Expansion2D
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

            const SpatialDomains::Geometry2DSharedPtr& GetGeom2D() const
            {
                return m_geom;
            }

            const SpatialDomains::GeometrySharedPtr GetGeom() const
            {
                return m_geom;
            }

            void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true);

            //----------------------------
            // Integration Methods
            //----------------------------

            /// \brief Integrate the physical point list \a inarray over region
            NekDouble Integral(const Array<OneD, const NekDouble> &inarray);

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
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                if(m_base[0]->Collocation() && m_base[1]->Collocation())
                {  
                    MultiplyByQuadratureMetric(inarray,outarray);
                }
                else
                {
                    bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eQuadExp, NekOptimize::eIProductWRTBase>::
                        DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                    
                    if(doMatOp)
                    {
                        QuadExp::IProductWRTBase_MatOp(inarray,outarray);
                    }
                    else
                    {
                        QuadExp::IProductWRTBase_SumFac(inarray,outarray);
                    }  
                }
            }

            void IProductWRTDerivBase(const int dir,
                                      const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray)
            {
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eQuadExp, NekOptimize::eIProductWRTDerivBase>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    QuadExp::IProductWRTDerivBase_MatOp(dir,inarray,outarray);
                }
                else
                {
                    QuadExp::IProductWRTDerivBase_SumFac(dir,inarray,outarray);
                }  
            }
        
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

            void FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray, 
                                         Array<OneD, NekDouble> &outarray);
        
            NekDouble PhysEvaluate(const Array<OneD, const NekDouble> &coord);        
        
            StdRegions::StdExpansion1DSharedPtr GetEdgeExp(int edge, bool SetUpNormal=true);

            void GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray);

            void GetEdgePhysVals(const int edge, const StdRegions::StdExpansion1DSharedPtr &EdgeExp, 
                                 const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray);
               
            void MassMatrixOp(const Array<OneD, const NekDouble> &inarray, 
                              Array<OneD,NekDouble> &outarray,
                              const StdRegions::StdMatrixKey &mkey)
            {              
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eQuadExp, NekOptimize::eMassMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    QuadExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
                }
            }

            void LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {           
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eQuadExp, NekOptimize::eLaplacianMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    QuadExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    QuadExp::LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
                }
            }

            void LaplacianMatrixOp(const int k1, const int k2, 
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {           
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eQuadExp, NekOptimize::eLaplacianMatrixIJOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    QuadExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    StdExpansion::LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,mkey);
                }
            }

            void WeakDerivMatrixOp(const int i,
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eQuadExp, NekOptimize::eWeakDerivMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    QuadExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
                }
            }
            
            void HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eQuadExp, NekOptimize::eHelmholtzMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    QuadExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    QuadExp::HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
                }
            }       

        protected:

            void GenMetricInfo();

            DNekMatSharedPtr GenMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekMatSharedPtr CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekScalMatSharedPtr  CreateMatrix(const MatrixKey &mkey);
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
            
            void LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,NekDouble> &outarray,
                                                 const StdRegions::StdMatrixKey &mkey);
            void HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD,NekDouble> &outarray,
                                                 const StdRegions::StdMatrixKey &mkey);   

        private:
            SpatialDomains::Geometry2DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            QuadExp();

            virtual int v_GetNumPoints(const int dir) const 
            {
                return GetNumPoints(dir);
            }

            virtual int v_GetNcoeffs() const
            {
                return m_ncoeffs;
            }
                
            virtual int v_GetNedges() const
            {
                return 4;
            }

            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                return GetEdgeNcoeffs(i);
            }

            virtual int v_GetEdgeNumPoints(const int i) const
            {
                return GetEdgeNumPoints(i);
            }

            virtual bool v_IsBoundaryInteriorExpansion()
            {
                return StdQuadExp::IsBoundaryInteriorExpansion();
            }

            virtual int v_NumBndryCoeffs() const
            {
                ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                         GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                         GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");

                return 4 + 2*(GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
            } 

            virtual int v_NumDGBndryCoeffs() const
            {
                ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                         GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                         GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");

                return  2*GetBasisNumModes(0) + 2*GetBasisNumModes(1);
            } 
            virtual StdRegions::ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            }
            
            virtual void v_GetEdgeToElementMap(const int eid, const StdRegions::EdgeOrientation edgeOrient, Array<OneD, unsigned int> &maparray, Array<OneD, int> &signarray)
            {
                StdQuadExp::GetEdgeToElementMap(eid,edgeOrient,maparray,signarray);
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
    
            /// Virtual call to QuadExp::FwdTrans
            virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray,outarray);
            }

            virtual void v_FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray, 
                                                   Array<OneD, NekDouble> &outarray)
            {
                FwdTrans_BndConstrained(inarray, outarray); 
            }

            /// Virtual call to StdQuadExp::BwdTrans
            virtual void v_BwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray,outarray);
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
        
            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const StdRegions::MatrixType mtype, 
                                                         NekDouble lambdaval = NekConstants::kNekUnsetDouble, 
                                                         NekDouble tau = NekConstants::kNekUnsetDouble)
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

            virtual void v_AddEdgeNormBoundaryInt(const int edge, 
                                                  StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                  const Array<OneD, const NekDouble> &Fx,  
                                                  const Array<OneD, const NekDouble> &Fy,  
                                                  Array<OneD, NekDouble> &outarray)
            {
                Expansion2D::AddEdgeNormBoundaryInt(edge,EdgeExp,Fx,Fy,outarray);
                
            }


            virtual void v_AddEdgeNormBoundaryInt(const int edge, 
                                                  StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                  const Array<OneD, const NekDouble> &Fn,  
                                                  Array<OneD, NekDouble> &outarray)
            {
                Expansion2D::AddEdgeNormBoundaryInt(edge,EdgeExp,Fn,outarray);
                
            }


            virtual void v_AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp, 
                                                     Array <OneD,NekDouble > &outarray)
            {
                Expansion2D::AddHDGHelmholtzTraceTerms(tau,inarray,EdgeExp,outarray);
            }
            


            virtual StdRegions::StdExpansion1DSharedPtr v_GetEdgeExp(const int edge, bool SetUpNormals=true)
            {
                return GetEdgeExp(edge,true);
            }

            virtual void v_GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                GetEdgePhysVals(edge,inarray,outarray);
            }
            
            virtual void v_GetEdgePhysVals(const int edge, 
                                           const StdRegions::StdExpansion1DSharedPtr &EdgeExp, 
                                           const Array<OneD, const NekDouble> &inarray, 
                                           Array<OneD,NekDouble> &outarray)
            {
                GetEdgePhysVals(edge,EdgeExp,inarray,outarray);
            }
            
            virtual void v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                BwdTrans_SumFac(inarray,outarray);
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
            
            virtual void v_LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                           Array<OneD,NekDouble> &outarray,
                                                           const StdRegions::StdMatrixKey &mkey)
            {
                LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
            }
            
            virtual void v_HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                           Array<OneD,NekDouble> &outarray,
                                                           const StdRegions::StdMatrixKey &mkey)
            {
                HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
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
 *    Revision 1.48  2009/04/20 16:12:28  sherwin
 *    Updates related to output format and optimising DG solver
 *
 *    Revision 1.47  2009/01/21 16:59:57  pvos
 *    Added additional geometric factors to improve efficiency
 *
 *    Revision 1.46  2008/12/18 14:08:24  pvos
 *    NekConstants update
 *
 *    Revision 1.45  2008/11/24 10:31:14  pvos
 *    Changed name from _PartitionedOp to _MatFree
 *
 *    Revision 1.44  2008/11/05 16:08:15  pvos
 *    Added elemental optimisation functionality
 *
 *    Revision 1.43  2008/10/04 19:35:07  sherwin
 *    Added Upwind method which takes normal rather than component fluxes
 *
 *    Revision 1.42  2008/09/23 18:20:25  pvos
 *    Updates for working ProjectContField3D demo
 *
 *    Revision 1.41  2008/08/27 16:35:13  pvos
 *    Small efficiency update
 *
 *    Revision 1.40  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.39  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *    Revision 1.38  2008/07/31 11:13:22  sherwin
 *    Depracated GetEdgeBasis and replaced with DetEdgeBasisKey
 *
 *    Revision 1.37  2008/07/29 22:25:34  sherwin
 *    general update for DG Advection including separation of GetGeom() into GetGeom1D,2D,3D()
 *
 *    Revision 1.36  2008/07/19 21:15:38  sherwin
 *    Removed MapTo function, made orientation anticlockwise, changed enum from BndSys to BndLam
 *
 *    Revision 1.35  2008/07/16 22:20:54  sherwin
 *    Added AddEdgeNormBoundaryInt
 *
 *    Revision 1.34  2008/07/12 19:08:29  sherwin
 *    Modifications for DG advection routines
 *
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
