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
#include <LocalRegions/SegExp.h>
//#include <LocalRegions/GenSegExp.h>

#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
    namespace LocalRegions
    {
        
        class TriExp: public StdRegions::StdTriExp, public Expansion2D
        {

        public:

            /** \brief Constructor using BasisKey class for quadrature
                points and order definition */
            TriExp(const LibUtilities::BasisKey &Ba,
                   const LibUtilities::BasisKey &Bb,
                   const SpatialDomains::TriGeomSharedPtr &geom);

            /// Copy Constructor
            TriExp(const TriExp &T);

            /// Destructor
            ~TriExp();

            void GetCoords(Array<OneD,NekDouble> &coords_1,
                           Array<OneD,NekDouble> &coords_2, 
                           Array<OneD,NekDouble> &coords_3 = NullNekDouble1DArray);
            void GetCoord(const Array<OneD, const NekDouble>& Lcoords, 
                          Array<OneD,NekDouble> &coords);

            void GetSurfaceNormal(Array<OneD,NekDouble> &SurfaceNormal,
                                  const int k);

            const SpatialDomains::GeometrySharedPtr GetGeom() const
            {
                return m_geom;
            }

            const SpatialDomains::Geometry2DSharedPtr& GetGeom2D() const
            {
                return m_geom;
            }

            void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v");

            //----------------------------
            // Integration Methods
            //----------------------------

            /// \brief Integrate the physical point list \a inarray over region
            NekDouble Integral(const Array<OneD, const NekDouble> &inarray);

            /** 
                \brief Calculate the inner product of inarray with respect to
                the basis B=base0*base1 and put into outarray:
            
                \f$ 
                \begin{array}{rcl}
                I_{pq} = (\phi_q \phi_q, u) & = & \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1}
                \phi_p(\xi_{0,i}) \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} \xi_{1,j}) 
                J_{i,j}\\
                & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
                \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}  J_{i,j}
                \end{array}
                \f$ 
            
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
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eIProductWRTBase>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::IProductWRTBase_MatOp(inarray,outarray);
                }
                else
                {
                    TriExp::IProductWRTBase_SumFac(inarray,outarray);
                }  
            }  

            void IProductWRTDerivBase(const int dir,
                                      const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray)
            {
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eIProductWRTDerivBase>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::IProductWRTDerivBase_MatOp(dir,inarray,outarray);
                }
                else
                {
                    TriExp::IProductWRTDerivBase_SumFac(dir,inarray,outarray);
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

            void PhysDirectionalDeriv(const Array<OneD, const NekDouble> &inarray, 
                                      const Array<OneD, const NekDouble>& direction,
                                      Array<OneD, NekDouble> &out);  

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
        
            /** \brief Extract the physical values along edge \a edge
                from \a inarray into \a outarray following the local
                edge orientation and point distribution defined by
                defined in \a EdgeExp.
                
                Note: this function will check to see if points
                distribution along the Tri expansion \a edge is the
                same as the local edge definition and if not
                interpolate. If they are the same no interpolation
                will be performed as can be seen in the functino
                LibUtilities::Interp1D
                
                \param edge the edge id which is to be extracted

                \param EdgeExp The Edge Expansion defining the
                orientation and point distrubution points are to be
                interpolated. 

                \param inarray is the 2D physical point set from which
                data is to be extracted.

                \param outarray is the output data
            */
            void GetEdgePhysVals(const int edge, 
                                 const StdRegions::StdExpansion1DSharedPtr &EdgeExp, 
                                 const Array<OneD, const NekDouble> &inarray, 
                                 Array<OneD,NekDouble> &outarray);

            StdRegions::StdExpansion1DSharedPtr GetEdgeExp(int edge, bool SetUpNormals=true);

            void MassMatrixOp(const Array<OneD, const NekDouble> &inarray, 
                              Array<OneD,NekDouble> &outarray,
                              const StdRegions::StdMatrixKey &mkey)
            {              
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eMassMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
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
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eLaplacianMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    TriExp::LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
                }
            }

            void LaplacianMatrixOp(const int k1, const int k2, 
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {           
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eLaplacianMatrixIJOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
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
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eWeakDerivMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
                }
            }

            void WeakDirectionalDerivMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                              Array<OneD,NekDouble> &outarray,
                                              const StdRegions::StdMatrixKey &mkey)
            {
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eWeakDerivMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    StdExpansion::WeakDirectionalDerivMatrixOp_MatFree(inarray,outarray,mkey);
                }
            }
            
            void MassLevelCurvatureMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                              Array<OneD,NekDouble> &outarray,
                                              const StdRegions::StdMatrixKey &mkey)
            {
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eMassMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    StdExpansion::MassLevelCurvatureMatrixOp_MatFree(inarray,outarray,mkey);
                }
            }

            void HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdRegions::StdMatrixKey &mkey)
            {
                bool doMatOp = NekOptimize::ElementalOptimization<StdRegions::eTriExp, NekOptimize::eHelmholtzMatrixOp>::
                    DoMatOp(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                
                if(doMatOp)
                {
                    TriExp::GeneralMatrixOp_MatOp(inarray,outarray,mkey);
                }
                else
                {
                    TriExp::HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
                }
            }       


        protected:
            DNekMatSharedPtr     GenMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekMatSharedPtr     CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);
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


            TriExp();

            virtual StdRegions::ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            }
	     
	    virtual const LibUtilities::BasisSharedPtr& v_GetBasis(int dir) const
            {
	      ASSERTL1(dir >= 0 &&dir <= 1,"input dir is out of range");
	      return m_base[dir];
	    }

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
                return 3;
            }

            virtual int v_NumBndryCoeffs() const
            {
                ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B,
                         "BasisType is not a boundary interior form");

                return 3 + (GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
            } 

            virtual int v_NumDGBndryCoeffs() const
            {
                ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B,
                         "BasisType is not a boundary interior form");

                return GetBasisNumModes(0) + 2*GetBasisNumModes(1);
            } 


            virtual bool v_IsBoundaryInteriorExpansion()
            {
                return StdTriExp::IsBoundaryInteriorExpansion();
            }

            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                return GetEdgeNcoeffs(i);
            }

            virtual int v_GetEdgeNumPoints(const int i) const
            {
                return GetEdgeNumPoints(i);
            }

            virtual void v_GetEdgeToElementMap(const int eid, const StdRegions::EdgeOrientation edgeOrient, 
                                               Array<OneD, unsigned int> &maparray, Array<OneD, int> &signarray)
            {
                StdTriExp::GetEdgeToElementMap(eid,edgeOrient,maparray,signarray);
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

            virtual void v_GetSurfaceNormal(Array<OneD, NekDouble> &SurfaceNormal,
                                            const int k)
            {
                GetSurfaceNormal(SurfaceNormal, k);
            }

            virtual int v_GetCoordim()
            {
                return m_geom->GetCoordim();
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
        
            virtual void v_PhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1,
                                     Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
            }

            virtual void v_PhysDeriv(const int dir, 
                                     const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &outarray)
            {
                PhysDeriv(dir,inarray,outarray);
            }

            virtual void v_PhysDirectionalDeriv(const Array<OneD, const NekDouble>& inarray,
                                                const Array<OneD, const NekDouble>& direction,
                                                Array<OneD, NekDouble> &outarray)
            {
                PhysDirectionalDeriv(inarray,direction,outarray);
            }
        
            /// Virtual call to TriExp::FwdTrans
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

            /// Virtual call to TriExp::BwdTrans
            virtual void v_BwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray,outarray);
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


            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const StdRegions::MatrixType mtype,
							 const Array<OneD, NekDouble> &dir1Forcing,
                                                         const int matrixid = 0,
                                                         NekDouble lambdaval = NekConstants::kNekUnsetDouble, 
                                                         NekDouble tau = NekConstants::kNekUnsetDouble)
            {
	         MatrixKey mkey(mtype,DetExpansionType(),*this,lambdaval,tau,dir1Forcing,matrixid);
		 return m_matrixManager[mkey];
            }

            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const StdRegions::MatrixType mtype,
							 const Array<OneD, Array<OneD, NekDouble> >& dirForcing,
                                                         const int matrixid = 0,
                                                         NekDouble lambdaval = NekConstants::kNekUnsetDouble, 
                                                         NekDouble tau = NekConstants::kNekUnsetDouble)
            {
	        MatrixKey mkey(mtype,DetExpansionType(),*this,lambdaval,tau,dirForcing,matrixid);
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

            virtual void v_AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
						     const Array<OneD, Array<OneD, NekDouble> > &dirForcing, 
                                                     Array <OneD,NekDouble > &outarray,
                                                     const int matrixid)
            {
                Expansion2D::AddHDGHelmholtzTraceTerms(tau,inarray,EdgeExp,dirForcing,outarray,matrixid);
            }
            
            virtual void v_DGDeriv(int dir, 
                                   const Array<OneD, const NekDouble>&incoeffs,
                                   Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                   Array<OneD, NekDouble> &out_d)
            {
              Expansion2D::DGDeriv(dir,incoeffs,EdgeExp,out_d);
            }

            virtual StdRegions::StdExpansion1DSharedPtr v_GetEdgeExp(const int edge, bool SetUpNormals=true)
            {
                return GetEdgeExp(edge,SetUpNormals);
            }

            virtual void v_GetEdgePhysVals(const int edge, const StdRegions::StdExpansion1DSharedPtr &EdgeExp, 
                                           const Array<OneD, const NekDouble> &inarray, 
                                           Array<OneD,NekDouble> &outarray)
            {
                GetEdgePhysVals(edge,EdgeExp,inarray,outarray);
            }

            
            virtual void v_AddEdgeNormBoundaryInt(const int edge, 
                                                  StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                  const Array<OneD, const NekDouble> &Fx,  
                                                  const Array<OneD, const NekDouble> &Fy,  
                                                  Array<OneD, NekDouble> &outarray)
            {
                Expansion2D::AddEdgeNormBoundaryInt(edge,EdgeExp,Fx,Fy,outarray);
                
            }

            virtual void v_AddEdgeNormBoundaryBiInt(const int edge, 
                                                  StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                  const Array<OneD, const NekDouble> &Fwd,  
                                                  const Array<OneD, const NekDouble> &Bwd,  
                                                  Array<OneD, NekDouble> &outarray)
            {
                Expansion2D::AddEdgeNormBoundaryBiInt(edge,EdgeExp,Fwd,Bwd,outarray);
                
            }

            virtual void v_AddEdgeNormBoundaryInt(const int edge, 
                                                  StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                  const Array<OneD, const NekDouble> &Fn,  
                                                  Array<OneD, NekDouble> &outarray)
            {
                Expansion2D::AddEdgeNormBoundaryInt(edge,EdgeExp,Fn,outarray);
                
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

            virtual void v_WeakDirectionalDerivMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD,NekDouble> &outarray,
                                                        const StdRegions::StdMatrixKey &mkey)
            {
                WeakDirectionalDerivMatrixOp(inarray,outarray,mkey);
            }
            
            virtual void v_MassLevelCurvatureMatrixOp(const Array<OneD, const NekDouble> &inarray, 
                                                      Array<OneD,NekDouble> &outarray,
                                                      const StdRegions::StdMatrixKey &mkey)
            {
                MassLevelCurvatureMatrixOp(inarray,outarray,mkey);
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
        
        // type defines for use of TriExp in a boost vector
        typedef boost::shared_ptr<TriExp> TriExpSharedPtr;
        typedef std::vector< TriExpSharedPtr > TriExpVector;
        typedef std::vector< TriExpSharedPtr >::iterator TriExpVectorIter;

    } //end of namespace
} //end of namespace

#endif // TRIEXP_H

/**
 *    $Log: TriExp.h,v $
 *    Revision 1.56  2009/11/10 19:04:24  sehunchun
 *    Variable coefficients for HDG2D Solver
 *
 *    Revision 1.55  2009/11/09 15:43:50  sehunchun
 *    HDG2DManifold Solver with Variable coefficients
 *
 *    Revision 1.54  2009/11/06 21:43:56  sherwin
 *    DGDeriv function
 *
 *    Revision 1.53  2009/11/02 19:15:42  cantwell
 *    Moved ContField1D to inherit from DisContField1D.
 *    Moved ContField3D to inherit from DisContField3D.
 *    Incorporated GenExpList1D functionality into ExpList1D.
 *    Tidied up and added documentation to various classes.
 *    Moved Namespace documentation and introductions to separate files along with
 *    doxygen configuration.
 *    Added option to use system ZLIB library instead of libboost_zlib on UNIX.
 *    Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
 *    Updated Linux compiling instructions.
 *    Updated regDemo to use Helmholtz2D-g when built as debug.
 *
 *    Revision 1.52  2009/10/06 09:43:57  cbiotto
 *    Adding virtual function GetBasis
 *
 *    Revision 1.51  2009/09/06 22:24:00  sherwin
 *    Updates for Navier-Stokes solver
 *
 *    Revision 1.50  2009/08/14 09:20:08  cbiotto
 *    Fixed bug in v_GetEdgeExp
 *
 *    Revision 1.49  2009/07/08 17:19:48  sehunchun
 *    Deleting GetTanBasis
 *
 *    Revision 1.48  2009/07/08 11:11:24  sehunchun
 *    Adding GetSurfaceNormal Function
 *
 *    Revision 1.47  2009/07/07 16:31:47  sehunchun
 *    Adding AddEdgeBoundaryBiInt to line integrate depending on Fwd and Bwd
 *
 *    Revision 1.46  2009/07/03 15:34:52  sehunchun
 *    Adding GetTanBasis function
 *
 *    Revision 1.45  2009/04/27 21:34:07  sherwin
 *    Updated WriteToField
 *
 *    Revision 1.44  2009/04/27 09:38:22  pvos
 *    Fixed some bugs
 *
 *    Revision 1.43  2009/04/20 16:12:28  sherwin
 *    Updates related to output format and optimising DG solver
 *
 *    Revision 1.42  2009/01/21 16:59:57  pvos
 *    Added additional geometric factors to improve efficiency
 *
 *    Revision 1.41  2008/12/18 14:08:24  pvos
 *    NekConstants update
 *
 *    Revision 1.40  2008/11/24 10:31:14  pvos
 *    Changed name from _PartitionedOp to _MatFree
 *
 *    Revision 1.39  2008/11/05 16:08:15  pvos
 *    Added elemental optimisation functionality
 *
 *    Revision 1.38  2008/10/04 19:35:07  sherwin
 *    Added Upwind method which takes normal rather than component fluxes
 *
 *    Revision 1.37  2008/09/23 18:20:25  pvos
 *    Updates for working ProjectContField3D demo
 *
 *    Revision 1.36  2008/08/27 16:35:13  pvos
 *    Small efficiency update
 *
 *    Revision 1.35  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.34  2008/08/14 22:12:57  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *    Revision 1.33  2008/07/31 21:25:13  sherwin
 *    Mods for DG Advection
 *
 *    Revision 1.32  2008/07/31 11:13:22  sherwin
 *    Depracated GetEdgeBasis and replaced with DetEdgeBasisKey
 *
 *    Revision 1.31  2008/07/29 22:25:35  sherwin
 *    general update for DG Advection including separation of GetGeom() into GetGeom1D,2D,3D()
 *
 *    Revision 1.30  2008/07/12 17:27:07  sherwin
 *    Update for AddBoundaryInt and moved various members to be private rather than protected
 *
 *    Revision 1.29  2008/07/04 10:19:05  pvos
 *    Some updates
 *
 *    Revision 1.28  2008/07/02 14:09:18  pvos
 *    Implementation of HelmholtzMatOp and LapMatOp on shape level
 *
 *    Revision 1.27  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.26  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.25  2008/05/10 18:27:33  sherwin
 *    Modifications necessary for QuadExp Unified DG Solver
 *
 *    Revision 1.24  2008/04/06 05:59:05  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.23  2008/04/02 22:19:26  pvos
 *    Update for 2D local to global mapping
 *
 *    Revision 1.22  2008/03/12 15:24:29  pvos
 *    Clean up of the code
 *
 *    Revision 1.21  2008/02/28 10:04:11  sherwin
 *    Modes for UDG codes
 *
 *    Revision 1.20  2007/11/08 16:54:27  pvos
 *    Updates towards 2D helmholtz solver
 *
 *    Revision 1.19  2007/07/28 05:09:33  sherwin
 *    Fixed version with updated MemoryManager
 *
 *    Revision 1.18  2007/07/22 23:04:19  bnelson
 *    Backed out Nektar::ptr.
 *
 *    Revision 1.17  2007/07/20 00:45:52  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.16  2007/07/13 09:02:23  sherwin
 *    Mods for Helmholtz solver
 *
 *    Revision 1.15  2007/07/11 19:26:04  sherwin
 *    update for new Manager structure
 *
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
