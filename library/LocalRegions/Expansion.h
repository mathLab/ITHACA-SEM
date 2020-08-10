///////////////////////////////////////////////////////////////////////////////
//
// File Expansion.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Header file for Expansion routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION_H
#define EXPANSION_H

#include <StdRegions/StdExpansion.h>
#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/GeomFactors.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/IndexMapKey.h>
#include <memory>
#include <vector>
#include <map>

namespace Nektar
{
    namespace LocalRegions 
    {

        class Expansion;
        class MatrixKey;
        
        typedef Array<OneD, Array<OneD, NekDouble> > NormalVector;
        
        enum MetricType
        {
            eMetricLaplacian00,
            eMetricLaplacian01,
            eMetricLaplacian02,
            eMetricLaplacian11,
            eMetricLaplacian12,
            eMetricLaplacian22,
            eMetricQuadrature
        };

        typedef std::shared_ptr<Expansion> ExpansionSharedPtr;
        typedef std::weak_ptr<Expansion> ExpansionWeakPtr;
        typedef std::vector< ExpansionSharedPtr > ExpansionVector;
        typedef std::map<MetricType, Array<OneD, NekDouble> > MetricMap;
        
        class Expansion : virtual public StdRegions::StdExpansion
        {
        public:
            LOCAL_REGIONS_EXPORT Expansion(SpatialDomains::GeometrySharedPtr pGeom); // default constructor.
            LOCAL_REGIONS_EXPORT Expansion(const Expansion &pSrc); // copy constructor.
            LOCAL_REGIONS_EXPORT virtual ~Expansion();
            
            LOCAL_REGIONS_EXPORT void SetTraceExp(const int traceid, ExpansionSharedPtr &f);                
            LOCAL_REGIONS_EXPORT ExpansionSharedPtr GetTraceExp(const int traceid);            
            
            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr
                          GetLocMatrix(const LocalRegions::MatrixKey &mkey);
            
            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr GetLocMatrix
                (const StdRegions::MatrixType mtype,
                 const StdRegions::ConstFactorMap &factors = StdRegions::NullConstFactorMap,
                 const StdRegions::VarCoeffMap &varcoeffs = StdRegions::NullVarCoeffMap);

            LOCAL_REGIONS_EXPORT SpatialDomains::GeometrySharedPtr GetGeom() const;
            
            LOCAL_REGIONS_EXPORT void Reset();
            
            LOCAL_REGIONS_EXPORT IndexMapValuesSharedPtr CreateIndexMap(const IndexMapKey &ikey);

            LOCAL_REGIONS_EXPORT const SpatialDomains::GeomFactorsSharedPtr& GetMetricInfo() const;

            LOCAL_REGIONS_EXPORT DNekMatSharedPtr BuildTransformationMatrix(
                          const DNekScalMatSharedPtr &r_bnd, 
                          const StdRegions::MatrixType matrixType);
            
            LOCAL_REGIONS_EXPORT DNekMatSharedPtr BuildVertexMatrix(
                          const DNekScalMatSharedPtr &r_bnd);

            LOCAL_REGIONS_EXPORT void ExtractDataToCoeffs(
                          const NekDouble *data,
                          const std::vector<unsigned int > &nummodes,
                          const int nmodes_offset,
                          NekDouble *coeffs,
                          std::vector<LibUtilities::BasisType> &fromType);

            LOCAL_REGIONS_EXPORT void AddEdgeNormBoundaryInt(
                    const int                           edge,
                    const std::shared_ptr<Expansion>   &EdgeExp,
                    const Array<OneD, const NekDouble> &Fx,
                    const Array<OneD, const NekDouble> &Fy,
                          Array<OneD,       NekDouble> &outarray);
            LOCAL_REGIONS_EXPORT void AddEdgeNormBoundaryInt(
                    const int                           edge,
                    const std::shared_ptr<Expansion>   &EdgeExp,
                    const Array<OneD, const NekDouble> &Fn,
                          Array<OneD,       NekDouble> &outarray);
            LOCAL_REGIONS_EXPORT void AddFaceNormBoundaryInt(
                    const int                           face,
                    const std::shared_ptr<Expansion>   &FaceExp,
                    const Array<OneD, const NekDouble> &Fn,
                          Array<OneD,       NekDouble> &outarray);
            LOCAL_REGIONS_EXPORT void DGDeriv(
                    const int                                   dir,
                    const Array<OneD, const NekDouble>&         inarray,
                          Array<OneD, ExpansionSharedPtr>      &EdgeExp,
                          Array<OneD, Array<OneD, NekDouble> > &coeffs,
                          Array<OneD,             NekDouble>   &outarray);
            LOCAL_REGIONS_EXPORT NekDouble VectorFlux(
                    const Array<OneD, Array<OneD, NekDouble > > &vec);

            inline IndexMapValuesSharedPtr GetIndexMap(const IndexMapKey &ikey)
            {
                return m_indexMapManager[ikey];
            }

            inline ExpansionSharedPtr GetLeftAdjacentElementExp() const;

            inline ExpansionSharedPtr GetRightAdjacentElementExp() const;
            
            inline int GetLeftAdjacentElementTrace() const;
            
            inline int GetRightAdjacentElementTrace() const;

            inline void SetAdjacentElementExp(
                    int                traceid,
                    ExpansionSharedPtr &e);
                
            inline StdRegions::Orientation GetTraceOrient(int trace)
            {
                return v_GetTraceOrient(trace);
            }
            
            inline void SetCoeffsToOrientation(
                StdRegions::Orientation dir,
                Array<OneD, const NekDouble> &inarray,
                Array<OneD, NekDouble> &outarray)
            {
                v_SetCoeffsToOrientation(dir,inarray,outarray);
            }

            /// Divided by the metric jacobi and quadrature weights
            inline void DivideByQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray)
            {
                v_DivideByQuadratureMetric(inarray, outarray);
            }
            
            /**
             * @brief Extract the metric factors to compute the contravariant
             * fluxes along edge \a edge and stores them into \a outarray
             * following the local edge orientation (i.e. anticlockwise
             * convention).
             */
            inline void GetTraceQFactors(const int trace,
                                     Array<OneD, NekDouble> &outarray)
            {
                v_GetTraceQFactors(trace, outarray);
            }

            inline void GetTracePhysVals(const int trace,
                  const StdRegions::StdExpansionSharedPtr &TraceExp,
                  const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
                  StdRegions::Orientation orient = StdRegions::eNoOrientation)
            {
                v_GetTracePhysVals(trace,TraceExp,inarray,outarray,orient);
            }
            
            inline void GetTracePhysMap(
                    const int           edge,
                    Array<OneD, int>   &outarray)
            {
                v_GetTracePhysMap(edge, outarray);
            }
        
            inline void ReOrientTracePhysMap(
                            const StdRegions::Orientation orient,
                            Array<OneD, int> &idmap,
                            const int nq0,
                            const int nq1)
            {
                v_ReOrientTracePhysMap(orient,idmap,nq0,nq1);
            }
            
            inline const NormalVector &GetTraceNormal(const int id) const
            {
                return v_GetTraceNormal(id);
            }
            
            inline void ComputeTraceNormal(const int id)
            {
                v_ComputeTraceNormal(id);
            }
        
            inline const Array<OneD, const NekDouble>& GetPhysNormals(void)
            {
                return v_GetPhysNormals();
            }
            
            inline void SetPhysNormals(Array<OneD, const NekDouble> &normal)
            {
                v_SetPhysNormals(normal);
            }
            
            inline void SetUpPhysNormals(const int trace)
            {
                v_SetUpPhysNormals(trace);
            }
            
            inline void AddRobinMassMatrix(const int traceid,
                          const Array<OneD, const NekDouble > &primCoeffs,
                                           DNekMatSharedPtr &inoutmat)
            {
                v_AddRobinMassMatrix(traceid,primCoeffs,inoutmat);
            }
            

            virtual void AddRobinTraceContribution(
                const int traceid,
                const Array<OneD, const NekDouble> &primCoeffs,
                const Array<OneD, NekDouble> &incoeffs,
                Array<OneD, NekDouble> &coeffs)
            {
                v_AddRobinTraceContribution(traceid,primCoeffs,incoeffs,coeffs);
            }

            LOCAL_REGIONS_EXPORT const Array<OneD, const NekDouble > 
                        &GetElmtBndNormDirElmtLen(const int nbnd) const;
            
        protected:
	    LibUtilities::NekManager<IndexMapKey,
                      IndexMapValues, IndexMapKey::opLess> m_indexMapManager;

            std::vector<ExpansionWeakPtr>        m_traceExp;
            SpatialDomains::GeometrySharedPtr    m_geom;
            SpatialDomains::GeomFactorsSharedPtr m_metricinfo;
            MetricMap        m_metrics;
            ExpansionWeakPtr m_elementLeft;
            ExpansionWeakPtr m_elementRight;
            int              m_elementTraceLeft  = -1;
            int              m_elementTraceRight = -1;

            /// the element length in each element boundary(Vertex, edge
            /// or face) normal direction calculated based on the local
            /// m_metricinfo times the standard element length (which is
            /// 2.0)
            std::map<int, Array<OneD, NekDouble>> m_elmtBndNormDirElmtLen;
            
            void ComputeLaplacianMetric();
            void ComputeQuadratureMetric();
            void ComputeGmatcdotMF(const Array<TwoD, const NekDouble> &df,
                                   const Array<OneD, const NekDouble> &direction,
                                   Array<OneD, Array<OneD, NekDouble> > &dfdir);
            
            virtual void v_MultiplyByQuadratureMetric(
                     const Array<OneD, const NekDouble> &inarray,
                     Array<OneD,       NekDouble> &outarray);

            virtual void v_DivideByQuadratureMetric(
                     const Array<OneD, 
                     const NekDouble>& inarray,
                     Array<OneD, NekDouble> &outarray);

            virtual void v_ComputeLaplacianMetric()
                         {}
            
            virtual void v_GetCoords(Array<OneD,NekDouble> &coords_1,
                                     Array<OneD,NekDouble> &coords_2,
                                     Array<OneD,NekDouble> &coords_3);
            
            Array<OneD, NekDouble> v_GetMF(
                    const int dir,
                    const int shapedim,
                    const StdRegions::VarCoeffMap &varcoeffs);

            Array<OneD, NekDouble> v_GetMFDiv(
                    const int dir,
                    const StdRegions::VarCoeffMap &varcoeffs);

            Array<OneD, NekDouble> v_GetMFMag(
                    const int dir,
                    const StdRegions::VarCoeffMap &varcoeffs);

            virtual DNekScalMatSharedPtr v_GetLocMatrix(
                    const LocalRegions::MatrixKey &mkey);

            virtual DNekMatSharedPtr v_BuildTransformationMatrix(
                    const DNekScalMatSharedPtr &r_bnd,
                    const StdRegions::MatrixType matrixType);

            virtual DNekMatSharedPtr v_BuildVertexMatrix(
                    const DNekScalMatSharedPtr &r_bnd);

            virtual void v_ExtractDataToCoeffs(
                    const NekDouble *data,
                    const std::vector<unsigned int > &nummodes,
                    const int nmodes_offset,
                    NekDouble *coeffs,
                    std::vector<LibUtilities::BasisType> &fromType);

            virtual void v_AddEdgeNormBoundaryInt(
                    const int                           edge,
                    const std::shared_ptr<Expansion>   &EdgeExp,
                    const Array<OneD, const NekDouble> &Fx,
                    const Array<OneD, const NekDouble> &Fy,
                          Array<OneD,       NekDouble> &outarray);
            virtual void v_AddEdgeNormBoundaryInt(
                    const int                           edge,
                    const std::shared_ptr<Expansion>   &EdgeExp,
                    const Array<OneD, const NekDouble> &Fn,
                          Array<OneD,       NekDouble> &outarray);
            virtual void v_AddFaceNormBoundaryInt(
                    const int                           face,
                    const std::shared_ptr<Expansion>   &FaceExp,
                    const Array<OneD, const NekDouble> &Fn,
                          Array<OneD,       NekDouble> &outarray);
            virtual void v_DGDeriv(
                    const int                                   dir,
                    const Array<OneD, const NekDouble>&         inarray,
                          Array<OneD, ExpansionSharedPtr>      &EdgeExp,
                          Array<OneD, Array<OneD, NekDouble> > &coeffs,
                          Array<OneD,             NekDouble>   &outarray);
            virtual NekDouble v_VectorFlux(
                    const Array<OneD, Array<OneD, NekDouble > > &vec);


            virtual StdRegions::Orientation v_GetTraceOrient(int trace);

            virtual void v_SetCoeffsToOrientation(StdRegions::Orientation dir,
                                                  Array<OneD, const NekDouble> &inarray,
                                                  Array<OneD, NekDouble> &outarray);

            virtual void v_GetTraceQFactors(const int trace,
                                           Array<OneD, NekDouble> &outarray);

            virtual void v_GetTracePhysVals(const int trace,
                     const StdRegions::StdExpansionSharedPtr &TraceExp,
                     const Array<OneD, const NekDouble> &inarray,
                     Array<OneD,       NekDouble> &outarray,
                     StdRegions::Orientation  orient);
            
            virtual void v_GetTracePhysMap( const int       edge,
                                            Array<OneD,int> &outarray);

            virtual void v_ReOrientTracePhysMap(
                                const StdRegions::Orientation orient,
                                 Array<OneD, int> &idmap,
                                 const int nq0,
                                 const int nq1 = -1);
            
            virtual const NormalVector & v_GetTraceNormal(const int id) const;

            virtual void v_ComputeTraceNormal(const int id);

            virtual const Array<OneD, const NekDouble>& v_GetPhysNormals(void);

            virtual void v_SetPhysNormals(Array<OneD, const NekDouble> &normal);
            
            virtual void v_SetUpPhysNormals(const int id);
                
            virtual void v_AddRobinMassMatrix(
                const int                           face, 
                const Array<OneD, const NekDouble> &primCoeffs, 
                DNekMatSharedPtr                   &inoutmat);

            virtual void v_AddRobinTraceContribution(
                const int traceid,
                const Array<OneD, const NekDouble> &primCoeffs,
                const Array<OneD, NekDouble> &incoeffs,
                Array<OneD, NekDouble> &coeffs);

        private:

        };

        inline ExpansionSharedPtr Expansion::GetTraceExp(int  traceid)
        {
            ASSERTL1(traceid < GetNtraces(), "Trace is out of range.");
            return m_traceExp[traceid].lock();
        }

        inline void Expansion::SetTraceExp(
            const int           traceid,
            ExpansionSharedPtr &exp)
        {
            int nTraces = GetNtraces();
            ASSERTL1(traceid < nTraces, "Edge out of range.");
            if (m_traceExp.size() < nTraces)
            {
                m_traceExp.resize(nTraces);
            }

            m_traceExp[traceid] = exp;
        }
        
        inline ExpansionSharedPtr Expansion::
            GetLeftAdjacentElementExp() const
        {
            ASSERTL1(m_elementLeft.lock().get(),
                     "Left adjacent element not set.");
            return m_elementLeft.lock();
        }

        inline ExpansionSharedPtr Expansion::
            GetRightAdjacentElementExp() const
        {
            ASSERTL1(m_elementLeft.lock().get(),
                     "Right adjacent element not set.");
            
            return m_elementRight.lock();
        }

        inline int Expansion::GetLeftAdjacentElementTrace() const
        {
            return m_elementTraceLeft;
        }

        inline int Expansion::GetRightAdjacentElementTrace() const
        {
            return m_elementTraceRight;
        }

        inline void Expansion::SetAdjacentElementExp(
            int                 traceid,
            ExpansionSharedPtr &exp)
        {
            if (m_elementLeft.lock().get())
            {
                m_elementRight      = exp;
                m_elementTraceRight = traceid;
            }
            else
            {
                m_elementLeft      = exp;
                m_elementTraceLeft = traceid;
            }
        }

        
        
    } //end of namespace
} //end of namespace

#endif
