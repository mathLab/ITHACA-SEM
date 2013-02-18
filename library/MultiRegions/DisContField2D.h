///////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.h
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
// Description: Field definition in two-dimensions for a discontinuous LDG-H
// expansion.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SpatialDomains/Conditions.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField2D : public ExpList2D
        {
        public:
            MULTI_REGIONS_EXPORT DisContField2D();

            MULTI_REGIONS_EXPORT DisContField2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr   &graph2D,
                const std::string                          &variable,
                const bool SetUpJustDG            = true,
                const bool DeclareCoeffPhysArrays = true);
            
            MULTI_REGIONS_EXPORT DisContField2D(
                const DisContField2D                     &In,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const std::string                        &variable,
                const bool SetUpJustDG            = false,
                const bool DeclareCoeffPhysArrays = true);
            
            MULTI_REGIONS_EXPORT DisContField2D(
                const DisContField2D &In,
                const bool DeclareCoeffPhysArrays = true);

            MULTI_REGIONS_EXPORT virtual ~DisContField2D();
            
            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr GetGlobalBndLinSys(
                const GlobalLinSysKey &mkey);

            MULTI_REGIONS_EXPORT NekDouble L2_DGDeriv(
                const int                           dir,
                const Array<OneD, const NekDouble> &soln);

            MULTI_REGIONS_EXPORT void EvaluateHDGPostProcessing(
                Array<OneD, NekDouble> &outarray);

        protected:
            /**
             * @brief An object which contains the discretised boundary
             * conditions.
             *
             * It is an array of size equal to the number of boundary regions
             * and consists of entries of the type MultiRegions#ExpList1D. Every
             * entry corresponds to the one-dimensional spectral/hp expansion on
             * a single boundary region.  The values of the boundary conditions
             * are stored as the coefficients of the one-dimensional expansion.
             */
            Array<OneD,MultiRegions::ExpListSharedPtr> m_bndCondExpansions;
            
            /**
             * @brief An array which contains the information about
             * the boundary condition on the different boundary regions.
             */
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            GlobalLinSysMapShPtr        m_globalBndMat;
            ExpListSharedPtr            m_trace;
            AssemblyMapDGSharedPtr      m_traceMap;
            
            Array<OneD, Array<OneD, unsigned int> >     m_mapEdgeToElmn;
            Array<OneD, Array<OneD, unsigned int> >     m_signEdgeToElmn;
            Array<OneD,StdRegions::Orientation>    m_edgedir;

            /**
             * @brief A set storing the global IDs of any boundary edges.
             */
            std::set<int> m_boundaryEdges;
            
            /**
             * @brief A map which identifies pairs of periodic edges.
             */
            map<int,int> m_periodicEdges;
	    
            /**
             * @brief A map identifying pairs of periodic vertices.
             */
            vector<map<int,int> > m_periodicVertices;
            
            /**
             * @brief Auxiliary map for periodic boundary conditions.
             * 
             * Takes geometry IDs of periodic edges to a pair (n,e), where n
             * is the expansion containing the edge and e the local edge number.
             */
            boost::unordered_map<int,pair<int,int> > m_perEdgeToExpMap;

            /*
             * @brief A map identifying which edges are left- and right-adjacent
             * for DG.
             */
            vector<bool> m_leftAdjacentEdges;

            void SetUpDG();
            bool SameTypeOfBoundaryConditions(const DisContField2D &In);
            void GenerateBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string                        &variable,
                const bool DeclareCoeffPhysArrays = true);
            void FindPeriodicEdges(
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string                        &variable);
            
            bool IsLeftAdjacentEdge(const int n, const int e);

            virtual void v_GetFwdBwdTracePhys(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD,       NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_GetFwdBwdTracePhys(
                      Array<OneD,       NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fx,
                const Array<OneD, const NekDouble> &Fy,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fn,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_AddFwdBwdTraceIntegral(
                const Array<OneD, const NekDouble> &Fwd, 
                const Array<OneD, const NekDouble> &Bwd, 
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_ExtractTracePhys(
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_ExtractTracePhys(
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const FlagList &flags,
                const StdRegions::ConstFactorMap   &factors,
                const StdRegions::VarCoeffMap      &varcoeff,
                const Array<OneD, const NekDouble> &dirForcing);
            virtual void v_GeneralMatrixOp(
                const GlobalMatrixKey             &gkey,
                const Array<OneD,const NekDouble> &inarray,
                      Array<OneD,      NekDouble> &outarray,
                CoeffState                         coeffstate = eLocal);
            virtual void v_GetBoundaryToElmtMap(
                Array<OneD, int> &ElmtID,
                Array<OneD, int> &EdgeID);
            
            /**
             * @brief Obtain a copy of the periodic edges and vertices for this
             * field.
             */
            virtual void v_GetPeriodicEdges(
                vector<map<int,int> > &periodicVertices,
                map<int,int>          &periodicEdges)
            {
                periodicVertices = m_periodicVertices;
                periodicEdges    = m_periodicEdges;
            }
            
            virtual ExpListSharedPtr &v_GetTrace()
            {
                if(m_trace == NullExpListSharedPtr)
                {
                    SetUpDG();
                }

                return m_trace;
            }

            virtual AssemblyMapDGSharedPtr &v_GetTraceMap()
            {
                return m_traceMap;
            }

            virtual const Array<OneD,const MultiRegions::ExpListSharedPtr>
                &v_GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }

            virtual const 
                Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                &v_GetBndConditions()
            {
                return m_bndConditions;
            }

            virtual MultiRegions::ExpListSharedPtr 
                &v_UpdateBndCondExpansion(int i)
            {
                return m_bndCondExpansions[i];
            }

            virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr> 
                &v_UpdateBndConditions()
            {
                return m_bndConditions;
            }

            virtual void v_EvaluateBoundaryConditions(
                const NekDouble time = 0.0,
                const NekDouble x2_in = NekConstants::kNekUnsetDouble,
                const NekDouble x3_in = NekConstants::kNekUnsetDouble);

            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo();
        };
        
        typedef boost::shared_ptr<DisContField2D>   DisContField2DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
