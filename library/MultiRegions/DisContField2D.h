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
// Description: Field definition in two-dimensions for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/Conditions.h>
#include <SpatialDomains/SegGeom.h>

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
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const std::string &variable,
                const bool SetUpJustDG = true,
                const bool DeclareCoeffPhysArrays = true);
            
            MULTI_REGIONS_EXPORT DisContField2D(
                const DisContField2D &In,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const std::string &variable,
                const bool SetUpJustDG = false,
                const bool DeclareCoeffPhysArrays = true);
            
            MULTI_REGIONS_EXPORT DisContField2D(
                const DisContField2D &In,
                const bool DeclareCoeffPhysArrays = true);

            MULTI_REGIONS_EXPORT virtual ~DisContField2D();

            /// Determines if another ContField2D shares the same boundary
            /// conditions as this field.
            MULTI_REGIONS_EXPORT bool SameTypeOfBoundaryConditions(const DisContField2D &In);

            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr GetGlobalBndLinSys(const GlobalLinSysKey &mkey);

            MULTI_REGIONS_EXPORT NekDouble L2_DGDeriv(const int dir,
                                 const Array<OneD, const NekDouble> &soln);

            /// \brief Set up an stl map containing the information
            /// for a robin aboundary condition in the location of the
            /// element id
            MULTI_REGIONS_EXPORT void EvaluateHDGPostProcessing(Array<OneD, NekDouble> &outarray);

            /// Generates a map of periodic edges in the mesh.
            MULTI_REGIONS_EXPORT void GetPeriodicEdges(
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string &variable,
                vector<map<int,int> > & periodicVertices,
                map<int,int>& periodicEdges);
            
        protected:
            /**
             * \brief An object which contains the discretised
             * boundary conditions.
             *
             * It is an array of size equal to the number of boundary regions
             * and consists of entries of the type
             * MultiRegions#ExpList1D. Every entry corresponds to the
             * one-dimensional spectral/hp expansion on a single boundary
             * region.  The values of the boundary conditions are stored as
             * the coefficients of the one-dimensional expansion.
             */
            Array<OneD,MultiRegions::ExpListSharedPtr>       m_bndCondExpansions;              
            

            /**
             * \brief An array which contains the information about
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
             * @brief Auxiliary map for periodic boundary conditions.
             * 
             * Takes geometry IDs of periodic edges to a pair (n,e), where n
             * is the expansion containing the edge and e the local edge number.
             */
            boost::unordered_map<int,pair<int,int> > m_perEdgeToExpMap;

            /**
             * \brief This function discretises the boundary conditions by setting up
             * a list of one-dimensional boundary expansions.
             *
             * According to their boundary region, the separate segmental boundary
             * expansions are bundled together in an object of the class
             * MultiRegions#ExpList1D.
             * The list of expansions of the Dirichlet boundary regions are listed
             * first in the array #m_bndCondExpansions.
             *
             * \param graph2D A mesh, containing information about the domain and
             * the spectral/hp element expansion.
             * \param bcs An entity containing information about the boundaries and
             * boundary conditions.
             * \param variable An optional parameter to indicate for which variable
             * the boundary conditions should be discretised.
             */
            void GenerateBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string &variable,
                const bool DeclareCoeffPhysArrays = true);

            virtual ExpListSharedPtr &v_GetTrace()
            {
                return m_trace;
            }

            virtual AssemblyMapDGSharedPtr &v_GetTraceMap()
            {
                return m_traceMap;
            }
            
            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fx,
                const Array<OneD, const NekDouble> &Fy,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fn,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_AddTraceBiIntegral(
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                      Array<OneD,       NekDouble> &outarray);

            /**
             * \brief This method extracts the "forward" and "backward" trace
             * data from the array \a field and puts the data into output
             * vectors \a Fwd and \a Bwd.
             * 
             * We first define the convention which defines "forwards" and
             * "backwards". First an association is made between the edge of
             * each element and its corresponding edge in the trace space
             * using the mapping #m_traceMap. The element can either be
             * left-adjacent or right-adjacent to this trace edge (see
             * Expansion1D::GetLeftAdjacentElementExp). Boundary edges are
             * always left-adjacent since left-adjacency is populated first.
             * 
             * If the element is left-adjacent we extract the edge trace data
             * from \a field into the forward trace space \a Fwd; otherwise,
             * we place it in the backwards trace space \a Bwd. In this way,
             * we form a unique set of trace normals since these are always
             * extracted from left-adjacent elements.
             *
             * \param field is a NekDouble array which contains the 2D data
             * from which we wish to extract the backward and forward
             * orientated trace/edge arrays.
             *
             * \return Updates a NekDouble array \a Fwd and \a Bwd
             */
            virtual void v_GetFwdBwdTracePhys(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD,       NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_GetFwdBwdTracePhys(
                      Array<OneD,       NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);

            /**
             * \brief This method extracts the trace (edges in 2D)
             * from the field \a inarray and puts the values in \a outarray.
             *
             * It assumes the field is C0 continuous so that it can overwrite
             * the edge data when visited by the two adjacent elements.
             *
             * \param inarray is a NekDouble array which contains the 2D
             * data from which we wish to extract the edge data
             *
             * \return Updates a NekDouble array \a outarray which
             * contains the edge information
             */
            virtual void v_ExtractTracePhys(
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_ExtractTracePhys(
                      Array<OneD,       NekDouble> &outarray);

            virtual void v_GetPeriodicEdges(const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                            const SpatialDomains::BoundaryConditions &bcs,
                                            const std::string &variable,
                                            vector<map<int,int> > & periodicVertices,
                                            map<int,int>& periodicEdges)
            {
                GetPeriodicEdges(graph2D,bcs,variable,periodicVertices,periodicEdges);
            }

            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const FlagList &flags,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const Array<OneD, const NekDouble> &dirForcing);

            /// Calculates the result of the multiplication of a global
            /// matrix of type specified by \a mkey with a vector given by \a
            /// inarray.
            virtual void v_GeneralMatrixOp(
                   const GlobalMatrixKey             &gkey,
                   const Array<OneD,const NekDouble> &inarray,
                   Array<OneD,      NekDouble> &outarray,
                   bool UseContCoeffs = false);

            virtual void v_GetBoundaryToElmtMap(
                Array<OneD, int> &ElmtID,
                Array<OneD, int> &EdgeID);

            virtual const Array<OneD,const MultiRegions::ExpListSharedPtr>
                &v_GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }

            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                &v_GetBndConditions()
            {
                return m_bndConditions;
            }

            virtual MultiRegions::ExpListSharedPtr &v_UpdateBndCondExpansion(int i)
            {
                return m_bndCondExpansions[i];
            }

            virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr> &v_UpdateBndConditions()
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
