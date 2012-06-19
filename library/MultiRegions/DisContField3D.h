///////////////////////////////////////////////////////////////////////////////
//
// File DisContField3D.h
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
// Description: Field definition in three-dimensions for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/Conditions.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * Abstraction of a global discontinuous three-dimensional spectral/hp
         * element expansion which approximates the solution of a set of
         * partial differential equations.
         */
        class DisContField3D : public ExpList3D
        {
        public:
            /**
             * Default constructor
             */
            MULTI_REGIONS_EXPORT DisContField3D();

            /**
             * Constructs a global discontinuous field based on an input mesh
             * with boundary conditions.
             */
            MULTI_REGIONS_EXPORT DisContField3D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph3D,
                const std::string &variable,
                const bool SetUpJustDG = true);

            MULTI_REGIONS_EXPORT DisContField3D(
                const DisContField3D &In,
                const SpatialDomains::MeshGraphSharedPtr &graph3D,
                const std::string &variable,
                const bool SetUpJustDG = false);
            
            /// Constructs a global discontinuous field based on another
            /// discontinuous field.
            MULTI_REGIONS_EXPORT DisContField3D(const DisContField3D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~DisContField3D();

        protected:
            /**
             * An array of size equal to the number of boundary regions and
             * consists of entries of the type MultiRegions#ExpList1D. Every
             * entry corresponds to the one-dimensional spectral/hp expansion
             * on a single boundary region.  The values of the boundary
             * conditions are stored as the coefficients of the
             * one-dimensional expansion.
             */
            Array<OneD,MultiRegions::ExpListSharedPtr>        m_bndCondExpansions;

            /// An array which contains the information about the boundary
            /// condition on the different boundary regions.
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            GlobalLinSysMapShPtr        m_globalBndMat;
            ExpListSharedPtr            m_trace;
            AssemblyMapDGSharedPtr      m_traceMap;

            /// This function discretises the boundary conditions by setting up
            /// a list of one-dimensional boundary expansions.
            void GenerateBoundaryConditionExpansion(
                    const SpatialDomains::MeshGraphSharedPtr &graph3D,
                    const SpatialDomains::BoundaryConditions &bcs,
                    const std::string &variable);

            bool SameTypeOfBoundaryConditions(const DisContField3D &In);

            /// Populates the list of boundary condition expansions.
            void SetBoundaryConditionExpansion(
                    const SpatialDomains::MeshGraphSharedPtr &graph3D,
                    const SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable,
                    Array<OneD, ExpListSharedPtr> &bndCondExpansions,
                    Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                    &bndConditions);

            /// Generates a map of periodic faces in the mesh.
            void GetPeriodicFaces(const SpatialDomains::MeshGraphSharedPtr &graph3D,
                                  const SpatialDomains::BoundaryConditions &bcs,
                                  const std::string &variable,
                                  map<int,int>& periodicVertices,
                                  map<int,int>& periodicEdges,
                                  map<int,int>& periodicFaces);

            virtual void v_EvaluateBoundaryConditions(
                const NekDouble time = 0.0,
                const NekDouble x2_in = NekConstants::kNekUnsetDouble,
                const NekDouble x3_in = NekConstants::kNekUnsetDouble);

            /**
             * \brief This method extracts the "forward" and "backward" trace
             * data from the array \a field and puts the data into output
             * vectors \a Fwd and \a Bwd.
             * 
             * We first define the convention which defines "forwards" and
             * "backwards". First an association is made between the face of
             * each element and its corresponding face in the trace space
             * using the mapping #m_traceMap. The element can either be
             * left-adjacent or right-adjacent to this trace face (see
             * Expansion2D::GetLeftAdjacentElementExp). Boundary faces are
             * always left-adjacent since left-adjacency is populated first.
             * 
             * If the element is left-adjacent we extract the face trace data
             * from \a field into the forward trace space \a Fwd; otherwise,
             * we place it in the backwards trace space \a Bwd. In this way,
             * we form a unique set of trace normals since these are always
             * extracted from left-adjacent elements.
             *
             * \param field is a NekDouble array which contains the 3D data
             * from which we wish to extract the backward and forward
             * orientated trace/face arrays.
             *
             * \return Updates a NekDouble array \a Fwd and \a Bwd
             */
            virtual void v_GetFwdBwdTracePhys(
                Array<OneD,NekDouble> &Fwd,
                Array<OneD,NekDouble> &Bwd);
            virtual void v_GetFwdBwdTracePhys(
                const Array<OneD,const NekDouble> &field,
                      Array<OneD,      NekDouble> &Fwd,
                      Array<OneD,      NekDouble> &Bwd);
            virtual void v_ExtractTracePhys(
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_ExtractTracePhys(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fn,
                      Array<OneD,       NekDouble> &outarray);

            virtual const Array<OneD, const MultiRegions::ExpListSharedPtr> 
                &v_GetBndCondExpansions();
            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr> 
                &v_GetBndConditions();

            /// \brief Set up a list of element ids and edge ids the link to the
            /// boundary conditions
            virtual void v_GetBoundaryToElmtMap(Array<OneD,int> &ElmtID,
                                                Array<OneD,int> &FaceID);

            virtual ExpListSharedPtr &v_GetTrace(void)
            {
                return m_trace;
            }
            
            virtual AssemblyMapDGSharedPtr &v_GetTraceMap()
            {
                return m_traceMap;
            }

            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo();
        };

        typedef boost::shared_ptr<DisContField3D> DisContField3DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3D_H
