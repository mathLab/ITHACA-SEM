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
#include <MultiRegions/LocalToGlobalDGMap.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/Conditions.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Abstraction of a global discontinuous three-dimensional spectral/hp
        /// element expansion which approximates the solution of a set of
        /// partial differential equations.
        class DisContField3D: public ExpList3D
        {
        public:
            /// Default constructor
            MULTI_REGIONS_EXPORT DisContField3D();

            /// Constructs a global discontinuous field based on an input mesh.
            MULTI_REGIONS_EXPORT DisContField3D(
                           LibUtilities::CommSharedPtr &pComm,
                           SpatialDomains::MeshGraph3D &graph3D,
                           const GlobalSysSolnType solnType = eDirectStaticCond,
                           bool SetUpJustDG = true);

            /// Constructs a global discontinuous field based on an input mesh
            /// with boundary conditions.
            MULTI_REGIONS_EXPORT DisContField3D(
                           LibUtilities::CommSharedPtr &pComm,
                           SpatialDomains::MeshGraph3D &graph3D,
                           SpatialDomains::BoundaryConditions &bcs,
                           const int bc_loc,
                           const GlobalSysSolnType solnType = eDirectStaticCond,
                           bool SetUpJustDG = true);

            /// Constructs a global discontinuous field based on an input mesh
            /// with boundary conditions.
            MULTI_REGIONS_EXPORT DisContField3D(
                           LibUtilities::CommSharedPtr &pComm,
                           SpatialDomains::MeshGraph3D &graph3D,
                           SpatialDomains::BoundaryConditions &bcs,
                           const std::string variable,
                           const GlobalSysSolnType solnType = eDirectStaticCond,
                           bool SetUpJustDG = true);

            /// Constructs a global discontinuous field based on an input mesh.
/*            DisContField3D(SpatialDomains::MeshGraph3D &graph3D,
                          SpatialDomains::BoundaryConditions &bcs,
                          const GlobalSysSolnType solnType = eDirectStaticCond,
                          const bool constructMap = true);
*/

            /// Constructs a global discontinuous field based on another
            /// discontinuous field.
            MULTI_REGIONS_EXPORT DisContField3D(const DisContField3D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT ~DisContField3D();

            inline const
                    Array<OneD,const MultiRegions::ExpListSharedPtr>&
                    GetBndCondExpansions();

            inline const
                    Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                    GetBndConditions();

            /// \brief Set up a list of element ids and edge ids the link to the
            /// boundary conditions
            MULTI_REGIONS_EXPORT void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                      Array<OneD,int> &FaceID);

        protected:

            /// An object which contains the discretised boundary conditions.
            /**
             * It is an array of size equal to the number of boundary
             * regions and consists of entries of the type
             * MultiRegions#ExpList1D. Every entry corresponds to the
             * one-dimensional spectral/hp expansion on a single
             * boundary region.  The values of the boundary conditions
             * are stored as the coefficients of the one-dimensional
             * expansion.
             */
            Array<OneD,MultiRegions::ExpListSharedPtr>       m_bndCondExpansions;

            /// An array which contains the information about the boundary
            /// condition on the different boundary regions.
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            /// This function discretises the boundary conditions by setting up
            /// a list of one-dimensional boundary expansions.
            void GenerateBoundaryConditionExpansion(
                    SpatialDomains::MeshGraph3D &graph3D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable);

            /// Populates the list of boundary condition expansions.
            void SetBoundaryConditionExpansion(
                    SpatialDomains::MeshGraph3D &graph3D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable,
                    Array<OneD, ExpListSharedPtr> &bndCondExpansions,
                    Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                    &bndConditions);

            /// Generates a map of periodic faces in the mesh.
            void GetPeriodicFaces(SpatialDomains::MeshGraph3D &graph3D,
                                  SpatialDomains::BoundaryConditions &bcs,
                                  const std::string variable,
                                  map<int,int>& periodicVertices,
                                  map<int,int>& periodicEdges,
                                  map<int,int>& periodicFaces);

            /// \brief Set up an stl map containing the information
            /// for a robin aboundary condition in the location of the
            /// element id
            map<int, RobinBCInfoSharedPtr> GetRobinBCInfo(void);

        private:
            GlobalLinSysMapShPtr                                m_globalBndMat;
            ExpList2DSharedPtr                                  m_trace;
            LocalToGlobalDGMapSharedPtr                         m_traceMap;

            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0,
													  const NekDouble x2_in = NekConstants::kNekUnsetDouble,
													  const NekDouble x3_in = NekConstants::kNekUnsetDouble);


            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& v_GetBndConditions()
            {
                return m_bndConditions;
            }

            virtual void v_GetBoundaryToElmtMap(Array<OneD,int> &ElmtID,
                                                Array<OneD,int> &FaceID);

            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo()
            {
                return GetRobinBCInfo();
            }
        };

        typedef boost::shared_ptr<DisContField3D>   DisContField3DSharedPtr;

        inline const Array<OneD,const MultiRegions::ExpListSharedPtr>&
            DisContField3D::GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }

        inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
            DisContField3D::GetBndConditions()
        {
            return m_bndConditions;
        }

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3D_H
