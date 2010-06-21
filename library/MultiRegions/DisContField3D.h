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

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/LocalToGlobalDGMap.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/BoundaryConditions.h>
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
            DisContField3D();

            /// Constructs a global discontinuous field based on an input mesh.
            DisContField3D(SpatialDomains::MeshGraph3D &graph3D,
                           const GlobalSysSolnType solnType = eDirectStaticCond,
                           bool SetUpJustDG = true);

            /// Constructs a global discontinuous field based on an input mesh
            /// with boundary conditions.
            DisContField3D(SpatialDomains::MeshGraph3D &graph3D,
                           SpatialDomains::BoundaryConditions &bcs,
                           const int bc_loc,
                           const GlobalSysSolnType solnType = eDirectStaticCond,
                           bool SetUpJustDG = true);

            /// Constructs a global discontinuous field based on an input mesh
            /// with boundary conditions.
            DisContField3D(SpatialDomains::MeshGraph3D &graph3D,
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
            DisContField3D(const DisContField3D &In);

            /// Destructor.
            ~DisContField3D();

            inline const
                    Array<OneD,const MultiRegions::ExpList2DSharedPtr>&
                    GetBndCondExpansions();

            inline const
                    Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                    GetBndConditions();

            /// \brief Set up a list of element ids and edge ids the link to the
            /// boundary conditions
            void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                      Array<OneD,int> &FaceID);

        protected:
            /// The number of boundary segments on which Dirichlet boundary
            /// conditions are imposed.
            int m_numDirBndCondExpansions;

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
            Array<OneD,MultiRegions::ExpList2DSharedPtr>       m_bndCondExpansions;

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
//            void SetBoundaryConditionExpansion(
//                    SpatialDomains::MeshGraph3D &graph3D,
//                    SpatialDomains::BoundaryConditions &bcs,
//                    const std::string variable,
//                    Array<OneD, ExpList2DSharedPtr> &bndCondExpansions,
//                    Array<OneD, SpatialDomains::BoundaryConditionShPtr>
//                                                                &bndConditions);

            /// Evaluates boundary conditions.
            inline void EvaluateBoundaryConditions(
                    const NekDouble time = 0.0);


            virtual void v_GetBoundaryToElmtMap(Array<OneD,int> &ElmtID,
                                                Array<OneD,int> &FaceID)
            {
                GetBoundaryToElmtMap(ElmtID,FaceID);
            }

            /// \brief Set up an stl map containing the information
            /// for a robin aboundary condition in the location of the
            /// element id
            map<int, RobinBCInfoSharedPtr> GetRobinBCInfo(void);
        private:
            GlobalLinSysMapShPtr                                m_globalBndMat;
            ExpList2DSharedPtr                                  m_trace;
            LocalToGlobalDGMapSharedPtr                         m_traceMap;

            virtual void v_EvaluateBoundaryConditions(
                    const NekDouble time = 0.0);

            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& v_GetBndConditions()
            {
                return m_bndConditions;
            }

            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo()
            {
                return GetRobinBCInfo();
            }
        };

        typedef boost::shared_ptr<DisContField3D>   DisContField3DSharedPtr;

        inline const Array<OneD,const MultiRegions::ExpList2DSharedPtr>&
                    DisContField3D::GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }

        /**
         * \brief This function evaluates the boundary conditions at a certain
         * time-level.
         *
         * Based on the boundary condition \f$g(\boldsymbol{x},t)\f$ evaluated
         * at a given time-level \a t, this function transforms the boundary
         * conditions onto the coefficients of the (one-dimensional) boundary
         * expansion. Depending on the type of boundary conditions, these
         * expansion coefficients are calculated in different ways:
         * - <b>Dirichlet boundary conditions</b><BR>
         *   In order to ensure global \f$C^0\f$ continuity of the spectral/hp
         *   approximation, the Dirichlet boundary conditions are projected onto
         *   the boundary expansion by means of a modified \f$C^0\f$ continuous
         *   Galerkin projection. This projection can be viewed as a collocation
         *   projection at the vertices, followed by an \f$L^2\f$ projection on
         *   the interior modes of the edges. The resulting coefficients
         *   \f$\boldsymbol{\hat{u}}^{\mathcal{D}}\f$ will be stored for the
         *   boundary expansion.
         * - <b>Neumann boundary conditions</b>
         *   In the discrete Galerkin formulation of the problem to be solved,
         *   the Neumann boundary conditions appear as the set of surface
         *   integrals: \f[\boldsymbol{\hat{g}}=\int_{\Gamma}
         *   \phi^e_n(\boldsymbol{x})g(\boldsymbol{x})d(\boldsymbol{x})\quad
         *   \forall n \f]
         *   As a result, it are the coefficients \f$\boldsymbol{\hat{g}}\f$
         *   that will be stored in the boundary expansion
         *
         * \param time The time at which the boundary conditions should be
         * evaluated
         */
        inline void DisContField3D::EvaluateBoundaryConditions(const NekDouble time)
        {
            ExpList3D::EvaluateBoundaryConditions(time,m_bndCondExpansions,m_bndConditions);
        }

        inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                    DisContField3D::GetBndConditions()
        {
            return m_bndConditions;
        }

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3D_H
