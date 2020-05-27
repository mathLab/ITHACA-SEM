///////////////////////////////////////////////////////////////////////////////
//
// File DisContField3DHomogeneous1D.h
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
// Description: Field definition in three-dimensions for a discontinuous
// LDG-H expansion with a homogeneous direction in 1D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3DHOMO1D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3DHOMO1D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SpatialDomains/Conditions.h>
#include <MultiRegions/GlobalLinSys.h>


namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField3DHomogeneous1D: public ExpList3DHomogeneous1D
        {
        public:
            MULTI_REGIONS_EXPORT DisContField3DHomogeneous1D();

            MULTI_REGIONS_EXPORT DisContField3DHomogeneous1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const LibUtilities::BasisKey &HomoBasis,
                const NekDouble lhom,
                const bool useFFT,
                const bool dealiasing);

            MULTI_REGIONS_EXPORT DisContField3DHomogeneous1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const LibUtilities::BasisKey               &HomoBasis,
                const NekDouble                             lhom,
                const bool                                  useFFT,
                const bool                                  dealiasing,
                const SpatialDomains::MeshGraphSharedPtr   &graph2D,
                const std::string                          &variable,
                const Collections::ImplementationType ImpType
                = Collections::eNoImpType);
                                                       

            /// Copy constructor.
            MULTI_REGIONS_EXPORT DisContField3DHomogeneous1D(
                const DisContField3DHomogeneous1D &In,
                const bool DeclarePlanesSetCoeffPhys = true);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~DisContField3DHomogeneous1D();

            MULTI_REGIONS_EXPORT void SetupBoundaryConditions(
                const LibUtilities::BasisKey       &HomoBasis,
                const NekDouble                     lhom,
                SpatialDomains::BoundaryConditions &bcs,
                const std::string                   variable);

            /**
             * \brief This function evaluates the boundary conditions
             * at a certaintime-level.
             *
             * Based on the boundary condition \f$g(\boldsymbol{x},t)\f$
             * evaluated at a given time-level \a t, this function transforms
             * the boundary conditions onto the coefficients of the
             * (one-dimensional) boundary expansion.
             * Depending on the type of boundary conditions, these expansion
             * coefficients are calculated in different ways:
             * - <b>Dirichlet boundary conditions</b><BR>
             *   In order to ensure global \f$C^0\f$ continuity of the
             *   spectral/hp approximation, the Dirichlet boundary conditions
             *   are projected onto the boundary expansion by means of a
             *   modified \f$C^0\f$ continuous Galerkin projection.
             *   This projection can be viewed as a collocation projection at the
             *   vertices, followed by an \f$L^2\f$ projection on the interior
             *   modes of the edges. The resulting coefficients
             *   \f$\boldsymbol{\hat{u}}^{\mathcal{D}}\f$ will be stored for the
             *   boundary expansion.
             * - <b>Neumann boundary conditions</b>
             *   In the discrete Galerkin formulation of the problem to be
             *   solved, the Neumann boundary conditions appear as the set of
             *   surface integrals: \f[\boldsymbol{\hat{g}}=\int_{\Gamma}
             *   \phi^e_n(\boldsymbol{x})g(\boldsymbol{x})d(\boldsymbol{x})\quad
             *   \forall n \f]
             *   As a result, it are the coefficients \f$\boldsymbol{\hat{g}}\f$
             *   that will be stored in the boundary expansion
             *
             * \param time The time at which the boundary conditions should be
             * evaluated
             */
            MULTI_REGIONS_EXPORT void EvaluateBoundaryConditions(
                const NekDouble   time    = 0.0,
                const std::string varName = "");

            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>
                &GetBndCondExpansions();

            inline const Array<OneD,const SpatialDomains::
                BoundaryConditionShPtr> &GetBndConditions();

            inline std::shared_ptr<ExpList> &UpdateBndCondExpansion(int i);

            inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>&
                UpdateBndConditions();

            /// \brief Set up a list of element ids and edge ids the link to the
            /// boundary conditions
            MULTI_REGIONS_EXPORT void GetBoundaryToElmtMap(
                Array<OneD, int> &ElmtID,
                Array<OneD,int> &EdgeID);
            
            virtual void v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays);

            /// This funtion extract form a vector containing a full
            /// 3D-homogenous-1D field the value associated with a
            /// specific boundary conditions.
            /// TotField is the full field contaning all the physical values
            /// BndVals is the vector where the boundary physical values are
            /// stored BndID is the identifier of the boundary region
            MULTI_REGIONS_EXPORT void GetBCValues(
                      Array<OneD, NekDouble> &BndVals,
                const Array<OneD, NekDouble> &TotField,
                int                           BndID);

            /// This function calculate the inner product of two vectors
            /// (V1 and V2)  respect to the basis along a boundary region.
            /// outarray is the inner product result multiplied by the normal to
            /// the edge (specified by the BndID)
            MULTI_REGIONS_EXPORT void NormVectorIProductWRTBase(
                Array<OneD, const NekDouble> &V1,
                Array<OneD, const NekDouble> &V2,
                Array<OneD,       NekDouble> &outarray,
                int                           BndID);

            /// Storage space for the boundary to element and boundary to trace
            /// map. This member variable is really allocated just in case
            /// a boundary expansion recasting is required at the solver level.
            /// Otherwise is the 2 vectors are not filled up.
            /// If is needed all the funcitons whihc require to use this map
            /// do not have to recalculate it anymore.
            Array<OneD, int> m_BCtoElmMap;
            Array<OneD, int> m_BCtoEdgMap;

        protected:
            /**
             * \brief An object which contains the discretised
             * boundary conditions.
             *
             * It is an array of size equal to the number of boundary
             * regions and consists of entries of the type
             * MultiRegions#ExpList1D. Every entry corresponds to the
             * one-dimensional spectral/hp expansion on a single
             * boundary region.  The values of the boundary conditions
             * are stored as the coefficients of the one-dimensional
             * expansion.
             */

            Array<OneD, MultiRegions::ExpListSharedPtr>  m_bndCondExpansions;

            Array<OneD, NekDouble >                      m_bndCondBndWeight;

            ExpListSharedPtr m_trace;

            Array<OneD, int> m_traceBndMap;

            /**
             * \brief An array which contains the information about
             * the boundary condition on the different boundary
             * regions.
             */
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            virtual void v_GetBoundaryToElmtMap(
                Array<OneD,int> &ElmtID,
                Array<OneD,int> &EdgeID)
            {
                GetBoundaryToElmtMap(ElmtID, EdgeID);
            }

            virtual void v_GetBCValues(
                      Array<OneD, NekDouble> &BndVals,
                const Array<OneD, NekDouble> &TotField,
                int                           BndID)
            {
                GetBCValues(BndVals, TotField, BndID);
            }

            virtual void v_NormVectorIProductWRTBase(
                Array<OneD, const NekDouble> &V1,
                Array<OneD, const NekDouble> &V2,
                Array<OneD,       NekDouble> &outarray,
                int                           BndID)
            {
                NormVectorIProductWRTBase(V1,V2,outarray,BndID);
            }

            /// Set up all DG member variables and maps
            MULTI_REGIONS_EXPORT void SetUpDG();

            /// @todo Fix in another way considering all the planes
            virtual ExpListSharedPtr &v_GetTrace()
            {
                return m_trace;
            }

            /// @todo Fix in another way considering all the planes
            virtual AssemblyMapDGSharedPtr &v_GetTraceMap()
            {
                return m_planes[0]->GetTraceMap();
            }

            virtual const Array<OneD,const MultiRegions::ExpListSharedPtr>
                &v_GetBndCondExpansions(void)
            {
                return GetBndCondExpansions();
            }

            virtual const
                Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                &v_GetBndConditions()
           {
               return GetBndConditions();
           }

           /// @todo Fix Robin BCs for homogeneous case
           virtual std::map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo()
           {
               return std::map<int, RobinBCInfoSharedPtr>();
           }

           virtual void v_ExtractTracePhys(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

           virtual void v_ExtractTracePhys(
                      Array<OneD,       NekDouble> &outarray);
           
           virtual void v_GetBoundaryNormals(int i,
                            Array<OneD, Array<OneD, NekDouble> > &normals);

        private:
            // virtual functions
            virtual void v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const StdRegions::ConstFactorMap   &factors,
                const StdRegions::VarCoeffMap      &varcoeff,
                const MultiRegions::VarFactorsMap &varfactors,
                const Array<OneD, const NekDouble> &dirForcing,
                const bool PhysSpaceForcing);

            virtual void v_EvaluateBoundaryConditions(
                const NekDouble   time    = 0.0,
                const std::string varName = "",
                const NekDouble   x2_in   = NekConstants::kNekUnsetDouble,
                const NekDouble   x3_in   = NekConstants::kNekUnsetDouble);

            virtual std::shared_ptr<ExpList> &v_UpdateBndCondExpansion(int i);

            virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                &v_UpdateBndConditions();

            virtual const Array<OneD, const int> &v_GetTraceBndMap()
            {
                return m_traceBndMap;
            }
            inline virtual void v_SetBndCondBwdWeight(
                const int index, 
                const NekDouble value);
        };

        typedef std::shared_ptr<DisContField3DHomogeneous1D>
            DisContField3DHomogeneous1DSharedPtr;

        inline const Array<OneD,const MultiRegions::ExpListSharedPtr>
            &DisContField3DHomogeneous1D::GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }

        inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
            &DisContField3DHomogeneous1D::GetBndConditions()
        {
            return m_bndConditions;
        }

        inline MultiRegions::ExpListSharedPtr &DisContField3DHomogeneous1D::
            UpdateBndCondExpansion(int i)
        {
            return m_bndCondExpansions[i];
        }

        inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>
            &DisContField3DHomogeneous1D::UpdateBndConditions()
        {
            return m_bndConditions;
        }
        inline void DisContField3DHomogeneous1D::v_SetBndCondBwdWeight(
            const int index, 
            const NekDouble value)
        {
            m_bndCondBndWeight[index]   =   value;
        }
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD3DHOMO1D_H
