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

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
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
                const bool DeclareCoeffPhysArrays = true,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);

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

            virtual ExpListSharedPtr &v_GetTrace()
            {
                if(m_trace == NullExpListSharedPtr)
                {
                    SetUpDG();
                }

                return m_trace;
            }

            Array<OneD, int> m_BCtoElmMap;
            Array<OneD, int> m_BCtoEdgMap;

        protected:

            Array<OneD, LibUtilities::BasisSharedPtr> m_base; /**< Bases needed for the expansion */

            /** \brief This function gets the shared point to basis
             *
             *  \return returns the shared pointer to the bases
             */
            inline const Array<OneD, const LibUtilities::BasisSharedPtr>& GetBase() const
            {
                return(m_base);
            }

            /** \brief This function returns the type of basis used in the \a dir
             *  direction
             *
             *  The different types of bases implemented in the code are defined
             *  in the LibUtilities::BasisType enumeration list. As a result, the
             *  function will return one of the types of this enumeration list.
             *
             *  \param dir the direction
             *  \return returns the type of basis used in the \a dir direction
             */
            inline  LibUtilities::BasisType GetBasisType(const int dir) const
            {
                ASSERTL1(dir < m_base.size(), "dir is larger than m_numbases");
                return(m_base[dir]->GetBasisType());
            }

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
             * @brief A Double to determine the weight of Bwd trace value when averaging 
             *        the Fwd and Bwd in diffusion flux calculation.
             *
             * It is necessary because for some boundaries the Bwd itself is the target value,
             * while for other boundaries some average of Fwd and Bwd is the target value.
             * Using this parameter they can be treated uniformly.
             * This array stores a Double value for each variable for each ExpListSharedPtr, which 
             * assumes that for each ExpListSharedPtr the same kind of boundary treatment is used.
             *  
             */
            Array<OneD,NekDouble> m_BndCondBwdWeight;

            /**
             * @brief An array which contains the information about
             * the boundary condition on the different boundary regions.
             */
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            GlobalLinSysMapShPtr   m_globalBndMat;
            ExpListSharedPtr       m_trace;
            AssemblyMapDGSharedPtr m_traceMap;

            /**
             * Map of local trace (the points at the face of
             * the element) to the trace space discretisation
             */
            LocTraceToTraceMapSharedPtr m_locTraceToTraceMap;

            Array<OneD, Array<OneD, unsigned int> > m_mapEdgeToElmn;
            Array<OneD, Array<OneD, unsigned int> > m_signEdgeToElmn;
            Array<OneD,StdRegions::Orientation>     m_edgedir;

            /**
             * @brief A set storing the global IDs of any boundary edges.
             */
            std::set<int> m_boundaryEdges;

            /**
             * @brief A map which identifies groups of periodic vertices.
             */
            PeriodicMap m_periodicVerts;

            /**
             * @brief A map which identifies pairs of periodic edges.
             */
            PeriodicMap m_periodicEdges;


            /**
             * @brief A vector indicating degress of freedom which need to be
             * copied from forwards to backwards space in case of a periodic
             * boundary condition.
             */
            std::vector<int> m_periodicFwdCopy;
            std::vector<int> m_periodicBwdCopy;

            /*
             * @brief A map identifying which edges are left- and right-adjacent
             * for DG.
             */
            std::vector<bool> m_leftAdjacentEdges;

            void SetUpDG(const std::string  = "DefaultVar");
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
                const MultiRegions::VarFactorsMap &varfactors,
                const Array<OneD, const NekDouble> &dirForcing,
                const bool  PhysSpaceForcing);
            virtual void v_GeneralMatrixOp(
                const GlobalMatrixKey             &gkey,
                const Array<OneD,const NekDouble> &inarray,
                      Array<OneD,      NekDouble> &outarray,
                CoeffState                         coeffstate = eLocal);
            virtual void v_GetBoundaryToElmtMap(
                Array<OneD, int> &ElmtID,
                Array<OneD, int> &EdgeID);
            virtual void v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays);
            virtual void v_Reset();

            /**
             * @brief Obtain a copy of the periodic edges and vertices for this
             * field.
             */
            virtual void v_GetPeriodicEntities(
                PeriodicMap &periodicVerts,
                PeriodicMap &periodicEdges,
                PeriodicMap &periodicFaces)
            {
                boost::ignore_unused(periodicFaces);
                periodicVerts = m_periodicVerts;
                periodicEdges = m_periodicEdges;
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
                const NekDouble   time    = 0.0,
                const std::string varName = "",
                const NekDouble   x2_in   = NekConstants::kNekUnsetDouble,
                const NekDouble   x3_in   = NekConstants::kNekUnsetDouble);

            virtual std::map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo();

            virtual void v_FillBwdWithBound(
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_FillBwdWithBoundDeriv(
                const int                          Dir,
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_FillBwdWITHBwdWeight(
                    Array<OneD,       NekDouble> &weightave,
                    Array<OneD,       NekDouble> &weightjmp);

            virtual void v_GetFwdBwdTracePhysInterior(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD,       NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            
            virtual void v_AddTraceQuadPhysToField(
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                    Array<OneD,       NekDouble> &field);
            
            inline virtual void v_GetFwdBwdTracePhys(
                Array<OneD, NekDouble> &Fwd,
                Array<OneD, NekDouble> &Bwd);
            
            inline virtual void v_GetFwdBwdTracePhys(
                const Array<OneD, const NekDouble> &field,
                    Array<OneD,       NekDouble> &Fwd,
                    Array<OneD,       NekDouble> &Bwd);

            inline virtual void v_GetFwdBwdTracePhysNoBndFill(
                const Array<OneD, const NekDouble> &field,
                    Array<OneD,       NekDouble> &Fwd,
                    Array<OneD,       NekDouble> &Bwd);

            inline virtual void v_GetFwdBwdTracePhysDeriv(
                const int                          Dir,
                const Array<OneD, const NekDouble> &field,
                    Array<OneD,       NekDouble> &Fwd,
                    Array<OneD,       NekDouble> &Bwd);

            inline virtual void v_GetFwdBwdTracePhysDeriv_serial(
                const int                          Dir,
                const Array<OneD, const NekDouble> &field,
                    Array<OneD,       NekDouble> &Fwd,
                    Array<OneD,       NekDouble> &Bwd);

            inline virtual void v_GetFwdBwdTracePhys_serial(
                const Array<OneD, const NekDouble> &field,
                    Array<OneD,       NekDouble> &Fwd,
                    Array<OneD,       NekDouble> &Bwd);

            inline virtual void v_PeriodicBwdCopy(
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);

            inline virtual const Array<OneD,const NekDouble>
                &v_GetBndCondBwdWeight();

            inline virtual void v_SetBndCondBwdWeight(
                const int index, 
                const NekDouble value);

            inline virtual const LocTraceToTraceMapSharedPtr 
                    &v_GetlocTraceToTraceMap() const;
        };

        void DisContField2D::v_GetFwdBwdTracePhys(
            Array<OneD, NekDouble> &Fwd,
            Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(m_phys, Fwd, Bwd);
        }
        
        void DisContField2D::v_GetFwdBwdTracePhys(
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys_serial(field, Fwd, Bwd);
            m_traceMap->UniversalTraceAssemble(Fwd);
            m_traceMap->UniversalTraceAssemble(Bwd);
        }

        void DisContField2D::v_GetFwdBwdTracePhysNoBndFill(
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysInterior(field, Fwd, Bwd);
            m_traceMap->UniversalTraceAssemble(Fwd);
            m_traceMap->UniversalTraceAssemble(Bwd);
        }

        void DisContField2D::v_GetFwdBwdTracePhysDeriv(
            const int                          Dir,
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysDeriv_serial(Dir,field, Fwd, Bwd);

            m_traceMap->UniversalTraceAssemble(Fwd);
            m_traceMap->UniversalTraceAssemble(Bwd);
        }

        void DisContField2D::v_GetFwdBwdTracePhysDeriv_serial(
            const int                          Dir,
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysInterior(field, Fwd, Bwd);
            v_FillBwdWithBoundDeriv(Dir, Fwd, Bwd);
        }

        void DisContField2D::v_GetFwdBwdTracePhys_serial(
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysInterior(field, Fwd, Bwd);
            v_FillBwdWithBound(Fwd, Bwd);
        }

        void DisContField2D::v_PeriodicBwdCopy(
            const Array<OneD, const NekDouble> &Fwd,
                    Array<OneD,       NekDouble> &Bwd)
        {
            for (int n = 0; n < m_periodicFwdCopy.size(); ++n)
            {
                Bwd[m_periodicBwdCopy[n]] = Fwd[m_periodicFwdCopy[n]];
            }
        }

        const Array<OneD,const NekDouble>
            &DisContField2D::v_GetBndCondBwdWeight()
        {
            return m_BndCondBwdWeight;
        }

        void DisContField2D::v_SetBndCondBwdWeight(
            const int index, 
            const NekDouble value)
        {
            m_BndCondBwdWeight[index]   =   value;
        }

        const LocTraceToTraceMapSharedPtr 
                &DisContField2D::v_GetlocTraceToTraceMap() const
        {
            return m_locTraceToTraceMap;
        }

        typedef std::shared_ptr<DisContField2D>   DisContField2DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
