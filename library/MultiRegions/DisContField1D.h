///////////////////////////////////////////////////////////////////////////////
//
// File DisContField1D.h
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
// Description: Field definition in one-dimension for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD1D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>
#include <SpatialDomains/Conditions.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>


namespace Nektar
{
    namespace MultiRegions
    {
        /// This class is the abstraction of a global discontinuous two-
        /// dimensional spectral/hp element expansion which approximates the
        /// solution of a set of partial differential equations.
        class DisContField1D: public ExpList1D
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT DisContField1D();

            /// Constructs a 1D discontinuous field based on a mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT DisContField1D(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const std::string &variable,
                const bool SetUpJustDG  = true,
                const Collections::ImplementationType ImpType
                = Collections::eNoImpType);
            
            /// Constructor for a DisContField1D from a List of subdomains
            /// New Constructor for arterial network 
            MULTI_REGIONS_EXPORT DisContField1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const SpatialDomains::CompositeMap& domain,
                const SpatialDomains::BoundaryConditions &Allbcs, 
                const std::string &variable,
                bool SetToOneSpaceDimensions = false,
                const Collections::ImplementationType ImpType
                = Collections::eNoImpType);

            /// Constructs a 1D discontinuous field based on an existing field.
            MULTI_REGIONS_EXPORT DisContField1D(const DisContField1D &In);
            
            /// Constructs a 1D discontinuous field based on an existing field.
	    /// (needed in order to use ContField( const ExpList1D &In) constructor
            MULTI_REGIONS_EXPORT DisContField1D(const ExpList1D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~DisContField1D();
            
            /// For a given key, returns the associated global linear system.
            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr GetGlobalBndLinSys(
                const GlobalLinSysKey &mkey);


            // Return the internal vector which directs whether the normal flux
            // at the trace defined by Left and Right Adjacent elements
            // is negated with respect to the segment normal
            MULTI_REGIONS_EXPORT std::vector<bool> &GetNegatedFluxNormal(void);

        protected:
            /// The number of boundary segments on which Dirichlet boundary
            /// conditions are imposed.
            int m_numDirBndCondExpansions;

            /// Discretised boundary conditions.
            /**
             * It is an array of size equal to the number of boundary points
             * and consists of entries of the type LocalRegions#PointExp. Every
             * entry corresponds to a point on a single boundary region.
             */
            Array<OneD,MultiRegions::ExpListSharedPtr>      m_bndCondExpansions;

            Array<OneD,NekDouble> m_bndCondBndWeight;
            
            /// An array which contains the information about the boundary
            /// condition on the different boundary regions.
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            /// Global boundary matrix.
            GlobalLinSysMapShPtr                               m_globalBndMat;
            
            /// Trace space storage for points between elements.
            ExpListSharedPtr                                   m_trace;

            /// Local to global DG mapping for trace space.
            AssemblyMapDGSharedPtr                        m_traceMap;

            /**
             * @brief A set storing the global IDs of any boundary edges.
             */
            std::set<int> m_boundaryVerts;

            
            /**
             * @brief A map which identifies groups of periodic vertices.
             */
            PeriodicMap m_periodicVerts;

            
            /**
             * @brief A vector indicating degress of freedom which need to be
             * copied from forwards to backwards space in case of a periodic
             * boundary condition.
             */
            std::vector<int> m_periodicFwdCopy;
            std::vector<int> m_periodicBwdCopy;


            /*
             * @brief A map identifying which verts are left- and right-adjacent
             * for DG.
             */
            std::vector<bool> m_leftAdjacentVerts;


            /// Discretises the boundary conditions.
            void GenerateBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string variable);
            
            
            /// Generate a associative map of periodic vertices in a mesh.
            void FindPeriodicVertices(const SpatialDomains::BoundaryConditions &bcs,
                                      const std::string variable);
            
            virtual ExpListSharedPtr &v_GetTrace()
            {
                return m_trace;
            }
            
            virtual AssemblyMapDGSharedPtr &v_GetTraceMap(void)
            {
                return m_traceMap;
            }
            
            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fn,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_AddTraceQuadPhysToField(
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                    Array<OneD,       NekDouble> &field);
            virtual void v_ExtractTracePhys(
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_ExtractTracePhys(
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,       NekDouble> &outarray);

            /// Populates the list of boundary condition expansions.
            void SetBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string variable,
                Array<OneD, MultiRegions::ExpListSharedPtr>
                    &bndCondExpansions,
                Array<OneD, SpatialDomains
                    ::BoundaryConditionShPtr> &bndConditions);
            
            /// Populates the list of boundary condition expansions in multidomain case.
            void SetMultiDomainBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string variable,
                Array<OneD, MultiRegions::ExpListSharedPtr>
                    &bndCondExpansions,
                Array<OneD, SpatialDomains
                    ::BoundaryConditionShPtr> &bndConditions,
                int subdomain);
            
            void GenerateFieldBnd1D(
                SpatialDomains::BoundaryConditions &bcs,
                const std::string variable);
            
            virtual std::map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo();
            
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

            virtual void v_GetBoundaryToElmtMap(
                Array<OneD,int> &ElmtID, Array<OneD,int> &VertID);
            virtual void v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays);
            virtual void v_Reset();

            /// Evaluate all boundary conditions at a given time..
            virtual void v_EvaluateBoundaryConditions(
                const NekDouble   time    = 0.0,
                const std::string varName = "",
                const NekDouble   x2_in   = NekConstants::kNekUnsetDouble,
                const NekDouble   x3_in   = NekConstants::kNekUnsetDouble);
            
            /// Solve the Helmholtz equation.
            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const MultiRegions::VarFactorsMap &varfactors,
                    const Array<OneD, const NekDouble> &dirForcing,
                    const bool PhysSpaceForcing);
            
            virtual void v_FillBwdWithBound(
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_FillBwdWithBoundDeriv(
                const int                          Dir,
                const Array<OneD, const NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_FillBwdWithBwdWeight(
                    Array<OneD,       NekDouble> &weightave,
                    Array<OneD,       NekDouble> &weightjmp);

            virtual void v_GetFwdBwdTracePhysInterior(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD,       NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            
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

            inline virtual void v_GetFwdBwdTracePhysDerivSerial(
                const int                          Dir,
                const Array<OneD, const NekDouble> &field,
                    Array<OneD,       NekDouble> &Fwd,
                    Array<OneD,       NekDouble> &Bwd);

            inline virtual void v_GetFwdBwdTracePhysSerial(
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

        private:
            void SetUpDG(const std::string &variable);
            
            bool IsLeftAdjacentVertex(const int n, const int e);

            std::vector<bool> m_negatedFluxNormal;

            SpatialDomains::BoundaryConditionsSharedPtr GetDomainBCs(const SpatialDomains::CompositeMap &domain,
                                                                     const SpatialDomains::BoundaryConditions &Allbcs,
                                                                     const std::string &variable);
        };

        void DisContField1D::v_GetFwdBwdTracePhys(
            Array<OneD, NekDouble> &Fwd,
            Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(m_phys, Fwd, Bwd);
        }
        
        void DisContField1D::v_GetFwdBwdTracePhys(
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysSerial(field, Fwd, Bwd);
            m_traceMap->GetAssemblyCommDG()->PerformExchange(Fwd, Bwd);
        }

        void DisContField1D::v_GetFwdBwdTracePhysNoBndFill(
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysInterior(field, Fwd, Bwd);
            m_traceMap->GetAssemblyCommDG()->PerformExchange(Fwd, Bwd);
        }

        void DisContField1D::v_GetFwdBwdTracePhysDeriv(
            const int                          Dir,
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysDerivSerial(Dir,field, Fwd, Bwd);
            m_traceMap->GetAssemblyCommDG()->PerformExchange(Fwd, Bwd);
        }

        void DisContField1D::v_GetFwdBwdTracePhysDerivSerial(
            const int                          Dir,
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysInterior(field, Fwd, Bwd);
            v_FillBwdWithBoundDeriv(Dir, Fwd, Bwd);
        }

        void DisContField1D::v_GetFwdBwdTracePhysSerial(
            const Array<OneD, const NekDouble> &field,
                Array<OneD,       NekDouble> &Fwd,
                Array<OneD,       NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhysInterior(field, Fwd, Bwd);
            v_FillBwdWithBound(Fwd, Bwd);
        }

        void DisContField1D::v_PeriodicBwdCopy(
            const Array<OneD, const NekDouble> &Fwd,
                    Array<OneD,       NekDouble> &Bwd)
        {
            for (int n = 0; n < m_periodicFwdCopy.size(); ++n)
            {
                Bwd[m_periodicBwdCopy[n]] = Fwd[m_periodicFwdCopy[n]];
            }
        }

        const Array<OneD,const NekDouble>
            &DisContField1D::v_GetBndCondBwdWeight()
        {
            return m_bndCondBndWeight;
        }

        void DisContField1D::v_SetBndCondBwdWeight(
            const int index, 
            const NekDouble value)
        {
            m_bndCondBndWeight[index]   =   value;
        }

        typedef std::shared_ptr<DisContField1D>   DisContField1DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD1D_H
