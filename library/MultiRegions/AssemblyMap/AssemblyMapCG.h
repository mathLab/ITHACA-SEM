///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapCG.h
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
// Description: C0-continuous Local to Global mapping routines, base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_ASSEMBLYMAPCG_H
#define MULTIREGIONS_ASSEMBLYMAPCG_H

#include <tuple>

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMap.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
    namespace MultiRegions
    {
        static std::map<int,int> NullIntIntMap;
        const static std::vector<std::map<int,int> > NullVecIntIntMap;

        class AssemblyMapCG;
        typedef std::shared_ptr<AssemblyMapCG>  AssemblyMapCGSharedPtr;
        typedef std::tuple<int, int, NekDouble> ExtraDirDof;

        typedef std::vector<std::map<int, int> > DofGraph;

        MULTI_REGIONS_EXPORT
        std::pair<int, StdRegions::Orientation> DeterminePeriodicEdgeOrientId(
            int                           meshEdgeId,
            StdRegions::Orientation       edgeOrient,
            const std::vector<PeriodicEntity> &periodicEdges);

        MULTI_REGIONS_EXPORT
        StdRegions::Orientation  DeterminePeriodicFaceOrient(
            StdRegions::Orientation   faceOrient1,
            StdRegions::Orientation   faceOrient2);


        /// Constructs mappings for the C0 scalar continuous Galerkin formulation.
        class AssemblyMapCG: public AssemblyMap
        {
            typedef Array<OneD, const ExpListSharedPtr> BndCondExp;
            typedef Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                BndCond;

        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT AssemblyMapCG(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::string variable = "DefaultVar");

            /// General constructor for expansions of all dimensions without
            /// boundary conditions.
            MULTI_REGIONS_EXPORT AssemblyMapCG(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int                                   numLocalCoeffs,
                const ExpList                              &locExp,
                const BndCondExp                           &bndCondExp
                                                    = NullExpListSharedPtrArray,
                const BndCond                              &bndConditions
                              = SpatialDomains::NullBoundaryConditionShPtrArray,
                const bool                                  checkIfSingular
                                                                        = false,
                const std::string                           variable
                                                                 = "defaultVar",
                const PeriodicMap                          &periodicVerts
                                                              = NullPeriodicMap,
                const PeriodicMap                          &periodicEdges
                                                              = NullPeriodicMap,
                const PeriodicMap                          &periodicFaces
                                                             = NullPeriodicMap);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~AssemblyMapCG();

            MULTI_REGIONS_EXPORT std::set<ExtraDirDof> &GetCopyLocalDirDofs()
            {
                return m_copyLocalDirDofs;
            }

            MULTI_REGIONS_EXPORT std::set<int> &GetParallelDirBndSign()
            {
                return m_parallelDirBndSign;
            }
            
        protected:
            /// Integer map of local coeffs to global space
            Array<OneD,int> m_localToGlobalMap;
            /// Integer sign of local coeffs to global space
            Array<OneD,NekDouble> m_localToGlobalSign;
            /// Bandwith of the full matrix system (no static condensation).
            int m_fullSystemBandWidth;
            /// Integer map of process coeffs to universal space
            Array<OneD,int> m_globalToUniversalMap;
            /// Integer map of unique process coeffs to universal space (signed)
            Array<OneD,int> m_globalToUniversalMapUnique;
            /// Number of non Dirichlet vertex modes
            int m_numNonDirVertexModes;
            /// Number of non Dirichlet edge modes
            int m_numNonDirEdgeModes;
            /// Number of non Dirichlet face modes
            int m_numNonDirFaceModes;
            /// Number of Dirichlet edges
            int m_numDirEdges;
            /// Number of Dirichlet faces
            int m_numDirFaces;
            /// Number of Dirichlet edges
            int m_numNonDirEdges;
            /// Number of Dirichlet faces
            int m_numNonDirFaces;
            /// Number of local boundary condition coefficients
            int m_numLocalBndCondCoeffs;
            /// Extra dirichlet edges in parallel
            Array<OneD, int> m_extraDirEdges;
            /// Number of local boundary condition degrees of freedom.
            int m_numLocDirBndCondDofs;
            /// Maximum static condensation level.
            int m_maxStaticCondLevel;
            /// Set indicating degrees of freedom which are Dirichlet but whose
            /// value is stored on another processor.
            std::set<ExtraDirDof> m_copyLocalDirDofs;
            /// Set indicating the local coeffs just touching parallel
            /// dirichlet boundary that have a sign change
            std::set<int> m_parallelDirBndSign;

            MULTI_REGIONS_EXPORT int CreateGraph(
                const ExpList                       &locExp,
                const BndCondExp                    &bndCondExp,
                const Array<OneD, const BndCond>    &bndConditions,
                const bool                           checkIfSystemSingular,
                const PeriodicMap                   &periodicVerts,
                const PeriodicMap                   &periodicEdges,
                const PeriodicMap                   &periodicFaces,
                DofGraph                            &graph,
                BottomUpSubStructuredGraphSharedPtr &bottomUpGraph,
                std::set<int>                       &extraDirVerts,
                std::set<int>                       &extraDirEdges,
                int                                 &firstNonDirGraphVertId,
                int                                 &nExtraDirichlet,
                int                                  mdswitch = 1);

            void SetUpUniversalC0ContMap(
                const ExpList     &locExp,
                const PeriodicMap &perVerts = NullPeriodicMap,
                const PeriodicMap &perEdges = NullPeriodicMap,
                const PeriodicMap &perFaces = NullPeriodicMap);

            /// Calculate the bandwith of the full matrix system.
            void CalculateFullSystemBandWidth();

            MULTI_REGIONS_EXPORT virtual int v_GetLocalToGlobalMap(const int i) const;

            MULTI_REGIONS_EXPORT virtual int v_GetGlobalToUniversalMap(const int i) const;

            MULTI_REGIONS_EXPORT virtual int v_GetGlobalToUniversalMapUnique(const int i) const;

            MULTI_REGIONS_EXPORT virtual const Array<OneD,const int>&  v_GetLocalToGlobalMap();

            MULTI_REGIONS_EXPORT virtual const Array<OneD, const int>& v_GetGlobalToUniversalMap();

            MULTI_REGIONS_EXPORT virtual const Array<OneD, const int>& v_GetGlobalToUniversalMapUnique();

            MULTI_REGIONS_EXPORT virtual NekDouble v_GetLocalToGlobalSign(const int i) const;

            MULTI_REGIONS_EXPORT virtual const Array<OneD, NekDouble>& v_GetLocalToGlobalSign() const;

            MULTI_REGIONS_EXPORT virtual void v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,       NekDouble>& global,
                    bool useComm) const;

            MULTI_REGIONS_EXPORT virtual void v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                    NekVector<      NekDouble>& global,
                    bool useComm) const;

            MULTI_REGIONS_EXPORT virtual void v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT virtual void v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT virtual void v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            MULTI_REGIONS_EXPORT virtual void v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            MULTI_REGIONS_EXPORT virtual void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT virtual void v_UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT virtual void v_UniversalAssemble(
                Array<OneD,     NekDouble>& pGlobal,
                int offset) const;

            MULTI_REGIONS_EXPORT virtual int v_GetFullSystemBandWidth() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumNonDirVertexModes() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumNonDirEdgeModes() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumNonDirFaceModes() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumDirEdges() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumDirFaces() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumNonDirEdges() const;

            MULTI_REGIONS_EXPORT virtual int v_GetNumNonDirFaces() const;

            MULTI_REGIONS_EXPORT virtual const Array<OneD, const int>& v_GetExtraDirEdges();

            MULTI_REGIONS_EXPORT virtual AssemblyMapSharedPtr v_LinearSpaceMap(
                const ExpList &locexp, GlobalSysSolnType solnType);
        };


    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_ASSEMBLYMAPCG_H

