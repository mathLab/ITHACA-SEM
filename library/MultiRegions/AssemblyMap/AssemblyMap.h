///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMap.h
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
// Description: Assembly (e.g. local to global) base mapping routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_ASSEMBLY_MAP_H
#define MULTIREGIONS_ASSEMBLY_MAP_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/SubStructuredGraph.h>
#include <vector>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Communication/GsLib.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class AssemblyMap;
        class ExpList;
        typedef std::shared_ptr<AssemblyMap>  AssemblyMapSharedPtr;
        static AssemblyMapSharedPtr NullAssemblyMapSharedPtr;

        /// Base class for constructing local to global mapping of degrees of
        /// freedom.
        class AssemblyMap
        {
        public:
        	/// Default constructor.
            MULTI_REGIONS_EXPORT AssemblyMap();
            /// Constructor with a communicator
            MULTI_REGIONS_EXPORT AssemblyMap(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::string variable = "DefaultVar");

            /// Constructor for next level in multi-level static condensation.
            MULTI_REGIONS_EXPORT AssemblyMap(AssemblyMap* oldLevelMap,
                    const BottomUpSubStructuredGraphSharedPtr& multiLevelGraph);
            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~AssemblyMap();

            /// Retrieves the communicator
            MULTI_REGIONS_EXPORT LibUtilities::CommSharedPtr GetComm();

            /// Retrieves the hash of this map
            MULTI_REGIONS_EXPORT size_t GetHash() const;

            MULTI_REGIONS_EXPORT int GetLocalToGlobalMap(const int i) const;

            MULTI_REGIONS_EXPORT int GetGlobalToUniversalMap(const int i) const;

            MULTI_REGIONS_EXPORT int GetGlobalToUniversalMapUnique(const int i) const;

            MULTI_REGIONS_EXPORT const Array<OneD,const int>&  GetLocalToGlobalMap();

            MULTI_REGIONS_EXPORT const Array<OneD, const int>& GetGlobalToUniversalMap();

            MULTI_REGIONS_EXPORT const Array<OneD, const int>& GetGlobalToUniversalMapUnique();

            MULTI_REGIONS_EXPORT NekDouble GetLocalToGlobalSign(const int i) const;

            MULTI_REGIONS_EXPORT const Array<OneD, NekDouble>& GetLocalToGlobalSign() const;

            MULTI_REGIONS_EXPORT void LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,       NekDouble>& global,
                    bool useComm = true) const;

            MULTI_REGIONS_EXPORT void LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                    NekVector<      NekDouble>& global,
                    bool useComm = true) const;

            MULTI_REGIONS_EXPORT void GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT void GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT void Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            MULTI_REGIONS_EXPORT void Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            MULTI_REGIONS_EXPORT void UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT void UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT void UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal,
                          int                         offset) const;

            MULTI_REGIONS_EXPORT void UniversalAbsMaxBnd(
                           Array<OneD, NekDouble> &bndvals);

            /// Retrieve the global index of a given local boundary mode.
            MULTI_REGIONS_EXPORT int GetLocalToGlobalBndMap(const int i) const;
            /// Retrieve the global indices of the local boundary modes.
            MULTI_REGIONS_EXPORT const Array<OneD,const int>&  GetLocalToGlobalBndMap();

            MULTI_REGIONS_EXPORT const Array<OneD, const int>& GetGlobalToUniversalBndMap();

            MULTI_REGIONS_EXPORT const Array<OneD, const int>& GetGlobalToUniversalBndMapUnique();

            /// Returns true if using a modal expansion requiring a change of
            /// sign of some modes.
            MULTI_REGIONS_EXPORT bool GetSignChange();

            /// Retrieve the sign change of a given local boundary mode.
            MULTI_REGIONS_EXPORT NekDouble GetLocalToGlobalBndSign(const int i) const;
            /// Retrieve the sign change for all local boundary modes.
            MULTI_REGIONS_EXPORT Array<OneD, const NekDouble> GetLocalToGlobalBndSign() const;
            /// Retrieves the local indices corresponding to the
            /// boundary expansion modes.
            MULTI_REGIONS_EXPORT const Array<OneD,const int>
                                              &GetBndCondCoeffsToLocalCoeffsMap();
            /// Returns the modal sign associated with a given
            /// boundary expansion mode.
            MULTI_REGIONS_EXPORT const Array<OneD, NekDouble>
                                              &GetBndCondCoeffsToLocalCoeffsSign();

            /// Retrieves the local indices corresponding to the
            /// boundary expansion modes to global trace
            MULTI_REGIONS_EXPORT const Array<OneD,const int>
                                              &GetBndCondCoeffsToLocalTraceMap();

            /// Returns the global index of the boundary trace giving the
            /// index on the boundary expansion
            MULTI_REGIONS_EXPORT int GetBndCondIDToGlobalTraceID(const int i);
            MULTI_REGIONS_EXPORT const Array<OneD, const int>
                &GetBndCondIDToGlobalTraceID();
 
            /// Returns the number of global Dirichlet boundary coefficients.
            MULTI_REGIONS_EXPORT int GetNumGlobalDirBndCoeffs() const;
            /// Returns the number of local Dirichlet boundary coefficients.
            MULTI_REGIONS_EXPORT int GetNumLocalDirBndCoeffs() const;
            /// Returns the total number of global boundary coefficients.
            MULTI_REGIONS_EXPORT int GetNumGlobalBndCoeffs() const;
            /// Returns the total number of local boundary coefficients.
            MULTI_REGIONS_EXPORT int GetNumLocalBndCoeffs() const;
            /// Returns the total number of local coefficients.
            MULTI_REGIONS_EXPORT int GetNumLocalCoeffs() const;
            /// Returns the total number of global coefficients.
            MULTI_REGIONS_EXPORT int GetNumGlobalCoeffs() const;
            /// Retrieves if the system is singular (true) or not (false)
            MULTI_REGIONS_EXPORT bool GetSingularSystem() const;

            ///
            MULTI_REGIONS_EXPORT void GlobalToLocalBnd(
                    const NekVector<NekDouble>& global,
                    NekVector<NekDouble>& loc,
                    int offset) const;

            MULTI_REGIONS_EXPORT void GlobalToLocalBnd(
                    const NekVector<NekDouble>& global,
                    NekVector<NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT void GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc,
                    int offset) const;

            MULTI_REGIONS_EXPORT void GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT void LocalBndToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,NekDouble>& global,
                    int offset, bool UseComm = true) const;

            MULTI_REGIONS_EXPORT void LocalBndToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,NekDouble>& global, bool UseComm = true)  const;

            MULTI_REGIONS_EXPORT void LocalToLocalBnd(
                    const Array<OneD, const NekDouble>& local,
                    Array<OneD,NekDouble>& locbnd)  const;

            MULTI_REGIONS_EXPORT void LocalToLocalInt(
                    const Array<OneD, const NekDouble>& local,
                    Array<OneD,NekDouble>& locint)  const;

            MULTI_REGIONS_EXPORT void LocalBndToLocal(
                    const Array<OneD, const NekDouble>& locbnd,
                    Array<OneD,NekDouble>& local)  const;

            MULTI_REGIONS_EXPORT void LocalIntToLocal(
                    const Array<OneD, const NekDouble>& locbnd,
                    Array<OneD,NekDouble>& local)  const;

            MULTI_REGIONS_EXPORT void AssembleBnd(const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global, int offset) const;

            MULTI_REGIONS_EXPORT void AssembleBnd(const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global) const;

            MULTI_REGIONS_EXPORT void AssembleBnd(const Array<OneD,const NekDouble>& loc,
                    Array<OneD, NekDouble>& global, int offset) const;

            MULTI_REGIONS_EXPORT void AssembleBnd(const Array<OneD, const NekDouble>& loc,
                    Array<OneD, NekDouble>& global) const;

            MULTI_REGIONS_EXPORT void UniversalAssembleBnd(
                          Array<OneD,     NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT void UniversalAssembleBnd(
                          NekVector<      NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT void UniversalAssembleBnd(
                          Array<OneD,     NekDouble>& pGlobal,
                          int                         offset) const;

            MULTI_REGIONS_EXPORT int GetFullSystemBandWidth() const;

            MULTI_REGIONS_EXPORT int GetNumNonDirVertexModes() const;

            MULTI_REGIONS_EXPORT int GetNumNonDirEdgeModes() const;

            MULTI_REGIONS_EXPORT int GetNumNonDirFaceModes() const;

            MULTI_REGIONS_EXPORT int GetNumDirEdges() const;

            MULTI_REGIONS_EXPORT int GetNumDirFaces() const;

            MULTI_REGIONS_EXPORT int GetNumNonDirEdges() const;

            MULTI_REGIONS_EXPORT int GetNumNonDirFaces() const;

            MULTI_REGIONS_EXPORT void PrintStats(std::ostream &out, std::string variable, bool printHeader = true) const;

            MULTI_REGIONS_EXPORT const Array<OneD, const int>& 
                GetExtraDirEdges();

            MULTI_REGIONS_EXPORT std::shared_ptr<AssemblyMap> LinearSpaceMap(const ExpList &locexp, GlobalSysSolnType solnType);

            /// Returns the bandwidth of the boundary system.
            MULTI_REGIONS_EXPORT int GetBndSystemBandWidth() const;
            /// Returns the level of static condensation for this map.
            MULTI_REGIONS_EXPORT int GetStaticCondLevel() const;
            /// Returns the number of patches in this static condensation level.
            MULTI_REGIONS_EXPORT int GetNumPatches() const;
            /// Returns the number of local boundary coefficients in each patch.
            MULTI_REGIONS_EXPORT const Array<OneD,const unsigned int>&
                    GetNumLocalBndCoeffsPerPatch();
            /// Returns the number of local interior coefficients in each patch.
            MULTI_REGIONS_EXPORT const Array<OneD,const unsigned int>&
                    GetNumLocalIntCoeffsPerPatch();
            /// Returns the local to global mapping for the next level in the
            /// multi-level static condensation.
            MULTI_REGIONS_EXPORT const AssemblyMapSharedPtr
                    GetNextLevelLocalToGlobalMap() const;

            MULTI_REGIONS_EXPORT void SetNextLevelLocalToGlobalMap( AssemblyMapSharedPtr  pNextLevelLocalToGlobalMap);

            /// Returns the patch map from the previous level 
            /// of the multi-level static condensation.
            MULTI_REGIONS_EXPORT const PatchMapSharedPtr&
                GetPatchMapFromPrevLevel(void) const;

            /// Returns true if this is the last level in the multi-level
            /// static condensation.
            MULTI_REGIONS_EXPORT bool AtLastLevel() const;
            /// Returns the method of solving global systems.
            MULTI_REGIONS_EXPORT GlobalSysSolnType GetGlobalSysSolnType() const;
            MULTI_REGIONS_EXPORT PreconditionerType GetPreconType() const;
            MULTI_REGIONS_EXPORT NekDouble GetIterativeTolerance() const;
            MULTI_REGIONS_EXPORT int GetMaxIterations() const;
            MULTI_REGIONS_EXPORT int GetSuccessiveRHS() const;

            MULTI_REGIONS_EXPORT int GetLowestStaticCondLevel() const
            {
                return m_lowestStaticCondLevel;
            }

            MULTI_REGIONS_EXPORT void PatchLocalToGlobal(
                               const Array<OneD, const NekDouble>& loc,
                               Array<OneD,       NekDouble>& global) const;

            MULTI_REGIONS_EXPORT void PatchGlobalToLocal(
                               const Array<OneD, const NekDouble>& global,
                               Array<OneD,       NekDouble>& loc) const;
                        
            MULTI_REGIONS_EXPORT void PatchAssemble(
                    const Array<OneD, const NekDouble> &loc,
                    Array<OneD,       NekDouble> &global) const;
        protected:
            /// Session object
            LibUtilities::SessionReaderSharedPtr m_session;

            /// Communicator
            LibUtilities::CommSharedPtr m_comm;

            /// Hash for map
            size_t m_hash;

            /// Number of local boundary coefficients
            int m_numLocalBndCoeffs;
            /// Total number of global boundary coefficients
            int m_numGlobalBndCoeffs;
            /// Number of Local Dirichlet Boundary Coefficients
            int m_numLocalDirBndCoeffs;
            /// Number of Global Dirichlet Boundary Coefficients
            int m_numGlobalDirBndCoeffs;
            /// Flag indicating if the system is singular or not
            bool m_systemSingular;

            /// Total number of local coefficients
            /** This corresponds to the number of total number of coefficients
             *  - For CG this corresponds to the total of bnd + int DOFs
             *  - For DG this corresponds to the number of bnd DOFs.
             *    This means that #m_numLocalCoeffs = #m_numLocalBndCoeffs
             *    This way, we can consider the trace-system solve as a
             *    statically condensed solve without interior DOFs. This allows
             *    us to use the same global system solver for both cases.
             */
            int m_numLocalCoeffs;

            /// Total number of global coefficients
            /** This corresponds to the number of total number of coefficients
             *  - For CG this corresponds to the total of bnd + int DOFs.
             *  - For DG this corresponds to the number of bnd DOFs.
             *    This means that #m_numGlobalCoeffs = #m_numGlobalBndCoeffs
             *    This way, we can consider the trace-system solve as a
             *    statically condensed solve without interior DOFs. This allows
             *    us to use the same global system solver for both cases.
             */
            int m_numGlobalCoeffs;

            /// Flag indicating if modes require sign reversal.
            bool m_signChange;

            /// Integer map of local coeffs to global Boundary Dofs
            Array<OneD,int>       m_localToGlobalBndMap; 
            /// Integer sign of local boundary coeffs to global space
            Array<OneD,NekDouble> m_localToGlobalBndSign;
            /// Integer map of local boundary coeffs to local boundary system numbering
            Array<OneD,int>       m_localToLocalBndMap;
            /// Integer map of local boundary coeffs to local interior system numbering
            Array<OneD,int>       m_localToLocalIntMap;
            /// Integer map of bnd cond coeffs to local coefficients
            Array<OneD,int>       m_bndCondCoeffsToLocalCoeffsMap;
            /// Integer map of sign of bnd cond coeffs to local coefficients
            Array<OneD,NekDouble> m_bndCondCoeffsToLocalCoeffsSign;
            /// Integer map of bnd cond coeff to local trace coeff
            Array<OneD,int>       m_bndCondCoeffsToLocalTraceMap;
            /// Integer map of bnd cond trace number to global trace number
            Array<OneD,int>       m_bndCondIDToGlobalTraceID;
            /// Integer map of process coeffs to universal space
            Array<OneD,int>       m_globalToUniversalBndMap;
            /// Integer map of unique process coeffs to universal space (signed)
            Array<OneD,int>       m_globalToUniversalBndMapUnique;

            /// The solution type of the global system
            GlobalSysSolnType m_solnType;
            /// The bandwith of the global bnd system
            int m_bndSystemBandWidth;

            /// Type type of preconditioner to use in iterative solver.
            PreconditionerType m_preconType;

            /// Maximum iterations for iterative solver
            int m_maxIterations;

            /// Tolerance for iterative solver
            NekDouble  m_iterativeTolerance;

            /// sucessive RHS  for iterative solver
            int  m_successiveRHS;

            Gs::gs_data * m_gsh;
            Gs::gs_data * m_bndGsh;
            /// gs gather communication to impose Dirhichlet BCs. 
            Gs::gs_data * m_dirBndGsh; 
            
            /// The level of recursion in the case of multi-level static
            /// condensation.
            int m_staticCondLevel;
            /// The number of patches (~elements) in the current level
            int m_numPatches;
            /// The number of bnd dofs per patch
            Array<OneD, unsigned int> m_numLocalBndCoeffsPerPatch;
            /// The number of int dofs per patch
            Array<OneD, unsigned int> m_numLocalIntCoeffsPerPatch;
            /// Map from the patches of the previous level to the patches of
            /// the current level

            /// The local to global mapping of the next level of recursion
            AssemblyMapSharedPtr m_nextLevelLocalToGlobalMap;
            /// Lowest static condensation level.
            int m_lowestStaticCondLevel;
            
            /// Calculates the bandwidth of the boundary system.
            void CalculateBndSystemBandWidth();

            void GlobalToLocalBndWithoutSign(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc);

        private:
            /// Mapping information for previous level in MultiLevel Solver
            PatchMapSharedPtr m_patchMapFromPrevLevel;

            virtual int v_GetLocalToGlobalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMap(const int i) const;

            virtual int v_GetGlobalToUniversalMapUnique(const int i) const;

            virtual const Array<OneD,const int>&  v_GetLocalToGlobalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMap();

            virtual const Array<OneD, const int>& v_GetGlobalToUniversalMapUnique();

            virtual NekDouble v_GetLocalToGlobalSign(const int i) const;

            virtual const Array<OneD, NekDouble>& v_GetLocalToGlobalSign() const;

            virtual void v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,       NekDouble>& global,
                    bool useComm) const;

            virtual void v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                    NekVector<      NekDouble>& global,
                    bool useComm) const;

            virtual void v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            virtual void v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            virtual void v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            virtual void v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            virtual void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            virtual void v_UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            virtual void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal,
                          int                         offset) const;

            virtual int v_GetFullSystemBandWidth() const;

            virtual int v_GetNumNonDirVertexModes() const;

            virtual int v_GetNumNonDirEdgeModes() const;

            virtual int v_GetNumNonDirFaceModes() const;

            virtual int v_GetNumDirEdges() const;

            virtual int v_GetNumDirFaces() const;

            virtual int v_GetNumNonDirEdges() const;

            virtual int v_GetNumNonDirFaces() const;
            
            virtual const Array<OneD, const int>& 
                v_GetExtraDirEdges();
            
            /// Generate a linear space mapping from existing mapping 
            virtual std::shared_ptr<AssemblyMap> v_LinearSpaceMap(
                const ExpList &locexp, GlobalSysSolnType solnType);
        };


    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_ASSEMBLY_MAP_H


