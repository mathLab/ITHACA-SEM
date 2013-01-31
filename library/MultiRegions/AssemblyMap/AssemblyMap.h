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
        class AssemblyMap;
        typedef boost::shared_ptr<AssemblyMap>  AssemblyMapSharedPtr;
        static AssemblyMapSharedPtr NullAssemblyMapSharedPtr;



        /// Base class for constructing local to global mapping of degrees of
        /// freedom.
        class AssemblyMap
        {
        public:
        	/// Default constructor.
            MULTI_REGIONS_EXPORT AssemblyMap();
            /// Constructor with a communicator
            MULTI_REGIONS_EXPORT AssemblyMap(const LibUtilities::SessionReaderSharedPtr &pSession);

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

            MULTI_REGIONS_EXPORT const void LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const;

            MULTI_REGIONS_EXPORT const void LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            MULTI_REGIONS_EXPORT const void GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT const void GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            MULTI_REGIONS_EXPORT const void Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            MULTI_REGIONS_EXPORT const void Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            MULTI_REGIONS_EXPORT const void UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT const void UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT const void UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal,
                          int                         offset) const;

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
            NekDouble GetLocalToGlobalBndSign(const int i) const;
            /// Retrieve the sign change for all local boundary modes.
            Array<OneD, const NekDouble> GetLocalToGlobalBndSign() const;
            /// Retrieves the global index corresponding to a boundary expansion
            /// mode.
            MULTI_REGIONS_EXPORT int GetBndCondCoeffsToGlobalCoeffsMap(const int i);
            /// Retrieves the global indices corresponding to the boundary
            /// expansion modes.
            MULTI_REGIONS_EXPORT const Array<OneD,const int>&
                    GetBndCondCoeffsToGlobalCoeffsMap();
            /// Returns the modal sign associated with a given boundary
            /// expansion mode.
            MULTI_REGIONS_EXPORT NekDouble GetBndCondCoeffsToGlobalCoeffsSign(const int i);

            /// Returns the global index of the boundary trace giving the
            /// index on the boundary  expansion
            MULTI_REGIONS_EXPORT int GetBndCondTraceToGlobalTraceMap(const int i);
 
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
                    const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global,
                    int offset) const;

            MULTI_REGIONS_EXPORT void LocalBndToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,NekDouble>& global,
                    int offset) const;

            MULTI_REGIONS_EXPORT void AssembleBnd(const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global, int offset) const;

            MULTI_REGIONS_EXPORT void AssembleBnd(const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global) const;

            MULTI_REGIONS_EXPORT void AssembleBnd(const Array<OneD,const NekDouble>& loc,
                    Array<OneD, NekDouble>& global, int offset) const;

            MULTI_REGIONS_EXPORT void AssembleBnd(const Array<OneD, const NekDouble>& loc,
                    Array<OneD, NekDouble>& global) const;

            MULTI_REGIONS_EXPORT const void UniversalAssembleBnd(
                          Array<OneD,     NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT const void UniversalAssembleBnd(
                          NekVector<      NekDouble>& pGlobal) const;

            MULTI_REGIONS_EXPORT const void UniversalAssembleBnd(
                          Array<OneD,     NekDouble>& pGlobal,
                          int                         offset) const;

            MULTI_REGIONS_EXPORT const int GetFullSystemBandWidth() const;

            MULTI_REGIONS_EXPORT int GetNumNonDirVertexModes() const;

            MULTI_REGIONS_EXPORT int GetNumNonDirEdgeModes() const;

            MULTI_REGIONS_EXPORT int GetNumNonDirFaceModes() const;

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
            MULTI_REGIONS_EXPORT const GlobalSysSolnType  GetGlobalSysSolnType() const;
            MULTI_REGIONS_EXPORT const PreconditionerType  GetPreconType() const;

            MULTI_REGIONS_EXPORT const int GetLowestStaticCondLevel() const
            {
                return m_lowestStaticCondLevel;
            }

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

            /// Integer map of local boundary coeffs to global space
            Array<OneD,int>       m_localToGlobalBndMap;
            /// Integer sign of local boundary coeffs to global space
            Array<OneD,NekDouble> m_localToGlobalBndSign;
            /// Integer map of bnd cond coeffs to global coefficients
            Array<OneD,int>       m_bndCondCoeffsToGlobalCoeffsMap;
            /// Integer map of bnd cond coeffs to global coefficients
            Array<OneD,NekDouble> m_bndCondCoeffsToGlobalCoeffsSign;
            /// Integer map of bnd cond trace number to global trace number
            Array<OneD,int>       m_bndCondTraceToGlobalTraceMap;
            /// Integer map of process coeffs to universal space
            Array<OneD,int>       m_globalToUniversalBndMap;
            /// Integer map of unique process coeffs to universal space (signed)
            Array<OneD,int>       m_globalToUniversalBndMapUnique;

            /// The solution type of the global system
            GlobalSysSolnType m_solnType;
            /// The bandwith of the global bnd system
            int m_bndSystemBandWidth;

            PreconditionerType m_preconType;

            Gs::gs_data * m_gsh;
            Gs::gs_data * m_bndGsh;

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

            virtual const void v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const;

            virtual const void v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            virtual const void v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            virtual const void v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            virtual const void v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            virtual const void v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            virtual const void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            virtual const void v_UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;

            virtual const void v_UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal,
                          int                         offset) const;

            virtual const int v_GetFullSystemBandWidth() const;

            virtual int v_GetNumNonDirVertexModes() const;

            virtual int v_GetNumNonDirEdgeModes() const;

            virtual int v_GetNumNonDirFaceModes() const;

        };


    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_ASSEMBLY_MAP_H


