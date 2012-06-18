///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapBase.h
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

#ifndef MULTIREGIONS_ASSEMBLY_MAP_BASE_H
#define MULTIREGIONS_ASSEMBLY_MAP_BASE_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/SubStructuredGraph.h>
#include <vector>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Communication/GsLib.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        class AssemblyMapBase;
        typedef boost::shared_ptr<AssemblyMapBase>  AssemblyMapBaseSharedPtr;
        static AssemblyMapBaseSharedPtr NullAssemblyMapBaseSharedPtr;

        /// Base class for constructing local to global mapping of degrees of
        /// freedom.
        class AssemblyMapBase
        {
        public:
        	/// Default constructor.
            MULTI_REGIONS_EXPORT AssemblyMapBase();
            /// Constructor with a communicator
            MULTI_REGIONS_EXPORT AssemblyMapBase(const LibUtilities::SessionReaderSharedPtr &pSession);

            /// Constructor for next level in multi-level static condensation.
            MULTI_REGIONS_EXPORT AssemblyMapBase(AssemblyMapBase* oldLevelMap,
                    const BottomUpSubStructuredGraphSharedPtr& multiLevelGraph);
            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~AssemblyMapBase();

            /// Retrieves the communicator
            inline LibUtilities::CommSharedPtr GetComm();

            /// Retrieves the hash of this map
            inline size_t GetHash() const;

            inline int GetLocalToGlobalMap(const int i) const;

            inline int GetGlobalToUniversalMap(const int i) const;

            inline int GetGlobalToUniversalMapUnique(const int i) const;

            inline const Array<OneD,const int>&  GetLocalToGlobalMap();

            inline const Array<OneD, const int>& GetGlobalToUniversalMap();

            inline const Array<OneD, const int>& GetGlobalToUniversalMapUnique();

            inline NekDouble GetLocalToGlobalSign(const int i) const;

            inline const Array<OneD, NekDouble>& GetLocalToGlobalSign() const;

            inline const void LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const;

            inline const void LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            inline const void GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const;

            inline const void GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const;

            inline const void Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const;

            inline const void Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const;

            inline const void UniversalAssemble(
                          Array<OneD,     NekDouble>& pGlobal) const;

            inline const void UniversalAssemble(
                          NekVector<      NekDouble>& pGlobal) const;


            /// Retrieve the global index of a given local boundary mode.
            inline int GetLocalToGlobalBndMap(const int i) const;
            /// Retrieve the global indices of the local boundary modes.
            inline const Array<OneD,const int>&  GetLocalToGlobalBndMap();

            inline const Array<OneD, const int>& GetGlobalToUniversalBndMap();

            inline const Array<OneD, const int>& GetGlobalToUniversalBndMapUnique();

            /// Returns true if using a modal expansion requiring a change of
            /// sign of some modes.
            inline bool GetSignChange();            

            /// Retrieve the sign change of a given local boundary mode.
            inline NekDouble GetLocalToGlobalBndSign(const int i) const;
            /// Retrieve the sign change for all local boundary modes.
            inline Array<OneD, const NekDouble> GetLocalToGlobalBndSign() const;
            /// Retrieves the global index corresponding to a boundary expansion
            /// mode.
            inline int GetBndCondCoeffsToGlobalCoeffsMap(const int i);
            /// Retrieves the global indices corresponding to the boundary
            /// expansion modes.
            inline const Array<OneD,const int>&
                    GetBndCondCoeffsToGlobalCoeffsMap();
            /// Returns the modal sign associated with a given boundary
            /// expansion mode.
            inline NekDouble GetBndCondCoeffsToGlobalCoeffsSign(const int i);

            /// Returns the global index of the boundary trace giving the
            /// index on the boundary  expansion
            inline int GetBndCondTraceToGlobalTraceMap(const int i);
 
            /// Returns the number of global Dirichlet boundary coefficients.
            inline int GetNumGlobalDirBndCoeffs() const;
            /// Returns the number of local Dirichlet boundary coefficients.
            inline int GetNumLocalDirBndCoeffs() const;
            /// Returns the total number of global boundary coefficients.
            inline int GetNumGlobalBndCoeffs() const;
            /// Returns the total number of local boundary coefficients.
            inline int GetNumLocalBndCoeffs() const;
            /// Returns the total number of local coefficients.
            inline int GetNumLocalCoeffs() const;
            /// Returns the total number of global coefficients.
            inline int GetNumGlobalCoeffs() const;
            /// Retrieves if the system is singular (true) or not (false)
            inline bool GetSingularSystem() const;

            ///
            inline void GlobalToLocalBnd(
                    const NekVector<NekDouble>& global,
                    NekVector<NekDouble>& loc,
                    int offset) const;

            inline void GlobalToLocalBnd(
                    const NekVector<NekDouble>& global,
                    NekVector<NekDouble>& loc) const;

            inline void GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc,
                    int offset) const;

            inline void GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc) const;

            inline void LocalBndToGlobal(
                    const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global,
                    int offset) const;

            inline void LocalBndToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,NekDouble>& global,
                    int offset) const;

            inline void AssembleBnd(const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global, int offset) const;

            inline void AssembleBnd(const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global) const;

            inline void AssembleBnd(const Array<OneD,const NekDouble>& loc,
                    Array<OneD, NekDouble>& global, int offset) const;

            inline void AssembleBnd(const Array<OneD, const NekDouble>& loc,
                    Array<OneD, NekDouble>& global) const;

            inline const void UniversalAssembleBnd(
                          Array<OneD,     NekDouble>& pGlobal) const;

            inline const void UniversalAssembleBnd(
                          NekVector<      NekDouble>& pGlobal) const;

            inline const int GetFullSystemBandWidth() const;

            inline int GetNumNonDirVertexModes() const;

            inline int GetNumNonDirEdgeModes() const;

            inline int GetNumNonDirFaceModes() const;

            /// Returns the bandwidth of the boundary system.
            inline int GetBndSystemBandWidth() const;
            /// Returns the level of static condensation for this map.
            inline int GetStaticCondLevel() const;
            /// Returns the number of patches in this static condensation level.
            inline int GetNumPatches() const;
            /// Returns the number of local boundary coefficients in each patch.
            inline const Array<OneD,const unsigned int>&
                    GetNumLocalBndCoeffsPerPatch();
            /// Returns the number of local interior coefficients in each patch.
            inline const Array<OneD,const unsigned int>&
                    GetNumLocalIntCoeffsPerPatch();
            /// Returns the local to global mapping for the next level in the
            /// multi-level static condensation.
            inline const AssemblyMapBaseSharedPtr
                    GetNextLevelLocalToGlobalMap() const;

            inline void SetNextLevelLocalToGlobalMap( AssemblyMapBaseSharedPtr  pNextLevelLocalToGlobalMap);

            /// Returns the patch map from the previous level 
            /// of the multi-level static condensation.
            inline const PatchMapSharedPtr&
                GetPatchMapFromPrevLevel(void) const;

            /// Returns true if this is the last level in the multi-level
            /// static condensation.
            inline bool AtLastLevel() const;
            /// Returns the method of solving global systems.
            inline const GlobalSysSolnType  GetGlobalSysSolnType() const;
            inline const PreconditionerType  GetPreconType() const;
            
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
            AssemblyMapBaseSharedPtr m_nextLevelLocalToGlobalMap;

            /// Calculates the bandwidth of the boundary system.
            void CalculateBndSystemBandWidth();

            inline void GlobalToLocalBndWithoutSign(
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

            virtual const int v_GetFullSystemBandWidth() const;

            virtual int v_GetNumNonDirVertexModes() const;

            virtual int v_GetNumNonDirEdgeModes() const;

            virtual int v_GetNumNonDirFaceModes() const;

        };


    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_ASSEMBLY_MAP_BASE_H


