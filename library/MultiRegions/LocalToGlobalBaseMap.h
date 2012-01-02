///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalBaseMap.h
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
// Description: Local to Global base mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_LOCALTOGLOBALBASEMAP_H
#define MULTIREGIONS_LOCALTOGLOBALBASEMAP_H
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
        class LocalToGlobalBaseMap;
        typedef boost::shared_ptr<LocalToGlobalBaseMap>  LocalToGlobalBaseMapSharedPtr;
        static LocalToGlobalBaseMapSharedPtr NullLocalToGlobalBaseMapSharedPtr;

        /// Base class for constructing local to global mapping of degrees of
        /// freedom.
        class LocalToGlobalBaseMap
        {
        public:
        	/// Default constructor.
            MULTI_REGIONS_EXPORT LocalToGlobalBaseMap();
            /// Constructor with a communicator
            MULTI_REGIONS_EXPORT LocalToGlobalBaseMap(const LibUtilities::SessionReaderSharedPtr &pSession);

            /// Constructor for next level in multi-level static condensation.
            MULTI_REGIONS_EXPORT LocalToGlobalBaseMap(LocalToGlobalBaseMap* oldLevelMap,
                    const BottomUpSubStructuredGraphSharedPtr& multiLevelGraph);
            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~LocalToGlobalBaseMap();

            /// Retrieves the communicator
            inline LibUtilities::CommSharedPtr GetComm();

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
            inline const LocalToGlobalBaseMapSharedPtr
                    GetNextLevelLocalToGlobalMap() const;

            inline void SetNextLevelLocalToGlobalMap( LocalToGlobalBaseMapSharedPtr  pNextLevelLocalToGlobalMap);

            /// Returns the patch map from the previous level 
            /// of the multi-level static condensation.
            inline const PatchMapSharedPtr&
                GetPatchMapFromPrevLevel(void) const;

            /// Returns true if this is the last level in the multi-level
            /// static condensation.
            inline bool AtLastLevel() const;
            /// Returns the method of solving global systems.
            inline const GlobalSysSolnType  GetGlobalSysSolnType() const;
            
        protected:
            /// Session object
            LibUtilities::SessionReaderSharedPtr m_session;

            /// Communicator
            LibUtilities::CommSharedPtr m_comm;

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
            LocalToGlobalBaseMapSharedPtr m_nextLevelLocalToGlobalMap;

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

        };


        inline LibUtilities::CommSharedPtr LocalToGlobalBaseMap::GetComm()
        {
            return m_comm;
        }

        inline int LocalToGlobalBaseMap::GetLocalToGlobalMap(const int i) const
        {
            return v_GetLocalToGlobalMap(i);
        }

        inline int LocalToGlobalBaseMap::GetGlobalToUniversalMap(const int i) const
        {
            return v_GetGlobalToUniversalMap(i);
        }

        inline int LocalToGlobalBaseMap::GetGlobalToUniversalMapUnique(const int i) const
        {
            return v_GetGlobalToUniversalMapUnique(i);
        }

        inline const Array<OneD,const int>&  LocalToGlobalBaseMap::GetLocalToGlobalMap()
        {
            return v_GetLocalToGlobalMap();
        }

        inline const Array<OneD, const int>& LocalToGlobalBaseMap::GetGlobalToUniversalMap()
        {
            return v_GetGlobalToUniversalMap();
        }

        inline const Array<OneD, const int>& LocalToGlobalBaseMap::GetGlobalToUniversalMapUnique()
        {
            return v_GetGlobalToUniversalMapUnique();
        }

        inline NekDouble LocalToGlobalBaseMap::GetLocalToGlobalSign(const int i) const
        {
            return v_GetLocalToGlobalSign(i);
        }

        inline const Array<OneD, NekDouble>& LocalToGlobalBaseMap::GetLocalToGlobalSign() const
        {
            return v_GetLocalToGlobalSign();
        }

        inline const void LocalToGlobalBaseMap::LocalToGlobal(
                const Array<OneD, const NekDouble>& loc,
                      Array<OneD,       NekDouble>& global) const
        {
            v_LocalToGlobal(loc,global);
        }

        inline const void LocalToGlobalBaseMap::LocalToGlobal(
                const NekVector<NekDouble>& loc,
                      NekVector<      NekDouble>& global) const
        {
            v_LocalToGlobal(loc,global);
        }

        inline const void LocalToGlobalBaseMap::GlobalToLocal(
                const Array<OneD, const NekDouble>& global,
                      Array<OneD,       NekDouble>& loc) const
        {
            v_GlobalToLocal(global,loc);
        }

        inline const void LocalToGlobalBaseMap::GlobalToLocal(
                const NekVector<NekDouble>& global,
                      NekVector<      NekDouble>& loc) const
        {
            v_GlobalToLocal(global,loc);
        }

        inline const void LocalToGlobalBaseMap::Assemble(
                const Array<OneD, const NekDouble> &loc,
                      Array<OneD,       NekDouble> &global) const
        {
            v_Assemble(loc,global);
        }

        inline const void LocalToGlobalBaseMap::Assemble(
                const NekVector<NekDouble>& loc,
                      NekVector<      NekDouble>& global) const
        {
            v_Assemble(loc,global);
        }

        inline const void LocalToGlobalBaseMap::UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            v_UniversalAssemble(pGlobal);
        }

        inline const void LocalToGlobalBaseMap::UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            v_UniversalAssemble(pGlobal);
        }

        inline const int LocalToGlobalBaseMap::GetFullSystemBandWidth() const
        {
            return v_GetFullSystemBandWidth();
        }

        inline int LocalToGlobalBaseMap::GetLocalToGlobalBndMap(const int i)
                                                                        const
        {
            return m_localToGlobalBndMap[i];
        }


        inline const Array<OneD,const int>&
                    LocalToGlobalBaseMap::GetLocalToGlobalBndMap(void)
        {
            return m_localToGlobalBndMap;
        }

        inline bool LocalToGlobalBaseMap::GetSignChange()
        {
            return m_signChange;
        }


        inline Array<OneD, const NekDouble>
                    LocalToGlobalBaseMap::GetLocalToGlobalBndSign(void) const
        {
            return m_localToGlobalBndSign;
        }

        inline const Array<OneD, const int>& LocalToGlobalBaseMap::GetGlobalToUniversalBndMap()
        {
            return m_globalToUniversalBndMap;
        }

        inline const Array<OneD, const int>& LocalToGlobalBaseMap::GetGlobalToUniversalBndMapUnique()
        {
            return m_globalToUniversalBndMapUnique;
        }

        inline NekDouble LocalToGlobalBaseMap::GetLocalToGlobalBndSign(
                    const int i) const
        {
            if(m_signChange)
            {
                return m_localToGlobalBndSign[i];
            }
            else
            {
                return 1.0;
            }
        }


        inline int LocalToGlobalBaseMap::GetBndCondCoeffsToGlobalCoeffsMap(
                    const int i)
        {
            return m_bndCondCoeffsToGlobalCoeffsMap[i];
        }


        inline int LocalToGlobalBaseMap::GetBndCondTraceToGlobalTraceMap(
                    const int i)
        {
            return m_bndCondTraceToGlobalTraceMap[i];
        }


        inline NekDouble
                    LocalToGlobalBaseMap::GetBndCondCoeffsToGlobalCoeffsSign(
                    const int i)
        {
            if(m_signChange)
            {
                return m_bndCondCoeffsToGlobalCoeffsSign[i];
            }
            else
            {
                return 1.0;
            }
        }

        inline const Array<OneD,const int>&
                    LocalToGlobalBaseMap::GetBndCondCoeffsToGlobalCoeffsMap()
        {
            return m_bndCondCoeffsToGlobalCoeffsMap;
        }


        inline int LocalToGlobalBaseMap::GetNumGlobalDirBndCoeffs() const
        {
            return m_numGlobalDirBndCoeffs;
        }


        inline int LocalToGlobalBaseMap::GetNumLocalDirBndCoeffs() const
        {
            return m_numLocalDirBndCoeffs;
        }

        inline int LocalToGlobalBaseMap::GetNumLocalBndCoeffs() const
        {
            return m_numLocalBndCoeffs;
        }

        inline int LocalToGlobalBaseMap::GetNumGlobalBndCoeffs() const
        {
            return m_numGlobalBndCoeffs;
        }

        inline int LocalToGlobalBaseMap::GetNumLocalCoeffs() const
        {
            return m_numLocalCoeffs;
        }

        inline int LocalToGlobalBaseMap::GetNumGlobalCoeffs() const
        {
            return m_numGlobalCoeffs;
        }

        inline bool LocalToGlobalBaseMap::GetSingularSystem() const
        {
            return m_systemSingular;
        }

        inline void LocalToGlobalBaseMap::GlobalToLocalBnd(
                    const NekVector<NekDouble>& global,
                    NekVector<NekDouble>& loc,
                    int offset) const
        {
            GlobalToLocalBnd(global.GetPtr(), loc.GetPtr(), offset);
        }


        inline void LocalToGlobalBaseMap::GlobalToLocalBnd(
                    const NekVector<NekDouble>& global,
                    NekVector<NekDouble>& loc) const
        {
            GlobalToLocalBnd(global.GetPtr(), loc.GetPtr());
        }


        inline void LocalToGlobalBaseMap::GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc, int offset) const
        {
            ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
            ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs-offset,"Global vector is not of correct dimension");

            // offset input data by length "offset" for Dirichlet boundary conditions.
            Array<OneD,NekDouble> tmp(m_numGlobalBndCoeffs,0.0);
            Vmath::Vcopy(m_numGlobalBndCoeffs-offset, global.get(), 1, tmp.get() + offset, 1);

            if(m_signChange)
            {
                Vmath::Gathr(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), tmp.get(), m_localToGlobalBndMap.get(), loc.get());
            }
            else
            {
                Vmath::Gathr(m_numLocalBndCoeffs, tmp.get(), m_localToGlobalBndMap.get(), loc.get());
            }
        }


        inline void LocalToGlobalBaseMap::GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc) const
        {
            ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
            ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

            if(m_signChange)
            {
                Vmath::Gathr(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), global.get(), m_localToGlobalBndMap.get(), loc.get());
            }
            else
            {
                Vmath::Gathr(m_numLocalBndCoeffs, global.get(), m_localToGlobalBndMap.get(), loc.get());
            }
        }


        inline void LocalToGlobalBaseMap::AssembleBnd(
                    const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global, int offset) const
        {
            AssembleBnd(loc.GetPtr(), global.GetPtr(), offset);
        }


        inline void LocalToGlobalBaseMap::AssembleBnd(
                    const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global) const
        {
            AssembleBnd(loc.GetPtr(), global.GetPtr());
        }


        inline void LocalToGlobalBaseMap::AssembleBnd(
                    const Array<OneD,const NekDouble>& loc,
                    Array<OneD, NekDouble>& global, int offset) const
        {
            ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local array is not of correct dimension");
            ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs-offset,"Global array is not of correct dimension");
            Array<OneD,NekDouble> tmp(m_numGlobalBndCoeffs,0.0);

            if(m_signChange)
            {
                Vmath::Assmb(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(),loc.get(), m_localToGlobalBndMap.get(), tmp.get());
            }
            else
            {
                Vmath::Assmb(m_numLocalBndCoeffs,loc.get(), m_localToGlobalBndMap.get(), tmp.get());
            }
            UniversalAssembleBnd(tmp);
            Vmath::Vcopy(m_numGlobalBndCoeffs-offset, tmp.get() + offset, 1, global.get(), 1);
        }


        inline void LocalToGlobalBaseMap::AssembleBnd(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD, NekDouble>& global) const
        {
            ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
            ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

            Vmath::Zero(m_numGlobalBndCoeffs, global.get(), 1);

            if(m_signChange)
            {
                Vmath::Assmb(m_numLocalBndCoeffs,m_localToGlobalBndSign.get(),
                             loc.get(), m_localToGlobalBndMap.get(), global.get());
            }
            else
            {
                Vmath::Assmb(m_numLocalBndCoeffs,loc.get(), m_localToGlobalBndMap.get(), global.get());
            }
            UniversalAssembleBnd(global);
        }

        inline const void LocalToGlobalBaseMap::UniversalAssembleBnd(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            ASSERTL1(pGlobal.num_elements() == m_numGlobalBndCoeffs,
                    "Wrong size.");
            Gs::Gather(pGlobal, Gs::gs_add, m_bndGsh);
        }

        inline const void LocalToGlobalBaseMap::UniversalAssembleBnd(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssembleBnd(pGlobal.GetPtr());
        }

        inline int LocalToGlobalBaseMap::GetBndSystemBandWidth() const
        {
            return m_bndSystemBandWidth;
        }

        inline int LocalToGlobalBaseMap::GetStaticCondLevel() const
        {
            return m_staticCondLevel;
        }
        
        inline int LocalToGlobalBaseMap::GetNumPatches() const
        {
            return m_numPatches;
        }

        inline const Array<OneD,const unsigned int>&
                    LocalToGlobalBaseMap::GetNumLocalBndCoeffsPerPatch()
        {
            return m_numLocalBndCoeffsPerPatch;
        }


        inline const Array<OneD,const unsigned int>&
            LocalToGlobalBaseMap::GetNumLocalIntCoeffsPerPatch()
        {
            return m_numLocalIntCoeffsPerPatch;
        }

        inline const LocalToGlobalBaseMapSharedPtr
                    LocalToGlobalBaseMap::GetNextLevelLocalToGlobalMap() const
        {
            return  m_nextLevelLocalToGlobalMap;
        }

        inline const PatchMapSharedPtr&
                    LocalToGlobalBaseMap::GetPatchMapFromPrevLevel(void)
                                                                        const
        {
            return m_patchMapFromPrevLevel;
        }

        inline bool LocalToGlobalBaseMap::AtLastLevel() const
        {
            return !( (bool) m_nextLevelLocalToGlobalMap.get() );
        }


        inline const GlobalSysSolnType
                    LocalToGlobalBaseMap::GetGlobalSysSolnType() const
        {
            return m_solnType;
        }

        inline void LocalToGlobalBaseMap::GlobalToLocalBndWithoutSign(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc)
        {
            ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
            ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

            Vmath::Gathr(m_numLocalBndCoeffs, global.get(), m_localToGlobalBndMap.get(), loc.get());
        }
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP_H


/**
 $Log: LocalToGlobalBaseMap.h,v $
 Revision 1.16  2009/10/30 14:02:55  pvos
 Multi-level static condensation updates

 Revision 1.15  2009/09/06 22:28:45  sherwin
 Updates for Navier-Stokes solver

 Revision 1.14  2009/04/27 11:42:14  sherwin
 Updated to make DG helmsolve efficient

 Revision 1.13  2009/04/08 06:38:55  sherwin
 Put eigensolve into NekMatrix. Some bandwidths mods

 Revision 1.12  2009/04/02 13:06:42  sherwin
 Modified to take symmetric banded system for HDH solver

 Revision 1.11  2009/03/04 14:17:38  pvos
 Removed all methods that take and Expansion as argument

 Revision 1.10  2009/02/08 09:10:47  sherwin
 Added NUllLocalToGlobalBaseMapSharedPtr definition

 Revision 1.9  2008/12/19 15:12:50  pvos
 Updates for precomputed dirichlet forcing functionality

 Revision 1.8  2008/12/17 17:08:22  pvos
 Performance updates

 Revision 1.7  2008/12/16 14:08:43  pvos
 Performance updates

 Revision 1.6  2008/11/01 22:36:06  bnelson
 Removed uneeded files.

 Revision 1.5  2008/11/01 22:07:46  bnelson
 Fixed compiler warning

 Revision 1.4  2008/09/23 18:21:00  pvos
 Updates for working ProjectContField3D demo

 Revision 1.3  2008/09/17 13:46:40  pvos
 Added LocalToGlobalC0ContMap for 3D expansions

 Revision 1.2  2008/09/16 13:36:06  pvos
 Restructured the LocalToGlobalMap classes

 Revision 1.1  2008/08/18 08:16:23  sherwin
 First version of this new class container for mappings

 */

