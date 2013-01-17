///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMap.cpp
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

#include <MultiRegions/AssemblyMap/AssemblyMap.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class AssemblyMap
         * This class acts as a base class for constructing mappings between
         * local, global and boundary degrees of freedom. It holds the storage
         * for the maps and provides the accessors needed to retrieve them.
         *
         * There are two derived classes: AssemblyMapCG and
         * AssemblyMapDG. These perform the actual construction of the
         * maps within their specific contexts.
         *
         */

        /// Rounds a double precision number to an integer.
        int RoundNekDoubleToInt(NekDouble x)
        {
            return int(x > 0.0 ? x + 0.5 : x - 0.5);
        }

        /// Rounds an array of double precision numbers to integers.
        void RoundNekDoubleToInt(const Array<OneD,const NekDouble> inarray, Array<OneD,int> outarray)
        {
            int size = inarray.num_elements();
            ASSERTL1(outarray.num_elements()>=size,"Array sizes not compatible");

            NekDouble x;
            for(int i = 0; i < size; i++)
            {
                x = inarray[i];
                outarray[i] =  int(x > 0.0 ? x + 0.5 : x - 0.5);
            }
        }

        /**
         * Initialises an empty mapping.
         */
        AssemblyMap::AssemblyMap():
            m_session(),
            m_comm(),
            m_hash(0),
            m_solnType(eNoSolnType),
            m_numLocalBndCoeffs(0),
            m_numGlobalBndCoeffs(0),
            m_numLocalDirBndCoeffs(0),
            m_numGlobalDirBndCoeffs(0),
            m_bndSystemBandWidth(0),
            m_gsh(0),
            m_bndGsh(0)
        {
        }

        AssemblyMap::AssemblyMap(const LibUtilities::SessionReaderSharedPtr &pSession):
            m_session(pSession),
            m_comm(pSession->GetComm()),
            m_hash(0),
            m_solnType(pSession->GetSolverInfoAsEnum<GlobalSysSolnType>("GlobalSysSoln")),
            m_preconType(pSession->GetSolverInfoAsEnum<PreconditionerType>("Preconditioner")),
            m_numLocalBndCoeffs(0),
            m_numGlobalBndCoeffs(0),
            m_numLocalDirBndCoeffs(0),
            m_numGlobalDirBndCoeffs(0),
            m_bndSystemBandWidth(0),
            m_gsh(0),
            m_bndGsh(0)
        {
        }


        /** 
         * Create a new level of mapping using the information in
         * multiLevelGraph and performing the following steps:
         */
        AssemblyMap::AssemblyMap(
                    AssemblyMap* oldLevelMap,
                    const BottomUpSubStructuredGraphSharedPtr& multiLevelGraph):
            m_session(oldLevelMap->m_session),
            m_comm(oldLevelMap->GetComm()),
            m_hash(0),
            m_solnType(oldLevelMap->m_solnType),
            m_globalToUniversalBndMap(oldLevelMap->GetGlobalToUniversalBndMap()),
            m_globalToUniversalBndMapUnique(oldLevelMap->GetGlobalToUniversalBndMapUnique()),
            m_gsh(oldLevelMap->m_gsh),
            m_bndGsh(oldLevelMap->m_bndGsh),
            m_preconType(oldLevelMap->m_preconType),
            m_lowestStaticCondLevel(oldLevelMap->m_lowestStaticCondLevel)
        {
            int i;
            int j;
            int cnt;

            //--------------------------------------------------------------
            // -- Extract information from the input argument
            int numGlobalBndCoeffsOld    = oldLevelMap->GetNumGlobalBndCoeffs();
            int numGlobalDirBndCoeffsOld = oldLevelMap->GetNumGlobalDirBndCoeffs();
            int numGlobalHomBndCoeffsOld = numGlobalBndCoeffsOld - numGlobalDirBndCoeffsOld;
            int numLocalBndCoeffsOld     = oldLevelMap->GetNumLocalBndCoeffs();
            int numLocalDirBndCoeffsOld  = oldLevelMap->GetNumLocalDirBndCoeffs();
            bool signChangeOld           = oldLevelMap->GetSignChange();

            int staticCondLevelOld       = oldLevelMap->GetStaticCondLevel();
            int numPatchesOld            = oldLevelMap->GetNumPatches();
            GlobalSysSolnType solnTypeOld = oldLevelMap->GetGlobalSysSolnType();
            const Array<OneD, const unsigned int>& numLocalBndCoeffsPerPatchOld = oldLevelMap->GetNumLocalBndCoeffsPerPatch();
            //--------------------------------------------------------------

            //--------------------------------------------------------------
            int newLevel = staticCondLevelOld+1;
            /** - STEP 1: setup a mask array to determine to which patch
             *          of the new level every patch of the current
             *          level belongs.  To do so we make four arrays,
             *          #gloPatchMask, #globHomPatchMask,
             *          #locPatchMask_NekDouble and #locPatchMask.
             *          These arrays are then used to check which local
             *          dofs of the old level belong to which patch of
             *          the new level
             */
            Array<OneD, NekDouble> globPatchMask         (numGlobalBndCoeffsOld,-1.0);
            Array<OneD, NekDouble> globHomPatchMask      (globPatchMask+numGlobalDirBndCoeffsOld);
            Array<OneD, NekDouble> locPatchMask_NekDouble(numLocalBndCoeffsOld,-3.0);
            Array<OneD, int>       locPatchMask          (numLocalBndCoeffsOld);

            // Fill the array globPatchMask as follows:
            // - The first part (i.e. the glob bnd dofs) is filled with the
            //   value -1
            // - The second part (i.e. the glob interior dofs) is numbered
            //   according to the patch it belongs to (i.e. dofs in first block
            //   all are numbered 0, the second block numbered are 1, etc...)
            multiLevelGraph->MaskPatches(newLevel,globHomPatchMask);
            
            // Map from Global Dofs to Local Dofs
            // As a result, we know for each local dof whether
            // it is mapped to the boundary of the next level, or to which
            // patch. Based upon this, we can than later associate every patch
            // of the current level with a patch in the next level.
            oldLevelMap->GlobalToLocalBndWithoutSign(globPatchMask,locPatchMask_NekDouble);

            // Convert the result to an array of integers rather than doubles
            RoundNekDoubleToInt(locPatchMask_NekDouble,locPatchMask);

            /** - STEP 2: We calculate how many local bnd dofs of the
             *  old level belong to the boundaries of each patch at
             *  the new level. We need this information to set up the
             *  mapping between different levels.
             */

            // Retrieve the number of patches at the next level
            int numPatchesWithIntNew = multiLevelGraph->GetNpatchesWithInterior(newLevel);
            int numPatchesNew        = numPatchesWithIntNew;
            
            // Allocate memory to store the number of local dofs associated to
            // each of elemental boundaries of these patches
            map<int, int> numLocalBndCoeffsPerPatchNew;
            for(int i = 0; i < numPatchesNew; i++)
            {
                numLocalBndCoeffsPerPatchNew[i] = 0;
            }

            int minval;
            int maxval;
            int curPatch;
            for(i = cnt = 0; i < numPatchesOld; i++)
            {
                // For every patch at the current level, the mask array
                // locPatchMask should be filled with
                // - the same (positive) number for each entry (which will
                //   correspond to the patch at the next level it belongs to)
                // - the same (positive) number for each entry, except some
                //   entries that are -1 (the enties correspond to -1, will be
                //   mapped to the local boundary of the next level patch given
                //   by the positive number)
                // - -1 for all entries. In this case, we will make an
                //   additional patch only consisting of boundaries at the next
                //   level
                minval = *min_element(&locPatchMask[cnt],
                                      &locPatchMask[cnt]+numLocalBndCoeffsPerPatchOld[i]);
                maxval = *max_element(&locPatchMask[cnt],
                                      &locPatchMask[cnt]+numLocalBndCoeffsPerPatchOld[i]);
                ASSERTL0((minval==maxval)||(minval==-1),"These values should never be the same");
                
                if(maxval == -1)
                {
                    curPatch = numPatchesNew;
                    numLocalBndCoeffsPerPatchNew[curPatch] = 0;
                    numPatchesNew++;
                }
                else
                {
                    curPatch = maxval;
                }
                
                for(j = 0; j < numLocalBndCoeffsPerPatchOld[i]; j++ )
                {
                    ASSERTL0((locPatchMask[cnt]==maxval)||(locPatchMask[cnt]==minval),
                             "These values should never be the same");
                    if(locPatchMask[cnt] == -1)
                    {
                        ++numLocalBndCoeffsPerPatchNew[curPatch];
                    }
                    cnt++;
                }
            }

            // Count how many local dofs of the old level are mapped
            // to the local boundary dofs of the new level
            m_numLocalBndCoeffs  = 0;
            m_numPatches         = numLocalBndCoeffsPerPatchNew.size();
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches,0u);
            for(int i = 0; i < m_numPatches; i++)
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) numLocalBndCoeffsPerPatchNew[i];
                m_numLocalBndCoeffs           += numLocalBndCoeffsPerPatchNew[i];
            }
            multiLevelGraph->GetNintDofsPerPatch(newLevel,m_numLocalIntCoeffsPerPatch);

            // Also initialise some more data members
            m_solnType              = solnTypeOld;
            ASSERTL1(m_solnType==eDirectMultiLevelStaticCond
                    ||m_solnType==eIterativeMultiLevelStaticCond
                    ||m_solnType==eXxtMultiLevelStaticCond,
                     "This method should only be called for in "
                     "case of multi-level static condensation.");
            m_staticCondLevel       = newLevel;
            m_signChange            = signChangeOld;
            m_numLocalDirBndCoeffs  = numLocalDirBndCoeffsOld;
            m_numGlobalDirBndCoeffs = numGlobalDirBndCoeffsOld;
            m_numGlobalBndCoeffs    = multiLevelGraph->GetInteriorOffset(newLevel) +
                m_numGlobalDirBndCoeffs;
            m_numGlobalCoeffs       = multiLevelGraph->GetNumGlobalDofs(newLevel) +
                m_numGlobalDirBndCoeffs;
            m_localToGlobalBndMap   = Array<OneD,int>(m_numLocalBndCoeffs);
            if(m_signChange)
            {
                m_localToGlobalBndSign    = Array<OneD,NekDouble>(m_numLocalBndCoeffs);
            }
            
            m_patchMapFromPrevLevel = MemoryManager<PatchMap>::AllocateSharedPtr(numLocalBndCoeffsOld);

            // Set up an offset array that denotes the offset of the local
            // boundary degrees of freedom of the next level
            Array<OneD, int> numLocalBndCoeffsPerPatchOffset(m_numPatches+1,0);
            for(int i = 1; i < m_numPatches+1; i++)
            {
                numLocalBndCoeffsPerPatchOffset[i] += numLocalBndCoeffsPerPatchOffset[i-1] + numLocalBndCoeffsPerPatchNew[i-1];
            }

            int additionalPatchCnt = numPatchesWithIntNew;
            int newid;
            int blockid;
            bool isBndDof;
            NekDouble sign;
            Array<OneD, int> bndDofPerPatchCnt(m_numPatches,0);
            for(i = cnt = 0; i < numPatchesOld; i++)
            {
                minval = *min_element(&locPatchMask[cnt],&locPatchMask[cnt]+numLocalBndCoeffsPerPatchOld[i]);
                maxval = *max_element(&locPatchMask[cnt],&locPatchMask[cnt]+numLocalBndCoeffsPerPatchOld[i]);
                ASSERTL0((minval==maxval)||(minval==-1),"These values should never be the same");

                if(maxval == -1)
                {
                    curPatch = additionalPatchCnt;
                    additionalPatchCnt++;
                }
                else
                {
                    curPatch = maxval;
                }

                for(j = 0; j < numLocalBndCoeffsPerPatchOld[i]; j++ )
                {
                    ASSERTL0((locPatchMask[cnt]==maxval)||(locPatchMask[cnt]==minval),
                             "These values should never be the same");

                    sign = oldLevelMap->GetLocalToGlobalBndSign(cnt);

                    if(locPatchMask[cnt] == -1)
                    {
                        newid = numLocalBndCoeffsPerPatchOffset[curPatch];

                        m_localToGlobalBndMap[newid] = oldLevelMap->GetLocalToGlobalBndMap(cnt);
                        if(m_signChange)
                        {
                            m_localToGlobalBndSign[ newid ] = sign;
                        }

                        blockid = bndDofPerPatchCnt[curPatch];
                        isBndDof = true;


                        numLocalBndCoeffsPerPatchOffset[curPatch]++;
                        bndDofPerPatchCnt[curPatch]++;
                    }
                    else
                    {
                        newid = oldLevelMap->GetLocalToGlobalBndMap(cnt) -
                            m_numGlobalBndCoeffs+m_numLocalBndCoeffs;

                        blockid = oldLevelMap->GetLocalToGlobalBndMap(cnt)-
                            m_numGlobalDirBndCoeffs - multiLevelGraph->GetInteriorOffset(newLevel,curPatch);
                        isBndDof = false;
                    }

                    sign = isBndDof?1.0:sign;

                    m_patchMapFromPrevLevel->SetPatchMap(cnt,curPatch,blockid,isBndDof,sign);

                    cnt++;
                }
            }
            CalculateBndSystemBandWidth();

            // Postprocess the computed information - Update the old
            // level with the mapping to new evel
            // oldLevelMap->SetLocalBndToNextLevelMap(oldLocalBndToNewLevelMap,oldLocalBndToNewLevelSign);
            // - Construct the next level mapping object
            if(m_staticCondLevel < (multiLevelGraph->GetNlevels()-1))
            {
                m_nextLevelLocalToGlobalMap = MemoryManager<AssemblyMap>::AllocateSharedPtr(this,multiLevelGraph);
            }
        }


        AssemblyMap::~AssemblyMap(void)
        {
        }


        /**
         * The bandwidth calculated corresponds to what is referred to as
         * half-bandwidth.  If the elements of the matrix are designated as
         * a_ij, it corresponds to the maximum value of |i-j| for non-zero
         * a_ij.  As a result, the value also corresponds to the number of
         * sub- or super-diagonals.
         *
         * The bandwith can be calculated elementally as it corresponds to the
         * maximal elemental bandwith (i.e. the maximal difference in global
         * DOF index for every element).
         *
         * We here calculate the bandwith of the global boundary system (as
         * used for static condensation).
         */
        void AssemblyMap::CalculateBndSystemBandWidth()
        {
            int i,j;
            int cnt = 0;
            int locSize;
            int maxId;
            int minId;
            int bwidth = -1;
            for(i = 0; i < m_numPatches; ++i)
            {
                locSize = m_numLocalBndCoeffsPerPatch[i];
                maxId = -1;
                minId = m_numLocalBndCoeffs+1;
                for(j = 0; j < locSize; j++)
                {
                    if(m_localToGlobalBndMap[cnt+j] >= m_numGlobalDirBndCoeffs)
                    {
                        if(m_localToGlobalBndMap[cnt+j] > maxId)
                        {
                            maxId = m_localToGlobalBndMap[cnt+j];
                        }

                        if(m_localToGlobalBndMap[cnt+j] < minId)
                        {
                            minId = m_localToGlobalBndMap[cnt+j];
                        }
                    }
                }
                bwidth = (bwidth>(maxId-minId))?bwidth:(maxId-minId);

                cnt+=locSize;
            }

            m_bndSystemBandWidth = bwidth;
        }


        int AssemblyMap::v_GetLocalToGlobalMap(const int i) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            return 0;
        }

        int AssemblyMap::v_GetGlobalToUniversalMap(const int i) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            return 0;
        }

        int AssemblyMap::v_GetGlobalToUniversalMapUnique(const int i) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            return 0;
        }

        const Array<OneD,const int>&  AssemblyMap::v_GetLocalToGlobalMap()
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            static Array<OneD,const int> result;
            return result;
        }

        const Array<OneD, const int>& AssemblyMap::v_GetGlobalToUniversalMap()
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            static Array<OneD, const int> result;
            return result;
        }

        const Array<OneD, const int>& AssemblyMap::v_GetGlobalToUniversalMapUnique()
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            static Array<OneD, const int> result;
            return result;
        }

        NekDouble AssemblyMap::v_GetLocalToGlobalSign(const int i) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            return 0.0;
        }

        const Array<OneD, NekDouble>& AssemblyMap::v_GetLocalToGlobalSign() const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            static Array<OneD, NekDouble> result;
            return result;
        }

        const void AssemblyMap::v_LocalToGlobal(
                const Array<OneD, const NekDouble>& loc,
                      Array<OneD,       NekDouble>& global) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
        }

        const void AssemblyMap::v_LocalToGlobal(
                const NekVector<NekDouble>& loc,
                      NekVector<      NekDouble>& global) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
        }

        const void AssemblyMap::v_GlobalToLocal(
                const Array<OneD, const NekDouble>& global,
                      Array<OneD,       NekDouble>& loc) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
        }

        const void AssemblyMap::v_GlobalToLocal(
                const NekVector<NekDouble>& global,
                      NekVector<      NekDouble>& loc) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
        }

        const void AssemblyMap::v_Assemble(
                const Array<OneD, const NekDouble> &loc,
                      Array<OneD,       NekDouble> &global) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
        }

        const void AssemblyMap::v_Assemble(
                const NekVector<NekDouble>& loc,
                      NekVector<      NekDouble>& global) const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
        }

        const void AssemblyMap::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            // Do nothing here since multi-level static condensation uses a
            // AssemblyMap and thus will call this routine in serial.
        }

        const void AssemblyMap::v_UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            // Do nothing here since multi-level static condensation uses a
            // AssemblyMap and thus will call this routine in serial.
        }

        const void AssemblyMap::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal,
                      int                         offset) const
        {
            // Do nothing here since multi-level static condensation uses a
            // AssemblyMap and thus will call this routine in serial.
        }

        const int AssemblyMap::v_GetFullSystemBandWidth() const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            return 0;
        }

        int AssemblyMap::v_GetNumNonDirVertexModes() const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            return 0;
        }

        int AssemblyMap::v_GetNumNonDirEdgeModes() const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            return 0;
        }

        int AssemblyMap::v_GetNumNonDirFaceModes() const
        {
            ASSERTL0(false, "Not defined for this type of mapping.");
            return 0;
        }



        LibUtilities::CommSharedPtr AssemblyMap::GetComm()
        {
            return m_comm;
        }

        size_t AssemblyMap::GetHash() const
        {
            return m_hash;
        }

        int AssemblyMap::GetLocalToGlobalMap(const int i) const
        {
            return v_GetLocalToGlobalMap(i);
        }

        int AssemblyMap::GetGlobalToUniversalMap(const int i) const
        {
            return v_GetGlobalToUniversalMap(i);
        }

        int AssemblyMap::GetGlobalToUniversalMapUnique(const int i) const
        {
            return v_GetGlobalToUniversalMapUnique(i);
        }

        const Array<OneD,const int>&  AssemblyMap::GetLocalToGlobalMap()
        {
            return v_GetLocalToGlobalMap();
        }

        const Array<OneD, const int>& AssemblyMap::GetGlobalToUniversalMap()
        {
            return v_GetGlobalToUniversalMap();
        }

        const Array<OneD, const int>& AssemblyMap::GetGlobalToUniversalMapUnique()
        {
            return v_GetGlobalToUniversalMapUnique();
        }

        NekDouble AssemblyMap::GetLocalToGlobalSign(const int i) const
        {
            return v_GetLocalToGlobalSign(i);
        }

        const Array<OneD, NekDouble>& AssemblyMap::GetLocalToGlobalSign() const
        {
            return v_GetLocalToGlobalSign();
        }

        const void AssemblyMap::LocalToGlobal(
                const Array<OneD, const NekDouble>& loc,
                      Array<OneD,       NekDouble>& global) const
        {
            v_LocalToGlobal(loc,global);
        }

        const void AssemblyMap::LocalToGlobal(
                const NekVector<NekDouble>& loc,
                      NekVector<      NekDouble>& global) const
        {
            v_LocalToGlobal(loc,global);
        }

        const void AssemblyMap::GlobalToLocal(
                const Array<OneD, const NekDouble>& global,
                      Array<OneD,       NekDouble>& loc) const
        {
            v_GlobalToLocal(global,loc);
        }

        const void AssemblyMap::GlobalToLocal(
                const NekVector<NekDouble>& global,
                      NekVector<      NekDouble>& loc) const
        {
            v_GlobalToLocal(global,loc);
        }

        const void AssemblyMap::Assemble(
                const Array<OneD, const NekDouble> &loc,
                      Array<OneD,       NekDouble> &global) const
        {
            v_Assemble(loc,global);
        }

        const void AssemblyMap::Assemble(
                const NekVector<NekDouble>& loc,
                      NekVector<      NekDouble>& global) const
        {
            v_Assemble(loc,global);
        }

        const void AssemblyMap::UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            v_UniversalAssemble(pGlobal);
        }

        const void AssemblyMap::UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            v_UniversalAssemble(pGlobal);
        }

        const void AssemblyMap::UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal,
                      int                         offset) const
        {
            v_UniversalAssemble(pGlobal, offset);
        }

        const int AssemblyMap::GetFullSystemBandWidth() const
        {
            return v_GetFullSystemBandWidth();
        }

        int AssemblyMap::GetNumNonDirVertexModes() const
        {
            return v_GetNumNonDirVertexModes();
        }

        int AssemblyMap::GetNumNonDirEdgeModes() const
        {
            return v_GetNumNonDirEdgeModes();
        }

        int AssemblyMap::GetNumNonDirFaceModes() const
        {
            return v_GetNumNonDirFaceModes();
        }

        int AssemblyMap::GetLocalToGlobalBndMap(const int i) const
        {
            return m_localToGlobalBndMap[i];
        }


        const Array<OneD,const int>&
                    AssemblyMap::GetLocalToGlobalBndMap(void)
        {
            return m_localToGlobalBndMap;
        }

        bool AssemblyMap::GetSignChange()
        {
            return m_signChange;
        }


        Array<OneD, const NekDouble>
                    AssemblyMap::GetLocalToGlobalBndSign(void) const
        {
            return m_localToGlobalBndSign;
        }

        const Array<OneD, const int>& AssemblyMap::GetGlobalToUniversalBndMap()
        {
            return m_globalToUniversalBndMap;
        }

        const Array<OneD, const int>& AssemblyMap::GetGlobalToUniversalBndMapUnique()
        {
            return m_globalToUniversalBndMapUnique;
        }

        NekDouble AssemblyMap::GetLocalToGlobalBndSign(const int i) const
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


        int AssemblyMap::GetBndCondCoeffsToGlobalCoeffsMap(
                    const int i)
        {
            return m_bndCondCoeffsToGlobalCoeffsMap[i];
        }


        int AssemblyMap::GetBndCondTraceToGlobalTraceMap(
                    const int i)
        {
            return m_bndCondTraceToGlobalTraceMap[i];
        }


        NekDouble AssemblyMap::GetBndCondCoeffsToGlobalCoeffsSign(const int i)
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

        const Array<OneD,const int>&
                    AssemblyMap::GetBndCondCoeffsToGlobalCoeffsMap()
        {
            return m_bndCondCoeffsToGlobalCoeffsMap;
        }


        int AssemblyMap::GetNumGlobalDirBndCoeffs() const
        {
            return m_numGlobalDirBndCoeffs;
        }


        int AssemblyMap::GetNumLocalDirBndCoeffs() const
        {
            return m_numLocalDirBndCoeffs;
        }

        int AssemblyMap::GetNumLocalBndCoeffs() const
        {
            return m_numLocalBndCoeffs;
        }

        int AssemblyMap::GetNumGlobalBndCoeffs() const
        {
            return m_numGlobalBndCoeffs;
        }

        int AssemblyMap::GetNumLocalCoeffs() const
        {
            return m_numLocalCoeffs;
        }

        int AssemblyMap::GetNumGlobalCoeffs() const
        {
            return m_numGlobalCoeffs;
        }

        bool AssemblyMap::GetSingularSystem() const
        {
            return m_systemSingular;
        }

        void AssemblyMap::GlobalToLocalBnd(
                    const NekVector<NekDouble>& global,
                    NekVector<NekDouble>& loc,
                    int offset) const
        {
            GlobalToLocalBnd(global.GetPtr(), loc.GetPtr(), offset);
        }


        void AssemblyMap::GlobalToLocalBnd(
                    const NekVector<NekDouble>& global,
                    NekVector<NekDouble>& loc) const
        {
            GlobalToLocalBnd(global.GetPtr(), loc.GetPtr());
        }


        void AssemblyMap::GlobalToLocalBnd(
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


        void AssemblyMap::GlobalToLocalBnd(
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

        void AssemblyMap::LocalBndToGlobal(
                    const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global,
                    int offset) const
        {
            LocalBndToGlobal(loc.GetPtr(), global.GetPtr(), offset);
        }


        void AssemblyMap::LocalBndToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,NekDouble>& global,
                    int offset) const
        {
            ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
            ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs-offset,"Global vector is not of correct dimension");

            // offset input data by length "offset" for Dirichlet boundary conditions.
            Array<OneD,NekDouble> tmp(m_numGlobalBndCoeffs,0.0);

            if(m_signChange)
            {
                Vmath::Scatr(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), loc.get(), m_localToGlobalBndMap.get(), tmp.get());
            }
            else
            {
                Vmath::Scatr(m_numLocalBndCoeffs, loc.get(), m_localToGlobalBndMap.get(), tmp.get());
            }
            Vmath::Vcopy(m_numGlobalBndCoeffs-offset, tmp.get()+offset, 1, global.get(), 1);
        }

        void AssemblyMap::AssembleBnd(
                    const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global, int offset) const
        {
            AssembleBnd(loc.GetPtr(), global.GetPtr(), offset);
        }


        void AssemblyMap::AssembleBnd(
                    const NekVector<NekDouble>& loc,
                    NekVector<NekDouble>& global) const
        {
            AssembleBnd(loc.GetPtr(), global.GetPtr());
        }


        void AssemblyMap::AssembleBnd(
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


        void AssemblyMap::AssembleBnd(
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

        const void AssemblyMap::UniversalAssembleBnd(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            ASSERTL1(pGlobal.num_elements() == m_numGlobalBndCoeffs,
                     "Wrong size.");
            Gs::Gather(pGlobal, Gs::gs_add, m_bndGsh);
        }

        const void AssemblyMap::UniversalAssembleBnd(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssembleBnd(pGlobal.GetPtr());
        }

        const void AssemblyMap::UniversalAssembleBnd(
                      Array<OneD,     NekDouble>& pGlobal,
                      int                         offset) const
        {
            Array<OneD, NekDouble> tmp(offset);
            if (offset > 0)  Vmath::Vcopy(offset, pGlobal, 1, tmp, 1);
            UniversalAssembleBnd(pGlobal);
            if (offset > 0)  Vmath::Vcopy(offset, tmp, 1, pGlobal, 1);
        }

        int AssemblyMap::GetBndSystemBandWidth() const
        {
            return m_bndSystemBandWidth;
        }

        int AssemblyMap::GetStaticCondLevel() const
        {
            return m_staticCondLevel;
        }
        
        int AssemblyMap::GetNumPatches() const
        {
            return m_numPatches;
        }

        const Array<OneD,const unsigned int>&
                    AssemblyMap::GetNumLocalBndCoeffsPerPatch()
        {
            return m_numLocalBndCoeffsPerPatch;
        }


        const Array<OneD,const unsigned int>&
            AssemblyMap::GetNumLocalIntCoeffsPerPatch()
        {
            return m_numLocalIntCoeffsPerPatch;
        }

        const AssemblyMapSharedPtr
                    AssemblyMap::GetNextLevelLocalToGlobalMap() const
        {
            return  m_nextLevelLocalToGlobalMap;
        }

        const PatchMapSharedPtr&
                    AssemblyMap::GetPatchMapFromPrevLevel(void)
                                                                        const
        {
            return m_patchMapFromPrevLevel;
        }

        bool AssemblyMap::AtLastLevel() const
        {
            return !( m_nextLevelLocalToGlobalMap.get() );
        }


        const GlobalSysSolnType
                    AssemblyMap::GetGlobalSysSolnType() const
        {
            return m_solnType;
        }

        const PreconditionerType
                    AssemblyMap::GetPreconType() const
        {
            return m_preconType;
        }

        void AssemblyMap::GlobalToLocalBndWithoutSign(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc)
        {
            ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
            ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

            Vmath::Gathr(m_numLocalBndCoeffs, global.get(), m_localToGlobalBndMap.get(), loc.get());
        }
    } // namespace
} // namespace
