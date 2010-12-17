///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.cpp
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
// Description: GlobalLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysDirectStaticCond.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>
#include <MultiRegions/LocalMatrixSystem.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysDirect
         *
         * Solves a linear system using a full matrix.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysDirectStaticCond::className
                = GlobalLinSysFactory::RegisterCreatorFunction(
                    "DirectStaticCond",
                    GlobalLinSysDirectStaticCond::create,
                    "Direct (possibly multi-level) static condensation.");

        /**
         * For a matrix system of the form @f[
         * \left[ \begin{array}{cc}
         * \boldsymbol{A} & \boldsymbol{B}\\
         * \boldsymbol{C} & \boldsymbol{D}
         * \end{array} \right]
         * \left[ \begin{array}{c} \boldsymbol{x_1}\\ \boldsymbol{x_2}
         * \end{array}\right]
         * = \left[ \begin{array}{c} \boldsymbol{y_1}\\ \boldsymbol{y_2}
         * \end{array}\right],
         * @f]
         * where @f$\boldsymbol{D}@f$ and
         * @f$(\boldsymbol{A-BD^{-1}C})@f$ are invertible, store and assemble
         * a static condensation system, according to a given local to global
         * mapping. #m_linSys is constructed by AssembleSchurComplement().
         * @param   mKey        Associated matrix key.
         * @param   pLocMatSys  LocalMatrixSystem
         * @param   locToGloMap Local to global mapping.
         */
        GlobalLinSysDirectStaticCond::GlobalLinSysDirectStaticCond(const GlobalLinSysKey &mkey,
                     const boost::shared_ptr<LocalMatrixSystem> &pLocMatSys,
                     const boost::shared_ptr<LocalToGlobalBaseMap>
                                                            &locToGloMap)
                : GlobalLinSysDirect(mkey, pLocMatSys, locToGloMap)
        {
            ASSERTL1((mkey.GetGlobalSysSolnType()==eDirectStaticCond)||
                     (mkey.GetGlobalSysSolnType()==eDirectMultiLevelStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(mkey.GetGlobalSysSolnType()
                        == locToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");

            m_locMatSys = pLocMatSys;

            if(locToGloMap->AtLastLevel())
            {
                AssembleSchurComplement(locToGloMap);
            }
            else
            {
                ConstructNextLevelCondensedSystem(
                        locToGloMap->GetNextLevelLocalToGlobalMap());
            }
        }

        GlobalLinSysDirectStaticCond::~GlobalLinSysDirectStaticCond()
        {

        }


        void GlobalLinSysDirectStaticCond::Solve( const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out,
                    const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                    const Array<OneD, const NekDouble> &dirForcing)
        {
            bool dirForcCalculated = (bool) dirForcing.num_elements();
            bool atLastLevel       = locToGloMap->AtLastLevel();

            int nGlobDofs          = locToGloMap->GetNumGlobalCoeffs();
            int nGlobBndDofs       = locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = locToGloMap->GetNumLocalBndCoeffs();
            int nIntDofs           = locToGloMap->GetNumGlobalCoeffs()
                                                                - nGlobBndDofs;

            Array<OneD, NekDouble> F(nGlobDofs);
            if(nDirBndDofs && dirForcCalculated)
            {
                Vmath::Vsub(nGlobDofs,in.get(),1,dirForcing.get(),1,F.get(),1);
            }
            else
            {
                Vmath::Vcopy(nGlobDofs,in.get(),1,F.get(),1);
            }

            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,F+nDirBndDofs,
                                          eWrapper);
            NekVector<NekDouble> F_Int(nIntDofs,F+nGlobBndDofs,eWrapper);

            NekVector<NekDouble> V_GlobBnd(nGlobBndDofs,out,eWrapper);
            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,out+nDirBndDofs,
                                              eWrapper);
            NekVector<NekDouble> V_Int(nIntDofs,out+nGlobBndDofs,eWrapper);
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,0.0);

            NekVector<NekDouble> V_GlobHomBndTmp(nGlobHomBndDofs,0.0);

            if(nGlobHomBndDofs)
            {
                if(nIntDofs || ((nDirBndDofs) && (!dirForcCalculated)
                                              && (atLastLevel)) )
                {
                    // construct boundary forcing
                    if( nIntDofs  && ((nDirBndDofs) && (!dirForcCalculated)
                                                    && (atLastLevel)) )
                    {
                        //include dirichlet boundary forcing
                        DNekScalBlkMat &BinvD      = *m_locMatSys->GetLocalSystem()[1];
                        DNekScalBlkMat &SchurCompl = *m_locMatSys->GetLocalSystem()[0];
                        locToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                        V_LocBnd = BinvD*F_Int + SchurCompl*V_LocBnd;
                    }
                    else if((nDirBndDofs) && (!dirForcCalculated)
                                          && (atLastLevel))
                    {
                        //include dirichlet boundary forcing
                        DNekScalBlkMat &SchurCompl = *m_locMatSys->GetLocalSystem()[0];
                        locToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                        V_LocBnd = SchurCompl*V_LocBnd;
                    }
                    else
                    {
                        DNekScalBlkMat &BinvD      = *m_locMatSys->GetLocalSystem()[1];
                        V_LocBnd = BinvD*F_Int;
                    }
                    locToGloMap->AssembleBnd(V_LocBnd,V_GlobHomBndTmp,
                                             nDirBndDofs);
                    F_HomBnd = F_HomBnd - V_GlobHomBndTmp;
                }

                // solve boundary system
                if(atLastLevel)
                {
                    m_linSys->Solve(F_HomBnd,V_GlobHomBnd);
                }
                else
                {
                    m_recursiveSchurCompl->Solve(F,
                                V_GlobBnd.GetPtr(),
                                locToGloMap->GetNextLevelLocalToGlobalMap());
                }
            }

            // solve interior system
            if(nIntDofs)
            {
                DNekScalBlkMat &invD  = *m_locMatSys->GetLocalSystem()[3];

                if(nGlobHomBndDofs || nDirBndDofs)
                {
                    DNekScalBlkMat &C     = *m_locMatSys->GetLocalSystem()[2];

                    if(dirForcCalculated && nDirBndDofs)
                    {
                        locToGloMap->GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd,
                                                      nDirBndDofs);
                    }
                    else
                    {
                        locToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    }
                    F_Int = F_Int - C*V_LocBnd;
                }

                V_Int = invD*F_Int;
            }
        }



        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysDirectStaticCond::AssembleSchurComplement(
                    const LocalToGlobalBaseMapSharedPtr &locToGloMap)
        {
            int i,j,n,cnt,gid1,gid2;
            NekDouble sign1,sign2,value;

            int nBndDofs  = locToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = locToGloMap->GetNumGlobalDirBndCoeffs();

            DNekScalBlkMatSharedPtr SchurCompl = m_locMatSys->GetLocalSystem()[0];
            DNekScalBlkMatSharedPtr BinvD      = m_locMatSys->GetLocalSystem()[1];
            DNekScalBlkMatSharedPtr C          = m_locMatSys->GetLocalSystem()[2];
            DNekScalBlkMatSharedPtr invD       = m_locMatSys->GetLocalSystem()[3];

            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;
            NekDouble zero = 0.0;

            DNekMatSharedPtr Gmat;
            int bwidth = locToGloMap->GetBndSystemBandWidth();
            MatrixStorage matStorage;

            switch(m_linSysKey.GetMatrixType())
            {
                // case for all symmetric matices
            case StdRegions::eMass:
            case StdRegions::eLaplacian:
            case StdRegions::eHelmholtz:
            case StdRegions::eHybridDGHelmBndLam:
                {
                    if( (2*(bwidth+1)) < rows)
                    {
                        try {
                            matStorage = ePOSITIVE_DEFINITE_SYMMETRIC_BANDED;
                            Gmat = MemoryManager<DNekMat>
                                ::AllocateSharedPtr(rows, cols, zero,
                                                    matStorage,
                                                    bwidth, bwidth);
                        }
                        catch (...) {
                            NEKERROR(ErrorUtil::efatal,
                                     "Insufficient memory for GlobalLinSys.");
                        }
                    }
                    else
                    {
                        matStorage = ePOSITIVE_DEFINITE_SYMMETRIC;
                        Gmat = MemoryManager<DNekMat>
                                        ::AllocateSharedPtr(rows, cols, zero,
                                                            matStorage);
                    }
                }
                break;
            case StdRegions::eLinearAdvectionReaction:
            case StdRegions::eLinearAdvectionDiffusionReaction:
                {
                    // Current inversion techniques do not seem to
                    // allow banded matrices to be used as a linear
                    // system
                    matStorage = eFULL;
                    Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols, zero,
                                                                     matStorage);

                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal, "Add MatrixType to switch "
                             "statement");
                }
            }

            // fill global matrix
            DNekScalMatSharedPtr loc_mat;
            int loc_lda;
            for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
            {
                loc_mat = SchurCompl->GetBlock(n,n);
                loc_lda = loc_mat->GetRows();

                // Set up  Matrix;
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1  = locToGloMap->GetLocalToGlobalBndMap (cnt + i)
                                                                    - NumDirBCs;
                    sign1 = locToGloMap->GetLocalToGlobalBndSign(cnt + i);

                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2  = locToGloMap->GetLocalToGlobalBndMap(cnt+j)
                                                                 - NumDirBCs;
                            sign2 = locToGloMap->GetLocalToGlobalBndSign(cnt+j);

                            if(gid2 >= 0)
                            {
                                // As the global matrix should be
                                // symmetric, only add the value for
                                // the upper triangular part in order
                                // to avoid entries to be entered
                                // twice
                                if((matStorage == eFULL)||(gid2 >= gid1))
                                {
                                    value = Gmat->GetValue(gid1,gid2)
                                                + sign1*sign2*(*loc_mat)(i,j);
                                    Gmat->SetValue(gid1,gid2,value);
                                }
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            if(rows)
            {
                PointerWrapper w = eWrapper;
                m_linSys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,w);
            }
        }


        /**
         *
         */
        void GlobalLinSysDirectStaticCond::ConstructNextLevelCondensedSystem(
                        const LocalToGlobalBaseMapSharedPtr& locToGloMap)
        {
            int i,j,n,cnt;
            NekDouble one  = 1.0;
            NekDouble zero = 0.0;
            LocalMatrixSystemSharedPtr pNextLevelLocMatSys = MemoryManager<LocalMatrixSystemStaticCond>::AllocateSharedPtr(locToGloMap);
            DNekScalBlkMatSharedPtr blkMatrices[4];

            // Create temporary matrices within an inner-local scope to ensure
            // any references to the intermediate storage is lost before
            // the recursive step, rather than at the end of the routine.
            // This allows the schur complement matrix from this level to be
            // disposed of in the next level after use without being retained
            // due to lingering shared pointers.
            {

                const Array<OneD,const unsigned int>& nBndDofsPerPatch = locToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nIntDofsPerPatch = locToGloMap->GetNumLocalIntCoeffsPerPatch();

                // STEP 1:
                // Based upon the schur complement of the the current level we will
                // substructure this matrix in the form
                //      --     --
                //      | A   B |
                //      | C   D |
                //      --     --
                // All matrices A,B,C and D are (diagonal) blockmatrices. However,
                // as a start we will use an array of DNekMatrices as it is too hard
                // to change the individual entries of a DNekScalBlkMatSharedPtr.

                // In addition, we will also try to ensure that the memory of the
                // blockmatrices will be contiguous. This will probably enhance
                // the efficiency
                // - Calculate the total number of entries in the blockmatrices
                int nPatches  = locToGloMap->GetNumPatches();
                int nEntriesA = 0; int nEntriesB = 0;
                int nEntriesC = 0; int nEntriesD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    nEntriesA += nBndDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesB += nBndDofsPerPatch[i]*nIntDofsPerPatch[i];
                    nEntriesC += nIntDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesD += nIntDofsPerPatch[i]*nIntDofsPerPatch[i];
                }

                // Now create the DNekMatrices and link them to the memory allocated above
                Array<OneD, DNekMatSharedPtr> substructeredMat[4] = {Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix A
                                                                     Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix B
                                                                     Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix C
                                                                     Array<OneD, DNekMatSharedPtr>(nPatches)}; //Matrix D

                // Initialise storage for the matrices. We do this separately
                // for each matrix so the matrices may be independently
                // deallocated when no longer required.
                Array<OneD, NekDouble> storageA(nEntriesA,0.0);
                Array<OneD, NekDouble> storageB(nEntriesB,0.0);
                Array<OneD, NekDouble> storageC(nEntriesC,0.0);
                Array<OneD, NekDouble> storageD(nEntriesD,0.0);

                Array<OneD, NekDouble> tmparray;
                PointerWrapper wType = eWrapper;
                int cntA = 0;
                int cntB = 0;
                int cntC = 0;
                int cntD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    // Matrix A
                    tmparray = storageA+cntA;
                    substructeredMat[0][i] = MemoryManager<DNekMat>::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                                                       nBndDofsPerPatch[i],
                                                                                       tmparray,wType);
                    // Matrix B
                    tmparray = storageB+cntB;
                    substructeredMat[1][i] = MemoryManager<DNekMat>::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                                                       nIntDofsPerPatch[i],
                                                                                       tmparray,wType);
                    // Matrix C
                    tmparray = storageC+cntC;
                    substructeredMat[2][i] = MemoryManager<DNekMat>::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                                                       nBndDofsPerPatch[i],
                                                                                       tmparray,wType);
                    // Matrix D
                    tmparray = storageD+cntD;
                    substructeredMat[3][i] = MemoryManager<DNekMat>::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                                                       nIntDofsPerPatch[i],
                                                                                       tmparray,wType);

                    cntA += nBndDofsPerPatch[i]*nBndDofsPerPatch[i];
                    cntB += nBndDofsPerPatch[i]*nIntDofsPerPatch[i];
                    cntC += nIntDofsPerPatch[i]*nBndDofsPerPatch[i];
                    cntD += nIntDofsPerPatch[i]*nIntDofsPerPatch[i];
                }

                // Then, project SchurComplement onto
                // the substructured matrices of the next level
                DNekScalBlkMatSharedPtr SchurCompl  = m_locMatSys->GetLocalSystem()[0];
                DNekScalMatSharedPtr schurComplSubMat;
                int       schurComplSubMatnRows;
                int       patchId_i ,patchId_j;
                int       dofId_i   ,dofId_j;
                bool      isBndDof_i,isBndDof_j;
                NekDouble sign_i    ,sign_j;
                NekDouble value;
                int       ABCorD;
                for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
                {
                    schurComplSubMat      = SchurCompl->GetBlock(n,n);
                    schurComplSubMatnRows = schurComplSubMat->GetRows();

                    // Set up  Matrix;
                    for(i = 0; i < schurComplSubMatnRows; ++i)
                    {
                        patchId_i  = locToGloMap->GetPatchMapFromPrevLevel(cnt+i)->GetPatchId();
                        dofId_i    = locToGloMap->GetPatchMapFromPrevLevel(cnt+i)->GetDofId();
                        isBndDof_i = locToGloMap->GetPatchMapFromPrevLevel(cnt+i)->IsBndDof();
                        sign_i     = locToGloMap->GetPatchMapFromPrevLevel(cnt+i)->GetSign();

                        for(j = 0; j < schurComplSubMatnRows; ++j)
                        {
                            patchId_j  = locToGloMap->GetPatchMapFromPrevLevel(cnt+j)->GetPatchId();
                            dofId_j    = locToGloMap->GetPatchMapFromPrevLevel(cnt+j)->GetDofId();
                            isBndDof_j = locToGloMap->GetPatchMapFromPrevLevel(cnt+j)->IsBndDof();
                            sign_j     = locToGloMap->GetPatchMapFromPrevLevel(cnt+j)->GetSign();

                            ASSERTL0(patchId_i==patchId_j,"These values should be equal");

                            ABCorD = 2*(isBndDof_i?0:1)+(isBndDof_j?0:1);
                            value = substructeredMat[ABCorD][patchId_i]->GetValue(dofId_i,dofId_j) +
                                sign_i*sign_j*(*schurComplSubMat)(i,j);

                            substructeredMat[ABCorD][patchId_i]->SetValue(dofId_i,dofId_j,value);
                        }
                    }
                    cnt += schurComplSubMatnRows;
                }

                // STEP 2: condense the system
                // This can be done elementally (i.e. patch per patch)
                for(i = 0; i < nPatches; i++)
                {
                    if(nIntDofsPerPatch[i])
                    {
                        // 1. D -> InvD
                        substructeredMat[3][i]->Invert();
                        // 2. B -> BInvD
                        (*substructeredMat[1][i]) = (*substructeredMat[1][i])*(*substructeredMat[3][i]);
                        // 3. A -> A - BInvD*C (= schurcomplement)
                        (*substructeredMat[0][i]) = (*substructeredMat[0][i]) -
                            (*substructeredMat[1][i])*(*substructeredMat[2][i]);
                    }
                }

                // STEP 3: fill the blockmatrices
                // however, do note that we first have to convert them to
                // a DNekScalMat in order to be compatible with the first
                // level of static condensation

                blkMatrices[0] = pNextLevelLocMatSys->GetLocalSystem()[0];
                blkMatrices[1] = pNextLevelLocMatSys->GetLocalSystem()[1];
                blkMatrices[2] = pNextLevelLocMatSys->GetLocalSystem()[2];
                blkMatrices[3] = pNextLevelLocMatSys->GetLocalSystem()[3];

                DNekScalMatSharedPtr tmpscalmat;
                for(i = 0; i < nPatches; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        tmpscalmat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,substructeredMat[j][i]);
                        blkMatrices[j]->SetBlock(i,i,tmpscalmat);
                    }
                }
            }

            // We've finished with the Schur complement matrix passed to this
            // level, so return the memory to the system.
            // The Schur complement matrix need only be retained at the last
            // level. Save the other matrices at this level though.
            m_locMatSys->DeallocateSystem(0);

            m_recursiveSchurCompl = MemoryManager<GlobalLinSysDirectStaticCond>::
                AllocateSharedPtr(m_linSysKey,pNextLevelLocMatSys,locToGloMap);
        }

    }
}
