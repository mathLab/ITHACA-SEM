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

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysDirect
         *
         * Solves a linear system using single- or multi-level static
         * condensation.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysDirectStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "DirectStaticCond",
                    GlobalLinSysDirectStaticCond::create,
                    "Direct static condensation.");

        string GlobalLinSysDirectStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "DirectMultiLevelStaticCond",
                    GlobalLinSysDirectStaticCond::create,
                    "Direct multi-level static condensation.");

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
        GlobalLinSysDirectStaticCond::GlobalLinSysDirectStaticCond(
                     const GlobalLinSysKey &pKey,
                     const boost::weak_ptr<ExpList> &pExpList,
                     const boost::shared_ptr<AssemblyMap>
                     &pLocToGloMap)
                : GlobalLinSysDirect(pKey, pExpList, pLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==eDirectStaticCond)||
                     (pKey.GetGlobalSysSolnType()==eDirectMultiLevelStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(pKey.GetGlobalSysSolnType()
                        == pLocToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");

            // Allocate memory for top-level structure
            SetupTopLevel(pLocToGloMap);

            // Construct this level
            Initialise(pLocToGloMap,DetermineMatrixStorage(pLocToGloMap));
        }


//        GlobalLinSysDirectStaticCond::GlobalLinSysDirectStaticCond(
//                     const DNekScalBlkMatSharedPtr pSchurCompl,
//                     const DNekScalBlkMatSharedPtr pBinvD,
//                     const DNekScalBlkMatSharedPtr pC,
//                     const DNekScalBlkMatSharedPtr pInvD,
//                     const boost::shared_ptr<AssemblyMap>
//                     &pLocToGloMap)
//        {
//            m_schurCompl = pSchurCompl;
//            m_BinvD      = pBinvD;
//            m_C          = pC;
//            m_invD       = pInvD;
//
//            // Construct this level
//            Initialise(pLocToGloMap, eFULL);
//        }
        

        /**
         *
         */
        GlobalLinSysDirectStaticCond::GlobalLinSysDirectStaticCond(
                     const GlobalLinSysKey &pKey,
                     const boost::weak_ptr<ExpList> &pExpList,
                     const DNekScalBlkMatSharedPtr pSchurCompl,
                     const DNekScalBlkMatSharedPtr pBinvD,
                     const DNekScalBlkMatSharedPtr pC,
                     const DNekScalBlkMatSharedPtr pInvD,
                     const boost::shared_ptr<AssemblyMap>
                                                            &pLocToGloMap)
                : GlobalLinSysDirect(pKey, pExpList, pLocToGloMap),
                  m_schurCompl ( pSchurCompl ),
                  m_BinvD      ( pBinvD ),
                  m_C          ( pC ),
                  m_invD       ( pInvD )
        {
            // Construct this level
            Initialise(pLocToGloMap,DetermineMatrixStorage(pLocToGloMap));
        }


        /**
         *
         */
        GlobalLinSysDirectStaticCond::~GlobalLinSysDirectStaticCond()
        {

        }


        /**
         *
         */
        void GlobalLinSysDirectStaticCond::v_Solve(
                    const Array<OneD, const NekDouble>  &in,
                          Array<OneD,       NekDouble>  &out,
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &dirForcing)
        {
            bool dirForcCalculated = (bool) dirForcing.num_elements();
            bool atLastLevel       = pLocToGloMap->AtLastLevel();

            int nGlobDofs          = pLocToGloMap->GetNumGlobalCoeffs();
            int nGlobBndDofs       = pLocToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = pLocToGloMap->GetNumLocalBndCoeffs();
            int nIntDofs           = pLocToGloMap->GetNumGlobalCoeffs()
                                                                - nGlobBndDofs;

            Array<OneD, NekDouble> F(nGlobDofs);
            Array<OneD, NekDouble> tmp;

            if(nDirBndDofs && dirForcCalculated)
            {
                Vmath::Vsub(nGlobDofs,in.get(),1,dirForcing.get(),1,F.get(),1);
            }
            else
            {
                Vmath::Vcopy(nGlobDofs,in.get(),1,F.get(),1);
            }

            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,tmp=F+nDirBndDofs,
                                          eWrapper);
            NekVector<NekDouble> F_Int(nIntDofs,tmp=F+nGlobBndDofs,eWrapper);

            NekVector<NekDouble> V_GlobBnd(nGlobBndDofs,out,eWrapper);
            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,
                                              tmp=out+nDirBndDofs,
                                              eWrapper);
            NekVector<NekDouble> V_Int(nIntDofs,tmp=out+nGlobBndDofs,eWrapper);
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
                        DNekScalBlkMat &BinvD      = *m_BinvD;
                        DNekScalBlkMat &SchurCompl = *m_schurCompl;
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);

                        V_LocBnd = BinvD*F_Int + SchurCompl*V_LocBnd;

                    }
                    else if((nDirBndDofs) && (!dirForcCalculated)
                                          && (atLastLevel))
                    {
                        //include dirichlet boundary forcing
                        DNekScalBlkMat &SchurCompl = *m_schurCompl;
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                        V_LocBnd = SchurCompl*V_LocBnd;
                    }
                    else
                    {
                        DNekScalBlkMat &BinvD      = *m_BinvD;
                        V_LocBnd = BinvD*F_Int;
                    }
                    pLocToGloMap->AssembleBnd(V_LocBnd,V_GlobHomBndTmp,
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
                                pLocToGloMap->GetNextLevelLocalToGlobalMap());
                }
            }

            // solve interior system
            if(nIntDofs)
            {
                DNekScalBlkMat &invD  = *m_invD;

                if(nGlobHomBndDofs || nDirBndDofs)
                {
                    DNekScalBlkMat &C     = *m_C;

                    if(dirForcCalculated && nDirBndDofs)
                    {
                        pLocToGloMap->GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd,
                                                      nDirBndDofs);
                    }
                    else
                    {
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    }
                    F_Int = F_Int - C*V_LocBnd;
                }

                V_Int = invD*F_Int;
            }
        }



        /**
         * If at the last level of recursion (or the only level in the case of
         * single-level static condensation), assemble the Schur complement.
         * For other levels, in the case of multi-level static condensation,
         * the next level of the condensed system is computed.
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysDirectStaticCond::Initialise(
           const boost::shared_ptr<AssemblyMap>& pLocToGloMap,
           const MatrixStorage matStorage)
        {
            if(pLocToGloMap->AtLastLevel())
            {
                AssembleSchurComplement(pLocToGloMap,matStorage);
            }
            else
            {
                ConstructNextLevelCondensedSystem(
                        pLocToGloMap->GetNextLevelLocalToGlobalMap());
            }
        }


        /**
         * For the first level in multi-level static condensation, or the only
         * level in the case of single-level static condensation, allocate the
         * condensed matrices and populate them with the local matrices
         * retrieved from the expansion list.
         * @param
         */
        void GlobalLinSysDirectStaticCond::SetupTopLevel(
                const boost::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            int n;
            int n_exp = m_expList.lock()->GetNumElmts();

            const Array<OneD,const unsigned int>& nbdry_size
                    = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
            const Array<OneD,const unsigned int>& nint_size
                    = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

            // Setup Block Matrix systems
            MatrixStorage blkmatStorage = eDIAGONAL;
            m_schurCompl = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
            m_BinvD      = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
            m_C          = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
            m_invD       = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

            for(n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() == StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr loc_mat = GetBlock(m_expList.lock()->GetOffset_Elmt_Id(n));
                    m_schurCompl->SetBlock(n,n,loc_mat);
                }
                else
                {
                    DNekScalBlkMatSharedPtr loc_mat = GetStaticCondBlock(m_expList.lock()->GetOffset_Elmt_Id(n));
                    DNekScalMatSharedPtr    tmp_mat;
                    m_schurCompl->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,0));
                    m_BinvD     ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,1));
                    m_C         ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,0));
                    m_invD      ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,1));
                }
            }
        }

        MatrixStorage GlobalLinSysDirectStaticCond::DetermineMatrixStorage(const AssemblyMapSharedPtr &pLocToGloMap)
        {
            int nBndDofs  = pLocToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int rows = nBndDofs - NumDirBCs;
            int bwidth = pLocToGloMap->GetBndSystemBandWidth();

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
                        matStorage = ePOSITIVE_DEFINITE_SYMMETRIC_BANDED; 
                    }
                    else
                    {
                        matStorage = ePOSITIVE_DEFINITE_SYMMETRIC;
                    }
                }
                break;
            case StdRegions::eLinearAdvectionReaction:
            case StdRegions::eLinearAdvectionDiffusionReaction:
            default:
                {
                    // Current inversion techniques do not seem to
                    // allow banded matrices to be used as a linear
                    // system
                    matStorage = eFULL;
                }
                break;
            }

            return matStorage;
        }

        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysDirectStaticCond::AssembleSchurComplement(
                      const AssemblyMapSharedPtr &pLocToGloMap,
                      MatrixStorage matStorage )
        {
            int i,j,n,cnt,gid1,gid2;
            NekDouble sign1,sign2,value;

            int nBndDofs  = pLocToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            DNekScalBlkMatSharedPtr SchurCompl = m_schurCompl;
            DNekScalBlkMatSharedPtr BinvD      = m_BinvD;
            DNekScalBlkMatSharedPtr C          = m_C;
            DNekScalBlkMatSharedPtr invD       = m_invD;

            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;
            NekDouble zero = 0.0;

            DNekMatSharedPtr Gmat;
            int bwidth = pLocToGloMap->GetBndSystemBandWidth();
         
            switch(matStorage)
            {
            case ePOSITIVE_DEFINITE_SYMMETRIC_BANDED:
                {
                    try {
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
                break;
            case ePOSITIVE_DEFINITE_SYMMETRIC:
            case eFULL:
                {
                    Gmat = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(rows, cols, zero, matStorage);
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal, "Unknown matrix storage type of type not set up");
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
                    gid1  = pLocToGloMap->GetLocalToGlobalBndMap (cnt + i)
                                                                    - NumDirBCs;
                    sign1 = pLocToGloMap->GetLocalToGlobalBndSign(cnt + i);

                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2  = pLocToGloMap->GetLocalToGlobalBndMap(cnt+j)
                                                                 - NumDirBCs;
                            sign2 = pLocToGloMap->GetLocalToGlobalBndSign(cnt+j);

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
                        const AssemblyMapSharedPtr& pLocToGloMap)
        {
            int i,j,n,cnt;
            DNekScalBlkMatSharedPtr blkMatrices[4];

            
            // Create temporary matrices within an inner-local scope to ensure
            // any references to the intermediate storage is lost before
            // the recursive step, rather than at the end of the routine.
            // This allows the schur complement matrix from this level to be
            // disposed of in the next level after use without being retained
            // due to lingering shared pointers.
            {

                const Array<OneD,const unsigned int>& nBndDofsPerPatch
                                = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nIntDofsPerPatch
                                = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

                // STEP 1:
                // Based upon the schur complement of the the current level we
                // will substructure this matrix in the form
                //      --     --
                //      | A   B |
                //      | C   D |
                //      --     --
                // All matrices A,B,C and D are (diagonal) blockmatrices.
                // However, as a start we will use an array of DNekMatrices as
                // it is too hard to change the individual entries of a
                // DNekScalBlkMatSharedPtr.

                // In addition, we will also try to ensure that the memory of
                // the blockmatrices will be contiguous. This will probably
                // enhance the efficiency
                // - Calculate the total number of entries in the blockmatrices
                int nPatches  = pLocToGloMap->GetNumPatches();
                int nEntriesA = 0; int nEntriesB = 0;
                int nEntriesC = 0; int nEntriesD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    nEntriesA += nBndDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesB += nBndDofsPerPatch[i]*nIntDofsPerPatch[i];
                    nEntriesC += nIntDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesD += nIntDofsPerPatch[i]*nIntDofsPerPatch[i];
                }

                // Now create the DNekMatrices and link them to the memory
                // allocated above
                Array<OneD, DNekMatSharedPtr> substructuredMat[4]
                    = {Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix A
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
                    substructuredMat[0][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                        nBndDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix B
                    tmparray = storageB+cntB;
                    substructuredMat[1][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                        nIntDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix C
                    tmparray = storageC+cntC;
                    substructuredMat[2][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                        nBndDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix D
                    tmparray = storageD+cntD;
                    substructuredMat[3][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                        nIntDofsPerPatch[i],
                                                        tmparray, wType);

                    cntA += nBndDofsPerPatch[i] * nBndDofsPerPatch[i];
                    cntB += nBndDofsPerPatch[i] * nIntDofsPerPatch[i];
                    cntC += nIntDofsPerPatch[i] * nBndDofsPerPatch[i];
                    cntD += nIntDofsPerPatch[i] * nIntDofsPerPatch[i];
                }

                // Then, project SchurComplement onto
                // the substructured matrices of the next level
                DNekScalBlkMatSharedPtr SchurCompl  = m_schurCompl;
                DNekScalMatSharedPtr schurComplSubMat;
                int       schurComplSubMatnRows;
                Array<OneD, const int>       patchId, dofId;
                Array<OneD, const unsigned int>      isBndDof;
                Array<OneD, const NekDouble> sign;
                NekDouble scale;

                for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
                {
                    schurComplSubMat      = SchurCompl->GetBlock(n,n);
                    schurComplSubMatnRows = schurComplSubMat->GetRows();
                    
                    scale = SchurCompl->GetBlock(n,n)->Scale();
                    const Array<OneD, const NekDouble> schurSubMat = SchurCompl->GetBlock(n,n)->GetOwnedMatrix()->GetPtr();
                    

                    patchId  = pLocToGloMap->GetPatchMapFromPrevLevel()->GetPatchId()+ cnt;
                    dofId    = pLocToGloMap->GetPatchMapFromPrevLevel()->GetDofId()+ cnt;
                    isBndDof = pLocToGloMap->GetPatchMapFromPrevLevel()->IsBndDof() + cnt;
                    sign     = pLocToGloMap->GetPatchMapFromPrevLevel()->GetSign() + cnt;

                    // Set up  Matrix;
                    for(i = 0; i < schurComplSubMatnRows; ++i)
                    {
                        Array<OneD, NekDouble> subMat0 = substructuredMat[0][patchId[i]]->GetPtr();
                        int subMat0rows = substructuredMat[0][patchId[i]]->GetRows();
                        Array<OneD, NekDouble> subMat1 = substructuredMat[1][patchId[i]]->GetPtr();
                        int subMat1rows = substructuredMat[1][patchId[i]]->GetRows();
                        Array<OneD, NekDouble> subMat2 = substructuredMat[2][patchId[i]]->GetPtr();
                        int subMat2rows = substructuredMat[2][patchId[i]]->GetRows();
                        Array<OneD, NekDouble> subMat3 = substructuredMat[3][patchId[i]]->GetPtr();
                        int subMat3rows = substructuredMat[3][patchId[i]]->GetRows();
                        
                        if(isBndDof[i])
                        {
                            for(j = 0; j < schurComplSubMatnRows; ++j)
                            {
                                ASSERTL0(patchId[i]==patchId[j],"These values should be equal");
                                
                                if(isBndDof[j])
                                {
                                    subMat0[dofId[i]+dofId[j]*subMat0rows] += 
                                        sign[i]*sign[j]*(scale*schurSubMat[i+j*schurComplSubMatnRows]);
                                }
                                else
                                {
                                    subMat1[dofId[i]+dofId[j]*subMat1rows] += 
                                        sign[i]*sign[j]*(scale*schurSubMat[i+j*schurComplSubMatnRows]);
                                }
                            }
                        }
                        else
                        {
                            for(j = 0; j < schurComplSubMatnRows; ++j)
                            {
                                ASSERTL0(patchId[i]==patchId[j],"These values should be equal");
                                
                                if(isBndDof[j])
                                {
                                    subMat2[dofId[i]+dofId[j]*subMat2rows] += 
                                        sign[i]*sign[j]*(scale*schurSubMat[i+j*schurComplSubMatnRows]);
                                }
                                else
                                {
                                    subMat3[dofId[i]+dofId[j]*subMat3rows] += 
                                        sign[i]*sign[j]*(scale*schurSubMat[i+j*schurComplSubMatnRows]);
                                }
                            }
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
                        Array<OneD, NekDouble> subMat0 = substructuredMat[0][i]->GetPtr();
                        int subMat0rows = substructuredMat[0][i]->GetRows();
                        Array<OneD, NekDouble> subMat1 = substructuredMat[1][i]->GetPtr();
                        int subMat1rows = substructuredMat[1][i]->GetRows();
                        Array<OneD, NekDouble> subMat2 = substructuredMat[2][i]->GetPtr();
                        int subMat2rows = substructuredMat[2][i]->GetRows();
                        int subMat2cols = substructuredMat[2][i]->GetColumns();

                        // 1. D -> InvD
                        substructuredMat[3][i]->Invert();
                        // 2. B -> BInvD
                        (*substructuredMat[1][i]) = (*substructuredMat[1][i])*
                            (*substructuredMat[3][i]);

                        // 3. A -> A - BInvD*C (= schurcomplement)
                        //(*substructuredMat[0][i]) = (*substructuredMat[0][i]) -
                        //(*substructuredMat[1][i])*(*substructuredMat[2][i]);
                        // Note: faster to use blas directly
                        Blas::Dgemm('N','N', subMat1rows, subMat2cols, subMat2rows, -1.0,
                                    &subMat1[0], subMat1rows, &subMat2[0], subMat2rows, 
                                    1.0, &subMat0[0], subMat0rows);
                    }
                }

                // STEP 3: fill the blockmatrices
                // however, do note that we first have to convert them to
                // a DNekScalMat in order to be compatible with the first
                // level of static condensation

                const Array<OneD,const unsigned int>& nbdry_size = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nint_size  = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();
                MatrixStorage blkmatStorage = eDIAGONAL;

                blkMatrices[0] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
                blkMatrices[1] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
                blkMatrices[2] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
                blkMatrices[3] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

                DNekScalMatSharedPtr tmpscalmat;
                for(i = 0; i < nPatches; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        tmpscalmat = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0,substructuredMat[j][i]);
                        blkMatrices[j]->SetBlock(i,i,tmpscalmat);
                    }
                }
            }

            // We've finished with the Schur complement matrix passed to this
            // level, so return the memory to the system.
            // The Schur complement matrix need only be retained at the last
            // level. Save the other matrices at this level though.
            m_schurCompl.reset();

            m_recursiveSchurCompl = MemoryManager<GlobalLinSysDirectStaticCond>::
                AllocateSharedPtr(m_linSysKey,m_expList,blkMatrices[0],blkMatrices[1],blkMatrices[2],blkMatrices[3],pLocToGloMap);
            
        }

    }
}
