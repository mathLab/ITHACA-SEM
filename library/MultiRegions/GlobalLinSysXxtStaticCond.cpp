///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterativeStaticCond.cpp
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
// Description: GlobalLinSysIterativeStaticCond definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/Xxt.hpp>
#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/GlobalLinSysXxtStaticCond.h>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterativeStaticCond
         *
         * Solves a linear system iteratively using single- or multi-level
         * static condensation.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysXxtStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "XxtStaticCond",
                    GlobalLinSysXxtStaticCond::create,
                    "Iterative static condensation.");

        string GlobalLinSysXxtStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "XxtMultiLevelStaticCond",
                    GlobalLinSysXxtStaticCond::create,
                    "Xxt multi-level static condensation.");

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
        GlobalLinSysXxtStaticCond::GlobalLinSysXxtStaticCond(
                     const GlobalLinSysKey &pKey,
                     const boost::weak_ptr<ExpList> &pExpList,
                     const boost::shared_ptr<AssemblyMap>
                                                            &pLocToGloMap)
                : GlobalLinSysXxt(pKey, pExpList, pLocToGloMap),
                  m_locToGloMap (pLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==eXxtStaticCond)||
                     (pKey.GetGlobalSysSolnType()==eXxtMultiLevelStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(pKey.GetGlobalSysSolnType()
                        == pLocToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");

            // Allocate memory for top-level structure
            SetupTopLevel(pLocToGloMap);

            // Construct this level
            Initialise(pLocToGloMap);
        }


        /**
         *
         */
        GlobalLinSysXxtStaticCond::GlobalLinSysXxtStaticCond(
                     const GlobalLinSysKey &pKey,
                     const boost::weak_ptr<ExpList> &pExpList,
                     const DNekScalBlkMatSharedPtr pSchurCompl,
                     const DNekScalBlkMatSharedPtr pBinvD,
                     const DNekScalBlkMatSharedPtr pC,
                     const DNekScalBlkMatSharedPtr pInvD,
                     const boost::shared_ptr<AssemblyMap>
                                                            &pLocToGloMap)
                : GlobalLinSysXxt(pKey, pExpList, pLocToGloMap),
                  m_schurCompl ( pSchurCompl ),
                  m_BinvD      ( pBinvD ),
                  m_C          ( pC ),
                  m_invD       ( pInvD ),
                  m_locToGloMap( pLocToGloMap )
        {
            // Construct this level
            Initialise(pLocToGloMap);
        }


        /**
         *
         */
        GlobalLinSysXxtStaticCond::~GlobalLinSysXxtStaticCond()
        {

        }


        /**
         *
         */
        void GlobalLinSysXxtStaticCond::v_Solve(
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

            Array<OneD, NekDouble> F = m_wsp + nLocBndDofs;
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
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,m_wsp,eWrapper);

            NekVector<NekDouble> V_GlobHomBndTmp(nGlobHomBndDofs,0.0);

            if(nGlobHomBndDofs)
            {
                // construct boundary forcing
                if( nIntDofs  && ((!dirForcCalculated) && (atLastLevel)) )
                {
                    //include dirichlet boundary forcing
                    DNekScalBlkMat &BinvD      = *m_BinvD;
                    DNekScalBlkMat &SchurCompl = *m_schurCompl;
                    pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    V_LocBnd = BinvD*F_Int + SchurCompl*V_LocBnd;
                }
                else if((!dirForcCalculated) && (atLastLevel))
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

                // solve boundary system
                if(atLastLevel)
                {
                    Array<OneD, NekDouble> tmp(nGlobBndDofs, 0.0);
                    SolveLinearSystem(pLocToGloMap->GetNumLocalBndCoeffs(),
                                F, tmp, pLocToGloMap, nDirBndDofs);

                    // Enforce the Dirichlet boundary conditions on the solution
                    // array.
                    Vmath::Vadd(nGlobBndDofs, out, 1,
                                           tmp,    1,
                                           out, 1);
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
        void GlobalLinSysXxtStaticCond::Initialise(
                const boost::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            if(pLocToGloMap->AtLastLevel())
            {
                int nLocalBnd = m_locToGloMap->GetNumLocalBndCoeffs();
                int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
                m_wsp = Array<OneD, NekDouble>(nLocalBnd + nGlobal);

                CreateMap(pLocToGloMap);
                AssembleSchurComplementMatrixArrays(pLocToGloMap);
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
        void GlobalLinSysXxtStaticCond::SetupTopLevel(
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
                    DNekScalMatSharedPtr tmp_mat;
                    m_schurCompl->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,0));
                    m_BinvD     ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,1));
                    m_C         ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,0));
                    m_invD      ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,1));
                }
            }
        }


        /**
         * Create the inverse multiplicity map.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysXxtStaticCond::CreateMap(
                    const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
        {
            const Array<OneD, const int> &vMap
                                    = pLocToGloMap->GetLocalToGlobalBndMap();
            unsigned int nGlo       = pLocToGloMap->GetNumGlobalBndCoeffs();
            unsigned int nEntries   = pLocToGloMap->GetNumLocalBndCoeffs();
            unsigned int i;

            // Count the multiplicity of each global DOF on this process
            Array<OneD, NekDouble> vCounts(nGlo, 0.0);
            for (i = 0; i < nEntries; ++i)
            {
                vCounts[vMap[i]] += 1.0;
            }

            // Get universal multiplicity by globally assembling counts
            pLocToGloMap->UniversalAssembleBnd(vCounts);

            // Construct a map of 1/multiplicity for use in XXT solve
            m_locToGloSignMult = Array<OneD, NekDouble>(nEntries);
            for (i = 0; i < nEntries; ++i)
            {
                m_locToGloSignMult[i] = 1.0/vCounts[vMap[i]];
            }

            m_map = pLocToGloMap->GetLocalToGlobalBndMap();
        }

        /**
         * Construct the local matrix row index, column index and value index
         * arrays and initialize the XXT data structure with this information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysXxtStaticCond::AssembleSchurComplementMatrixArrays(
                        const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
        {
            ExpListSharedPtr vExp = m_expList.lock();
            unsigned int nElmt = m_schurCompl->GetNumberOfBlockRows();
            DNekScalMatSharedPtr loc_mat;
            unsigned int iCount     = 0;
            unsigned int rCount     = 0;
            unsigned int nRows      = 0;
            unsigned int nEntries   = 0;
            unsigned int numDirBnd  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int nLocal     = pLocToGloMap->GetNumLocalBndCoeffs();
            const Array<OneD, NekDouble> &vMapSign
                                    = pLocToGloMap->GetLocalToGlobalBndSign();
            bool doSign = pLocToGloMap->GetSignChange();
            unsigned int i = 0, j = 0, k = 0, n = 0;
            int gid1;
            Array<OneD, unsigned int> vSizes(nElmt);

            // First construct a map of the number of local DOFs in each block
            // and the number of matrix entries for each block
            for (n = 0; n < nElmt; ++n)
            {
                loc_mat = m_schurCompl->GetBlock(n,n);
                vSizes[n] = loc_mat->GetRows();
                nEntries += vSizes[n]*vSizes[n];
            }

            // Set up i-index, j-index and value arrays
            m_Ai = Array<OneD, unsigned int>(nEntries);
            m_Aj = Array<OneD, unsigned int>(nEntries);
            m_Ar = Array<OneD, double>(nEntries, 0.0);

            // Set up the universal ID array for XXT
            Array<OneD, unsigned long> vId(nLocal);

            // Loop over each elemental block, extract matrix indices and value
            // and set the universal ID array
            for(n = iCount = 0; n < nElmt; ++n)
            {
                loc_mat = m_schurCompl->GetBlock(n,n);
                nRows = loc_mat->GetRows();

                for(i = 0; i < nRows; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalBndMap(iCount + i);
                    for(j = 0; j < nRows; ++j)
                    {
                            k = rCount + i*vSizes[n] + j;
                            m_Ai[k] = iCount + i;
                            m_Aj[k] = iCount + j;
                            m_Ar[k] = (*loc_mat)(i,j);
                            if (doSign)
                            {
                                m_Ar[k] *= vMapSign[iCount+i]*vMapSign[iCount+j];
                            }
                    }

                    // Dirichlet DOFs are not included in the solve, so we set
                    // these to the special XXT id=0.
                    if (gid1 < numDirBnd)
                    {
                        vId[iCount + i] = 0;
                    }
                    else
                    {
                        vId[iCount + i]
                            = pLocToGloMap->GetGlobalToUniversalBndMap()[gid1];
                    }
                }
                iCount += vSizes[n];
                rCount += vSizes[n]*vSizes[n];
            }

            // Set up XXT and output some stats
            LibUtilities::CommSharedPtr vComm = pLocToGloMap->GetComm();
            m_crsData = Xxt::Init(nLocal, vId, m_Ai, m_Aj, m_Ar, vComm);
            Xxt::nektar_crs_stats(m_crsData);
        }



        /**
         *
         */
        void GlobalLinSysXxtStaticCond::ConstructNextLevelCondensedSystem(
                        const AssemblyMapSharedPtr& pLocToGloMap)
        {
            int i,j,n,cnt;
            NekDouble one  = 1.0;
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
#if 0 
                Array<OneD, const bool>      isBndDof;
#else
                Array<OneD, const unsigned int>      isBndDof;
#endif
                Array<OneD, const NekDouble> sign;
                NekDouble scale;

                for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
                {
                    schurComplSubMat      = SchurCompl->GetBlock(n,n);
                    schurComplSubMatnRows = schurComplSubMat->GetRows();
                    
                    scale = SchurCompl->GetBlock(n,n)->Scale();
                    Array<OneD, NekDouble> schurSubMat = SchurCompl->GetBlock(n,n)->GetOwnedMatrix()->GetPtr();
                    
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
                        (*substructuredMat[1][i]) = (*substructuredMat[1][i])*(*substructuredMat[3][i]);
                        // 3. A -> A - BInvD*C (= schurcomplement)
                        // (*substructuredMat[0][i]) = (*substructuredMat[0][i]) -
                        // (*substructuredMat[1][i])*(*substructuredMat[2][i]);
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
                        tmpscalmat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,substructuredMat[j][i]);
                        blkMatrices[j]->SetBlock(i,i,tmpscalmat);
                    }
                }
            }

            // We've finished with the Schur complement matrix passed to this
            // level, so return the memory to the system.
            // The Schur complement matrix need only be retained at the last
            // level. Save the other matrices at this level though.
            m_schurCompl.reset();

            m_recursiveSchurCompl = MemoryManager<GlobalLinSysXxtStaticCond>::
                AllocateSharedPtr(m_linSysKey,m_expList,blkMatrices[0],blkMatrices[1],blkMatrices[2],blkMatrices[3],pLocToGloMap);
        }


        /**
         *
         */
        void GlobalLinSysXxtStaticCond::v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            int nLocal = m_locToGloMap->GetNumLocalBndCoeffs();
            int nDir = m_locToGloMap->GetNumGlobalDirBndCoeffs();

            if (m_globalSchurCompl)
            {
                // Do matrix multiply globally

                Array<OneD, NekDouble> in  = pInput + nDir;
                Array<OneD, NekDouble> out = pOutput+ nDir;

                m_globalSchurCompl->Multiply(in,out);
            }
            else
            {
                // Do matrix multiply locally

                NekVector<NekDouble> loc(nLocal, m_wsp, eWrapper);
                m_locToGloMap->GlobalToLocalBnd(pInput, m_wsp);
                loc = (*m_schurCompl)*loc;
                m_locToGloMap->AssembleBnd(m_wsp, pOutput);
            }
        }


        /**
         * Diagonal preconditioner computed by evaluating the local matrix
         * acting on each basis vector (0,...,0,1,0,...,0). (deprecated)
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysXxtStaticCond::v_ComputePreconditioner()
        {
        /*
            ASSERTL1(m_gmat.get(),
                     "Matrix must be defined to compute preconditioner.");
            ASSERTL1(!m_preconditioner.get(),
                     "Preconditioner has already been defined.");

            int nGlobalBnd = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int n = m_gmat->GetRows();
            m_map = m_locToGloMap->GetGlobalToUniversalBndMapUnique();
            MatrixStorage storage = eDIAGONAL;
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(n, n, storage);
            DNekMat &M = (*m_preconditioner);

            // Extract diagonal contributions
            Array<OneD, NekDouble> vOutput(nGlobalBnd,0.0);
            for (unsigned int i = 0; i < n; ++i)
            {
                vOutput[nDirBnd + i] = (*m_gmat)(i,i);
            }

            // Assemble diagonal contributions across processes
            m_locToGloMap->UniversalAssembleBnd(vOutput);

            // Populate preconditioner matrix
            for (unsigned int i = 0; i < n; ++i)
            {
                M.SetValue(i,i,1.0/vOutput[nDirBnd + i]);
            }
        */
        }

        void GlobalLinSysXxtStaticCond::v_UniqueMap()
        {
            //m_map = m_locToGloMap->GetGlobalToUniversalBndMapUnique();
        }

    }
}
