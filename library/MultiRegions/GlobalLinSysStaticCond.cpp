///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysStaticCond.cpp
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
// Description: Implementation to linear solver using single-
//              or multi-level static condensation
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysStaticCond.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>
#include <LibUtilities/LinearAlgebra/SparseDiagBlkMatrix.hpp>
#include <LibUtilities/LinearAlgebra/SparseUtils.hpp>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysStaticCond
         *
         * Solves a linear system using single- or multi-level static
         * condensation.
         */

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
        GlobalLinSysStaticCond::GlobalLinSysStaticCond(
            const GlobalLinSysKey                &pKey,
            const std::weak_ptr<ExpList>         &pExpList,
            const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
                : GlobalLinSys(pKey, pExpList, pLocToGloMap),
                  m_locToGloMap (pLocToGloMap)
        {
        }

        void GlobalLinSysStaticCond::v_InitObject()
        {
            // Allocate memory for top-level structure
            SetupTopLevel(m_locToGloMap.lock());

            // Construct this level
            Initialise(m_locToGloMap.lock());
        }

        /**
         *
         */
        GlobalLinSysStaticCond::~GlobalLinSysStaticCond()
        {

        }


        /**
         *
         */
        void GlobalLinSysStaticCond::v_Solve(
            const Array<OneD, const NekDouble> &pLocInput,
                  Array<OneD,       NekDouble> &pLocOutput,
            const AssemblyMapSharedPtr         &pLocToGloMap,
            const Array<OneD, const NekDouble> &dirForcing)
        {
            boost::ignore_unused(dirForcing);
            ASSERTL1( dirForcing.size() == 0,
                      "GlobalLinSysStaticCond: Not setup for dirForcing");

            bool atLastLevel       = pLocToGloMap->AtLastLevel();
            int  scLevel           = pLocToGloMap->GetStaticCondLevel();
            
            int nGlobDofs    = pLocToGloMap->GetNumGlobalCoeffs();
            int nLocBndDofs  = pLocToGloMap->GetNumLocalBndCoeffs();
            int nGlobBndDofs = pLocToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nIntDofs     = nGlobDofs - nGlobBndDofs;

            if((nGlobDofs-nDirBndDofs) == 0)
            {
                return; //nothing to solve; 
            }

            Array<OneD, NekDouble> F_bnd, F_bnd1, F_int, V_bnd; 
            Array<OneD, NekDouble> tmp;

            F_bnd  = m_wsp;
            F_bnd1 = m_wsp + nLocBndDofs;
            V_bnd  = m_wsp + 2*nLocBndDofs;
            F_int  = m_wsp + 3*nLocBndDofs; 

            pLocToGloMap->LocalToLocalBnd(pLocOutput,V_bnd);

            NekVector<NekDouble> F_Int(nIntDofs, F_int,eWrapper);
            NekVector<NekDouble> V_Bnd(nLocBndDofs,V_bnd,eWrapper);

            if(nIntDofs)
            {
                m_locToGloMap.lock()->LocalToLocalInt(pLocInput,F_int);
            }

            // Boundary system solution
            if(nGlobBndDofs-nDirBndDofs)
            {
                pLocToGloMap->LocalToLocalBnd(pLocInput,F_bnd);


                // set up normalisation factor for right hand side on first SC level
                v_PreSolve(scLevel, F_bnd);

                // Gather boundary expansison into locbnd 
                NekVector<NekDouble> F_Bnd(nLocBndDofs,F_bnd1,eWrapper);
                
                // construct boundary forcing
                if(nIntDofs)
                {
                    DNekScalBlkMat &BinvD      = *m_BinvD;
                    
                    F_Bnd = BinvD*F_Int; 

                    Vmath::Vsub(nLocBndDofs, F_bnd,1, F_bnd1,1, F_bnd,1);
               }

                if(atLastLevel)
                {
                    // Transform to new basis if required 
                    v_BasisFwdTransform(F_bnd);

                    DNekScalBlkMat &SchurCompl = *m_schurCompl;

                    v_CoeffsFwdTransform(V_bnd,V_bnd);
                        
                    // subtract dirichlet boundary forcing
                    F_Bnd = SchurCompl*V_Bnd;

                    Vmath::Vsub(nLocBndDofs, F_bnd,1, F_bnd1, 1, F_bnd,1);

                    Array<OneD, NekDouble> F_hom, pert(nGlobBndDofs,0.0);
                    
                    pLocToGloMap->AssembleBnd(F_bnd, F_bnd1);
                    
                    // Solve for difference from initial solution given inout;
                    SolveLinearSystem(nGlobBndDofs, F_bnd1, pert, pLocToGloMap,
                                      nDirBndDofs);
                    
                    Array<OneD, NekDouble> outloc = F_bnd; 
                    pLocToGloMap->GlobalToLocalBnd(pert,outloc);
                    
                    // Add back initial conditions onto difference
                    Vmath::Vadd(nLocBndDofs, V_bnd, 1, outloc, 1, V_bnd,1);

                    // Transform back to original basis
                    v_CoeffsBwdTransform(V_bnd);

                    // put final bnd solution back in output array
                    m_locToGloMap.lock()->LocalBndToLocal(V_bnd,pLocOutput);
                }
                else // Process next level of recursion for multi level SC
                {
                    AssemblyMapSharedPtr nextLevLocToGloMap = pLocToGloMap->
                        GetNextLevelLocalToGlobalMap();
                    
                    // partially assemble F for next level and
                    // reshuffle V_bnd
                    nextLevLocToGloMap->PatchAssemble     (F_bnd,F_bnd);
                    nextLevLocToGloMap->PatchLocalToGlobal(V_bnd,V_bnd);
                    
                    m_recursiveSchurCompl->Solve(F_bnd,V_bnd, nextLevLocToGloMap);
                    
                    // unpack V_bnd
                    nextLevLocToGloMap->PatchGlobalToLocal(V_bnd,V_bnd);
                    
                    // place V_bnd in output array
                    m_locToGloMap.lock()->LocalBndToLocal(V_bnd, pLocOutput);
                }
            }

            // solve interior system
            if(nIntDofs)
            {
                Array<OneD, NekDouble> V_int(nIntDofs);
                NekVector<NekDouble>   V_Int(nIntDofs, V_int ,eWrapper);

                // get array of local solutions
                DNekScalBlkMat &invD  = *m_invD;
                DNekScalBlkMat &C     = *m_C;
                    
                F_Int = F_Int - C*V_Bnd;
                
                Multiply(V_Int, invD, F_Int);
                
                m_locToGloMap.lock()->LocalIntToLocal(V_int, pLocOutput);
            }
        }

        /**
         * If at the last level of recursion (or the only level in the case of
         * single-level static condensation), assemble the Schur complement.
         * For other levels, in the case of multi-level static condensation,
         * the next level of the condensed system is computed.
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysStaticCond::v_Initialise(
                const std::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            int nLocalBnd = m_locToGloMap.lock()->GetNumLocalBndCoeffs();
            int nIntDofs = m_locToGloMap.lock()->
                GetNumLocalCoeffs() - nLocalBnd; 
            m_wsp = Array<OneD, NekDouble>
                    (3*nLocalBnd+nIntDofs, 0.0);

            if (pLocToGloMap->AtLastLevel())
            {
                v_AssembleSchurComplement(pLocToGloMap);
            }
            else
            {
                ConstructNextLevelCondensedSystem(
                        pLocToGloMap->GetNextLevelLocalToGlobalMap());
            }
        }

        int GlobalLinSysStaticCond::v_GetNumBlocks()
        {
            return m_schurCompl->GetNumberOfBlockRows();
        }

        /**
         * For the first level in multi-level static condensation, or the only
         * level in the case of single-level static condensation, allocate the
         * condensed matrices and populate them with the local matrices
         * retrieved from the expansion list.
         * @param
         */
        void GlobalLinSysStaticCond::SetupTopLevel(
                const std::shared_ptr<AssemblyMap>& pLocToGloMap)
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
                if (m_linSysKey.GetMatrixType() ==
                        StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr loc_mat
                        = GlobalLinSys::v_GetBlock(n);
                    m_schurCompl->SetBlock(n,n,loc_mat);
                }
                else
                {
                    DNekScalBlkMatSharedPtr loc_schur
                        = GlobalLinSys::v_GetStaticCondBlock(n);
                    DNekScalMatSharedPtr t;
                    m_schurCompl->SetBlock(n, n, t = loc_schur->GetBlock(0,0));
                    m_BinvD     ->SetBlock(n, n, t = loc_schur->GetBlock(0,1));
                    m_C         ->SetBlock(n, n, t = loc_schur->GetBlock(1,0));
                    m_invD      ->SetBlock(n, n, t = loc_schur->GetBlock(1,1));
                }
            }
        }

        /**
         *
         */
        void GlobalLinSysStaticCond::ConstructNextLevelCondensedSystem(
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

                // Use symmetric storage for invD if possible
                MatrixStorage storageTypeD = eFULL;
                if ( (m_linSysKey.GetMatrixType() == StdRegions::eMass)      ||
                     (m_linSysKey.GetMatrixType() == StdRegions::eLaplacian) ||
                     (m_linSysKey.GetMatrixType() == StdRegions::eHelmholtz))
                {
                    storageTypeD = eSYMMETRIC;
                }

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
                                                        tmparray, wType,
                                                        storageTypeD);

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

                for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
                {
                    schurComplSubMat      = SchurCompl->GetBlock(n,n);
                    schurComplSubMatnRows = schurComplSubMat->GetRows();

                    patchId  = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetPatchId() + cnt;
                    dofId    = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetDofId()   + cnt;
                    isBndDof = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->IsBndDof() + cnt;
                    sign     = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetSign() + cnt;

                    // Set up  Matrix;
                    for(i = 0; i < schurComplSubMatnRows; ++i)
                    {
                        int pId = patchId[i];
                        Array<OneD, NekDouble> subMat0
                            = substructuredMat[0][pId]->GetPtr();
                        Array<OneD, NekDouble> subMat1
                            = substructuredMat[1][patchId[i]]->GetPtr();
                        Array<OneD, NekDouble> subMat2
                            = substructuredMat[2][patchId[i]]->GetPtr();
                        DNekMatSharedPtr subMat3
                            = substructuredMat[3][patchId[i]];
                        int subMat0rows = substructuredMat[0][pId]->GetRows();
                        int subMat1rows = substructuredMat[1][pId]->GetRows();
                        int subMat2rows = substructuredMat[2][pId]->GetRows();

                        if(isBndDof[i])
                        {
                            for(j = 0; j < schurComplSubMatnRows; ++j)
                            {
                                ASSERTL0(patchId[i]==patchId[j],
                                         "These values should be equal");

                                if(isBndDof[j])
                                {
                                    subMat0[dofId[i]+dofId[j]*subMat0rows] +=
                                        sign[i]*sign[j]*
                                            (*schurComplSubMat)(i,j);
                                }
                                else
                                {
                                    subMat1[dofId[i]+dofId[j]*subMat1rows] +=
                                        sign[i]*sign[j]*
                                            (*schurComplSubMat)(i,j);
                                }
                            }
                        }
                        else
                        {
                            for(j = 0; j < schurComplSubMatnRows; ++j)
                            {
                                ASSERTL0(patchId[i]==patchId[j],
                                         "These values should be equal");

                                if(isBndDof[j])
                                {
                                    subMat2[dofId[i]+dofId[j]*subMat2rows] +=
                                        sign[i]*sign[j]*
                                            (*schurComplSubMat)(i,j);
                                }
                                else
                                {
                                    if (storageTypeD == eSYMMETRIC)
                                    {
                                        if (dofId[i] <= dofId[j])
                                        {
                                            (*subMat3)(dofId[i],dofId[j]) +=
                                                sign[i]*sign[j]*
                                                (*schurComplSubMat)(i,j);
                                        }
                                    }
                                    else
                                    {
                                        (*subMat3)(dofId[i],dofId[j]) +=
                                            sign[i]*sign[j]*
                                            (*schurComplSubMat)(i,j);
                                    }
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
                        Array<OneD, NekDouble> subMat0
                            = substructuredMat[0][i]->GetPtr();
                        Array<OneD, NekDouble> subMat1
                            = substructuredMat[1][i]->GetPtr();
                        Array<OneD, NekDouble> subMat2
                            = substructuredMat[2][i]->GetPtr();
                        int subMat0rows = substructuredMat[0][i]->GetRows();
                        int subMat1rows = substructuredMat[1][i]->GetRows();
                        int subMat2rows = substructuredMat[2][i]->GetRows();
                        int subMat2cols = substructuredMat[2][i]->GetColumns();

                        // 1. D -> InvD
                        substructuredMat[3][i]->Invert();
                        // 2. B -> BInvD
                        (*substructuredMat[1][i]) = (*substructuredMat[1][i])*
                            (*substructuredMat[3][i]);
                        // 3. A -> A - BInvD*C (= schurcomplement)
                        // (*substructuredMat[0][i]) =(*substructuredMat[0][i])-
                        // (*substructuredMat[1][i])*(*substructuredMat[2][i]);
                        // Note: faster to use blas directly
                        Blas::Dgemm('N','N', subMat1rows, subMat2cols,
                                    subMat2rows, -1.0, &subMat1[0], subMat1rows,
                                    &subMat2[0], subMat2rows, 1.0, &subMat0[0],
                                    subMat0rows);
                    }
                }

                // STEP 3: fill the block matrices. However, do note that we
                // first have to convert them to a DNekScalMat in order to be
                // compatible with the first level of static condensation
                const Array<OneD,const unsigned int>& nbdry_size
                    = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nint_size
                    = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();
                MatrixStorage blkmatStorage = eDIAGONAL;

                blkMatrices[0] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
                blkMatrices[1] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
                blkMatrices[2] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
                blkMatrices[3] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

                DNekScalMatSharedPtr tmpscalmat;
                for(i = 0; i < nPatches; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        tmpscalmat= MemoryManager<DNekScalMat>
                            ::AllocateSharedPtr(1.0,substructuredMat[j][i]);
                        blkMatrices[j]->SetBlock(i,i,tmpscalmat);
                    }
                }
            }

            // We've finished with the Schur complement matrix passed to this
            // level, so return the memory to the system. The Schur complement
            // matrix need only be retained at the last level. Save the other
            // matrices at this level though.
            m_schurCompl.reset();

            m_recursiveSchurCompl = v_Recurse(
                m_linSysKey, m_expList, blkMatrices[0], blkMatrices[1],
                blkMatrices[2], blkMatrices[3], pLocToGloMap);
        }
    }
}
