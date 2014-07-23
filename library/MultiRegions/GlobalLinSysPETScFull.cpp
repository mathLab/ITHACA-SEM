///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysPETScFull.cpp
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
// Description: GlobalLinSysPETScFull definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysPETScFull.h>
#include <MultiRegions/ExpList.h>

#include "petscao.h"
#include "petscis.h"

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysPETScFull
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysPETScFull::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "PETScFull",
                    GlobalLinSysPETScFull::create,
                    "PETSc Full Matrix.");


        /// Constructor for full direct matrix solve.
        GlobalLinSysPETScFull::GlobalLinSysPETScFull(
            const GlobalLinSysKey                &pLinSysKey,
            const boost::weak_ptr<ExpList>       &pExp,
            const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
                : GlobalLinSysPETSc(pLinSysKey, pExp, pLocToGloMap)
        {
            ASSERTL1(m_linSysKey.GetGlobalSysSolnType() == ePETScFullMatrix,
                     "This routine should only be used when using a Full PETSc"
                     " matrix solve");

            int i, j, n, cnt, gid1, gid2;
            NekDouble sign1, sign2, value;

            int nDofs    = pLocToGloMap->GetNumGlobalCoeffs();
            int nDirDofs = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nHomDofs = nDofs - nDirDofs;
            
            // fill global matrix
            DNekScalMatSharedPtr loc_mat;

            LibUtilities::CommSharedPtr vComm
                = pExp.lock()->GetSession()->GetComm();
            const Array<OneD, const int>& glo2unique
                = pLocToGloMap->GetGlobalToUniversalMapUnique();
            int nLocal = Vmath::Vsum(nHomDofs, glo2unique + nDirDofs, 1);

            int nProc = vComm->GetSize();
            int rank  = vComm->GetRank();

            m_reorderedMap.resize(nHomDofs);

            Array<OneD, int> localCounts(nProc, 0);
            Array<OneD, int> localOffset(nProc, 0);
            localCounts[rank] = nHomDofs;
            vComm->AllReduce(localCounts, LibUtilities::ReduceSum);

            for (i = 1; i < nProc; ++i)
            {
                localOffset[i] = localOffset[i-1] + localCounts[i-1];
            }

            int totHomDofs = Vmath::Vsum(nProc, localCounts, 1);
            vector<unsigned int> allUniIds(totHomDofs, 0);

            for (n = cnt = 0; n < m_expList.lock()->GetNumElmts(); ++n)
            {
                loc_mat = GetBlock(m_expList.lock()->GetOffset_Elmt_Id(n));
                int loc_lda = loc_mat->GetRows();

                for (i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt+i) - nDirDofs;

                    if (gid1 >= 0)
                    {
                        allUniIds[gid1+localOffset[rank]] =
                            pLocToGloMap->GetGlobalToUniversalMap(gid1 + nDirDofs);
                    }
                }
                cnt += loc_lda;
            }

            vComm->AllReduce(allUniIds, LibUtilities::ReduceSum);
            std::sort(allUniIds.begin(), allUniIds.end());
            map<int,int> uniIdReorder;

            for (cnt = n = 0; n < allUniIds.size(); ++n)
            {
                if (uniIdReorder.count(allUniIds[n]) > 0)
                {
                    continue;
                }

                uniIdReorder[allUniIds[n]] = cnt++;
            }

            map<int,int>::iterator mIt;

            for (mIt = uniIdReorder.begin(); mIt != uniIdReorder.end(); ++mIt)
            {
                cout << "RANK " << vComm->GetRank() << ": " << mIt->first
                     << " -> " << mIt->second << endl;
            }

            for (n = cnt = 0; n < m_expList.lock()->GetNumElmts(); ++n)
            {
                loc_mat = GetBlock(m_expList.lock()->GetOffset_Elmt_Id(n));
                int loc_lda = loc_mat->GetRows();

                for (i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt+i) - nDirDofs;

                    if (gid1 >= 0)
                    {
                        int uniid = pLocToGloMap->GetGlobalToUniversalMap(gid1 + nDirDofs);
                        ASSERTL0(uniIdReorder.count(uniid) > 0, "wat");
                        m_reorderedMap[gid1] = uniIdReorder[uniid];
                    }
                }
                cnt += loc_lda;
            }

            // CREATE VECTORS
            VecCreate        (PETSC_COMM_WORLD, &m_x);
            VecSetSizes      (m_x, nLocal, PETSC_DECIDE);
            VecSetFromOptions(m_x);
            VecDuplicate     (m_x, &m_b);

            // CREATE MATRICES
            MatCreate        (PETSC_COMM_WORLD, &m_matrix);
            MatSetType       (m_matrix, MATMPIAIJ);
            MatSetSizes      (m_matrix, nLocal, nLocal,
                              PETSC_DETERMINE, PETSC_DETERMINE);
            MatSetFromOptions(m_matrix);
            MatSetUp         (m_matrix);

#if 0
            // preallocate
            int loc_lda;
            Array<OneD, int> nnz(nDofs-NumDirBCs, 0);
            vector<set<int> > hack(rows);

            for(n = cnt = 0; n < m_expList.lock()->GetNumElmts(); ++n)
            {
                loc_mat = GetBlock(m_expList.lock()->GetOffset_Elmt_Id(n));
                loc_lda = loc_mat->GetRows();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt + i)-NumDirBCs;
                    if(gid1 < 0)
                    {
                        continue;
                    }

                    for(j = 0; j < loc_lda; ++j)
                    {
                        gid2 = pLocToGloMap->GetLocalToGlobalMap(cnt + j)
                            - NumDirBCs;
                        if(gid2 < 0)
                        {
                            continue;
                        }

                        if (hack[gid1].count(gid2) > 0)
                        {
                            continue;
                        }

                        hack[gid1].insert(gid2);
                        nnz[gid1]++;
                    }
                }
                cnt += loc_lda;
            }

            MatSeqAIJSetPreallocation(m_matrix, 0, nnz.get());
#endif

            int loc_lda;

            for(n = cnt = 0; n < m_expList.lock()->GetNumElmts(); ++n)
            {
                loc_mat = GetBlock(m_expList.lock()->GetOffset_Elmt_Id(n));
                loc_lda = loc_mat->GetRows();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt+i) - nDirDofs;
                    sign1 = pLocToGloMap->GetLocalToGlobalSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        int gid1ro = m_reorderedMap[gid1];
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = pLocToGloMap->GetLocalToGlobalMap(cnt + j)
                                                                     - nDirDofs;
                            sign2 = pLocToGloMap->GetLocalToGlobalSign(cnt + j);
                            if(gid2 >= 0)
                            {
                                int gid2ro = m_reorderedMap[gid2];
                                value = sign1*sign2*(*loc_mat)(i,j);
                                MatSetValue(
                                    m_matrix, gid1ro, gid2ro, value, ADD_VALUES);
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            // ASSEMBLE MATRIX
            MatAssemblyBegin(m_matrix,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(m_matrix,MAT_FINAL_ASSEMBLY);

            MatView(m_matrix,  PETSC_VIEWER_STDOUT_WORLD);

            // CONSTRUCT KSP OBJECT
            KSPCreate(PETSC_COMM_WORLD, &m_ksp);
            KSPSetTolerances(
                m_ksp, pLocToGloMap->GetIterativeTolerance(),
                PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
            KSPSetFromOptions(m_ksp);
            KSPSetOperators(m_ksp, m_matrix, m_matrix);
        }


        GlobalLinSysPETScFull::~GlobalLinSysPETScFull()
        {

        }


        /**
         * Solve the linear system using a full global matrix system.
         */
        void GlobalLinSysPETScFull::v_Solve(
                    const Array<OneD, const NekDouble>  &pInput,
                          Array<OneD,       NekDouble>  &pOutput,
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &pDirForcing)
        {
            bool dirForcCalculated = (bool) pDirForcing.num_elements();
            int nDirDofs  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobDofs = pLocToGloMap->GetNumGlobalCoeffs();
            Array<OneD, NekDouble> tmp(nGlobDofs), tmp2;
            
            if(nDirDofs)
            {
                // calculate the dirichlet forcing
                if(dirForcCalculated)
                {
                    Vmath::Vsub(nGlobDofs,
                                pInput.get(),      1,
                                pDirForcing.get(), 1,
                                tmp.get(),         1);
                }
                else
                {
                    // Calculate the dirichlet forcing and substract it
                    // from the rhs
                    m_expList.lock()->GeneralMatrixOp(
                        m_linSysKey, pOutput, tmp, eGlobal);
                    
                    Vmath::Vsub(nGlobDofs, 
                                pInput.get(), 1,
                                tmp.get(),    1,
                                tmp.get(),    1);
                }

                SolveLinearSystem(nGlobDofs, tmp,
                                  tmp2 = pOutput + nDirDofs,
                                  pLocToGloMap, nDirDofs);
            }
            else
            {
                SolveLinearSystem(nDirDofs, pInput, pOutput, pLocToGloMap);
            }
        }
    }
}
