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

            int nDofs     = pLocToGloMap->GetNumGlobalCoeffs();
            int NumDirBCs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            unsigned int rows = nDofs - NumDirBCs;
            unsigned int cols = nDofs - NumDirBCs;

            // fill global matrix
            DNekScalMatSharedPtr loc_mat;

            MatCreate(PETSC_COMM_WORLD, &m_matrix);
            MatSetType(m_matrix, MATSEQAIJ);
            MatSetSizes(m_matrix, rows, cols, PETSC_DETERMINE, PETSC_DETERMINE);
            MatSetFromOptions(m_matrix);
            MatSeqAIJSetPreallocation(m_matrix, 2500, NULL);
            //MatSetUp(m_matrix);
            MatGetVecs(m_matrix, &m_x, &m_b);

            int loc_lda;
            for(n = cnt = 0; n < m_expList.lock()->GetNumElmts(); ++n)
            {
                cout << n << endl;
                loc_mat = GetBlock(m_expList.lock()->GetOffset_Elmt_Id(n));
                loc_lda = loc_mat->GetRows();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt + i)-NumDirBCs;
                    sign1 = pLocToGloMap->GetLocalToGlobalSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = pLocToGloMap->GetLocalToGlobalMap(cnt + j)
                                                                    - NumDirBCs;
                            sign2 = pLocToGloMap->GetLocalToGlobalSign(cnt + j);
                            if(gid2 >= 0)
                            {
                                value = sign1*sign2*(*loc_mat)(i,j);
                                MatSetValue(
                                    m_matrix, gid1, gid2, value, ADD_VALUES);
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            MatAssemblyBegin(m_matrix,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(m_matrix,MAT_FINAL_ASSEMBLY);

            KSPCreate(PETSC_COMM_WORLD, &m_ksp);

            PC pc;
            //KSPSetType(m_ksp, KSPCG);
            KSPSetTolerances(
                m_ksp, pLocToGloMap->GetIterativeTolerance(),
                PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
            //KSPGetPC(m_ksp, &pc);
            //PCSetType(pc, PCGAMG);
            KSPSetFromOptions(m_ksp);
            KSPSetOperators(m_ksp, m_matrix, m_matrix, DIFFERENT_NONZERO_PATTERN);
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
