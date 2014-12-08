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
            : GlobalLinSys     (pLinSysKey, pExp, pLocToGloMap),
              GlobalLinSysPETSc(pLinSysKey, pExp, pLocToGloMap)
        {
            ASSERTL1(m_linSysKey.GetGlobalSysSolnType() == ePETScFullMatrix,
                     "This routine should only be used when using a Full PETSc"
                     " matrix solve");

            const int nDirDofs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            int i, j, n, cnt, gid1, gid2, loc_lda;
            NekDouble sign1, sign2, value;
            DNekScalMatSharedPtr loc_mat;

            // CALCULATE REORDERING MAPPING
            CalculateReordering(pLocToGloMap->GetGlobalToUniversalMap(),
                                pLocToGloMap->GetGlobalToUniversalMapUnique(),
                                pLocToGloMap);

            // SET UP VECTORS AND MATRIX
            SetUpMatVec();

            // POPULATE MATRIX
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
            MatAssemblyBegin(m_matrix, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd  (m_matrix, MAT_FINAL_ASSEMBLY);

            // SET UP SCATTER OBJECTS
            SetUpScatter();

            // CONSTRUCT KSP OBJECT
            SetUpSolver(pLocToGloMap->GetIterativeTolerance());
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

            int nDirTotal = nDirDofs;
            m_expList.lock()->GetComm()->GetRowComm()->AllReduce(
                nDirTotal, LibUtilities::ReduceSum);
            
            if(nDirTotal)
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

                Array<OneD, NekDouble> out(nGlobDofs,0.0);
                SolveLinearSystem(nGlobDofs, tmp, out, pLocToGloMap, nDirDofs);
                Vmath::Vadd(nGlobDofs-nDirDofs,    &out    [nDirDofs], 1,
                            &pOutput[nDirDofs], 1, &pOutput[nDirDofs], 1);
            }
            else
            {
                SolveLinearSystem(nDirDofs, pInput, pOutput, pLocToGloMap);
            }
        }
    }
}
