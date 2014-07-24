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

#include <MultiRegions/GlobalLinSysPETScStaticCond.h>

#include <petscsys.h>
#include <petscksp.h>
#include <petscmat.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysPETSc
         *
         * Solves a linear system using single- or multi-level static
         * condensation.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysPETScStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "PETScStaticCond",
                    GlobalLinSysPETScStaticCond::create,
                    "PETSc static condensation.");

        string GlobalLinSysPETScStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "PETScMultiLevelStaticCond",
                    GlobalLinSysPETScStaticCond::create,
                    "PETSc multi-level static condensation.");

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
        GlobalLinSysPETScStaticCond::GlobalLinSysPETScStaticCond(
                     const GlobalLinSysKey &pKey,
                     const boost::weak_ptr<ExpList> &pExpList,
                     const boost::shared_ptr<AssemblyMap>
                     &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysPETSc     (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==ePETScStaticCond)||
                     (pKey.GetGlobalSysSolnType()==ePETScMultiLevelStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(pKey.GetGlobalSysSolnType()
                        == pLocToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");
        }

        /**
         *
         */
        GlobalLinSysPETScStaticCond::GlobalLinSysPETScStaticCond(
                     const GlobalLinSysKey &pKey,
                     const boost::weak_ptr<ExpList> &pExpList,
                     const DNekScalBlkMatSharedPtr pSchurCompl,
                     const DNekScalBlkMatSharedPtr pBinvD,
                     const DNekScalBlkMatSharedPtr pC,
                     const DNekScalBlkMatSharedPtr pInvD,
                     const boost::shared_ptr<AssemblyMap>
                                                            &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysPETSc     (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            m_schurCompl = pSchurCompl;
            m_BinvD      = pBinvD;
            m_C          = pC;
            m_invD       = pInvD;
        }

        /**
         *
         */
        GlobalLinSysPETScStaticCond::~GlobalLinSysPETScStaticCond()
        {

        }

        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysPETScStaticCond::v_AssembleSchurComplement(
            AssemblyMapSharedPtr pLocToGloMap)
        {
            int i, j, n, cnt, gid1, gid2, loc_lda;
            NekDouble sign1, sign2, value;

            const int nDirDofs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            DNekScalBlkMatSharedPtr SchurCompl = m_schurCompl;
            DNekScalBlkMatSharedPtr BinvD      = m_BinvD;
            DNekScalBlkMatSharedPtr C          = m_C;
            DNekScalBlkMatSharedPtr invD       = m_invD;
            DNekScalMatSharedPtr    loc_mat;

            // CALCULATE REORDERING MAPPING
            CalculateReordering(pLocToGloMap->GetGlobalToUniversalBndMap(),
                                pLocToGloMap->GetGlobalToUniversalBndMapUnique(),
                                pLocToGloMap);

            // SET UP VECTORS AND MATRIX
            SetUpMatVec();

            for(n = cnt = 0; n < m_schurCompl->GetNumberOfBlockRows(); ++n)
            {
                loc_mat = m_schurCompl->GetBlock(n,n);
                loc_lda = loc_mat->GetRows();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalBndMap(cnt + i)-nDirDofs;
                    sign1 = pLocToGloMap->GetLocalToGlobalBndSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        int gid1ro = m_reorderedMap[gid1];
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = pLocToGloMap->GetLocalToGlobalBndMap(cnt + j)
                                                                    - nDirDofs;
                            sign2 = pLocToGloMap->GetLocalToGlobalBndSign(cnt + j);
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
                cnt   += loc_lda;
            }

            // ASSEMBLE MATRIX
            MatAssemblyBegin(m_matrix, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd  (m_matrix, MAT_FINAL_ASSEMBLY);

            // SET UP SCATTER OBJECTS
            SetUpScatter();

            // CONSTRUCT KSP OBJECT
            SetUpSolver(pLocToGloMap->GetIterativeTolerance());
        }

        GlobalLinSysStaticCondSharedPtr GlobalLinSysPETScStaticCond::v_Recurse(
            const GlobalLinSysKey                &mkey,
            const boost::weak_ptr<ExpList>       &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const boost::shared_ptr<AssemblyMap> &l2gMap)
        {
            GlobalLinSysPETScStaticCondSharedPtr sys = MemoryManager<
                GlobalLinSysPETScStaticCond>::AllocateSharedPtr(
                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap);
            sys->Initialise(l2gMap);
            return sys;
        }
    }
}
