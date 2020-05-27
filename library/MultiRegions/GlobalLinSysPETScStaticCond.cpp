///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysPETScStaticCond.cpp
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
// Description: GlobalLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/GlobalLinSysPETScStaticCond.h>

#include <petscsys.h>
#include <petscksp.h>
#include <petscmat.h>

using namespace std;

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
                     const GlobalLinSysKey                &pKey,
                     const std::weak_ptr<ExpList>         &pExpList,
                     const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
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
                     const GlobalLinSysKey                &pKey,
                     const std::weak_ptr<ExpList>         &pExpList,
                     const DNekScalBlkMatSharedPtr         pSchurCompl,
                     const DNekScalBlkMatSharedPtr         pBinvD,
                     const DNekScalBlkMatSharedPtr         pC,
                     const DNekScalBlkMatSharedPtr         pInvD,
                     const std::shared_ptr<AssemblyMap>   &pLocToGloMap,
                     const PreconditionerSharedPtr         pPrecon)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysPETSc     (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            m_schurCompl = pSchurCompl;
            m_BinvD      = pBinvD;
            m_C          = pC;
            m_invD       = pInvD;
            m_precon     = pPrecon;
        }

        /**
         *
         */
        GlobalLinSysPETScStaticCond::~GlobalLinSysPETScStaticCond()
        {

        }

        void GlobalLinSysPETScStaticCond::v_InitObject()
        {
            auto asmMap = m_locToGloMap.lock();

            m_precon = CreatePrecon(asmMap);

            // Allocate memory for top-level structure
            SetupTopLevel(asmMap);

            // Setup Block Matrix systems
            int n, n_exp = m_expList.lock()->GetNumElmts();

            // Build preconditioner
            m_precon->BuildPreconditioner();

            // Do transform of Schur complement matrix
            int cnt = 0;
            for (n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() !=
                        StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr mat = m_schurCompl->GetBlock(n, n);
                    DNekScalMatSharedPtr t = m_precon->TransformedSchurCompl(
                                                             n, cnt, mat);
                    m_schurCompl->SetBlock(n, n, t);
                    cnt += mat->GetRows();
                }
            }

            // Construct this level
            Initialise(asmMap);
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

            // Build precon again if we in multi-level static condensation (a
            // bit of a hack)
            if (m_linSysKey.GetGlobalSysSolnType() ==
                    ePETScMultiLevelStaticCond)
            {
                m_precon = CreatePrecon(m_locToGloMap.lock());
                m_precon->BuildPreconditioner();
            }

            // CALCULATE REORDERING MAPPING
            CalculateReordering(pLocToGloMap->GetGlobalToUniversalBndMap(),
                                pLocToGloMap->GetGlobalToUniversalBndMapUnique(),
                                pLocToGloMap);

            // SET UP VECTORS AND MATRIX
            SetUpMatVec(pLocToGloMap->GetNumGlobalBndCoeffs(), nDirDofs);

            // SET UP SCATTER OBJECTS
            SetUpScatter();

            // CONSTRUCT KSP OBJECT
            SetUpSolver(pLocToGloMap->GetIterativeTolerance());

            // If we are using the matrix multiplication shell don't try to
            // populate the matrix.
            if (m_matMult == ePETScMatMultShell)
            {
                return;
            }

            // POPULATE MATRIX
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
        }

        DNekScalBlkMatSharedPtr GlobalLinSysPETScStaticCond::
            v_GetStaticCondBlock(unsigned int n)
        {
            DNekScalBlkMatSharedPtr schurComplBlock;
            DNekScalMatSharedPtr    localMat = m_schurCompl->GetBlock(n,n);
            unsigned int nbdry    = localMat->GetRows();
            unsigned int nblks    = 1;
            unsigned int esize[1] = {nbdry};

            schurComplBlock = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nblks, nblks, esize, esize);
            schurComplBlock->SetBlock(0, 0, localMat);

            return schurComplBlock;
        }

        void GlobalLinSysPETScStaticCond::v_PreSolve(
            int                     scLevel,
            Array<OneD, NekDouble>   &F_bnd)
        {
            boost::ignore_unused(F_bnd);

            if (scLevel == 0)
            {
                // When matrices are supplied to the constructor at the top
                // level, the preconditioner is never set up.
                if (!m_precon)
                {
                    m_precon = CreatePrecon(m_locToGloMap.lock());
                    m_precon->BuildPreconditioner();
                }
            }


        }
        
        void GlobalLinSysPETScStaticCond::v_BasisFwdTransform(
                                     Array<OneD, NekDouble>& pInOut)
        {
            m_precon->DoTransformBasisToLowEnergy(pInOut);            
        }

        void GlobalLinSysPETScStaticCond::v_CoeffsBwdTransform(
            Array<OneD, NekDouble>& pInOut)
        {
	    m_precon->DoTransformCoeffsFromLowEnergy(pInOut);
        }

        void GlobalLinSysPETScStaticCond::v_CoeffsFwdTransform(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput)
        {
	    m_precon->DoTransformCoeffsToLowEnergy(pInput,pOutput);
        }

        /**
         * @brief Apply matrix-vector multiplication using local approach and
         * the assembly map.
         *
         * @param input   Vector input.
         * @param output  Result of multiplication.
         *
         * @todo This can possibly be made faster by using the sparse
         *       block-matrix multiplication code from the iterative elastic
         *       systems.
         */
        void GlobalLinSysPETScStaticCond::v_DoMatrixMultiply(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output)
        {
            auto asmMap = m_locToGloMap.lock();

            int nLocBndDofs = asmMap->GetNumLocalBndCoeffs();
            int nDirDofs    = asmMap->GetNumGlobalDirBndCoeffs();

            NekVector<NekDouble> in(nLocBndDofs), out(nLocBndDofs);
            asmMap->GlobalToLocalBnd(input, in.GetPtr(), nDirDofs);
            out = (*m_schurCompl) * in;
            asmMap->AssembleBnd(out.GetPtr(), output, nDirDofs);
        }

        GlobalLinSysStaticCondSharedPtr GlobalLinSysPETScStaticCond::v_Recurse(
            const GlobalLinSysKey                &mkey,
            const std::weak_ptr<ExpList>         &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const std::shared_ptr<AssemblyMap>   &l2gMap)
        {
            GlobalLinSysPETScStaticCondSharedPtr sys = MemoryManager<
                GlobalLinSysPETScStaticCond>::AllocateSharedPtr(
                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap,
                    m_precon);
            sys->Initialise(l2gMap);
            return sys;
        }
    }
}
