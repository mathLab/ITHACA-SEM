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

#include <MultiRegions/GlobalLinSysPETSc.h>
#include <MultiRegions/Preconditioner.h>

#include "petscis.h"
#include "petscversion.h"

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        std::string GlobalLinSysPETSc::matMult =
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                "PETScMatMult", "Sparse");
        std::string GlobalLinSysPETSc::matMultIds[] = {
            LibUtilities::SessionReader::RegisterEnumValue(
                "PETScMatMult", "Sparse",
                MultiRegions::ePETScMatMultSparse),
            LibUtilities::SessionReader::RegisterEnumValue(
                "PETScMatMult", "Shell",
                MultiRegions::ePETScMatMultShell)
        };

        /**
         * @class GlobalLinSysPETSc
         *
         * Solves a linear system using PETSc.
         */
        GlobalLinSysPETSc::GlobalLinSysPETSc(
            const GlobalLinSysKey                &pKey,
            const std::weak_ptr<ExpList>         &pExp,
            const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
            : GlobalLinSys(pKey, pExp, pLocToGloMap)
        {
            // Determine whether to use standard sparse matrix approach or
            // shell.
            m_matMult = pExp.lock()->GetSession()
                ->GetSolverInfoAsEnum<PETScMatMult>(
                    "PETScMatMult");

            // Check PETSc is initialized. For some reason, this is needed on
            // OS X as logging is not initialized properly in the call within
            // CommMpi.
            PetscBool isInitialized;
            PetscInitialized(&isInitialized);
            if (!isInitialized)
            {
#ifdef NEKTAR_USE_MPI
                std::string commType =
                    m_expList.lock()->GetSession()->GetComm()->GetType();
                if (commType.find("MPI") != std::string::npos)
                {
                    LibUtilities::CommMpiSharedPtr comm =
                        std::static_pointer_cast<LibUtilities::CommMpi>(
                            m_expList.lock()->GetSession()->GetComm());
                    PETSC_COMM_WORLD = comm->GetComm();
                }
#endif
                PetscInitializeNoArguments();
            }

            // Create matrix
            MatCreate(PETSC_COMM_WORLD, &m_matrix);
        }

        /**
         * @brief Clean up PETSc objects.
         *
         * Note that if SessionReader::Finalize is called before the end of the
         * program, PETSc may have been finalized already, at which point we
         * cannot deallocate our objects. If that's the case we do nothing and
         * let the kernel clear up after us.
         */
        GlobalLinSysPETSc::~GlobalLinSysPETSc()
        {
            PetscBool isFinalized;
            PetscFinalized(&isFinalized);

            // Sometimes, PetscFinalized returns false when (in fact) CommMpi's
            // Finalise routine has been called. We therefore also need to check
            // whether MPI has been finalised. This might arise from the
            // additional call to PetscInitializeNoArguments in the constructor
            // above.
#ifdef NEKTAR_USE_MPI
            int mpiFinal = 0;
            MPI_Finalized(&mpiFinal);
            isFinalized = isFinalized || mpiFinal ? PETSC_TRUE : PETSC_FALSE;
#endif

            if (!isFinalized)
            {
                KSPDestroy(&m_ksp);
                PCDestroy (&m_pc);
                MatDestroy(&m_matrix);
                VecDestroy(&m_x);
                VecDestroy(&m_b);
                VecDestroy(&m_locVec);
            }
        }

        /**
         * @brief Solve linear system using PETSc.
         *
         * The general strategy being a PETSc solve is to:
         *
         * - Copy values into the PETSc vector #m_b
         * - Solve the system #m_ksp and place result into #m_x.
         * - Scatter results back into #m_locVec using #m_ctx scatter object.
         * - Copy from #m_locVec to output array #pOutput.
         */
        void GlobalLinSysPETSc::v_SolveLinearSystem(
            const int                          pNumRows,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr        &locToGloMap,
            const int                          pNumDir)
        {
            const int nHomDofs = pNumRows - pNumDir;

            if (!m_precon && m_matMult == ePETScMatMultShell)
            {
                m_precon = CreatePrecon(locToGloMap);
                m_precon->BuildPreconditioner();
            }

            // Populate RHS vector from input
            VecSetValues(m_b, nHomDofs, &m_reorderedMap[0],
                         &pInput[pNumDir], INSERT_VALUES);

            // Assemble RHS vector
            VecAssemblyBegin(m_b);
            VecAssemblyEnd  (m_b);

            // Do system solve
            KSPSolve(m_ksp, m_b, m_x);

            KSPConvergedReason reason;
            KSPGetConvergedReason(m_ksp, &reason);
            ASSERTL0(reason > 0,
                    "PETSc solver diverged, reason is: " +
                        std::string(KSPConvergedReasons[reason]));


            // Scatter results to local vector
            VecScatterBegin(m_ctx, m_x, m_locVec,
                            INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd  (m_ctx, m_x, m_locVec,
                            INSERT_VALUES, SCATTER_FORWARD);

            // Copy results into output vector
            PetscScalar *tmp;
            VecGetArray    (m_locVec, &tmp);
            Vmath::Vcopy   (nHomDofs, tmp, 1, &pOutput[pNumDir], 1);
            VecRestoreArray(m_locVec, &tmp);
        }

        /**
         * @brief Set up PETSc local (equivalent to Nektar++ global) and global
         * (equivalent to universal) scatter maps.
         *
         * These maps are used in GlobalLinSysPETSc::v_SolveLinearSystem to
         * scatter the solution vector back to each process.
         */
        void GlobalLinSysPETSc::SetUpScatter()
        {
            const int nHomDofs = m_reorderedMap.size();

            // Create local and global numbering systems for vector
            IS isGlobal, isLocal;
            ISCreateGeneral  (PETSC_COMM_SELF, nHomDofs, &m_reorderedMap[0],
                              PETSC_COPY_VALUES, &isGlobal);
            ISCreateStride   (PETSC_COMM_SELF, nHomDofs, 0, 1, &isLocal);

            // Create local vector for output
            VecCreate        (PETSC_COMM_SELF, &m_locVec);
            VecSetSizes      (m_locVec, nHomDofs, PETSC_DECIDE);
            VecSetFromOptions(m_locVec);

            // Create scatter context
            VecScatterCreate (m_x, isGlobal, m_locVec, isLocal, &m_ctx);

            // Clean up
            ISDestroy(&isGlobal);
            ISDestroy(&isLocal);
        }

        /**
         * @brief Calculate a reordering of universal IDs for PETSc.
         *
         * PETSc requires a unique, contiguous index of all global and universal
         * degrees of freedom which represents its position inside the
         * matrix. Presently Gs does not guarantee this, so this routine
         * constructs a new universal mapping.
         *
         * @param glo2uniMap    Global to universal map
         * @param glo2unique    Global to unique map
         * @param pLocToGloMap  Assembly map for this system
         */
        void GlobalLinSysPETSc::CalculateReordering(
            const Array<OneD, const int> &glo2uniMap,
            const Array<OneD, const int> &glo2unique,
            const AssemblyMapSharedPtr   &pLocToGloMap)
        {
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetSession()->GetComm();

            const int nDirDofs = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            const int nHomDofs = glo2uniMap.size() - nDirDofs;
            const int nProc    = vComm->GetSize();
            const int rank     = vComm->GetRank();

            int n, cnt;

            // Count number of unique degrees of freedom on each process.
            m_nLocal = Vmath::Vsum(nHomDofs, glo2unique + nDirDofs, 1);
            m_reorderedMap.resize(nHomDofs);

            // Reduce coefficient counts across all processors.
            Array<OneD, int> localCounts(nProc, 0), localOffset(nProc, 0);
            localCounts[rank] = nHomDofs;
            vComm->AllReduce(localCounts, LibUtilities::ReduceSum);

            for (n = 1; n < nProc; ++n)
            {
                localOffset[n] = localOffset[n-1] + localCounts[n-1];
            }

            int totHomDofs = Vmath::Vsum(nProc, localCounts, 1);
            vector<unsigned int> allUniIds(totHomDofs, 0);

            // Assemble list of universal IDs
            for (n = 0; n < nHomDofs; ++n)
            {
                int gid = n + nDirDofs;
                allUniIds[n + localOffset[rank]] = glo2uniMap[gid];
            }

            // Reduce this across processors so that each process has a list of
            // all universal IDs.
            vComm->AllReduce(allUniIds, LibUtilities::ReduceSum);
            std::sort(allUniIds.begin(), allUniIds.end());
            map<int,int> uniIdReorder;

            // Renumber starting from 0.
            for (cnt = n = 0; n < allUniIds.size(); ++n)
            {
                if (uniIdReorder.count(allUniIds[n]) > 0)
                {
                    continue;
                }

                uniIdReorder[allUniIds[n]] = cnt++;
            }

            // Populate reordering map.
            for (n = 0; n < nHomDofs; ++n)
            {
                int gid = n + nDirDofs;
                int uniId = glo2uniMap[gid];
                ASSERTL0(uniIdReorder.count(uniId) > 0, "Error in ordering");
                m_reorderedMap[n] = uniIdReorder[uniId];
            }
        }

        /**
         * @brief Construct PETSc matrix and vector handles.
         *
         * @todo Preallocation should be done at this point, since presently
         *       matrix allocation takes a significant amount of time.
         *
         * @param nGlobal  Number of global degrees of freedom in the system (on
         *                 this processor)
         * @param nDir     Number of Dirichlet degrees of freedom (on this
         *                 processor).
         */
        void GlobalLinSysPETSc::SetUpMatVec(int nGlobal, int nDir)
        {
            // CREATE VECTORS
            VecCreate        (PETSC_COMM_WORLD, &m_x);
            VecSetSizes      (m_x, m_nLocal, PETSC_DECIDE);
            VecSetFromOptions(m_x);
            VecDuplicate     (m_x, &m_b);

            // CREATE MATRICES
            if (m_matMult == ePETScMatMultShell)
            {
                // Create ShellCtx context object which will store the matrix
                // size and a pointer to the linear system. We do this so that
                // we can call a member function to the matrix-vector and
                // preconditioning multiplication in a subclass.
                ShellCtx *ctx1 = new ShellCtx(), *ctx2 = new ShellCtx();
                ctx1->nGlobal = ctx2->nGlobal = nGlobal;
                ctx1->nDir    = ctx2->nDir    = nDir;
                ctx1->linSys  = ctx2->linSys  = this;

                // Set up MatShell object.
                MatCreateShell      (PETSC_COMM_WORLD, m_nLocal, m_nLocal,
                                     PETSC_DETERMINE, PETSC_DETERMINE,
                                     (void *)ctx1, &m_matrix);
                MatShellSetOperation(m_matrix, MATOP_MULT,
                                     (void(*)(void))DoMatrixMultiply);
                MatShellSetOperation(m_matrix, MATOP_DESTROY,
                                     (void(*)(void))DoDestroyMatCtx);

                // Create a PCShell to go alongside the MatShell.
                PCCreate         (PETSC_COMM_WORLD, &m_pc);
#if PETSC_VERSION_GE(3,5,0)
                PCSetOperators   (m_pc, m_matrix, m_matrix);
#else
                PCSetOperators   (m_pc, m_matrix, m_matrix, SAME_NONZERO_PATTERN);
#endif
                PCSetType        (m_pc, PCSHELL);
                PCShellSetApply  (m_pc, DoPreconditioner);
                PCShellSetDestroy(m_pc, DoDestroyPCCtx);
                PCShellSetContext(m_pc, ctx2);
            }
            else
            {
                // Otherwise we create a PETSc matrix and use MatSetFromOptions
                // so that we can set various options on the command line.
                MatCreate        (PETSC_COMM_WORLD, &m_matrix);
                MatSetType       (m_matrix, MATAIJ);
                MatSetSizes      (m_matrix, m_nLocal, m_nLocal,
                                  PETSC_DETERMINE, PETSC_DETERMINE);
                MatSetFromOptions(m_matrix);
                MatSetUp         (m_matrix);
            }
        }

        /**
         * @brief Perform either matrix multiplication or preconditioning using
         * Nektar++ routines.
         *
         * This static function uses Nektar++ routines to calculate the
         * matrix-vector product of @p M with @p in, storing the output in @p
         * out.
         *
         * @todo There's a lot of scatters and copies that might possibly be
         *       eliminated to make this more efficient.
         *
         * @param in      Input vector.
         * @param out     Output vector.
         * @param ctx     ShellCtx object that points to our instance of
         *                GlobalLinSysPETSc.
         * @param precon  If true, we apply a preconditioner, if false, we
         *                perform a matrix multiplication.
         */
        void GlobalLinSysPETSc::DoNekppOperation(
            Vec &in, Vec &out, ShellCtx *ctx, bool precon)
        {
            const int          nGlobal  = ctx->nGlobal;
            const int          nDir     = ctx->nDir;
            const int          nHomDofs = nGlobal - nDir;
            GlobalLinSysPETSc *linSys   = ctx->linSys;

            // Scatter from PETSc ordering to our local ordering. It's actually
            // unclear whether this step might also do some communication in
            // parallel, which is probably not ideal.
            VecScatterBegin(linSys->m_ctx, in, linSys->m_locVec,
                            INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd  (linSys->m_ctx, in, linSys->m_locVec,
                            INSERT_VALUES, SCATTER_FORWARD);

            // Temporary storage to pass to Nektar++
            Array<OneD, NekDouble> tmpIn(nHomDofs), tmpOut(nHomDofs);

            // Get values from input vector and copy to tmpIn.
            PetscScalar *tmpLocIn;
            VecGetArray    (linSys->m_locVec, &tmpLocIn);
            Vmath::Vcopy   (nHomDofs, tmpLocIn, 1, &tmpIn[0], 1);
            VecRestoreArray(linSys->m_locVec, &tmpLocIn);

            // Do matrix multiply in Nektar++, store in tmpOut.
            if (precon)
            {
                linSys->m_precon->DoPreconditioner(tmpIn, tmpOut);
            }
            else
            {
                linSys->v_DoMatrixMultiply(tmpIn, tmpOut);
            }

            // Scatter back to PETSc ordering and put in out.
            VecSetValues(out, nHomDofs, &linSys->m_reorderedMap[0],
                         &tmpOut[0], INSERT_VALUES);
            VecAssemblyBegin(out);
            VecAssemblyEnd  (out);
        }

        /**
         * @brief Perform matrix multiplication using Nektar++ routines.
         *
         * This static function uses Nektar++ routines to calculate the
         * matrix-vector product of @p M with @p in, storing the output in @p
         * out.
         *
         * @param M    Original MatShell matrix, which stores the ShellCtx
         *             object.
         * @param in   Input vector.
         * @param out  Output vector.
         */
        PetscErrorCode GlobalLinSysPETSc::DoMatrixMultiply(
            Mat M, Vec in, Vec out)
        {
            // Grab our shell context from M.
            void *ptr;
            MatShellGetContext(M, &ptr);
            ShellCtx *ctx = (ShellCtx *)ptr;

            DoNekppOperation(in, out, ctx, false);

            // Must return 0, otherwise PETSc complains.
            return 0;
        }

        /**
         * @brief Apply preconditioning using Nektar++ routines.
         *
         * This static function uses Nektar++ routines to apply the
         * preconditioner stored in GlobalLinSysPETSc::m_precon from the context
         * of @p pc to the vector @p in, storing the output in @p out.
         *
         * @param pc   Preconditioner object that stores the ShellCtx.
         * @param in   Input vector.
         * @param out  Output vector.
         */
        PetscErrorCode GlobalLinSysPETSc::DoPreconditioner(
            PC pc, Vec in, Vec out)
        {
            // Grab our PCShell context from pc.
            void *ptr;
            PCShellGetContext(pc, &ptr);
            ShellCtx *ctx = (ShellCtx *)ptr;

            DoNekppOperation(in, out, ctx, true);

            // Must return 0, otherwise PETSc complains.
            return 0;
        }

        /**
         * @brief Destroy matrix shell context object.
         *
         * Note the matrix shell and preconditioner share a common context so
         * this might have already been deallocated below, in which case we do
         * nothing.
         *
         * @param M  Matrix shell object
         */
        PetscErrorCode GlobalLinSysPETSc::DoDestroyMatCtx(Mat M)
        {
            void *ptr;
            MatShellGetContext(M, &ptr);
            ShellCtx *ctx = (ShellCtx *)ptr;
            delete ctx;
            return 0;
        }

        /**
         * @brief Destroy preconditioner context object.
         *
         * Note the matrix shell and preconditioner share a common context so
         * this might have already been deallocated above, in which case we do
         * nothing.
         *
         * @param pc  Preconditioner object
         */
        PetscErrorCode GlobalLinSysPETSc::DoDestroyPCCtx(PC pc)
        {
            void *ptr;
            PCShellGetContext(pc, &ptr);
            ShellCtx *ctx = (ShellCtx *)ptr;
            delete ctx;
            return 0;
        }

        /**
         * @brief Set up KSP solver object.
         *
         * This is reasonably generic setup -- most solver types can be changed
         * using the .petscrc file.
         *
         * @param tolerance  Residual tolerance to converge to.
         */
        void GlobalLinSysPETSc::SetUpSolver(NekDouble tolerance)
        {
            KSPCreate(PETSC_COMM_WORLD, &m_ksp);
            KSPSetTolerances(
                m_ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
            KSPSetFromOptions(m_ksp);
#if PETSC_VERSION_GE(3,5,0)
            KSPSetOperators(m_ksp, m_matrix, m_matrix);
#else
            KSPSetOperators(m_ksp, m_matrix, m_matrix, SAME_NONZERO_PATTERN);
#endif

            if (m_matMult == ePETScMatMultShell)
            {
                KSPSetPC(m_ksp, m_pc);
            }
        }
    }
}
