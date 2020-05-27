///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.h
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
// Description: GlobalLinSysPETSc header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSPETSC_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSPETSC_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSys.h>

#include <petscmat.h>
#include <petscksp.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;

        /// Enumerator
        enum PETScMatMult
        {
            ePETScMatMultSparse,
            ePETScMatMultShell
        };

        /// A PETSc global linear system.
        class GlobalLinSysPETSc : virtual public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysPETSc(
                const GlobalLinSysKey                &pKey,
                const std::weak_ptr<ExpList>         &pExp,
                const std::shared_ptr<AssemblyMap>   &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysPETSc();

            virtual void v_SolveLinearSystem(
                const int                          pNumRows,
                const Array<OneD,const NekDouble> &pInput,
                      Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr        &locToGloMap,
                const int                          pNumDir);

        protected:
            /// PETSc matrix object.
            Mat               m_matrix;
            /// PETSc vector objects used for local storage.
            Vec               m_x, m_b, m_locVec;
            /// KSP object that represents solver system.
            KSP               m_ksp;
            /// PCShell for preconditioner.
            PC                m_pc;
            /// Enumerator to select matrix multiplication type.
            PETScMatMult      m_matMult;
            /// Reordering that takes universal IDs to a unique row in the PETSc
            /// matrix. @see GlobalLinSysPETSc::CalculateReordering
            std::vector<int>  m_reorderedMap;
            /// PETSc scatter context that takes us between Nektar++ global
            /// ordering and PETSc vector ordering.
            VecScatter        m_ctx;
            /// Number of unique degrees of freedom on this process.
            int               m_nLocal;

            PreconditionerSharedPtr m_precon;

            /**
             * @brief Internal struct for MatShell and PCShell calls to store
             * current context for callback.
             *
             * To use the MatShell/PCShell representation inside PETSc KSP and
             * PC objects (so that we can use the local spectral element
             * approach) requires the use of a callback function, which must be
             * static. This is a lightweight wrapper allowing us to call a
             * virtual function so that we can handle the static
             * condensation/full variants of the global system.
             *
             * @see GlobalLinSysPETSc::DoMatrixMultiply
             */
            struct ShellCtx
            {
                /// Number of global degrees of freedom.
                int nGlobal;
                /// Number of Dirichlet degrees of freedom.
                int nDir;
                /// Pointer to the original calling object.
                GlobalLinSysPETSc *linSys;
            };

            void SetUpScatter();
            void SetUpMatVec(int nGlobal, int nDir);
            void SetUpSolver(NekDouble tolerance);
            void CalculateReordering(
                const Array<OneD, const int> &glo2uniMap,
                const Array<OneD, const int> &glo2unique,
                const AssemblyMapSharedPtr   &pLocToGloMap);

            virtual void v_DoMatrixMultiply(
                const Array<OneD, const NekDouble>& pInput,
                      Array<OneD,       NekDouble>& pOutput) = 0;
        private:
            static std::string matMult;
            static std::string matMultIds[];

            static PetscErrorCode DoMatrixMultiply(Mat M, Vec in, Vec out);
            static PetscErrorCode DoPreconditioner(PC pc, Vec in, Vec out);
            static void DoNekppOperation(
                Vec &in, Vec &out, ShellCtx *ctx, bool precon);
            static PetscErrorCode DoDestroyMatCtx(Mat M);
            static PetscErrorCode DoDestroyPCCtx (PC pc);
        };
    }
}

#endif
