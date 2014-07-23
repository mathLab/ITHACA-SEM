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

#include <MultiRegions/GlobalLinSysPETSc.h>

#include "petscis.h"

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysPETSc
         *
         * Solves a linear system using direct methods.
         */


        /// Constructor for full direct matrix solve.
        GlobalLinSysPETSc::GlobalLinSysPETSc(
            const GlobalLinSysKey                &pKey,
            const boost::weak_ptr<ExpList>       &pExp,
            const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
            : GlobalLinSys(pKey, pExp, pLocToGloMap)
        {
            // Initialise PETSc
            PetscInitialize(0, NULL, NULL, NULL);

            // Create matrix
            MatCreate(PETSC_COMM_WORLD, &m_matrix);
        }

        GlobalLinSysPETSc::~GlobalLinSysPETSc()
        {
        }

        void GlobalLinSysPETSc::v_SolveLinearSystem(
            const int                          pNumRows,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr        &locToGloMap,
            const int                          pNumDir)
        {
            const int nHomDofs = pNumRows - pNumDir;
            int i;
            VecSetValues(m_b, nHomDofs, &m_reorderedMap[0], &pInput[pNumDir], INSERT_VALUES);

            VecAssemblyBegin(m_b);
            VecAssemblyEnd  (m_b);

            PetscErrorCode ierr = KSPSolve(m_ksp, m_b, m_x);

            PetscInt its;
            KSPGetIterationNumber(m_ksp,&its);
            cout << "iteration = " << its << endl;

            IS isGlobal, isLocal;
            ISCreateGeneral(PETSC_COMM_SELF, nHomDofs, &m_reorderedMap[0], PETSC_COPY_VALUES, &isGlobal);
            ISCreateStride(PETSC_COMM_SELF, nHomDofs, 0, 1, &isLocal);

            Vec locVec;
            VecCreate        (PETSC_COMM_SELF, &locVec);
            VecSetSizes      (locVec, nHomDofs, PETSC_DECIDE);
            VecSetFromOptions(locVec);

            VecScatter ctx;
            VecScatterCreate(m_x, isGlobal, locVec, isLocal, &ctx);
            VecScatterBegin(ctx, m_x, locVec, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(ctx, m_x, locVec, INSERT_VALUES, SCATTER_FORWARD);

            PetscScalar *avec;
            VecGetArray(locVec, &avec);
            Vmath::Vcopy(nHomDofs, avec, 1, &pOutput[0], 1);
            VecRestoreArray(locVec, &avec);
        }
    }
}
