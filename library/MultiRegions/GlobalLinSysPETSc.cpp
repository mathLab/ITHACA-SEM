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
            int i;
            for (i = pNumDir; i < pNumRows; ++i)
            {
                VecSetValue(m_b, i-pNumDir, pInput[i], INSERT_VALUES);
            }

            PetscErrorCode ierr = KSPSolve(m_ksp, m_b, m_x);
            PetscScalar   *avec;
            VecGetArray(m_x, &avec);

            PetscInt its;
            KSPGetIterationNumber(m_ksp,&its);
            cout << "iteration = " << its << endl;

            for (i = 0; i < pNumRows - pNumDir; ++i)
            {
                pOutput[i] = avec[i];
            }
        }
    }
}
