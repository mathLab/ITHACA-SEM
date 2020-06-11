///////////////////////////////////////////////////////////////////////////////
//
// File Preconditioner.cpp
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
// Description: Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/PreconditionerDiagonal.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <cmath>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * Registers the class with the Factory.
         */
        string PreconditionerDiagonal::className
                = GetPreconFactory().RegisterCreatorFunction(
                    "Diagonal",
                    PreconditionerDiagonal::create,
                    "Diagonal Preconditioning");
        /**
         * @class Preconditioner
         *
         * This class implements diagonal preconditioning for the conjugate
	 * gradient matrix solver.
	 */

         PreconditionerDiagonal::PreconditionerDiagonal(
                         const std::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap)
           : Preconditioner(plinsys, pLocToGloMap)
         {
	 }

        void PreconditionerDiagonal::v_InitObject()
        {
	}

        void PreconditionerDiagonal::v_BuildPreconditioner()
        {
            GlobalSysSolnType solvertype =
                m_locToGloMap.lock()->GetGlobalSysSolnType();
            if (solvertype == eIterativeFull)
            {
                DiagonalPreconditionerSum();
            }
            else if(solvertype == eIterativeStaticCond ||
                    solvertype == eIterativeMultiLevelStaticCond ||
                    solvertype == ePETScStaticCond ||
                    solvertype == ePETScMultiLevelStaticCond)
            {
                StaticCondDiagonalPreconditionerSum();
            }
            else
            {
                ASSERTL0(0,"Unsupported solver type");
            }
        }

        /**
         * Diagonal preconditioner computed by summing the relevant elements of
         * the local matrix system.
         */
         void PreconditionerDiagonal::DiagonalPreconditionerSum()
         {
             std::shared_ptr<MultiRegions::ExpList> expList =
                 ((m_linsys.lock())->GetLocMat()).lock();

             LocalRegions::ExpansionSharedPtr locExpansion;

             auto asmMap = m_locToGloMap.lock();

             int i,j,n,cnt,gid1,gid2;
             NekDouble sign1,sign2,value;
             int nGlobal = asmMap->GetNumGlobalCoeffs();
             int nDir    = asmMap->GetNumGlobalDirBndCoeffs();
             int nInt    = nGlobal - nDir;

             // fill global matrix
             DNekScalMatSharedPtr loc_mat;
             Array<OneD, NekDouble> vOutput(nGlobal,0.0);

             int loc_row;
             int nElmt = expList->GetNumElmts();
             for(n = cnt = 0; n < nElmt; ++n)
             {
                 loc_mat = (m_linsys.lock())->GetBlock(n);
                 loc_row = loc_mat->GetRows();

                 for(i = 0; i < loc_row; ++i)
                 {
                     gid1  = asmMap->GetLocalToGlobalMap(cnt+i)-nDir;
                     sign1 = asmMap->GetLocalToGlobalSign(cnt+i);

                     if(gid1 >= 0)
                     {
                         for(j = 0; j < loc_row; ++j)
                         {
                             gid2  = asmMap->GetLocalToGlobalMap(cnt+j)
                                 - nDir;
                             sign2 = asmMap->GetLocalToGlobalSign(cnt+j);
                             if(gid2 == gid1)
                             {
                                 // When global matrix is symmetric,
                                 // only add the value for the upper
                                 // triangular part in order to avoid
                                 // entries to be entered twice
                                 value = vOutput[gid1 + nDir]
                                            + sign1*sign2*(*loc_mat)(i,j);
                                 vOutput[gid1 + nDir] = value;
                             }
                         }
                     }
                 }
                 cnt   += loc_row;
             }

             // Assemble diagonal contributions across processes
             asmMap->UniversalAssemble(vOutput);

             m_diagonals = Array<OneD, NekDouble> (nInt);
             Vmath::Sdiv(nInt, 1.0, &vOutput[nDir], 1, &m_diagonals[0], 1);
         }

        /**
         * Diagonal preconditioner defined as the inverse of the main
	 * diagonal of the Schur complement
	 *
         */
        void PreconditionerDiagonal::StaticCondDiagonalPreconditionerSum()
        {
            auto asmMap = m_locToGloMap.lock();

            int nGlobalBnd = asmMap->GetNumGlobalBndCoeffs();
            int nDirBnd = asmMap->GetNumGlobalDirBndCoeffs();
            int rows = nGlobalBnd - nDirBnd;

            Array<OneD, NekDouble> vOutput(nGlobalBnd,0.0);

            // Extract diagonal contributions
            Array<OneD, NekDouble> diagonals = AssembleStaticCondGlobalDiagonals();
            for (unsigned int i = 0; i < rows; ++i)
            {
                vOutput[nDirBnd + i] = diagonals[i];
            }

            // Assemble diagonal contributions across processes
            asmMap->UniversalAssembleBnd(vOutput);

            m_diagonals = Array<OneD, NekDouble> (rows);
            Vmath::Sdiv(rows, 1.0, &vOutput[nDirBnd], 1, &m_diagonals[0], 1);
        }

        /**
         *
         */
        void PreconditionerDiagonal::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            auto asmMap = m_locToGloMap.lock();

            GlobalSysSolnType solvertype = asmMap->GetGlobalSysSolnType();

            int nGlobal = solvertype == eIterativeFull ?
                asmMap->GetNumGlobalCoeffs() :
                asmMap->GetNumGlobalBndCoeffs();
            int nDir    = asmMap->GetNumGlobalDirBndCoeffs();
            int nNonDir = nGlobal-nDir;
            Vmath::Vmul(nNonDir, &pInput[0], 1, &m_diagonals[0], 1, &pOutput[0], 1);
	}

        string PreconditionerNull::className
        = GetPreconFactory().RegisterCreatorFunction(
            "Null",
            PreconditionerNull::create,
            "No Preconditioning");

        /**
         * @class Null Preconditioner
         *
         * This class implements no preconditioning for the conjugate
	 * gradient matrix solver.
	 */
         PreconditionerNull::PreconditionerNull(
                         const std::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap)
           : Preconditioner(plinsys, pLocToGloMap)
         {
	 }

        /**
         *
         */
        void PreconditionerNull::v_InitObject()
        {
	}

        /**
         *
         */
        void PreconditionerNull::v_BuildPreconditioner()
        {
        }

        /**
         *
         */
        void PreconditionerNull::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            Vmath::Vcopy(pInput.size(), pInput, 1, pOutput, 1);
	}


    }
}






