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
// Description: Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/PreconditionerDiagonal.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <math.h>

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

        string PreconditionerDiagonal::className1
                = GetPreconFactory().RegisterCreatorFunction(
                    "Null",
                    PreconditionerDiagonal::create,
                    "No Preconditioning");

        /**
         * @class Preconditioner
         *
         * This class implements preconditioning for the conjugate 
	 * gradient matrix solver.
	 */

         PreconditionerDiagonal::PreconditionerDiagonal(
                         const boost::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap)
           : Preconditioner(plinsys, pLocToGloMap),
	   m_linsys(plinsys),
           m_locToGloMap(pLocToGloMap),
           m_preconType(pLocToGloMap->GetPreconType())
         {
	 }

        void PreconditionerDiagonal::v_InitObject()
        {
	    if(m_preconType == MultiRegions::eDiagonal)
	    {
 	        GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
                if (solvertype == eIterativeFull)
                {
                    DiagonalPreconditionerSum();
                }
                else if(solvertype == eIterativeStaticCond ||
                      solvertype == eIterativeMultiLevelStaticCond)
                {
                    StaticCondDiagonalPreconditionerSum();
                }
                else
                {
                    ASSERTL0(0,"Unsupported solver type");
                }
	    }
	}


        /**
         * Populates preconditioner matrix with the identity i.e no
         * preconditioning.
         * @param   pLocToGloMap    Local to Global mapping.
         */
        void Preconditioner::NullPreconditioner()
        {
            int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
            int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt = nGlobal - nDir;
            MatrixStorage storage = eDIAGONAL;
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nInt, storage);
            DNekMat &M = (*m_preconditioner);

            for (unsigned int i = 0; i < nInt; ++i)
            {
                M.SetValue(i,i,1.0);
            }
        }


        /**
         * Diagonal preconditioner computed by summing the relevant elements of
         * the local matrix system.
         */
         void PreconditionerDiagonal::DiagonalPreconditionerSum()
         {
             boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();

             const StdRegions::StdExpansionVector &locExpVector = *(expList->GetExp());
             StdRegions::StdExpansionSharedPtr locExpansion;

             int i,j,n,cnt,gid1,gid2;
             NekDouble sign1,sign2,value;
             int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
             int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
             int nInt    = nGlobal - nDir;

             NekDouble zero = 0.0;

             // fill global matrix
             DNekScalMatSharedPtr loc_mat;
             Array<OneD, NekDouble> vOutput(nGlobal,0.0);
             MatrixStorage storage = eDIAGONAL;
             m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nInt, storage);
             DNekMat &M = (*m_preconditioner);

             int loc_lda;
             for(n = cnt = 0; n < expList->GetNumElmts(); ++n)
             {
                 loc_mat = (m_linsys.lock())->GetBlock(expList->GetOffset_Elmt_Id(n));
                 loc_lda = loc_mat->GetRows();

                 for(i = 0; i < loc_lda; ++i)
                 {
                     gid1 = m_locToGloMap->GetLocalToGlobalMap(cnt + i) - nDir;
                     sign1 =  m_locToGloMap->GetLocalToGlobalSign(cnt + i);
                     if(gid1 >= 0)
                     {
                         for(j = 0; j < loc_lda; ++j)
                         {
                             gid2 = m_locToGloMap->GetLocalToGlobalMap(cnt + j)
                                                                    - nDir;
                             sign2 = m_locToGloMap->GetLocalToGlobalSign(cnt + j);
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
                 cnt   += loc_lda;
             }

             // Assemble diagonal contributions across processes
             m_locToGloMap->UniversalAssemble(vOutput);

             // Populate preconditioner with reciprocal of diagonal elements
             for (unsigned int i = 0; i < nInt; ++i)
             {
                 M.SetValue(i,i,1.0/vOutput[i + nDir]);
             }
         }





        /**
         * Diagonal preconditioner defined as the inverse of the main
	 * diagonal of the Schur complement
	 *
         */
        void PreconditionerDiagonal::StaticCondDiagonalPreconditionerSum()
        {
            int nGlobalBnd = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int rows = nGlobalBnd - nDirBnd;

            MatrixStorage storage = eDIAGONAL;
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(rows, rows, storage);
            DNekMat &M = (*m_preconditioner);

            Array<OneD, NekDouble> vOutput(nGlobalBnd,0.0);

            // Extract diagonal contributions
            Array<OneD, NekDouble> diagonals = AssembleStaticCondGlobalDiagonals();
            for (unsigned int i = 0; i < rows; ++i)
            {
                vOutput[nDirBnd + i] = diagonals[i];
            }

            // Assemble diagonal contributions across processes
            m_locToGloMap->UniversalAssembleBnd(vOutput);

            // Populate preconditioner matrix
            for (unsigned int i = 0; i < rows; ++i)
            {
                M.SetValue(i,i,1.0/vOutput[nDirBnd + i]);
            }
        }

        /**
         *
         */
        void PreconditionerDiagonal::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            switch(m_preconType)
            {
            case MultiRegions::eDiagonal:
                 {
                     if (solvertype == eIterativeFull)
                     {
                         int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
                         int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                         int nNonDir = nGlobal-nDir;
                         DNekMat &M = (*m_preconditioner);

                         NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                         NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                         z = M * r;
		     }
                     else if(solvertype == eIterativeStaticCond)
                     {
                         int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                         int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
                         int nNonDir = nGlobal-nDir;
                         DNekMat &M = (*m_preconditioner);

                         NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                         NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                         z = M * r;
		     }
		     else
		     {
                         ASSERTL0(0,"Unsupported solver type");
		     }
		 }
                 break;
            case MultiRegions::eNull:
	         {
                     Vmath::Vcopy(pInput.num_elements(), pInput, 1, pOutput, 1);
		 }
		 break;
            default:
            ASSERTL0(0,"Unknown preconditioner");
            break;
	    }
	}
    }
}






