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

#include <map>

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/GlobalLinSysIterativeFull.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>
#include <MultiRegions/LocalToGlobalDGMap.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterativeCG
         *
         * This class implements a conjugate gradient matrix solver.
         * Preconditioning is implemented using a Jacobi (diagonal)
         * preconditioner.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysIterativeFull::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "IterativeFull",
                    GlobalLinSysIterativeFull::create,
                    "Iterative solver for full matrix system.");


        /**
         * Constructor for full direct matrix solve.
         * @param   pKey        Key specifying matrix to solve.
         * @param   pExp        Shared pointer to expansion list for applying
         *                      matrix evaluations.
         * @param   pLocToGloMap Local to global mapping.
         */
        GlobalLinSysIterativeFull::GlobalLinSysIterativeFull(
                    const GlobalLinSysKey &pKey,
                    const boost::shared_ptr<ExpList> &pExp,
                    const boost::shared_ptr<LocalToGlobalBaseMap> &pLocToGloMap)
                : GlobalLinSysIterative(pKey, pExp, pLocToGloMap)
        {
            ASSERTL1(m_linSysKey.GetGlobalSysSolnType()==eIterativeFull,
                     "This routine should only be used when using an Iterative "
                     "conjugate gradient matrix solve.");
        }


        /**
         *
         */
        GlobalLinSysIterativeFull::~GlobalLinSysIterativeFull()
        {

        }


        /**
         * Solve a global linear system with Dirichlet forcing using a
         * conjugate gradient method. This routine performs handling of the
         * Dirichlet forcing terms and wraps the underlying iterative solver
         * used for the remaining degrees of freedom.
         *
         * Consider solving for \f$x\f$, the matrix system \f$Ax=b\f$, where
         * \f$b\f$ is known. To enforce the Dirichlet terms we instead solve
         * \f[A(x-x_0) = b - Ax_0 \f]
         * where \f$x_0\f$ is the Dirichlet forcing.
         *
         * @param           pInput      RHS of linear system, \f$b\f$.
         * @param           pOutput     On input, values of dirichlet degrees
         *                              of freedom. On output, the solution
         *                              \f$x\f$.
         * @param           pLocToGloMap    Local to global mapping.
         * @param           pDirForcing Precalculated Dirichlet forcing.
         */
        void GlobalLinSysIterativeFull::v_Solve(
                    const Array<OneD, const NekDouble>  &pInput,
                          Array<OneD,       NekDouble>  &pOutput,
                    const LocalToGlobalBaseMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &pDirForcing)
        {
            bool vCG;
            if (m_locToGloMap
                = boost::dynamic_pointer_cast<LocalToGlobalC0ContMap>(
                                                                pLocToGloMap))
            {
                vCG = true;
            }
            else if (m_locToGloMap
                = boost::dynamic_pointer_cast<LocalToGlobalDGMap>(pLocToGloMap))
            {
                vCG = false;
            }

            bool dirForcCalculated = (bool) pDirForcing.num_elements();
            int nDirDofs  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobDofs = pLocToGloMap->GetNumGlobalCoeffs();
            int nLocDofs  = pLocToGloMap->GetNumLocalCoeffs();

            int nDirTotal = nDirDofs;
            m_expList->GetComm()->AllReduce(nDirTotal, LibUtilities::ReduceSum);

            if(nDirTotal)
            {
                // calculate the Dirichlet forcing
                Array<OneD, NekDouble> global_tmp(nGlobDofs);
                Array<OneD, NekDouble> offsetarray;

                if(dirForcCalculated)
                {
                    Vmath::Vsub(nGlobDofs, pInput.get(), 1,
                                pDirForcing.get(), 1,
                                global_tmp.get(), 1);
                }
                else
                {
                    // Calculate the dirichlet forcing B_b (== X_b) and
                    // substract it from the rhs
                    m_expList->GeneralMatrixOp(
                                    *m_linSysKey.GetGlobalMatrixKey(),
                                    pOutput, global_tmp, true);

                    Vmath::Vsub(nGlobDofs,  pInput.get(), 1,
                                            global_tmp.get(), 1,
                                            global_tmp.get(), 1);
                }
                if (vCG)
                {
                    SolveLinearSystem(nGlobDofs, global_tmp, pOutput, nDirDofs);
                }
                else
                {
                    ASSERTL0(false, "Need DG solve if using Dir BCs");
                }
            }
            else
            {
                SolveLinearSystem(nGlobDofs, pInput, pOutput);
            }
        }


        /**
         * Populates preconditioner with the identity to apply no
         * preconditioning.
         * @param   pLocToGloMap    Local to Global mapping.
         */
        void GlobalLinSysIterativeFull::ComputeNullPreconditioner(
                const boost::shared_ptr<LocalToGlobalBaseMap> &pLocToGloMap)
        {
            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
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
         * Diagonal preconditioner computed by evaluating the local matrix
         * acting on each basis vector (0,...,0,1,0,...,0). (deprecated)
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysIterativeFull::ComputeDiagonalPreconditioner(
                const boost::shared_ptr<LocalToGlobalBaseMap> &pLocToGloMap)
        {
            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nLocal  = pLocToGloMap->GetNumLocalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt = nGlobal - nDir;
            MatrixStorage storage = eDIAGONAL;
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nInt, storage);
            DNekMat &M = (*m_preconditioner);

            Array<OneD, int> vMap = pLocToGloMap->GetLocalToGlobalMap();

            for (unsigned int i = 0; i < nInt; ++i)
            {
                Array<OneD, NekDouble> test(nGlobal, 0.0);
                Array<OneD, NekDouble> test_local(nLocal, 0.0);
                test[i+nDir] = 1.0;
                m_expList->GeneralMatrixOp(*m_linSysKey.GetGlobalMatrixKey(),
                                test, test, true);

                M.SetValue(i,i,1.0/test[i+nDir]);
            }
        }


        /**
         * Diagonal preconditioner computed by summing the relevant elements of
         * the local matrix system.
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysIterativeFull::ComputeDiagonalPreconditionerSum(
                const boost::shared_ptr<LocalToGlobalBaseMap> &pLocToGloMap)
        {
            int i,j,n,cnt,gid1,gid2;
            NekDouble sign1,sign2,value;
            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt    = nGlobal - nDir;

            NekDouble zero = 0.0;

            // fill global matrix
            DNekScalMatSharedPtr loc_mat;
            Array<OneD, NekDouble> vOutput(nGlobal,0.0);
            MatrixStorage storage = eDIAGONAL;
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nInt, storage);
            DNekMat &M = (*m_preconditioner);

            int loc_lda;
            for(n = cnt = 0; n < m_expList->GetNumElmts(); ++n)
            {
                //loc_mat = vBlkMat->GetBlock(n,n);
                loc_mat = GetBlock(n);
                loc_lda = loc_mat->GetRows();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt + i) - nDir;
                    sign1 =  pLocToGloMap->GetLocalToGlobalSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = pLocToGloMap->GetLocalToGlobalMap(cnt + j)
                                                                    - nDir;
                            sign2 = pLocToGloMap->GetLocalToGlobalSign(cnt + j);
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
            pLocToGloMap->UniversalAssemble(vOutput);

            // Populate preconditioner with reciprocal of diagonal elements
            for (unsigned int i = 0; i < nInt; ++i)
            {
                M.SetValue(i,i,1.0/vOutput[i + nDir]);
            }
        }


        /**
         *
         */
        void GlobalLinSysIterativeFull::v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            // Perform matrix-vector operation A*d_i
            m_expList->GeneralMatrixOp(*m_linSysKey.GetGlobalMatrixKey(),
                                        pInput, pOutput, true);

            // retrieve robin boundary condition information and apply robin
            // boundary conditions to the solution.
            const std::map<int, RobinBCInfoSharedPtr> vRobinBCInfo
                                                = m_expList->GetRobinBCInfo();
            if(vRobinBCInfo.size() > 0)
            {
                ASSERTL0(false,
                        "Robin boundaries not set up in IterativeFull solver.");
                int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
                int nLocal  = m_locToGloMap->GetNumLocalCoeffs();
                int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                int nNonDir = nGlobal - nDir;
                Array<OneD, NekDouble> robin_A(nGlobal, 0.0);
                Array<OneD, NekDouble> robin_l(nLocal,  0.0);
                NekVector<NekDouble> robin(nNonDir,robin_A + nDir, eWrapper);

                // Operation: p_A = A * d_A
                // First map d_A to local solution
                m_locToGloMap->GlobalToLocal(pInput, robin_l);

                // Iterate over all the elements computing Robin BCs where
                // necessary
                for (int n = 0; n < m_expList->GetNumElmts(); ++n)
                {
                    int nel = m_expList->GetOffset_Elmt_Id(n);
                    int offset = m_expList->GetCoeff_Offset(n);
                    int ncoeffs = m_expList->GetExp(nel)->GetNcoeffs();

                    if(vRobinBCInfo.count(nel) != 0) // add robin mass matrix
                    {
                        RobinBCInfoSharedPtr rBC;
                        Array<OneD, NekDouble> tmp;
                        StdRegions::StdExpansionSharedPtr vExp = m_expList->GetExp(nel);

                        // add local matrix contribution
                        for(rBC = vRobinBCInfo.find(nel)->second;rBC; rBC = rBC->next)
                        {
                            vExp->AddRobinEdgeContribution(rBC->m_robinID,rBC->m_robinPrimitiveCoeffs, tmp = robin_l + offset);
                        }
                    }
                    else
                    {
                        Vmath::Zero(ncoeffs, &robin_l[offset], 1);
                    }
                }

                // Map local Robin contribution back to global coefficients
                m_locToGloMap->LocalToGlobal(robin_l, robin_A);
                // Add them to the output of the GeneralMatrixOp
                Vmath::Vadd(nGlobal, pOutput, 1, robin_A, 1, pOutput, 1);
            }

        }


        /**
         *
         */
        void GlobalLinSysIterativeFull::v_ComputePreconditioner()
        {
            ASSERTL1(!m_preconditioner.get(),
                     "Preconditioner has already been defined.");

            ComputeDiagonalPreconditionerSum(m_locToGloMap);
            m_map = m_locToGloMap->GetGlobalToUniversalMapUnique();
        }
    }
}
