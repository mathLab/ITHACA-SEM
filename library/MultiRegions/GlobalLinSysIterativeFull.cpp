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
#include <MultiRegions/GlobalLinSysIterativeFull.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

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
                    const boost::weak_ptr<ExpList> &pExp,
                    const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
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
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &pDirForcing)
        {
            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            bool vCG;
            if ((m_locToGloMap = boost::dynamic_pointer_cast<AssemblyMapCG>(
                     pLocToGloMap)))
            {
                vCG = true;
            }
            else if ((m_locToGloMap = boost::dynamic_pointer_cast<
                          AssemblyMapDG>(pLocToGloMap)))
            {
                vCG = false;
            }

            bool dirForcCalculated = (bool) pDirForcing.num_elements();
            int nDirDofs  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobDofs = pLocToGloMap->GetNumGlobalCoeffs();
            int nDirTotal = nDirDofs;
            
            expList->GetComm()->AllReduce(nDirTotal, LibUtilities::ReduceSum);
            
            Array<OneD, NekDouble> tmp(nGlobDofs);

            if(nDirTotal)
            {
                // calculate the Dirichlet forcing
                if(dirForcCalculated)
                {
                    Vmath::Vsub(nGlobDofs, pInput.get(), 1,
                                pDirForcing.get(), 1,
                                tmp.get(), 1);
                }
                else
                {
                    // Calculate the dirichlet forcing B_b (== X_b) and
                    // substract it from the rhs
                    expList->GeneralMatrixOp(
                        m_linSysKey, pOutput, tmp, eGlobal);

                    Vmath::Vsub(nGlobDofs, pInput.get(), 1,
                                           tmp.get(),    1,
                                           tmp.get(),    1);
                }
                if (vCG)
                {
                    SolveLinearSystem(
                        nGlobDofs, tmp, pOutput, pLocToGloMap, nDirDofs);
                }
                else
                {
                    ASSERTL0(false, "Need DG solve if using Dir BCs");
                }
            }
            else
            {
                Vmath::Vcopy(nGlobDofs, pInput, 1, tmp, 1);
                SolveLinearSystem(nGlobDofs, tmp, pOutput, pLocToGloMap);
            }
        }


        /**
         *
         */
        void GlobalLinSysIterativeFull::v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            // Perform matrix-vector operation A*d_i
            expList->GeneralMatrixOp(m_linSysKey,
                                     pInput, pOutput, eGlobal);

            // retrieve robin boundary condition information and apply robin
            // boundary conditions to the solution.
            const std::map<int, RobinBCInfoSharedPtr> vRobinBCInfo
                                                = expList->GetRobinBCInfo();
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
                Array<OneD, NekDouble> tmp;
                NekVector<NekDouble> robin(nNonDir,
                                           tmp = robin_A + nDir, eWrapper);

                // Operation: p_A = A * d_A
                // First map d_A to local solution
                m_locToGloMap->GlobalToLocal(pInput, robin_l);

                // Iterate over all the elements computing Robin BCs where
                // necessary
                for (int n = 0; n < expList->GetNumElmts(); ++n)
                {
                    int nel = expList->GetOffset_Elmt_Id(n);
                    int offset = expList->GetCoeff_Offset(n);
                    int ncoeffs = expList->GetExp(nel)->GetNcoeffs();

                    if(vRobinBCInfo.count(nel) != 0) // add robin mass matrix
                    {
                        RobinBCInfoSharedPtr rBC;
                        Array<OneD, NekDouble> tmp;
                        StdRegions::StdExpansionSharedPtr vExp = expList->GetExp(nel);

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
        void GlobalLinSysIterativeFull::v_UniqueMap()
        {
            m_map = m_locToGloMap->GetGlobalToUniversalMapUnique();
        }

    }
}
