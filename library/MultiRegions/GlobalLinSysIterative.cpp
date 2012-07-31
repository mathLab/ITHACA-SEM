///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterative.cpp
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
// Description: GlobalLinSysIterative definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/GlobalLinSysIterative.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterative
         *
         * Solves a linear system using direct methods.
         */

        /// Constructor for full direct matrix solve.
        GlobalLinSysIterative::GlobalLinSysIterative(
                const GlobalLinSysKey &pKey,
                const boost::weak_ptr<ExpList> &pExpList,
                const boost::shared_ptr<AssemblyMap>
                                                       &pLocToGloMap)
                : GlobalLinSys(pKey, pExpList, pLocToGloMap)
        {
            LibUtilities::SessionReaderSharedPtr vSession
                                            = pExpList.lock()->GetSession();
            vSession->LoadParameter("IterativeSolverTolerance",
                                    m_tolerance,
                                    NekConstants::kNekIterativeTol);

        }

        GlobalLinSysIterative::~GlobalLinSysIterative()
        {
        }

        /**
         * Solve a global linear system using the conjugate gradient method.
         * We solve only for the non-Dirichlet modes. The operator is evaluated
         * using an auxiliary function v_DoMatrixMultiply defined by the
         * specific solver. Distributed math routines are used to support
         * parallel execution of the solver.
         *
         * The implemented algorithm uses a reduced-communication reordering of
         * the standard PCG method (Amin, Sadayappan, Gudavalli, IEEE 1994)
         *
         * @param       pInput      Input vector of all DOFs.
         * @param       pOutput     Solution vector of all DOFs.
         */
        void GlobalLinSysIterative::v_SolveLinearSystem(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &plocToGloMap,
                    const int nDir)
        {
            // Check if preconditioner has been computed and compute if needed.
            if (!m_precon)
            {
                v_UniqueMap();
                m_precon = MemoryManager<Preconditioner>::AllocateSharedPtr(
                                            GetSharedThisPtr(),plocToGloMap);
            }

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                                = m_expList.lock()->GetComm()->GetRowComm();

            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Allocate array storage
            Array<OneD, NekDouble> p_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> q_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> y_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> r_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> z_A    (nNonDir, 0.0);

            // Create NekVector wrappers for linear algebra operations
            NekVector<NekDouble> in (nNonDir,pInput + nDir, eWrapper);
            NekVector<NekDouble> out(nNonDir,pOutput + nDir,eWrapper);
            NekVector<NekDouble> p  (nNonDir,p_A + nDir,    eWrapper);
            NekVector<NekDouble> q  (nNonDir,q_A + nDir,    eWrapper);
            NekVector<NekDouble> y  (nNonDir,y_A,           eWrapper);
            NekVector<NekDouble> r  (nNonDir,r_A,           eWrapper);
            NekVector<NekDouble> z  (nNonDir,z_A,           eWrapper);

            int k;
            NekDouble alpha, beta, r_dot_d, b_dot_b, min_resid;
            Array<OneD, NekDouble> vExchange(3);

            // Initialise with zero as the initial guess.
            r = in;
            m_precon->DoPreconditioner(r_A,z_A);
            p = z;
            k = 0;

            vExchange[0] = Vmath::Dot2(nNonDir, r_A, r_A, m_map + nDir);
            vExchange[1] = Vmath::Dot2(nNonDir, r_A, z_A, m_map + nDir);
            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

            // If input vector is zero, set zero output and skip solve.
            if (vExchange[0] < NekConstants::kNekZeroTol)
            {
                Vmath::Zero(nGlobal, pOutput, 1);
                return;
            }
            b_dot_b = vExchange[0];
            r_dot_d = vExchange[1];
            min_resid = b_dot_b;

            // Continue until convergence
            while (true)
            {
                ASSERTL0(k < 20000,
                         "Exceeded maximum number of iterations (20000)");

                ASSERTL0(vExchange[0] <= b_dot_b || k < 10,
                         "Conjugate gradient diverged. Tolerance too small?"
                         "Minimum residual achieved: "
                         + boost::lexical_cast<std::string>(sqrt(min_resid)));

                // Perform the method-specific matrix-vector multiply operation.
                v_DoMatrixMultiply(p_A, q_A);

                // Apply preconditioner
                m_precon->DoPreconditioner(q_A + nDir, y_A);

                // <r_k, r_k>
                vExchange[0] = Vmath::Dot2(nNonDir,
                                        r_A,
                                        r_A,
                                        m_map + nDir);
                // <p_k, q_k>
                vExchange[1] = Vmath::Dot2(nNonDir,
                                        p_A + nDir,
                                        q_A + nDir,
                                        m_map + nDir);

                // <q_k, y_k>
                vExchange[2] = Vmath::Dot2(nNonDir,
                                        q_A + nDir,
                                        y_A,
                                        m_map + nDir);

                // Perform exchange
                vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

                min_resid = min(min_resid, vExchange[0]);

                // test if norm is within tolerance
                if (vExchange[0] < m_tolerance * m_tolerance)
                {
                    break;
                }

                // compute alpha
                alpha   =  r_dot_d / vExchange[1];
                beta    =  alpha * vExchange[2] / vExchange[1] - 1.0;
                r_dot_d *= beta;

                // Update residual vector
                r = r   - alpha * q;
                z = z   - alpha * y;

                // Update solution
                out   = out + alpha * p;

                // Compute new search direction
                p = z + beta * p;

                k++;
            }
        }

    }
}
