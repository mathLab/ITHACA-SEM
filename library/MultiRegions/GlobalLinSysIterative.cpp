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
         * the standard PCG method (Demmel, Heath and Vorst, 1993)
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
            Array<OneD, NekDouble> w_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> s_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> p_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> r_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> q_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> tmp;

            // Create NekVector wrappers for linear algebra operations
            NekVector<NekDouble> in (nNonDir,pInput  + nDir,eWrapper);
            NekVector<NekDouble> out(nNonDir,pOutput + nDir,eWrapper);
            NekVector<NekDouble> w  (nNonDir,w_A + nDir,    eWrapper);
            NekVector<NekDouble> s  (nNonDir,s_A + nDir,    eWrapper);
            NekVector<NekDouble> p  (nNonDir,p_A,           eWrapper);
            NekVector<NekDouble> r  (nNonDir,r_A,           eWrapper);
            NekVector<NekDouble> q  (nNonDir,q_A,           eWrapper);

            int k;
            NekDouble alpha, beta, rho, rho_new, mu, eps, bb_inv, min_resid;
            Array<OneD, NekDouble> vExchange(3);

            // Initialise with zero as the initial guess.
            r = in;
            m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);
            v_DoMatrixMultiply(w_A, s_A);
            k = 0;

            vExchange[0] = Vmath::Dot2(nNonDir,
                                       r_A,
                                       w_A + nDir,
                                       m_map + nDir);

            vExchange[1] = Vmath::Dot2(nNonDir,
                                       s_A + nDir,
                                       w_A + nDir,
                                       m_map + nDir);

            vExchange[2] = Vmath::Dot2(nNonDir,
                                       r_A,
                                       r_A,
                                       m_map + nDir);

            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

            // If input vector is zero, set zero output and skip solve.
            if (vExchange[0] < NekConstants::kNekZeroTol)
            {
                Vmath::Zero(nNonDir, tmp = pOutput+nDir, 1);
                return;
            }

            rho       = vExchange[0];
            mu        = vExchange[1];
            beta      = 0.0;
            alpha     = rho/mu;
            eps       = 0.0;
            bb_inv    = 1.0/vExchange[2];
            min_resid = bb_inv;

            // Continue until convergence
            while (true)
            {
                ASSERTL0(k < 20000,
                         "Exceeded maximum number of iterations (20000)");

                ASSERTL0(eps*bb_inv <= 1.0 || k < 10,
                         "Conjugate gradient diverged. Tolerance too small?"
                         "Minimum residual achieved: "
                         + boost::lexical_cast<std::string>(sqrt(min_resid)));

                // Compute new search direction p_k, q_k
                p   = w   + beta  * p;
                q   = s   + beta  * q;

                // Update solution x_{k+1}
                out = out + alpha * p;

                // Update residual vector r_{k+1}
                r   = r   - alpha * q;

                // Apply preconditioner
                m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);

                // Perform the method-specific matrix-vector multiply operation.
                v_DoMatrixMultiply(w_A, s_A);

                // <r_{k+1}, w_{k+1}>
                vExchange[0] = Vmath::Dot2(nNonDir,
                                        r_A,
                                        w_A + nDir,
                                        m_map + nDir);
                // <s_{k+1}, w_{k+1}>
                vExchange[1] = Vmath::Dot2(nNonDir,
                                        s_A + nDir,
                                        w_A + nDir,
                                        m_map + nDir);

                // <r_{k+1}, r_{k+1}>
                vExchange[2] = Vmath::Dot2(nNonDir,
                                        r_A,
                                        r_A,
                                        m_map + nDir);

                // Perform inner-product exchanges
                vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

                rho_new = vExchange[0];
                mu      = vExchange[1];
                eps     = vExchange[2];

                // test if norm is within tolerance
                if (eps*bb_inv < m_tolerance * m_tolerance)
                {
                    break;
                }
                min_resid = min(min_resid, eps);

                // Compute search direction and solution coefficients
                beta  = rho_new/rho;
                alpha = rho_new/(mu - rho_new*beta/alpha);
                rho   = rho_new;

                k++;
            }
        }

#if 0
        /**
         * Solve a global linear system using the conjugate gradient method.
         * We solve only for the non-Dirichlet modes. The operator is evaluated
         * using the local-matrix representation. Distributed math routines are
         * used to support parallel execution of the solver.
         *
         * Standard CG algorithm.
         *
         * @param       pInput      Input vector of non-Dirichlet DOFs.
         * @param       pOutput     Solution vector of non-Dirichlet DOFs.
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
                m_precon = MemoryManager<Preconditioner>::AllocateSharedPtr(GetSharedThisPtr(),plocToGloMap);
            }

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm = m_expList.lock()->GetComm()->GetRowComm();

            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Allocate array storage
            Array<OneD, NekDouble> d_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> p_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> z_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> z_new_A(nNonDir, 0.0);
            Array<OneD, NekDouble> r_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> r_new_A(nNonDir, 0.0);

            // Create NekVector wrappers for linear algebra operations
            NekVector<NekDouble> in(nNonDir,pInput + nDir,eWrapper);
            NekVector<NekDouble> out(nNonDir,pOutput + nDir,eWrapper);
            NekVector<NekDouble> r(nNonDir,r_A,eWrapper);
            NekVector<NekDouble> r_new(nNonDir,r_new_A,eWrapper);
            NekVector<NekDouble> z(nNonDir,z_A,eWrapper);
            NekVector<NekDouble> z_new(nNonDir,z_new_A,eWrapper);
            NekVector<NekDouble> d(nNonDir,d_A + nDir, eWrapper);
            NekVector<NekDouble> p(nNonDir,p_A + nDir,eWrapper);

            int k;
            NekDouble alpha, beta, normsq, r_dot_z_old;
            Array<OneD, NekDouble> vExchange(2);

            // INVERSE of preconditioner matrix.
            //const DNekMat &M = (*m_preconditioner);

            // Initialise with zero as the initial guess.
            r = in;
            m_precon->DoPreconditioner(r_A,z_A);
            d = z;
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
            r_dot_z_old = vExchange[1];

            // Continue until convergence
            while (true)
            {
                // Perform the method-specific matrix-vector multiply operation.
                v_DoMatrixMultiply(d_A, p_A);

                // compute step length alpha
                // alpha denominator
                vExchange[0] = Vmath::Dot2(nNonDir,
                                        d_A + nDir,
                                        p_A + nDir,
                                        m_map + nDir);
                // perform exchange
                vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

                // compute alpha
                alpha = r_dot_z_old/vExchange[0];

                // approximate solution
                out   = out + alpha*d;

                // compute residual
                r_new = r   - alpha*p;

                // Apply preconditioner to new residual
                m_precon->DoPreconditioner(r_new_A,z_new_A);

                // beta
                vExchange[0] = Vmath::Dot2(nNonDir,
                                        r_new_A,
                                        z_new_A,
                                        m_map + nDir);

                // residual
                vExchange[1] = Vmath::Dot2(nNonDir,
                                        r_new_A,
                                        r_new_A,
                                        m_map + nDir);

                // perform exchange
                vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

                // extract values for beta and norm
                beta        = vExchange[0]/r_dot_z_old;
                r_dot_z_old = vExchange[0];
                normsq      = vExchange[1];

                // test if norm is within tolerance
                if (normsq < m_tolerance * m_tolerance)
                {
                    break;
                }

                // Compute new search direction
                d = z_new + beta*d;

                // Next step
                r = r_new;
                z = z_new;
                k++;

                ASSERTL0(k < 20000,
                         "Exceeded maximum number of iterations (20000)");
            }
        }
#endif
    }
}
