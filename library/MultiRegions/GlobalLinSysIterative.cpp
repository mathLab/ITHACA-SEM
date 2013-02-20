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

#include <MultiRegions/GlobalLinSysIterative.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterative
         *
         * Solves a linear system using iterative methods.
         */

        /// Constructor for full direct matrix solve.
        GlobalLinSysIterative::GlobalLinSysIterative(
                const GlobalLinSysKey &pKey,
                const boost::weak_ptr<ExpList> &pExpList,
                const boost::shared_ptr<AssemblyMap>
                                                       &pLocToGloMap)
                : GlobalLinSys(pKey, pExpList, pLocToGloMap),
                  m_totalIterations(0),
                  m_useProjection(false),
                  m_numPrevSols(0)
        {
            LibUtilities::SessionReaderSharedPtr vSession
                                            = pExpList.lock()->GetSession();
            vSession->LoadParameter("IterativeSolverTolerance",
                                    m_tolerance,
                                    NekConstants::kNekIterativeTol);

            LibUtilities::CommSharedPtr vComm = m_expList.lock()->GetComm()->GetRowComm();
            m_root = (vComm->GetRank())? false : true;

            m_verbose = (vSession->DefinesCmdLineArgument("verbose"))? true :false;

            std::string successiveRhs;
            vSession->LoadSolverInfo("SuccessiveRHS",  successiveRhs );
            try
            {
                int solutionsToStore = boost::lexical_cast<int>(successiveRhs);
                m_prevLinSol.set_capacity(solutionsToStore);
                m_useProjection = true;
                std::cout << "Using successive rhs projection with " << solutionsToStore << " solutions to be stored" << std::endl;
            }
            catch(...)
            {
                // lexical_cast will throw bad_lexical_cast if successiveRhs is not integer-valued
            }
        }

        GlobalLinSysIterative::~GlobalLinSysIterative()
        {
        }


        /**
         * 
         */
        void GlobalLinSysIterative::v_SolveLinearSystem(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &plocToGloMap,
                    const int nDir)
        {
            // Get vector sizes
            int nNonDir = nGlobal - nDir;
            Array<OneD, NekDouble> tmp;

            if (m_useProjection)
            {
                DoAconjugateProjection(nGlobal, pInput, pOutput, plocToGloMap, nDir);

                if (0)
                {
                    // check correctness: solve the same system with plain CG and compare
                    Array<OneD, NekDouble> cg_s   (nGlobal, 0.0);
                    NekVector<NekDouble>   cg     (nNonDir, tmp = cg_s + nDir, eWrapper);
                    NekVector<NekDouble>   x      (nNonDir, tmp = pOutput + nDir, eWrapper);

                    DoConjugateGradient(nGlobal, pInput, cg_s, plocToGloMap, nDir);

                    cg -= x;
                    NekDouble norm = CalculateAnorm(nGlobal, cg_s, nDir);
                    std::cout << "norm of solutions difference is = " << norm << std::endl;
                }
            }
            else
            {
                // applying plain Conjugate Gradient
                DoConjugateGradient(nGlobal, pInput, pOutput, plocToGloMap, nDir);
            }
        }


        /**
         * This method implements A-conjugate projection technique
         * in order to speed up successive linear solves with
         * right-hand sides arising from time-dependent discretisations.
         * (P.F.Fischer, Comput. Methods Appl. Mech. Engrg. 163, 1998)
         *
         */
        void GlobalLinSysIterative::DoAconjugateProjection(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &plocToGloMap,
                    const int nDir)
        {
            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                                = m_expList.lock()->GetComm()->GetRowComm();

            // Get vector sizes
            int nNonDir = nGlobal - nDir;
            Array<OneD, NekDouble> tmp;

            if (0 == m_numPrevSols)
            {
                // no previous solutions found, call CG

                DoConjugateGradient(nGlobal, pInput, pOutput, plocToGloMap, nDir);

                UpdateKnownSolutions(nGlobal, pOutput, nDir);
            }
            else
            {
                // Create NekVector wrappers for linear algebra operations
                NekVector<NekDouble> b     (nNonDir, pInput  + nDir, eWrapper);
                NekVector<NekDouble> x     (nNonDir, tmp = pOutput + nDir, eWrapper);

                // check the input vector (rhs) is not zero

                NekDouble rhsNorm = Vmath::Dot2(nNonDir,
                                                pInput + nDir,
                                                pInput + nDir,
                                                m_map + nDir);

                vComm->AllReduce(rhsNorm, Nektar::LibUtilities::ReduceSum);

                if (rhsNorm < NekConstants::kNekZeroTol)
                {
                    Array<OneD, NekDouble> tmp = pOutput+nDir;
                    Vmath::Zero(nNonDir, tmp, 1);
                    return;
                }

                // Allocate array storage
                Array<OneD, NekDouble> px_s       (nGlobal, 0.0);
                Array<OneD, NekDouble> pb_s       (nGlobal, 0.0);
                Array<OneD, NekDouble> tmpAx_s    (nGlobal, 0.0);
                Array<OneD, NekDouble> tmpx_s     (nGlobal, 0.0);

                NekVector<NekDouble> pb    (nNonDir, tmp = pb_s    + nDir, eWrapper);
                NekVector<NekDouble> px    (nNonDir, tmp = px_s    + nDir, eWrapper);
                NekVector<NekDouble> tmpAx (nNonDir, tmp = tmpAx_s + nDir, eWrapper);
                NekVector<NekDouble> tmpx  (nNonDir, tmp = tmpx_s  + nDir, eWrapper);


                // notation follows the paper cited:
                // \alpha_i = \tilda{x_i}^T b^n
                // projected x, px = \sum \alpha_i \tilda{x_i}

                Array<OneD, NekDouble> alpha     (m_prevLinSol.size(), 0.0);
                for (int i = 0; i < m_prevLinSol.size(); i++)
                {
                    alpha[i] = Vmath::Dot2(nNonDir,
                                           m_prevLinSol[i],
                                           pInput + nDir,
                                           m_map + nDir);
                }
                vComm->AllReduce(alpha, Nektar::LibUtilities::ReduceSum);

                for (int i = 0; i < m_prevLinSol.size(); i++)
                {
                    if (alpha[i] < NekConstants::kNekZeroTol)
                    {
                        continue;
                    }

                    NekVector<NekDouble> xi (nNonDir, m_prevLinSol[i], eWrapper);
                    px += alpha[i] * xi;
                }

                // pb = b^n - A px
                Vmath::Vcopy(nNonDir,
                             pInput.get() + nDir, 1,
                             pb_s.get()   + nDir, 1);

                v_DoMatrixMultiply(px_s, tmpAx_s);

                pb -= tmpAx;


                // solve the system with projected rhs
                DoConjugateGradient(nGlobal, pb_s, tmpx_s, plocToGloMap, nDir);


                // remainder solution + projection of previous solutions
                x = tmpx + px;

                // save the auxiliary solution to prev. known solutions
                UpdateKnownSolutions(nGlobal, tmpx_s, nDir);
            }
        }


        /**
         * Calculating A-norm of an input vector,
         * A-norm(x) := sqrt( < x, Ax > )
         */
        NekDouble GlobalLinSysIterative::CalculateAnorm(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &in,
                    const int nDir)
        {
            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                                = m_expList.lock()->GetComm()->GetRowComm();

            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Allocate array storage
            Array<OneD, NekDouble> tmpAx_s    (nGlobal, 0.0);

            v_DoMatrixMultiply(in, tmpAx_s);

            NekDouble anorm_sq = Vmath::Dot2(nNonDir,
                                             in      + nDir,
                                             tmpAx_s + nDir,
                                             m_map   + nDir);
            vComm->AllReduce(anorm_sq, Nektar::LibUtilities::ReduceSum);
            return std::sqrt(anorm_sq);
        }

        /**
         * Updates the storage of previously known solutions.
         * Performs normalisation of input vector wrt A-norm.
         */
        void GlobalLinSysIterative::UpdateKnownSolutions(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &newX,
                    const int nDir)
        {
            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                                = m_expList.lock()->GetComm()->GetRowComm();

            // Check the solution is non-zero
            NekDouble solNorm = Vmath::Dot2(nNonDir,
                                            newX + nDir,
                                            newX + nDir,
                                            m_map + nDir);
            vComm->AllReduce(solNorm, Nektar::LibUtilities::ReduceSum);

            if (solNorm < NekConstants::kNekZeroTol)
            {
                return;
            }


            // Allocate array storage
            Array<OneD, NekDouble> tmpAx_s    (nGlobal, 0.0);
            Array<OneD, NekDouble> px_s       (nGlobal, 0.0);
            Array<OneD, NekDouble> tmp1, tmp2;

            // Create NekVector wrappers for linear algebra operations
            NekVector<NekDouble> px           (nNonDir, tmp1 = px_s    + nDir, eWrapper);
            NekVector<NekDouble> tmpAx        (nNonDir, tmp2 = tmpAx_s + nDir, eWrapper);


            // calculating \tilda{x} - sum \alpha_i\tilda{x}_i

            Vmath::Vcopy(nNonDir,
                         tmp1 = newX + nDir, 1,
                         tmp2 = px_s + nDir, 1);

            if (m_prevLinSol.size() > 0)
            {
                v_DoMatrixMultiply(newX, tmpAx_s);
            }

            Array<OneD, NekDouble> alpha (m_prevLinSol.size(), 0.0);
            for (int i = 0; i < m_prevLinSol.size(); i++)
            {
                alpha[i] = Vmath::Dot2(nNonDir,
                                       m_prevLinSol[i],
                                       tmpAx_s + nDir,
                                       m_map + nDir);
            }
            vComm->AllReduce(alpha, Nektar::LibUtilities::ReduceSum);

            for (int i = 0; i < m_prevLinSol.size(); i++)
            {
                if (alpha[i] < NekConstants::kNekZeroTol)
                {
                    continue;
                }

                NekVector<NekDouble> xi (nNonDir, m_prevLinSol[i], eWrapper);
                px -= alpha[i] * xi;
            }


            // Some solutions generated by CG are identical zeros, see
            // solutions generated for Test_Tet_equitri.xml (IncNavierStokesSolver).
            // Not going to store identically zero solutions.

            NekDouble anorm = CalculateAnorm(nGlobal, px_s, nDir);
            if (anorm < NekConstants::kNekZeroTol)
            {
                return;
            }

            // normalisation of new solution
            Vmath::Smul(nNonDir, 1.0/anorm, px_s.get() + nDir, 1, px_s.get() + nDir, 1);

            // updating storage with non-Dirichlet-dof part of new solution vector
            m_prevLinSol.push_back(px_s + nDir);
            m_numPrevSols++;
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
        void GlobalLinSysIterative::DoConjugateGradient(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &plocToGloMap,
                    const int nDir)
        {
            // Check if preconditioner has been computed and compute if needed.
            if (!m_precon)
            {
                MultiRegions::PreconditionerType pType = plocToGloMap->GetPreconType();
                
                std::string PreconType = MultiRegions::PreconditionerTypeMap[pType];
                
                v_UniqueMap();
                m_precon = GetPreconFactory().CreateInstance(PreconType,GetSharedThisPtr(),plocToGloMap);
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
            NekVector<NekDouble> in (nNonDir,pInput  + nDir,      eWrapper);
            NekVector<NekDouble> out(nNonDir,tmp = pOutput + nDir,eWrapper);
            NekVector<NekDouble> w  (nNonDir,tmp = w_A + nDir,    eWrapper);
            NekVector<NekDouble> s  (nNonDir,tmp = s_A + nDir,    eWrapper);
            NekVector<NekDouble> p  (nNonDir,p_A,                 eWrapper);
            NekVector<NekDouble> r  (nNonDir,r_A,                 eWrapper);
            NekVector<NekDouble> q  (nNonDir,q_A,                 eWrapper);

            int k;
            NekDouble alpha, beta, rho, rho_new, mu, eps, bb_inv, min_resid;
            Array<OneD, NekDouble> vExchange(3);

            // Initialise with input initial guess.
            r = in;
            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            Vmath::Zero(nNonDir,tmp = pOutput + nDir,1);

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
            
            m_totalIterations = 0;
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
                ASSERTL0(k < 5000,
                         "Exceeded maximum number of iterations (5000)");

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
                    if (m_verbose && m_root)
                    {
                        cout << "CG iterations made = " << m_totalIterations 
                             << " using tolerance of "  << m_tolerance 
                             << " (eps = " << sqrt(eps) << ")" << endl;
                    }
                    break;
                }
                min_resid = min(min_resid, eps);

                // Compute search direction and solution coefficients
                beta  = rho_new/rho;
                alpha = rho_new/(mu - rho_new*beta/alpha);
                rho   = rho_new;

                k++;
                m_totalIterations++;
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
