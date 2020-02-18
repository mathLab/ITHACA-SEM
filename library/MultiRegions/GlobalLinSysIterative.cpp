///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterative.cpp
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
// Description: GlobalLinSysIterative definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterative.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        std::string GlobalLinSysIterative::lookupIds[2] =
        {
            LibUtilities::SessionReader::RegisterEnumValue(
                "IterativeMethod", "ConjugateGradient", eConjugateGradient),
        };
        std::string GlobalLinSysIterative::def =
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                "IterativeMethod", "ConjugateGradient");

        /**
         * @class GlobalLinSysIterative
         *
         * Solves a linear system using iterative methods.
         */

        /// Constructor for full direct matrix solve.
        GlobalLinSysIterative::GlobalLinSysIterative(
            const GlobalLinSysKey &pKey,
            const std::weak_ptr<ExpList> &pExpList,
            const std::shared_ptr<AssemblyMap>
            &pLocToGloMap)
            : GlobalLinSys(pKey, pExpList, pLocToGloMap),
              m_rhs_magnitude(NekConstants::kNekUnsetDouble),
              m_prec_factor(NekConstants::kNekUnsetDouble),
              m_rhs_mag_sm(0.9),
              m_precon(NullPreconditionerSharedPtr),
              m_totalIterations(0),
              m_useProjection(false),
              m_numPrevSols(0)
        {
            m_tolerance = pLocToGloMap->GetIterativeTolerance();
            m_maxiter   = pLocToGloMap->GetMaxIterations();
           
            m_maxstorage = pLocToGloMap->GetMaxStorage();
            m_maxhesband= pLocToGloMap->GetMaxHesband();

            LibUtilities::CommSharedPtr vComm = m_expList.lock()->GetComm()->GetRowComm();
            m_root    = (vComm->GetRank()) ? false : true;

            int successiveRHS;

            if((successiveRHS = pLocToGloMap->GetSuccessiveRHS()))
            {
                m_prevLinSol.set_capacity(successiveRHS);
                m_useProjection = true;
            }
            else
            {
                m_useProjection = false;
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
            const Array<OneD, const NekDouble> &pInput,
            Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr &plocToGloMap,
            const int nDir)
        {
            IterativeMethodType pType = plocToGloMap->GetIteraterType();
            switch(pType)
            {
            case eConjugateGradient:
                if (m_useProjection)
                {
                    DoAconjugateProjection(nGlobal, pInput, pOutput, plocToGloMap, nDir);
                }
                else
                {
                    // applying plain Conjugate Gradient
                    DoConjugateGradient(nGlobal, pInput, pOutput, plocToGloMap, nDir);
                }
                break;
            default:
                ASSERTL0(false, "IterativeMethodType NOT CORRECT.");
                break;
            }
        }


        /**
         * This method implements A-conjugate projection technique
         * in order to speed up successive linear solves with
         * right-hand sides arising from time-dependent discretisations.
         * (P.F.Fischer, Comput. Methods Appl. Mech. Engrg. 163, 1998)
         */
        void GlobalLinSysIterative::DoAconjugateProjection(
            const int nGlobal,
            const Array<OneD, const NekDouble> &pInput,
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
                    Array<OneD, NekDouble> tmp = pOutput + nDir;
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
            const Array<OneD, const NekDouble> &in,
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
            const Array<OneD, const NekDouble> &newX,
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
            Vmath::Smul(nNonDir, 1.0 / anorm, px_s.get() + nDir, 1, px_s.get() + nDir, 1);

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
         * @param       pInput      Input residual  of all DOFs.  
         * @param       pOutput     Solution vector of all DOFs.  
         */
        void GlobalLinSysIterative::DoConjugateGradient(
            const int                          nGlobal,
            const Array<OneD, const NekDouble> &pInput,
            Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr        &plocToGloMap,
            const int                          nDir)
        {
            if (!m_precon)
            {
                v_UniqueMap();
                m_precon = CreatePrecon(plocToGloMap);
                m_precon->BuildPreconditioner();
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
            NekVector<NekDouble> in (nNonDir, pInput  + nDir,      eWrapper);
            NekVector<NekDouble> out(nNonDir, tmp = pOutput + nDir, eWrapper);
            NekVector<NekDouble> w  (nNonDir, tmp = w_A + nDir,    eWrapper);
            NekVector<NekDouble> s  (nNonDir, tmp = s_A + nDir,    eWrapper);
            NekVector<NekDouble> p  (nNonDir, p_A,                 eWrapper);
            NekVector<NekDouble> r  (nNonDir, r_A,                 eWrapper);
            NekVector<NekDouble> q  (nNonDir, q_A,                 eWrapper);

            int k;
            NekDouble alpha, beta, rho, rho_new, mu, eps,  min_resid;
            Array<OneD, NekDouble> vExchange(3, 0.0);

            // Copy initial residual from input
            r = in;
            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            Vmath::Zero(nNonDir, tmp = pOutput + nDir, 1);


            // evaluate initial residual error for exit check
            vExchange[2] = Vmath::Dot2(nNonDir,
                                       r_A,
                                       r_A,
                                       m_map + nDir);

            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

            eps       = vExchange[2];

            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                NekVector<NekDouble> inGlob (nGlobal, pInput, eWrapper);
                Set_Rhs_Magnitude(inGlob);
            }

            m_totalIterations = 0;

            // If input residual is less than tolerance skip solve.
            if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
            {
                if (m_verbose && m_root)
                {
                    cout << "CG iterations made = " << m_totalIterations 
                         << " using tolerance of "  << m_tolerance 
                         << " (error = " << sqrt(eps/m_rhs_magnitude) 
                         << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")" 
                         << endl;
                }
                return;
            }

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

            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

            rho               = vExchange[0];
            mu                = vExchange[1];
            min_resid         = m_rhs_magnitude;
            beta              = 0.0;
            alpha             = rho / mu;
            m_totalIterations = 1;

            // Continue until convergence
            while (true)
            {
                if(k >= m_maxiter)
                {
                    if (m_root)
                    {
                        cout << "CG iterations made = " << m_totalIterations
                             << " using tolerance of "  << m_tolerance
                             << " (error = " << sqrt(eps / m_rhs_magnitude)
                             << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                             << endl;
                    }
                    ROOTONLY_NEKERROR(ErrorUtil::efatal,
                                      "Exceeded maximum number of iterations");
                }

                // Compute new search direction p_k, q_k
                Vmath::Svtvp(nNonDir, beta, &p_A[0], 1, &w_A[nDir], 1, &p_A[0], 1);
                Vmath::Svtvp(nNonDir, beta, &q_A[0], 1, &s_A[nDir], 1, &q_A[0], 1);

                // Update solution x_{k+1}
                Vmath::Svtvp(nNonDir, alpha, &p_A[0], 1, &pOutput[nDir], 1, &pOutput[nDir], 1);

                // Update residual vector r_{k+1}
                Vmath::Svtvp(nNonDir, -alpha, &q_A[0], 1, &r_A[0], 1, &r_A[0], 1);

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

                m_totalIterations++;

                // test if norm is within tolerance
                if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
                {
                    if (m_verbose && m_root)
                    {
                        cout << "CG iterations made = " << m_totalIterations
                             << " using tolerance of "  << m_tolerance
                             << " (error = " << sqrt(eps / m_rhs_magnitude)
                             << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                             << endl;
                    }
                    break;
                }
                min_resid = min(min_resid, eps);

                // Compute search direction and solution coefficients
                beta  = rho_new / rho;
                alpha = rho_new / (mu - rho_new * beta / alpha);
                rho   = rho_new;
                k++;
            }
        }

        // Arnoldi Subroutine
        void GlobalLinSysIterative::DoArnoldi(
            const int starttem,
            const int endtem,
            const int nGlobal,
            const int nDir,
            // V_total(:,1:nd)
            Array<OneD, Array<OneD,  NekDouble> > &V_local,
            // V[nd]
            Array<OneD, NekDouble> &Vsingle1,
            // V[nd+1]
            Array<OneD, NekDouble> &Vsingle2,
            //h
            Array<OneD, NekDouble> &hsingle)
        {
            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();

            // To notice, V_local's order not certainly equal to starttem:endtem
            // starttem:endtem is the entry position in Hessenburg matrix
            NekDouble alpha, beta;
            Array<OneD, NekDouble> tmp1, tmp2;
            int numbertem;
            int nNonDir = nGlobal - nDir;
            //later for parallel
            NekDouble vExchange = 0.0;
            // w=AV(:,nd)
            Array<OneD, NekDouble> w(nGlobal, 0.0);

            v_DoMatrixMultiply(Vsingle1, w);

            tmp1 = w + nDir;
            tmp2 = w + nDir;
            m_precon->DoPreconditioner(tmp1, tmp2);

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Modified Gram-Schmidt
            // The pointer not certainly equal to starttem.
            // Like initially, Gmres-deep need to use numbertem=0
            numbertem = starttem;
            for(int i = starttem; i < endtem; ++i)
            {
                vExchange = Vmath::Dot2(nNonDir,
                                        &w[0] + nDir,
                                        &V_local[numbertem][0] + nDir,
                                        &m_map[0] + nDir);

                vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

                hsingle[i] = vExchange;

                beta = -1.0 * vExchange;
                Vmath::Svtvp(nNonDir, beta, &V_local[numbertem][0] + nDir, 1, &w[0] + nDir, 1, &w[0] + nDir, 1);

                numbertem = numbertem + 1;
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Classical Gram-Schmidt
            // int number = endtem - starttem;
            // Array<OneD,NekDouble> vExchangeArray(number, 0.0);
            // numbertem=0;
            // for (int i = starttem; i < endtem; ++i)
            // {
            //     vExchangeArray[numbertem] =
            //         Vmath::Dot2(nNonDir, &w[0] + nDir, &V_local[i][0] + nDir,
            //                     &m_map[0] + nDir);
            //     numbertem = numbertem + 1;
            // }

            // vComm->AllReduce(vExchangeArray, Nektar::LibUtilities::ReduceSum);

            // numbertem = 0;
            // for (int i = starttem; i < endtem; ++i)
            // {
            //     hsingle[i] = vExchangeArray[numbertem];
            //     beta       = -1.0 * vExchangeArray[numbertem];
            //     Vmath::Svtvp(nNonDir, beta, &V_local[i][0] + nDir, 1,
            //                  &w[0] + nDir, 1, &w[0] + nDir, 1);
            //     numbertem = numbertem + 1;
            // }
            // end of Classical Gram-Schmidt
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // calculate the L2 norm and normalize
            vExchange = Vmath::Dot2(nNonDir,
                                    &w[0] + nDir,
                                    &w[0] + nDir,
                                    &m_map[0] + nDir);

            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

            hsingle[endtem] = sqrt(vExchange);

            alpha = 1.0 / hsingle[endtem];
            Vmath::Smul(nNonDir, alpha, &w[0] + nDir, 1, &Vsingle2[0] + nDir, 1);
        }

        // QR factorization through Givens rotation
        void GlobalLinSysIterative::DoGivensRotation(
            const int starttem,
            const int endtem,
            const int nGlobal,
            const int nDir,
            Array<OneD, NekDouble> &c,
            Array<OneD, NekDouble> &s,
            Array<OneD, NekDouble> &hsingle,
            Array<OneD, NekDouble> &eta)
        {
            boost::ignore_unused(nGlobal,nDir);
            NekDouble temp_dbl, dd, hh;
            int idtem = endtem - 1;
            // The starttem and endtem are beginning and ending order of Givens rotation
            // They usually equal to the beginning position and ending position of Hessenburg matrix
            // But sometimes starttem will change, like if it is initial 0 and becomes nonzero because previous Givens rotation
            // See Yu Pan's User Guide
            for(int i = starttem; i < idtem; ++i)
            {
                temp_dbl = c[i] * hsingle[i] - s[i] * hsingle[i + 1];
                hsingle[i + 1] = s[i] * hsingle[i] + c[i] * hsingle[i + 1];
                hsingle[i] = temp_dbl;
            }
            dd = hsingle[idtem];
            hh = hsingle[endtem];
            if(hh == 0.0)
            {
                c[idtem] = 1.0;
                s[idtem] = 0.0;
            }
            else if (abs(hh) > abs(dd))
            {
                temp_dbl = -dd / hh;
                s[idtem] = 1.0 / sqrt(1.0 + temp_dbl * temp_dbl);
                c[idtem] = temp_dbl * s[idtem];
            }
            else
            {
                temp_dbl = -hh / dd;
                c[idtem] = 1.0 / sqrt(1.0 + temp_dbl * temp_dbl);
                s[idtem] = temp_dbl * c[idtem];
            }

            hsingle[idtem] = c[idtem] * hsingle[idtem] - s[idtem] * hsingle[endtem];
            hsingle[endtem] = 0.0;

            temp_dbl = c[idtem] * eta[idtem] - s[idtem] * eta[endtem];
            eta[endtem]       = s[idtem] * eta[idtem] + c[idtem] * eta[endtem];
            eta[idtem] = temp_dbl;
        }

        // Backward calculation
        // to notice, Hesssenburg matrix's column and row changes due to use Array<OneD,Array<OneD,NekDouble>>
        void GlobalLinSysIterative::DoBackward(
            const int  number,
            Array<OneD, Array<OneD, NekDouble> > &A,
            const Array<OneD, const NekDouble> &b,
            Array <OneD, NekDouble> &y)
        {
            // number is the entry number, but C++'s order need to be one smaller
            int maxid = number - 1;
            NekDouble sum;
            y[maxid] = b[maxid] / A[maxid][maxid];

            for (int i = maxid - 1; i > -1; --i)
            {
                sum = b[i];

                for (int j = i + 1; j < number; ++j)
                {
                    // i and j changes due to use Array<OneD,Array<OneD,NekDouble>>
                    sum = sum - y[j] * A[j][i];
                }
                y[i] = sum / A[i][i];
            }
        }

        NekDouble GlobalLinSysIterative::DoGmresRestart(
            const bool                         restarted,
            const bool                         truncted,
            const int                          nGlobal,
            const Array<OneD, const NekDouble> &pInput,
            Array<OneD,      NekDouble> &pOutput,
            const int                          nDir)
        {
            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();

            int nNonDir = nGlobal - nDir;

            // Allocate array storage of coefficients
            // Hessenburg matrix
            Array<OneD, Array<OneD, NekDouble> > hes    (m_maxstorage);
            for (int i = 0; i < m_maxstorage; ++i)
            {
                hes[i] = Array<OneD, NekDouble>(m_maxstorage + 1, 0.0);
            }
            // Hesseburg matrix after rotation
            Array<OneD, Array<OneD, NekDouble> >  Upper  (m_maxstorage);
            for (int i = 0; i < m_maxstorage; ++i)
            {
                Upper[i] = Array<OneD, NekDouble>(m_maxstorage + 1, 0.0);
            }
            // Total search directions
            Array<OneD, Array<OneD, NekDouble> >  V_total(m_maxstorage + 1);
            for (int i = 0; i < m_maxstorage + 1; ++i)
            {
                V_total[i] = Array<OneD, NekDouble>(nGlobal, 0.0);
            }
            //Residual
            Array<OneD, NekDouble> eta    (m_maxstorage + 1, 0.0);
            //Givens rotation c
            Array<OneD, NekDouble> cs     (m_maxstorage, 0.0);
            //Givens rotation s
            Array<OneD, NekDouble> sn     (m_maxstorage, 0.0);
            //Total coefficients, just for check
            Array<OneD, NekDouble> y_total     (m_maxstorage, 0.0);
            // Residual
            NekDouble eps;
            //Search direction order
            Array<OneD, int> id (m_maxstorage, 0);
            Array<OneD, int> id_start (m_maxstorage, 0);
            Array<OneD, int> id_end (m_maxstorage, 0);
            // temporary variables
            int idtem, starttem, endtem;
            NekDouble beta, alpha;
            NekDouble vExchange = 0;
            // temporary Array
            Array<OneD, NekDouble> r0(nGlobal, 0.0);
            Array<OneD, NekDouble> tmp1;
            Array<OneD, NekDouble> tmp2;
            Array<OneD, NekDouble> Vsingle1;
            Array<OneD, NekDouble> Vsingle2;
            Array<OneD, NekDouble> hsingle1;
            Array<OneD, NekDouble> hsingle2;

            ///////////////////////////////////////////////////////////////////////////////
            // // tmp2 for preconditioner multiplication, later consider it
            if(restarted)
            {
                // tmp2=Ax
                v_DoMatrixMultiply(pOutput, r0);

                //The first search direction
                beta = -1.0;
                //PYT: r0=b-AX
                Vmath::Svtvp(nNonDir, beta, &r0[0] + nDir, 1, &pInput[0] + nDir, 1, &r0[0] + nDir, 1);
            }
            else
            {
                // If not restarted, x0 should be zero
                Vmath::Vcopy(nNonDir, &pInput[0] + nDir, 1, &r0[0] + nDir, 1);

            }

            tmp1 = r0 + nDir;
            tmp2 = r0 + nDir;
            m_precon->DoPreconditioner(tmp1, tmp2);

            // norm of (r0)
            // m_map tells how to connect
            vExchange    = Vmath::Dot2(nNonDir,
                                       &r0[0] + nDir,
                                       &r0[0] + nDir,
                                       &m_map[0] + nDir);

            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
            eps = vExchange;

            eta[0] = sqrt(eps);

            if(!restarted)
            {
                if(m_prec_factor == NekConstants::kNekUnsetDouble)
                {
                    // evaluate initial residual error for exit check
                    vExchange    = Vmath::Dot2(nNonDir,
                                               &pInput[0] + nDir,
                                               &pInput[0] + nDir,
                                               &m_map[0] + nDir);
                    vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
                    m_prec_factor = vExchange / eps;
                }
            }

            // If input residual is less than tolerance skip solve.
            if (eps * m_prec_factor < m_tolerance * m_tolerance * m_rhs_magnitude)
            {
                m_converged = true;
                return eps;
            }

            // Give an order for the entries in Hessenburg matrix
            for(int nd = 0; nd < m_maxstorage; ++nd)
            {
                id[nd] = nd;
                id_end[nd] = nd + 1;
                starttem = id_end[nd] - m_maxhesband;
                if(truncted && (starttem) > 0)
                {
                    id_start[nd] = starttem;
                }
                else
                {
                    id_start[nd] = 0;
                }
            }

            //Normlized by r0 norm V(:,1)=r0/norm(r0)
            alpha = 1.0 / eta[0];
            //Scalar multiplication
            Vmath::Smul(nNonDir, alpha, &r0[0] + nDir, 1, &V_total[0][0] + nDir, 1);

            // restarted Gmres(m) process
            int nswp = 0;

            for (int nd = 0; nd < m_maxstorage; ++nd)
            {
                Vsingle1 = V_total[nd];
                Vsingle2 = V_total[nd + 1];
                hsingle1 = hes[nd];

                // w here is no need to add nDir due to temporary Array
                idtem = id[nd];
                starttem = id_start[idtem];
                endtem = id_end[idtem];

                DoArnoldi(starttem, endtem, nGlobal, nDir, V_total, Vsingle1, Vsingle2, hsingle1);

                if(starttem > 0)
                {
                    starttem = starttem - 1;
                }

                hsingle2 = Upper[nd];
                Vmath::Vcopy(m_maxstorage + 1, &hsingle1[0], 1, &hsingle2[0], 1);
                DoGivensRotation(starttem, endtem, nGlobal, nDir, cs, sn, hsingle2, eta);

                // This Gmres merge truncted Gmres to accelerate.
                // If truncted, cannot jump out because the last term of eta is not residual
                if((!truncted) || (nd < m_maxhesband))
                {
                    eps = eta[nd + 1] * eta[nd + 1];

                    if (eps * m_prec_factor < m_tolerance * m_tolerance * m_rhs_magnitude)
                    {
                        m_converged = true;
                    }

                }
                nswp++;
                m_totalIterations++;

                if(m_converged)
                {
                    break;
                }
            }

            DoBackward(nswp, Upper, eta, y_total);

            // calculate output x = x + y_total*V_total
            for(int i = 0; i < nswp; ++i)
            {
                beta = y_total[i];
                // y_total[i] = beta;

                Vmath::Svtvp(nNonDir, beta, &V_total[i][0] + nDir, 1, &pOutput[0] + nDir, 1, &pOutput[0] + nDir, 1);
            }
            return eps;
        }

        void GlobalLinSysIterative::Set_Rhs_Magnitude(
            const NekVector<NekDouble> &pIn)
        {
            Array<OneD, NekDouble> vExchange(1, 0.0);
            if (m_map.num_elements() > 0)
            {
                vExchange[0] = Vmath::Dot2(pIn.GetDimension(),
                                           &pIn[0], &pIn[0], &m_map[0]);
            }

            m_expList.lock()->GetComm()->GetRowComm()->AllReduce(
                vExchange, Nektar::LibUtilities::ReduceSum);

            // To ensure that very different rhs values are not being
            // used in subsequent solvers such as the velocit solve in
            // INC NS. If this works we then need to work out a better
            // way to control this.
            NekDouble new_rhs_mag = (vExchange[0] > 1e-6) ? vExchange[0] : 1.0;

            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                m_rhs_magnitude = new_rhs_mag;
            }
            else
            {
                m_rhs_magnitude = (m_rhs_mag_sm * (m_rhs_magnitude) +
                                   (1.0 - m_rhs_mag_sm) * new_rhs_mag);
            }
        }
        
    }
}
