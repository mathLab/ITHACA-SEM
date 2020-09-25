///////////////////////////////////////////////////////////////////////////////
//
// File:  NekLinSysIteratCG.cpp
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
// Description:  NekLinSysIteratCG definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekLinSysIteratCG.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{      
    namespace LibUtilities
    {  
        /**
         * @class  NekLinSysIteratCG
         *
         * Solves a linear system using iterative methods.
         */
        string NekLinSysIteratCG::className =
        LibUtilities::GetNekLinSysIteratFactory().RegisterCreatorFunction(
            "ConjugateGradient", NekLinSysIteratCG::create,
            "NekLinSysIteratCG solver.");

        /// Constructor for full direct matrix solve.
        NekLinSysIteratCG::NekLinSysIteratCG(
            const LibUtilities::SessionReaderSharedPtr  &pSession,
            const LibUtilities::CommSharedPtr           &vComm,
            const int                                   nDimen)
            : NekLinSysIterat(pSession, vComm, nDimen)
        {
            std::vector<std::string>  variables(1);
            variables[0] =  pSession->GetVariable(0);
            string variable = variables[0];

            if(pSession->DefinesGlobalSysSolnInfo(variable, "SuccessiveRHS"))
            {
                m_successiveRHS = boost::lexical_cast<int>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "SuccessiveRHS").c_str());
            }
            else
            {
                pSession->LoadParameter("SuccessiveRHS",
                                        m_successiveRHS, 0);
            }

            int successiveRHS;
            if((successiveRHS = m_successiveRHS))
            {
                m_prevLinSol.set_capacity(successiveRHS);
                m_useProjection = true;
            }
            else
            {
                m_useProjection = false;
            }

        }

        void NekLinSysIteratCG::v_InitObject()
        {
            NekLinSysIterat::v_InitObject();
        }


        NekLinSysIteratCG::~NekLinSysIteratCG()
        {
        }

        /**
         *
         */
        int NekLinSysIteratCG::v_SolveSystem(
            const int                           nGlobal,
            const Array<OneD, const NekDouble>  &pInput,
            Array<OneD,      NekDouble>         &pOutput,
            const int                           nDir,
            const NekDouble                     tol,
            const NekDouble                     factor)
        {
            boost::ignore_unused(tol);

            m_tolerance = max(tol, 1.0E-16);
            m_prec_factor = factor;
            if (m_useProjection)
            {
                DoAconjugateProjection(nGlobal, pInput, pOutput, nDir);
            }
            else
            {
                // applying plain Conjugate Gradient
                DoConjugateGradient(nGlobal, pInput, pOutput, nDir);
            }

            return 0;
        }

        /**
         * This method implements A-conjugate projection technique
         * in order to speed up successive linear solves with
         * right-hand sides arising from time-dependent discretisations.
         * (P.F.Fischer, Comput. Methods Appl. Mech. Engrg. 163, 1998)
         */
        void NekLinSysIteratCG::DoAconjugateProjection(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const int nDir)
        {

            // Get vector sizes
            int nNonDir = nGlobal - nDir;
            Array<OneD, NekDouble> tmp;

            if (0 == m_numPrevSols)
            {
                // no previous solutions found, call CG

                DoConjugateGradient(nGlobal, pInput, pOutput, nDir);

                UpdateKnownSolutions(nGlobal, pOutput, nDir);
            }
            else
            {
                // Create NekVector wrappers for linear algebra operations
                NekVector<NekDouble> b     (nNonDir, pInput  + nDir, eWrapper);
                NekVector<NekDouble> x     (nNonDir, tmp = pOutput + nDir,
                                            eWrapper);

                // check the input vector (rhs) is not zero

                NekDouble rhsNorm = Vmath::Dot2(nNonDir,
                                                pInput + nDir,
                                                pInput + nDir,
                                                m_map + nDir);

                m_Comm->AllReduce(rhsNorm, Nektar::LibUtilities::ReduceSum);

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

                NekVector<NekDouble> pb(nNonDir, tmp = pb_s + nDir, eWrapper);
                NekVector<NekDouble> px(nNonDir, tmp = px_s + nDir, eWrapper);
                NekVector<NekDouble> tmpAx(nNonDir, tmp = tmpAx_s + nDir, 
                                           eWrapper);
                NekVector<NekDouble> tmpx(nNonDir, tmp = tmpx_s  + nDir,
                                     eWrapper);


                // notation follows the paper cited:
                // \alpha_i = \tilda{x_i}^T b^n
                // projected x, px = \sum \alpha_i \tilda{x_i}

                Array<OneD, NekDouble> alpha(m_prevLinSol.size(), 0.0);
                for (int i = 0; i < m_prevLinSol.size(); ++i)
                {
                    alpha[i] = Vmath::Dot2(nNonDir,
                                           m_prevLinSol[i],
                                           pInput + nDir,
                                           m_map + nDir);
                }
                m_Comm->AllReduce(alpha, Nektar::LibUtilities::ReduceSum);

                for (int i = 0; i < m_prevLinSol.size(); ++i)
                {
                    if (alpha[i] < NekConstants::kNekZeroTol)
                    {
                        continue;
                    }

                    NekVector<NekDouble> xi(nNonDir, m_prevLinSol[i], eWrapper);
                    px += alpha[i] * xi;
                }

                // pb = b^n - A px
                Vmath::Vcopy(nNonDir,
                             pInput.get() + nDir, 1,
                             pb_s.get()   + nDir, 1);

                m_operator.DoNonlinLinSysLhsEval(px_s, tmpAx_s);

                pb -= tmpAx;


                // solve the system with projected rhs
                DoConjugateGradient(nGlobal, pb_s, tmpx_s, nDir);


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
        NekDouble NekLinSysIteratCG::CalculateAnorm(
                                        const int nGlobal,
                                        const Array<OneD,const NekDouble> &in,
                                        const int nDir)
        {
            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Allocate array storage
            Array<OneD, NekDouble> tmpAx_s(nGlobal, 0.0);

            m_operator.DoNonlinLinSysLhsEval(in, tmpAx_s);

            NekDouble anorm_sq = Vmath::Dot2(nNonDir,
                                             in      + nDir,
                                             tmpAx_s + nDir,
                                             m_map   + nDir);
            m_Comm->AllReduce(anorm_sq, Nektar::LibUtilities::ReduceSum);
            return std::sqrt(anorm_sq);
        }

        /**
         * Updates the storage of previously known solutions.
         * Performs normalisation of input vector wrt A-norm.
         */
        void NekLinSysIteratCG::UpdateKnownSolutions(
                                        const int nGlobal,
                                        const Array<OneD,const NekDouble> &newX,
                                        const int nDir)
        {
            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Check the solution is non-zero
            NekDouble solNorm = Vmath::Dot2(nNonDir,
                                            newX + nDir,
                                            newX + nDir,
                                            m_map + nDir);
            m_Comm->AllReduce(solNorm, Nektar::LibUtilities::ReduceSum);

            if (solNorm < NekConstants::kNekZeroTol)
            {
                return;
            }

            // Allocate array storage
            Array<OneD, NekDouble> tmpAx_s    (nGlobal, 0.0);
            Array<OneD, NekDouble> px_s       (nGlobal, 0.0);
            Array<OneD, NekDouble> tmp1, tmp2;

            // Create NekVector wrappers for linear algebra operations
            NekVector<NekDouble>px   (nNonDir, tmp1 = px_s    + nDir, eWrapper);
            NekVector<NekDouble>tmpAx(nNonDir, tmp2 = tmpAx_s + nDir, eWrapper);


            // calculating \tilda{x} - sum \alpha_i\tilda{x}_i

            Vmath::Vcopy(nNonDir,
                         tmp1 = newX + nDir, 1,
                         tmp2 = px_s + nDir, 1);

            if (m_prevLinSol.size() > 0)
            {
                m_operator.DoNonlinLinSysLhsEval(newX, tmpAx_s);
            }

            Array<OneD, NekDouble> alpha (m_prevLinSol.size(), 0.0);
            for (int i = 0; i < m_prevLinSol.size(); i++)
            {
                alpha[i] = Vmath::Dot2(nNonDir,
                                       m_prevLinSol[i],
                                       tmpAx_s + nDir,
                                       m_map + nDir);
            }
            m_Comm->AllReduce(alpha, Nektar::LibUtilities::ReduceSum);

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
            // solutions generated for 
            // Test_Tet_equitri.xml (IncNavierStokesSolver).
            // Not going to store identically zero solutions.

            NekDouble anorm = CalculateAnorm(nGlobal, px_s, nDir);
            if (anorm < NekConstants::kNekZeroTol)
            {
                return;
            }

            // normalisation of new solution
            Vmath::Smul(nNonDir, 1.0/anorm, 
                        px_s.get() + nDir, 1, px_s.get() + nDir, 1);

            // updating storage with non-Dirichlet-dof part of new solution vector
            m_prevLinSol.push_back(px_s + nDir);
            m_numPrevSols++;
        }



        /**  
         * Solve a global linear system using the conjugate gradient method.  
         * We solve only for the non-Dirichlet modes. The operator is evaluated  
         * using an auxiliary function m_operator.DoNonlinLinSysLhsEval defined by the  
         * specific solver. Distributed math routines are used to support  
         * parallel execution of the solver.  
         *  
         * The implemented algorithm uses a reduced-communication reordering of  
         * the standard PCG method (Demmel, Heath and Vorst, 1993)  
         *  
         * @param       pInput      Input residual  of all DOFs.  
         * @param       pOutput     Solution vector of all DOFs.  
         */
        void NekLinSysIteratCG::DoConjugateGradient(
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const int                          nDir)
        {
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
            NekDouble alpha, beta, rho, rho_new, mu, eps,  min_resid;
            Array<OneD, NekDouble> vExchange(3,0.0);

            // Copy initial residual from input
            r = in;
            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            Vmath::Zero(nNonDir,tmp = pOutput + nDir,1);


            // evaluate initial residual error for exit check
            vExchange[2] = Vmath::Dot2(nNonDir,
                                       r_A,
                                       r_A,
                                       m_map + nDir);

            m_Comm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
            
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

            m_operator.DoNonlinLinPrecond(r_A, tmp = w_A + nDir);

            m_operator.DoNonlinLinSysLhsEval(w_A, s_A);

            k = 0;

            vExchange[0] = Vmath::Dot2(nNonDir,
                                       r_A,
                                       w_A + nDir,
                                       m_map + nDir);

            vExchange[1] = Vmath::Dot2(nNonDir,
                                       s_A + nDir,
                                       w_A + nDir,
                                       m_map + nDir);

            m_Comm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

            rho               = vExchange[0];
            mu                = vExchange[1];
            min_resid         = m_rhs_magnitude;
            beta              = 0.0;
            alpha             = rho/mu;
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
                             << " (error = " << sqrt(eps/m_rhs_magnitude)
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
                m_operator.DoNonlinLinPrecond(r_A, tmp = w_A + nDir);

                // Perform the method-specific matrix-vector multiply operation.
                m_operator.DoNonlinLinSysLhsEval(w_A, s_A);

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
                m_Comm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

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
                             << " (error = " << sqrt(eps/m_rhs_magnitude)
                             << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                             << endl;
                    }
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
    }
}