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
//#define GMRES_DEBUG
//#define GMRES_PRECDTN
//#define ITERATIVE_DEBUG_OUTPUT        

#include <MultiRegions/GlobalLinSysIterative.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        std::string GlobalLinSysIterative::lookupIds[2] = {
            LibUtilities::SessionReader::RegisterEnumValue(
                "IterativeMethod", "ConjugateGradient", eConjugateGradient),
            LibUtilities::SessionReader::RegisterEnumValue(
                "IterativeMethod", "GMRES", eGMRES),
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

            LibUtilities::CommSharedPtr vComm = m_expList.lock()->GetComm()->GetRowComm();
            m_root    = (vComm->GetRank())? false : true;

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
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &plocToGloMap,
                    const int nDir)
        {
            
            IterativeMethodType pType = plocToGloMap->GetIteraterType();
            switch(pType)
            {
            case eGMRES:
                DoGMRES(nGlobal, pInput, pOutput, plocToGloMap, nDir);
                break;
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
         * @param       pInput      Input residual  of all DOFs.  
         * @param       pOutput     Solution vector of all DOFs.  
         */
        void GlobalLinSysIterative::DoConjugateGradient(
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
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

        void GlobalLinSysIterative::DoGMRES(
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr &plocToGloMap,
            const int                          nDir)
        {
            
            m_prec_factor = NekConstants::kNekUnsetDouble;
            if (!m_precon)
            {
                v_UniqueMap();
                m_precon = CreatePrecon(plocToGloMap);
                m_precon->BuildPreconditioner();
            }

              
            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                NekVector<NekDouble> inGlob (nGlobal, pInput, eWrapper);
                Set_Rhs_Magnitude(inGlob);
            }


            // Get vector sizes
            NekDouble eps = 0.0;
            int nNonDir = nGlobal - nDir;

            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            Vmath::Zero(nNonDir, &pOutput[nDir],1);
            //Vmath::Zero(nGlobal, &qk_a[iqk0],1);


            m_totalIterations = 0;
            m_converged       = false;
            
            bool restarted = false;

            m_maxdirction = 50;
            m_maxrestart = m_maxiter/m_maxdirction+1;


            for(int nrestart=0;nrestart<m_maxrestart;++nrestart)
            {
                eps = DoGmresRestart(restarted, nGlobal,pInput,pOutput,nDir);


                if(m_converged)
                {   
                    if (m_verbose && m_root)
                    {
                        cout << "GMRES iterations made = " << m_totalIterations 
                             << " using tolerance of "  << m_tolerance 
                             << " (error = " << sqrt(eps*m_prec_factor/m_rhs_magnitude) 
                             << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")" 
                             << endl;
                    }

                    return;

                }
                restarted = true;
            }
        
            
            if(m_root)
            {
                cout << "GMRES iterations made = " << m_totalIterations 
                     << " using tolerance of "  << m_tolerance 
                     << " (error = " << sqrt(eps*m_prec_factor/m_rhs_magnitude)
                     << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                     << endl;
            }
            ROOTONLY_NEKERROR(ErrorUtil::efatal,
                              "Exceeded maximum number of iterations");

        }


        /**  
         * Solve a global linear system(Ax=f, r =f-Ax) using the conjugate gradient method.  
         * We solve only for the non-Dirichlet modes. The operator is evaluated  
         * using an auxiliary function v_DoMatrixMultiply defined by the  
         * specific solver. Distributed math routines are used to support  
         * parallel execution of the solver.  
         *  
         * The implemented algorithm uses a reduced-communication reordering of  
         * the standard PCG method (Demmel, Heath and Vorst, 1993)  
         *  
         * @param       pInput      Input residual(f)  of all DOFs.  
         * @param       pOutput     Solution vector(x) of all DOFs.  
         */
        NekDouble GlobalLinSysIterative::DoGmresRestart(
            const bool                         rested,
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const int                          nDir)
        {

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();
            
            int nNonDir = nGlobal - nDir;


            // Allocate array storage of coefficients
            Array<TwoD, NekDouble> hes    (m_maxdirction+1,m_maxdirction, 0.0);
            Array<OneD, NekDouble> eta    (m_maxdirction+1, 0.0);
            Array<OneD, NekDouble> cs     (m_maxdirction, 0.0);
            Array<OneD, NekDouble> sn     (m_maxdirction, 0.0);
            Array<OneD, NekDouble> yk     (m_maxdirction, 0.0);
            Array<OneD, NekDouble> tmp0;
            Array<OneD, NekDouble> tmp1;

            int nqk =   (m_maxdirction+1)*nGlobal;
            int iqk0 =0;
            int iqk1 =0;
            // Allocate array storage of search directions
            Array<OneD, NekDouble> qk_a   (nqk, 0.0);

            // Allocate array storage of temperary search directions
            // theoretically not needed but v_DoMatrixMultiply uses the size of its second para-
            // meter, which leads to error when using tmp1 = qk_a+iqk1
            // todo modify v_DoMatrixMultiply to eliminate the need of these two temperary arrays
            Array<OneD, NekDouble> tm_a   (nGlobal, 0.0);
            Array<OneD, NekDouble> tm_b   (nGlobal, 0.0);


            NekDouble alpha, beta, eps, dd, hh,temp_dbl;
            NekDouble   vExchange=0.0;
            Array<OneD, NekDouble> vExchange_a(m_maxdirction,0.0);
            for(int nd=0;nd<m_maxdirction+1;++nd)
            {
                iqk0 = nd*nGlobal;
                Vmath::Zero(nGlobal,&qk_a[iqk0],1);
            }

            iqk0 = 0;
            if(rested)
            {
                // qk_a[0] = A*x0 if not restarted x0 should be zero!!
                v_DoMatrixMultiply(pOutput, tm_b);
                tmp0 = qk_a+iqk0;
                for(int k =0;k<nGlobal;++k)
                {
                    tmp0[k] = tm_b[k];
                }
                // q_k[0] = f-A*x0
                beta = -1.0;
                Vmath::Svtvp(nNonDir, beta, &qk_a[iqk0]+nDir, 1, &pInput[0]+nDir, 1, &qk_a[iqk0]+nDir, 1);
            }
            else
            {
                Vmath::Vcopy(nNonDir, &pInput[0]+nDir, 1, &qk_a[iqk0]+nDir, 1);
            }

            // if null preconditioner is used, vcopy is wrong without the following copy process 
            iqk0 = 0;
            for(int k =0;k<nGlobal;++k)
            {
                tm_b[k] = qk_a[iqk0+k];
            }
            tmp0 = tm_b + nDir;
            tmp1 = tm_a + nDir;
            m_precon->DoPreconditioner(tmp0, tmp1);
            for(int k =0;k<nGlobal;++k)
            {
                qk_a[iqk0+k] = tm_a[k];
            }

            // evaluate initial residual error for exit check
            vExchange    = Vmath::Dot2(nNonDir,
                                       &qk_a[iqk0]+nDir,
                                       &qk_a[iqk0]+nDir,
                                       &m_map[0] + nDir);
            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

            eps          = vExchange;
            eta[0]       = sqrt(eps);

            if(!rested)
            {
                if(m_prec_factor == NekConstants::kNekUnsetDouble)
                {

                    // evaluate initial residual error for exit check
                    vExchange    = Vmath::Dot2(nNonDir,
                                               &pInput[0]+nDir,
                                               &pInput[0]+nDir,
                                               &m_map[0] + nDir);
                    vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
                    m_prec_factor = vExchange/eps;
                    
                }
            }
            

            // If input residual is less than tolerance skip solve.
            if (eps*m_prec_factor < m_tolerance * m_tolerance * m_rhs_magnitude)
            {
                m_converged = true;
                return eps;
            }

            // normalize the first search vector
            alpha        = 1.0/eta[0] ;
            iqk0 = 0;
            Vmath::Smul(nNonDir,alpha,&qk_a[iqk0]+nDir,1,&qk_a[iqk0]+nDir,1);

            // to store the total search direction in the following for loop
            int nswp = 0;
            // iterative loop 
            for(int nd=0; nd<(m_maxdirction); ++nd)
            {
                iqk0 = nd*nGlobal;
                iqk1 = (nd+1)*nGlobal;

                tmp0 = qk_a+iqk0 ;
                v_DoMatrixMultiply(tmp0, tm_b);
                


                tmp0 = tm_b + nDir;
                tmp1 = tm_a + nDir;
                m_precon->DoPreconditioner(tmp0, tmp1);

                for(int k =0;k<nNonDir;++k)
                {
                    qk_a[iqk1+nDir+k] = tmp1[k];
                }

                // orthogonalize new search vector with previous ones(3 steps)
                // 1st step: calculate coefficient of Hessenberg matrix()
                for(int i=0;i<nd+1;++i)
                {
                    int iqki = i*nGlobal;
                    vExchange    = Vmath::Dot2(nNonDir,
                                               &qk_a[iqki]+nDir,
                                               &qk_a[iqk1]+nDir,
                                               &m_map[0] + nDir);
                    vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
                    hes[i][nd] = vExchange;
                
                    beta = -1.0*vExchange;
                    Vmath::Svtvp(nNonDir, beta, &qk_a[iqki]+nDir, 1, &qk_a[iqk1]+nDir, 1, &qk_a[iqk1]+nDir, 1);

                }

                // 1st step: calculate coefficient of Hessenberg matrix()
                // for(int i=0;i<nd+1;++i)
                // {
                //     int iqki = i*nGlobal;
                //     vExchange    = Vmath::Dot2(nNonDir,
                //                                &qk_a[iqki]+nDir,
                //                                &qk_a[iqk1]+nDir,
                //                                &m_map[0] + nDir);
                //     hes[i][nd] = vExchange;
                // }

                // for(int i=0;i<nd+1;++i)
                // {
                //     vExchange_a[i]    = hes[i][nd];
                // }
                // vComm->AllReduce(vExchange_a, Nektar::LibUtilities::ReduceSum);

                // for(int i=0;i<nd+1;++i)
                // {
                //     hes[i][nd] = vExchange_a[i];
                // }
                
                

                // for(int i=0;i<nd+1;++i)
                // {
                //     int iqki = i*nGlobal;
                //     vExchange   = hes[i][nd];
                //     beta = -1.0*vExchange;
                //     Vmath::Svtvp(nNonDir, beta, &qk_a[iqki]+nDir, 1, &qk_a[iqk1]+nDir, 1, &qk_a[iqk1]+nDir, 1);
                // }
                
                // calculate the L2 norm and normalize
                vExchange    = Vmath::Dot2(nNonDir,
                                           &qk_a[iqk1]+nDir,
                                           &qk_a[iqk1]+nDir,
                                           &m_map[0] + nDir);
                vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
                hes[nd+1][nd] = sqrt(vExchange);
                alpha        = 1.0/hes[nd+1][nd] ;
                Vmath::Smul(nNonDir,alpha,&qk_a[iqk1]+nDir,1,&qk_a[iqk1]+nDir,1);

                // multiply rotate matrix to QR factorization and calculate error
                // R matrix is stored in hes
                for(int i=0;i<nd;++i)
                {
                    temp_dbl        = cs[i]*hes[i][nd] - sn[i]*hes[i+1][nd];
                    hes[i+1][nd]    = sn[i]*hes[i][nd] + cs[i]*hes[i+1][nd];
                    hes[i][nd]      = temp_dbl;
                }
                dd = hes[nd][nd];
                hh = hes[nd+1][nd];
                if(hh==0.0)
                {
                    cs[nd] = 1.0;
                    sn[nd] = 0.0;
                }
                else if (abs(hh) > abs(dd))
                {
                    temp_dbl = -dd/hh;
                    sn[nd] = 1.0 / sqrt(1.0 + temp_dbl*temp_dbl);
                    cs[nd] = temp_dbl * sn[nd];
                }
                else
                {
                    temp_dbl = -hh/dd;
                    cs[nd] = 1.0 / sqrt(1.0 + temp_dbl*temp_dbl);
                    sn[nd] = temp_dbl * cs[nd];
                }
                hes[nd][nd] = cs[nd]*hes[nd][nd]-sn[nd]*hes[nd+1][nd];
                hes[nd+1][nd] = 0.0;

                temp_dbl        = cs[nd]*eta[nd] - sn[nd]*eta[nd+1];
                eta[nd+1]       = sn[nd]*eta[nd] + cs[nd]*eta[nd+1];
                eta[nd]         = temp_dbl;

                eps          = eta[nd+1]*eta[nd+1];

                nswp++;
                m_totalIterations++;

                // If input residual is less than tolerance break the for loop.
                if (eps*m_prec_factor < m_tolerance * m_tolerance * m_rhs_magnitude)
                {
                    m_converged = true;
                    break;
                }
            }

            // calculate yk using R matrix and eta.
            for(int i=0;i<nswp;++i)
            {
                yk[i] = eta[i];
            }

            for(int i=nswp-1;i>-1;--i)
            {
                yk[i] = yk[i]/hes[i][i];
                for(int j=0;j<i;++j)
                {
                    yk[j] = yk[j]-hes[j][i]*yk[i];
                }
            }

            // calculate output x = x + yk*qk_a
            for(int i=0;i<nswp;++i)
            {
                beta = yk[i];
                int iqki = i*nGlobal;
                Vmath::Svtvp(nNonDir, beta, &qk_a[iqki]+nDir, 1, &pOutput[0]+nDir, 1, &pOutput[0]+nDir, 1);
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
                                        &pIn[0],&pIn[0],&m_map[0]);
            }

            m_expList.lock()->GetComm()->GetRowComm()->AllReduce(
                vExchange, Nektar::LibUtilities::ReduceSum);

            // To ensure that very different rhs values are not being
            // used in subsequent solvers such as the velocit solve in
            // INC NS. If this works we then need to work out a better
            // way to control this.
            NekDouble new_rhs_mag = (vExchange[0] > 1e-6)? vExchange[0] : 1.0;

            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                m_rhs_magnitude = new_rhs_mag;
            }
            else
            {
                m_rhs_magnitude = (m_rhs_mag_sm*(m_rhs_magnitude) + 
                                   (1.0-m_rhs_mag_sm)*new_rhs_mag); 
            }
        }

    }
}
