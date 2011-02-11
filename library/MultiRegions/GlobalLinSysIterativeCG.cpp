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

#include <MultiRegions/GlobalLinSysIterativeCG.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>

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
        string GlobalLinSysIterativeCG::className
                = GlobalLinSysFactory::RegisterCreatorFunction(
                    "IterativeCG",
                    GlobalLinSysIterativeCG::create,
                    "Iterative conjugate gradient solver.");


        /**
         * Constructor for full direct matrix solve.
         * @param   pKey        Key specifying matrix to solve.
         * @param   pExp        Shared pointer to expansion list for applying
         *                      matrix evaluations.
         * @param   pLocToGloMap Local to global mapping.
         */
        GlobalLinSysIterativeCG::GlobalLinSysIterativeCG(
                    const GlobalLinSysKey &pKey,
                    const boost::shared_ptr<ExpList> &pExp,
                    const boost::shared_ptr<LocalToGlobalBaseMap> &pLocToGloMap)
                : GlobalLinSysIterative(pKey, pExp, pLocToGloMap)
        {
            ASSERTL1(m_linSysKey.GetGlobalSysSolnType()==eIterativeCG,
                     "This routine should only be used when using an Iterative "
                     "conjugate gradient matrix solve.");
            m_expList = pExp;

            // Initialise diagonal preconditioner
            LocalToGlobalC0ContMapSharedPtr vLocToGloMap
                = boost::dynamic_pointer_cast<LocalToGlobalC0ContMap>(
                                                                pLocToGloMap);

            //ComputeDiagonalPreconditionerSum(vLocToGloMap);
            ComputeNullPreconditioner(vLocToGloMap);
        }


        /**
         *
         */
        GlobalLinSysIterativeCG::~GlobalLinSysIterativeCG()
        {

        }


        /**
         * Solve a global linear system using the conjugate gradient method.
         * We solve only for the non-Dirichlet modes. The operator is evaluated
         * using the local-matrix representation.
         * @param       pInput      Input vector of non-Dirichlet DOFs.
         * @param       pOutput     Solution vector of non-Dirichlet DOFs.
         */
        void GlobalLinSysIterativeCG::Solve(
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput)
        {
            int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
            int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nLocal  = m_locToGloMap->GetNumLocalCoeffs();
            int nNonDir = nGlobal - nDir;
            Array<OneD, NekDouble> p_global(nGlobal, 0.0);
            Array<OneD, NekDouble> tmp_global(nGlobal, 0.0);
            NekVector<NekDouble> in(nNonDir,pInput,eWrapper);
            NekVector<NekDouble> out(nNonDir,pOutput,eWrapper);
            NekVector<NekDouble> r(nNonDir);
            NekVector<NekDouble> r_new(nNonDir);
            NekVector<NekDouble> z(nNonDir);
            NekVector<NekDouble> z_new(nNonDir);
            NekVector<NekDouble> d_g(nGlobal,p_global,eWrapper);
            NekVector<NekDouble> d(nNonDir,p_global + nDir, eWrapper);
            NekVector<NekDouble> tmp_g(nGlobal,tmp_global,eWrapper);
            NekVector<NekDouble> tmp(nNonDir,tmp_global + nDir,eWrapper);
            NekVector<NekDouble> local_tmp(nLocal);
            int k;
            NekDouble alpha;
            NekDouble beta;
            //DNekScalBlkMat &A = *m_locMatSys->GetLocalSystem()[0];
            DNekMat &M = m_preconditioner; // M inverse in algorithm

            // Initialise with zero as the initial guess.
            r = in;
            z = M * r;
            d = z;
            k = 0;

            // Continue until convergence
            while (true)
            {
                // Perform matrix-vector operation A*d_i
                m_expList->GeneralMatrixOp(*m_linSysKey.GetGlobalMatrixKey(),
                                            p_global, tmp_global, true);

                // compute step length
                alpha = d.Dot(tmp);
                alpha = z.Dot(r)/alpha;

                // approximate solution
                out   = out + alpha*d;

                // compute residual
                r_new = r   - alpha*tmp;

                // Test if residual is small enough
                if (r_new.L2Norm() < NekConstants::kNekIterativeTol)
                {
                    break;
                }

                // Apply preconditioner to new residual
                z_new = M * r_new;

                // Improvement achieved
                beta = r_new.Dot(z_new) / r.Dot(z);

                // Compute new search direction
                d = z_new + beta*d;

                // Next step
                r = r_new;
                z = z_new;
                k++;
            }
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
        void GlobalLinSysIterativeCG::Solve(
                    const Array<OneD, const NekDouble>  &pInput,
                          Array<OneD,       NekDouble>  &pOutput,
                    const LocalToGlobalBaseMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &pDirForcing)
        {
            LocalToGlobalC0ContMapSharedPtr vLocToGloMap
                = boost::dynamic_pointer_cast<LocalToGlobalC0ContMap>(
                                                                pLocToGloMap);
            m_locToGloMap = vLocToGloMap;

            bool dirForcCalculated = (bool) pDirForcing.num_elements();
            int nDirDofs  = vLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobDofs = vLocToGloMap->GetNumGlobalCoeffs();
            int nLocDofs  = vLocToGloMap->GetNumLocalCoeffs();
            if(nDirDofs)
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
                Solve(global_tmp+nDirDofs,offsetarray = pOutput+nDirDofs);
            }
            else
            {
                Solve(pInput, pOutput);
            }
        }


        /**
         * Populates preconditioner with the identity to apply no
         * preconditioning.
         * @param   pLocToGloMap    Local to Global mapping.
         */
        void GlobalLinSysIterativeCG::ComputeNullPreconditioner(
                const boost::shared_ptr<LocalToGlobalC0ContMap> &pLocToGloMap)
        {
            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt = nGlobal - nDir;
            m_preconditioner = DNekMat(nInt, nInt, eDIAGONAL);

            for (unsigned int i = 0; i < nInt; ++i)
            {
                m_preconditioner.SetValue(i,i,1.0);
            }
        }


        /**
         * Diagonal preconditioner computed by evaluating the local matrix
         * acting on each basis vector (0,...,0,1,0,...,0). (deprecated)
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysIterativeCG::ComputeDiagonalPreconditioner(
                const boost::shared_ptr<LocalToGlobalC0ContMap> &pLocToGloMap)
        {
            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nLocal  = pLocToGloMap->GetNumLocalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt = nGlobal - nDir;
            m_preconditioner = DNekMat(nInt, nInt, eDIAGONAL);

            Array<OneD, int> vMap = pLocToGloMap->GetLocalToGlobalMap();

            for (unsigned int i = 0; i < nInt; ++i)
            {
                Array<OneD, NekDouble> test(nGlobal, 0.0);
                Array<OneD, NekDouble> test_local(nLocal, 0.0);
                test[i+nDir] = 1.0;
                m_expList->GeneralMatrixOp(*m_linSysKey.GetGlobalMatrixKey(),
                                test, test, true);

                m_preconditioner.SetValue(i,i,1.0/test[i+nDir]);
            }
        }


        /**
         * Diagonal preconditioner computed by summing the relevant elements of
         * the local matrix system.
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysIterativeCG::ComputeDiagonalPreconditionerSum(
                const boost::shared_ptr<LocalToGlobalC0ContMap> &pLocToGloMap)
        {
            int i,j,n,cnt,gid1,gid2;
            NekDouble sign1,sign2,value;
            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt    = nGlobal - nDir;

            NekDouble zero = 0.0;

            // fill global matrix
            DNekScalMatSharedPtr loc_mat;
            Array<OneD, NekDouble> vOutput(nInt,0.0);
            m_preconditioner = DNekMat(nInt, nInt, eDIAGONAL);

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
                                value = vOutput[gid1]
                                            + sign1*sign2*(*loc_mat)(i,j);
                                vOutput[gid1] = value;
                            }
                        }
                    }
                }
                cnt   += loc_lda;
            }

            // Populate preconditioner
            for (unsigned int i = 0; i < nInt; ++i)
            {
                m_preconditioner.SetValue(i,i,1.0/vOutput[i]);
            }
        }
    }
}
