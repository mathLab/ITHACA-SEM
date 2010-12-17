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
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysIterativeCG::className
                = GlobalLinSysFactory::RegisterCreatorFunction(
                    "IterativeCG",
                    GlobalLinSysIterativeCG::create,
                    "Iterative conjugate gradient solver.");


        /// Constructor for full direct matrix solve.
        GlobalLinSysIterativeCG::GlobalLinSysIterativeCG(const GlobalLinSysKey &pLinSysKey,
                     const boost::shared_ptr<LocalMatrixSystem> &pLocMatSys,
                     const boost::shared_ptr<LocalToGlobalBaseMap>
                                                            &pLocToGloMap)
                : GlobalLinSysIterative(pLinSysKey, pLocMatSys, pLocToGloMap)
        {
            ASSERTL1(m_linSysKey.GetGlobalSysSolnType()==eIterativeCG,
                     "This routine should only be used when using an Iterative CG"
                     " matrix solve");
            m_locMatSys = pLocMatSys;

            // Initialise diagonal preconditioner
            LocalToGlobalC0ContMapSharedPtr vLocToGloMap
                = boost::dynamic_pointer_cast<LocalToGlobalC0ContMap>(pLocToGloMap);

            ComputeDiagonalPreconditioner(vLocToGloMap);
        }

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
        void GlobalLinSysIterativeCG::Solve( const Array<OneD,const NekDouble> &pInput,
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
            NekVector<NekDouble> p_g(nGlobal,p_global,eWrapper);
            NekVector<NekDouble> p(nNonDir,p_global + nDir, eWrapper);
            NekVector<NekDouble> tmp_g(nGlobal,tmp_global,eWrapper);
            NekVector<NekDouble> tmp(nNonDir,tmp_global + nDir,eWrapper);
            NekVector<NekDouble> local_tmp(nLocal);
            int k;
            NekDouble alpha;
            NekDouble beta;
            DNekScalBlkMat &A = *m_locMatSys->GetLocalSystem()[0];
            DNekMat &M = m_preconditioner;

            // Initialise with zero as the initial guess.
            r = in;
            z = M * r;
            p = r;
            k = 0;

            // Continue until convergence
            while (true)
            {
                // Perform matrix-vector operation in local space
                m_locToGloMap->GlobalToLocal(p_g,local_tmp);
                local_tmp = A*local_tmp;
                m_locToGloMap->Assemble(local_tmp,tmp_g);

                alpha = p.Dot(tmp);
                alpha = r.Dot(r)/alpha;
                out   = out + alpha*p;
                r_new = r   - alpha*tmp;

                // Test if residual is small enough
                if (r_new.L2Norm() < 1e-8)
                {
                    break;
                }

                // Update
                z_new = M * r_new;
                beta = r_new.Dot(z_new) / r.Dot(z);
                p = z_new + beta*p;
                r = r_new;
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
        void GlobalLinSysIterativeCG::Solve( const Array<OneD, const NekDouble> &pInput,
                          Array<OneD,       NekDouble> &pOutput,
                    const LocalToGlobalBaseMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble> &pDirForcing)
        {
            LocalToGlobalC0ContMapSharedPtr vLocToGloMap
                = boost::dynamic_pointer_cast<LocalToGlobalC0ContMap>(pLocToGloMap);
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
                    DNekScalBlkMat &Mat = *m_locMatSys->GetLocalSystem()[0];
                    Array<OneD, NekDouble> local_tmp(nLocDofs);
                    NekVector<NekDouble> V_glob(nGlobDofs,global_tmp,eWrapper);
                    NekVector<NekDouble> V_loc(nLocDofs,local_tmp,eWrapper);
                    NekVector<NekDouble> V_out(nGlobDofs,pOutput,eWrapper);

                    vLocToGloMap->GlobalToLocal(V_out,V_loc);
                    V_loc = Mat*V_loc;
                    vLocToGloMap->Assemble(V_loc,V_glob);

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

        void GlobalLinSysIterativeCG::ComputeDiagonalPreconditioner(
                const boost::shared_ptr<LocalToGlobalC0ContMap> &pLocToGloMap)
        {
            DNekScalBlkMat &Mat = *m_locMatSys->GetLocalSystem()[0];

            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nLocal  = pLocToGloMap->GetNumLocalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt = nGlobal - nDir;
            m_preconditioner = DNekMat(nInt, nInt, eDIAGONAL);

            Array<OneD, int> vMap = pLocToGloMap->GetLocalToGlobalMap();

            NekVector<NekDouble> test_local(nLocal, 0.0);
            for (unsigned int i = nDir; i < nGlobal; ++i)
            {
                NekVector<NekDouble> test(nGlobal, 0.0);
                test[i] = 1.0;
                pLocToGloMap->GlobalToLocal(test, test_local);
                test_local = Mat * test_local;
                pLocToGloMap->Assemble(test_local, test);
                m_preconditioner.SetValue(i-nDir,i-nDir,1.0/test[i]);
            }
        }

        void GlobalLinSysIterativeCG::ComputeDiagonalPreconditionerSum(
                const boost::shared_ptr<LocalToGlobalC0ContMap> &pLocToGloMap)
        {
            // TODO: Some entries are being computed incorrectly.
            // Need to fix this.
            DNekScalBlkMat &Mat = *m_locMatSys->GetLocalSystem()[0];

            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nLocal  = pLocToGloMap->GetNumLocalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt = nGlobal - nDir;
            m_preconditioner = DNekMat(nInt, nInt, eDIAGONAL);

            Array<OneD, int> vMap = pLocToGloMap->GetLocalToGlobalMap();

            // Alternative implementation through summing matrix entries.
            Array<OneD, NekDouble> vOutput(nInt,0.0);
            // Scan through the local to global map
            for (unsigned int i = 0; i < nLocal; ++i)
            {
                // Retrieve the corresponding global ID and check if it's been
                // evaluated already.
                int vGlobalIdx = vMap[i];
                if (vGlobalIdx < nDir)
                {
                    continue;
                }

                vOutput[vGlobalIdx-nDir] = 0.0;
                // Scan across the corresponding row of the matrix summing
                // those entries referenced by the local to global map.
                for (unsigned int p = 0; p < nLocal; p++)
                {
                    if (vMap[p] == vGlobalIdx)
                    {
                        vOutput[vGlobalIdx-nDir] += Mat(i,p);
                    }
                }
            }

            // Finally form the preconditioner
            for (unsigned int i = 0; i < nInt; ++i)
            {
                m_preconditioner.SetValue(i,i,1.0/vOutput[i]);
            }
        }
    }
}
