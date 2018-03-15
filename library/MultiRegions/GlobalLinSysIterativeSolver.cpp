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
// Description: GlobalLinSysIterativeSolver definition
//
///////////////////////////////////////////////////////////////////////////////

#include <map>
#include <MultiRegions/GlobalLinSysIterativeSolver.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

using namespace std;

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
        // string GlobalLinSysIterativeSolver::className
        //         = GetGlobalLinSysFactory().RegisterCreatorFunction(
        //             "IterativeSolver",
        //             GlobalLinSysIterativeSolver::create,
        //             "Iterative solver for full matrix system.");


        GlobalLinSysIterativeSolver::GlobalLinSysIterativeSolver(
                    const GlobalLinSysKey &pKey,
                    const std::weak_ptr<ExpList> &pExp,
                    const std::shared_ptr<AssemblyMap> &pLocToGloMap)
            : GlobalLinSys         (pKey, pExp, pLocToGloMap),
              GlobalLinSysIterative(pKey, pExp, pLocToGloMap)
        {
            
        }


        /**
         * Constructor for full direct matrix solve.
         * @param   pKey        Key specifying matrix to solve.
         * @param   pExp        Shared pointer to expansion list for applying
         *                      matrix evaluations.
         * @param   pLocToGloMap Local to global mapping.
         */
        void GlobalLinSysIterativeSolver::initializeLinSys(const LibUtilities::SessionReaderSharedPtr &session)
        {
            
            m_nlinsys                   = session->GetParameter("LinSysDimens");
            m_maxdirction               = session->GetParameter("MaxDirction" );
            m_maxrestart                = session->GetParameter("MaxRestart" );
            m_tolerance                 = session->GetParameter("Tolerance" );
            
            m_rhs       =  Array<OneD, NekDouble>(m_nlinsys);
            m_mat       =  Array<OneD, NekDouble>(m_nlinsys*m_nlinsys);
            Vmath::Fill(m_nlinsys,1.0,&m_rhs[0],1);

            // initial the matrix(A) in Ax = f
            for(int i=0; i<m_nlinsys*m_nlinsys; ++i)
            {
                m_mat[i]    =   1.0*i;
            }
            for(int i=0; i<m_nlinsys; ++i)
            {
                m_mat[i*m_nlinsys+i]    =   m_mat[i*m_nlinsys+i]*10.0;
            }
            return;
        }


        /**
         *
         */
        GlobalLinSysIterativeSolver::~GlobalLinSysIterativeSolver()
        {

        }




        /**
         *
         */
        void GlobalLinSysIterativeSolver::v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            //ASSERTL1(m_nlinsys > pInput.num_elements()+pInput.GetOffset(),"Array out of bounds");
            //ASSERTL1(m_nlinsys > pInput.num_elements()+pInput.GetOffset(),"Array out of bounds");
            for(int i=0;i<m_nlinsys;++i)
            {
                pOutput[i]  =   Vmath::Dot(m_nlinsys,&m_mat[i*m_nlinsys],&pInput[0]);
                pOutput[i]  =   sqrt(pOutput[i]);
            }
            return;
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
         *                              of freedom with initial guess on other values.
         *                              On output, the solution \f$x\f$.
         * @param           pLocToGloMap    Local to global mapping.
         * @param           pDirForcing Precalculated Dirichlet forcing.
         */
        void GlobalLinSysIterativeSolver::v_Solve(
                    const Array<OneD, const NekDouble>  &pInput,
                          Array<OneD,       NekDouble>  &pOutput,
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &pDirForcing)
        {
            std::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            bool vCG;
            if ((m_locToGloMap = std::dynamic_pointer_cast<AssemblyMapCG>(
                     pLocToGloMap)))
            {
                vCG = true;
            }
            else if ((m_locToGloMap = std::dynamic_pointer_cast<
                          AssemblyMapDG>(pLocToGloMap)))
            {
                vCG = false;
            }
            else
            {
                ASSERTL0(false, "Unknown map type");
            }

            bool dirForcCalculated = (bool) pDirForcing.num_elements();
            int nDirDofs  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobDofs = pLocToGloMap->GetNumGlobalCoeffs();
            int nDirTotal = nDirDofs;
            
            expList->GetComm()->GetRowComm()
                   ->AllReduce(nDirTotal, LibUtilities::ReduceSum);
            
            Array<OneD, NekDouble> tmp(nGlobDofs), tmp2;

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
                    Array<OneD, NekDouble> out(nGlobDofs,0.0);

                    // solve for perturbation from intiial guess in pOutput
                    SolveLinearSystem(
                        nGlobDofs, tmp, out, pLocToGloMap, nDirDofs);
                    Vmath::Vadd(nGlobDofs-nDirDofs,    &out    [nDirDofs], 1,
                                &pOutput[nDirDofs], 1, &pOutput[nDirDofs], 1);
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
        void GlobalLinSysIterativeSolver::v_UniqueMap()
        {
            m_map = m_locToGloMap->GetGlobalToUniversalMapUnique();
        }



    }
}
