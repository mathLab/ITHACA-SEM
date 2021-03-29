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

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        std::string GlobalLinSysIterative::IteratSolverlookupIds[2] =
        {
            LibUtilities::SessionReader::RegisterEnumValue(
                "LinSysIterSolver", "ConjugateGradient", MultiRegions::
                eConjugateGradient),
            LibUtilities::SessionReader::RegisterEnumValue(
                "LinSysIterSolver", "GMRES", MultiRegions::eGMRES),
        };

        std::string GlobalLinSysIterative::IteratSolverdef =
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                "LinSysIterSolver", "ConjugateGradient");

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
                  m_rhs_mag_sm(0.9),
                  m_precon(NullPreconditionerSharedPtr),
                  m_totalIterations(0),
                  m_useProjection(false),
                  m_numPrevSols(0)
        {
            m_tolerance = pLocToGloMap->GetIterativeTolerance();
            m_maxiter   = pLocToGloMap->GetMaxIterations();
            m_linSysIterSolver = pLocToGloMap->GetLinSysIterSolver();

            LibUtilities::CommSharedPtr vComm = m_expList.lock()->GetComm()->GetRowComm();
            m_root    = (vComm->GetRank())? false : true;

            m_numSuccessiveRHS = pLocToGloMap->GetSuccessiveRHS();
            m_isAconjugate = m_numSuccessiveRHS > 0;
            m_numSuccessiveRHS = std::abs(m_numSuccessiveRHS);
            m_useProjection = m_numSuccessiveRHS > 0;

            if(m_isAconjugate && 0 == m_linSysIterSolver.compare("GMRES"))
            {
                WARNINGL0(false, "To use A-conjugate projection, the matrix should be symmetric positive definite.");
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
            if (!m_linsol)
            {
                LibUtilities::CommSharedPtr v_Comm = 
                m_expList.lock()->GetComm()->GetRowComm();
                LibUtilities::SessionReaderSharedPtr pSession =
                m_expList.lock()->GetSession();
                
                // Check such a module exists for this equation.
                ASSERTL0(LibUtilities::GetNekLinSysIterFactory().
                ModuleExists(m_linSysIterSolver),
                    "NekLinSysIter '" + m_linSysIterSolver +
                    "' is not defined.\n");
                m_linsol = LibUtilities::GetNekLinSysIterFactory().
                           CreateInstance(m_linSysIterSolver, pSession,
                            v_Comm, nGlobal - nDir, LibUtilities::NekSysKey());

                m_NekSysOp.DefineNekSysLhsEval(
                    &GlobalLinSysIterative::DoMatrixMultiplyFlag, this);
                m_NekSysOp.DefineNekSysPrecon(
                    &GlobalLinSysIterative::DoPreconditionerFlag, this);
                m_linsol->SetSysOperators(m_NekSysOp);
                    v_UniqueMap();
                m_linsol->setUniversalUniqueMap(m_map);
            }
            
            if (!m_precon)
            {
                m_precon = CreatePrecon(plocToGloMap);
                m_precon->BuildPreconditioner();
            }

            m_linsol->setRhsMagnitude(m_rhs_magnitude);
            if (m_useProjection)
            {
                DoProjection(nGlobal, pInput, pOutput, nDir, m_tolerance, m_isAconjugate);
            }
            else
            {
                m_linsol->SolveSystem(nGlobal, pInput, pOutput, nDir, m_tolerance);
            }
        }

        /**
         * This method implements projection techniques
         * in order to speed up successive linear solves with
         * right-hand sides arising from time-dependent discretisations.
         * (P.F.Fischer, Comput. Methods Appl. Mech. Engrg. 163, 1998)
         */
        void GlobalLinSysIterative::DoProjection(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const int nDir,
                    const NekDouble tol,
                    const bool isAconjugate)
        {
            int numIterations = 0;
            if (0 == m_numPrevSols)
            {
                // no previous solutions found
                numIterations = m_linsol->SolveSystem(nGlobal, pInput, pOutput, nDir, tol);
            }
            else
            {
                // Get the communicator for performing data exchanges
                LibUtilities::CommSharedPtr vComm
                    = m_expList.lock()->GetComm()->GetRowComm();

                // Get vector sizes
                int nNonDir = nGlobal - nDir;

                // check the input vector (rhs) is not zero
                Array<OneD, NekDouble> tmp;

                NekDouble rhsNorm = Vmath::Dot2(nNonDir,
                                                pInput + nDir,
                                                pInput + nDir,
                                                m_map + nDir);

                vComm->AllReduce(rhsNorm, Nektar::LibUtilities::ReduceSum);

                if (rhsNorm < tol * tol * m_rhs_magnitude)
                {
                    Vmath::Zero(nNonDir, tmp = pOutput+nDir, 1);
                    if(m_verbose && m_root)
                    {
                        cout << "No iterations made"
                        << " using tolerance of " << tol
                        << " (error = " << sqrt(rhsNorm / m_rhs_magnitude)
                        << ", rhs_mag = " << sqrt(m_rhs_magnitude) << ")" << endl;
                    }
                    return;
                }

                // Create NekVector wrappers for linear algebra operations
                NekVector<NekDouble> b     (nNonDir, pInput  + nDir, eWrapper);
                NekVector<NekDouble> x     (nNonDir, tmp = pOutput + nDir, eWrapper);
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

                Array<OneD, NekDouble> alpha(m_prevBasis.size(), 0.0);
                Array<OneD, NekDouble> alphaback(m_prevBasis.size(), 0.0);
                for (int i = 0; i < m_prevBasis.size(); i++)
                {
                    alpha[i] = Vmath::Dot2(nNonDir,
                                           m_prevBasis[i],
                                           pInput + nDir,
                                           m_map + nDir);
                }
                vComm->AllReduce(alpha, Nektar::LibUtilities::ReduceSum);
                int n = m_prevBasis.size(), info = -1;
                Vmath::Vcopy(m_prevBasis.size(), alpha, 1, alphaback, 1);
                Lapack::Dsptrs('U', n, 1, m_coeffMatrixFactor.get(), m_ipivot.get(), alpha.get(), n, info);
                if(info!=0)
                {
                    // Dsptrs fails, only keep the latest solution
                    int latest = ResetKnownSolutionsToLatestOne();
                    alpha[0] = alphaback[latest];
                }
                for (int i = 0; i < m_prevBasis.size(); ++i)
                {
                    NekVector<NekDouble> xi (nNonDir, m_prevLinSol[i], eWrapper);
                    px += alpha[i] * xi;
                }

                // pb = b^n - A px
                Vmath::Vcopy(nNonDir,
                             pInput.get() + nDir, 1,
                             pb_s.get()   + nDir, 1);

                DoMatrixMultiplyFlag(px_s, tmpAx_s, false);

                pb -= tmpAx;

                if (m_verbose)
                {
                    if(m_root)
                        cout << "SuccessiveRHS: " << m_prevBasis.size() << "-bases projection reduces L2-norm of RHS from " << std::sqrt(rhsNorm) << " to ";
                    NekDouble tmprhsNorm = Vmath::Dot2(nNonDir,
                                                       pb_s + nDir,
                                                       pb_s + nDir,
                                                       m_map + nDir);
                    vComm->AllReduce(tmprhsNorm, Nektar::LibUtilities::ReduceSum);
                    if(m_root)
                        cout << std::sqrt(tmprhsNorm) << endl;
                }

                // solve the system with projected rhs
                numIterations = m_linsol->SolveSystem(nGlobal, pb_s, tmpx_s, nDir, tol);

                // remainder solution + projection of previous solutions
                x = tmpx + px;
            }
            // save the auxiliary solution to prev. known solutions
            if(numIterations)
            {
                UpdateKnownSolutions(nGlobal, pOutput, nDir, isAconjugate);
            }
        }

        int GlobalLinSysIterative::ResetKnownSolutionsToLatestOne()
        {
            if(m_numPrevSols==0)
            {
                return -1;
            }
            int latest = (m_numPrevSols -1 + m_numSuccessiveRHS) % m_numSuccessiveRHS;
            Array<OneD, NekDouble> b = m_prevBasis[latest];
            Array<OneD, NekDouble> x = m_prevLinSol[latest];
            m_prevBasis.clear();
            m_prevLinSol.clear();
            m_prevBasis.push_back(b);
            m_prevLinSol.push_back(x);
            m_numPrevSols = 1;
            return latest;
        }

        /**
         * Updates the storage of previously known solutions.
         * Performs normalisation of input vector wrt A-norm.
         */
        void GlobalLinSysIterative::UpdateKnownSolutions(
                                                         const int nGlobal,
                                                         const Array<OneD,const NekDouble> &newX,
                                                         const int nDir,
                                                         const bool isAconjugate)
        {
            // Get vector sizes
            int nNonDir = nGlobal - nDir;
            int insertLocation = m_numPrevSols % m_numSuccessiveRHS;
            int fullbuffer = (m_prevBasis.size() == m_numSuccessiveRHS);

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();

            Array<OneD, NekDouble> tmpAx_s (nGlobal, 0.0);
            Array<OneD, NekDouble> y_s     (m_prevBasis.size() - fullbuffer, 0.0);
            Array<OneD, NekDouble> invMy_s (y_s.size(), 0.0);
            Array<OneD, int>       ipivot  (m_numSuccessiveRHS);
            Array<OneD, NekDouble> tmp, newBasis;

            DoMatrixMultiplyFlag(newX, tmpAx_s, false);

            if(isAconjugate)
            {
                newBasis = newX + nDir;
            }
            else
            {
                newBasis = tmpAx_s + nDir;
            }

            // Check the solution is non-zero
            NekDouble solNorm = Vmath::Dot2(nNonDir,
                                            newBasis,
                                            tmpAx_s + nDir,
                                            m_map + nDir);
            vComm->AllReduce(solNorm, Nektar::LibUtilities::ReduceSum);

            if (solNorm < 22.2 * NekConstants::kNekSparseNonZeroTol)
            {
                return;
            }

            // normalisation of A x
            Vmath::Smul(nNonDir, 1.0/sqrt(solNorm),
                        tmpAx_s + nDir, 1,
                        tmp = tmpAx_s + nDir, 1);

            for (int i = 0; i < m_prevBasis.size(); ++i)
            {
                if(i == insertLocation) continue;
                int skip = i > insertLocation;
                y_s[i-skip] = Vmath::Dot2(nNonDir,
                                         m_prevBasis[i],
                                         tmpAx_s + nDir,
                                         m_map + nDir);
            }
            vComm->AllReduce(y_s, Nektar::LibUtilities::ReduceSum);

            //check if linearly dependent
            DNekMatSharedPtr tilCoeffMatrix;
            if(fullbuffer && m_numSuccessiveRHS>1)
            {
                tilCoeffMatrix = MemoryManager<DNekMat>::
                            AllocateSharedPtr(m_numSuccessiveRHS-1,m_numSuccessiveRHS-1,0.0,eSYMMETRIC);
                for(int i=0; i<m_numSuccessiveRHS; ++i)
                {
                    if(i == insertLocation) continue;
                    int iskip = i >  insertLocation;
                    for(int j=i; j<m_numSuccessiveRHS; ++j)
                    {
                        if(j == insertLocation) continue;
                        int jskip = j >  insertLocation;
                        tilCoeffMatrix->SetValue(i-iskip, j-jskip, m_coeffMatrix->GetValue(i, j));
                    }
                }
            }
            else if(!fullbuffer && m_prevBasis.size())
            {
                tilCoeffMatrix = MemoryManager<DNekMat>::AllocateSharedPtr(m_prevBasis.size(),m_prevBasis.size(),0.0,eSYMMETRIC);
                Vmath::Vcopy(tilCoeffMatrix->GetStorageSize(), m_coeffMatrix->GetPtr(), 1, tilCoeffMatrix->GetPtr(), 1);
            }

            int n, info1 = 0, info2 = 0, info3 = 0;
            if(y_s.size())
            {
                n = tilCoeffMatrix->GetRows();
                Array<OneD, NekDouble> tilCoeffMatrixFactor(tilCoeffMatrix->GetStorageSize());
                Vmath::Vcopy(tilCoeffMatrix->GetStorageSize(), tilCoeffMatrix->GetPtr(), 1, tilCoeffMatrixFactor, 1);
                Lapack::Dsptrf('U', n, tilCoeffMatrixFactor.get(), ipivot.get(), info1);
                if(info1==0)
                {
                    Vmath::Vcopy(n, y_s, 1, invMy_s, 1);
                    Lapack::Dsptrs('U', n, 1, tilCoeffMatrixFactor.get(), ipivot.get(), invMy_s.get(), n, info2);
                }
            }
            if(info1 || info2)
            {
                int latest = ResetKnownSolutionsToLatestOne();
                y_s[0]    = y_s[latest - (latest > insertLocation)];
                invMy_s[0]  = y_s[0];
                insertLocation = m_numPrevSols % m_numSuccessiveRHS;
                fullbuffer = (m_prevBasis.size() == m_numSuccessiveRHS);
            }
            NekDouble residual = 1.;
            NekDouble epsilon = 10. * NekConstants::kNekZeroTol;
            for (int i = 0; i < m_prevBasis.size()-fullbuffer; i++)
            {
                residual  -= y_s[i] * invMy_s[i];
            }
            if(m_verbose && m_root) cout << "SuccessiveRHS: residual " << residual;
            if (residual < epsilon)
            {
                if(m_verbose && m_root) cout << " < " << epsilon << ", reject" << endl;
                return;
            }

            //calculate new coefficient matrix and its factor
            DNekMatSharedPtr newCoeffMatrix;
            if(fullbuffer)
            {
                newCoeffMatrix = MemoryManager<DNekMat>::AllocateSharedPtr(m_numSuccessiveRHS,m_numSuccessiveRHS,0.0,eSYMMETRIC);
                Vmath::Vcopy(m_coeffMatrix->GetStorageSize(), m_coeffMatrix->GetPtr(), 1, newCoeffMatrix->GetPtr(), 1);
                newCoeffMatrix->SetValue(insertLocation, insertLocation, 1.);
                for(int i=0; i<m_numSuccessiveRHS; ++i)
                {
                    if(i == insertLocation) continue;
                    int iskip = i >  insertLocation;
                    newCoeffMatrix->SetValue(insertLocation, i, y_s[i-iskip]);
                }
            }
            else
            {
                newCoeffMatrix = MemoryManager<DNekMat>::AllocateSharedPtr(m_prevBasis.size()+1,m_prevBasis.size()+1,0.0,eSYMMETRIC);
                newCoeffMatrix->SetValue(insertLocation, insertLocation, 1.);
                for(int i=0; i<m_prevBasis.size(); ++i)
                {
                    newCoeffMatrix->SetValue(insertLocation, i, y_s[i]);
                    for(int j=i; j<m_prevBasis.size(); ++j)
                    {
                        newCoeffMatrix->SetValue(i, j, m_coeffMatrix->GetValue(i, j));
                    }
                }
            }
            n = newCoeffMatrix->GetRows();
            Array<OneD, NekDouble> coeffMatrixFactor(newCoeffMatrix->GetStorageSize());
            Vmath::Vcopy(newCoeffMatrix->GetStorageSize(), newCoeffMatrix->GetPtr(), 1, coeffMatrixFactor, 1);
            Lapack::Dsptrf('U', n, coeffMatrixFactor.get(), ipivot.get(), info3);
            if(info3)
            {
                if(m_verbose && m_root) cout << " >= " << epsilon << ", reject (Dsptrf fails)" << endl;
                return;
            }
            if(m_verbose && m_root) cout << " >= " << epsilon << ", accept" << endl;

            //if success, update basis, rhs, coefficient matrix, and its factor
            if(m_prevBasis.size() < m_numSuccessiveRHS)
            {
                m_prevBasis.push_back(tmp = tmpAx_s + nDir);
                if(isAconjugate)
                {
                    m_prevLinSol.push_back(tmp);
                }
                else
                {
                    Array<OneD, NekDouble> solution(nNonDir, 0.0);
                    m_prevLinSol.push_back(solution);
                }
            }
            Vmath::Smul(nNonDir, 1./sqrt(solNorm),
                        tmp = newX + nDir, 1,
                        m_prevLinSol[insertLocation], 1);
            if(!isAconjugate)
            {
                m_prevBasis[insertLocation] = tmpAx_s + nDir;
            }
            m_coeffMatrix = newCoeffMatrix;
            m_coeffMatrixFactor = coeffMatrixFactor;
            m_ipivot = ipivot;
            ++m_numPrevSols;
        }

        void GlobalLinSysIterative::Set_Rhs_Magnitude(
            const NekVector<NekDouble> &pIn)
        {
            Array<OneD, NekDouble> vExchange(1, 0.0);
            if (m_map.size() > 0)
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
