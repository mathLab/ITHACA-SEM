///////////////////////////////////////////////////////////////////////////////
//
// File:  NekLinSysIterFixedpointJacobi.cpp
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
// Description:  NekLinSysIterFixedpointJacobi definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekLinSysIterFixedpointJacobi.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{      
    namespace LibUtilities
    {  
        /**
         * @class  NekLinSysIterFixedpointJacobi
         *
         * Solves a linear system using iterative methods.
         */
        string NekLinSysIterFixedpointJacobi::className =
        LibUtilities::GetNekLinSysIterFactory().RegisterCreatorFunction(
            "FixedpointJacobi", NekLinSysIterFixedpointJacobi::create,
            "NekLinSysIterFixedpointJacobi solver.");

        /// Constructor for full direct matrix solve.
        NekLinSysIterFixedpointJacobi::NekLinSysIterFixedpointJacobi(
            const LibUtilities::SessionReaderSharedPtr  &pSession,
            const LibUtilities::CommSharedPtr           &vComm,
            const int                                   nDimen)
            : NekLinSysIter(pSession, vComm, nDimen)
        {
        }

        void NekLinSysIterFixedpointJacobi::v_InitObject()
        {
            NekLinSysIter::v_InitObject();
        }


        NekLinSysIterFixedpointJacobi::~NekLinSysIterFixedpointJacobi()
        {
        }

        /**
         *
         */
        int NekLinSysIterFixedpointJacobi::v_SolveSystem(
            const int                           nGlobal,
            const Array<OneD, const NekDouble>  &pRhs,
            Array<OneD,      NekDouble>         &pSolution,
            const int                           nDir,
            const NekDouble                     tol,
            const NekDouble                     factor)
        {
            boost::ignore_unused(tol, nDir);

            int niterations = 0;
            m_tolerance     = max(tol, 1.0E-16);
            m_prec_factor   = factor;

            Array<OneD, NekDouble> pSol0 (nGlobal);
            Vmath::Vcopy(nGlobal, pSolution, 1, pSol0, 1);
            for (int i = 0; i < m_maxiter; ++i)
            {
                m_operator.DoNonlinLinFixPointIte(pRhs, pSol0, pSolution);
                Vmath::Vsub(nGlobal, pSolution, 1, pSol0, 1, pSol0, 1);
                m_converged = ConvergenceCheck(i, pSol0, m_tolerance);
                Vmath::Vcopy(nGlobal, pSolution,1, pSol0, 1);
                niterations++;
            }

            return niterations;
        }
    }
}

