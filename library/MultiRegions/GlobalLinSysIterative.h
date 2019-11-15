///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterative.h
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
// Description: GlobalLinSysIterative header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVE_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVE_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>

#include <boost/circular_buffer.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;

        /// A global linear system.
        class GlobalLinSysIterative : virtual public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterative(
                const GlobalLinSysKey                &pKey,
                const std::weak_ptr<ExpList>         &pExpList,
                const std::shared_ptr<AssemblyMap>   &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysIterative();

            /// Actual iterative solve-CG
            void DoConjugateGradient(
                const int pNumRows,
                const Array<OneD, const NekDouble> &pInput,
                Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr &locToGloMap,
                const int pNumDir);

        protected:
            /// Global to universal unique map
            Array<OneD, int>                            m_map;

            /// maximum gmres restart iteration
            int                                         m_maxrestart;

            /// maximum gmres search directions for one restart(determines the max storage usage)
            int                                         m_maxstorage;

            /// maximum bandwidth of Hessenburg matrix if use truncted Gmres(m)
            int                                         m_maxhesband;
            /// maximum iterations (for gmres m_maxiter =  m_maxrestart*m_maxdirction)
            int                                         m_maxiter;

            /// Tolerance of iterative solver.
            NekDouble                                   m_tolerance;

            /// dot product of rhs to normalise stopping criterion
            NekDouble                                   m_rhs_magnitude;

            /// dot product of rhs to normalise stopping criterion
            NekDouble                                   m_prec_factor;

            /// cnt to how many times rhs_magnitude is called
            NekDouble                                   m_rhs_mag_sm;

            PreconditionerSharedPtr                     m_precon;

            MultiRegions::PreconditionerType            m_precontype;

            MultiRegions::IterativeMethodType            m_IteraterType;

            int                                         m_totalIterations;

            /// Whether to apply projection technique
            bool                                        m_useProjection;

            /// Whether the iteration has been converged
            bool                                        m_converged;

            /// Root if parallel
            bool                                        m_root;

            /// Storage for solutions to previous linear problems
            boost::circular_buffer<Array<OneD, NekDouble> > m_prevLinSol;

            /// Total counter of previous solutions
            int m_numPrevSols;

            /// A-conjugate projection technique
            void DoAconjugateProjection(
                const int pNumRows,
                const Array<OneD, const NekDouble> &pInput,
                Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr &locToGloMap,
                const int pNumDir);

            void Set_Rhs_Magnitude(const NekVector<NekDouble> &pIn);

            virtual void v_UniqueMap() = 0;

        private:
            void UpdateKnownSolutions(
                const int pGlobalBndDofs,
                const Array<OneD, const NekDouble> &pSolution,
                const int pNumDirBndDofs);

            NekDouble CalculateAnorm(
                const int nGlobal,
                const Array<OneD, const NekDouble> &in,
                const int nDir);

            /// Solve the matrix system
            virtual void v_SolveLinearSystem(
                const int pNumRows,
                const Array<OneD, const NekDouble> &pInput,
                Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr &locToGloMap,
                const int pNumDir);

            virtual void v_DoMatrixMultiply(
                const Array<OneD, NekDouble> &pInput,
                Array<OneD, NekDouble> &pOutput) = 0;

            /// Actual iterative gmres solver for one restart
            NekDouble DoGmresRestart(
                const bool                         restarted,
                const bool                         truncted,
                const int                          nGlobal,
                const Array<OneD, const NekDouble> &pInput,
                Array<OneD,      NekDouble> &pOutput,
                const int                          nDir);
            
            // Arnoldi process
            void DoArnoldi(
                const int starttem,
                const int endtem,
                const int nGlobal,
                const int nDir,
                // V_total(:,1:nd) total search directions
                Array<OneD, Array<OneD,  NekDouble> > &V_local,
                // V[nd] current search direction
                Array<OneD,  NekDouble> &Vsingle1,
                // V[nd+1] new search direction
                Array<OneD, NekDouble> &Vsingle2,
                // One line of Hessenburg matrix
                Array<OneD, NekDouble> &hsingle
            );

            // QR fatorization through Givens rotation
            void DoGivensRotation(
                const int starttem,
                const int endtem,
                const int nGlobal,
                const int nDir,
                Array<OneD, NekDouble> &c,
                Array<OneD, NekDouble> &s,
                Array<OneD, NekDouble> &hsingle,
                Array<OneD, NekDouble> &eta
            );

            // Backward calculation to calculate coeficients of least square problem
            // To notice, Hessenburg's columnns and rows are reverse
            void DoBackward(
                const int  number,
                Array<OneD, Array<OneD, NekDouble> > &A,
                const Array<OneD, const NekDouble> &b,
                Array <OneD, NekDouble> &y
            );
            static std::string lookupIds[];
            static std::string def;
        };
    }
}

#endif
