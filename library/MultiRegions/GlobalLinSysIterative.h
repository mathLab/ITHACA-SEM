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
        class GlobalLinSysIterative : public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterative(
                    const GlobalLinSysKey &pKey,
                    const boost::weak_ptr<ExpList> &pExpList,
                    const boost::shared_ptr<AssemblyMap>
                                                           &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysIterative();

        protected:
            /// Global to universal unique map
            Array<OneD, int>                            m_map;

            /// Tolerance of iterative solver.
            NekDouble                                   m_tolerance;

            PreconditionerSharedPtr                     m_precon;

            MultiRegions::PreconditionerType            m_precontype;


            int                                         m_totalIterations;

            /// Whether to apply projection technique
            bool                                        m_useProjection;

            /// Provide verbose output and root if parallel. 
            bool                                        m_root;
            bool                                        m_verbose;

            /// Storage for solutions to previous linear problems
            boost::circular_buffer<Array<OneD, NekDouble> > m_prevLinSol;

            /// Total counter of previous solutions
            int m_numPrevSols;


            /// A-conjugate projection technique
            void DoAconjugateProjection(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const int pNumDir);

            /// Actual iterative solve
            void DoConjugateGradient(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const int pNumDir);


        private:

            void printArray(
                    const std::string& msg,
                    const Array<OneD, const NekDouble>  &in,
                    const int len,
                    const int offset);


            void UpdateKnownSolutions(
                    const int pGlobalBndDofs,
                    const Array<OneD,const NekDouble> &pSolution,
                    const int pNumDirBndDofs);

            NekDouble CalculateAnorm(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &in,
                    const int nDir);


            /// Solve the matrix system
            virtual void v_SolveLinearSystem(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const int pNumDir);

            virtual void v_DoMatrixMultiply(
                    const Array<OneD, NekDouble>& pInput,
                          Array<OneD, NekDouble>& pOutput) = 0;

            virtual void v_UniqueMap() = 0;
        };
    }
}

#endif
