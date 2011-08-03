///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterativeStaticCond.h
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
// Description: GlobalLinSysIterativeStaticCond header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVESTATICCOND_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVESTATICCOND_H

#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/GlobalLinSysIterative.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class LocalToGlobalC0ContMap;
        class ExpList;
        class GlobalLinSysIterativeStaticCond;

        typedef boost::shared_ptr<GlobalLinSysIterativeStaticCond>
                                        GlobalLinSysIterativeStaticCondSharedPtr;

        /// A global linear system.
        class GlobalLinSysIterativeStaticCond : public GlobalLinSysIterative
        {
        public:
            /// Creates an instance of this class
            static GlobalLinSysSharedPtr create(
                        const GlobalLinSysKey &pLinSysKey,
                        const boost::shared_ptr<ExpList> &pExpList,
                        const boost::shared_ptr<LocalToGlobalBaseMap>
                                                               &pLocToGloMap)
            {
                return MemoryManager<GlobalLinSysIterativeStaticCond>
                    ::AllocateSharedPtr(pLinSysKey, pExpList, pLocToGloMap);
            }

            /// Name of class
            static std::string className;
            static std::string className2;

            /// Constructor for full direct matrix solve.
            GlobalLinSysIterativeStaticCond(
                        const GlobalLinSysKey &mkey,
                        const boost::shared_ptr<ExpList> &pExpList,
                        const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &locToGloMap);

            /// Constructor for full direct matrix solve.
            GlobalLinSysIterativeStaticCond(
                        const GlobalLinSysKey &mkey,
                        const boost::shared_ptr<ExpList> &pExpList,
                        const DNekScalBlkMatSharedPtr pSchurCompl,
                        const DNekScalBlkMatSharedPtr pBinvD,
                        const DNekScalBlkMatSharedPtr pC,
                        const DNekScalBlkMatSharedPtr pInvD,
                        const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &locToGloMap);

            virtual ~GlobalLinSysIterativeStaticCond();

        private:
            /// Schur complement for Direct Static Condensation.
            GlobalLinSysIterativeStaticCondSharedPtr m_recursiveSchurCompl;

            /// Block matrices at this level
            DNekScalBlkMatSharedPtr m_schurCompl;
            DNekScalBlkMatSharedPtr m_BinvD;
            DNekScalBlkMatSharedPtr m_C;
            DNekScalBlkMatSharedPtr m_invD;

            /// Global matrix on which the iterative solver acts.
            DNekMatSharedPtr m_gmat;

            // Local to global map.
            boost::shared_ptr<LocalToGlobalBaseMap>     m_locToGloMap;

            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            virtual void v_Solve(
                        const Array<OneD, const NekDouble>  &in,
                              Array<OneD,       NekDouble>  &out,
                        const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                        const Array<OneD, const NekDouble>  &dirForcing
                                                        = NullNekDouble1DArray);

            /// Initialise this object
            void Initialise(
                    const boost::shared_ptr<LocalToGlobalBaseMap>& locToGloMap);

            /// Set up the storage for the Schur complement or the top level
            /// of the multi-level Schur complement.
            void SetupTopLevel(
                    const boost::shared_ptr<LocalToGlobalBaseMap>& locToGloMap);

            /// Assemble the Schur complement matrix.
            void AssembleSchurComplement(
                    const boost::shared_ptr<LocalToGlobalBaseMap>& locToGloMap);

            ///
            void ConstructNextLevelCondensedSystem(
                    const boost::shared_ptr<LocalToGlobalBaseMap>& locToGloMap);

            /// Compute a diagonal preconditioner of the Shur-complement matrix.
            void ComputeDiagonalPreconditioner(
                    const boost::shared_ptr<LocalToGlobalBaseMap> &pLocToGloMap);

            /// Perform a Shur-complement matrix multiply operation.
            virtual void v_DoMatrixMultiply(
                    const Array<OneD, NekDouble>& pInput,
                          Array<OneD, NekDouble>& pOutput);

            /// Compute the preconditioner for the Shur-complement matrix.
            virtual void v_ComputePreconditioner();
         };
    }
}

#endif
