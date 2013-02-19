///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterativeStaticCond.h
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

#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalLinSysIterative.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSysIterativeStaticCond;

        typedef boost::shared_ptr<GlobalLinSysIterativeStaticCond>
            GlobalLinSysIterativeStaticCondSharedPtr;

        enum LocalMatrixStorageStrategy
        {
            eContiguous,
            eNonContiguous,
        };

        const char* const LocalMatrixStorageStrategyMap[] =
        {
            "Contiguous",
            "Non-contiguous",
        };


        /// A global linear system.
        class GlobalLinSysIterativeStaticCond : public GlobalLinSysIterative
        {
        public:
            /// Creates an instance of this class
            static GlobalLinSysSharedPtr create(
                const GlobalLinSysKey                &pLinSysKey,
                const boost::weak_ptr<ExpList>       &pExpList,
                const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
            {
                GlobalLinSysSharedPtr p = MemoryManager<
                    GlobalLinSysIterativeStaticCond>::
                    AllocateSharedPtr(pLinSysKey, pExpList, pLocToGloMap);
                p->InitObject();
                return p;
            }

            /// Name of class
            static std::string className;
            static std::string className2;

            /// Constructor for full direct matrix solve.
            GlobalLinSysIterativeStaticCond(
                const GlobalLinSysKey                &mkey,
                const boost::weak_ptr<ExpList>       &pExpList,
                const boost::shared_ptr<AssemblyMap> &locToGloMap);

            /// Constructor for full direct matrix solve.
            GlobalLinSysIterativeStaticCond(
                const GlobalLinSysKey                &mkey,
                const boost::weak_ptr<ExpList>       &pExpList,
                const DNekScalBlkMatSharedPtr         pSchurCompl,
                const DNekScalBlkMatSharedPtr         pBinvD,
                const DNekScalBlkMatSharedPtr         pC,
                const DNekScalBlkMatSharedPtr         pInvD,
                const boost::shared_ptr<AssemblyMap> &locToGloMap);

            virtual ~GlobalLinSysIterativeStaticCond();

        protected:
            virtual int v_GetNumBlocks();
            virtual DNekScalBlkMatSharedPtr v_GetStaticCondBlock(unsigned int n);

        private:
            /// Schur complement for Direct Static Condensation.
            GlobalLinSysIterativeStaticCondSharedPtr m_recursiveSchurCompl;
            /// Block Schur complement matrix.
            DNekScalBlkMatSharedPtr                  m_schurCompl;
            /// Dense storage for block Schur complement matrix.
            std::vector<double>                      m_storage;
            /// Vector of pointers to local matrix data
            std::vector<const double*>               m_denseBlocks;
            /// Ranks of local matrices
            Array<OneD, unsigned int>                m_rows;
            /// Scaling factors for local matrices
            Array<OneD, NekDouble>                   m_scale;
            /// Block \f$ BD^{-1} \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_BinvD;
            /// Block \f$ C \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_C;
            /// Block \f$ D^{-1} \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_invD;
            /// Block matrices for low energy
            DNekScalBlkMatSharedPtr                  m_RBlk;
            DNekScalBlkMatSharedPtr                  m_RTBlk;
            DNekScalBlkMatSharedPtr                  m_S1Blk;
            /// Globally assembled Schur complement matrix at this level
            GlobalMatrixSharedPtr                    m_globalSchurCompl;
            /// Local to global map.
            boost::shared_ptr<AssemblyMap>           m_locToGloMap;
            /// Workspace array for matrix multiplication
            Array<OneD, NekDouble>                   m_wsp;
            /// Preconditioner object.
            PreconditionerSharedPtr                  m_precon;

            /// Utility strings
            static std::string                       storagedef;
            static std::string                       storagelookupIds[];


            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            virtual void v_Solve(
                const Array<OneD, const NekDouble> &in,
                      Array<OneD,       NekDouble> &out,
                const AssemblyMapSharedPtr         &locToGloMap,
                const Array<OneD, const NekDouble>  &dirForcing
                    = NullNekDouble1DArray);

            virtual void v_InitObject();

            /// Initialise this object
            void Initialise(
                    const boost::shared_ptr<AssemblyMap>& locToGloMap);

            /// Set up the storage for the Schur complement or the top level
            /// of the multi-level Schur complement.
            void SetupTopLevel(
                    const boost::shared_ptr<AssemblyMap>& locToGloMap);

            void SetupLowEnergyTopLevel(
                    const boost::shared_ptr<AssemblyMap>& locToGloMap);

            /// Assemble the Schur complement matrix.
            void AssembleSchurComplement(
                    const boost::shared_ptr<AssemblyMap>& locToGloMap);

            ///
            void ConstructNextLevelCondensedSystem(
                    const boost::shared_ptr<AssemblyMap>& locToGloMap);

            /// Perform a Shur-complement matrix multiply operation.
            virtual void v_DoMatrixMultiply(
                    const Array<OneD, NekDouble>& pInput,
                          Array<OneD, NekDouble>& pOutput);

            virtual void v_UniqueMap();
        };
    }
}

#endif
