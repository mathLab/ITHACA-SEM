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
#include <MultiRegions/GlobalLinSysStaticCond.h>
#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSysIterativeStaticCond;

        typedef std::shared_ptr<GlobalLinSysIterativeStaticCond>
            GlobalLinSysIterativeStaticCondSharedPtr;

        enum LocalMatrixStorageStrategy
        {
            eNoStrategy,
            eContiguous,
            eNonContiguous,
            eSparse
        };

        const char* const LocalMatrixStorageStrategyMap[] =
        {
            "Contiguous",
            "Non-contiguous",
            "Sparse"
        };


        /// A global linear system.
        class GlobalLinSysIterativeStaticCond : virtual public GlobalLinSysIterative,
                                                virtual public GlobalLinSysStaticCond
        {
        public:
            typedef NekSparseDiagBlkMatrix<StorageSmvBsr<NekDouble> >
                                            DNekSmvBsrDiagBlkMat;
            typedef std::shared_ptr<DNekSmvBsrDiagBlkMat>
                                            DNekSmvBsrDiagBlkMatSharedPtr;

            /// Creates an instance of this class
            static GlobalLinSysSharedPtr create(
                const GlobalLinSysKey                &pLinSysKey,
                const std::weak_ptr<ExpList>         &pExpList,
                const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
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
            MULTI_REGIONS_EXPORT GlobalLinSysIterativeStaticCond(
                const GlobalLinSysKey                &mkey,
                const std::weak_ptr<ExpList>         &pExpList,
                const std::shared_ptr<AssemblyMap>   &locToGloMap);

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterativeStaticCond(
                const GlobalLinSysKey                &mkey,
                const std::weak_ptr<ExpList>         &pExpList,
                const DNekScalBlkMatSharedPtr         pSchurCompl,
                const DNekScalBlkMatSharedPtr         pBinvD,
                const DNekScalBlkMatSharedPtr         pC,
                const DNekScalBlkMatSharedPtr         pInvD,
                const std::shared_ptr<AssemblyMap>   &locToGloMap,
                const PreconditionerSharedPtr         pPrecon);

            virtual ~GlobalLinSysIterativeStaticCond();

        protected:
            virtual DNekScalBlkMatSharedPtr v_GetStaticCondBlock(unsigned int n);
            virtual GlobalLinSysStaticCondSharedPtr v_Recurse(
                const GlobalLinSysKey                &mkey,
                const std::weak_ptr<ExpList>         &pExpList,
                const DNekScalBlkMatSharedPtr         pSchurCompl,
                const DNekScalBlkMatSharedPtr         pBinvD,
                const DNekScalBlkMatSharedPtr         pC,
                const DNekScalBlkMatSharedPtr         pInvD,
                const std::shared_ptr<AssemblyMap>   &locToGloMap);

            void v_PreSolve(int scLevel,
                            Array<OneD, NekDouble>  &F_bnd);
            virtual void v_BasisFwdTransform(
                Array<OneD, NekDouble>& pInOut);
            virtual void v_CoeffsBwdTransform(
                Array<OneD, NekDouble>& pInOut);
            virtual void v_CoeffsFwdTransform(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);

        private:
            /// Dense storage for block Schur complement matrix
            std::vector<double>                      m_storage;
            /// Vector of pointers to local matrix data
            std::vector<const double*>               m_denseBlocks;
            /// Ranks of local matrices
            Array<OneD, unsigned int>                m_rows;
            /// Scaling factors for local matrices
            Array<OneD, NekDouble>                   m_scale;
            /// Sparse representation of Schur complement matrix at this level
            DNekSmvBsrDiagBlkMatSharedPtr            m_sparseSchurCompl;
            /// Utility strings
            static std::string                       storagedef;
            static std::string                       storagelookupIds[];

            virtual void v_InitObject();

            /// Assemble the Schur complement matrix.
            void v_AssembleSchurComplement(
                const std::shared_ptr<AssemblyMap> locToGloMap);

            /// Prepares local representation of Schur complement
            /// stored as a sparse block-diagonal matrix.
            void PrepareLocalSchurComplement();

            /// Perform a Shur-complement matrix multiply operation.
            virtual void v_DoMatrixMultiply(
                    const Array<OneD, NekDouble>& pInput,
                          Array<OneD, NekDouble>& pOutput);

            virtual void v_UniqueMap();
        };
    }
}

#endif
