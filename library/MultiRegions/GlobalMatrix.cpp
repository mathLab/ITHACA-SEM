///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalMatrix.cpp
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
// Description: GlobalMatrix definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalMatrix.h>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>
#include <LibUtilities/LinearAlgebra/SparseMatrix.hpp>
#include <LibUtilities/LinearAlgebra/SparseUtils.hpp>

#include <iomanip>
#include <fstream>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        std::string GlobalMatrix::def = LibUtilities::SessionReader::
            RegisterDefaultSolverInfo("GlobalMatrixStorageType","SmvBSR");
        std::string GlobalMatrix::lookupIds[3] = {
            LibUtilities::SessionReader::RegisterEnumValue(
                "GlobalMatrixStorageType", "SmvBSR", MultiRegions::eSmvBSR)
        };


        /**
         * @class GlobalMatrix
         * This matrix is essentially a wrapper around a DNekSparseMat.
         */

        /**
         * Allocates a new DNekSparseMat object from the given specification.
         * @param   rows        Number of rows in matrix.
         * @param   columns     Number of columns in matrix.
         * @param   cooMat      ?
         */
        GlobalMatrix::GlobalMatrix(
                                   const LibUtilities::SessionReaderSharedPtr& pSession,
                                   unsigned int rows, 
                                   unsigned int columns,
                                   const COOMatType &cooMat,
                                   const MatrixStorage& matStorage):
            m_smvbsrmatrix(),
            m_rows(rows),
            m_mulCallsCounter(0)
        {
            MatrixStorageType storageType = pSession->
                GetSolverInfoAsEnum<MatrixStorageType>("GlobalMatrixStorageType");

            unsigned int brows, bcols;

            // Size of dense matrix sub-blocks
            int block_size = 1;

            BCOMatType bcoMat;

            // assuming current sparse format allows
            // block-sparse data representation

            if(pSession->DefinesParameter("SparseBlockSize"))
            {
                pSession->LoadParameter("SparseBlockSize", block_size);
                ASSERTL1(block_size > 0,"SparseBlockSize parameter must to be positive");
            }

            brows = rows / block_size + (rows % block_size > 0);
            bcols = columns / block_size + (columns % block_size > 0);

            if (rows % block_size > 0)  m_copyOp = true;

            if (m_copyOp)
            {
                m_tmpin  = Array<OneD, NekDouble> (brows*block_size, 0.0);
                m_tmpout = Array<OneD, NekDouble> (brows*block_size, 0.0);
            }

            convertCooToBco(block_size, cooMat, bcoMat);

            size_t matBytes = 0;
            switch(storageType)
            {
                case eSmvBSR:
                    {

                    // Create zero-based Smv-multiply BSR sparse storage holder
                    DNekSmvBsrMat::SparseStorageSharedPtr sparseStorage =
                            MemoryManager<DNekSmvBsrMat::StorageType>::
                                    AllocateSharedPtr(
                                        brows, bcols, block_size, bcoMat, matStorage );

                    // Create sparse matrix
                    m_smvbsrmatrix = MemoryManager<DNekSmvBsrMat>::
                                            AllocateSharedPtr( sparseStorage );

                    matBytes = m_smvbsrmatrix->GetMemoryFootprint();

                    }
                    break;

                default:
                    NEKERROR(ErrorUtil::efatal,"Unsupported sparse storage type chosen");
            }

            cout << "Global matrix storage type: " 
                    << MatrixStorageTypeMap[storageType] << endl;
            std::cout << "Global matrix memory, bytes = " << matBytes;
            if (matBytes/(1024*1024) > 0)
            {
                std::cout << " ("<< matBytes/(1024*1024) <<" MB)" << std::endl;
            }
            else
            {
                std::cout << " ("<< matBytes/1024 <<" KB)" << std::endl;
            }
            std::cout << "Sparse storage block size = " << block_size << std::endl;
        }

        /**
         * Performs a matrix-vector multiply using the sparse format-specific
         * multiply routine.
         * @param   in          Input vector.
         * @param   out         Output vector.
         */
        void GlobalMatrix::Multiply(const Array<OneD,const NekDouble> &in, 
                                          Array<OneD,      NekDouble> &out)
        {
            if (!m_copyOp)
            {
                if (m_smvbsrmatrix)  m_smvbsrmatrix->Multiply(in,out);
            }
            else
            {
                // if block size makes the last row/column bigger, one needs
                // using temporary storage for rhs and result vectors.
                Vmath::Vcopy(m_rows, &in[0], 1, &m_tmpin[0], 1);

                if (m_smvbsrmatrix)  m_smvbsrmatrix->Multiply(m_tmpin,m_tmpout);

                Vmath::Vcopy(m_rows, &m_tmpout[0], 1, &out[0], 1);
            }

            m_mulCallsCounter++;
        }

        unsigned long GlobalMatrix::GetMulCallsCounter() const 
        {
            if (m_smvbsrmatrix)  return m_smvbsrmatrix->GetMulCallsCounter();
            return -1;
        }

        unsigned int GlobalMatrix::GetNumNonZeroEntries() const
        {
            if (m_smvbsrmatrix)  return m_smvbsrmatrix->GetNumNonZeroEntries();
            return -1;
        }


    } //end of namespace
} //end of namespace

