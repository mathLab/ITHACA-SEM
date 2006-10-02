///////////////////////////////////////////////////////////////////////////////
//
// File: blas.cpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/blas.h>

namespace Nektar
{
#if defined(NEKTAR_USING_MKL) || defined(NEKTAR_USING_ATLAS)
    CBLAS_ORDER OrderMapping[] = { CblasRowMajor, CblasColMajor };
    CBLAS_TRANSPOSE TransposeMapping[] = { CblasNoTrans, CblasTrans, CblasConjTrans };
#endif

   void dgemm(const MatrixOrder order, const Transpose MatrixATranspose, const Transpose MatrixBTranspose,
              const int M, const int N, const int K,
              const double alpha, const double* A, const int lda, const double* B,
              const int ldb, const double beta, double* C, const int ldc)
    {
#ifdef NEKTAR_USING_MKL
        cblas_dgemm(OrderMapping[order], TransposeMapping[MatrixATranspose], TransposeMapping[MatrixBTranspose],
                    M, N, K,
                    alpha, A, lda, B, ldb, beta, C, ldc);
#endif

#ifdef NEKTAR_USING_ATLAS
        cblas_dgemm(OrderMapping[order], TransposeMapping[MatrixATranspose], TransposeMapping[MatrixBTranspose],
           M, N, K,
           alpha, A, lda, B, ldb, beta, C, ldc);
#endif

    }
    
    void dgemm(const int rowsInA, const int columnsInA, const int columnsInB, 
               const double* A, const double* B, double* result)
    {
        dgemm(eROW_MAJOR, eNO_TRANSPOSE, eNO_TRANSPOSE, rowsInA, columnsInB, columnsInA,
              1.0, A, columnsInA, B, columnsInB, 0.0, result, columnsInB);
    }
}


/**
    $Log: $

 **/
 
