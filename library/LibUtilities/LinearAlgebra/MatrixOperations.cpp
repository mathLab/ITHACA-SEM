///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixOperations.cpp
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
// Description: Defines the global functions needed for matrix operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/MatrixOperations.hpp>

namespace Nektar
{    
#ifdef NEKTAR_USING_BLAS    
    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                            const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& lhs,
                            const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& rhs)
    {
        unsigned int M = lhs.GetRows();
        unsigned int N = rhs.GetColumns();
        unsigned int K = lhs.GetColumns();

        unsigned int LDA = M;
        if( lhs.GetTransposeFlag() == 'T' )
        {
            LDA = K;
        }

        unsigned int LDB = K;
        if( rhs.GetTransposeFlag() == 'T' )
        {
            LDB = N;
        }

        Blas::Dgemm(lhs.GetTransposeFlag(), rhs.GetTransposeFlag(), M, N, K,
            1.0, lhs.GetRawPtr(), LDA, rhs.GetRawPtr(), LDB, 0.0,
            result.GetRawPtr(), lhs.GetRows());
    }
    
    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& lhs,
                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& rhs)
    {
        unsigned int M = lhs.GetRows();
        unsigned int N = rhs.GetColumns();
        unsigned int K = lhs.GetColumns();

        unsigned int LDA = M;
        if( lhs.GetTransposeFlag() == 'T' )
        {
            LDA = K;
        }

        unsigned int LDB = K;
        if( rhs.GetTransposeFlag() == 'T' )
        {
            LDB = N;
        }

        Blas::Dgemm(lhs.GetTransposeFlag(), rhs.GetTransposeFlag(), M, N, K,
            lhs.Scale()*rhs.Scale(), lhs.GetOwnedMatrix()->GetRawPtr(), LDA,
            rhs.GetOwnedMatrix()->GetRawPtr(), LDB, 0.0,
            result.GetRawPtr(), lhs.GetRows());
    }
    
    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& lhs,
                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& rhs)
    {
        unsigned int M = lhs.GetRows();
        unsigned int N = rhs.GetColumns();
        unsigned int K = lhs.GetColumns();

        unsigned int LDA = M;
        if( lhs.GetTransposeFlag() == 'T' )
        {
            LDA = K;
        }

        unsigned int LDB = K;
        if( rhs.GetTransposeFlag() == 'T' )
        {
            LDB = N;
        }

        Blas::Dgemm(lhs.GetTransposeFlag(), rhs.GetTransposeFlag(), M, N, K,
            rhs.Scale(), lhs.GetRawPtr(), LDA, 
            rhs.GetOwnedMatrix()->GetRawPtr(), LDB, 0.0,
            result.GetRawPtr(), lhs.GetRows());
    }

    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& lhs,
                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& rhs)
    {
        unsigned int M = lhs.GetRows();
        unsigned int N = rhs.GetColumns();
        unsigned int K = lhs.GetColumns();

        unsigned int LDA = M;
        if( lhs.GetTransposeFlag() == 'T' )
        {
            LDA = K;
        }

        unsigned int LDB = K;
        if( rhs.GetTransposeFlag() == 'T' )
        {
            LDB = N;
        }

        Blas::Dgemm(lhs.GetTransposeFlag(), rhs.GetTransposeFlag(), M, N, K,
            lhs.Scale(), lhs.GetOwnedMatrix()->GetRawPtr(), LDA, 
            rhs.GetRawPtr(), LDB, 0.0,
            result.GetRawPtr(), lhs.GetRows());
    }
    
    void NekMultiplyEqual(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                          const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& rhs)
    {
        unsigned int M = result.GetRows();
        unsigned int N = rhs.GetColumns();
        unsigned int K = result.GetColumns();

        unsigned int LDA = M;
        if( result.GetTransposeFlag() == 'T' )
        {
            LDA = K;
        }

        unsigned int LDB = K;
        if( rhs.GetTransposeFlag() == 'T' )
        {
            LDB = N;
        }

        Array<OneD, double> buf(result.GetRows()*result.GetColumns());
        Blas::Dgemm(result.GetTransposeFlag(), rhs.GetTransposeFlag(), M, N, K,
            1.0, result.GetRawPtr(), LDA, rhs.GetRawPtr(), LDB, 0.0,
            buf.data(), result.GetRows());
        result = NekMatrix<double, FullMatrixTag, StandardMatrixTag>(result.GetRows(), result.GetColumns(), buf, eWrapper);
    }
#endif //NEKTAR_USING_BLAS
}


