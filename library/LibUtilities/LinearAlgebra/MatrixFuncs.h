///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixFuncs.h
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
// Description: Matrix functions that depend on storage policy.  Putting
// methods in these separate classes makes it easier to use them from
// normal, scaled, and block matrices.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_FUNCS_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_FUNCS_H

#include <limits>
#include <tuple>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
    struct LIB_UTILITIES_EXPORT BandedMatrixFuncs
    {
        /// \brief Calculates and returns the storage size required.
        ///
        /// This method assumes that the matrix will be used with LU factorizationa and
        /// allocates additional storage as appropriate.
        static unsigned int GetRequiredStorageSize(unsigned int totalRows, unsigned int totalColumns,
                                                   unsigned int subDiags, unsigned int superDiags);

        static unsigned int CalculateNumberOfDiags(unsigned int totalRows, unsigned int diags);

        static unsigned int CalculateNumberOfRows(unsigned int totalRows, unsigned int subDiags, unsigned int superDiags);

        static unsigned int CalculateIndex(unsigned int totalRows,
                                                            unsigned int totalColumns,
                                                            unsigned int row, unsigned int column,
                                                            unsigned int sub, unsigned int super);


        static std::tuple<unsigned int, unsigned int>
        Advance(const unsigned int totalRows, const unsigned int totalColumns,
                const unsigned int curRow, const unsigned int curColumn);
    };

    struct LIB_UTILITIES_EXPORT FullMatrixFuncs
    {

        static unsigned int GetRequiredStorageSize(unsigned int rows, unsigned int columns);
        static unsigned int CalculateIndex(unsigned int totalRows, unsigned int totalColumns, unsigned int curRow, unsigned int curColumn);


        static std::tuple<unsigned int, unsigned int>
        Advance(const unsigned int totalRows, const unsigned int totalColumns,
                const unsigned int curRow, const unsigned int curColumn);

        template<typename DataType>
        static void Invert(unsigned int rows, unsigned int columns,
                           Array<OneD, DataType>& data,
                           const char transpose)
        {
                ASSERTL0(rows==columns, "Only square matrices can be inverted.");
                ASSERTL0(transpose=='N', "Only untransposed matrices may be inverted.");

                int m = rows;
                int n = columns;
                int info = 0;
                Array<OneD, int> ipivot(n);
                Array<OneD, DataType> work(n);

                Lapack::Dgetrf(m, n, data.get(), m, ipivot.get(), info);

                if( info < 0 )
                {
                    std::string message = "ERROR: The " + std::to_string(-info) + "th parameter had an illegal parameter for dgetrf";
                    ASSERTL0(false, message.c_str());
                }
                else if( info > 0 )
                {
                    std::string message = "ERROR: Element u_" + std::to_string(info) +   std::to_string(info) + " is 0 from dgetrf";
                    ASSERTL0(false, message.c_str());
                }

                Lapack::Dgetri(n, data.get(), n, ipivot.get(),
                               work.get(), n, info);

                if( info < 0 )
                {
                    std::string message = "ERROR: The " + std::to_string(-info) + "th parameter had an illegal parameter for dgetri";
                    ASSERTL0(false, message.c_str());
                }
                else if( info > 0 )
                {
                    std::string message = "ERROR: Element u_" + std::to_string(info) +   std::to_string(info) + " is 0 from dgetri";
                    ASSERTL0(false, message.c_str());
                }
        }

        static void EigenSolve(unsigned int n,
                               const Array<OneD, const double>& A,
                               Array<OneD, NekDouble> &EigValReal,
                               Array<OneD, NekDouble> &EigValImag,
                               Array<OneD, NekDouble> &EigVecs = NullNekDouble1DArray)
        {
            int lda = n,info = 0;
            NekDouble dum;
            char uplo = 'N';

            if(EigVecs != NullNekDouble1DArray) // calculate Right Eigen Vectors
            {
                int lwork = 4*lda;
                Array<OneD,NekDouble> work(4*lda);
                char lrev = 'V';
                Lapack::Dgeev(uplo,lrev,lda, A.get(),lda,
                              EigValReal.get(),
                              EigValImag.get(),
                              &dum,1,
                              EigVecs.get(),lda,
                              &work[0],lwork,info);
            }
            else
            {
                int lwork = 3*lda;
                Array<OneD,NekDouble> work(3*lda);
                char lrev = 'N';
                Lapack::Dgeev(uplo,lrev,lda,
                              A.get(),lda,
                              EigValReal.get(),
                              EigValImag.get(),
                              &dum,1,&dum,1,
                              &work[0],lwork,info);
            }
            ASSERTL0(info == 0,"Info is not zero");

        }
    };

    struct LIB_UTILITIES_EXPORT TriangularMatrixFuncs
    {
        static unsigned int GetRequiredStorageSize(unsigned int rows, unsigned int columns);
    };

    struct LIB_UTILITIES_EXPORT UpperTriangularMatrixFuncs : public TriangularMatrixFuncs
    {
        static unsigned int CalculateIndex(unsigned int curRow, unsigned int curColumn);

        static std::tuple<unsigned int, unsigned int>
        Advance(const unsigned int totalRows, const unsigned int totalColumns,
                const unsigned int curRow, const unsigned int curColumn);
    };


    struct LIB_UTILITIES_EXPORT LowerTriangularMatrixFuncs : public TriangularMatrixFuncs
    {
        static unsigned int CalculateIndex(unsigned int totalColumns, unsigned int curRow, unsigned int curColumn);

        static std::tuple<unsigned int, unsigned int>
        Advance(const unsigned int totalRows, const unsigned int totalColumns,
                const unsigned int curRow, const unsigned int curColumn,
                char transpose = 'N');
    };

        /// \internal
    /// Symmetric matrices use upper triangular packed storage.
    struct LIB_UTILITIES_EXPORT SymmetricMatrixFuncs : private TriangularMatrixFuncs
    {
        using TriangularMatrixFuncs::GetRequiredStorageSize;

        static unsigned int CalculateIndex(unsigned int curRow, unsigned int curColumn);

        template<typename DataType>
        static void Invert(unsigned int rows, unsigned int columns,
                           Array<OneD, DataType>& data)
        {
            ASSERTL0(rows==columns, "Only square matrices can be inverted.");

            int n = columns;
            int info = 0;
            Array<OneD, int> ipivot(n);
            Array<OneD, DataType> work(n);

            Lapack::Dsptrf('U', n, data.get(), ipivot.get(), info);

            if( info < 0 )
            {
                std::string message = "ERROR: The " + std::to_string(-info) + "th parameter had an illegal parameter for dsptrf";
                ASSERTL0(false, message.c_str());
            }
            else if( info > 0 )
            {
                std::string message = "ERROR: Element u_" + std::to_string(info) +   std::to_string(info) + " is 0 from dsptrf";
                ASSERTL0(false, message.c_str());
            }

            Lapack::Dsptri('U', n, data.get(), ipivot.get(),
                           work.get(), info);

            if( info < 0 )
            {
                std::string message = "ERROR: The " + std::to_string(-info) + "th parameter had an illegal parameter for dsptri";
                ASSERTL0(false, message.c_str());
            }
            else if( info > 0 )
            {
                std::string message = "ERROR: Element u_" + std::to_string(info) +   std::to_string(info) + " is 0 from dsptri";
                ASSERTL0(false, message.c_str());
            }
        }

        static std::tuple<unsigned int, unsigned int>
        Advance(const unsigned int totalRows, const unsigned int totalColumns,
                const unsigned int curRow, const unsigned int curColumn);
    };

    struct LIB_UTILITIES_EXPORT DiagonalMatrixFuncs
    {
        static std::tuple<unsigned int, unsigned int>
        Advance(const unsigned int totalRows, const unsigned int totalColumns,
                const unsigned int curRow, const unsigned int curColumn);

        template<typename DataType>
        static void Invert(unsigned int rows, unsigned int columns,
                           Array<OneD, DataType>& data)
        {
            ASSERTL0(rows==columns, "Only square matrices can be inverted.");
            for(unsigned int i = 0; i < rows; ++i)
            {
                data[i] = 1.0/data[i];
            }
        }

        static unsigned int GetRequiredStorageSize(unsigned int rows, unsigned int columns);

        static unsigned int CalculateIndex(unsigned int row, unsigned int col);
    };


    struct LIB_UTILITIES_EXPORT TriangularBandedMatrixFuncs
    {
        static unsigned int GetRequiredStorageSize(unsigned int rows, unsigned int columns,
                                                   unsigned int nSubSuperDiags);
    };

    struct LIB_UTILITIES_EXPORT UpperTriangularBandedMatrixFuncs : public TriangularBandedMatrixFuncs
    {
    };

    struct LIB_UTILITIES_EXPORT LowerTriangularBandedMatrixFuncs : public TriangularBandedMatrixFuncs
    {
    };

    /// \internal
    /// Symmetric banded matrices use upper triangular banded packed storage.
    struct LIB_UTILITIES_EXPORT SymmetricBandedMatrixFuncs : private TriangularBandedMatrixFuncs
    {
        using TriangularBandedMatrixFuncs::GetRequiredStorageSize;

        static unsigned int CalculateIndex(unsigned int curRow, unsigned int curColumn,
                                           unsigned int nSuperDiags);
    };

}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_FUNCS_H
