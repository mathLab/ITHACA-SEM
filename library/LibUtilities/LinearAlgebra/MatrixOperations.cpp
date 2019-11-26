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
#include <LibUtilities/LinearAlgebra/ExplicitInstantiation.h>

namespace Nektar
{     
    template<typename ResultDataType, typename LhsDataType, typename LhsMatrixType>
    void Multiply(NekMatrix<ResultDataType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                     const ResultDataType& rhs)
    {
        // TODO - optimize for the different matrix types.
        int n = lhs.GetRows();
        int m = lhs.GetColumns();
        for(unsigned int i = 0; i < n; ++i)
        {
            for(unsigned int j = 0; j < m; ++j)
            {
                result(i,j) = lhs(i,j)*rhs;
            }
        }
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (1, (const NekDouble&)))

    template<typename DataType, typename LhsDataType, typename LhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag>
    Multiply(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const DataType& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType ResultDataType;
        NekMatrix<ResultDataType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        Multiply(result, lhs, rhs);
        return result;
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_ALL_MATRIX_TYPES, (1, (DNekMat)), (0, ()), (1, (const NekDouble&)))

    template<typename RhsDataType, typename RhsMatrixType, typename ResultDataType>
    void Multiply(NekMatrix<ResultDataType, StandardMatrixTag>& result,
                     const ResultDataType& lhs,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs)

    {
        Multiply(result, rhs, lhs);
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (2, (DNekMat&, const NekDouble&)), (0, ()))

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<RhsDataType, RhsMatrixType>::NumberType, StandardMatrixTag>
    Multiply(const DataType& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)

    {
        return Multiply(rhs, lhs);
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_ALL_MATRIX_TYPES, (1, (DNekMat)), (1, (const NekDouble&)), (0, ()))

    template<typename LhsDataType>
    void MultiplyEqual(NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                       typename boost::call_traits<LhsDataType>::const_reference rhs)
    {
        typedef typename NekMatrix<LhsDataType, StandardMatrixTag>::iterator iterator;
        for( iterator iter = lhs.begin(); iter != lhs.end(); ++iter)
        {
            *iter *= rhs;
        }
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(MultiplyEqual, (1, (DNekMat&)), (1, (void)), (0, ()), (1, (const NekDouble&)))

    template<typename ResultType, typename LhsDataType, typename RhsDataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void NekMultiplyDefaultImpl(NekMatrix<ResultType, StandardMatrixTag>& result,
                                         const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                                         const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            std::to_string(lhs.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            std::to_string(rhs.GetRows()) + std::string(" can't be multiplied."));

        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < result.GetColumns(); ++j)
            {
                ResultType t = ResultType(0);

                // Set the result(i,j) element.
                for(unsigned int k = 0; k < lhs.GetColumns(); ++k)
                {
                    t += lhs(i,k)*rhs(k,j);
                }
                result(i,j) = t;
            }
        }
    }
     
    template<typename ResultType, typename LhsDataType, typename RhsDataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void NekMultiplyFullMatrixFullMatrix(NekMatrix<ResultType, StandardMatrixTag>& result,
                                         const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                                         const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        NekMultiplyDefaultImpl(result, lhs, rhs);
    }

    template<typename LhsDataType, typename RhsDataType, typename DataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void Multiply(NekMatrix<DataType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            std::to_string(lhs.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            std::to_string(rhs.GetRows()) + std::string(" can't be multiplied."));

        result.SetSize(lhs.GetRows(), rhs.GetColumns());
        if( lhs.GetType() == eFULL && rhs.GetType() == eFULL)
        {
            NekMultiplyFullMatrixFullMatrix(result, lhs, rhs);
        }
        else
        {
            NekMultiplyDefaultImpl(result, lhs, rhs);
        }
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES(Multiply, NEKTAR_ALL_MATRIX_TYPES, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void AddEqual(NekMatrix<DataType, StandardMatrixTag>& result,
                  const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(result.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            std::to_string(result.GetRows()) + std::string(" and ") +
            std::to_string(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(result.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            std::to_string(result.GetColumns()) + std::string(" and ") +
            std::to_string(rhs.GetColumns()) + std::string(" can't be added."));

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                result(i,j) += rhs(i,j);
            }
        }
    }

    
    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void AddEqualNegatedLhs(NekMatrix<DataType, StandardMatrixTag>& result,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(result.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            std::to_string(result.GetRows()) + std::string(" and ") +
            std::to_string(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(result.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            std::to_string(result.GetColumns()) + std::string(" and ") +
            std::to_string(rhs.GetColumns()) + std::string(" can't be added."));

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                result(i,j) = -result(i,j) + rhs(i,j);
            }
        }
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(AddEqual, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))
    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(AddEqualNegatedLhs, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))


    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void Add(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            std::to_string(lhs.GetRows()) + std::string(" and ") +
            std::to_string(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            std::to_string(lhs.GetColumns()) + std::string(" and ") +
            std::to_string(rhs.GetColumns()) + std::string(" can't be added."));

        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = lhs(i,j) + rhs(i,j);
            }
        }
    }

    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void AddNegatedLhs(NekMatrix<DataType, StandardMatrixTag>& result,
                       const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                       const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            std::to_string(lhs.GetRows()) + std::string(" and ") +
            std::to_string(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            std::to_string(lhs.GetColumns()) + std::string(" and ") +
            std::to_string(rhs.GetColumns()) + std::string(" can't be added."));

        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = -lhs(i,j) + rhs(i,j);
            }
        }
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES(Add, NEKTAR_ALL_MATRIX_TYPES, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))
    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES(AddNegatedLhs, NEKTAR_ALL_MATRIX_TYPES, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))
    
            
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag> 
    Add(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
        const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        Add(result, lhs, rhs);
        return result;
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES(Add, NEKTAR_ALL_MATRIX_TYPES, NEKTAR_ALL_MATRIX_TYPES, (1, (DNekMat)), (0, ()), (0, ()))

    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void Subtract(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            std::to_string(lhs.GetRows()) + std::string(" and ") +
            std::to_string(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            std::to_string(lhs.GetColumns()) + std::string(" and ") +
            std::to_string(rhs.GetColumns()) + std::string(" can't be subtracted."));

        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = lhs(i,j) - rhs(i,j);
            }
        }
    }

    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void SubtractNegatedLhs(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            std::to_string(lhs.GetRows()) + std::string(" and ") +
            std::to_string(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            std::to_string(lhs.GetColumns()) + std::string(" and ") +
            std::to_string(rhs.GetColumns()) + std::string(" can't be subtracted."));

        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = -lhs(i,j) - rhs(i,j);
            }
        }
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES(Subtract, NEKTAR_ALL_MATRIX_TYPES, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))
    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES(SubtractNegatedLhs, NEKTAR_ALL_MATRIX_TYPES, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))

    

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void SubtractEqual(NekMatrix<DataType, StandardMatrixTag>& result,
                          const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(result.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            std::to_string(result.GetRows()) + std::string(" and ") +
            std::to_string(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(result.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            std::to_string(result.GetColumns()) + std::string(" and ") +
            std::to_string(rhs.GetColumns()) + std::string(" can't be subtracted."));

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                result(i,j) -= rhs(i,j);
            }
        }
    }

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void SubtractEqualNegatedLhs(NekMatrix<DataType, StandardMatrixTag>& result,
                          const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(result.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            std::to_string(result.GetRows()) + std::string(" and ") +
            std::to_string(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(result.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            std::to_string(result.GetColumns()) + std::string(" and ") +
            std::to_string(rhs.GetColumns()) + std::string(" can't be subtracted."));

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                result(i,j) = -result(i,j) - rhs(i,j);
            }
        }
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(SubtractEqual, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))
    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(SubtractEqualNegatedLhs, NEKTAR_ALL_MATRIX_TYPES, (1, (void)), (1, (DNekMat&)), (0, ()))
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag> 
    Subtract(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
             const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        Subtract(result, lhs, rhs);
        return result;
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES(Subtract, NEKTAR_ALL_MATRIX_TYPES, NEKTAR_ALL_MATRIX_TYPES, (1, (DNekMat)), (0, ()), (0, ()))
}


