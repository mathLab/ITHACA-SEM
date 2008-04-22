///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixOperations.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

// Since this file defines all of the operations for all combination of matrix types, 
// we have to include all matrix specializations first.
#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <LibUtilities/BasicUtils/RawType.hpp>
#include <boost/utility/enable_if.hpp>

#include <string>

namespace Nektar
{
    ////////////////////////////////////////////////////////////////////////////////////
    // Multiplication
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename ResultDataType, typename LhsDataType, typename LhsStorageType, typename LhsMatrixType>
    void NekMultiply(NekMatrix<ResultDataType, LhsStorageType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>& lhs,
                     const ResultDataType& rhs)
    {
        // TODO - optimize for the different matrix types.
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = lhs(i,j)*rhs;
            }
        }
    }
    
    template<typename DataType, typename LhsDataType, typename LhsStorageType, typename LhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>::NumberType, LhsStorageType, StandardMatrixTag> 
    NekMultiply(const NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>& lhs,
                const DataType& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>::NumberType ResultDataType;
        NekMatrix<ResultDataType, LhsStorageType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekMultiply(result, lhs, rhs);
        return result;
    }

    template<typename RhsDataType, typename RhsStorageType, typename RhsMatrixType, typename ResultDataType>
    void NekMultiply(NekMatrix<ResultDataType, RhsStorageType, StandardMatrixTag>& result,
                     const ResultDataType& lhs,
                     const NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>& rhs)
                     
    {
        NekMultiply(result, rhs, lhs);
    }
    
    template<typename DataType, typename RhsDataType, typename RhsStorageType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>::NumberType, RhsStorageType, StandardMatrixTag> 
    NekMultiply(const DataType& lhs,
                const NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>& rhs)
                     
    {
        return NekMultiply(rhs, lhs);
    }
    
   
    template<typename LhsDataType, typename RhsDataType, typename DataType,
             typename LhsMatrixType, typename RhsMatrixType,
             typename LhsStorageType, typename RhsStorageType>
    void NekMultiply(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>& lhs,
                     const NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>& rhs
                     #ifdef NEKTAR_USING_BLAS
                     // This is necessary to force double matrices to use the blas specific versions below.
                     // Without this check, NekMatrix<double> would not use the versions below because
                     // this is a better match than NekMatrix<const double>
                     , typename boost::disable_if
                     <
                        boost::mpl::and_
                        <
                            boost::is_same<double, typename RawType<typename NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>::NumberType>::type>,
                            boost::is_same<double, typename RawType<typename NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>::NumberType>::type>
                        >
                     >::type* p = 0
                     #endif //NEKTAR_USING_BLAS
                     )
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < result.GetColumns(); ++j)
            {
                DataType t = DataType(0);

                // Set the result(i,j) element.
                for(unsigned int k = 0; k < lhs.GetColumns(); ++k)
                {
                    t += lhs(i,k)*rhs(k,j);
                }
                result(i,j) = t;
            }
        }
    }
        
        
    template<typename LhsDataType, typename MatrixType>
    void NekMultiplyEqual(NekMatrix<LhsDataType, MatrixType, StandardMatrixTag>& lhs,
                          typename boost::call_traits<LhsDataType>::const_reference rhs)
    {
        typedef typename NekMatrix<LhsDataType, MatrixType, StandardMatrixTag>::iterator iterator;
        for( iterator iter = lhs.begin(); iter != lhs.end(); ++iter)
        {
            *iter *= rhs;
        }
    }
    #ifdef NEKTAR_USING_BLAS    
    /// \brief Floating point specialization when blas is in use.
    ///
    /// The signature of this function should match the unspecialized NekMultiply.
    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& lhs,
                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& rhs);

    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& lhs,
                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& rhs);

    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& lhs,
                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& rhs);

    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& lhs,
                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& rhs);
                     
    void NekMultiplyEqual(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                          const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& rhs);
                          
                          
                     
    #endif //NEKTAR_USING_BLAS

	template<typename LhsDataType, typename RhsDataType,
             typename LhsMatrixType, typename RhsMatrixType,
             typename LhsStorageType, typename RhsStorageType>
    NekMatrix<typename boost::remove_const<typename NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>::NumberType>::type, FullMatrixTag, StandardMatrixTag> 
    NekMultiply(const NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>& rhs)
    {
        typedef typename boost::remove_const<typename NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>::NumberType>::type NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), rhs.GetColumns());
        NekMultiply(result, lhs, rhs);
        return result;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Matrix-Vector multiplication
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, FullMatrixTag, MatrixType>& lhs,
                    const NekVector<const DataType, dim, space>& rhs
                    #ifdef NEKTAR_USING_BLAS
                     , typename boost::disable_if
                     <
                        boost::is_same<double, typename RawType<typename NekMatrix<LhsDataType, FullMatrixTag, MatrixType>::NumberType>::type>
                     >::type* p = 0
                     #endif //NEKTAR_USING_BLAS
                     )
    {
        typedef NekMatrix<LhsDataType, FullMatrixTag, MatrixType> LhsType;
        typedef NekVector<const DataType, dim, space> RhsType;
        
//        BOOST_MPL_ASSERT(( boost::mpl::or_
//                           <
//                                boost::is_same<double, typename RawType<typename LhsType::NumberType>::type>,
//                                boost::is_same<BlockMatrixTag, MatrixType> 
//                           > ));
                           
       ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           DataType accum = DataType(0);
           for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
           {
               accum += lhs(i,j)*rhs(j);
           }
           result[i] = accum;
       }
    }

    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, BandedMatrixTag, MatrixType>& lhs,
                    const NekVector<const DataType, dim, space>& rhs
                    #ifdef NEKTAR_USING_BLAS
                     , typename boost::disable_if
                     <
                        boost::is_same<double, typename RawType<typename NekMatrix<LhsDataType, BandedMatrixTag, MatrixType>::NumberType>::type>
                     >::type* p = 0
                     #endif //NEKTAR_USING_BLAS
                     )
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        const typename NekMatrix<LhsDataType, BandedMatrixTag, MatrixType>::PolicySpecificDataHolderType& 
           dataHolder = lhs.GetPolicySpecificDataHolderType();
        unsigned int subDiagonals = dataHolder.GetNumberOfSubDiagonals(lhs.GetRows());
        unsigned int superDiagonals = dataHolder.GetNumberOfSuperDiagonals(lhs.GetRows());
        unsigned int packedRows = subDiagonals+superDiagonals+1;

        const DataType* rawData = lhs.GetRawPtr();
        for(unsigned int i = 0; i < result.GetDimension(); ++i)
        {
            result[i] = DataType(0);
        }
        
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {          
            unsigned int start = i <= superDiagonals ? superDiagonals-i : 0;
            unsigned int end = packedRows;
            if( lhs.GetColumns() - i <= subDiagonals )
            {
                end -= subDiagonals - (lhs.GetColumns() - i) + 1;
            }
            
            unsigned int resultOffset = (i <= superDiagonals) ? 0 : i - superDiagonals;

            for(unsigned int j = start; j < end; ++j, ++resultOffset)
            {
               result[resultOffset] += rawData[i*packedRows+j]*rhs[i];
            }
        }
    }

    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, UpperTriangularMatrixTag, MatrixType>& lhs,
                    const NekVector<const DataType, dim, space>& rhs
                    #ifdef NEKTAR_USING_BLAS
                     , typename boost::disable_if
                     <
                        boost::is_same<double, typename RawType<typename NekMatrix<LhsDataType, UpperTriangularMatrixTag, MatrixType>::NumberType>::type>
                     >::type* p = 0
                     #endif //NEKTAR_USING_BLAS
                    )
    {
       ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           DataType accum = DataType(0);
           for(unsigned int j = i; j < lhs.GetColumns(); ++j)
           {
               accum += lhs(i,j)*rhs(j);
           }
           result[i] = accum;
       }
    }

    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, LowerTriangularMatrixTag, MatrixType>& lhs,
                    const NekVector<const DataType, dim, space>& rhs
                    #ifdef NEKTAR_USING_BLAS
                     , typename boost::disable_if
                     <
                        boost::is_same<double, typename RawType<typename NekMatrix<LhsDataType, LowerTriangularMatrixTag, MatrixType>::NumberType>::type>
                     >::type* p = 0
                     #endif //NEKTAR_USING_BLAS
                     )
    {
       ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           DataType accum = DataType(0);
           for(unsigned int j = 0; j <= i; ++j)
           {
               accum += lhs(i,j)*rhs(j);
           }
           result[i] = accum;
       }
    }

    template<typename DataType, typename LhsInnerMatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                     const NekMatrix<LhsInnerMatrixType, FullMatrixTag, BlockMatrixTag>& lhs,
                     const NekVector<const DataType, dim, space>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
        unsigned int numberOfBlockColumns = lhs.GetNumberOfBlockColumns();

        for(unsigned int i = 0; i < result.GetDimension(); ++i)
        {
            result[i] = DataType(0);
        }
        
        unsigned int curResultRow = 0;
        for(unsigned int blockRow = 0; blockRow < numberOfBlockRows; ++blockRow)
        {
            unsigned int rowsInBlock = lhs.GetNumberOfRowsInBlockRow(blockRow);

            if( blockRow != 0 )
            {
                curResultRow += lhs.GetNumberOfRowsInBlockRow(blockRow-1);
            }

            NekVector<DataType, VariableSizedVector, space> resultWrapper(rowsInBlock, result.GetPtr() + curResultRow, eWrapper);

            unsigned int curWrapperRow = 0;
            for(unsigned int blockColumn = 0; blockColumn < numberOfBlockColumns; ++blockColumn)
            {
                if( blockColumn != 0 )
                {
                    curWrapperRow += lhs.GetNumberOfColumnsInBlockColumn(blockColumn-1);
                }

                const boost::shared_ptr<const LhsInnerMatrixType>& block = lhs.GetBlock(blockRow, blockColumn);
                if( !block )
                {
                    continue;
                }

                unsigned int columnsInBlock = lhs.GetNumberOfColumnsInBlockColumn(blockColumn);
                NekVector<const DataType, VariableSizedVector, space> rhsWrapper(columnsInBlock, rhs.GetPtr() + curWrapperRow, eWrapper);
                
                resultWrapper += (*block)*rhsWrapper;
            }
        }
    }

    #ifdef NEKTAR_USING_BLAS  
        /// \brief Floating point specialization when blas is in use.
        template<typename dim, typename space>
        void NekMultiply(NekVector<double, VariableSizedVector, space>& result,
                        const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& lhs,
                        const NekVector<const double, dim, space>& rhs)
        {
            ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
               boost::lexical_cast<std::string>(lhs.GetColumns()) + 
               std::string(" and a right side vector with row count ") + 
               boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

            int m = lhs.GetRows();
            int n = lhs.GetColumns();

            char t = lhs.GetTransposeFlag();
            if( t == 'T' )
            {
                std::swap(m, n);
            }

            double alpha = 1.0;
            const double* a = lhs.GetRawPtr();
            int lda = m;
            const double* x = rhs.GetRawPtr();
            int incx = 1;
            double beta = 0.0;
            double* y = result.GetRawPtr();
            int incy = 1;
            
            Blas::Dgemv(t, m, n, alpha, a, lda, x, incx, beta, y, incy);
        }
               
        template<typename dim, typename space>
        void NekMultiply(NekVector<double, VariableSizedVector, space>& result,
                        const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& lhs,
                        const NekVector<const double, dim, space>& rhs)
        {
            ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
               boost::lexical_cast<std::string>(lhs.GetColumns()) + 
               std::string(" and a right side vector with row count ") + 
               boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

            char transpose = lhs.GetTransposeFlag();
            int m = lhs.GetRows();
            int n = lhs.GetColumns();

            if( transpose == 'T' )
            {
                std::swap(m, n);
            }
            
            double alpha = lhs.Scale();
            const double* a = lhs.GetOwnedMatrix()->GetRawPtr();
            int lda = m;
            const double* x = rhs.GetRawPtr();
            int incx = 1;
            double beta = 0.0;
            double* y = result.GetRawPtr();
            int incy = 1;
            
            Blas::Dgemv(transpose, m, n, alpha, a, lda, x, incx, beta, y, incy);
        }

        template<typename dim, typename space>
        void NekMultiply(NekVector<double, VariableSizedVector, space>& result,
                         const NekMatrix<const double, UpperTriangularMatrixTag, StandardMatrixTag>& lhs,
                         const NekVector<const double, dim, space>& rhs)
        {
            ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
               boost::lexical_cast<std::string>(lhs.GetColumns()) + 
               std::string(" and a right side vector with row count ") + 
               boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

            result = rhs;
            int n = lhs.GetRows();
            const double* a = lhs.GetRawPtr();
            double* x = result.GetRawPtr();
            int incx = 1;
            
            Blas::Dtpmv('U', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
        }

        template<typename dim, typename space>
        void NekMultiply(NekVector<double, VariableSizedVector, space>& result,
                         const NekMatrix<const double, LowerTriangularMatrixTag, StandardMatrixTag>& lhs,
                         const NekVector<const double, dim, space>& rhs)
        {
            ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
               boost::lexical_cast<std::string>(lhs.GetColumns()) + 
               std::string(" and a right side vector with row count ") + 
               boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

            result = rhs;
            int n = lhs.GetRows();
            const double* a = lhs.GetRawPtr();
            double* x = result.GetRawPtr();
            int incx = 1;
            
            Blas::Dtpmv('L', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
        }

        template<typename LhsDataType, typename MatrixType, typename dim, typename space>
        void NekMultiply(NekVector<double, VariableSizedVector, space>& result,
                        const NekMatrix<LhsDataType, BandedMatrixTag, MatrixType>& lhs,
                        const NekVector<const double, dim, space>& rhs)
        {
            ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
               boost::lexical_cast<std::string>(lhs.GetColumns()) + 
               std::string(" and a right side vector with row count ") + 
               boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

            const typename NekMatrix<LhsDataType, BandedMatrixTag, MatrixType>::PolicySpecificDataHolderType& 
                dataHolder = lhs.GetPolicySpecificDataHolderType();
       
            int m = lhs.GetRows();
            int n = lhs.GetColumns();
            int kl = dataHolder.GetNumberOfSubDiagonals(lhs.GetRows());
            int ku = dataHolder.GetNumberOfSuperDiagonals(lhs.GetRows());
            double alpha = 1.0;
            const double* a = lhs.GetRawPtr();
            int lda = kl + ku + 1;
            const double* x = rhs.GetRawPtr();
            int incx = 1;
            double beta = 0.0;
            double* y = result.GetRawPtr();
            int incy = 1;
            Blas::Dgbmv('N', m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);

        }
    #endif //NEKTAR_USING_BLAS

    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, DiagonalMatrixTag, MatrixType>& lhs,
                    const NekVector<DataType, dim, space>& rhs)
    {
       ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           result[i] = lhs(i,i)*rhs(i);
       }
    }
   
    /// We need two version of NekMultiply, one for constant sized vectors and one for 
    /// variable sized vectors because the constructors to initialize the vector to 0 are 
    /// different in each case.
    ///
    /// Note that I can't make a compile time decision about the size of the result vector
    /// because the matrix dimensions are runtime only.
    template<typename DataType, typename LhsDataType, typename StorageType, typename MatrixType, typename dim, typename space>
    NekVector<typename boost::remove_const<DataType>::type, VariableSizedVector, space> 
    NekMultiply(const NekMatrix<LhsDataType, StorageType, MatrixType>& lhs,
                const NekVector<DataType, dim, space>& rhs)
    {
       NekVector<typename boost::remove_const<DataType>::type, VariableSizedVector, space> result(lhs.GetRows(), DataType(0));
       NekMultiply(result, lhs, rhs);
       return result;
    }

    template<typename DataType, typename LhsDataType, typename StorageType, typename MatrixType, typename space>
    NekVector<typename boost::remove_const<DataType>::type, VariableSizedVector, space> 
    NekMultiply(const NekMatrix<LhsDataType, StorageType, MatrixType>& lhs,
                const NekVector<DataType, VariableSizedVector, space>& rhs)
    {
       NekVector<typename boost::remove_const<DataType>::type, VariableSizedVector, space> result(lhs.GetRows(), DataType(0));
       NekMultiply(result, lhs, rhs);
       return result;
    }

    ///////////////////////////////////////////////////////////////////
    // Addition
    ///////////////////////////////////////////////////////////////////
    template<typename ResultType, typename RhsType>
    void AddEqualMatricesWithSameStorageType(ResultType& result, const RhsType& rhs)
    {
        //typename ResultType::NumberType* r = result.GetRawPtr();
        //const typename RhsType::NumberType* rhs_buf = rhs.GetRawPtr();
        //unsigned int numElements = result.GetStorageSize();
        //for(unsigned int i = 0; i < numElements; ++i)
        //{
        //    r[i] += rhs_buf[i];
        //}
        typename ResultType::iterator result_iter = result.begin();
        typename RhsType::const_iterator rhs_iter = rhs.begin();
        while( result_iter != result.end() )
        {
            (*result_iter) += (*rhs_iter);
            ++result_iter;
            ++rhs_iter;
        }
    }

    
    template<typename ResultType, typename LhsType, typename RhsType>
    void AddMatricesWithSameStorageType(ResultType& result, const LhsType& lhs, const RhsType& rhs)
    {
        //if( lhs.GetTransposeFlag() == rhs.GetTransposeFlag() )
        //{
        //    // This is much more efficient if we don't need to worry about stepping through 
        //    // the data in different ways for each operand.
        //    typename ResultType::NumberType* r = result.GetRawPtr();
        //    const typename LhsType::NumberType* lhs_buf = lhs.GetRawPtr();
        //    const typename RhsType::NumberType* rhs_buf = rhs.GetRawPtr();
        //    unsigned int end = result.GetStorageSize();
        //    for(unsigned int i = 0; i < end; ++i)
        //    {
        //        r[i] = lhs_buf[i] + rhs_buf[i];
        //    }
        //}
        //else
        //{
            typename ResultType::iterator result_iter = result.begin();
            typename RhsType::const_iterator rhs_iter = rhs.begin();
            typename LhsType::const_iterator lhs_iter = lhs.begin();
            while( result_iter != result.end() )
            {
                (*result_iter) = (*lhs_iter) + (*rhs_iter);
                ++result_iter;
                ++rhs_iter;
                ++lhs_iter;
            } 
        //}
    }


    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void NekAddEqual(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(result.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(result.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(result.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(result.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));

        return AddEqualMatricesWithSameStorageType(result, rhs);
    }



    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));

        AddMatricesWithSameStorageType(result, lhs, rhs);
    }

    
            
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
    NekAdd(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
           const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekAdd(result, lhs, rhs);
        return result;
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
        
        result = lhs;
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) += rhs(i,i);
        }
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
    NekAdd(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
           const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekAdd(result, lhs, rhs);
        return result;
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));

        result = rhs;
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) += lhs(i,i);
        }
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
    NekAdd(const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
           const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekAdd(result, lhs, rhs);
        return result;
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));

        result = lhs;
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = i; j < result.GetColumns(); ++j)
            {
                result(i,i) += lhs(i,i);
            }
        }
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
    NekAdd(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
           const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekAdd(result, lhs, rhs);
        return result;
    }

    // UT = UT + UT
    // Line 22
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, UpperTriangularMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));

        return AddMatricesWithSameStorageType(result, lhs, rhs);
    }

    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>::NumberType, UpperTriangularMatrixTag, StandardMatrixTag> 
    NekAdd(const NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>& lhs,
           const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, UpperTriangularMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekAdd(result, lhs, rhs);
        return result;
    }

    // LT = LT + LT
    // Line 32
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, LowerTriangularMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, LowerTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
            
        return AddMatricesWithSameStorageType(result, lhs, rhs);
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>::NumberType, LowerTriangularMatrixTag, StandardMatrixTag> 
    NekAdd(const NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>& lhs,
           const NekMatrix<RhsDataType, LowerTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, LowerTriangularMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekAdd(result, lhs, rhs);
        return result;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // Subtraction
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename ResultType, typename RhsType>
    void SubtractEqualMatricesWithSameStorageType(ResultType& result, const RhsType& rhs)
    {
        typename ResultType::iterator result_iter = result.begin();
        typename RhsType::const_iterator rhs_iter = rhs.begin();
        while( result_iter != result.end() )
        {
            (*result_iter) -= (*rhs_iter);
            ++result_iter;
            ++rhs_iter;
        }
    }


    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));

        typename NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>::iterator result_iter = result.begin();
        typename NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>::const_iterator rhs_iter = rhs.begin();
        typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::const_iterator lhs_iter = lhs.begin();
        while( result_iter != result.end() )
        {
            (*result_iter) = (*lhs_iter) - (*rhs_iter);
            ++result_iter;
            ++rhs_iter;
            ++lhs_iter;
        }
    }

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtractEqual(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                          const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(result.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(result.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(result.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(result.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));

        SubtractEqualMatricesWithSameStorageType(result, rhs);
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
    NekSubtract(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekSubtract(result, lhs, rhs);
        return result;
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));
        
        result = lhs;
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) -= rhs(i,i);
        }
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
    NekSubtract(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekSubtract(result, lhs, rhs);
        return result;
    }

    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));
        
        result = rhs;
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) -= lhs(i,i);
        }
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
    NekSubtract(const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColuomns());
        NekSubtract(result, lhs, rhs);
        return result;
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));

        result = rhs;
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = i; j < result.GetColumns(); ++j)
            {
                result(i,i) -= lhs(i,i);
            }
        }
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
    NekSubtract(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColumns());
        NekSubtract(result, lhs, rhs);
        return result;
    }


    
    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 3, NekMatrix, 3);
    
//    // TODO - Any update possible for this case?  The constant value must be the same as the matrix 
//    // number type.  Seems pretty custom to me.
//    template<typename DataType, typename LhsDataType, typename StorageType, typename MatrixType>
//    typename MultiplicationTraits<NekMatrix<LhsDataType, StorageType, MatrixType>, DataType>::ResultType
//    operator*(const NekMatrix<LhsDataType, StorageType, MatrixType>& lhs, const DataType& rhs)
//    {
//        return NekMultiply(lhs, rhs);
//    }
//                        
//    template<typename DataType, typename RhsDataType, typename StorageType, typename MatrixType>
//    typename MultiplicationTraits<DataType, NekMatrix<RhsDataType, StorageType, MatrixType> >::ResultType
//    operator*(const DataType& lhs, const NekMatrix<RhsDataType, StorageType, MatrixType>& rhs)
//    {
//        return NekMultiply(lhs, rhs);
//    }
    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 3, NekDouble, 0);
    GENERATE_MULTIPLICATION_OPERATOR(NekDouble, 0, NekMatrix, 3);
    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 3, NekVector, 3);
    
    GENERATE_DIVISION_OPERATOR(NekMatrix, 3, NekMatrix, 3);
    GENERATE_ADDITION_OPERATOR(NekMatrix, 3, NekMatrix, 3);
    GENERATE_SUBTRACTION_OPERATOR(NekMatrix, 3, NekMatrix, 3);
    
    // TODO - Either update the GENERATE macros to allow non-type template parameters,
    // or put an if/else here for expression templates.
//    template<typename MatrixDataType, typename StorageType, typename MatrixType, 
//             typename VectorDataType, typename dim, typename space>
//    typename MultiplicationTraits<NekMatrix<MatrixDataType, StorageType, MatrixType>, NekVector<VectorDataType, dim, space> >::ResultType
//    operator*(const NekMatrix<MatrixDataType, StorageType, MatrixType>& lhs,
//              const NekVector<VectorDataType, dim, space>& rhs)
//    {
//        return NekMultiply(lhs, rhs);
//    }          
}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

