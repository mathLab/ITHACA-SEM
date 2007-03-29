///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrix.hpp
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
// Description: Generic Matrix
//
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekDiagonalMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekFullMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekBlockDiagonalMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekBlockFullMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

namespace Nektar
{
    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space, typename enabled>
    class NekMatrix
    {
        public:
            typedef NekMatrix<DataType, form, BlockType, space, void> ThisType;
            typedef NekMatrixStoragePolicy<DataType, form> StoragePolicy;
            typedef NekMatrixArithmeticPolicy<DataType, form> ArtihmeticPolicy;
            typedef NekMatrixAssignmentPolicy<DataType, form> AssignmentPolicy;

        public:
            /// \brief Creates an empty, 0x0 matrix.  
            NekMatrix() :
                m_rows(0),
                m_columns(0),
                m_data()
            {
            }
    
            /// \brief Create a matrix with the given size, initialized to the default value of DataType.
            ///
            /// It is not possible to create a block matrix with this constructor since the block size is not specified.
            NekMatrix(unsigned int rows, unsigned int columns) :
                m_rows(rows),
                m_columns(columns),
                m_data()
            {
                // If you get a compiler error here, then you are trying to use this constructor with a block matrix.
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                DataType d(0);
                Initialize(rows, columns, d);
            }

            NekMatrix(unsigned int rows, unsigned int columns, DataType* ptr, PointerWrapper t = eCopy) :
                m_rows(rows),
                m_columns(columns),
                m_data()
            {
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                Initialize(rows, columns, ptr, t);
            }

            NekMatrix(unsigned int rows, unsigned int columns, const DataType* const ptr) :
                m_rows(rows),
                m_columns(columns),
                m_data()
            {
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                Initialize(rows, columns, ptr);
            }

            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d) :
                m_rows(),
                m_columns(),
                m_data()
            {
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                Initialize(rows, columns, d);
            }

            NekMatrix(const ThisType& rhs) :
                m_rows(rhs.m_rows),
                m_columns(rhs.m_columns),
                m_data()
            {
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                Initialize(m_rows, m_columns, rhs.m_data.get());
            }

            ThisType& operator=(const ThisType& rhs)
            {
                ThisType temp(rhs);
                Swap(temp);
                return *this;
            }

            template<NekMatrixForm rhs_form>
            ThisType& operator=(const NekMatrix<DataType, rhs_form, BlockType, space>& rhs)
            {
                AssignmentPolicy::Assign(*this, rhs);
                return *this;
            }

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

            template<typename ExpressionPolicyType>
            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                m_rows(rhs.GetMetadata().Rows),
                m_columns(rhs.GetMetadata().Columns),
                m_data()
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, ThisType> ));
                Initialize(m_rows, m_columns);
                rhs.Apply(*this);
            }

            template<typename ExpressionPolicyType>
            NekMatrix& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, ThisType> ));
                m_rows = rhs.GetMetadata().Rows;
                m_columns = rhs.GetMetadata().Columns;
                Initialize(m_rows, m_columns);
                rhs.Apply(*this);
                return *this;
            }
#endif

            ~NekMatrix() 
            {
                m_rows = 0;
                m_columns = 0;
                m_data.reset();
            }


            unsigned int GetRows() const { return m_rows; }
            unsigned int GetColumns() const { return m_columns; }

            typename boost::call_traits<DataType>::reference operator()(unsigned int rowNumber, unsigned int colNumber)
            {
                ASSERTL2(rowNumber < m_rows, "Invalid row number to NekMatrix::operator()");
                ASSERTL2(colNumber < m_columns, "Invalid column number to NekMatrix::operator()");
                return StoragePolicy::GetData(rowNumber, colNumber, m_rows, m_columns, m_data);
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int rowNumber, unsigned int colNumber) const
            {
                ASSERTL2(rowNumber < m_rows, "Invalid row number to NekMatrix::operator()");
                ASSERTL2(colNumber < m_columns, "Invalid column number to NekMatrix::operator()");
                return StoragePolicy::GetConstData(rowNumber, colNumber, m_rows, m_columns, m_data);
            }

            //DataType* GetPtr(unsigned int rowNumber, unsigned int colNumber)
            //{
            //    return m_data.get() + rowNumber*m_columns + colNumber;
            //}

            SharedArray<DataType>& GetPtr()
            {
                return m_data;
            }

            const SharedArray<DataType>& GetPtr() const
            {
                return m_data;
            }

            typedef DataType* iterator;
            typedef const DataType* const_iterator;

            iterator begin() 
            { 
                return m_data.get();
            }

            iterator end() 
            { 
                return m_data.get() + StoragePolicy::NumStorageElements(m_rows, m_columns);
            }

            const_iterator begin() const 
            { 
                return m_data.get();
            }

            const_iterator end() const 
            { 
                return m_data.get() + StoragePolicy::NumStorageElements(m_rows, m_columns);
            }

            void Negate()
            {
                for(iterator iter = begin(); iter != end(); ++iter)
                {
                    *iter = -*iter;
                }
            }
            
            void Transpose()
            {
                StoragePolicy::Transpose(m_data, m_rows, m_columns);
            }

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<ThisType> >, expt::NegateOp> > operator-() const
            {
                return expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<ThisType> >, expt::NegateOp> >(
                        expt::Expression<expt::ConstantExpressionPolicy<ThisType> >(*this));
            }
#endif
            template<NekMatrixForm rhs_form>
            ThisType& operator+=(const NekMatrix<DataType, rhs_form, BlockType, space>& rhs)
            {
                ArtihmeticPolicy::PlusEqual(*this, rhs);
                return *this;
            }

            //// This is wrong as well.  What if this is diagonal?
            //// enable if on the output.
            //NekMatrix<DataType, eFull, eNormal, space> operator+=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            //{
            //    ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");

            //    for(unsigned int i = 0; i < rhs.GetRows(); ++i)
            //    {
            //        (*this)(i,i) += rhs(i, i);
            //    }

            //    return *this;
            //}

            ThisType& operator-=(const ThisType& rhs)
            {
                ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator-=");
                DataType* lhs_data = begin();
                const DataType* rhs_data = rhs.begin();

                for( ; lhs_data < end(); ++lhs_data, ++rhs_data )
                {
                    *lhs_data -= *rhs_data;
                }

                return *this;
            }

            //// This is wrong as well.  What if this is diagonal?
            //// enable if on the output.
            //NekMatrix<DataType, eFull, eNormal, space> operator-=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            //{
            //    ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator-=");

            //    for(unsigned int i = 0; i < rhs.GetRows(); ++i)
            //    {
            //        (*this)(i,i) -= rhs(i, i);
            //    }

            //    return *this;
            //}

            ThisType& operator*=(const ThisType& rhs)
            {
                ASSERTL0(GetColumns() == rhs.GetRows(), "Invalid matrix dimensions in operator*");

                NekMatrix<DataType, eFull, eNormal, space> result(GetRows(), rhs.GetColumns());

                for(unsigned int i = 0; i < result.GetRows(); ++i)
                {
                    for(unsigned int j = 0; j < result.GetColumns(); ++j)
                    {
                        DataType t = DataType(0);

                        // Set the result(i,j) element.
                        for(unsigned int k = 0; k < GetColumns(); ++k)
                        {
                            t += (*this)(i,k)*rhs(k,j);
                        }
                        result(i,j) = t;
                    }
                }

                Swap(result);

                return *this;
            }

            

        private:
            void Initialize(unsigned int rows, unsigned int columns)
            {
                m_rows = rows;
                m_columns = columns;
                unsigned int storageSize = StoragePolicy::NumStorageElements(m_rows, m_columns);
                m_data = MemoryManager::AllocateSharedArray<DataType>(storageSize);
            }
            
            void Initialize(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d)
            {
                Initialize(rows, columns);
                std::fill(begin(), end(), d);
            }

            void Initialize(unsigned int rows, unsigned int columns, DataType* ptr, 
                                   PointerWrapper t)
            {
                m_rows = rows;
                m_columns = columns;
                if( t == eCopy )
                {
                    Initialize(rows, columns, ptr);
                }
                else
                {
                    m_data = SharedArray<DataType>(ptr, StoragePolicy::NumStorageElements(rows, columns),
                        DeleteNothing<DataType>());
                }
            }

            void Initialize(unsigned int rows, unsigned int columns, const DataType* const ptr)
            {
                Initialize(rows, columns);
                std::copy(ptr, ptr + StoragePolicy::NumStorageElements(rows, columns), begin());
            }

            void Swap(NekMatrix<DataType, eFull, eNormal, space>& rhs)
            {
                std::swap(m_rows, rhs.m_rows);
                std::swap(m_columns, rhs.m_columns);
                std::swap(m_data, rhs.m_data);
            }

            unsigned int m_rows;
            unsigned int m_columns;
            SharedArray<DataType> m_data;
    };
    
    
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
    // All of the expression interfaces for NekMatrix should go here.
    namespace expt
    {
        template<typename DataType, Nektar::NekMatrixForm Form, MatrixBlockType BlockType, unsigned int space>
        class ConstantExpressionTraits<Nektar::NekMatrix<DataType, Form, BlockType, space> >
        {
            public:
                typedef Nektar::NekMatrix<DataType, Form, BlockType, space> result_type;
                typedef NekMatrixConstantMetadata MetadataType;
        };
                
        template<typename DataType, Nektar::NekMatrixForm lhsForm, Nektar::NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
        class AdditionTraits<Nektar::NekMatrix<DataType, lhsForm, BlockType, space>, Nektar::NekMatrix<DataType, rhsForm, BlockType, space> >
        {
            public:
                typedef Nektar::NekMatrix<DataType, lhsForm, BlockType, space> LhsType;
                typedef Nektar::NekMatrix<DataType, rhsForm, BlockType, space> RhsType;
                typedef Nektar::NekMatrixAdditionAndSubtractionMetadata MetadataType;
                typedef Nektar::NekMatrix<DataType, Nektar::eFull, BlockType, space> result_type;
                static const bool HasOpEqual = true;
                static const bool HasOpLeftEqual = false;
                
                static void Add(result_type& result, const LhsType& lhs, const RhsType& rhs)
                {
                    result = lhs;
                    result += rhs;
                }
        };
        
        template<typename DataType, Nektar::NekMatrixForm lhsForm, Nektar::NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
        class SubtractionTraits<Nektar::NekMatrix<DataType, lhsForm, BlockType, space>, Nektar::NekMatrix<DataType, rhsForm, BlockType, space> >
        {
            public:
                typedef Nektar::NekMatrix<DataType, lhsForm, BlockType, space> LhsType;
                typedef Nektar::NekMatrix<DataType, rhsForm, BlockType, space> RhsType;
                typedef Nektar::NekMatrixAdditionAndSubtractionMetadata MetadataType;
                typedef Nektar::NekMatrix<DataType, Nektar::eFull, BlockType, space> result_type;
                static const bool HasOpEqual = true;
                static const bool HasOpLeftEqual = false;
                
                static void Subtract(result_type& result, const LhsType& lhs, const RhsType& rhs)
                {
                    result = lhs;
                    result -= rhs;
                }
        };
        
        template<typename DataType, Nektar::NekMatrixForm lhsForm, Nektar::NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
        class MultiplicationTraits<Nektar::NekMatrix<DataType, lhsForm, BlockType, space>, Nektar::NekMatrix<DataType, rhsForm, BlockType, space> >
        {
            public:
                typedef Nektar::NekMatrix<DataType, lhsForm, BlockType, space> LhsType;
                typedef Nektar::NekMatrix<DataType, rhsForm, BlockType, space> RhsType;
                typedef Nektar::NekMatrixMultiplicationMetadata MetadataType;
                typedef Nektar::NekMatrix<DataType, Nektar::eFull, BlockType, space> result_type;
                static const bool HasOpEqual = true;
                static const bool HasOpLeftEqual = false;
                
                static void Multiply(result_type& result, const LhsType& lhs, const RhsType& rhs)
                {
                    result = lhs;
                    result *= rhs;
                }

        };
        
       
        template<typename FirstDataType, typename SecondDataType, 
                 Nektar::NekMatrixForm firstForm, Nektar::NekMatrixForm secondForm,
                 MatrixBlockType firstBlockType, MatrixBlockType secondBlockType, 
                 unsigned int space>
        class CommutativeTraits<NekMatrix<FirstDataType, firstForm, firstBlockType, space>, MultiplyOp,
                                NekMatrix<SecondDataType, secondForm, secondBlockType, space> >
        {
            public:
                static const bool IsCommutative = false;
        };    
    }
#endif
    
    // Now define general purpose operators.
    
    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    void negate(NekMatrix<DataType, form, BlockType, space>& rhs)
    {
        rhs.Negate();
    }
        
    



#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    typename expt::BinaryExpressionType<NekMatrix<DataType, lhsForm, BlockType, space>,
                                         expt::AddOp,
                                         NekMatrix<DataType, rhsForm, BlockType, space> >::Type
    operator+(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
    {
        return expt::CreateBinaryExpression<expt::AddOp>(lhs, rhs);
    }

#else
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::AdditionResultType 
    operator+(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
    {
        typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::AdditionResultType result(lhs);
        result += rhs;
        return result;
    }
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES


    /////////////////////////////////////
    // Subtraction
    ////////////////////////////////////

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    typename expt::BinaryExpressionType<NekMatrix<DataType, lhsForm, BlockType, space>,
                                         expt::SubtractOp,
                                         NekMatrix<DataType, rhsForm, BlockType, space> >::Type
    operator-(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
    {
        return expt::CreateBinaryExpression<expt::SubtractOp>(lhs, rhs);
    }
#else
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::SubtractionResultType 
    operator+(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
    {
        typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::SubtractionResultType result(lhs);
        result -= rhs;
        return result;
    }
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

    /////////////////////////////////////////
    // Multiplication
    /////////////////////////////////////////
//     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
//     void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
//                 const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
//                 typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::MultiplicationResultType& result,
//                 typename boost::enable_if<boost::is_same<DataType, double> >::type* p = NULL )
//     {
// #ifdef NEKTAR_USING_BLAS
//         dgemm(lhs.GetRows(), lhs.GetColumns(), rhs.GetColumns(), lhs.begin(), rhs.begin(), result.begin());
// #else
//         result = lhs;
//         result *= rhs;
// #endif
//     }

//     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
//     void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
//                   const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
//                   typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::MultiplicationResultType& result,
//                 typename boost::disable_if<boost::is_same<DataType, double> >::type* p = NULL )
//     {
//         result = lhs;
//         result *= rhs;
//     }

//     template<typename DataType, NekMatrixForm lhsForm, unsigned int vectorDim, MatrixBlockType BlockType, unsigned int space>
//     void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekVector<DataType, vectorDim, space>& rhs,
//                   typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekVector<DataType, vectorDim, space> >::MultiplicationResultType& result)
//     {
//         ASSERTL0(lhs.GetColumns() == rhs.GetDimension(), "Invalid matrix dimensions in operator*");
//         for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
//         {
//             DataType t = DataType(0);
//             for(unsigned int j = 0; j < rhs.GetDimension(); ++j)
//             {
//                 t += lhs(i,j)*rhs(j);
//             }
//             result(i) = t;
//         }
//     }

               
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    typename expt::BinaryExpressionType<NekMatrix<DataType, lhsForm, BlockType, space>,
                                         expt::MultiplyOp,
                                         NekMatrix<DataType, rhsForm, BlockType, space> >::Type
    operator*(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
    {
        return expt::CreateBinaryExpression<expt::MultiplyOp>(lhs, rhs);
    }

    template<typename DataType, NekMatrixForm lhsForm, unsigned int vectorDim, MatrixBlockType BlockType, unsigned int space>
    typename expt::BinaryExpressionType<NekMatrix<DataType, lhsForm, BlockType, space>,
                                         expt::MultiplyOp,
                                         NekVector<DataType, vectorDim, space> >::Type
    operator*(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekVector<DataType, vectorDim, space>& rhs)
    {
        return expt::CreateBinaryExpression<expt::MultiplyOp>(lhs, rhs);
    }
#else

    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    NekMatrix<DataType, form, BlockType, space> operator*(const NekMatrix<DataType, form, BlockType, space>& lhs,
                                               const NekMatrix<DataType, form, BlockType, space>& rhs)
    {
        NekMatrix<DataType, form, BlockType, space> result(lhs.GetRows(), rhs.GetColumns());
        multiply(lhs, rhs, result);

        return result;
    }

    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    NekVector<DataType, 0, space> operator*(
        const NekMatrix<DataType, form, space>& lhs,
        const NekVector<DataType, 0, space>& rhs)
    {
        ASSERTL0(lhs.GetColumns() == rhs.GetDimension(), "Invalid matrix dimensions in operator*");

        NekVector<DataType, 0, space> result(rhs.GetDimension(), DataType(0));

        for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
        {
            DataType t = DataType(0);
            for(unsigned int j = 0; j < rhs.GetDimension(); ++j)
            {
                t += lhs(i,j)*rhs(j);
            }
            result(i) = t;
        }

        return result;
    }
#endif


    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    bool operator==(const NekMatrix<DataType, form, BlockType, space>& lhs,
                    const NekMatrix<DataType, form, BlockType, space>& rhs)
    {
        if( lhs.GetRows() != rhs.GetRows() )
        {
            return false;
        }

        if( lhs.GetColumns() != rhs.GetColumns() )
        {
            return false;
        }

        typename NekMatrix<DataType, form, BlockType, space>::const_iterator lhs_iter = lhs.begin();
        typename NekMatrix<DataType, form, BlockType, space>::const_iterator rhs_iter = rhs.begin();

        for( ; lhs_iter != lhs.end(); ++lhs_iter, ++rhs_iter )
        {
            if( *lhs_iter != *rhs_iter )
            {
                return false;
            }
        }

        return true;
    }
    
    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    bool operator!=(const NekMatrix<DataType, form, BlockType, space>& lhs,
                     const NekMatrix<DataType, form, BlockType, space>& rhs)
    {
        return !(lhs == rhs);
    }

    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    std::ostream& operator<<(std::ostream& os, const NekMatrix<DataType, form, BlockType, space>& rhs)
    {
        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            os << "[";
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                os << rhs(i,j);
                if( j != rhs.GetColumns() - 1 )
                {
                    os << ", ";
                }
            }
            os << "]";
            if( i != rhs.GetRows()-1 )
            {
                os << std::endl;
            }
        }
        return os;
    }

}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

/**
    $Log: NekMatrix.hpp,v $
    Revision 1.21  2007/01/29 01:31:07  bnelson
    *** empty log message ***

    Revision 1.20  2007/01/23 03:12:50  jfrazier
    Added more conditional compilation directives for expression templates.

    Revision 1.19  2007/01/17 01:11:21  bnelson
    Removed old code.

    Revision 1.18  2007/01/16 05:30:33  bnelson
    Major improvements for expression templates.

    Revision 1.17  2006/11/08 04:16:14  bnelson
    Added subtraction operators.

    Revision 1.16  2006/11/06 17:09:10  bnelson
    *** empty log message ***

    Revision 1.15  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

    Revision 1.14  2006/10/04 03:07:59  bnelson
    Added a check around blas.h to not break existing code.

    Revision 1.13  2006/10/04 03:02:36  bnelson
    Fixed a conflict problem from the previous commit.

    Revision 1.12  2006/10/02 01:16:14  bnelson
    Started working on adding BLAS and LAPACK

    Revision 1.11  2006/09/30 15:18:37  bnelson
    no message

    Revision 1.10  2006/09/21 01:03:31  bnelson
    Added addition and subtraction expression templates.

    Revision 1.9  2006/09/16 23:53:35  bnelson
    Modified the negation operation to reflect changes in the unary expression templates.

    Revision 1.8  2006/09/14 02:06:16  bnelson
    Fixed gcc compiler errors.

    Revision 1.7  2006/09/11 03:26:26  bnelson
    Updated to use new policy based expression templates.

    Revision 1.6  2006/08/28 02:40:21  bnelson
    *** empty log message ***

    Revision 1.5  2006/08/25 03:05:16  bnelson
    Fixed gcc compile errors.

    Revision 1.4  2006/08/25 01:28:53  bnelson
    Changed the way specialized matrices are handled.

    Added support for a matrix wrapper around a pre-allocated array of data.

    Added the NekMemoryManager for allocations.

    Revision 1.3  2006/08/14 02:29:49  bnelson
    Updated points, vectors, and matrix classes to work with ElVis.  Added a variety of methods to all of these classes.

    Revision 1.2  2006/06/01 13:44:28  kirby
    *** empty log message ***

    Revision 1.1  2006/06/01 09:12:41  kirby
    *** empty log message ***

    Revision 1.11  2006/05/31 23:24:21  bnelson
    Moved the diagonal matrix to NekMatrix.hpp.

    Updated method names for the coding standard.

    Revision 1.10  2006/05/31 04:20:16  bnelson
    Changed matrix implementation so the form is a template parameter.

    Revision 1.9  2006/05/29 04:32:18  bnelson
    Changed the data holder to boost::shared_array.

    Revision 1.8  2006/05/29 03:45:04  bnelson
    Updated operator+= to be more efficient using iterators.

    Revision 1.7  2006/05/29 03:40:12  bnelson
    Changed the implementation from individual NekMatrixImpl objects for each type of matrix to an enumeration to store the type.

    Revision 1.6  2006/05/25 03:02:40  bnelson
    Added Matrix/Vector multiplication.

    Revision 1.5  2006/05/18 04:21:06  bnelson
    Removed the ability to specify the rows and columns in the template parameter list.  If this behavior is desired we'll need to create a fixed size array class.

    Added multiplication to the arrays.

    Revision 1.4  2006/05/15 05:06:55  bnelson
    Removed use of MemoryManager pending review of some memory problems.

    Revision 1.3  2006/05/15 04:13:36  bnelson
    no message

    Revision 1.2  2006/05/14 21:32:03  bnelson
 *** empty log message ***

    Revision 1.1  2006/05/04 18:57:43  kirby
 *** empty log message ***

    Revision 1.1  2006/04/11 02:00:43  bnelson
    Initial Revision


 **/


