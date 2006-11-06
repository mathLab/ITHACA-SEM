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
    // All of the expression interfaces for NekMatrix should go here.
    namespace expt
    {
        template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
        class ExpressionTraits<NekMatrix<DataType, form, BlockType, space> >
        {
            public:
                typedef NekMatrixMetadata MetadataType;
        };
        
        template<typename DataType, NekMatrixForm form, unsigned int space>
        class ExpressionTraits<NekMatrix<DataType, form, eBlock, space> >
        {
            public:
                typedef NekBlockMatrixMetadata MetadataType;
        };
    
        template<typename DataType, NekMatrixForm form, unsigned int space>
        class ExpressionTraits<NekMatrix<DataType, form, ePointerBlock, space> >
        {
            public:
                typedef NekBlockMatrixMetadata MetadataType;
        };
        
        // Binary expression specializations for NekMatrix.
        template<typename DataType, Nektar::NekMatrixForm lhsForm, Nektar::NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
        class BinaryExpressionTraits<Nektar::NekMatrix<DataType, lhsForm, BlockType, space>, Nektar::NekMatrix<DataType, rhsForm, BlockType, space> >
        {
            public:
                typedef Nektar::NekMatrix<DataType, Nektar::eFull, BlockType, space> AdditionResultType;
                typedef Nektar::NekMatrix<DataType, Nektar::eFull, BlockType, space> SubtractionResultType;
                typedef Nektar::NekMatrix<DataType, Nektar::eFull, BlockType, space> DivisionResultType;
                typedef Nektar::NekMatrix<DataType, Nektar::eFull, BlockType, space> MultiplicationResultType;
                
                static const bool AdditionIsAssociative = true;
        };
    
        //// Binary expression specializations for NekMatrix.
        //template<typename DataType, Nektar::LibUtilities::NekMatrixForm rhsForm, unsigned int space>
        //class BinaryExpressionTraits<Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space>, Nektar::LibUtilities::NekMatrix<DataType, rhsForm, space> >
        //{
        //    public:
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> AdditionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> SubtractionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> DivisionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> MultiplicationResultType;
        //};
    
        //// Binary expression specializations for NekMatrix.
        //template<typename DataType, Nektar::LibUtilities::NekMatrixForm lhsForm, unsigned int space>
        //class BinaryExpressionTraits<Nektar::LibUtilities::NekMatrix<DataType, lhsForm, space>, Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> >
        //{
        //    public:
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> AdditionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> SubtractionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> DivisionResultType;
        //        typedef Nektar::LibUtilities::NekMatrix<DataType, Nektar::LibUtilities::eFull, space> MultiplicationResultType;
        //};
    
        // Binary expression specializations for NekMatrix.
        template<typename DataType, MatrixBlockType BlockType, unsigned int space>
        class BinaryExpressionTraits<Nektar::NekMatrix<DataType, Nektar::eDiagonal, BlockType, space>, Nektar::NekMatrix<DataType, Nektar::eDiagonal, BlockType, space> >
        {
            public:
                typedef Nektar::NekMatrix<DataType, Nektar::eDiagonal, BlockType, space> AdditionResultType;
                typedef Nektar::NekMatrix<DataType, Nektar::eDiagonal, BlockType, space> SubtractionResultType;
                typedef Nektar::NekMatrix<DataType, Nektar::eDiagonal, BlockType, space> DivisionResultType;
                typedef Nektar::NekMatrix<DataType, Nektar::eDiagonal, BlockType, space> MultiplicationResultType;
        };
            
            
        template<typename DataType, Nektar::NekMatrixForm lhsForm, unsigned int vectorDim, MatrixBlockType BlockType, unsigned int space>
        class BinaryExpressionTraits<Nektar::NekMatrix<DataType, lhsForm, BlockType, space>, Nektar::NekVector<DataType, vectorDim, space> >
        {
            public:
                typedef Nektar::NekVector<DataType, vectorDim, space> MultiplicationResultType;
        };
    
    }
    
    
    // Now define general purpose operators.
    
    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    void negate(NekMatrix<DataType, form, BlockType, space>& rhs)
    {
        rhs.Negate();
    }
        
    
    /// \brief Defines an addition operator for two matrices.
    ///
    /// Expression template can't use operator+, so we define the add method for those.
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    void add(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
             typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::AdditionResultType& result)
    {
        result = lhs;
        result += rhs;
    }


#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    expt::Expression<expt::BinaryExpressionPolicy<
        expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > >,
        expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, BlockType, space> > >, expt::AddOp > > 
        operator+(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
    {
        typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > > LhsExpressionType;
        typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, BlockType, space> > > RhsExpressionType;

        return expt::Expression<expt::BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, expt::AddOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
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
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    void subtract(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
                const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
                typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::SubtractionResultType& result)
    {
        result = lhs;
        result -= rhs;
    }


    /////////////////////////////////////////
    // Multiplication
    /////////////////////////////////////////
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
                const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
                typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::MultiplicationResultType& result,
                typename boost::enable_if<boost::is_same<DataType, double> >::type* p = NULL )
    {
#ifdef NEKTAR_USING_BLAS
        dgemm(lhs.GetRows(), lhs.GetColumns(), rhs.GetColumns(), lhs.begin(), rhs.begin(), result.begin());
#else
        result = lhs;
        result *= rhs;
#endif
    }

    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
                  const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
                  typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::MultiplicationResultType& result,
                typename boost::disable_if<boost::is_same<DataType, double> >::type* p = NULL )
    {
        result = lhs;
        result *= rhs;
    }

    template<typename DataType, NekMatrixForm lhsForm, unsigned int vectorDim, MatrixBlockType BlockType, unsigned int space>
    void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekVector<DataType, vectorDim, space>& rhs,
                  typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekVector<DataType, vectorDim, space> >::MultiplicationResultType& result)
    {
        ASSERTL0(lhs.GetColumns() == rhs.GetDimension(), "Invalid matrix dimensions in operator*");
        for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
        {
            DataType t = DataType(0);
            for(unsigned int j = 0; j < rhs.GetDimension(); ++j)
            {
                t += lhs(i,j)*rhs(j);
            }
            result(i) = t;
        }
    }

               
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
    template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
    expt::Expression<expt::BinaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > >, expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, BlockType, space> > >, expt::MultiplyOp > > operator*(
            const NekMatrix<DataType, lhsForm, BlockType, space>& lhs,
            const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
    {
        typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > > LhsExpressionType;
        typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, BlockType, space> > > RhsExpressionType;

        return expt::Expression<expt::BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, expt::MultiplyOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
    }

    template<typename DataType, NekMatrixForm lhsForm, unsigned int vectorDim, MatrixBlockType BlockType, unsigned int space>
            expt::Expression<expt::BinaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > >, expt::Expression<expt::ConstantExpressionPolicy<NekVector<DataType, vectorDim, space> > >, expt::MultiplyOp > > operator*(
            const NekMatrix<DataType, lhsForm, BlockType, space>& lhs,
            const NekVector<DataType, vectorDim, space>& rhs)
    {
        typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > > LhsExpressionType;
        typedef expt::Expression<expt::ConstantExpressionPolicy<NekVector<DataType, vectorDim, space> > > RhsExpressionType;

        return expt::Expression<expt::BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, expt::MultiplyOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
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
    NekVector<DataType, 0, BlockType, space> operator*(
        const NekMatrix<DataType, form, BlockType, space>& lhs,
        const NekVector<DataType, 0, BlockType, space>& rhs)
    {
        ASSERTL0(lhs.GetColumns() == rhs.GetDimension(), "Invalid matrix dimensions in operator*");

        NekVector<DataType, 0, BlockType, space> result(rhs.GetDimension(), DataType(0));

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

// #include <loki/Factory.h>
// #include <loki/Singleton.h>
// #include <loki/Typelist.h>
// #include <iostream>
// 
// #include <boost/shared_ptr.hpp>
// #include <boost/call_traits.hpp>
// #include <boost/concept_check.hpp>
// #include <boost/mpl/assert.hpp>
// #include <boost/static_assert.hpp>
// 
// #include <LibUtilities/Memory/NekMemoryManager.hpp>
// #include <LibUtilities/BasicUtils/ErrorUtil.hpp>
// #include <LibUtilities/LinearAlgebra/NekVector.hpp>
// 
// #include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
// #include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>
// #include <LibUtilities/LinearAlgebra/NekMatrixForm.h>
// #include <LibUtilities/LinearAlgebra/NekMatrixImpl.hpp>
// #include <LibUtilities/LinearAlgebra/MatrixBlockType.h>
// 
// // This check will eventually need to go away, but for now it allow people not using
// // cmake to use the NekMatrix header.
// #ifdef NEKTAR_USING_CMAKE
// #include <LibUtilities/LinearAlgebra/blas.h>
// #endif
// 
// namespace Nektar
// {
// 
//    
//     
// 
//     
//     /// NekMatrix Class
//     /// \param DataType The type of data to store in each element of the matrix.
//     /// \param form How is memory allocated for this matrix?  Full is the most straightforward, but can use a lot more
//     ///               memory than needed if most of the elements are 0.
//     /// \param implType Is this a normal or block matrix?
//     template<typename DataType, NekMatrixForm form = eFull, MatrixBlockType BlockType = eNormal, unsigned int space = 0>
//     class NekMatrix
//     {
//         public:
//             typedef NekMatrix<DataType, form, BlockType, space> ThisType;
//             typedef typename MatrixDataType<DataType, form, BlockType, space>::Type BlockDataType;
//             
//             /// The normal version is the friend because that is the matrix type that the block matrices need to initialize.
//             //template<typename A, NekMatrixForm B, MatrixBlockType C, unsigned int D>
//             //friend class NekMatrixImpl<A, B, C, D>;
//             //friend class NekMatrixImpl<DataType, form, eNormal, space>;
//             friend class NekMatrixImpl<DataType, form, eBlock, space>;
//             friend class NekMatrixImpl<DataType, form, eNormal, space>;
//             
//             /// Necessary to allow construction of default constructed matrices.
//             friend class MemoryManager;
//             static const MemoryManager::MemoryPoolEnabler MemoryPoolEnabled = MemoryManager::eEnabled;
//             
//         public:
//             
//             /// \brief Create a matrix with the given size, initialized to the default value of DataType.
//             ///
//             /// It is not possible to create a block matrix with this constructor since the block size is not specified.
//             NekMatrix(unsigned int rows, unsigned int columns);
//             
//             /// It is not possible to create a block matrix with this constructor since the block size is not specified.
//             NekMatrix(unsigned int rows, unsigned int columns, const DataType* const ptr);
//             
//             /// It is not possible to create a block matrix with this constructor since the block size is not specified.
//             NekMatrix(unsigned int rows, unsigned int columns, DataType* ptr, MatrixDataHolderType t = eCopy);
//             
//             /// It is not possible to create a block matrix with this constructor since the block size is not specified.
//             NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d);
//             
//             /// It is not possible to create a block matrix with this constructor since the block size is not specified.
//             NekMatrix(const NekMatrix<DataType, form, BlockType, space>& rhs);
// 
//             /// \brief Constructor for block matrices.
//             NekMatrix(unsigned int blockRows, unsigned int blockColumns, unsigned int* blockRowSizes, unsigned int* blockColumnSizes);
//             
//             template<NekMatrixForm rhsForm>
//             NekMatrix(const NekMatrix<DataType, rhsForm, BlockType, space>& rhs) :
//                 m_rows(rhs.m_rows),
//                 m_columns(rhs.m_columns),
//                 m_data(),
//                 m_dataIsDeletable(eDeletable)
//             {
//             
//                 m_data = NekMatrixImpl<DataType, form, BlockType, space>::CreateMatrixStorage(m_rows, m_columns);
//                 NekMatrixImpl<DataType, form, BlockType, space>::CopyMatrixValues(m_data, rhs, m_rows, m_columns);
//             }
// 
//             template<typename ExpressionPolicyType>
//             NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
//                 m_rows(rhs.GetMetadata().Rows),
//                 m_columns(rhs.GetMetadata().Columns),
//                 m_data(),
//                 m_dataIsDeletable(eDeletable)
//             {
//                 BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, form, BlockType, space> > ));
//                 m_data = NekMatrixImpl<DataType, form, BlockType, space>::CreateMatrixStorage(m_rows, m_columns);
//                 rhs.Apply(*this);
//             }
// 
//             NekMatrix<DataType, form, BlockType, space>& operator=(const NekMatrix<DataType, form, BlockType, space>& rhs);
// 
//             template<NekMatrixForm rhsForm>
//             NekMatrix<DataType, form, BlockType, space>& operator=(const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
//             {
//                 m_rows = rhs.GetRows();
//                 m_columns = rhs.GetColumns();
//                 m_data = NekMatrixImpl<DataType, form, BlockType, space>::CreateMatrixStorage(m_rows, m_columns);
//                 m_dataIsDeletable = eDeletable;
//                 NekMatrixImpl<DataType, form, BlockType, space>::CopyMatrixValues(m_data, rhs, m_rows, m_columns);
//                 return *this;
//             }
// 
//             template<typename ExpressionPolicyType>
//             NekMatrix<DataType, form, BlockType, space>& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
//             {
//                 BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, form, BlockType, space> > ));
//                 m_rows = rhs.GetMetadata().Rows;
//                 m_columns = rhs.GetMetadata().Columns;
//                 m_data = NekMatrixImpl<DataType, form, BlockType, space>::CreateMatrixStorage(m_rows, m_columns);
//                 m_dataIsDeletable = eDeletable;
//                 rhs.Apply(*this);
//                 return *this;
//             }
// 
//             ~NekMatrix();
// 
//             unsigned int GetRows() const;
//             unsigned int GetColumns() const;
// 
//             inline typename boost::call_traits<BlockDataType>::reference operator()
//                 (unsigned int rowNumber, unsigned int colNumber)
//             {
//                 return NekMatrixImpl<DataType, form, BlockType, space>::GetValue(m_data, m_columns, rowNumber, colNumber);
//             }
// 
//             inline typename boost::call_traits<BlockDataType>::const_reference operator()
//                     (unsigned int rowNumber, unsigned int colNumber) const
//             {
//                 return NekMatrixImpl<DataType, form, BlockType, space>::GetValue(m_data, m_columns, rowNumber, colNumber);
//             }
// 
// 
//             typename boost::call_traits<BlockDataType>::reference GetValue
//                     (unsigned int rowNumber, unsigned int colNumber)
//             {
//                 return NekMatrixImpl<DataType, form, BlockType, space>::GetValue(m_data, m_columns, rowNumber, colNumber);
//             }
// 
//             typename boost::call_traits<BlockDataType>::const_reference GetValue
//                     (unsigned int rowNumber, unsigned int colNumber) const
//             {
//                 return NekMatrixImpl<DataType, form, BlockType, space>::GetValue(m_data, m_columns, rowNumber, colNumber);
//             }
// 
//             void SetValue(unsigned int rowNumber, unsigned int colNumber, typename boost::call_traits<BlockDataType>::const_reference rhs)
//             {
//                 NekMatrixImpl<DataType, form, BlockType, space>::SetValue(m_data, m_rows, rowNumber, colNumber, rhs);
//             }
// 
//             BlockDataType* GetPtr(unsigned int rowNumber, unsigned int colNumber)
//             {
//                 return NekMatrixImpl<DataType, form, BlockType, space>::GetPtr(m_data, m_columns, rowNumber, colNumber);
//             }
//             
//             BlockDataType* GetPtr()
//             {
//                 return NekMatrixImpl<DataType, form, BlockType, space>::GetPtr(m_data, m_columns, 0, 0);
//             }
//             
//             const BlockDataType* GetPtr() const
//             {
//                 return NekMatrixImpl<DataType, form, BlockType, space>::GetPtr(m_data, m_columns, 0, 0);
//             }       
//                 
// 
//             typedef BlockDataType* iterator;
//             typedef const BlockDataType* const_iterator;
// 
//             iterator begin();
// 
//             iterator end();
// 
//             const_iterator begin() const;
// 
//             const_iterator end() const;
// 
//             NekMatrixForm GetForm() const;
// 
//             void negate();
//             void Transpose();
//             
//             expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, form, BlockType, space> > >, expt::NegateOp> > operator-() const
//             {
//                 return expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, form, BlockType, space> > >, expt::NegateOp> >(
//                         expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, form, BlockType, space> > >(*this));
//             }
// 
// 
//             NekMatrix<DataType, eFull, BlockType, space> operator+=(const NekMatrix<DataType, eFull, BlockType, space>& rhs)
//             {
//                 ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+");
//                 DataType* lhs_data = begin();
//                 const DataType* rhs_data = rhs.begin();
// 
//                 for(unsigned int i = 0; i < m_rows*m_columns; ++i)
//                 {
//                     lhs_data[i] += rhs_data[i];
//                 }
// 
//                 return *this;
//             }
// 
//             // This is wrong as well.  What if this is diagonal?
//             // enable if on the output.
//             NekMatrix<DataType, eFull, BlockType, space> operator+=(const NekMatrix<DataType, eDiagonal, BlockType, space>& rhs)
//             {
//                 ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+");
// 
//                 for(unsigned int i = 0; i < rhs.GetRows(); ++i)
//                 {
//                     (*this)(i,i) += rhs(i,i);
//                 }
// 
//                 return *this;
//             }
// 
//             // Full *= full = full
//             // Full *= diagonal = full
// 
//             // diag *= diag = diag
//             // diag *= full = full
// 
// 
//             template<NekMatrixForm rhsForm>
//             NekMatrix<DataType, eFull, BlockType, space>& operator*=(const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
//             {
//                 ASSERTL0(GetColumns() == rhs.GetRows(), "Invalid matrix dimensions in operator*");
// 
//                 NekMatrix<DataType, eFull, BlockType, space> result(GetRows(), rhs.GetColumns());
// 
//                 for(unsigned int i = 0; i < result.GetRows(); ++i)
//                 {
//                     for(unsigned int j = 0; j < result.GetColumns(); ++j)
//                     {
//                         DataType t = DataType(0);
// 
//                         // Set the result(i,j) element.
//                         for(unsigned int k = 0; k < GetColumns(); ++k)
//                         {
//                             t += (*this)(i,k)*rhs(k,j);
//                         }
//                         result(i,j) = t;
//                     }
//                 }
// 
//                 Swap(result);
// 
//                 return *this;
//             }
// 
//         private:
//             /// \brief Constructor used for block matrices.  Care must be used with this constructor
//             ///         to initialize the matrix before use.
//             NekMatrix();
//             
//             void Initialize(unsigned int rows, unsigned int columns);
//             
//             void Swap(NekMatrix<DataType, form, BlockType, space>& rhs)
//             {
//                 std::swap(m_rows, rhs.m_rows);
//                 std::swap(m_columns, rhs.m_columns);
//                 std::swap(m_data, rhs.m_data);
//                 std::swap(m_dataIsDeletable, rhs.m_dataIsDeletable);
//             }
// 
//             unsigned int m_rows;
//             unsigned int m_columns;
//             BlockDataType* m_data;
//             MatrixDataDeltetableType m_dataIsDeletable;
//     };
// 
// 
// 
// 
// 
// 
//     template<typename DataType, MatrixBlockType BlockType, unsigned int space>
//     const NekMatrix<DataType, eFull, BlockType, space>& convertToFull(const NekMatrix<DataType, eFull, BlockType, space>& m)
//     {
//         return m;
//     }
// 
//     template<typename DataType, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, eFull, BlockType, space> convertToFull(const NekMatrix<DataType, eDiagonal, BlockType, space>& m)
//     {
//         NekMatrix<DataType, eFull, space> result(m.GetRows(), m.GetColumns(), DataType());
//         for(unsigned int i = 0; i < m.GetRows(); ++i)
//         {
//             result(i,i) = m(i,i);
//         }
// 
//         return result;
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>::NekMatrix(unsigned int rows, unsigned int columns) :
//         m_rows(),
//         m_columns(),
//         m_data(),
//         m_dataIsDeletable(eDeletable)
//     {
//         BOOST_STATIC_ASSERT(BlockType == eNormal);
//         Initialize(rows, columns);
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>::NekMatrix(unsigned int rows, unsigned int columns, const DataType* const ptr) :
//         m_rows(),
//         m_columns(),
//         m_data(),
//         m_dataIsDeletable(eDeletable)
//     {
//         BOOST_STATIC_ASSERT(BlockType == eNormal);
//         Initialize(rows, columns);
//         NekMatrixImpl<DataType, form, BlockType, space>::CopyMatrixValues(m_data, ptr, m_rows, m_columns);
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>::NekMatrix(unsigned int rows, unsigned int columns, DataType* ptr, MatrixDataHolderType t) :
//         m_rows(rows),
//         m_columns(columns),
//         m_data(),
//         m_dataIsDeletable(eNotDeletable)
//     {
//         BOOST_STATIC_ASSERT(BlockType == eNormal);
//         if( t == eCopy )
//         {
//             Initialize(rows, columns);
//             NekMatrixImpl<DataType, form, BlockType, space>::CopyMatrixValues(m_data, ptr, m_rows, m_columns);
//             m_dataIsDeletable = eDeletable;
//         }
//         else
//         {
//             m_data = ptr;
//         }
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>::NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d) :
//         m_rows(),
//         m_columns(),
//         m_data(),
//         m_dataIsDeletable(eDeletable)
//     {
//         BOOST_STATIC_ASSERT(BlockType == eNormal);
//         Initialize(rows, columns);
//         std::fill(begin(), end(), d);
//     }
// 
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>::NekMatrix(const NekMatrix<DataType, form, BlockType, space>& rhs) :
//         m_rows(rhs.m_rows),
//         m_columns(rhs.m_columns),
//         m_data(),
//         m_dataIsDeletable(rhs.m_dataIsDeletable)
//     {
//         if( m_dataIsDeletable == eDeletable )
//         {
//             Initialize(m_rows, m_columns);
//             std::copy(rhs.begin(), rhs.end(), begin());
//         }
//         else
//         {
//             m_data = rhs.m_data;
//         }
//     }
//     
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>& NekMatrix<DataType, form, BlockType, space>::operator=(const NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         NekMatrix<DataType, form, BlockType, space> temp(rhs);
//         Swap(temp);
//         return *this;
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>::NekMatrix() :
//         m_rows(0),
//         m_columns(0),
//         m_data(),
//         m_dataIsDeletable(eDeletable)
//     {
//     }
//     
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>::NekMatrix(unsigned int blockRows, unsigned int blockColumns, unsigned int* blockRowSizes, unsigned int* blockColumnSizes) :
//         m_rows(blockRows),
//         m_columns(blockColumns),
//         m_data(NekMatrixImpl<DataType, form, BlockType, space>::CreateMatrixStorage(blockRows, blockColumns, blockRowSizes, blockColumnSizes)),
//         m_dataIsDeletable(eDeletable)
//     {
//         BOOST_STATIC_ASSERT(BlockType != eNormal);
//     }
//             
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     void NekMatrix<DataType, form, BlockType, space>::Initialize(unsigned int rows, unsigned int columns)
//     {
//         m_rows = rows;
//         m_columns = columns;
//         m_dataIsDeletable = eDeletable;
//         m_data = NekMatrixImpl<DataType, form, BlockType, space>::CreateMatrixStorage(m_rows, m_columns);
//                 
//     }
//     
//     
// 
//     
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space>::~NekMatrix()
//     {
//         if( m_dataIsDeletable == eDeletable )
//         {
//             NekMatrixImpl<DataType, form, BlockType, space>::DeallocateMatrixStorage(m_data, m_rows, m_columns);
//         }
//         m_data = NULL;
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     unsigned int NekMatrix<DataType, form, BlockType, space>::GetRows() const
//     {
//         return m_rows;
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     unsigned int NekMatrix<DataType, form, BlockType, space>::GetColumns() const
//     {
//         return m_columns;
//     }
//     
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     typename NekMatrix<DataType, form, BlockType, space>::iterator NekMatrix<DataType, form, BlockType, space>::begin()
//     {
//         return &m_data[0];
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     typename NekMatrix<DataType, form, BlockType, space>::iterator NekMatrix<DataType, form, BlockType, space>::end()
//     {
//         return NekMatrixImpl<DataType, form, BlockType, space>::end(m_data, m_rows, m_columns);
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     typename NekMatrix<DataType, form, BlockType, space>::const_iterator NekMatrix<DataType, form, BlockType, space>::begin() const
//     {
//         return &m_data[0];
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     typename NekMatrix<DataType, form, BlockType, space>::const_iterator NekMatrix<DataType, form, BlockType, space>::end() const
//     {
//         return NekMatrixImpl<DataType, form, BlockType, space>::end(m_data, m_rows, m_columns);
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrixForm NekMatrix<DataType, form, BlockType, space>::GetForm() const
//     {
//         return form;
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     void NekMatrix<DataType, form, BlockType, space>::negate()
//     {
//         for(iterator iter = begin(); iter != end(); ++iter)
//         {
//             *iter = -*iter;
//         }
//     }
//     
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     void NekMatrix<DataType, form, BlockType, space>::Transpose()
//     {
//         NekMatrixImpl<DataType, form, BlockType, space>::Transpose(m_data, m_rows, m_columns);
//     }
//     

// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     void negate(NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         rhs.negate();
//     }
//     
// 
//     ////////////////////////////////////////
//     // Addition
//     ////////////////////////////////////////
//     
//     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
//     void add(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
//                 const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
//                 typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::AdditionResultType& result)
//     {
//         result = lhs;
//         result += rhs;
//     }
// 
// #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
//     expt::Expression<expt::BinaryExpressionPolicy<
//         expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > >,
//         expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, BlockType, space> > >,
//         expt::AddOp > > operator+(
//         const NekMatrix<DataType, lhsForm, BlockType, space>& lhs,
//         const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
//     {
//         typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > > LhsExpressionType;
//         typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, BlockType, space> > > RhsExpressionType;
// 
//         return expt::Expression<expt::BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, expt::AddOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//     }
// #else
//     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, eFull, BlockType, space> operator+(
//         const NekMatrix<DataType, lhsForm, BlockType, space>& lhs,
//         const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
//     {
//         typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::AdditionResultType result(lhs);
//         result += rhs;
//         return result;
//     }
// 
//     template<typename DataType, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, eFull, BlockType, space> operator+(
//         const NekMatrix<DataType, eFull, BlockType, space>& lhs,
//         const NekMatrix<DataType, eDiagonal, BlockType, space>& rhs)
//     {
//         NekMatrix<DataType, eFull, BlockType, space> result = lhs;
//         result += rhs;
//         return result;
//     }
// 
//     template<typename DataType, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, eFull, BlockType, space> operator+(
//         const NekMatrix<DataType, eDiagonal, BlockType, space>& lhs,
//         const NekMatrix<DataType, eFull, BlockType, space>& rhs)
//     {
//         NekMatrix<DataType, eFull, BlockType, space> result = rhs;
//         result += lhs;
//         return result;
//     }
// #endif
// 
// 
//     /////////////////////////////////////
//     // Subtraction
//     ////////////////////////////////////
//     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
//     void subtract(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
//                 const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
//                 typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::SubtractionResultType& result)
//     {
//         result = lhs;
//         result -= rhs;
//     }
// 
// 
//     /////////////////////////////////////////
//     // Multiplication
//     /////////////////////////////////////////
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
// 
//     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
//     void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
//                   const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
//                   typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::MultiplicationResultType& result,
//                 typename boost::disable_if<boost::is_same<DataType, double> >::type* p = NULL )
//     {
//         result = lhs;
//         result *= rhs;
//     }
// 
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
// 
//                
// #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
//     expt::Expression<expt::BinaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > >, expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, BlockType, space> > >, expt::MultiplyOp > > operator*(
//             const NekMatrix<DataType, lhsForm, BlockType, space>& lhs,
//             const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
//     {
//         typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > > LhsExpressionType;
//         typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, rhsForm, BlockType, space> > > RhsExpressionType;
// 
//         return expt::Expression<expt::BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, expt::MultiplyOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//     }
// 
//     template<typename DataType, NekMatrixForm lhsForm, unsigned int vectorDim, MatrixBlockType BlockType, unsigned int space>
//             expt::Expression<expt::BinaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > >, expt::Expression<expt::ConstantExpressionPolicy<NekVector<DataType, vectorDim, space> > >, expt::MultiplyOp > > operator*(
//             const NekMatrix<DataType, lhsForm, BlockType, space>& lhs,
//             const NekVector<DataType, vectorDim, space>& rhs)
//     {
//         typedef expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, lhsForm, BlockType, space> > > LhsExpressionType;
//         typedef expt::Expression<expt::ConstantExpressionPolicy<NekVector<DataType, vectorDim, space> > > RhsExpressionType;
// 
//         return expt::Expression<expt::BinaryExpressionPolicy<LhsExpressionType, RhsExpressionType, expt::MultiplyOp> >(LhsExpressionType(lhs), RhsExpressionType(rhs));
//     }
// #else
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space> operator*(const NekMatrix<DataType, form, BlockType, space>& lhs,
//                                                const NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         NekMatrix<DataType, form, BlockType, space> result(lhs.GetRows(), rhs.GetColumns());
//         multiply(lhs, rhs, result);
// 
//         return result;
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekVector<DataType, 0, BlockType, space> operator*(
//         const NekMatrix<DataType, form, BlockType, space>& lhs,
//         const NekVector<DataType, 0, BlockType, space>& rhs)
//     {
//         ASSERTL0(lhs.GetColumns() == rhs.GetDimension(), "Invalid matrix dimensions in operator*");
// 
//         NekVector<DataType, 0, BlockType, space> result(rhs.GetDimension(), DataType(0));
// 
//         for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
//         {
//             DataType t = DataType(0);
//             for(unsigned int j = 0; j < rhs.GetDimension(); ++j)
//             {
//                 t += lhs(i,j)*rhs(j);
//             }
//             result(i) = t;
//         }
// 
//         return result;
//     }
// #endif
// 
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     bool operator==(const NekMatrix<DataType, form, BlockType, space>& lhs,
//                     const NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         if( lhs.GetRows() != rhs.GetRows() )
//         {
//             return false;
//         }
// 
//         if( lhs.GetColumns() != rhs.GetColumns() )
//         {
//             return false;
//         }
// 
//         typename NekMatrix<DataType, form, BlockType, space>::const_iterator lhs_iter = lhs.begin();
//         typename NekMatrix<DataType, form, BlockType, space>::const_iterator rhs_iter = rhs.begin();
// 
//         for( ; lhs_iter != lhs.end(); ++lhs_iter, ++rhs_iter )
//         {
//             if( *lhs_iter != *rhs_iter )
//             {
//                 return false;
//             }
//         }
// 
//         return true;
//     }
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     std::ostream& operator<<(std::ostream& os, const NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         for(unsigned int i = 0; i < rhs.GetRows(); ++i)
//         {
//             os << "[";
//             for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
//             {
//                 os << rhs(i,j);
//                 if( j != rhs.GetColumns() - 1 )
//                 {
//                     os << ", ";
//                 }
//             }
//             os << "]";
//             if( i != rhs.GetRows()-1 )
//             {
//                 os << std::endl;
//             }
//         }
//         return os;
//     }
// }

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

/**
    $Log: NekMatrix.hpp,v $
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


