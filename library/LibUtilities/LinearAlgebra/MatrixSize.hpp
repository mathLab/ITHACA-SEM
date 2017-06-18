///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_LINEAR_ALGEBRA_MATRIX_SIZE_HPP
#define NEKTAR_LIBUTILITIES_LINEAR_ALGEBRA_MATRIX_SIZE_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <type_traits>
#include <tuple>

#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/Operators.hpp>

namespace Nektar
{
    /// \brief Determines the size of a matrix from a matrix expression, as 
    ///         well as the required buffer size for the largest intermediate matrix
    ///         in the expression.
    template<typename NodeType, typename Indices, unsigned int index>
    struct MatrixSize;

    /// \brief Specialization for a constant matrix node.
    template<typename T, typename Indices, unsigned int index>
    struct MatrixSize<expt::Node<T, void, void>, Indices, index>
    {
        static const unsigned int MappedIndex = boost::mpl::at_c<Indices, index>::type::value;

        template<typename R>
        static unsigned int GetRequiredRowsFromMatrix(const R& matrix,
            typename std::enable_if<IsMatrix<R> >::type* dummy = 0)
        {
            return matrix.GetRows();
        }

        template<typename R>
        static unsigned int GetRequiredColumnsFromMatrix(const R& matrix,
            typename std::enable_if<IsMatrix<R> >::type* dummy = 0)
        {
            return matrix.GetColumns();
        }

        template<typename R>
        static unsigned int GetRequiredRowsFromMatrix(const R& matrix,
            typename std::enable_if<IsVector<R> >::type* dummy = 0)
        {
            return matrix.GetRows();
        }

        template<typename R>
        static unsigned int GetRequiredColumnsFromMatrix(const R& matrix,
            typename std::enable_if<IsVector<R> >::type* dummy = 0)
        {
            return 1;
        }

        template<typename R>
        static unsigned int GetRequiredRowsFromMatrix(
            const R& matrix,
            typename std::enable_if<!IsMatrix<R>::value && !IsVector<R>::value>::type* dummy = 0)
        {
            return 1;
        }


        template<typename R>
        static unsigned int GetRequiredColumnsFromMatrix(
            const R& matrix,
            typename std::enable_if<!IsMatrix<R>::value && !IsVector<R>::value>::type* dummy = 0)
        {
            return 1;
        }      
        
        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            unsigned int rows = GetRequiredRowsFromMatrix(boost::fusion::at_c<MappedIndex>(args));
            unsigned int columns = GetRequiredColumnsFromMatrix(boost::fusion::at_c<MappedIndex>(args));

            return std::make_tuple(rows, columns, rows*columns);
        }

        
    };


    template<typename ChildType,
             typename Op,
             typename Indices, unsigned int index>
    struct MatrixSize<expt::Node<ChildType, Op, void>, Indices, index>
    {
        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            return MatrixSize<ChildType, Indices, index>::GetRequiredSize(args);
        }
        
        template<typename ArgumentVectorType>
        static unsigned int
        GetRequiredRows(const ArgumentVectorType& args)
        {
            std::tuple<unsigned int, unsigned int, unsigned int> values = GetRequiredSize(args);
            return values.get<0>();
        }
    };

    
    /// \brief Calculates the size of the given matrix expression, as well as the 
    ///        largest intermediate buffer size, assuming multiplication.
    template<typename LeftNodeType, typename RightNodeType, typename Indices, unsigned int index>
    struct CalculateLargestRequiredSize
    {
        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            std::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
            MatrixSize<LeftNodeType, Indices, index>::GetRequiredSize(args);

            std::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
                MatrixSize<RightNodeType, Indices, index + LeftNodeType::TotalCount>::GetRequiredSize(args);

            unsigned int leftRows = lhsSizes.get<0>();
            unsigned int rightColumns = rhsSizes.get<1>();

            unsigned int matrixSize = leftRows*rightColumns;
            unsigned int bufferSize = std::max(std::max(lhsSizes.get<2>(), rhsSizes.get<2>()), matrixSize);
            
            return std::make_tuple(leftRows, rightColumns, bufferSize);
        }
    };


    template<typename LhsType, typename OpType, typename RhsType, typename Indices, unsigned int index, typename enabled=void>
    struct BinaryMatrixSizeEvaluator;

    // Addition or subtraction.
    template<typename LhsType, typename OpType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, OpType, RhsType, Indices, index,
                                     typename std::enable_if<std::is_same<OpType, expt::AddOp>::value ||
                                                             std::is_same<OpType, expt::SubtractOp>::value
                                                             >::type>
    {
        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            std::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
                MatrixSize<LhsType, Indices, index>::GetRequiredSize(args);
            return lhsSizes;
        }
    };

    // Multiplication with double lhs.
    template<typename LhsType, typename OpType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, OpType, RhsType, Indices, index,
        typename std::enable_if
        <
            expt::IsConstantNode<LhsType>::value &&
            std::is_same<typename LhsType::ResultType, double>::value
        >::type>
    {
        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            std::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
                MatrixSize<RhsType, Indices, index + LhsType::TotalCount>::GetRequiredSize(args);
            return rhsSizes;
        }
    };


    // Double rhs
    template<typename LhsType, typename OpType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, OpType, RhsType, Indices, index,
        typename std::enable_if
        <
            expt::IsConstantNode<RhsType>::value
            std::is_same<typename RhsType::ResultType, double>
        >::type>
    {
        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            std::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
                MatrixSize<LhsType, Indices, index>::GetRequiredSize(args);
            return lhsSizes;
        }
    };

    // Multiplication.
    template<typename LhsType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, expt::MultiplyOp, RhsType, Indices, index,
        typename std::enable_if
        <
            // Test for both nodes being constant matrices.
            // Negate to account for the next specialization.
            !(expt::IsConstantNode<RhsType>::value &&
              !std::is_same<typename RhsType::ResultType, double>::value &&
              expt::IsConstantNode<LhsType>::value &&
              !std::is_same<typename LhsType::ResultType, double>::value)
            && !(
                // Tests to make sure if either is constant, it is a matrix.
                expt::IsConstantNode<RhsType>::value &&
                std::is_same<typename RhsType::ResultType, double>::value)
            && !(
                expt::IsConstantNode<LhsType>::value &&
                std::is_same<typename LhsType::ResultType, double>::value)
            >::type>
    {
        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            return CalculateLargestRequiredSize<LhsType, RhsType, Indices, index>::GetRequiredSize(args);
        }
    };

    template<typename LhsType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, expt::MultiplyOp, RhsType, Indices, index,
        typename std::enable_if
        <
            // Test for both nodes being constant matrices.
            expt::IsConstantNode<RhsType>::value &&
            !std::is_same<typename RhsType::ResultType, double>::value &&
            expt::IsConstantNode<LhsType>::value &&
            !std::is_same<typename LhsType::ResultType, double>
        >::type>
    {
        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            std::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
                MatrixSize<LhsType, Indices, index>::GetRequiredSize(args);

            std::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
                MatrixSize<RhsType, Indices, index + LhsType::TotalCount>::GetRequiredSize(args);

            unsigned int leftRows = lhsSizes.get<0>();
            unsigned int rightColumns = rhsSizes.get<1>();

            unsigned int bufferSize = leftRows*rightColumns;
            
            return std::make_tuple(leftRows, rightColumns, bufferSize);
        }
    };

    template<typename L1, typename LOp, typename L2,
             typename Op,
             typename R1, typename ROp, typename R2,
             typename Indices, unsigned int index>
    struct MatrixSize<expt::Node<expt::Node<L1, LOp, L2>, Op, expt::Node<R1, ROp, R2> >, Indices, index>
    {
        typedef expt::Node<L1, LOp, L2> LhsType;
        typedef expt::Node<R1, ROp, R2> RhsType;

        template<typename ArgumentVectorType>
        static std::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            return BinaryMatrixSizeEvaluator<LhsType, Op, RhsType, Indices, index>::GetRequiredSize(args);
        }
        
        template<typename ArgumentVectorType>
        static unsigned int
        GetRequiredRows(const ArgumentVectorType& args)
        {
            std::tuple<unsigned int, unsigned int, unsigned int> values = GetRequiredSize(args);
            return values.get<0>();
        }
    };

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIBUTILITIES_LINEAR_ALGEBRA_MATRIX_SIZE_HPP

