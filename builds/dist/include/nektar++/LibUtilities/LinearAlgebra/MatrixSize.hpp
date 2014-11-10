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

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/tuple/tuple.hpp>

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
            typename boost::enable_if<IsMatrix<R> >::type* dummy = 0)
        {
            return matrix.GetRows();
        }

        template<typename R>
        static unsigned int GetRequiredColumnsFromMatrix(const R& matrix,
            typename boost::enable_if<IsMatrix<R> >::type* dummy = 0)
        {
            return matrix.GetColumns();
        }

        template<typename R>
        static unsigned int GetRequiredRowsFromMatrix(const R& matrix,
            typename boost::enable_if<IsVector<R> >::type* dummy = 0)
        {
            return matrix.GetRows();
        }

        template<typename R>
        static unsigned int GetRequiredColumnsFromMatrix(const R& matrix,
            typename boost::enable_if<IsVector<R> >::type* dummy = 0)
        {
            return 1;
        }

        template<typename R>
        static unsigned int GetRequiredRowsFromMatrix(const R& matrix,
            typename boost::disable_if
            <
                boost::mpl::or_
                <
                    IsMatrix<R> ,
                    IsVector<R> 
                >
            >::type* dummy = 0)
        {
            return 1;
        }


        template<typename R>
        static unsigned int GetRequiredColumnsFromMatrix(const R& matrix,
            typename boost::disable_if
            <
                boost::mpl::or_
                <
                    IsMatrix<R> ,
                    IsVector<R> 
                >
            >::type* dummy = 0)
        {
            return 1;
        }      
        
        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            unsigned int rows = GetRequiredRowsFromMatrix(boost::fusion::at_c<MappedIndex>(args));
            unsigned int columns = GetRequiredColumnsFromMatrix(boost::fusion::at_c<MappedIndex>(args));

            return boost::make_tuple(rows, columns, rows*columns);
        }

        
    };


    template<typename ChildType,
             typename Op,
             typename Indices, unsigned int index>
    struct MatrixSize<expt::Node<ChildType, Op, void>, Indices, index>
    {
        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            return MatrixSize<ChildType, Indices, index>::GetRequiredSize(args);
        }
        
        template<typename ArgumentVectorType>
        static unsigned int
        GetRequiredRows(const ArgumentVectorType& args)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> values = GetRequiredSize(args);
            return values.get<0>();
        }
    };

    
    /// \brief Calculates the size of the given matrix expression, as well as the 
    ///        largest intermediate buffer size, assuming multiplication.
    template<typename LeftNodeType, typename RightNodeType, typename Indices, unsigned int index>
    struct CalculateLargestRequiredSize
    {
        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
            MatrixSize<LeftNodeType, Indices, index>::GetRequiredSize(args);

            boost::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
                MatrixSize<RightNodeType, Indices, index + LeftNodeType::TotalCount>::GetRequiredSize(args);

            unsigned int leftRows = lhsSizes.get<0>();
            unsigned int rightColumns = rhsSizes.get<1>();

            unsigned int matrixSize = leftRows*rightColumns;
            unsigned int bufferSize = std::max(std::max(lhsSizes.get<2>(), rhsSizes.get<2>()), matrixSize);
            
            return boost::make_tuple(leftRows, rightColumns, bufferSize);
        }
    };


    template<typename LhsType, typename OpType, typename RhsType, typename Indices, unsigned int index, typename enabled=void>
    struct BinaryMatrixSizeEvaluator;

    // Addition or subtraction.
    template<typename LhsType, typename OpType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, OpType, RhsType, Indices, index,
        typename boost::enable_if
        <
            boost::mpl::or_
            <
                boost::is_same<OpType, expt::AddOp>,
                boost::is_same<OpType, expt::SubtractOp>
            >
        >::type>
    {
        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
                MatrixSize<LhsType, Indices, index>::GetRequiredSize(args);
            return lhsSizes;
        }
    };

    // Multiplication with double lhs.
    template<typename LhsType, typename OpType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, OpType, RhsType, Indices, index,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                expt::IsConstantNode<LhsType>,
                boost::is_same<typename LhsType::ResultType, double>
            >
        >::type>
    {
        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
                MatrixSize<RhsType, Indices, index + LhsType::TotalCount>::GetRequiredSize(args);
            return rhsSizes;
        }
    };


    // Double rhs
    template<typename LhsType, typename OpType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, OpType, RhsType, Indices, index,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                expt::IsConstantNode<RhsType>,
                boost::is_same<typename RhsType::ResultType, double>
            >
        >::type>
    {
        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
                MatrixSize<LhsType, Indices, index>::GetRequiredSize(args);
            return lhsSizes;
        }
    };

    // Multiplication.
    template<typename LhsType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, expt::MultiplyOp, RhsType, Indices, index,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                // Test for both nodes being constant matrices.
                // Negate to account for the next specialization.
                boost::mpl::not_<boost::mpl::and_
                <
                    boost::mpl::and_
                    <
                        expt::IsConstantNode<RhsType>,
                        boost::mpl::not_<boost::is_same<typename RhsType::ResultType, double> >
                    >,
                    boost::mpl::and_
                    <
                        expt::IsConstantNode<LhsType>,
                        boost::mpl::not_<boost::is_same<typename LhsType::ResultType, double> >
                    >
                > >,

                // Tests to make sure if either is constant, it is a matrix.
                boost::mpl::not_<boost::mpl::and_
                <
                    expt::IsConstantNode<RhsType>,
                    boost::is_same<typename RhsType::ResultType, double>
                > >,
                boost::mpl::not_<boost::mpl::and_
                <
                    expt::IsConstantNode<LhsType>,
                    boost::is_same<typename LhsType::ResultType, double>
                > >
            >
        >::type>
    {
        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            return CalculateLargestRequiredSize<LhsType, RhsType, Indices, index>::GetRequiredSize(args);
        }
    };

    template<typename LhsType, typename RhsType,typename Indices, unsigned int index>
    struct BinaryMatrixSizeEvaluator<LhsType, expt::MultiplyOp, RhsType, Indices, index,
        typename boost::enable_if
        <
            // Test for both nodes being constant matrices.
            boost::mpl::and_
            <
                boost::mpl::and_
                <
                    expt::IsConstantNode<RhsType>,
                    boost::mpl::not_<boost::is_same<typename RhsType::ResultType, double> >
                >,
                boost::mpl::and_
                <
                    expt::IsConstantNode<LhsType>,
                    boost::mpl::not_<boost::is_same<typename LhsType::ResultType, double> >
                >
            >
        >::type>
    {
        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
                MatrixSize<LhsType, Indices, index>::GetRequiredSize(args);

            boost::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
                MatrixSize<RhsType, Indices, index + LhsType::TotalCount>::GetRequiredSize(args);

            unsigned int leftRows = lhsSizes.get<0>();
            unsigned int rightColumns = rhsSizes.get<1>();

            unsigned int bufferSize = leftRows*rightColumns;
            
            return boost::make_tuple(leftRows, rightColumns, bufferSize);
        }
    };

    ///// \brief Implementation of the MatrixSize class for binary nodes.
    /////        The code is slightly cleaner by moving the evaluation to this class
    /////        instead of creating the specializations in the MatrixSize class.
    //template<typename L1, typename LOp, typename L2,
    //         typename Op,
    //         typename R1, typename ROp, typename R2,
    //        typename Indices, unsigned int index,
    //        typename enabled = void>
    //struct BinaryMatrixSizeEvaluator;

    ///// \brief Matrix sizes for two constant children with addition/subtraction.
    //template<typename L, typename Op, typename R, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<L, void, void, Op, R, void, void, Indices, index>
    //{
    //    typedef expt::Node<L, void, void> LeftNodeType;
    //    typedef expt::Node<R, void, void> RightNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
    //            MatrixSize<LeftNodeType, Indices, index>::GetRequiredSize(args);
    //        return lhsSizes;
    //    }
    //};
 
    ///// \brief Matrix sizes for two constant children with multiplication.
    //template<typename L, typename R, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<L, void, void, expt::MultiplyOp, R, void, void, Indices, index>
    //{
    //    typedef expt::Node<L, void, void> LeftNodeType;
    //    typedef expt::Node<R, void, void> RightNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
    //            MatrixSize<LeftNodeType, Indices, index>::GetRequiredSize(args);

    //        boost::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
    //            MatrixSize<RightNodeType, Indices, index + LeftNodeType::TotalCount>::GetRequiredSize(args);

    //        unsigned int leftRows = lhsSizes.get<0>();
    //        unsigned int rightColumns = rhsSizes.get<1>();

    //        unsigned int bufferSize = leftRows*rightColumns;
    //        
    //        return boost::make_tuple(leftRows, rightColumns, bufferSize);
    //    }
    //};

    //template<typename R, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<NekDouble, void, void, expt::MultiplyOp, R, void, void, Indices, index>
    //{
    //    typedef expt::Node<NekDouble, void, void> LeftNodeType;
    //    typedef expt::Node<R, void, void> RightNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        boost::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
    //            MatrixSize<RightNodeType, Indices, index + LeftNodeType::TotalCount>::GetRequiredSize(args);
    //        return rhsSizes;
    //    }
    //};

    //template<typename L, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<L, void, void, expt::MultiplyOp, NekDouble, void, void, Indices, index>
    //{
    //    typedef expt::Node<L, void, void> LeftNodeType;
    //    typedef expt::Node<NekDouble, void, void> RightNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
    //            MatrixSize<LeftNodeType, Indices, index>::GetRequiredSize(args);
    //        return lhsSizes;
    //    }
    //};


    //template<typename L1, typename LOp, typename L2, typename Op, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<L1, LOp, L2, Op, NekDouble, void, void, Indices, index,
    //    typename boost::enable_if
    //    <
    //        boost::mpl::and_
    //        <
    //            boost::mpl::not_<boost::is_same<LOp, void> >,
    //            boost::mpl::not_<boost::is_same<L1, void> >
    //        >
    //    >::type >
    //{
    //    typedef expt::Node<L1, LOp, L2> LeftNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
    //            MatrixSize<LeftNodeType, Indices, index>::GetRequiredSize(args);
    //        return lhsSizes;
    //    }
    //};

    //template<typename L1, typename LOp, typename L2, typename Op, typename R, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<L1, LOp, L2, Op, R, void, void, Indices, index,
    //    typename boost::enable_if
    //    <
    //        boost::mpl::and_
    //        <
    //            boost::mpl::not_<boost::is_same<LOp, void> >,
    //            boost::mpl::not_<boost::is_same<L1, void> >
    //        >
    //    >::type >
    //{
    //    typedef expt::Node<R, void, void> RightNodeType;
    //    typedef expt::Node<L1, LOp, L2> LeftNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        return CalculateLargestRequiredSize<LeftNodeType, RightNodeType, Indices, index>::GetRequiredSize(args);
    //    }
    //};

    //template<typename Op, typename R1, typename ROp, typename R2, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<NekDouble, void, void, Op, R1, ROp, R2, Indices, index,
    //    typename boost::enable_if
    //    <
    //        boost::mpl::and_
    //        <
    //            boost::mpl::not_<boost::is_same<ROp, void> >,
    //            boost::mpl::not_<boost::is_same<R1, void> >
    //        >
    //    >::type >
    //{
    //    typedef expt::Node<R1, ROp, R2> RightNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        boost::tuple<unsigned int, unsigned int, unsigned int> rhsSizes = 
    //            MatrixSize<RightNodeType, Indices, index+1>::GetRequiredSize(args);
    //        return rhsSizes;
    //    }
    //};

    //template<typename L, typename Op, typename R1, typename ROp, typename R2, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<L, void, void, Op, R1, ROp, R2, Indices, index,
    //    typename boost::enable_if
    //    <
    //        boost::mpl::and_
    //        <
    //            boost::mpl::not_<boost::is_same<ROp, void> >,
    //            boost::mpl::not_<boost::is_same<R1, void> >
    //        >
    //    >::type >
    //{
    //    typedef expt::Node<L, void, void> LeftNodeType;
    //    typedef expt::Node<R1, ROp, R2> RightNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        return CalculateLargestRequiredSize<LeftNodeType, RightNodeType, Indices, index>::GetRequiredSize(args);
    //    }
    //};

    //template<typename L1, typename LOp, typename L2, 
    //         typename Op,
    //         typename R1, typename ROp, typename R2, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<L1, LOp, L2, Op, R1, ROp, R2, Indices, index,
    //    typename boost::enable_if
    //    <
    //        boost::mpl::and_
    //        <
    //            boost::mpl::not_<boost::is_same<ROp, void> >,
    //            boost::mpl::not_<boost::is_same<R2, void> >,
    //            boost::mpl::not_<boost::is_same<LOp, void> >,
    //            boost::mpl::not_<boost::is_same<L2, void> >
    //        >
    //    >::type >
    //{
    //    typedef expt::Node<L1, LOp, L2> LeftNodeType;
    //    typedef expt::Node<R1, ROp, R2> RightNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        boost::tuple<unsigned int, unsigned int, unsigned int> lhsSizes = 
    //            MatrixSize<LeftNodeType, Indices, index>::GetRequiredSize(args);
    //        return lhsSizes;
    //    }
    //};

    //template<typename L1, typename LOp, typename L2, 
    //         typename R1, typename ROp, typename R2, typename Indices, unsigned int index>
    //struct BinaryMatrixSizeEvaluator<L1, LOp, L2, expt::MultiplyOp, R1, ROp, R2, Indices, index,
    //    typename boost::enable_if
    //    <
    //        boost::mpl::and_
    //        <
    //            boost::mpl::not_<boost::is_same<ROp, void> >,
    //            boost::mpl::not_<boost::is_same<R2, void> >,
    //            boost::mpl::not_<boost::is_same<LOp, void> >,
    //            boost::mpl::not_<boost::is_same<L2, void> >
    //        >
    //    >::type >
    //{
    //    typedef expt::Node<L1, LOp, L2> LeftNodeType;
    //    typedef expt::Node<R1, ROp, R2> RightNodeType;
    //    
    //    template<typename ArgumentVectorType>
    //    static boost::tuple<unsigned int, unsigned int, unsigned int>
    //    GetRequiredSize(const ArgumentVectorType& args)
    //    {
    //        return CalculateLargestRequiredSize<LeftNodeType, RightNodeType, Indices, index>::GetRequiredSize(args);
    //    }
    //};

    template<typename L1, typename LOp, typename L2,
             typename Op,
             typename R1, typename ROp, typename R2,
             typename Indices, unsigned int index>
    struct MatrixSize<expt::Node<expt::Node<L1, LOp, L2>, Op, expt::Node<R1, ROp, R2> >, Indices, index>
    {
        typedef expt::Node<L1, LOp, L2> LhsType;
        typedef expt::Node<R1, ROp, R2> RhsType;

        template<typename ArgumentVectorType>
        static boost::tuple<unsigned int, unsigned int, unsigned int>
        GetRequiredSize(const ArgumentVectorType& args)
        {
            return BinaryMatrixSizeEvaluator<LhsType, Op, RhsType, Indices, index>::GetRequiredSize(args);
        }
        
        template<typename ArgumentVectorType>
        static unsigned int
        GetRequiredRows(const ArgumentVectorType& args)
        {
            boost::tuple<unsigned int, unsigned int, unsigned int> values = GetRequiredSize(args);
            return values.get<0>();
        }
    };

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIBUTILITIES_LINEAR_ALGEBRA_MATRIX_SIZE_HPP

