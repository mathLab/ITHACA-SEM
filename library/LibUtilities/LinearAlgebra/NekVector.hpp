///////////////////////////////////////////////////////////////////////////////
//
// File: NekVector.hpp
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
// Description: Generic N-Dimensional Vector.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP
#define NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/NekPoint.hpp>

#include <LibUtilities/LinearAlgebra/NekVectorMetadata.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorConstantSized.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorVariableSized.hpp>

#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>

#include <functional>
#include <algorithm>
#include <math.h>

#include <boost/call_traits.hpp>
#include <boost/type_traits.hpp>
#include <boost/shared_array.hpp>


namespace Nektar
{
           
    template<typename DataType, typename dim, typename space>
    void Add(NekVector<DataType, dim, space>& result,
           const NekVector<const DataType, dim, space>& lhs,
           const NekVector<const DataType, dim, space>& rhs)
    {
        DataType* r_buf = result.GetRawPtr();
        const DataType* lhs_buf = lhs.GetRawPtr();
        const DataType* rhs_buf = rhs.GetRawPtr();
        for(int i = 0; i < lhs.GetDimension(); ++i)
        {
            r_buf[i] = lhs_buf[i] + rhs_buf[i];
        }
    }
    
    template<typename DataType, typename dim, typename space>
    void AddEqual(NekVector<DataType, dim, space>& result,
           const NekVector<const DataType, dim, space>& rhs)
    {
        DataType* r_buf = result.GetRawPtr();
        const DataType* rhs_buf = rhs.GetRawPtr();
        for(int i = 0; i < rhs.GetDimension(); ++i)
        {
            r_buf[i] += rhs_buf[i];
        }
    }
    
    template<typename DataType, typename dim, typename space>
    NekVector<DataType, dim, space> Add(const NekVector<DataType, dim, space>& lhs, 
                                           const NekVector<DataType, dim, space>& rhs)
    {
        NekVector<DataType, dim, space> result(lhs.GetDimension());
        Add(result, lhs, rhs);
        return result;
    }
    
    template<typename ResultDataType, typename InputDataType, typename dim, typename space>
    void Subtract(NekVector<ResultDataType, dim, space>& result,
           const NekVector<InputDataType, dim, space>& lhs,
           const NekVector<InputDataType, dim, space>& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        typename boost::add_const<InputDataType>::type* lhs_buf = lhs.GetRawPtr();
        typename boost::add_const<InputDataType>::type* rhs_buf = rhs.GetRawPtr();
        for(int i = 0; i < lhs.GetDimension(); ++i)
        {
            r_buf[i] = lhs_buf[i] - rhs_buf[i];
        }
    }
    
    template<typename ResultDataType, typename InputDataType, typename dim, typename space>
    void SubtractEqual(NekVector<ResultDataType, dim, space>& result,
           const NekVector<InputDataType, dim, space>& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        typename boost::add_const<InputDataType>::type* rhs_buf = rhs.GetRawPtr();
        for(int i = 0; i < rhs.GetDimension(); ++i)
        {
            r_buf[i] -= rhs_buf[i];
        }
    }
    
    template<typename DataType, typename dim, typename space>
    NekVector<typename boost::remove_const<DataType>::type, dim, space>
    Subtract(const NekVector<DataType, dim, space>& lhs,
                const NekVector<DataType, dim, space>& rhs)
    {
        NekVector<typename boost::remove_const<DataType>::type, dim, space> result(lhs.GetDimension());
        Subtract(result, lhs, rhs);
        return result;
    }





	template<typename ResultDataType, typename InputDataType, typename dim, typename space>
    void Divide(NekVector<ResultDataType, dim, space>& result,
           const NekVector<InputDataType, dim, space>& lhs,
           const NekDouble& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        typename boost::add_const<InputDataType>::type* lhs_buf = lhs.GetRawPtr();
        
        for(int i = 0; i < lhs.GetDimension(); ++i)
        {
            r_buf[i] = lhs_buf[i] / rhs;
        }
    }
    
    template<typename ResultDataType, typename dim, typename space>
    void DivideEqual(NekVector<ResultDataType, dim, space>& result,
           const NekDouble& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        for(int i = 0; i < result.GetDimension(); ++i)
        {
            r_buf[i] /= rhs;
        }
    }
    
    template<typename DataType, typename dim, typename space>
    NekVector<typename boost::remove_const<DataType>::type, dim, space>
    Divide(const NekVector<DataType, dim, space>& lhs,
                const NekDouble& rhs)
    {
        NekVector<typename boost::remove_const<DataType>::type, dim, space> result(lhs.GetDimension());
        Divide(result, lhs, rhs);
        return result;
    }


	template<typename ResultDataType, typename InputDataType, typename dim, typename space>
    void Multiply(NekVector<ResultDataType, dim, space>& result,
           const NekVector<InputDataType, dim, space>& lhs,
           const NekDouble& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        typename boost::add_const<InputDataType>::type* lhs_buf = lhs.GetRawPtr();
        
        for(int i = 0; i < lhs.GetDimension(); ++i)
        {
            r_buf[i] = lhs_buf[i] * rhs;
        }
    }
    
    template<typename ResultDataType, typename dim, typename space>
    void MultiplyEqual(NekVector<ResultDataType, dim, space>& result,
           const NekDouble& rhs)
    {
        ResultDataType* r_buf = result.GetRawPtr();
        for(int i = 0; i < result.GetDimension(); ++i)
        {
            r_buf[i] *= rhs;
        }
    }
    
    template<typename DataType, typename dim, typename space>
    NekVector<typename boost::remove_const<DataType>::type, dim, space>
    Multiply(const NekVector<DataType, dim, space>& lhs,
                const NekDouble& rhs)
    {
        NekVector<typename boost::remove_const<DataType>::type, dim, space> result(lhs.GetDimension());
        Multiply(result, lhs, rhs);
        return result;
    }

	template<typename ResultDataType, typename InputDataType, typename dim, typename space>
    void Multiply(NekVector<ResultDataType, dim, space>& result,
			const NekDouble& lhs,   
			const NekVector<InputDataType, dim, space>& rhs)
    {
		Multiply(result, rhs, lhs);
    }
        
    template<typename DataType, typename dim, typename space>
    NekVector<typename boost::remove_const<DataType>::type, dim, space>
	Multiply(const NekDouble& lhs,
			 const NekVector<DataType, dim, space>& rhs)
    {
		return Multiply(rhs, lhs);
    }

    GENERATE_MULTIPLICATION_OPERATOR(NekVector, 3, NekDouble, 0);
    GENERATE_MULTIPLICATION_OPERATOR(NekDouble, 0, NekVector, 3);
    
    GENERATE_DIVISION_OPERATOR(NekVector, 3, NekDouble, 0);
    GENERATE_ADDITION_OPERATOR(NekVector, 3, NekVector, 3);
    GENERATE_SUBTRACTION_OPERATOR(NekVector, 3, NekVector, 3);


    template<typename DataType, typename dim, typename space>
    std::ostream& operator<<(std::ostream& os, const NekVector<DataType, dim, space>& rhs)
    {
        os << rhs.AsString();
        return os;
    }

    template<typename DataType, typename dim, typename space>
    NekVector<DataType, dim, space> createVectorFromPoints(const NekPoint<DataType, dim, space>& source,
                                                           const NekPoint<DataType, dim, space>& dest)
    {
        NekVector<DataType, dim, space> result;
        for(unsigned int i = 0; i < dim::Value; ++i)
        {
            result[i] = dest[i]-source[i];
        }
        return result;
    }

    template<typename DataType, typename dim, typename space>
    NekPoint<DataType, dim, space> findPointAlongVector(const NekVector<DataType, dim, space>& lhs,
                                                        typename boost::call_traits<DataType>::const_reference t)
    {
        NekPoint<DataType, dim, space> result;
        for(unsigned int i = 0; i < dim::Value; ++i)
        {
            result[i] = lhs[i]*t;
        }

        return result;
    }

    template<typename DataType, typename lhsDim, typename rhsDim, typename space>
    bool operator==(const NekVector<const DataType, lhsDim, space>& lhs,
                    const NekVector<const DataType, rhsDim, space>& rhs)
    {
        if( lhs.GetDimension() != rhs.GetDimension() )
        {
            return false;
        }
        
        return std::equal(lhs.begin(), lhs.end(), rhs.begin());
    }

    template<typename DataType, typename lhsDim, typename rhsDim, typename space>
    bool operator!=(const NekVector<DataType, lhsDim, space>& lhs,
                    const NekVector<DataType, rhsDim, space>& rhs)
    {
        return !(lhs == rhs);
    }

    template<typename T>
    struct RemoveVectorConst;
    
    template<typename DataType, typename Dim, typename Space>
    struct RemoveVectorConst<NekVector<DataType, Dim, Space> >
    {
        typedef NekVector<DataType, Dim, Space> type;
    };
    
    template<typename DataType, typename Dim, typename Space>
    struct RemoveVectorConst<NekVector<const DataType, Dim, Space> >
    {
        typedef NekVector<DataType, Dim, Space> type;
    };
        

    #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        // Override default expression handling for addition/subtraction of vectors.
        // Optimal execution is obtained by loop unrolling.

        template<typename NodeType, typename enabled = void>
        struct NodeCanUnroll : public boost::false_type
        {
        };

        template<typename Type>
        struct NodeCanUnroll<Node<Type, void, void>,
            typename boost::enable_if
            <
                IsVector<typename Node<Type, void, void>::ResultType>
            >::type > : public boost::true_type
        {
        };
        
        template<typename LhsType, typename OpType, typename RhsType>
        struct NodeCanUnroll<Node<LhsType, OpType, RhsType>, 
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    IsVector<typename LhsType::ResultType>,
                    IsVector<typename RhsType::ResultType>,
                    NodeCanUnroll<LhsType>,
                    NodeCanUnroll<RhsType>,
                    boost::mpl::or_
                    <
                        boost::is_same<OpType, AddOp>,
                        boost::is_same<OpType, SubtractOp>
                    >
                >
            >::type >: public boost::true_type
        {
        };

        template<typename NodeType, typename IndicesType, unsigned int index>
        struct Accumulate;

        template<typename LhsType, typename IndicesType, unsigned int index>
        struct Accumulate<Node<LhsType, void, void>, IndicesType, index>
        {
            static const unsigned int MappedIndex = boost::mpl::at_c<IndicesType, index>::type::value;

            template<typename ResultType, typename ArgumentVectorType>
            static void Execute(ResultType& accumulator, const ArgumentVectorType& args, unsigned int i)
            {
                accumulator = boost::fusion::at_c<MappedIndex>(args)[i];
            }
        };

        template<typename LhsType, typename Op, typename RhsType, typename IndicesType, unsigned int index>
        struct Accumulate<Node<LhsType, Op, RhsType>, IndicesType, index>
        {
            static const int rhsNodeIndex = index + LhsType::TotalCount;

            template<typename ResultType, typename ArgumentVectorType>
            static void Execute(ResultType& accumulator, const ArgumentVectorType& args, unsigned int i)
            {
                Accumulate<LhsType, IndicesType, index>::Execute(accumulator, args, i);
                ResultType rhs;
                Accumulate<RhsType, IndicesType, rhsNodeIndex>::Execute(rhs, args, i);
                Op::OpEqual(accumulator, rhs);
            }
        };

        // Conditions
        // Lhs and Rhs must result in a vector.
        // Op must be Plus or Minus
        // Must apply recursively.
        template<typename LhsType, typename Op, typename RhsType, typename IndicesType, unsigned int index>
        struct BinaryBinaryEvaluateNodeOverride<LhsType, Op, RhsType, IndicesType, index,
            typename boost::enable_if
            <
                NodeCanUnroll<Node<LhsType, Op, RhsType> >
            >::type
         > : public boost::true_type 
        {
            template<typename ResultType, typename ArgumentVectorType>
            static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
            {
                static const int endIndex = index + LhsType::TotalCount + RhsType::TotalCount;
                typedef typename ResultType::DataType DataType;
                
                for(int i = 0; i < accumulator.GetRows(); ++i)
                {
                    DataType result = 0;
                    Accumulate<Node<LhsType, Op, RhsType>, IndicesType, index>::Execute(result, args, i);
                    accumulator[i] = result;
                }
            }
        };
    #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

}

#endif // NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

