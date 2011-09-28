///////////////////////////////////////////////////////////////////////////////
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


#ifndef EXPRESSION_TEMPLATES_NODE_HPP
#define EXPRESSION_TEMPLATES_NODE_HPP

#define FUSION_MAX_VECTOR_SIZE 50
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/include/size.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>
#include <boost/type_traits.hpp>
#include <boost/call_traits.hpp>
#include <iostream>
#include <ExpressionTemplates/Operators.hpp>

namespace expt
{
    // Creates a vector of constants from [0..n]
    template<typename T, template<int> class TypeConstructor, int n>
    struct CreateVectorC
    {
        typedef typename CreateVectorC<T, TypeConstructor, n-1>::type ListType;
        typedef typename boost::mpl::push_back<ListType, TypeConstructor<n> >::type type;
    };

    template<typename T, template<int> class TypeConstructor>
    struct CreateVectorC<T, TypeConstructor, 0>
    {
        typedef boost::mpl::vector_c<T, 0> type;
    };

    template<typename T1, typename Op = void, typename T2 = void>
    struct Node;


    // Binary node.
    template<typename L1, typename LOp, typename L2, typename Op, typename R1, typename ROp, typename R2>
    struct Node<Node<L1, LOp, L2>, Op, Node<R1, ROp, R2> >
    {
        public:
            typedef Node<L1, LOp, L2> Left;
            typedef Node<R1, ROp, R2> Right;
            typedef Op OpType;
            typedef Node<Left, Op, Right> ThisType;

            typedef typename Left::MPLVectorType LeftMPLVectorType;
            typedef typename boost::mpl::size<LeftMPLVectorType>::type LhsOperandCount;
            
            typedef typename Right::MPLVectorType RightMPLVectorType;
            typedef typename boost::mpl::size<RightMPLVectorType>::type RhsOperandCount;
            
            static const int TotalCount = LhsOperandCount::value + RhsOperandCount::value;
            
            typedef typename CreateVectorC<int, boost::mpl::int_, LhsOperandCount::value + RhsOperandCount::value - 1>::type Indices;
            
            typedef typename Left::VectorType LeftVectorType;        
            typedef typename Right::VectorType RightVectorType;
            
            typedef typename boost::mpl::insert_range
                             <
                                 LeftMPLVectorType,
                                 typename boost::mpl::end<LeftMPLVectorType>::type,
                                 RightMPLVectorType
                             >::type MPLVectorType;
            
            typedef typename boost::fusion::result_of::as_vector<MPLVectorType>::type VectorType;
            
            typedef typename OpType::template ResultType<typename Left::ResultType, typename Right::ResultType>::type ResultType;

        public:
                         
            Node(const Left& lhs, const Right& rhs) :
                m_data(CreateMData(lhs, rhs))
            {

            }
            
            Node(const ThisType& rhs) :
                m_data(rhs.m_data)
            {
            }

            ThisType& operator=(const ThisType& rhs)
            {
                m_data = rhs.m_data;
                return *this;
            }

            const VectorType& GetData() const { return m_data; }
            
        private:
        
            // Workaround because I can't seem to do this in the initializer list.
            VectorType CreateMData(const Left& lhs, const Right& rhs)
            {
                LeftVectorType l = lhs.GetData();
                RightVectorType r = rhs.GetData();
                boost::fusion::joint_view<LeftVectorType, RightVectorType> joint(l, r);
                return VectorType(joint);
            }
                        
            VectorType m_data;        
    };

    
    
    // Constant Node
    // T1 is a normal c++ data type and should not be a node object.
    template<typename T>
    struct Node<T, void, void>
    {
            typedef Node<T, void, void> ThisType;

            // ResultType should give a C++ type with all const/volatile removed.
            // We are assuming that we'll never be dealing with pointers.
            // Result type is only used in type deduction.
            typedef typename boost::remove_cv
                             <
                                typename boost::remove_reference<T>::type
                             >::type ResultType;
                             
            
            // We need to store a const& to the data.
            typedef typename boost::call_traits<T>::const_reference DataType;
            
            static const int TotalCount = 1;
            
            // I need MPL so I can get the types then convert to a 
            // fusion vector (need random access).  Using the fusion::typeof
            // with insert_range gets me a forward sequence, not random access.
            typedef boost::mpl::vector<DataType> MPLVectorType;
            typedef typename boost::fusion::result_of::as_vector<MPLVectorType>::type VectorType;
            typedef typename CreateVectorC<int, boost::mpl::int_, 1>::type Indices;
            typedef void OpType;
            
        public:
            explicit Node(DataType value) :
                m_data(value)
            {
            }
            
            Node(const ThisType& rhs) :
                m_data(rhs.m_data)
            {
            }
            
            ThisType& operator=(const ThisType& rhs)
            {
                m_data = rhs.m_data;
                return *this;
            }

            const VectorType& GetData() const { return m_data; }
            
        private:
            VectorType m_data;
    };

    
    template<typename NodeType, typename enabled = void>
    struct TemporaryCount;
    
    template<typename T>
    struct TemporaryCount<Node<T, void, void> >
    {
        static const unsigned int Value = 0;
    };

    template<typename T, typename OpType>
    struct TemporaryCount< Node<T, OpType, void> >
    {
        static const unsigned int Value = TemporaryCount<T>::Value;
    };

    template<typename L, typename OpType, typename R1, typename ROp, typename R2>
    struct TemporaryCount< Node<L, OpType, Node<R1, ROp, R2> >,
        typename boost::enable_if
        <
            boost::mpl::or_
            <
                boost::is_same<ROp, void>,
                boost::is_same<R2, void>
            >
        >::type>
    {
        static const unsigned int Value = TemporaryCount<L>::Value + TemporaryCount<Node<R1, ROp, R2> >::Value;
    };

    template<typename L, typename OpType, typename R1, typename ROp, typename R2>
    struct TemporaryCount< Node<L, OpType, Node<R1, ROp, R2> >, 
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                boost::mpl::not_<boost::is_same<ROp, void> >,
                boost::mpl::not_<boost::is_same<R2, void> >
            >
        >::type>
    {
        static const unsigned int Value = TemporaryCount<L>::Value + TemporaryCount<Node<R1, ROp, R2> >::Value + 1;
    };

    template<typename L, typename Op, typename R>
    std::ostream& operator<<(std::ostream& os, const Node<L, Op, R>& exp)
    {
        return os;
    }

    // Unary Node
    template<typename T1, typename OpType>
    struct Node<Node<T1>, OpType, void>
    {
        static const int i = 2;
        typedef T1 ResultType;
        typedef typename Node<T1>::Indices Indices;
    };


    template<typename L1, typename LOp, typename L2, typename RhsData>
    Node<Node<L1, LOp, L2>, expt::SubtractOp, Node<RhsData> >
    operator-(const Node<L1, LOp, L2>& lhsNode, const RhsData& rhs)
    {
        Node<RhsData> rhsNode(rhs);
        return Node<Node<L1, LOp, L2>, expt::SubtractOp, Node<RhsData> >(lhsNode, rhsNode);
    }

    template<typename L1, typename LOp, typename L2, typename RhsData>
    Node<Node<L1, LOp, L2>, expt::AddOp, Node<RhsData> >
    operator+(const Node<L1, LOp, L2>& lhsNode, const RhsData& rhs)
    {
        Node<RhsData> rhsNode(rhs);
        return Node<Node<L1, LOp, L2>, expt::AddOp, Node<RhsData> >(lhsNode, rhsNode);
    }

    template<typename L1, typename LOp, typename L2, typename RhsData>
    Node<Node<L1, LOp, L2>, expt::MultiplyOp, Node<RhsData> >
    operator*(const Node<L1, LOp, L2>& lhsNode, const RhsData& rhs)
    {
        Node<RhsData> rhsNode(rhs);
        return Node<Node<L1, LOp, L2>, expt::MultiplyOp, Node<RhsData> >(lhsNode, rhsNode);
    }

    template<typename L1, typename LOp, typename L2, typename RhsData>
    Node<Node<L1, LOp, L2>, expt::DivideOp, Node<RhsData> >
    operator/(const Node<L1, LOp, L2>& lhsNode, const RhsData& rhs)
    {
        Node<RhsData> rhsNode(rhs);
        return Node<Node<L1, LOp, L2>, expt::DivideOp, Node<RhsData> >(lhsNode, rhsNode);
    }


    template<typename LhsData, typename R1, typename ROp, typename R2>
    Node<Node<LhsData>, expt::SubtractOp, Node<R1, ROp, R2> >
    operator-(const LhsData& lhsData, const Node<R1, ROp, R2>& rhsNode)
    {
        Node<LhsData> lhsNode(lhsData);
        return Node<Node<LhsData>, expt::SubtractOp, Node<R1, ROp, R2> >(lhsNode, rhsNode);
    }

    template<typename LhsData, typename R1, typename ROp, typename R2>
    Node<Node<LhsData>, expt::AddOp, Node<R1, ROp, R2> >
    operator+(const LhsData& lhsData, const Node<R1, ROp, R2>& rhsNode)
    {
        Node<LhsData> lhsNode(lhsData);
        return Node<Node<LhsData>, expt::AddOp, Node<R1, ROp, R2> >(lhsNode, rhsNode);
    }

    template<typename LhsData, typename R1, typename ROp, typename R2>
    Node<Node<LhsData>, expt::MultiplyOp, Node<R1, ROp, R2> >
    operator*(const LhsData& lhsData, const Node<R1, ROp, R2>& rhsNode)
    {
        Node<LhsData> lhsNode(lhsData);
        return Node<Node<LhsData>, expt::MultiplyOp, Node<R1, ROp, R2> >(lhsNode, rhsNode);
    }

    template<typename LhsData, typename R1, typename ROp, typename R2>
    Node<Node<LhsData>, expt::DivideOp, Node<R1, ROp, R2> >
    operator/(const LhsData& lhsData, const Node<R1, ROp, R2>& rhsNode)
    {
        Node<LhsData> lhsNode(lhsData);
        return Node<Node<LhsData>, expt::DivideOp, Node<R1, ROp, R2> >(lhsNode, rhsNode);
    }

    template<typename L1, typename LOp, typename L2, typename R1, typename ROp, typename R2>
    Node<Node<L1, LOp, L2>, expt::AddOp, Node<R1, ROp, R2> >
    operator+(const Node<L1, LOp, L2>& lhsNode, const Node<R1, ROp, R2>& rhsNode)
    {
        return Node<Node<L1, LOp, L2>, expt::AddOp, Node<R1, ROp, R2> >(lhsNode, rhsNode);
    }

    template<typename L1, typename LOp, typename L2, typename R1, typename ROp, typename R2>
    Node<Node<L1, LOp, L2>, expt::SubtractOp, Node<R1, ROp, R2> >
    operator-(const Node<L1, LOp, L2>& lhsNode, const Node<R1, ROp, R2>& rhsNode)
    {
        return Node<Node<L1, LOp, L2>, expt::SubtractOp, Node<R1, ROp, R2> >(lhsNode, rhsNode);
    }

    template<typename L1, typename LOp, typename L2, typename R1, typename ROp, typename R2>
    Node<Node<L1, LOp, L2>, expt::MultiplyOp, Node<R1, ROp, R2> >
    operator*(const Node<L1, LOp, L2>& lhsNode, const Node<R1, ROp, R2>& rhsNode)
    {
        return Node<Node<L1, LOp, L2>, expt::MultiplyOp, Node<R1, ROp, R2> >(lhsNode, rhsNode);
    }

    template<typename L1, typename LOp, typename L2, typename R1, typename ROp, typename R2>
    Node<Node<L1, LOp, L2>, expt::DivideOp, Node<R1, ROp, R2> >
    operator/(const Node<L1, LOp, L2>& lhsNode, const Node<R1, ROp, R2>& rhsNode)
    {
        return Node<Node<L1, LOp, L2>, expt::DivideOp, Node<R1, ROp, R2> >(lhsNode, rhsNode);
    }

    template<typename Op, typename L, typename R>
    Node<Node<L>, Op, Node<R> > CreateBinaryExpression(const L& lhs, const R& rhs)
    {       
        Node<L> lhs_node(lhs);
        Node<R> rhs_node(rhs);
        return Node<Node<L>, Op, Node<R> >(lhs_node, rhs_node);
    }
}

#endif //EXPRESSION_TEMPLATES_NODE_HPP
