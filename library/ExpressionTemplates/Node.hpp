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

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

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

    //////////////////////////////////////////////////////////////////
    // Constant Node
    /////////////////////////////////////////////////////////////////

    /// \brief A traits class to deterimine how the data should be stored.
    template<typename T, typename enabled = void>
	struct ConstNodeDataAccessor
	{
		// ResultType should give a C++ type with all const/volatile removed.
        // We are assuming that we'll never be dealing with pointers.
        // Result type is only used in type deduction.
        typedef typename boost::remove_cv<typename boost::remove_reference<T>::type>::type ResultType;
        
        // We need to store a const& to the data.
        typedef typename boost::call_traits<T>::const_reference DataType;
	};
    
	template<typename T>
	struct ConstNodeDataAccessor<T, typename boost::enable_if<boost::is_arithmetic<T> >::type >
	{
		// ResultType should give a C++ type with all const/volatile removed.
        // We are assuming that we'll never be dealing with pointers.
        // Result type is only used in type deduction.
        typedef typename boost::remove_cv<typename boost::remove_reference<T>::type>::type ResultType;
                             
        // For numeric data, we can store the data by value.
        typedef typename boost::call_traits<T>::value_type DataType;
	};

    // Constant Node
    // T1 is a normal c++ data type and should not be a node object.
    template<typename T>
    struct Node<T, void, void>
    {
        public:
            typedef Node<T, void, void> ThisType;

            typedef typename ConstNodeDataAccessor<T>::ResultType ResultType;
			typedef typename ConstNodeDataAccessor<T>::DataType DataType;

            // Total number of nodes represented by the tree rooted at this node.
            static const int TotalCount = 1;
            
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

    // Unary Node
    template<typename T1, typename TOp, typename T2, typename Op>
    struct Node<Node<T1, TOp, T2>, Op, void>
    {
        public:
            typedef Op OpType;
            typedef Node<T1, TOp, T2> ChildNodeType;
            typedef Node<Node<T1, TOp, T2>, Op, void> ThisType;

            // Total number of leaf nodes represented by this node.
            static const int TotalCount = ChildNodeType::TotalCount;

            typedef typename OpType::template ResultType<typename ChildNodeType::ResultType>::type ResultType;
            typedef typename ChildNodeType::Indices Indices;
            typedef typename ChildNodeType::VectorType VectorType;
            typedef typename ChildNodeType::MPLVectorType MPLVectorType;

            explicit Node(const ChildNodeType& value) :
                    m_data(value.GetData())
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
            
            // Total number of leaf nodes represented by this node.
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

    
    // Utility metafunctions to help categorize nodes.

    template<typename NodeType>
    struct IsConstantNode : public boost::false_type {};

    template<typename DataType>
    struct IsConstantNode<Node<DataType, void, void> > : public boost::true_type {};

    template<typename NodeType>
    struct IsUnaryNode : public boost::false_type {};

    template<typename ChildNodeType, typename UnaryType>
    struct IsUnaryNode<Node<ChildNodeType, UnaryType, void> > : public boost::true_type {};

    template<typename NodeType>
    struct IsBinaryNode : public boost::true_type {};

    template<typename L1, typename OpType>
    struct IsBinaryNode<Node<L1, OpType, void> > : public boost::false_type {};

    template<typename L1>
    struct IsBinaryNode<Node<L1, void, void> > : public boost::false_type {};
   
    
    /// \brief Counts the number of temporaries required a tree.
    /// Value returns the number of temporaries required.
    /// AsRhs gives the number assuming that the node is on the rhs of a 
    /// binary node.
    template<typename NodeType, typename enabled = void>
    struct TemporaryCount;
    
    template<typename T>
    struct TemporaryCount<Node<T, void, void> >
    {
        static const unsigned int Value = 0;
        static const unsigned int AsRhs = 0;
    };

    template<typename T, typename OpType>
    struct TemporaryCount< Node<T, OpType, void> >
    {
        static const unsigned int Value = TemporaryCount<T>::Value;
        static const unsigned int AsRhs = Value + 1;
    };

    template<typename L, typename OpType, typename R>
    struct TemporaryCount<Node<L, OpType, R>, 
        typename boost::enable_if
        <
            IsConstantNode<R>
        >::type>
    {
        static const unsigned int Value = TemporaryCount<L>::Value;
        static const unsigned int AsRhs = Value + 1;
    };
    
    template<typename L, typename OpType, typename R>
    struct TemporaryCount<Node<L, OpType, R>, 
        typename boost::disable_if
        <
            IsConstantNode<R>
        >::type>
    {
        static const unsigned int Value = TemporaryCount<L>::Value + TemporaryCount<R>::Value + 1;
        static const unsigned int AsRhs = Value + 1;
    };


    template<typename L, typename Op, typename R>
    std::ostream& operator<<(std::ostream& os, const Node<L, Op, R>& exp)
    {
        return os;
    }

    // Operators to combine nodes.
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

    template<typename T1, typename TOp, typename T2>
    Node<Node<T1, TOp, T2>, NegateOp, void>
    operator-(const Node<T1, TOp, T2>& node)
    {
        return Node<Node<T1, TOp, T2>, NegateOp, void>(node);
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

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#endif //EXPRESSION_TEMPLATES_NODE_HPP
