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

#ifndef EXPRESSION_TEMPLATES_COMMUTATIVE_TRANSFORM_HPP
#define EXPRESSION_TEMPLATES_COMMUTATIVE_TRANSFORM_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/CommutativeTraits.hpp>

#include <boost/utility/enable_if.hpp>

namespace expt
{
    /// \brief Utility metafunction to swap indices.
    ///
    /// For a given sequence, this metafunction swaps the segment [start, partition) with 
    /// [partition, end).  
    ///
    /// For example, with the sequence 0, 1, 2, 3, 4, 5, start = 2, partition =4, end = 6, the 
    /// result of this metafunction is the sequence 0, 1, 4, 5, 2, 3.
    template<typename InputSequence, unsigned int start, unsigned int partition, unsigned int end>
    struct Swap
    {
        typedef typename boost::mpl::begin<InputSequence>::type Begin;
        typedef Begin PrefixStart;
        typedef typename boost::mpl::advance_c< Begin, start >::type PrefixEnd;
        
        typedef PrefixEnd FirstPartitionStart;
        typedef typename boost::mpl::advance_c< Begin, partition>::type FirstPartitionEnd;
        
        typedef FirstPartitionEnd SecondPartitionStart;
        typedef typename boost::mpl::advance_c<Begin, end>::type SecondPartitionEnd;
        
        typedef SecondPartitionEnd SuffixStart;
        typedef typename boost::mpl::end<InputSequence>::type SuffixEnd;
        
        typedef typename boost::mpl::iterator_range<PrefixStart, PrefixEnd>::type Prefix;
        typedef typename boost::mpl::iterator_range<FirstPartitionStart, FirstPartitionEnd>::type FirstPartition;
        typedef typename boost::mpl::iterator_range<SecondPartitionStart, SecondPartitionEnd>::type SecondPartition;
        typedef typename boost::mpl::iterator_range<SuffixStart, SuffixEnd>::type Suffix;
        
        typedef typename boost::mpl::joint_view<Prefix, SecondPartition>::type View1;
        typedef typename boost::mpl::joint_view<FirstPartition, Suffix>::type View2;
        typedef typename boost::mpl::joint_view<View1, View2>::type type;
        
    };

    /// \brief Performs a commutative transform on a node if the operator is commutative.
    template<typename NodeType, typename IndicesType, unsigned int IndexStart, typename enabled=void>
    struct CommutativeTransform
    {
        typedef NodeType TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };

    template<typename Left, typename Op, typename Right, typename IndicesType, unsigned int IndexStart>
    struct CommutativeTransform<Node<Left, Op, Right>, IndicesType, IndexStart,
        typename boost::enable_if<CommutativeTraits<typename Left::ResultType, Op, typename Right::ResultType> >::type>
    {
        typedef Node<Left, Op, Right> BaseNode;
        typedef Node<Right, Op, Left> TransformedNodeType;

        static const unsigned int partition = IndexStart + Left::TotalCount; 
        static const unsigned int end = IndexStart + Left::TotalCount + Right::TotalCount;
        typedef typename Swap<IndicesType, IndexStart, partition, end>::type TransformedIndicesType;    
    };

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_COMMUTATIVE_TRANSFORM_HPP
