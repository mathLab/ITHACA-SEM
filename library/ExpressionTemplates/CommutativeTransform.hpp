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

namespace expt
{
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

    template<typename NodeType, typename IndicesType, unsigned int IndexStart>
    struct CommutativeTransform;

    template<typename Left, typename Op, typename Right, typename IndicesType, unsigned int IndexStart>
    struct CommutativeTransform<Node<Left, Op, Right>, IndicesType, IndexStart >
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
