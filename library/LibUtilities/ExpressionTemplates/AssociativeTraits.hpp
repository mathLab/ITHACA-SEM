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

#ifndef NEKTAR_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
#define NEKTAR_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP

#include "Operators.hpp"
#include <boost/type_traits.hpp>

namespace Nektar
{

    /// \brief A traits class to indicate if an expression is associative or not.
    ///
    /// Given an expression T1 Op1 (R Op2 X), if AssociativeTraits inherits from 
    /// boost::true_type then this expression can be rewritten as 
    /// (T1 Op1 R) Op2 X, otherwise it can't be rewritten.
    ///
    /// For example, the matrix expression M1 + (M2 + M3) can be rewritten as 
    /// (M1 + M2) + M3.
    template<typename T1, typename Op1, typename R, typename Op2>
    struct AssociativeTraits : public boost::false_type
    {
    };

    /// \brief Specialization indicating multiplication is usually associative.
    template<typename T1, typename T2>
    struct AssociativeTraits<T1, MultiplyOp, T2, MultiplyOp> : public boost::true_type
    {
    };

    /// \brief Specialization indicating addition is usually associative.
    template<typename T1, typename T2>
    struct AssociativeTraits<T1, AddOp, T2, AddOp> : public boost::true_type
    {
    };


    template<typename NodeType>
    struct TreeIsAssociative;// : public boost::false_type {};

    template<typename L, typename RootOp, typename R1, typename ROp>
    struct TreeIsAssociative<Node<L, RootOp, Node<R1, ROp, void> > > : public boost::true_type
    {
    };

    template<typename L, typename RootOp, typename R1>
    struct TreeIsAssociative<Node<L, RootOp, Node<R1, void, void> > > : public boost::true_type
    {
    };

    template<typename L, typename RootOp, typename R1, typename ROp, typename R2>
    struct TreeIsAssociative<Node<L, RootOp, Node<R1, ROp, R2> > > :
        public boost::mpl::if_
        <
            boost::mpl::and_
            <
                AssociativeTraits<typename L::ResultType, RootOp, typename R1::ResultType, ROp>,
                TreeIsAssociative<Node<R1, ROp, R2> >
            >,
            boost::true_type,
            boost::false_type
        >::type
    {
    };
}

#endif //NEKTAR_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
