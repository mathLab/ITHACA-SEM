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

#ifndef EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
#define EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Operators.hpp>
#include <boost/type_traits.hpp>

namespace expt
{

    /// \brief Traits class indicating whether an operator is associative.
    ///
    /// Given an expression T1 Op1 (T2 Op2 X), where T1 and T2 are data types
    /// in an expression, Op1 and Op2 are the operators in the expression, and
    /// X is an arbitrary expression, AssociativeTraits initicates whether
    /// the expression can be rewritten using the associative property.  If
    /// it can, then the expression T1 Op1 (T2 Op2 X) can be rewritten as
    /// (T1 Op1 T2) Op2 X.
    ///
    /// This trait defaults to false.  The class should be specialized for
    /// all combinations of operators and data types that are associative.
    template<typename T1, typename Op1, typename T2, typename Op2>
    struct AssociativeTraits : public boost::false_type {};

    /// \brief Specialization indicating multiplication is usually associative.
    template<typename T>
    struct AssociativeTraits<T, MultiplyOp, T, MultiplyOp> : public boost::true_type {};

    /// \brief Specialization indicating addition is usually associative.
    template<typename T>
    struct AssociativeTraits<T, AddOp, T, AddOp> : public boost::true_type {};
}

#endif // NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
