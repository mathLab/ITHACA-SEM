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

#ifndef EXPRESSION_TEMPLATES_COMMUTATIVE_TRAITS_HPP
#define EXPRESSION_TEMPLATES_COMMUTATIVE_TRAITS_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Operators.hpp>
#include <boost/type_traits.hpp>

namespace expt
{
    /// \brief Traits class indicating whether an operator is commutative.
    ///
    /// Given an expression T1 Op T2, where T1 and T2 are data types
    /// in an expression, and Op is the operator, CommutativeTraits initicates whether
    /// the expression can be rewritten using the commutative property.  If
    /// it can, then the expression T1 Op T2 and be rewritten as T2 Op T1.
    ///
    /// This trait defaults to false.  The class should be specialized for
    /// all combinations of operators and data types that are commutative.
    template<typename L, typename Op, typename R, typename enabled = void>
    struct CommutativeTraits : public boost::false_type {};
    
    /// \brief Specialization indicating that addition is usually commutative.
    template<typename R, typename T>
    struct CommutativeTraits<R, AddOp, T> : public boost::true_type {};

    /// \brief Specialization indicating that multiplication is usually commutative.
    template<typename R, typename T>
    struct CommutativeTraits<R, MultiplyOp, T> : public boost::true_type {};
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_COMMUTATIVE_TRAITS_HPP
