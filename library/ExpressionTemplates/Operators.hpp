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


#ifndef EXPRESSION_TEMPLATES_OPERATORS_HPP
#define EXPRESSION_TEMPLATES_OPERATORS_HPP

#include <boost/typeof/typeof.hpp>

namespace expt
{
    struct AddOp 
    {
        template<typename L, typename R>
        struct ResultType
        {
            BOOST_TYPEOF_NESTED_TYPEDEF_TPL(nested, Add(L(), R()));
            typedef typename nested::type type;
        };
        
        template<typename L, typename R>
        static void OpEqual(L& accum, const R& rhs)
        {
            AddEqual(accum, rhs);
        }

        static void OpEqual(double& accum, const double& rhs)
        {
            accum += rhs;
        }

        template<typename L, typename R>
        static typename ResultType<L, R>::type Op(const L& lhs, const R& rhs)
        {
            return Add(lhs, rhs);
        }

        template<typename L, typename R>
        static void Op(typename ResultType<L, R>::type& accumulator, const L& lhs, const R& rhs)
        {
            return Add(accumulator, lhs, rhs);
        }
    };

    struct MultiplyOp 
    {
        template<typename L, typename R>
        class ResultType
        {
            public:
                BOOST_TYPEOF_NESTED_TYPEDEF_TPL(nested, Multiply(L(), R()));
                typedef typename nested::type type;
        };
        
        template<typename L, typename R>
        static void OpEqual(L& accum, const R& rhs)
        {
            MultiplyEqual(accum, rhs);
        }

        static void OpEqual(double& accum, const double& rhs)
        {
            accum *= rhs;
        }

        template<typename L, typename R>
        static typename ResultType<L, R>::type Op(const L& lhs, const R& rhs)
        {
            return Multiply(lhs, rhs);
        }

        template<typename L, typename R>
        static void Op(typename ResultType<L, R>::type& accumulator, const L& lhs, const R& rhs)
        {
            return Multiply(accumulator, lhs, rhs);
        }
    };

    struct DivideOp 
    {
        template<typename L, typename R>
        struct ResultType
        {
            BOOST_TYPEOF_NESTED_TYPEDEF_TPL(nested, Divide(L(), R()));
            typedef typename nested::type type;
        };

        template<typename L, typename R>
        static void OpEqual(L& accum, const R& rhs)
        {
            DivideEqual(accum, rhs);
        }

        static void OpEqual(double& accum, const double& rhs)
        {
            accum /= rhs;
        }

        template<typename L, typename R>
        static typename ResultType<L, R>::type Op(const L& lhs, const R& rhs)
        {
            return Divide(lhs, rhs);
        }

        template<typename L, typename R>
        static void Op(typename ResultType<L, R>::type& accumulator, const L& lhs, const R& rhs)
        {
            return Divide(accumulator, lhs, rhs);
        }
    };

    struct SubtractOp 
    {
        template<typename L, typename R>
        struct ResultType
        {
            BOOST_TYPEOF_NESTED_TYPEDEF_TPL(nested, Subtract(L(), R()));
            typedef typename nested::type type;
        };

        template<typename L, typename R>
        static void OpEqual(L& accum, const R& rhs)
        {
            SubtractEqual(accum, rhs);
        }

        static void OpEqual(double& accum, const double& rhs)
        {
            accum -= rhs;
        }

        template<typename L, typename R>
        static typename ResultType<L, R>::type Op(const L& lhs, const R& rhs)
        {
            return Subtract(lhs, rhs);
        }

        template<typename L, typename R>
        static void Op(typename ResultType<L, R>::type& accumulator, const L& lhs, const R& rhs)
        {
            return Subtract(accumulator, lhs, rhs);
        }

    };

    struct NegateOp
    {
        template<typename T>
        struct ResultType
        {
            typedef T type;
        };

        template<typename T>
        static void Op(T& accumulator)
        {
            NegateInPlace(accumulator);
        }
    };
}

#endif //NEKTAR_EXPRESSION_TEMPLATES_OPERATORS_HPP
