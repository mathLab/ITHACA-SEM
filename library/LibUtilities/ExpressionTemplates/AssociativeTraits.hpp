///////////////////////////////////////////////////////////////////////////////
//
// File: AssociativeTraits.hpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicyFwd.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionEvaluatorFwd.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

namespace Nektar
{  
    template<typename LhsResultType,
             template<typename, typename> class LhsOpType,
             typename RhsLhsPolicyType,
             template<typename, typename> class RhsOpType,
             typename RhsRhsPolicyType>
    struct ConstantBinaryCommonAssociativeTraits
    {
        typedef BinaryExpressionPolicy<ConstantExpressionPolicy<LhsResultType>, LhsOpType, RhsLhsPolicyType> LhsAsBinaryExpression;
        typedef typename LhsAsBinaryExpression::ResultType LhsAsBinaryExpressionResultType;
        typedef boost::true_type IsAssociative;
        typedef boost::true_type IsLeftAssociative;
        typedef boost::true_type IsRightAssociative;
    };
    
    template<typename LhsPolicy,
             template <typename, typename> class Op, 
             typename RhsPolicy>
    struct AssociativeTraits 
    {
        typedef typename LhsPolicy::ResultType LhsAsBinaryExpressionResultType;
        typedef boost::false_type IsAssociative;
        typedef boost::false_type OpEqualsAreDefined;
        typedef boost::false_type IsLeftAssociative;
        typedef boost::false_type IsRightAssociative;
    };
    
    template<typename LhsResultType,
             typename RhsLhsPolicyType,
             typename RhsRhsPolicyType>
    struct AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, AddOp, BinaryExpressionPolicy<RhsLhsPolicyType, AddOp, RhsRhsPolicyType> > 
        : public ConstantBinaryCommonAssociativeTraits<LhsResultType, AddOp, RhsLhsPolicyType, AddOp, RhsRhsPolicyType>
    {
                
        template<typename Lhs, typename Rhs, template <typename, typename> class ParentOpType = AddOp>
        struct LhsOpType
        {
            typedef AddOp<Lhs, Rhs> type;
        };
        
        template<typename Lhs, typename Rhs>
        struct LhsOpType<Lhs, Rhs, SubtractOp>
        {
            typedef SubtractOp<Lhs, Rhs> type;
        };
        
        template<typename Lhs, typename Rhs, template<typename, typename> class ParentOpType = AddOp>
        struct RhsOpType
        {
            typedef AddOp<Lhs, Rhs> type;
            
        };
        
        template<typename Lhs, typename Rhs>
        struct RhsOpType<Lhs, Rhs, SubtractOp>
        {
            typedef SubtractOp<Lhs, Rhs> type;
            
        };
        
        struct OpEqualsAreDefined : public boost::mpl::if_
                                           <
                                                boost::mpl::and_
                                                <
                                                    HasOpEqualTraits<LhsResultType, typename RhsLhsPolicyType::ResultType, AddOp>,
                                                    HasOpEqualTraits<LhsResultType, typename RhsRhsPolicyType::ResultType, AddOp>
                                                >,
                                                boost::true_type,
                                                boost::false_type
                                           >::type
        {
        };
        
        static void EvaluateNewRhs(const Expression<RhsRhsPolicyType>& exp, Accumulator<typename RhsRhsPolicyType::ResultType>& accum)
        {
            exp.template Evaluate<AddOp>(accum);
        }
    };
    
    template<typename LhsResultType,
             typename RhsLhsPolicyType,
             typename RhsRhsPolicyType>
    struct AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, AddOp, BinaryExpressionPolicy<RhsLhsPolicyType, SubtractOp, RhsRhsPolicyType> > 
            : public ConstantBinaryCommonAssociativeTraits<LhsResultType, AddOp, RhsLhsPolicyType, SubtractOp, RhsRhsPolicyType>
    {
        template<typename Lhs, typename Rhs, template <typename, typename> class ParentOpType = AddOp>
        struct LhsOpType
        {
            typedef AddOp<Lhs, Rhs> type;
        };
        
        template<typename Lhs, typename Rhs>
        struct LhsOpType<Lhs, Rhs, SubtractOp>
        {
            typedef SubtractOp<Lhs, Rhs> type;
        };
        
        template<typename Lhs, typename Rhs, template<typename, typename> class ParentOpType = AddOp>
        struct RhsOpType
        {
            typedef SubtractOp<Lhs, Rhs> type;
            
        };
        
        template<typename Lhs, typename Rhs>
        struct RhsOpType<Lhs, Rhs, SubtractOp>
        {
            typedef AddOp<Lhs, Rhs> type;
            
        };
        
        struct OpEqualsAreDefined : public boost::mpl::if_
                                           <
                                                boost::mpl::and_
                                                <
                                                    HasOpEqualTraits<LhsResultType, typename RhsLhsPolicyType::ResultType, AddOp>,
                                                    HasOpEqualTraits<LhsResultType, typename RhsRhsPolicyType::ResultType, SubtractOp>
                                                >,
                                                boost::true_type,
                                                boost::false_type
                                           >::type
        {
        };
        
        static void EvaluateNewRhs(const Expression<RhsRhsPolicyType>& exp, Accumulator<typename RhsRhsPolicyType::ResultType>& accum)
        {
            exp.template Evaluate<SubtractOp>(accum);
        }
    };
    
    template<typename LhsResultType,
             typename RhsLhsPolicyType,
             typename RhsRhsPolicyType>
    struct AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, SubtractOp, BinaryExpressionPolicy<RhsLhsPolicyType, AddOp, RhsRhsPolicyType> >
            : public ConstantBinaryCommonAssociativeTraits<LhsResultType, SubtractOp, RhsLhsPolicyType, AddOp, RhsRhsPolicyType>
    {
        
        template<typename Lhs, typename Rhs, template <typename, typename> class ParentOpType = AddOp>
        struct LhsOpType
        {
            typedef SubtractOp<Lhs, Rhs> type;
        };
        
        template<typename Lhs, typename Rhs>
        struct LhsOpType<Lhs, Rhs, SubtractOp>
        {
            typedef AddOp<Lhs, Rhs> type;
        };


        template<typename Lhs, typename Rhs, template <typename, typename> class ParentOpType = AddOp>
        struct RhsOpType
        {
            typedef SubtractOp<Lhs, Rhs> type;
            
        };
        
        template<typename Lhs, typename Rhs>
        struct RhsOpType<Lhs, Rhs, SubtractOp>
        {
            typedef AddOp<Lhs, Rhs> type;
            
        };
        
        struct OpEqualsAreDefined : public boost::mpl::if_
                                           <
                                                boost::mpl::and_
                                                <
                                                    HasOpEqualTraits<LhsResultType, typename RhsLhsPolicyType::ResultType, SubtractOp>,
                                                    HasOpEqualTraits<LhsResultType, typename RhsRhsPolicyType::ResultType, SubtractOp>
                                                >,
                                                boost::true_type,
                                                boost::false_type
                                           >::type
        {
        };
        
        static void EvaluateNewRhs(const Expression<RhsRhsPolicyType>& exp, Accumulator<typename RhsRhsPolicyType::ResultType>& accum)
        {
            exp.template Evaluate<SubtractOp>(accum);
        }
    };
    
    template<typename LhsResultType,
             typename RhsLhsPolicyType,
             typename RhsRhsPolicyType>
    struct AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, SubtractOp, BinaryExpressionPolicy<RhsLhsPolicyType, SubtractOp, RhsRhsPolicyType> >
            : public ConstantBinaryCommonAssociativeTraits<LhsResultType, SubtractOp, RhsLhsPolicyType, SubtractOp, RhsRhsPolicyType>
    {
        
        template<typename Lhs, typename Rhs, template <typename, typename> class ParentOpType = AddOp>
        struct LhsOpType
        {
            typedef SubtractOp<Lhs, Rhs> type;
        };
        
        template<typename Lhs, typename Rhs>
        struct LhsOpType<Lhs, Rhs, SubtractOp>
        {
            typedef AddOp<Lhs, Rhs> type;
        };


        template<typename Lhs, typename Rhs, template <typename, typename> class ParentOpType = AddOp>
        struct RhsOpType
        {
            typedef AddOp<Lhs, Rhs> type;
            
        };
        
        template<typename Lhs, typename Rhs>
        struct RhsOpType<Lhs, Rhs, SubtractOp>
        {
            typedef SubtractOp<Lhs, Rhs> type;
            
        };
        
        struct OpEqualsAreDefined : public boost::mpl::if_
                                           <
                                                boost::mpl::and_
                                                <
                                                    HasOpEqualTraits<LhsResultType, typename RhsLhsPolicyType::ResultType, SubtractOp>,
                                                    HasOpEqualTraits<LhsResultType, typename RhsRhsPolicyType::ResultType, AddOp>
                                                >,
                                                boost::true_type,
                                                boost::false_type
                                           >::type
        {
        };
        
        static void EvaluateNewRhs(const Expression<RhsRhsPolicyType>& exp, Accumulator<typename RhsRhsPolicyType::ResultType>& accum)
        {
            exp.template Evaluate<SubtractOp>(accum);
        }
    };
    
    template<typename LhsResultType,
             typename RhsLhsPolicyType,
             typename RhsRhsPolicyType>
    struct AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, MultiplyOp, BinaryExpressionPolicy<RhsLhsPolicyType, MultiplyOp, RhsRhsPolicyType> >
        : public ConstantBinaryCommonAssociativeTraits<LhsResultType, MultiplyOp, RhsLhsPolicyType, MultiplyOp, RhsRhsPolicyType>
    {
        
        template<typename Lhs, typename Rhs>
        struct RhsOpType
        {
            typedef MultiplyOp<Lhs, Rhs> type;
        };
        
        struct OpEqualsAreDefined : public boost::mpl::if_
                                           <
                                                boost::mpl::and_
                                                <
                                                    HasOpEqualTraits<LhsResultType, typename RhsLhsPolicyType::ResultType, MultiplyOp>,
                                                    HasOpEqualTraits<LhsResultType, typename RhsRhsPolicyType::ResultType, MultiplyOp>
                                                >,
                                                boost::true_type,
                                                boost::false_type
                                           >::type
        {
        };
        
        static void EvaluateNewRhs(const Expression<RhsRhsPolicyType>& exp, Accumulator<typename RhsRhsPolicyType::ResultType>& accum)
        {
            exp.template Evaluate<MultiplyOp>(accum);
        }
    };
    
    template<typename LhsResultType,
             typename RhsLhsPolicyType,
             typename RhsRhsPolicyType>
    struct AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, MultiplyOp, BinaryExpressionPolicy<RhsLhsPolicyType, DivideOp, RhsRhsPolicyType> >
        : public ConstantBinaryCommonAssociativeTraits<LhsResultType, MultiplyOp, RhsLhsPolicyType, DivideOp, RhsRhsPolicyType>
    {
        
        template<typename Lhs, typename Rhs>
        struct RhsOpType
        {
            typedef DivideOp<Lhs, Rhs> type;
        };
        
        struct OpEqualsAreDefined : public boost::mpl::if_
                                           <
                                                boost::mpl::and_
                                                <
                                                    HasOpEqualTraits<LhsResultType, typename RhsLhsPolicyType::ResultType, MultiplyOp>,
                                                    HasOpEqualTraits<LhsResultType, typename RhsRhsPolicyType::ResultType, MultiplyOp>
                                                >,
                                                boost::true_type,
                                                boost::false_type
                                           >::type
        {
        };
        
        static void EvaluateNewRhs(const Expression<RhsRhsPolicyType>& exp, Accumulator<typename RhsRhsPolicyType::ResultType>& accum)
        {
            exp.template Evaluate<DivideOp>(accum);
        }
    };
    
        template<typename LhsResultType,
             typename RhsLhsPolicyType,
             typename RhsRhsPolicyType>
    struct AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, DivideOp, BinaryExpressionPolicy<RhsLhsPolicyType, MultiplyOp, RhsRhsPolicyType> >
        : public ConstantBinaryCommonAssociativeTraits<LhsResultType, DivideOp, RhsLhsPolicyType, MultiplyOp, RhsRhsPolicyType>
    {
        
        template<typename Lhs, typename Rhs>
        struct RhsOpType
        {
            typedef DivideOp<Lhs, Rhs> type;
        };
        
        struct OpEqualsAreDefined : public boost::mpl::if_
                                           <
                                                boost::mpl::and_
                                                <
                                                    HasOpEqualTraits<LhsResultType, typename RhsLhsPolicyType::ResultType, MultiplyOp>,
                                                    HasOpEqualTraits<LhsResultType, typename RhsRhsPolicyType::ResultType, MultiplyOp>
                                                >,
                                                boost::true_type,
                                                boost::false_type
                                           >::type
        {
        };
        
        static void EvaluateNewRhs(const Expression<RhsRhsPolicyType>& exp, Accumulator<typename RhsRhsPolicyType::ResultType>& accum)
        {
            exp.template Evaluate<DivideOp>(accum);
        }
    };
    
    template<typename LhsResultType,
             typename RhsLhsPolicyType,
             typename RhsRhsPolicyType>
    struct AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, DivideOp, BinaryExpressionPolicy<RhsLhsPolicyType, DivideOp, RhsRhsPolicyType> >
        : public ConstantBinaryCommonAssociativeTraits<LhsResultType, DivideOp, RhsLhsPolicyType, DivideOp, RhsRhsPolicyType>
    {
        
        template<typename Lhs, typename Rhs>
        struct RhsOpType
        {
            typedef MultiplyOp<Lhs, Rhs> type;
        };
        
        struct OpEqualsAreDefined : public boost::mpl::if_
                                           <
                                                boost::mpl::and_
                                                <
                                                    HasOpEqualTraits<LhsResultType, typename RhsLhsPolicyType::ResultType, MultiplyOp>,
                                                    HasOpEqualTraits<LhsResultType, typename RhsRhsPolicyType::ResultType, MultiplyOp>
                                                >,
                                                boost::true_type,
                                                boost::false_type
                                           >::type
        {
        };
        
        static void EvaluateNewRhs(const Expression<RhsRhsPolicyType>& exp, Accumulator<typename RhsRhsPolicyType::ResultType>& accum)
        {
            exp.template Evaluate<MultiplyOp>(accum);
        }
    };
}

#endif // NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
