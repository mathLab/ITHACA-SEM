///////////////////////////////////////////////////////////////////////////////
//
// File: TestConstantBinarySpecialization0_ParentOpType.cpp
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


#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>

namespace Nektar
{
    namespace ConstantBinaryExpressionSpecialization_ParentOpType0
    {
//        template<typename LhsResultType, typename RhsLhsPolicyType, typename RhsRhsPolicyType,
//                 template <typename, typename> class RhsOpType,
//                 template <typename, typename> class OpType,
//                 template <typename, typename> class ParentOpType,
//                 typename ResultType>
//        struct BinaryExpressionEvaluator<ConstantExpressionPolicy<LhsResultType>,
//                                         BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>,
//                                         ResultType, OpType, ParentOpType,
//                                         typename boost::enable_if
//                                         <
//                                            boost::mpl::and_
//                                            <
//                                                typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative,
//                                                boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, ParentOpType<int, int> > >
//                                            >
//                                         >::type >
//        {
//            static void Eval(const Expression<ConstantExpressionPolicy<LhsResultType> >& lhs,
//                             const Expression<BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >& rhs, 
//                             Accumulator<ResultType>& accum)
//            {
//                typedef typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType RhsResultType;
//                typedef Expression<ConstantExpressionPolicy<LhsResultType> > LhsExpressionType;
//                typedef Expression<ConstantExpressionPolicy<RhsResultType> > RhsExpressionType;
//
//                ParentOpType<ResultType, LhsResultType>::ApplyEqual(accum, *lhs);
//                typedef AssociativeTraits<ConstantExpressionPolicy<ResultType>, ParentOpType,
//                                          BinaryExpressionPolicy<LhsExpressionType, OpType, RhsExpressionType> > TraitsType;
//
//                typedef typename TraitsType::template LhsOpType<ResultType, typename RhsLhsPolicyType::ResultType, ParentOpType>::type NewLhsOpType;
//                typedef typename TraitsType::template RhsOpType<ResultType, typename RhsRhsPolicyType::ResultType, ParentOpType>::type NewRhsOpType;
//                
//                NewLhsOpType::ApplyEqual(accum, (*rhs).first);
//                NewRhsOpType::ApplyEqual(accum, (*rhs).second);
//            }
//        };

    // Specialization 0 takes the form of A + (B + (C+D)).  The specialization 
    // deals with the evaluation of the (B + (C+D)) subexpression when the 
    // accumulator has already been initialized and set to the value of A.
    //
    // This specialization is called from constant-binary specialization 4
//                                            typename boost::enable_if
//                                        <
//                                            boost::mpl::and_
//                                            <
//                                                typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative,
//                                                boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType>,
//                                                typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::OpEqualsAreDefined
//                                            >
//                                         >::type
//                                         >

        class B
        {
            public:
                B() : Value(0) {}
                B(int v) : Value(v) {}
                int Value;
        };
        
        class C
        {
            public:
                C() : Value(0) {}
                C(int v) : Value(v) {}
                int Value;
        };
        
        class R
        {
            public:
                R() : Value(0) {}
                R(int v) : Value(v) {}
                int Value;
        };
        
        R NekAdd(const B& lhs, const C& rhs)
        {
            return R(lhs.Value*rhs.Value);
        }
        
        void NekAdd(R& result, const B& lhs, const C& rhs)
        {
            result.Value = lhs.Value+rhs.Value;
        }
        
        R NekAdd(const R& lhs, const B& rhs)
        {
            return R(lhs.Value + rhs.Value);
        }
        
        R NekAdd(const R& lhs, const C& rhs)
        {
            return R(lhs.Value + rhs.Value);
        }
        
        void NekAdd(R& result, const R& lhs, const B& rhs)
        {
            result.Value = lhs.Value + rhs.Value;
        }
        
        void NekAdd(R& result, const R& lhs, const C& rhs)
        {
            result.Value = lhs.Value + rhs.Value;
        }
        
        void NekAddEqual(R& result, const B& rhs)
        {
            result.Value += rhs.Value;
        }
        
        void NekAddEqual(R& result, const C& rhs)
        {
            result.Value += rhs.Value;
        }
            
        GENERATE_ADDITION_OPERATOR(B, 0, C, 0);     
        
        R NekAdd(const R& lhs, const R& rhs)
        {
            return R(lhs.Value + rhs.Value);
        }
        
        void NekAdd(R& result, const R& lhs, const R& rhs)
        {
            result.Value = lhs.Value + rhs.Value;
        }
        
        GENERATE_ADDITION_OPERATOR(R, 0, R, 0);
        GENERATE_ADDITION_OPERATOR(R, 0, B, 0);
        GENERATE_ADDITION_OPERATOR(R, 0, C, 0);
        
        R NekSubtract(const R& lhs, const R& rhs)
        {
            return R(lhs.Value - rhs.Value);
        }
        
        void NekSubtract(R& result, const R& lhs, const R& rhs)
        {
            result.Value = lhs.Value - rhs.Value;
        }
        
        GENERATE_SUBTRACTION_OPERATOR(R, 0, R, 0);
        
        
        R NekSubtract(const R& lhs, const B& rhs)
        {
            return R(lhs.Value - rhs.Value);
        }
        
        void NekSubtract(R& result, const R& lhs, const B& rhs)
        {
            result.Value = lhs.Value - rhs.Value;
        }
        
        void NekSubtractEqual(R& result, const B& rhs)
        {
            result.Value -= rhs.Value;
        }
        
        GENERATE_SUBTRACTION_OPERATOR(R, 0, B, 0);
        
        R NekSubtract(const R& lhs, const C& rhs)
        {
            return R(lhs.Value - rhs.Value);
        }
        
        void NekSubtract(R& result, const R& lhs, const C& rhs)
        {
            result.Value = lhs.Value - rhs.Value;
        }
        
        void NekSubtractEqual(R& result, const C& rhs)
        {
            result.Value -= rhs.Value;
        }
        
        GENERATE_SUBTRACTION_OPERATOR(R, 0, C, 0);
            
        BOOST_AUTO_TEST_CASE(TestNoOpChange)
        {
//            // R = R + (R + (B+C))
//            R obj1(8);
//            R obj2(-3);
//            B obj3(23);
//            C obj4(1);
//            
//            R result;
//            Assign(result, obj1 + (obj2 + (obj3+obj4)));
//            BOOST_CHECK_EQUAL(8 - 3 + 23 + 1, result.Value);
        }
    }
}
