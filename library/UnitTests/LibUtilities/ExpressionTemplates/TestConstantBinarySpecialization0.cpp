///////////////////////////////////////////////////////////////////////////////
//
// File: TestConstantBinarySpecialization0.cpp
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
// Description: Tests Constant-Binary specializations which require a temporary
// for each side.
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

    namespace TestConstantBinarySpecialization0_Case0
    {
//        boost::mpl::and_
//            <
//            boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
//            boost::mpl::not_<typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual>,
//            boost::mpl::not_<boost::is_same<LhsResultType, ResultType> >,
//            boost::mpl::not_<boost::is_same<ResultType, typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::LhsAsBinaryExpressionResultType> >
//        >,

    }
    
    namespace TestConstantBinarySpecialization0_Case1
    {
//        <
//        boost::mpl::not_< typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative>,
//        boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
//        boost::mpl::not_<typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual>,
//        boost::mpl::not_<boost::is_same<LhsResultType, ResultType> >
    
    }
    
    namespace TestConstantBinarySpecialization0_Case1
    {
//        boost::mpl::and_
//    <
//        boost::mpl::not_<boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType> >,
//        boost::mpl::not_<boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType> >
//    >

        // Need R = A + (B*B), where C=B*B
        class A
        {
            public:
                A() : Value(0) {}
                A(int v) : Value(v) {}
                int Value;
        };
        
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
        
        C NekMultiply(const B& lhs, const B& rhs)
        {
            return C(lhs.Value*rhs.Value);
        }
        
        void NekMultiply(C& result, const B& lhs, const B& rhs)
        {
            result.Value = lhs.Value*rhs.Value;
        }
        
        GENERATE_MULTIPLICATION_OPERATOR(B, B);     
        
        R NekAdd(const A& lhs, const C& rhs)
        {
            return R(lhs.Value + rhs.Value);
        }
        
        void NekAdd(R& result, const A& lhs, const C& rhs)
        {
            result.Value = lhs.Value + rhs.Value;
        }
        
        GENERATE_ADDITION_OPERATOR(A, C);
        
        BOOST_AUTO_TEST_CASE(TestConstantBinarySpecialization0_Case1)
        {
            A obj1(10);
            B obj2(2);
            B obj3(8);
            
            typedef BinaryExpressionPolicy<ConstantExpressionPolicy<B>,
                                           MultiplyOp,
                                           ConstantExpressionPolicy<B> > RhsExpressionType;
            typedef ConstantExpressionPolicy<A> LhsExpressionType;
            
            BOOST_STATIC_ASSERT(( BinaryExpressionEvaluator<LhsExpressionType, RhsExpressionType, R, AddOp, BinaryNullOp>::ClassNum == 1 ));
             
            Expression<BinaryExpressionPolicy<LhsExpressionType, AddOp, RhsExpressionType> > exp =
                obj1 + (obj2*obj3);

            R r;
            Assign(r, exp);
            BOOST_CHECK_EQUAL(r.Value, 10+2*8);
        }
    }
}