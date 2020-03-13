/////////////////////////////////////////////////////////////////////////////////
////
//// File: TestConstantBinarySpecialization0.cpp
////
//// For more information, please see: http://www.nektar.info
////
//// The MIT License
////
//// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//// Department of Aeronautics, Imperial College London (UK), and Scientific
//// Computing and Imaging Institute, University of Utah (USA).
////
//// License for the specific language governing rights and limitations under
//// Permission is hereby granted, free of charge, to any person obtaining a
//// copy of this software and associated documentation files (the "Software"),
//// to deal in the Software without restriction, including without limitation
//// the rights to use, copy, modify, merge, publish, distribute, sublicense,
//// and/or sell copies of the Software, and to permit persons to whom the
//// Software is furnished to do so, subject to the following conditions:
////
//// The above copyright notice and this permission notice shall be included
//// in all copies or substantial portions of the Software.
////
//// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//// DEALINGS IN THE SOFTWARE.
////
//// Description: Tests Constant-Binary specializations which require a temporary
//// for each side.
////
/////////////////////////////////////////////////////////////////////////////////
//
//#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
//#define NEKTAR_USE_EXPRESSION_TEMPLATES
//#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
//
//#include <boost/test/auto_unit_test.hpp>
//#include <boost/test/test_case_template.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <boost/test/unit_test.hpp>
//#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
//#include <UnitTests/CountedObject.h>
//
//namespace Nektar
//{
//
//    namespace TestConstantBinarySpecialization0_Case0
//    {
////        boost::mpl::and_
////            <
////            boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
////            boost::mpl::not_<typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual>,
////            boost::mpl::not_<boost::is_same<LhsResultType, ResultType> >,
////            boost::mpl::not_<boost::is_same<ResultType, typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::LhsAsBinaryExpressionResultType> >
////        >,
//
//    }
//    
//    namespace TestConstantBinarySpecialization0_Case1
//    {
////        <
////        boost::mpl::not_< typename AssociativeTraits<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> >::IsAssociative>,
////        boost::mpl::not_<IsCommutative<ConstantExpressionPolicy<LhsResultType>, OpType, BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType> > >,
////        boost::mpl::not_<typename OpType<LhsResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType>::HasOpLeftEqual>,
////        boost::mpl::not_<boost::is_same<LhsResultType, ResultType> >
//    
//    }
//    
//    //    .ilb LhsIsResultType RhsIsResultType IsCommutative IsAssociative
//    //    0-0-- 10000
//    //    00--- 10000
//    namespace TestConstantBinarySpecialization0_Case1
//    {
////        boost::mpl::and_
////    <
////        boost::mpl::not_<boost::is_same<ResultType, typename ConstantExpressionPolicy<LhsResultType>::ResultType> >,
////        boost::mpl::not_<boost::is_same<ResultType, typename BinaryExpressionPolicy<RhsLhsPolicyType, RhsOpType, RhsRhsPolicyType>::ResultType> >
////    >
//
//        // Need R = A + (B*B), where C=B*B
//        class A : public CountedObject<A>
//        {
//            public:
//                A() : CountedObject<A>() {}
//                A(int v) : CountedObject<A>(v)  {}
//                A(const A& rhs) : CountedObject<A>(rhs) {}
//                A& operator=(const A& rhs) { CountedObject<A>::operator=(rhs); return *this; }
//        };
//        
//        class B : public CountedObject<B>
//        {
//            public:
//                B() : CountedObject<B>() {}
//                B(int v) : CountedObject<B>(v) {}
//                B(const B& rhs) : CountedObject<B>(rhs) {}
//                B& operator=(const B& rhs) { CountedObject<B>::operator=(rhs);  return *this; }
//        };
//        
//        class C : public CountedObject<C>
//        {
//            public:
//                C() : CountedObject<C>() {}
//                C(int v) : CountedObject<C>(v) {}
//                C(const C& rhs) : CountedObject<C>(rhs) {}
//                C& operator=(const C& rhs) { CountedObject<C>::operator=(rhs); return *this; }
//        };
//        
//        class R : public CountedObject<R>
//        {
//            public:
//                R() : CountedObject<R>() {}
//                R(int v) : CountedObject<R>(v) {}
//                R(const R& rhs) : CountedObject<R>(rhs) {}
//                R& operator=(const R& rhs) { CountedObject<R>::operator=(rhs);  return *this; }                
//        };
//        
//        C Multiply(const B& lhs, const B& rhs)
//        {
//            return C(lhs.value*rhs.value);
//        }
//        
//        void Multiply(C& result, const B& lhs, const B& rhs)
//        {
//            result.value = lhs.value*rhs.value;
//        }
//        
//        GENERATE_MULTIPLICATION_OPERATOR(B, 0, B, 0);     
//        
//        R Add(const A& lhs, const C& rhs)
//        {
//            return R(lhs.value + rhs.value);
//        }
//        
//        void Add(R& result, const A& lhs, const C& rhs)
//        {
//            result.value = lhs.value + rhs.value;
//        }
//        
//        GENERATE_ADDITION_OPERATOR(A, 0, C, 0);
//        
//        BOOST_AUTO_TEST_CASE(TestConstantBinarySpecialization0_Case1)
//        {
////            A obj1(10);
////            B obj2(2);
////            B obj3(8);
////            
////            
////            typedef BinaryExpressionPolicy<ConstantExpressionPolicy<B>,
////                                           MultiplyOp,
////                                           ConstantExpressionPolicy<B> > RhsExpressionType;
////            typedef ConstantExpressionPolicy<A> LhsExpressionType;
////            
////            BOOST_STATIC_ASSERT(( BinaryExpressionEvaluator<LhsExpressionType, RhsExpressionType, R, AddOp, BinaryNullOp>::ClassNum == 301 ));
////             
////            Expression<BinaryExpressionPolicy<LhsExpressionType, AddOp, RhsExpressionType> > exp =
////                obj1 + (obj2*obj3);
////
////            R r;
////            
////            CountedObject<A>::ClearCounters();
////            CountedObject<B>::ClearCounters();
////            CountedObject<C>::ClearCounters();
////            
////            Assign(r, exp);
////            BOOST_CHECK_EQUAL(r.value, 10+2*8);
//            
////            CountedObject<A>::Check(1, 0, 1, 0, 0, 1);
////            CountedObject<B>::Check(0, 0, 0, 0, 0, 0);
////            CountedObject<C>::Check(1, 0, 1, 0, 0, 1);
////            CountedObject<R>::Check(0, 0, 0, 0, 0, 0);
//        }
//    }
//}
