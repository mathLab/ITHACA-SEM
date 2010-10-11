/////////////////////////////////////////////////////////////////////////////////
////
//// File: TestBinaryExpressionEvaluator_GeneralSpecialization_ParentsOpType1.cpp
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
//// Description:
////
/////////////////////////////////////////////////////////////////////////////////
//
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
//    namespace GeneralSpecializationParentOp1
//    {
//       class A : public CountedObject<A>
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
//        A Add(const A& lhs, const A& rhs)
//        {
//            return A(lhs.value + rhs.value);
//        }
//        
//        void Add(A& result, const A& lhs, const A& rhs)
//        {
//            result.value = lhs.value+rhs.value;
//        }
//        
//        void AddEqual(A& result, const A& rhs)
//        {
//            result.value += rhs.value;
//        }
//        
//        GENERATE_ADDITION_OPERATOR(A, 0, A, 0);     
//        
//        R Add(const R& lhs, const R& rhs)
//        {
//            return R(lhs.value + rhs.value);
//        }
//        
//        void Add(R& result, const R& lhs, const R& rhs)
//        {
//            result.value = lhs.value+rhs.value;
//        }
//        
//        void AddEqual(R& result, const R& rhs)
//        {
//            result.value += rhs.value;
//        }
//        
//        GENERATE_ADDITION_OPERATOR(R, 0, R, 0);    
//
//        R Add(const A& lhs, const R& rhs)
//        {
//            return R(lhs.value + rhs.value);
//        }
//        
//        void Add(R& result, const A& lhs, const R& rhs)
//        {
//            result.value = lhs.value + rhs.value;
//        }
//    
//        
//        GENERATE_ADDITION_OPERATOR(A, 0, R, 0);
//        
//        R Add(const R& lhs, const A& rhs)
//        {
//            return R(lhs.value + rhs.value);
//        }
//        
//        void Add(R& result, const R& lhs, const A& rhs)
//        {
//            result.value = lhs.value + rhs.value;
//        }
//        
//        void AddEqual(R& result, const A& rhs)
//        {
//            result.value += rhs.value;
//        }
//        
//        GENERATE_ADDITION_OPERATOR(R, 0, A, 0);
//
//        BOOST_AUTO_TEST_CASE(TestNoOpChange)
//        {
//#ifndef _WIN32
//            // Need an expression with non-constant leaves.
//            // Rhs has the same type as the result
//            // Lhs does not
//            // Rhs is associative.
//            
//            //R = R + ( (A+A) + (R-R) )
//            
//            A obj1(7);
//            A obj2(19);
//            
//            R obj3(-3);
//            R obj4(7);
//            
//            typedef BinaryExpressionPolicy<ConstantExpressionPolicy<A>, AddOp, ConstantExpressionPolicy<A> > LhsPolicy;
//            typedef BinaryExpressionPolicy<ConstantExpressionPolicy<R>, AddOp, ConstantExpressionPolicy<R> > RhsPolicy;
//            
//            Expression<LhsPolicy> lhsExp = obj1+obj2;
//            Expression<RhsPolicy> rhsExp = obj3 + obj4;
//                                    
//            typedef BinaryExpressionEvaluator<LhsPolicy, RhsPolicy, R, AddOp, AddOp> EvalType;
//            
//            BOOST_MPL_ASSERT(( boost::mpl::not_<boost::is_same<R, LhsPolicy::ResultType> >));
//            BOOST_MPL_ASSERT(( boost::is_same<R, RhsPolicy::ResultType> ));
//            BOOST_MPL_ASSERT(( AssociativeTraits<ConstantExpressionPolicy<R>, AddOp, RhsPolicy>::IsAssociative ));
//            BOOST_MPL_ASSERT(( boost::mpl::and_
//                                            <
//                                                IsBinaryExpressionPolicy<LhsPolicy>,
//                                                IsBinaryExpressionPolicy<RhsPolicy>
//                                            > ));
//            BOOST_MPL_ASSERT(( boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, AddOp<int, int> > > ));
//            BOOST_STATIC_ASSERT( EvalType::ClassNum == 201 );
//            //
//            BOOST_MPL_ASSERT(( boost::mpl::and_
//                                        <
//                                            boost::mpl::and_
//                                            <
//                                                IsBinaryExpressionPolicy<LhsPolicy>,
//                                                IsBinaryExpressionPolicy<RhsPolicy>
//                                            >,
//                                            boost::mpl::not_<boost::is_same<BinaryNullOp<int, int>, AddOp<int, int> > >,
//                                            boost::mpl::or_
//                                            <
//                                                boost::mpl::and_
//                                                <
//                                                    boost::mpl::not_<boost::is_same<R, LhsPolicy::ResultType> >,
//                                                    boost::is_same<R, RhsPolicy::ResultType>,
//                                                    AssociativeTraits<ConstantExpressionPolicy<R>, AddOp, RhsPolicy>::IsAssociative
//                                                >,
//                                                boost::mpl::and_
//                                                <
//                                                     boost::is_same<R, RhsPolicy::ResultType>,
//                                                     boost::mpl::not_<AssociativeTraits<ConstantExpressionPolicy<R>, AddOp, LhsPolicy>::IsAssociative>,
//                                                     AssociativeTraits<ConstantExpressionPolicy<R>, AddOp, RhsPolicy>::IsAssociative
//                                                >
//                                            >
//                                        > ));
//
//            R result;
//            Accumulator<R> accum(result);
//#endif
//            //EvalType::Eval(lhsExp, rhsExp, accum);
//            //
//            //BOOST_CHECK_EQUAL(result.value, 7+19-3+7);
//
//
//
//
//
//
//
////            // R = R + (R + (B+C))
////            R obj1(8);
////            R obj2(-3);
////            B obj3(23);
////            C obj4(1);
////            
////            R result;
////            Assign(result, obj1 + (obj2 + (obj3+obj4)));
////            BOOST_CHECK_EQUAL(8 - 3 + 23 + 1, result.Value);
//        }
//    }
//}
