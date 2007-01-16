///////////////////////////////////////////////////////////////////////////////
//
// File: testExpressionTemplates.cpp
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

#include <UnitTests/testExpressionTemplates.h>
#include <UnitTests/CountedObject.h>

#include <LibUtilities/LinearAlgebra/NekPoint.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>

#include <iostream>

using std::cout;
using std::endl;

namespace Nektar
{
    namespace UnitTests
    {
        using namespace Nektar;

        void testConstantExpressions()
        {
            using namespace Nektar;
            using namespace Nektar::expt;

            {
                int t = 7;
                Expression<ConstantExpressionPolicy<int> > e1(t);
                BOOST_CHECK_EQUAL(*e1, t);

                Expression<ConstantExpressionPolicy<int> > e2(e1);
                BOOST_CHECK_EQUAL(*e2, *e1);
            }

            typedef NekPoint<unsigned int, 3> Point;
            {
                Point p(1,2,3);
                Expression<ConstantExpressionPolicy<Point> > e1(p);

                Point p1(e1);
                BOOST_CHECK_EQUAL(p, p1);

                Point p2 = e1;
                BOOST_CHECK_EQUAL(p, p2);

                Point p3(9, 10, 11);
                BOOST_CHECK(p != p3);
                p3 = e1;
                BOOST_CHECK_EQUAL(p, p3);
            }

            {
                // TODO - Find a way to prevent temporaries (meaning that the parameter to
                // this call is temporary and that could cause problems).
                Nektar::expt::Expression<ConstantExpressionPolicy<Point> > e2(Point(1,2,3));
            }
        }

        void testUnaryExpressions()
        {
//             using namespace Nektar;
// 
//             NekPoint<double, 3> p(1,2,3);
//             NekPoint<double, 3> p1(-(-p));
// 
//             BOOST_CHECK_EQUAL(p, p1);
// 
//             NekMatrix<double, eFull> m(3,3);
//             NekMatrix<double, eFull> m1(-m);
        }

        void testNekMatrixMetadata()
        {
            using namespace Nektar;
            using namespace Nektar::expt;

            // Constant
            NekMatrix<double> m(3,3);
            ConstantExpressionTraits<NekMatrix<double> >::MetadataType t(m);
            expt::ExpressionMetadataChooser<ConstantExpressionTraits<NekMatrix<double> > >::MetadataType t10(m);
            ConstantExpressionPolicy<NekMatrix<double, eFull, eNormal, 0, void> >::MetadataType t1(m);
            Expression<ConstantExpressionPolicy<NekMatrix<double> > >::MetadataType t2(m);
            Expression<ConstantExpressionPolicy<NekMatrix<double> > > m_exp(m);
            unsigned int i = t.Rows;
            unsigned int j = t1.Rows;
            unsigned int k = t2.Rows;
            unsigned int l = t10.Rows;
            BOOST_CHECK_EQUAL(m.GetRows(), m_exp.GetMetadata().Rows);
            BOOST_CHECK_EQUAL(m.GetColumns(), m_exp.GetMetadata().Columns);

            // Unary
//             NekMatrix<double> m1(3,3);
//             Expression<UnaryExpressionPolicy<Expression<ConstantExpressionPolicy<NekMatrix<double> > >, NegateOp> > m1_exp = 
//                 Expression<UnaryExpressionPolicy<Expression<ConstantExpressionPolicy<NekMatrix<double> > >, NegateOp> >(
//                 Expression<ConstantExpressionPolicy<NekMatrix<double> > >(m1));
//             BOOST_CHECK_EQUAL(m1.GetRows(), m1_exp.GetMetadata().Rows);
//             BOOST_CHECK_EQUAL(m1.GetColumns(), m1_exp.GetMetadata().Columns);

            // Binary
            // TODO
         }
 
         void testBinaryExpressions()
         {
            using namespace Nektar;
            using namespace Nektar::expt;

            unsigned int buf1[] = {1, 2, 3, 4};
            unsigned int buf2[] = {3, 4, 5, 6};
// 
//             // The following crashes because of the temporary.  What to do here?
//             //Nektar::ConstantExpression<Nektar::LibUtilities::NekMatrix<unsigned int> > m1(Nektar::LibUtilities::NekMatrix<unsigned int>(2, 2, buf1));
//             //Nektar::ConstantExpression<Nektar::LibUtilities::NekMatrix<unsigned int> > m2(Nektar::LibUtilities::NekMatrix<unsigned int>(2, 2, buf2));
// 
            NekMatrix<unsigned int> lhs(2,2,buf1);
            NekMatrix<unsigned int> rhs(2,2,buf2);
            Expression<ConstantExpressionPolicy<Nektar::NekMatrix<unsigned int> > > m1(lhs);
            Expression<ConstantExpressionPolicy<Nektar::NekMatrix<unsigned int> > > m2(rhs);


            typedef ConstantExpressionPolicy<NekMatrix<unsigned int> > LhsPolicy;
            typedef ConstantExpressionPolicy<NekMatrix<unsigned int> > RhsPolicy;
            typedef BinaryExpressionPolicy<LhsPolicy, RhsPolicy, AddOp> BinaryPolicy;
            
            BinaryPolicy::DataType d(m1, m2);
            Expression<BinaryPolicy> bexp(d);

            unsigned int result_buf[] = {4, 6, 8, 10};
            NekMatrix<unsigned int> result(2, 2, result_buf);
            NekMatrix<unsigned int> evaluated_result(bexp);
            BOOST_CHECK_EQUAL(result, evaluated_result);

            

            {
                // Nested additions.
                NekMatrix<unsigned int> m1(2,2,buf1);
                NekMatrix<unsigned int> m2(2,2,buf2);
                NekMatrix<unsigned int> m3(m1);
                NekMatrix<unsigned int> m4(m2);

                NekMatrix<unsigned int> m5 = ((m1+m2)+m3)+m4;

                unsigned int buf3[] = { 8, 12, 16, 20 };
                NekMatrix<unsigned int> nested_add_result(2,2,buf3);

                BOOST_CHECK_EQUAL(nested_add_result, m5);
                
                NekMatrix<unsigned int> m6 = m1 + (m2 + (m3 + m4));
                BOOST_CHECK_EQUAL(nested_add_result, m6);
            }
        }
        
        void testNekMatrixMultiplication()
        {
            {
                unsigned int m1_buf[] = {1, 1, 1, 1};
                NekMatrix<unsigned int> m1(2,2,m1_buf);
                NekMatrix<unsigned int> m2(2,2,m1_buf);
                
                unsigned int result_buf[] = {2, 2, 2, 2};
                NekMatrix<unsigned int> result(2,2,result_buf);
                NekMatrix<unsigned int> m3 = m1*m2;
                BOOST_CHECK_EQUAL(m3, result);
            }
            
            double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
            double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
            double m3_buf[] = {77, 66, 54, -5, 99, -61, -50, -12, -18};

            NekMatrix<double> m1(3, 3, m1_buf);
            NekMatrix<double> m2(3, 3, m2_buf);
            NekMatrix<double> m3(3, 3, m3_buf);
            
            NekMatrix<double> result = m1*m2*m3;
            
            double result_buf[] = {-1456238, -1484136, -464512, 1026425, 505353, 583929, 1538925, 1557252, 504714};
            NekMatrix<double> expectedResult(3, 3, result_buf);
            
            double epsilon = 1e-11;
            for(unsigned int i = 0; i < 3; ++i)
            {
                for(unsigned int j = 0; j < 3; ++j)
                {
                    BOOST_CHECK_CLOSE(result(i,j), expectedResult(i,j), epsilon);
                }
            }
            
            NekMatrix<double> result1 = m1*m2*m3*m1*m2*m3*m1*m2*m3;
            double result_buf1[] = {223791291531519928.0, -139146145309301688.0, 241968403742002232.0, -81497861322837100.0, 109613922100149537.0, -116433655760219405.0, -233781330982473300.0, 141216567102193860.0, -250757429804037708.0};
            NekMatrix<double> expectedResult1(3,3,result_buf1);
            for(unsigned int i = 0; i < 3; ++i)
            {
                for(unsigned int j = 0; j < 3; ++j)
                {
                    BOOST_CHECK_CLOSE(result1(i,j), expectedResult1(i,j), epsilon);
                }
            }
        }
        
        void testNekMatrixSomewhatComplicatedExpression()
        {
            {
                double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
                double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
                double m3_buf[] = {77, 66, 54, -5, 99, -61, -50, -12, -18};
    
                NekMatrix<double> m1(3, 3, m1_buf);
                NekMatrix<double> m2(3, 3, m2_buf);
                NekMatrix<double> m3(3, 3, m3_buf);
                
                NekMatrix<double> result = (m1*m2) + m3;
                
                double result_buf[] = {-11157, -5930, 12478, 6755, -522, -10117, 11955, 6150, -12925};
                NekMatrix<double> expectedResult(3,3,result_buf);
                double epsilon = 1e-11;
                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK_CLOSE(result(i,j), expectedResult(i,j), epsilon);
                    }
                }
                
                NekMatrix<double> result1 = m3 + (m1*m2);
                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK_CLOSE(result1(i,j), expectedResult(i,j), epsilon);
                    }
                }
            }
        }
        
        void testNekMatrixComplicatedExpression()
        {
            {
                double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
                double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
                double m3_buf[] = {31, -26, -62, 1, -47, -91, -47, -61, 41};
    
                NekMatrix<double> m1(3, 3, m1_buf);
                NekMatrix<double> m2(3, 3, m2_buf);
                NekMatrix<double> m3(3, 3, m3_buf);
                
                NekMatrix<double> result = (((m3-m1)*m3) * (m1-m2)) + m1*m2*m3;
                
                double result_buf[] = {-1279130, -1162366, 243990, -1663904, 1197403, 2021293, 547959, 1802365, 1422677};
                NekMatrix<double> expectedResult(3,3,result_buf);
                
                double epsilon = 1e-11;
                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK_CLOSE(result(i,j), expectedResult(i,j), epsilon);
                    }
                }
            }
        }
        
        

/*        class TestTemporary : public CountedObject<TestTemporary>
        {
            public:
                TestTemporary() : CountedObject<TestTemporary>() {}
                TestTemporary(const TestTemporary& rhs) : CountedObject<TestTemporary>(rhs) {}
                virtual ~TestTemporary() {}
                TestTemporary& operator=(const TestTemporary& rhs) { CountedObject<TestTemporary>::operator=(rhs); return *this; }
                
                template<typename PolicyType>
                TestTemporary(const expt::Expression<PolicyType>& rhs) : CountedObject<TestTemporary>()
                {
                    //BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<PolicyType>::ResultType, TestTemporary ));
                    rhs.Apply(*this);
                }
                
                template<typename PolicyType>
                TestTemporary& operator=(const expt::Expression<PolicyType>& rhs)
                {
                    //BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<PolicyType>::ResultType, TestTemporary ));
                    rhs.Apply(*this);
                    return *this;
                }
                
                TestTemporary& operator+=(const TestTemporary& rhs)
                {
                    return *this;
                }
                
            private:
        };*/
        
/*        void add(const TestTemporary& lhs, const TestTemporary& rhs, TestTemporary& result)
        {
        }*/
        
        //ENABLE_EXPRESSION_TEMPLATE_OPERATORS(TestTemporary);

        // For these tests I will be testing items of the form A = B op C, where I vary
        // the data types for A, B, and C, as well as op types, to get all possible combinations.
        void testTemporaryGenerationFromSingleLevelBinaryExpressions()
        {
//             TestTemporary a;
//             TestTemporary b;
//             CountedObject<TestTemporary>::clearCounters();
//             TestTemporary c = a + b;
//             
//             CountedObject<TestTemporary>::check(1, 0, 0, 0, 0, 1);
        }
        
        
        class C;
        class B;
        class A
        {
            public:
                explicit A(int a) : val(a) {}
                
                A& operator=(const A& rhs);
                A& operator=(const C& rhs);
                A& operator=(const B& rhs);
                
                template<typename ExpressionPolicyType>
                A(const expt::Expression<ExpressionPolicyType>& rhs)
                {
                    rhs.Apply(*this);
                }

                A& operator+=(const A& rhs)
                {
                    val += rhs.val;
                    return *this;
                }
                
                A& operator-=(const A& rhs)
                {
                    val -= rhs.val;
                    return *this;
                }
                
                A& operator*=(const A& rhs)
                {
                    val *= rhs.val;
                    return *this;
                }
                
                A& operator/=(const A& rhs)
                {
                    val /= rhs.val;
                    return *this;
                }
                
                A& operator+=(const C& rhs);
                A& operator-=(const C& rhs);
                A& operator/=(const C& rhs);
                A& operator*=(const C& rhs);
                
                
                int val;
        };
        
        class B
        {
            public:
                explicit B(int a) : val(a) {}
                
                B& operator=(const B& rhs)
                {
                    val = rhs.val;
                    return *this;
                }
                
                B& operator=(const A& rhs)
                {
                    val = rhs.val;
                    return *this;
                }
                
                template<typename ExpressionPolicyType>
                B(const expt::Expression<ExpressionPolicyType>& rhs)
                {
                    rhs.Apply(*this);
                }
                
                B& operator+=(const A& rhs)
                {
                    val += rhs.val;
                    return *this;
                }
                
                B& operator-=(const A& rhs)
                {
                    val -= rhs.val;
                    return *this;
                }
                
                B& operator*=(const A& rhs)
                {
                    val *= rhs.val;
                    return *this;
                }
                
                B& operator/=(const A& rhs)
                {
                    val /= rhs.val;
                    return *this;
                }
                
                B& operator+=(const B& rhs)
                {
                    val += rhs.val;
                    return *this;
                }
                
                B& operator-=(const B& rhs)
                {
                    val -= rhs.val;
                    return *this;
                }
                
                B& operator*=(const B& rhs)
                {
                    val *= rhs.val;
                    return *this;
                }
                
                B& operator/=(const B& rhs)
                {
                    val /= rhs.val;
                    return *this;
                }
                
                int val;
        };
        
        class C
        {
            public:
                explicit C(int a) : val(a) {}
                
                C& operator=(const C& rhs)
                {
                    val = rhs.val;
                    return *this;
                }
                
                template<typename ExpressionPolicyType>
                C(const expt::Expression<ExpressionPolicyType>& rhs)
                {
                    rhs.Apply(*this);
                }
                
                int val;
        };
        
        A& A::operator+=(const C& rhs)
        {
            val += rhs.val;
            return *this;
        }
        
        A& A::operator-=(const C& rhs)
        {
            val -= rhs.val;
            return *this;
        }
        
        A& A::operator*=(const C& rhs)
        {
            val *= rhs.val;
            return *this;
        }
        
        A& A::operator/=(const C& rhs)
        {
            val /= rhs.val;
            return *this;
        }
        
        A& A::operator=(const A& rhs)
        {
            val = rhs.val;
            return *this;
        }
        
        A& A::operator=(const C& rhs)
        {
            val = rhs.val;
            return *this;
        }
        
        A& A::operator=(const B& rhs)
        {
            val = rhs.val;
            return *this;
        }
                
    }
    
    namespace expt
    {
        using Nektar::UnitTests::A;
        using Nektar::UnitTests::B;
        using Nektar::UnitTests::C;

        
        template<>
        class AdditionTraits<A, A>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Add(A& result, const A& lhs, const A& rhs)
                {
                    result.val = lhs.val + rhs.val;
                }
        };
        
        template<>
        class SubtractionTraits<A, A>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                static void Subtract(A& result, const A& lhs, const A& rhs)
                {
                    result.val = lhs.val - rhs.val;
                }
                
//                 static void SubtractEqual(const A& lhs, const A& rhs, A& result)
//                 {
//                     result.val -= lhs.val;
//                     result.val -= rhs.val;
//                 }
        };
        
        template<>
        class MultiplicationTraits<A, A>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Multiply(A& result, const A& lhs, const A& rhs)
                {
                    result.val = lhs.val * rhs.val;
                }
                
//                 static void MultiplyEqual(const A& lhs, const A& rhs, A& result)
//                 {
//                     result.val *= lhs.val;
//                     result.val *= rhs.val;
//                 }
        };
        
        template<>
        class DivisionTraits<A, A>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                static void Divide(A& result, const A& lhs, const A& rhs)
                {
                    result.val = lhs.val / rhs.val;
                }
                
//                 static void DivideEqual(const A& lhs, const A& rhs, A& result)
//                 {
//                     result.val /= lhs.val;
//                     result.val /= rhs.val;
//                 }
        };
        
        
        
        
        template<>
        class AdditionTraits<A, C>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Add(A& result, const A& lhs, const C& rhs)
                {
                    result.val = lhs.val + rhs.val;
                }
                
//                 static void AddEqual(const A& lhs, const C& rhs, A& result)
//                 {
//                     result.val += lhs.val;
//                     result.val += rhs.val;
//                 }
        };
        
        template<>
        class SubtractionTraits<A, C>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Subtract(A& result, const A& lhs, const C& rhs)
                {
                    result.val = lhs.val - rhs.val;
                }
                
//                 static void SubtractEqual(const A& lhs, const C& rhs, A& result)
//                 {
//                     result.val -= lhs.val;
//                     result.val -= rhs.val;
//                 }
        };
        
        template<>
        class MultiplicationTraits<A, C>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Multiply(A& result, const A& lhs, const C& rhs)
                {
                    result.val = lhs.val * rhs.val;
                }
                
//                 static void MultiplyEqual(const A& lhs, const C& rhs, A& result)
//                 {
//                     result.val *= lhs.val;
//                     result.val *= rhs.val;
//                 }
        };
        
        template<>
        class DivisionTraits<A, C>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Divide(A& result, const A& lhs, const C& rhs)
                {
                    result.val = lhs.val / rhs.val;
                }
                
//                 static void DivideEqual(const A& lhs, const C& rhs, A& result)
//                 {
//                     result.val /= lhs.val;
//                     result.val /= rhs.val;
//                 }
        };
        
        
        
        
    
        template<>
        class AdditionTraits<B, A>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Add(A& result, const B& lhs, const A& rhs)
                {
                    result.val = lhs.val + rhs.val;
                }
                
//                 static void AddEqual(const B& lhs, const A& rhs, A& result)
//                 {
//                     result.val += lhs.val;
//                     result.val += rhs.val;
//                 }
        };
        
        template<>
        class SubtractionTraits<B, A>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Subtract(A& result, const B& lhs, const A& rhs)
                {
                    result.val = lhs.val - rhs.val;
                }
                
//                 static void SubtractEqual(const B& lhs, const A& rhs, A& result)
//                 {
//                     result.val -= lhs.val;
//                     result.val -= rhs.val;
//                 }
        };
        
        template<>
        class MultiplicationTraits<B, A>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Multiply(A& result, const B& lhs, const A& rhs)
                {
                    result.val = lhs.val * rhs.val;
                }
                
//                 static void MultiplyEqual(const B& lhs, const A& rhs, A& result)
//                 {
//                     result.val *= lhs.val;
//                     result.val *= rhs.val;
//                 }
        };
        
        template<>
        class DivisionTraits<B, A>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                static void Divide(A& result, const B& lhs, const A& rhs)
                {
                    result.val = lhs.val / rhs.val;
                }
                
//                 static void DivideEqual(const B& lhs, const A& rhs, A& result)
//                 {
//                     result.val /= lhs.val;
//                     result.val /= rhs.val;
//                 }
        };
    
        template<>
        class AdditionTraits<B, C>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Add(A& result, const B& lhs, const C& rhs)
                {
                    result.val = lhs.val + rhs.val;
                }
                
//                 static void AddEqual(const B& lhs, const C& rhs, A& result)
//                 {
//                     result.val += lhs.val;
//                     result.val += rhs.val;
//                 }
        };
        
        template<>
        class SubtractionTraits<B, C>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Subtract(A& result, const B& lhs, const C& rhs)
                {
                    result.val = lhs.val - rhs.val;
                }
                
//                 static void SubtractEqual(const B& lhs, const C& rhs, A& result)
//                 {
//                     result.val -= lhs.val;
//                     result.val -= rhs.val;
//                 }
        };
        
        template<>
        class MultiplicationTraits<B, C>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Multiply(A& result, const B& lhs, const C& rhs)
                {
                    result.val = lhs.val * rhs.val;
                }
                
//                 static void MultiplyEqual(const B& lhs, const C& rhs, A& result)
//                 {
//                     result.val *= lhs.val;
//                     result.val *= rhs.val;
//                 }
        };
        
        template<>
        class DivisionTraits<B, C>
        {
            public:
                typedef A result_type;
                static const bool HasOpEqual = true;
                
                static void Divide(A& result, const B& lhs, const C& rhs)
                {
                    result.val = lhs.val / rhs.val;
                }
                
//                 static void DivideEqual(const B& lhs, const C& rhs, A& result)
//                 {
//                     result.val /= lhs.val;
//                     result.val /= rhs.val;
//                 }
        };
    }
    
    
    
    namespace UnitTests
    {
//         template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType, 
//         template <typename, typename> class OpType>
//                 class BinaryExpressionPolicy
//                 {
//                     public:
//                         typedef Expression<LhsExpressionPolicyType> LhsExpressionType;
//                         typedef Expression<RhsExpressionPolicyType> RhsExpressionType;
// 
//                         typedef typename LhsExpressionType::ResultType LhsResultType;
//                         typedef typename RhsExpressionType::ResultType RhsResultType;
// 
//                         typedef typename OpType<LhsResultType, RhsResultType>::TraitsType OpTraitsType;
//                 
//                         typedef typename OpType<LhsResultType, RhsResultType>::ResultType ResultType;
//                         typedef typename ExpressionMetadataChooser<OpTraitsType>::MetadataType MetadataType;
//                         typedef typename LhsExpressionType::MetadataType LhsMetadataType;
//                         typedef typename RhsExpressionType::MetadataType RhsMetadataType;
// 
//                         typedef std::pair<LhsExpressionType, RhsExpressionType> DataType;
        
        ENABLE_EXPRESSION_TEMPLATE_OPERATORS(A);
        ENABLE_EXPRESSION_TEMPLATE_OPERATORS2(B, A);
        ENABLE_EXPRESSION_TEMPLATE_OPERATORS2(A, C);
        ENABLE_EXPRESSION_TEMPLATE_OPERATORS2(B, C);

        std::ostream& operator<<(std::ostream& os, const A& rhs)
        {
            os << rhs.val;
        }
        
        std::ostream& operator<<(std::ostream& os, const B& rhs)
        {
            os << rhs.val;
        }
        
        std::ostream& operator<<(std::ostream& os, const C& rhs)
        {
            os << rhs.val;
        }
                
        typedef expt::BinaryExpressionPolicy<expt::ConstantExpressionPolicy<A>, 
                                 expt::ConstantExpressionPolicy<A>, expt::MultiplyOp > FirstType;
        typedef expt::BinaryExpressionPolicy<FirstType, expt::ConstantExpressionPolicy<A>, expt::AddOp> SecondType;
        typedef expt::BinaryExpressionPolicy<FirstType, FirstType, expt::AddOp> ThirdType;
               
//         expt::Expression<SecondType>
//         operator+(const expt::Expression<FirstType>& lhs, const A& rhs)
//         {
//             return lhs + expt::Expression<expt::ConstantExpressionPolicy<A> >(rhs);
//         }
        
//         template<typename LhsPolicyType>
//         typename expt::BinaryExpressionTraits<LhsPolicyType, expt::AddOp, expt::ConstantExpressionPolicy<A> >::ResultType
//         operator+(const expt::Expression<LhsPolicyType>& lhs, const A& rhs)
//         {
//             return lhs + expt::Expression<expt::ConstantExpressionPolicy<A> >(rhs);
//         }

        void testExhaustiveSingleLevelBinaryExpressions()
        {
            // These tests are for single A op B type expressions.
            {
                A lhs(10);
                A rhs(5);
                
                //cout << (lhs + lhs) + ((rhs*rhs)+lhs) << endl;
                //(rhs*rhs)+lhs;
                
                expt::Expression<FirstType> exp = rhs*rhs;
                
                
                expt::Expression<SecondType> exp1 = exp + lhs;
                //expt::Expression<ThirdType> exp3 = exp + exp;
                                 
                A result = lhs + rhs;
                BOOST_CHECK_EQUAL(result.val, 15);
            }
            
            {
                B lhs(10);
                A rhs(5);
                
                A result = lhs + rhs;
                BOOST_CHECK_EQUAL(result.val, 15);
            }
            
            {
                A lhs(10);
                C rhs(5);
                
                A result = lhs + rhs;
                BOOST_CHECK_EQUAL(result.val, 15);
            }
            
            {
                B lhs(10);
                C rhs(5);
                
                A result = lhs + rhs;
                BOOST_CHECK_EQUAL(result.val, 15);
            }
            
            ///////////////////////
            // Subtraction
            ///////////////////////
            {
                A lhs(10);
                A rhs(5);
                
                A result = lhs - rhs;
                BOOST_CHECK_EQUAL(result.val, 5);
            }
            
            {
                B lhs(10);
                A rhs(5);
                
                A result = lhs - rhs;
                BOOST_CHECK_EQUAL(result.val, 5);
            }
            
            {
                A lhs(10);
                C rhs(5);
                
                A result = lhs - rhs;
                BOOST_CHECK_EQUAL(result.val, 5);
            }
            
            {
                B lhs(10);
                C rhs(5);
                
                A result = lhs - rhs;
                BOOST_CHECK_EQUAL(result.val, 5);
            }
            
            ///////////////////////
            // Multiplication
            ///////////////////////
            {
                A lhs(10);
                A rhs(5);
                
                A result = lhs * rhs;
                BOOST_CHECK_EQUAL(result.val, 50);
            }
            
            {
                B lhs(10);
                A rhs(5);
                
                A result = lhs * rhs;
                BOOST_CHECK_EQUAL(result.val, 50);
            }
            
            {
                A lhs(10);
                C rhs(5);
                
                A result = lhs * rhs;
                BOOST_CHECK_EQUAL(result.val, 50);
            }
            
            {
                B lhs(10);
                C rhs(5);
                
                A result = lhs * rhs;
                BOOST_CHECK_EQUAL(result.val, 50);
            }
            
            //////////////////////
            // Division
            //////////////////////
            {
                A lhs(10);
                A rhs(5);
                
                A result = lhs / rhs;
                BOOST_CHECK_EQUAL(result.val, 2);
            }
            
            {
                B lhs(10);
                A rhs(5);
                
                A result = lhs / rhs;
                BOOST_CHECK_EQUAL(result.val, 2);
            }
            
            {
                A lhs(10);
                C rhs(5);
                
                A result = lhs / rhs;
                BOOST_CHECK_EQUAL(result.val, 2);
            }
            
            {
                B lhs(10);
                C rhs(5);
                
                A result = lhs / rhs;
                cout << lhs/rhs << endl;
                BOOST_CHECK_EQUAL(result.val, 2);
            }
        }
        
        void testExhaustive2OpBinaryExpressions()
        {
            // Tests all possible combinations of A op B op C,
            // with varying priorities among the operators.
        }
     }
}

/**
    $Log: testExpressionTemplates.cpp,v $
    Revision 1.11  2006/11/12 17:59:47  bnelson
    *** empty log message ***

    Revision 1.10  2006/11/11 01:32:52  bnelson
    *** empty log message ***

    Revision 1.9  2006/11/08 04:18:22  bnelson
    Added more expression template tests.

    Revision 1.8  2006/11/06 17:10:04  bnelson
    *** empty log message ***

    Revision 1.7  2006/09/30 15:38:29  bnelson
    no message

    Revision 1.6  2006/09/15 02:01:16  bnelson
    no message

    Revision 1.5  2006/09/11 03:28:41  bnelson
    no message

    Revision 1.4  2006/08/28 02:40:51  bnelson
    *** empty log message ***

    Revision 1.3  2006/08/27 02:14:09  bnelson
    Added support for negating an expression.

    Revision 1.2  2006/08/25 01:37:34  bnelson
    no message

    Revision 1.1  2006/08/25 01:36:25  bnelson
    no message


**/
