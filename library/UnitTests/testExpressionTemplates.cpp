/////////////////////////////////////////////////////////////////////////////////
////
//// File: testExpressionTemplates.cpp
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
//#include <UnitTests/testExpressionTemplates.h>
//#include <UnitTests/CountedObject.h>
//
//#include <LibUtilities/LinearAlgebra/NekPoint.hpp>
//#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
//
//
//#include <boost/test/auto_unit_test.hpp>
//#include <boost/test/test_case_template.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <boost/test/unit_test.hpp>
//#include <boost/progress.hpp>
//
//#include <iostream>
//
//using std::cout;
//using std::endl;
//
//namespace Nektar
//{
//    namespace UnitTests
//    {
//        #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//        using namespace Nektar;
//
//        BOOST_AUTO_TEST_CASE(testConstantExpressions)
//        {
//             using namespace Nektar;
// 
//             {
//                 int t = 7;
//                 Expression<ConstantExpressionPolicy<int> > e1(t);
//                 BOOST_CHECK_EQUAL(*e1, t);
// 
//                 Expression<ConstantExpressionPolicy<int> > e2(e1);
//                 BOOST_CHECK_EQUAL(*e2, *e1);
//             }
// 
//             typedef NekPoint<unsigned int, ThreeD> Point;
//             {
//                 Point p(1,2,3);
//                 Expression<ConstantExpressionPolicy<Point> > e1(p);
// 
//                 Point p1(e1);
//                 BOOST_CHECK_EQUAL(p, p1);
// 
//                 Point p2 = e1;
//                 BOOST_CHECK_EQUAL(p, p2);
// 
//                 Point p3(9, 10, 11);
//                 BOOST_CHECK(p != p3);
//                 p3 = e1;
//                 BOOST_CHECK_EQUAL(p, p3);
//             }
// 
//             {
//                 // TODO - Find a way to prevent temporaries (meaning that the parameter to
//                 // this call is temporary and that could cause problems).
//                 //Expression<ConstantExpressionPolicy<Point> > e2(Point(1,2,3));
//             }
//        }
//
//        BOOST_AUTO_TEST_CASE(testUnaryExpressions)
//        {
////             using namespace Nektar;
//// 
////             NekPoint<double, 3> p(1,2,3);
////             NekPoint<double, 3> p1(-(-p));
//// 
////             BOOST_CHECK_EQUAL(p, p1);
//// 
////             NekMatrix<double, FullMatrixTag> m(3,3);
////             NekMatrix<double, FullMatrixTag> m1(-m);
//        }
//
//        BOOST_AUTO_TEST_CASE(testNekMatrixMetadata)
//        {
////             using namespace Nektar;
//// 
////             // Constant
////             NekMatrix<double> m(3,3);
////             ConstantExpressionTraits<NekMatrix<double> >::MetadataType t(m);
////             expt::ExpressionMetadataChooser<ConstantExpressionTraits<NekMatrix<double> > >::MetadataType t10(m);
////             ConstantExpressionPolicy<NekMatrix<double, FullMatrixTag, StandardMatrixTag, 0, void> >::MetadataType t1(m);
////             Expression<ConstantExpressionPolicy<NekMatrix<double> > >::MetadataType t2(m);
////             Expression<ConstantExpressionPolicy<NekMatrix<double> > > m_exp(m);
////             unsigned int i = t.Rows;
////             unsigned int j = t1.Rows;
////             unsigned int k = t2.Rows;
////             unsigned int l = t10.Rows;
////             BOOST_CHECK_EQUAL(m.GetRows(), m_exp.GetMetadata().Rows);
////             BOOST_CHECK_EQUAL(m.GetColumns(), m_exp.GetMetadata().Columns);
//
//            // Unary
////             NekMatrix<double> m1(3,3);
////             Expression<UnaryExpressionPolicy<Expression<ConstantExpressionPolicy<NekMatrix<double> > >, expt::NegateOp> > m1_exp = 
////                 Expression<UnaryExpressionPolicy<Expression<ConstantExpressionPolicy<NekMatrix<double> > >, expt::NegateOp> >(
////                 Expression<ConstantExpressionPolicy<NekMatrix<double> > >(m1));
////             BOOST_CHECK_EQUAL(m1.GetRows(), m1_exp.GetMetadata().Rows);
////             BOOST_CHECK_EQUAL(m1.GetColumns(), m1_exp.GetMetadata().Columns);
//
//            // Binary
//            // TODO
//         }
// 
//         BOOST_AUTO_TEST_CASE(testBinaryExpressions)
//         {
//             using namespace Nektar;
// 
//             unsigned int b1[] = {1, 2, 3, 4};
//             unsigned int b2[] = {3, 4, 5, 6};
//             Array<OneD, unsigned int> buf1(4, b1);
//             Array<OneD, unsigned int> buf2(4, b2);
// // 
// //             // The following crashes because of the temporary.  What to do here?
// //             //Nektar::ConstantExpression<Nektar::LibUtilities::NekMatrix<unsigned int> > m1(Nektar::LibUtilities::NekMatrix<unsigned int>(2, 2, buf1));
// //             //Nektar::ConstantExpression<Nektar::LibUtilities::NekMatrix<unsigned int> > m2(Nektar::LibUtilities::NekMatrix<unsigned int>(2, 2, buf2));
// // 
//             NekMatrix<unsigned int> lhs(2,2,buf1);
//             NekMatrix<unsigned int> rhs(2,2,buf2);
//             Expression<ConstantExpressionPolicy<Nektar::NekMatrix<unsigned int> > > m1(lhs);
//             Expression<ConstantExpressionPolicy<Nektar::NekMatrix<unsigned int> > > m2(rhs);
// 
// 
//             typedef ConstantExpressionPolicy<NekMatrix<unsigned int> > LhsPolicy;
//             typedef ConstantExpressionPolicy<NekMatrix<unsigned int> > RhsPolicy;
//             typedef BinaryExpressionPolicy<LhsPolicy, expt::AddOp, RhsPolicy> BinaryPolicy;
//             
//             BinaryPolicy::DataType d(m1, m2);
//             Expression<BinaryPolicy> bexp(d);
// 
//             unsigned int result_buf[] = {4, 6, 8, 10};
//             NekMatrix<unsigned int> result(2, 2, result_buf);
//             NekMatrix<unsigned int> evaluated_result(bexp);
//             BOOST_CHECK_EQUAL(result, evaluated_result);
// 
//             
// 
//             {
//                 // Nested additions.
//                 NekMatrix<unsigned int> m1(2,2,buf1);
//                 NekMatrix<unsigned int> m2(2,2,buf2);
//                 NekMatrix<unsigned int> m3(m1);
//                 NekMatrix<unsigned int> m4(m2);
// 
//                 NekMatrix<unsigned int> m5 = ((m1+m2)+m3)+m4;
// 
//                 unsigned int buf3[] = { 8, 12, 16, 20 };
//                 NekMatrix<unsigned int> nested_add_result(2,2,buf3);
// 
//                 BOOST_CHECK_EQUAL(nested_add_result, m5);
//                 
//                 NekMatrix<unsigned int> m6 = m1 + (m2 + (m3 + m4));
//                 BOOST_CHECK_EQUAL(nested_add_result, m6);
//             }
//        }
//        
//        BOOST_AUTO_TEST_CASE(testNekMatrixMultiplication)
//        {
//             {
//                 unsigned int m1_buf[] = {1, 1, 1, 1};
//                 NekMatrix<unsigned int> m1(2,2,m1_buf);
//                 NekMatrix<unsigned int> m2(2,2,m1_buf);
//                 
//                 unsigned int result_buf[] = {2, 2, 2, 2};
//                 NekMatrix<unsigned int> result(2,2,result_buf);
//                 NekMatrix<unsigned int> m3 = m1*m2;
//                 BOOST_CHECK_EQUAL(m3, result);
//             }
//             
//             double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
//             double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
//             double m3_buf[] = {77, 66, 54, -5, 99, -61, -50, -12, -18};
// 
//             NekMatrix<double> m1(3, 3, m1_buf);
//             NekMatrix<double> m2(3, 3, m2_buf);
//             NekMatrix<double> m3(3, 3, m3_buf);
//             
//             NekMatrix<double> result = m1*m2*m3;
//             
//             double result_buf[] = { -2411761, -889268, -851464, -150650, -565220, -381906, 987268, 242006, 275320 };
//             NekMatrix<double> expectedResult(3, 3, result_buf);
//             
//             double epsilon = 1e-11;
//             for(unsigned int i = 0; i < 3; ++i)
//             {
//                 for(unsigned int j = 0; j < 3; ++j)
//                 {
//                     BOOST_CHECK_CLOSE(*result(i,j), *expectedResult(i,j), epsilon);
//                 }
//             }
//        }
//        
//        BOOST_AUTO_TEST_CASE(testNekMatrixSomewhatComplicatedExpression)
//        {
//             {
//                 double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
//                 double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
//                 double m3_buf[] = {77, 66, 54, -5, 99, -61, -50, -12, -18};
//     
//                 NekMatrix<double> m1(3, 3, m1_buf);
//                 NekMatrix<double> m2(3, 3, m2_buf);
//                 NekMatrix<double> m3(3, 3, m3_buf);
//                 
//                 NekMatrix<double> result = (m1*m2) + m3;
//                 
//                 double result_buf[] = { -11934, -1174, -2318, -10897, -8360, -6683, -14273, -4373, -4310 };
//                 NekMatrix<double> expectedResult(3,3,result_buf);
//                 double epsilon = 1e-11;
//                 for(unsigned int i = 0; i < 3; ++i)
//                 {
//                     for(unsigned int j = 0; j < 3; ++j)
//                     {
//                         BOOST_CHECK_CLOSE(*result(i,j), *expectedResult(i,j), epsilon);
//                     }
//                 }
//                 
//                 NekMatrix<double> result1 = m3 + (m1*m2);
//                 for(unsigned int i = 0; i < 3; ++i)
//                 {
//                     for(unsigned int j = 0; j < 3; ++j)
//                     {
//                         BOOST_CHECK_CLOSE(*result1(i,j), *expectedResult(i,j), epsilon);
//                     }
//                 }
//             }
//        }
//        
//        BOOST_AUTO_TEST_CASE(testNekMatrixComplicatedExpression)
//        {
//             {
//                 double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
//                 double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
//                 double m3_buf[] = {31, -26, -62, 1, -47, -91, -47, -61, 41};
//     
//                 NekMatrix<double> m1(3, 3, m1_buf);
//                 NekMatrix<double> m2(3, 3, m2_buf);
//                 NekMatrix<double> m3(3, 3, m3_buf);
//                 
//                 NekMatrix<double> result = (((m3-m1)*m3) * (m1-m2)) + m1*m2*m3;
//                 double result_buf[] = { -2146767, -3204880, -749636, 162494, 2048276, 2549188, -784134, 758853, 1439441 };
//
//                 NekMatrix<double> expectedResult(3,3,result_buf);
//                 
//                 double epsilon = 1e-11;
//                 for(unsigned int i = 0; i < 3; ++i)
//                 {
//                     for(unsigned int j = 0; j < 3; ++j)
//                     {
//                         BOOST_CHECK_CLOSE(*result(i,j), *expectedResult(i,j), epsilon);
//                     }
//                 }
//             }
//        }
//        
//        #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
//     }
//}
//
