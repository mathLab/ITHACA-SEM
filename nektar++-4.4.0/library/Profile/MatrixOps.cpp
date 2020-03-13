//
//#define NEKTAR_USE_EXPRESSION_TEMPLATES
//
//#include <LibUtilities/ExpressionTemplates/Expression.hpp>
//#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
//#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
//#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
//#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
//#include <LibUtilities/ExpressionTemplates/BinaryExpressionTraits.hpp>
//#include <iostream>
//
//#include "MatrixOpsExprTemp.h"
//#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
//
//namespace Nektar
//{
////    template<typename MatrixType>
////    class BinaryExpressionEvaluator<BinaryExpressionPolicy<ConstantExpressionPolicy<MatrixType>, AddOp, ConstantExpressionPolicy<MatrixType> >,
////                                    BinaryExpressionPolicy<ConstantExpressionPolicy<MatrixType>, AddOp, ConstantExpressionPolicy<MatrixType> >,
////                                    MatrixType, AddOp, BinaryNullOp>
////    {
////        public:
////            typedef BinaryExpressionPolicy<ConstantExpressionPolicy<MatrixType>, AddOp, ConstantExpressionPolicy<MatrixType> > PolicyType;
////            static void Eval(const Expression<PolicyType>& lhs, 
////                             const Expression<PolicyType>& rhs,
////                             Accumulator<MatrixType>& result)
////            {
////                const MatrixType& m1 = (*lhs).first;
////                const MatrixType& m2 = (*lhs).second;
////                const MatrixType& m3 = (*rhs).first;
////                const MatrixType& m4 = (*rhs).second;
////                
////                const typename MatrixType::NumberType* b1 = m1.GetRawPtr();
////                const typename MatrixType::NumberType* b2 = m2.GetRawPtr();
////                const typename MatrixType::NumberType* b3 = m3.GetRawPtr();
////                const typename MatrixType::NumberType* b4 = m4.GetRawPtr();
////                
////                unsigned int size = m1.GetStorageSize();
////
////                typename MatrixType::NumberType* r = (*result).GetRawPtr();
////                typename MatrixType::NumberType* end = r+size;
////                while( r < end )
////                {
////                    *r = *b1 + *b2 + *b3 + *b4;
////                    ++b1;
////                    ++b2;
////                    ++b3;
////                    ++b4;
////                    ++r;
////                } 
////            }
////    };
//}
//
//void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
//                const Nektar::NekMatrix<double>& v1,
//                const Nektar::NekMatrix<double>& v2)
//{
//    result = v1 + v2;
//}
//
//void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
//                const Nektar::NekMatrix<double>& v1,
//                const Nektar::NekMatrix<double>& v2,
//                const Nektar::NekMatrix<double>& v3)
//{
//    result = v1 + v2 + v3;
//}
//
//void AddMatricesExprTemp(Nektar::NekMatrix<double>& result,
//                const Nektar::NekMatrix<double>& v1,
//                const Nektar::NekMatrix<double>& v2,
//                const Nektar::NekMatrix<double>& v3,
//                const Nektar::NekMatrix<double>& v4)
//{
//    result = v1 + v2 + v3 + v4;
//}
//
//void AddMatricesExprTempResultAlloc(
//                const Nektar::NekMatrix<double>& v1,
//                const Nektar::NekMatrix<double>& v2)
//{
//    Nektar::NekMatrix<double> result = v1 + v2;
//}
//
//void AddMatricesExprTempResultAlloc(
//                const Nektar::NekMatrix<double>& v1,
//                const Nektar::NekMatrix<double>& v2,
//                const Nektar::NekMatrix<double>& v3)
//{
//    Nektar::NekMatrix<double> result = v1 + v2 + v3;
//}
//
//void AddMatricesExprTempResultAlloc(
//                const Nektar::NekMatrix<double>& v1,
//                const Nektar::NekMatrix<double>& v2,
//                const Nektar::NekMatrix<double>& v3,
//                const Nektar::NekMatrix<double>& v4)
//{
//    Nektar::NekMatrix<double> result = v1 + v2 + v3 + v4;
//}
//
//
//
//
//
//
//
////#include "IntWrapper.h"
////#include "IntWrapperExprTemp.h"
////#include <boost/lexical_cast.hpp>
////#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
////#include <Profile/StringConcat.h>
////#include <Profile/StringConcatExprTemp.h>
////#include <boost/timer.hpp>
////#include <boost/progress.hpp>
////#include "VectorOps.h"
////#include "VectorOpsExprTemp.h"
////#include "MatrixOps.h"
////#include "MatrixOpsExprTemp.h"
////
////
////#include <iostream>
////
////
////using namespace std;
////
////void Run2MatrixTests(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekMatrix<double> strVals[] = { Nektar::NekMatrix<double>(stringSize, stringSize), 
////        Nektar::NekMatrix<double>(stringSize, stringSize)};
////        
////    Nektar::NekMatrix<double> result(stringSize, stringSize);
////    
////    cout << "2 Matrix Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddMatrices(result, strVals[0], strVals[1]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesAccum(result, strVals[0], strVals[1]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesExprTemp(result, strVals[0], strVals[1]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesHandCoded(result, strVals[0], strVals[1]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Hand Coded - Total : " << elapsedTime << endl;
////    cout << "Hand Coded - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////
////}
////
////
////void Run3MatrixTests(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekMatrix<double> strVals[] = { Nektar::NekMatrix<double>(stringSize, stringSize), 
////        Nektar::NekMatrix<double>(stringSize, stringSize),
////        Nektar::NekMatrix<double>(stringSize, stringSize)};
////        
////    Nektar::NekMatrix<double> result(stringSize, stringSize);
////    
////    cout << "3 Matrix Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddMatrices(result, strVals[0], strVals[1], strVals[2]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesAccum(result, strVals[0], strVals[1], strVals[2]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesExprTemp(result, strVals[0], strVals[1], strVals[2]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesHandCoded(result, strVals[0], strVals[1], strVals[2]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Hand Coded - Total : " << elapsedTime << endl;
////    cout << "Hand Coded - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////
////}
////
////
////void Run4MatrixTests(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekMatrix<double> strVals[] = { Nektar::NekMatrix<double>(stringSize, stringSize), 
////        Nektar::NekMatrix<double>(stringSize, stringSize),
////        Nektar::NekMatrix<double>(stringSize, stringSize),
////        Nektar::NekMatrix<double>(stringSize, stringSize)};
////        
////    Nektar::NekMatrix<double> result(stringSize, stringSize);
////    
////    cout << "4 Matrix Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddMatrices(result, strVals[0], strVals[1], strVals[2], strVals[3]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesAccum(result, strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesExprTemp(result, strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesHandCoded(result, strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Hand Coded - Total : " << elapsedTime << endl;
////    cout << "Hand Coded - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////
////}
////
////
////void Run2MatrixTestsResultAlloc(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekMatrix<double> strVals[] = { Nektar::NekMatrix<double>(stringSize, stringSize), 
////        Nektar::NekMatrix<double>(stringSize, stringSize)};
////        
////    Nektar::NekMatrix<double> result(stringSize, stringSize);
////    
////    cout << "2 Matrix Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddMatricesResultAlloc(strVals[0], strVals[1]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesAccumResultAlloc(strVals[0], strVals[1]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesExprTempResultAlloc(strVals[0], strVals[1]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesHandCodedResultAlloc(strVals[0], strVals[1]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Hand Coded - Total : " << elapsedTime << endl;
////    cout << "Hand Coded - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////
////}
////
////
////void Run3MatrixTestsResultAlloc(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekMatrix<double> strVals[] = { Nektar::NekMatrix<double>(stringSize, stringSize), 
////        Nektar::NekMatrix<double>(stringSize, stringSize),
////        Nektar::NekMatrix<double>(stringSize, stringSize)};
////        
////    Nektar::NekMatrix<double> result(stringSize, stringSize);
////    
////    cout << "3 Matrix Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddMatricesResultAlloc(strVals[0], strVals[1], strVals[2]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesAccumResultAlloc(strVals[0], strVals[1], strVals[2]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesExprTempResultAlloc(strVals[0], strVals[1], strVals[2]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesHandCodedResultAlloc(strVals[0], strVals[1], strVals[2]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Hand Coded - Total : " << elapsedTime << endl;
////    cout << "Hand Coded - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////
////}
////
////
////void Run4MatrixTestsResultAlloc(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekMatrix<double> strVals[] = { Nektar::NekMatrix<double>(stringSize, stringSize), 
////        Nektar::NekMatrix<double>(stringSize, stringSize),
////        Nektar::NekMatrix<double>(stringSize, stringSize),
////        Nektar::NekMatrix<double>(stringSize, stringSize)};
////        
////    Nektar::NekMatrix<double> result(stringSize, stringSize);
////    
////    cout << "4 Matrix Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddMatricesResultAlloc(strVals[0], strVals[1], strVals[2], strVals[3]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesAccumResultAlloc(strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesExprTempResultAlloc(strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddMatricesHandCodedResultAlloc(strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Hand Coded - Total : " << elapsedTime << endl;
////    cout << "Hand Coded - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result(0,0) << endl;
////
////
////}
////
////
////void Run2VectorTests(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekVector<double> strVals[] = { Nektar::NekVector<double>(stringSize), 
////        Nektar::NekVector<double>(stringSize),
////        Nektar::NekVector<double>(stringSize),
////        Nektar::NekVector<double>(stringSize)};
////        
////    Nektar::NekVector<double> result(stringSize);
////    
////    cout << "2 Vector Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddVectors(result, strVals[0], strVals[1]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[0] << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddVectorsAccum(result, strVals[0], strVals[1]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[1] << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddVectorsExprTemp(result, strVals[0], strVals[1]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[2] << endl;
////    
////}
////
////void Run3VectorTests(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekVector<double> strVals[] = { Nektar::NekVector<double>(stringSize), 
////        Nektar::NekVector<double>(stringSize),
////        Nektar::NekVector<double>(stringSize),
////        Nektar::NekVector<double>(stringSize)};
////        
////    Nektar::NekVector<double> result(stringSize);
////    
////    cout << "3 Vector Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddVectors(result, strVals[0], strVals[1], strVals[2]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[0] << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddVectorsAccum(result, strVals[0], strVals[1], strVals[2]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[1] << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddVectorsExprTemp(result, strVals[0], strVals[1], strVals[2]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[2] << endl;
////    
////}
////
////void Run4VectorTests(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekVector<double> strVals[] = { Nektar::NekVector<double>(stringSize), 
////        Nektar::NekVector<double>(stringSize),
////        Nektar::NekVector<double>(stringSize),
////        Nektar::NekVector<double>(stringSize)};
////        
////    Nektar::NekVector<double> result(stringSize);
////    
////    cout << "4 Vector Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddVectors(result, strVals[0], strVals[1], strVals[2], strVals[3]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[0] << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddVectorsAccum(result, strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[1] << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddVectorsExprTemp(result, strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    cout << result[2] << endl;
////    
////}
////
////void Run2StringTests(unsigned int numTests, int stringSize)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    std::string strVals[] = { std::string(stringSize, 'a'), std::string(stringSize, 'a'),
////                              std::string(stringSize, 'a'), std::string(stringSize, 'a')};
////    std::string result;
////    
////    cout << "2 String Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddStrings(result, strVals[0], strVals[1], strVals[2], strVals[3]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddStringsAccum(result, strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////    cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////        AddStringsExprTemp(result, strVals[0], strVals[1], strVals[2], strVals[3]);
////    }
////    elapsedTime = t.elapsed();
////    cout << "Expression Templates - Total : " << elapsedTime << endl;
////    cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    
////    
////}
////
////void RunIntWrapperTests(unsigned int numTests)
////{
////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    {
////        norm::IntWrapper strVals[] = { norm::IntWrapper(1), norm::IntWrapper(2), norm::IntWrapper(3), norm::IntWrapper(4),
////                                 norm::IntWrapper(5), norm::IntWrapper(6) };
////        norm::IntWrapper result;
////        
////        cout << "Int Wrapper Tests" << endl;
////        cout << "---------------------" << endl;
////        t.restart();
////        for(unsigned int i = 0; i < numTests; ++i)
////        {
////           AddIntWrapper(result, strVals[0], strVals[1]);    
////        }
////        elapsedTime = t.elapsed();
////        cout << "Straight Addition - Total : " << elapsedTime << endl;
////        cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;
////        
////        t.restart();
////        for(unsigned int i = 0; i < numTests; ++i)
////        {
////            AddIntWrapperManualAccum(result, strVals[0], strVals[1]);
////        }
////        elapsedTime = t.elapsed();
////        cout << "Manual Accumulator - Total : " << elapsedTime << endl;
////        cout << "Manual Accumulator - PerOp : " << elapsedTime/(double)numTests << endl;
////    }
////    
////    {
////        expr::IntWrapper strVals[] = { expr::IntWrapper(1), expr::IntWrapper(2), expr::IntWrapper(3), expr::IntWrapper(4),
////                                 expr::IntWrapper(5), expr::IntWrapper(6) };
////        expr::IntWrapper result;
////        t.restart();
////        for(unsigned int i = 0; i < numTests; ++i)
////        {
////            AddIntWrapperExprTemp(result, strVals[0], strVals[1]);
////        }
////        elapsedTime = t.elapsed();
////        cout << "Expression Templates - Total : " << elapsedTime << endl;
////        cout << "Expression Templates - PerOp : " << elapsedTime/(double)numTests << endl;
////    }
////    
////    
////}
////
////int main(int argc, char** argv)
////{
////    if( argc != 3 )
////    {
////        std::cerr << "Usage: perf <ProblemSize> <NumTests>" << std::endl;
////        return 1;
////    }
////    
////    unsigned int n = boost::lexical_cast<unsigned int>(argv[1]);
////    unsigned int numTests = boost::lexical_cast<unsigned int>(argv[2]);
////    
////    //Run2StringTests(numTests, n);
////    //RunIntWrapperTests(numTests);
////    //Run4MatrixTestsResultAlloc(numTests, n);
////    Run4MatrixTests(numTests, n);
////    //Run2VectorTests(numTests, n);
////    //Run4VectorTests(numTests, n);
////    return 0;
////}
