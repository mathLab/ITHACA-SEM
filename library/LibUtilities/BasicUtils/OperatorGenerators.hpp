///////////////////////////////////////////////////////////////////////////////
//
// File: OperatorGenerators
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
// Description: Generates global operators for classes based upon traits.
//
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_OPERATOR_GENERATORS_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_OPERATOR_GENERATORS_HPP

#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/control/expr_if.hpp>
#include <boost/preprocessor/repetition/enum_trailing_binary_params.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/logical/and.hpp>
#include <boost/preprocessor/comparison/greater.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <LibUtilities/ExpressionTemplates/ArithmeticTraits.hpp>

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

namespace Nektar
{

    
    #define PP_ENUM_TWO_SETS_OF_BINARY_PARAMS(LhsTypeName, LhsVariableName, LhsNumber, RhsTypeName, RhsVariableName, RhsNumber) \
        BOOST_PP_ENUM_BINARY_PARAMS(LhsNumber, LhsTypeName, LhsVariableName) \
        BOOST_PP_COMMA_IF(BOOST_PP_AND(BOOST_PP_GREATER(LhsNumber,0), BOOST_PP_GREATER(RhsNumber,0))) \
        BOOST_PP_ENUM_BINARY_PARAMS(RhsNumber, RhsTypeName, RhsVariableName)
    
    #define PP_ENUM_TWO_SETS_OF_TYPES(LhsTypeName, LhsNumber, RhsTypeName, RhsNumber) \
        BOOST_PP_ENUM_PARAMS(LhsNumber, typename LhsTypeName) \
        BOOST_PP_COMMA_IF(BOOST_PP_AND(BOOST_PP_GREATER(LhsNumber,0), BOOST_PP_GREATER(RhsNumber,0))) \
        BOOST_PP_ENUM_PARAMS(RhsNumber, typename RhsTypeName)
    
    #define GET_TEMPLATED_TYPE(TypeName, TemplateTypeName, Number) TypeName<BOOST_PP_ENUM_PARAMS(Number, TemplateTypeName)>
    
 #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
     #define GENERATE_TEMPLATED_MULTIPLICATION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         template<PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams)> \
         Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams) >, ConstantExpressionPolicy<GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams) >,MultiplyOp > > \
         operator*(const GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)& rhs) \
         { \
             return CreateBinaryExpression<MultiplyOp>(lhs, rhs); \
         }
     #define GENERATE_MULTIPLICATION_OPERATOR(LeftType, RightType) \
         Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<LeftType>, ConstantExpressionPolicy<RightType>, MultiplyOp> > \
         operator*(const LeftType& lhs, const RightType& rhs) \
         { \
                return CreateBinaryExpression<MultiplyOp>(lhs, rhs); \
         } 
         
     #define GENERATE_DIVISION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         template<PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams)> \
         Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams) >, ConstantExpressionPolicy<GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams) >,DivideOp > > \
         operator/(const GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)& rhs) \
         { \
             return CreateBinaryExpression<DivideOp>(lhs, rhs); \
         }
         
     #define GENERATE_TEMPLATED_ADDITION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         template<PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams)> \
         Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams) >, ConstantExpressionPolicy<GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams) >,AddOp > > \
         operator+(const GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)& rhs) \
         { \
             return CreateBinaryExpression<AddOp>(lhs, rhs); \
         }
     
     #define GENERATE_ADDITION_OPERATOR(LeftType, RightType) \
         Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<LeftType>, ConstantExpressionPolicy<RightType>, AddOp> > \
         operator+(const LeftType& lhs, const RightType& rhs) \
         { \
                return CreateBinaryExpression<AddOp>(lhs, rhs); \
         } 
     
     #define GENERATE_SUBTRACTION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         template<PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams)> \
         Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams) >, ConstantExpressionPolicy<GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams) >,SubtractOp > > \
         operator-(const GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)& rhs) \
         { \
             return CreateBinaryExpression<SubtractOp>(lhs, rhs); \
         }
 #else //NEKTAR_USE_EXPRESSION_TEMPLATES
    #define GENERATE_TEMPLATED_MULTIPLICATION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
        template<PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams)> \
        typename MultiplicationTraits<GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams), GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)>::ResultType \
        operator*(const GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
              const GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)& rhs) \
        { \
            return NekMultiply(lhs, rhs); \
        }
        
    #define GENERATE_DIVISION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
        template<PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams)> \
        typename DivisionTraits<GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams), GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)>::ResultType \
        operator/(const GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
              const GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)& rhs) \
        { \
            return NekDivide(lhs, rhs); \
        }
        
    #define GENERATE_TEMPLATED_ADDITION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
        template<PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams)> \
        typename AdditionTraits<GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams), GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)>::ResultType \
        operator+(const GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
              const GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)& rhs) \
        { \
            return NekAdd(lhs, rhs); \
        }
        
    #define GENERATE_SUBTRACTION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
        template<PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams)> \
        typename SubtractionTraits<GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams), GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)>::ResultType \
        operator-(const GET_TEMPLATED_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
              const GET_TEMPLATED_TYPE(RightType, RhsType, NumRightParams)& rhs) \
        { \
            return NekSubtract(lhs, rhs); \
        }
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_OPERATOR_GENERATORS_HPP

