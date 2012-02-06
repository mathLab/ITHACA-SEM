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
#include <boost/preprocessor/logical/or.hpp>
#include <boost/preprocessor/comparison/greater.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/tuple/eat.hpp>

#include <ExpressionTemplates/Operators.hpp>

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
#include <ExpressionTemplates/Node.hpp>
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

namespace Nektar
{

    #define TEMPLATE() template<
    #define CLOSE_TEMPLATE() > 
    
    #define NEKTAR_TYPENAME() typename 
    #define SHOW_TYPENAME(LhsNumber, RhsNumber) \
        BOOST_PP_IF( BOOST_PP_OR(BOOST_PP_GREATER(LhsNumber,0), BOOST_PP_GREATER(RhsNumber,0)), NEKTAR_TYPENAME, BOOST_PP_EMPTY)()
        
    #define PP_ENUM_TWO_SETS_OF_TYPES(LhsTypeName, LhsNumber, RhsTypeName, RhsNumber) \
        BOOST_PP_IF( BOOST_PP_OR(BOOST_PP_GREATER(LhsNumber,0), BOOST_PP_GREATER(RhsNumber,0)), TEMPLATE, BOOST_PP_EMPTY)() \
        BOOST_PP_ENUM_PARAMS(LhsNumber, typename LhsTypeName) \
        BOOST_PP_COMMA_IF(BOOST_PP_AND(BOOST_PP_GREATER(LhsNumber,0), BOOST_PP_GREATER(RhsNumber,0))) \
        BOOST_PP_ENUM_PARAMS(RhsNumber, typename RhsTypeName) \
        BOOST_PP_IF( BOOST_PP_OR(BOOST_PP_GREATER(LhsNumber,0), BOOST_PP_GREATER(RhsNumber,0)), CLOSE_TEMPLATE, BOOST_PP_EMPTY)() 

    #define GET_TEMPLATED_TYPE(TypeName, TemplateTypeName, Number) TypeName<BOOST_PP_ENUM_PARAMS(Number, TemplateTypeName)>
    
    #define GET_TYPE(TypeName, TemplateTypeName, Number) \
        BOOST_PP_IF(BOOST_PP_GREATER(Number, 0), GET_TEMPLATED_TYPE, TypeName BOOST_PP_TUPLE_EAT(3))(TypeName, TemplateTypeName, Number)
    
    #define DECLARE_MULTIPLICATION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
         expt::Node<expt::Node<GET_TYPE(LeftType, LhsType, NumLeftParams) >, expt::MultiplyOp, expt::Node<GET_TYPE(RightType, RhsType, NumRightParams) > > \
         operator*(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TYPE(RightType, RhsType, NumRightParams)& rhs);
         
     #define DECLARE_DIVISION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
         expt::Node<expt::Node<GET_TYPE(LeftType, LhsType, NumLeftParams) >, expt::DivideOp, expt::Node<GET_TYPE(RightType, RhsType, NumRightParams) > >  \
         operator/(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TYPE(RightType, RhsType, NumRightParams)& rhs);
         
     #define DECLARE_ADDITION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
         expt::Node<expt::Node<GET_TYPE(LeftType, LhsType, NumLeftParams) >, expt::AddOp, expt::Node<GET_TYPE(RightType, RhsType, NumRightParams) > > \
         operator+(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TYPE(RightType, RhsType, NumRightParams)& rhs);
     
     #define DECLARE_SUBTRACTION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
         expt::Node<expt::Node<GET_TYPE(LeftType, LhsType, NumLeftParams) >, expt::SubtractOp, expt::Node<GET_TYPE(RightType, RhsType, NumRightParams) > >  \
         operator-(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TYPE(RightType, RhsType, NumRightParams)& rhs);

 #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
     #define GENERATE_MULTIPLICATION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
         expt::Node<expt::Node<GET_TYPE(LeftType, LhsType, NumLeftParams) >, expt::MultiplyOp, expt::Node<GET_TYPE(RightType, RhsType, NumRightParams) > > \
         operator*(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TYPE(RightType, RhsType, NumRightParams)& rhs) \
         { \
             return expt::CreateBinaryExpression<expt::MultiplyOp>(lhs, rhs); \
         }
         
     #define GENERATE_DIVISION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
         expt::Node<expt::Node<GET_TYPE(LeftType, LhsType, NumLeftParams) >, expt::DivideOp, expt::Node<GET_TYPE(RightType, RhsType, NumRightParams) > >  \
         operator/(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TYPE(RightType, RhsType, NumRightParams)& rhs) \
         { \
             return expt::CreateBinaryExpression<expt::DivideOp>(lhs, rhs); \
         }
         
     #define GENERATE_ADDITION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
         expt::Node<expt::Node<GET_TYPE(LeftType, LhsType, NumLeftParams) >, expt::AddOp, expt::Node<GET_TYPE(RightType, RhsType, NumRightParams) > > \
         operator+(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TYPE(RightType, RhsType, NumRightParams)& rhs) \
         { \
             return expt::CreateBinaryExpression<expt::AddOp>(lhs, rhs); \
         }
     
     #define GENERATE_SUBTRACTION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
         PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
         expt::Node<expt::Node<GET_TYPE(LeftType, LhsType, NumLeftParams) >, expt::SubtractOp, expt::Node<GET_TYPE(RightType, RhsType, NumRightParams) > >  \
         operator-(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
               const GET_TYPE(RightType, RhsType, NumRightParams)& rhs) \
         { \
             return expt::CreateBinaryExpression<expt::SubtractOp>(lhs, rhs); \
         }
         
 #else //NEKTAR_USE_EXPRESSION_TEMPLATES
 
    #define GENERATE_MULTIPLICATION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
        PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
        SHOW_TYPENAME(NumLeftParams, NumRightParams) expt::MultiplyOp::ResultType<GET_TYPE(LeftType, LhsType, NumLeftParams), GET_TYPE(RightType, RhsType, NumRightParams)>::type \
        operator*(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
              const GET_TYPE(RightType, RhsType, NumRightParams)& rhs) \
        { \
            return Multiply(lhs, rhs); \
        }
        
         
    #define GENERATE_DIVISION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
        PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
        SHOW_TYPENAME(NumLeftParams, NumRightParams) expt::DivideOp::ResultType<GET_TYPE(LeftType, LhsType, NumLeftParams), GET_TYPE(RightType, RhsType, NumRightParams)>::type \
        operator/(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
              const GET_TYPE(RightType, RhsType, NumRightParams)& rhs) \
        { \
            return Divide(lhs, rhs); \
        }
      
    #define GENERATE_ADDITION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
        PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
        SHOW_TYPENAME(NumLeftParams, NumRightParams) expt::AddOp::ResultType<GET_TYPE(LeftType, LhsType, NumLeftParams), GET_TYPE(RightType, RhsType, NumRightParams)>::type \
        operator+(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
              const GET_TYPE(RightType, RhsType, NumRightParams)& rhs) \
        { \
            return Add(lhs, rhs); \
        }
        

    #define GENERATE_SUBTRACTION_OPERATOR(LeftType, NumLeftParams, RightType, NumRightParams) \
        PP_ENUM_TWO_SETS_OF_TYPES(LhsType, NumLeftParams, RhsType, NumRightParams) \
        SHOW_TYPENAME(NumLeftParams, NumRightParams) expt::SubtractOp::ResultType<GET_TYPE(LeftType, LhsType, NumLeftParams), GET_TYPE(RightType, RhsType, NumRightParams)>::type \
        operator-(const GET_TYPE(LeftType, LhsType, NumLeftParams)& lhs, \
              const GET_TYPE(RightType, RhsType, NumRightParams)& rhs) \
        { \
            return Subtract(lhs, rhs); \
        }

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_OPERATOR_GENERATORS_HPP

