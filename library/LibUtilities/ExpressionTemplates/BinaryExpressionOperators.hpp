///////////////////////////////////////////////////////////////////////////////
//
// File: BinaryExpressionOperators.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_OPERATORS_HPP
#define NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_OPERATORS_HPP

#include <LibUtilities/ExpressionTemplates/BinaryExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>

#include <boost/call_traits.hpp>

namespace Nektar
{
    namespace expt
    {
        template<typename LhsType, typename RhsType>
        class BinaryOp
        {
            public:
                
        };

        /// \brief An expression to negate an object of ParameterType.
        /// Parameter type is the actual object, not the expression that may lead to it.
        template<typename LhsType, typename RhsType>
        class AddOp : public BinaryOp<LhsType, RhsType>
        {
            public:
                typedef typename BinaryExpressionTraits<LhsType, RhsType>::AdditionResultType ResultType;

                typedef typename ExpressionMetadataChooser<LhsType>::MetadataType LhsMetadataType;
                typedef typename ExpressionMetadataChooser<RhsType>::MetadataType RhsMetadataType;
                typedef typename ExpressionMetadataChooser<ResultType>::MetadataType ResultMetadataType;

                // This method is only valid if LhsType and ResultType are the same type.  Enforce it?
                static void ApplyEqual(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<RhsType>::const_reference rhs)
                {
                    result += rhs;
                }

                static void Apply(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<LhsType>::const_reference lhs, typename boost::call_traits<RhsType>::const_reference rhs)
                {
                    add(lhs, rhs, result);
                }

                static ResultMetadataType CreateBinaryMetadata(const LhsMetadataType& lhs, const RhsMetadataType& rhs)
                {
                    return ResultMetadataType::CreateForAddition(lhs, rhs);
                }

            private:

        };

        template<typename LhsType, typename RhsType>
        class MultiplyOp : public BinaryOp<LhsType, RhsType>
        {
            public:
                typedef typename BinaryExpressionTraits<LhsType, RhsType>::MultiplicationResultType ResultType;

                typedef typename ExpressionMetadataChooser<LhsType>::MetadataType LhsMetadataType;
                typedef typename ExpressionMetadataChooser<RhsType>::MetadataType RhsMetadataType;
                typedef typename ExpressionMetadataChooser<ResultType>::MetadataType ResultMetadataType;

                // This method is only valid if LhsType and ResultType are the same type.  Enforce it?
                static void ApplyEqual(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<RhsType>::const_reference rhs)
                {
                    result *= rhs;
                }

                static void ApplyLhsEqual(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<LhsType>::const_reference lhs)
                {
                    result.TimesEqualLhs(lhs);
                }
                
                static void Apply(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<LhsType>::const_reference lhs, typename boost::call_traits<RhsType>::const_reference rhs)
                {
                    multiply(lhs, rhs, result);
                }

                static ResultMetadataType CreateBinaryMetadata(const LhsMetadataType& lhs, const RhsMetadataType& rhs)
                {
                    return ResultMetadataType::CreateForMultiplication(lhs, rhs);
                }

                
            private:

        };

        template<typename LhsType, typename RhsType>
        class DivideOp : public BinaryOp<LhsType, RhsType>
        {
            public:
                typedef typename BinaryExpressionTraits<LhsType, RhsType>::DivideResultType ResultType;

                typedef typename ExpressionMetadataChooser<LhsType>::MetadataType LhsMetadataType;
                typedef typename ExpressionMetadataChooser<RhsType>::MetadataType RhsMetadataType;
                typedef typename ExpressionMetadataChooser<ResultType>::MetadataType ResultMetadataType;

                // This method is only valid if LhsType and ResultType are the same type.  Enforce it?
                static void ApplyEqual(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<RhsType>::const_reference rhs)
                {
                    result /= rhs;
                }

                static void Apply(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<LhsType>::const_reference lhs, typename boost::call_traits<RhsType>::const_reference rhs)
                {
                    divide(lhs, rhs, result);
                }

                static ResultMetadataType CreateBinaryMetadata(const LhsMetadataType& lhs, const RhsMetadataType& rhs)
                {
                    return ResultMetadataType::CreateForDivision(lhs, rhs);
                }

            private:

        };

        template<typename LhsType, typename RhsType>
        class SubtractOp : public BinaryOp<LhsType, RhsType>
        {
            public:
                typedef typename BinaryExpressionTraits<LhsType, RhsType>::MultiplyResultType ResultType;

                typedef typename ExpressionMetadataChooser<LhsType>::MetadataType LhsMetadataType;
                typedef typename ExpressionMetadataChooser<RhsType>::MetadataType RhsMetadataType;
                typedef typename ExpressionMetadataChooser<ResultType>::MetadataType ResultMetadataType;

                // This method is only valid if LhsType and ResultType are the same type.  Enforce it?
                static void ApplyEqual(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<RhsType>::const_reference rhs)
                {
                    result -= rhs;
                }

                static void Apply(typename boost::call_traits<ResultType>::reference result, typename boost::call_traits<LhsType>::const_reference lhs, typename boost::call_traits<RhsType>::const_reference rhs)
                {
                    subtract(lhs, rhs, result);
                }

                static ResultMetadataType CreateBinaryMetadata(const LhsMetadataType& lhs, const RhsMetadataType& rhs)
                {
                    return ResultMetadataType::CreateForSubtraction(lhs, rhs);
                }

            private:

        };
    }
}

#endif // NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_OPERATORS_HPP

/**
    $Log: BinaryExpressionOperators.hpp,v $
    Revision 1.3  2006/09/21 01:07:59  bnelson
    no message

    Revision 1.2  2006/09/14 02:08:59  bnelson
    Fixed gcc compile errors

    Revision 1.1  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

**/
