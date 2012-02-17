///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPRESSION_TEMPLATES_INVERSE_OPERATOR_TYPE_TRAITS_HPP
#define EXPRESSION_TEMPLATES_INVERSE_OPERATOR_TYPE_TRAITS_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Operators.hpp>

namespace expt
{    
    /// \brief Traits class indicating which operators have inverses and what they are.
    ///
    /// The ForwardInverseTransform and BackwardInverseTransform modify nodes such as 
    /// (A-B) to (A+(-B)).  To do this, they need to know if an operator is invertible 
    /// and what binary and unary operators are required to perform the inversion.  
    /// This information is specified in this traits class.
    ///
    /// This Traits class specified inverses fo BinaryOperators.  Unary operators 
    /// are handled by UnaryInverseOperatorTypeTraits below.
    template<typename OpType>
    struct InverseOperatorTypeTraits;

    template<>
    struct InverseOperatorTypeTraits<AddOp>
    {
        typedef SubtractOp InverseOpType;
        typedef NegateOp UnaryInversionType;
    };

    template<>
    struct InverseOperatorTypeTraits<SubtractOp>
    {
        typedef AddOp InverseOpType;
        typedef NegateOp UnaryInversionType;
    };

    template<>
    struct InverseOperatorTypeTraits<MultiplyOp>
    {
        typedef DivideOp InverseOpType;
        typedef InvertOp UnaryInversionType;;
    };

    template<>
    struct InverseOperatorTypeTraits<DivideOp>
    {
        typedef MultiplyOp InverseOpType;
        typedef InvertOp UnaryInversionType;;
    };

    template<typename OpType>
    struct UnaryInverseOperatorTypeTraits : public boost::false_type
    {
        typedef void UnaryType;

        template<typename TestType>
        struct Invert
        {
            typedef void Type;
        };
    };

    template<typename OpType>
    struct InverseOperatorHelper
    {
        typedef void Type;
    };

    template<>
    struct InverseOperatorHelper<MultiplyOp>
    {
        typedef DivideOp Type;
    };

    template<>
    struct InverseOperatorHelper<DivideOp>
    {
        typedef MultiplyOp Type;
    };

    template<>
    struct UnaryInverseOperatorTypeTraits<InvertOp> : public boost::true_type
    {
        template<typename OpType>
        struct Invert
        {
            typedef typename InverseOperatorHelper<OpType>::Type Type;
        };
    };

    template<typename OpType>
    struct NegateOperatorHelper
    {
        typedef void Type;
    };

    template<>
    struct NegateOperatorHelper<AddOp>
    {
        typedef SubtractOp Type;
    };

    template<>
    struct NegateOperatorHelper<SubtractOp>
    {
        typedef AddOp Type;
    };

    template<>
    struct UnaryInverseOperatorTypeTraits<NegateOp> : public boost::true_type
    {
        template<typename OpType>
        struct Invert
        {
            typedef typename NegateOperatorHelper<OpType>::Type Type;
        };
    };

    
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_INVERSE_OPERATOR_TYPE_TRAITS_HPP
