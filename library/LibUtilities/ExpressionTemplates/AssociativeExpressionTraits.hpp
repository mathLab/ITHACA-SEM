///////////////////////////////////////////////////////////////////////////////
//
// File: AssociativeExpressionTraits.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_EXPRESSION_TRAITS_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_EXPRESSION_TRAITS_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/AssociativeTraits.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicyFwd.hpp>

namespace Nektar
{
    namespace expt
    {
        // The idea here is to be able to give this object two expression, with an operator type between them, 
        // and have it figure out if it is associative
        //
        // A + (B*C) 
        // (B*C) + A - Don't even need to look since we'll never check in this scenario.
        // (A*B) + (B*C) - Look at the type of (A*B)
        //
        // Examples of change
        // A - (A+B) 
        
        // The class is responsible for
        // 1. Indicating yes or no, the expression combination is associative.
        //    This includes both true associativity and operator changed.
        // 2. Generate the appropriate return type for this expression.
        //

        
        template<typename LhsPolicy, template<typename, typename> class OpType, typename RhsPolicy>
        class AssociativeExpressionTraits
        {
            public:
                static const bool IsStrictlyAssociative = false;
                static const bool IsAssociativeWithOpChange = false;
                static const bool IsAssociative = false;
        };
        
        template<typename LhsPolicy, template<typename, typename> class OpType, 
                 typename RhsLhsPolicy, template<typename, typename> class RhsOpType,
                 typename RhsRhsPolicy>
        class AssociativeExpressionTraits<LhsPolicy, OpType, BinaryExpressionPolicy<RhsLhsPolicy, RhsRhsPolicy, RhsOpType> >
        {
            public:
                typedef AssociativeTraits<typename LhsPolicy::ResultType, OpType,
                                          typename RhsLhsPolicy::ResultType, RhsOpType,
                                          typename RhsRhsPolicy::ResultType> Traits;
                                          
                static const bool IsStrictlyAssociative = Traits::IsStrictlyAssociative;
                static const bool IsAssociativeWithOpChange = Traits::IsAssociativeWithOpChange;
                static const bool IsAssociative = Traits::IsAssociative;
        };
    }
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_EXPRESSION_TRAITS_HPP
