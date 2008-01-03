///////////////////////////////////////////////////////////////////////////////
//
// File: BinaryExpressionPolicy.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_POLICY_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_POLICY_HPP

#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicyFwd.hpp>
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/Accumulator.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionEvaluatorFwd.hpp>
#include <algorithm>

namespace Nektar
{
    template<typename LhsPolicy, typename RhsPolicy, template <typename, typename> class OpType>
    class BinaryExpressionPolicy
    {
        public:
            typedef Expression<LhsPolicy> LhsExpressionType;
            typedef Expression<RhsPolicy> RhsExpressionType;

            typedef typename LhsPolicy::ResultType LhsResultType;
            typedef typename RhsPolicy::ResultType RhsResultType;

            typedef typename OpType<LhsResultType, RhsResultType>::TraitsType OpTraitsType;
            
            typedef typename OpType<LhsResultType, RhsResultType>::ResultType ResultType;
            //typedef typename ExpressionMetadataChooser<OpTraitsType>::MetadataType MetadataType;
            typedef typename BinaryExpressionMetadataTraits<LhsResultType, RhsResultType, OpType>::MetadataType MetadataType;
            typedef typename LhsPolicy::MetadataType LhsMetadataType;
            typedef typename RhsPolicy::MetadataType RhsMetadataType;

            typedef std::pair<LhsExpressionType, RhsExpressionType> DataType;

        public:
        
            static void Evaluate(Accumulator<ResultType>& accum, const DataType& d)
            {
                Evaluate<BinaryNullOp>(accum, d);
            }
                
            template<template <typename, typename> class ParentOpType>
            static void Evaluate(Accumulator<ResultType>& accum, const DataType& d)
            {
                // Evaluation of a binary expression is separated into another class because
                // there are a wide variety of diffrent actions which must be taken based upon the 
                // policy types, operator type, and incoming parent operation type.
                BinaryExpressionEvaluator<LhsPolicy, RhsPolicy, 
                                         ResultType, OpType, 
                                         ParentOpType>::Eval(d.first, d.second, accum);
            }
            
            template<typename CheckType>
            static bool ContainsReference(const CheckType& result, const DataType& m_data)
            {
                return LhsPolicy::ContainsReference(result, *m_data.first) ||
                       RhsPolicy::ContainsReference(result, *m_data.second);
            }
            
            static void InitializeMetadata(const DataType& data, MetadataType& m)
            {
                m = MetadataType(data.first.GetMetadata(), data.second.GetMetadata());
            }
            
            static void Print(std::ostream& os, const DataType& data)
            {
                os << "(" << data.first << OpType<LhsResultType, RhsResultType>::AsString() << data.second << ")";
            }
                
    };
}

#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_POLICY_HPP
