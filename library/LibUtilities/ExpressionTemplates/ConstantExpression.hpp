///////////////////////////////////////////////////////////////////////////////
//
// File: ConstantExpression.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_CONSTANT_EXPRESSION_HPP
#define NEKTAR_LIB_UTILITIES_CONSTANT_EXPRESSION_HPP
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/call_traits.hpp>

#include <LibUtilities/ExpressionTemplates/ExpressionReturnType.hpp>
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
#include <LibUtilities/ExpressionTemplates/Accumulator.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/NullOp.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Nektar
{
    template<typename Type>
    class ConstantExpressionPolicy
    {
        public:
            typedef Type ResultType;
            typedef typename boost::call_traits<Type>::const_reference DataType;
            typedef ConstantNullOp NullOp;
            typedef typename ConstantExpressionTraits<Type>::MetadataType MetadataType;
            
            
        public:
            static void InitializeMetadata(typename boost::call_traits<DataType>::const_reference data, MetadataType& m)
            {
                m = MetadataType(data);
            }
            
            static void Evaluate(Accumulator<ResultType>& accum, typename boost::call_traits<DataType>::const_reference d)
            {
                *accum = d;
            }
            
            template<template <typename, typename> class ParentOpType>
            static void Evaluate(Accumulator<ResultType>& accum, typename boost::call_traits<DataType>::const_reference d)
            {
                ParentOpType<ResultType, DataType>::ApplyEqual(accum, d);
            }
            
            template<typename CheckType>
            static typename boost::enable_if<boost::is_same<CheckType, Type>, bool>::type
            ContainsReference(const CheckType& result, DataType m_data)
            {
                return &result == &m_data;
            }
            
            template<typename CheckType>
            static typename boost::disable_if<boost::is_same<CheckType, Type>, bool>::type
            ContainsReference(const CheckType& result, DataType m_data)
            {
                return false;
            }
            
            static void Print(std::ostream& os, typename boost::call_traits<DataType>::const_reference data)
            {
                os << data;
            }
    };
    
    template<typename DataType>
    Expression<ConstantExpressionPolicy<DataType> > MakeExpr(const DataType& rhs)
    {
        return Expression<ConstantExpressionPolicy<DataType> >(rhs);
    }
    
    template<typename T>
    struct IsConstantExpressionPolicy : public boost::false_type {};
    
    template<typename T>
    struct IsConstantExpressionPolicy<ConstantExpressionPolicy<T> > : public boost::true_type {};
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif // NEKTAR_LIB_UTILITIES_CONSTANT_EXPRESSION_HPP

/**
    $Log: ConstantExpression.hpp,v $
    Revision 1.15  2008/01/03 04:16:41  bnelson
    Changed method name in the expression library from Apply to Evaluate.

    Revision 1.14  2008/01/02 20:11:22  bnelson
    Fixed error detecting if incoming type was the same as the constant expression's type when looking for aliasing.

    Revision 1.13  2008/01/02 18:20:12  bnelson
    *** empty log message ***

    Revision 1.12  2007/12/19 05:09:21  bnelson
    First pass at detecting aliasing.  Still need to test performance implications.

    Revision 1.11  2007/11/13 18:05:27  bnelson
    Added MakeExpr helper function.

    Revision 1.10  2007/10/04 03:48:54  bnelson
    *** empty log message ***

    Revision 1.9  2007/08/16 02:14:21  bnelson
    Moved expression templates to the Nektar namespace.

    Revision 1.8  2007/01/30 23:37:16  bnelson
    *** empty log message ***

    Revision 1.7  2007/01/16 17:37:55  bnelson
    Wrapped everything with #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

    Revision 1.6  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.5  2006/09/14 02:08:59  bnelson
    Fixed gcc compile errors

    Revision 1.4  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

    Revision 1.3  2006/08/28 02:39:53  bnelson
    *** empty log message ***

    Revision 1.2  2006/08/25 01:33:47  bnelson
    no message

    Revision 1.1  2006/06/01 09:20:55  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:41  kirby
    *** empty log message ***

    Revision 1.2  2006/01/31 13:51:12  bnelson
    Updated for new configure.

    Revision 1.1  2006/01/10 14:50:30  bnelson
    Initial Revision.

**/

