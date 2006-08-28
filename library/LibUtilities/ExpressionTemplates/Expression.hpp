///////////////////////////////////////////////////////////////////////////////
//
// File: Expression.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_HPP

#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <algorithm>

namespace Nektar
{
    enum ExpressionOperation
    {
        eNEGATION,
        eSUBTRACTION,
        eADDITION,
        eDIVISION,
        eMULTIPLICATION
    };

    template<typename DataType>
    class DefaultExpressionMetadata
    {
        public:
            DefaultExpressionMetadata(typename boost::call_traits<DataType>::const_reference)
            {
            }

            static DefaultExpressionMetadata UpdateForNegation(const DefaultExpressionMetadata& rhs)
            {
                return DefaultExpressionMetadata(rhs);
            }
    };

    template<typename Type, typename enabled=void>
    class ExpressionMetadataChooser
    {
        public:
            typedef DefaultExpressionMetadata<Type> MetadataType;
    };

    template<typename Type>
    class ExpressionMetadataChooser<Type, typename boost::enable_if<
        boost::is_same<typename Type::ExpressionMetadataType, typename Type::ExpressionMetadataType>
        >::type >
    {
        public:
            typedef typename Type::ExpressionMetadataType MetadataType;
    };

//     template<typename ParameterType, unsigned int numberOfParameters, typename ParamType2 = void>
//     class Expression;

//     // An expression that, when fully evaulated, returns an object or value of type "Type".
//     template<typename Type>
//     class Expression
//     {
//         public:
//             typedef Type ResultType;
//             typedef typename ExpressionMetadataChooser<Type>::MetadataType MetadataType;
//
//         public:
//             Expression() :
//                 m_metadata()
//             {
//             }
//
//             Expression(const Expression<Type>& rhs) :
//                 m_metadata(rhs.m_metadata)
//             {
//             }
//
//             virtual ~Expression()
//             {
//             }
//
//             Expression<Type>& operator=(const Expression<Type>& rhs)
//             {
//                 Expression<Type> temp(rhs);
//                 Swap(rhs);
//                 return *this;
//             }
//
//             const MetadataType& GetMetadata() const
//             {
//                 return m_metadata;
//             }
//
//             void Apply(typename boost::call_traits<ResultType>::reference result) const
//             {
//                 return DoApply(result);
//             }
//
//         private:
//             virtual void DoApply(typename boost::call_traits<ResultType>::reference result) const = 0;
//
//             void Swap(Expression<Type>& rhs)
//             {
//                 std::swap(m_metadata, rhs.m_metadata);
//             }
//
//             MetadataType m_metadata;
//     };
}

#endif // NEKTAR_LIB_UTILITIES_EXPRESSION_HPP

/**
    $Log: Expression.hpp,v $
    Revision 1.1  2006/08/25 01:33:47  bnelson
    no message

**/
