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

#include <boost/call_traits.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionReturnType.hpp>
#include <LibUtilities/ExpressionTemplates/Expression.hpp>


namespace Nektar
{
    template<typename DataType>
    class ConstantExpressionPolicy
    {
    };

    template<typename DataType>
    class Expression<ConstantExpressionPolicy<DataType> >
    {
        public:
            typedef DataType ResultType;
            typedef typename ExpressionMetadataChooser<DataType>::MetadataType MetadataType;
            typedef Expression<ConstantExpressionPolicy<DataType> > ThisType;

        public:
            /// \note The parameter must not be a temporary or strange behavior will result.
            explicit Expression(typename boost::call_traits<ResultType>::const_reference value) :
                m_value(value),
                m_metadata(value)
            {
            }

            Expression(const ThisType& rhs) :
                m_value(rhs.m_value),
                m_metadata(rhs.m_metadata)
            {
            }

            ~Expression() {}

            typename boost::call_traits<ResultType>::const_reference operator*() const { return m_value; }
            typename boost::call_traits<ResultType>::const_reference GetValue() const { return m_value; }

            void Apply(typename boost::call_traits<ResultType>::reference result) const
            {
                result = m_value;
            }

            template<typename OpType>
            void ApplyEqual(typename boost::call_traits<ResultType>::reference result) const
            {
                OpType::ApplyEqual(result, m_value);
            }

            const MetadataType& GetMetadata() const
            {
                return m_metadata;
            }

        private:
            ThisType& operator=(const ThisType& rhs);

            typename boost::call_traits<ResultType>::const_reference m_value;
            MetadataType m_metadata;
    };
}

#endif // NEKTAR_LIB_UTILITIES_CONSTANT_EXPRESSION_HPP


/**
    $Log: ConstantExpression.hpp,v $
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

