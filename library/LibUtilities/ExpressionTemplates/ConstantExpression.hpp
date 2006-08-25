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
    // Expressions can support constant, unary, and binary expressions.
    template<typename ResultType>
    class ConstantExpression : public Expression<ResultType>
    {
        public:
            typedef ResultType result_type;

        public:
            /// \note The parameter must not be a temporary or strange behavior will result.
            explicit inline ConstantExpression(typename boost::call_traits<ResultType>::const_reference value) :
                m_value(value)
            {
            }

            inline ConstantExpression(const ConstantExpression& rhs) :
                m_value(rhs.m_value)
            {
            }

            inline ~ConstantExpression() {}

            typename boost::call_traits<ResultType>::const_reference operator*() const { return m_value; }
            typename boost::call_traits<ResultType>::const_reference GetValue() const { return m_value; }

        private:
            void DoApply(typename boost::call_traits<ResultType>::reference result) const
            {
                result = m_value;
            }

            ConstantExpression<ResultType>& operator=(const ConstantExpression<ResultType>& rhs);
            typename boost::call_traits<ResultType>::const_reference m_value;
    };

    template<typename DataType>
    class GetReturnType<ConstantExpression<DataType> >
    {
        public:
            typedef typename ConstantExpression<DataType>::result_type result_type;
    };
}

#endif // NEKTAR_LIB_UTILITIES_CONSTANT_EXPRESSION_HPP


/**
    $Log: ConstantExpression.hpp,v $
    Revision 1.1  2006/06/01 09:20:55  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:41  kirby
    *** empty log message ***

    Revision 1.2  2006/01/31 13:51:12  bnelson
    Updated for new configure.

    Revision 1.1  2006/01/10 14:50:30  bnelson
    Initial Revision.

**/

