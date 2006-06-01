////////////////////////////////////////////////////////////////////////////////
//
// ConstantExpression.hpp
// Blake Nelson
//
// Expression templates for constant expressions.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_CONSTANT_EXPRESSION_HPP
#define NEKTAR_LIB_UTILITIES_CONSTANT_EXPRESSION_HPP

#include <boost/call_traits.hpp>

namespace Nektar
{
    // Expressions can support constant, unary, and binary expressions.
    template<typename ResultType>
    class ConstantExpression
    {
        public:
            typedef ResultType result_type;

        public:
            explicit inline ConstantExpression(typename boost::call_traits<ResultType>::const_reference value) :
                m_value(value)
            {
            }

            inline ConstantExpression(const ConstantExpression& rhs) :
                m_value(rhs.m_value)
            {
            }

            inline ~ConstantExpression() {}

            inline ConstantExpression& operator=(const ConstantExpression& rhs)
            {
                m_value = rhs.m_value;
            }

            inline typename boost::call_traits<ResultType>::const_reference getValue() const { return m_value; }

        private:
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
    Revision 1.1  2006/05/04 18:57:41  kirby
    *** empty log message ***

    Revision 1.2  2006/01/31 13:51:12  bnelson
    Updated for new configure.

    Revision 1.1  2006/01/10 14:50:30  bnelson
    Initial Revision.

**/

