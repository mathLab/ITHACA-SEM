////////////////////////////////////////////////////////////////////////////////
// 
// exprt.hpp
// Blake Nelson
//
// Generic expression templates.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_TRAITS_H
#define NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_TRAITS_H

namespace Nektar
{
    namespace expt
    {
        /// \brief Defines the result type for various binary operations.
        ///
        /// Given a function or functor f and its arguments, there is 
        /// insufficient support
        /// in the C++ language to reliably determine the result of the the 
        /// expression (although some limited support can be found in the 
        /// boost::result_of library).
        ///
        /// This traits class is used to define the result type for the 
        /// operations +, -, /, * for an arbitrary combination of types.
        ///
        /// The result type of an operation is specified 
        /// TODO - finish docs.
        template<typename LhsType, typename RhsType>
        class BinaryExpressionTraits
        {
            public:
                typedef LhsType AdditionResultType;
                typedef LhsType SubtractionResultType;
                typedef LhsType DivisionResultType;
                typedef LhsType MultiplicationResultType;
                
                static const bool AdditionIsAssociative = true;
        };
    }
}

#endif // NEKTAR_LIB_UTILITIES_BINARY_EXPRESSION_TRAITS_H

/**
    $Log: BinaryExpressionTraits.hpp,v $
    Revision 1.4  2006/09/14 02:08:59  bnelson
    Fixed gcc compile errors

    Revision 1.3  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

    Revision 1.2  2006/08/25 01:33:47  bnelson
    no message

    Revision 1.1  2006/06/01 09:20:55  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:41  kirby
    *** empty log message ***

    Revision 1.1  2006/01/10 14:50:30  bnelson
    Initial Revision.

**/

