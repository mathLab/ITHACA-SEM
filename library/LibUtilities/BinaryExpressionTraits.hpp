////////////////////////////////////////////////////////////////////////////////
// 
// exprt.hpp
// Blake Nelson
//
// Generic expression templates.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_EXPRT_TRAITS_H
#define NEKTAR_LIB_UTILITIES_EXPRT_TRAITS_H

namespace Nektar
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
	};
}

#endif // NEKTAR_LIB_UTILITIES_EXPRT_TRAITS_H

/**
    $Log: BinaryExpressionTraits.hpp,v $
    Revision 1.1  2006/01/10 14:50:30  bnelson
    Initial Revision.

**/

