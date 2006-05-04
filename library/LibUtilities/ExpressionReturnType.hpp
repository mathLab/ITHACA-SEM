////////////////////////////////////////////////////////////////////////////////
//
// GetReturnType
// Blake Nelson
//
// Template to determine the return type of an expression.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_GET_RETURN_TYPE_H
#define NEKTAR_LIB_UTILITIES_GET_RETURN_TYPE_H

namespace Nektar
{
    template<typename DataType>
    class GetReturnType
    {
        public:
            typedef DataType result_type;
    };
}

#endif // NEKTAR_LIB_UTILITIES_GET_RETURN_TYPE_H

/**
    $Log: ExpressionReturnType.hpp,v $
    Revision 1.1  2006/01/10 14:50:31  bnelson
    Initial Revision.

**/

