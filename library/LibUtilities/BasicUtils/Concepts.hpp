///////////////////////////////////////////////////////////////////////////////
// 
// Concepts.hpp
//
// Defines concepts for the boost concept checking library.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_CONCEPTS_HPP
#define NEKTAR_LIB_UTILITIES_CONCEPTS_HPP

namespace Nektar
{
    template<typename DataType>
    class AssignableConcept
    {
        public:
            void constraints()
            {
                DataType* lhs = NULL;
                DataType* rhs = NULL;
                *lhs = *rhs;
                const_constraints(*rhs);
            }

            void const_constraints(const DataType& rhs)
            {
                DataType* lhs = NULL;
                (*lhs) = rhs;
            }
    };
}

#endif

/**
    $Log: Concepts.hpp,v $
    Revision 1.1  2006/05/04 18:57:41  kirby
    *** empty log message ***

    Revision 1.2  2006/02/12 15:06:12  sherwin

    Changed .h files to .hpp

    Revision 1.1  2006/01/31 13:51:12  bnelson
    Updated for new configure.

**/
