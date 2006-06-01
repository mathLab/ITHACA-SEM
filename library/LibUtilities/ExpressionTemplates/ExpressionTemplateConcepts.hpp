////////////////////////////////////////////////////////////////////////////////
// 
// ExpressionTemplateConcepts.hpp
// Blake Nelson
//
// Generic expression templates concepts.  These are to be used with the 
// boost concept checking classes to make it easier to diagnose problems
// with the expression templates.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_EXPRT_CONCEPTS_H
#define NEKTAR_LIB_UTILITIES_EXPRT_CONCEPTS_H

#include <LibUtilities/ExpressionTemplates/BinaryExpressionTraits.hpp>

namespace Nektar
{
    template<typename LhsType, typename RhsType>
    class OpAddEqualConcept
    {
        public:
            void constraints()
            {
                LhsType* lhs_obj;
                RhsType* rhs_obj;
                
                // All objects using expression templates must defined operator+=.
                *lhs_obj += *rhs_obj;
            }
    };

    template<typename LhsType, typename RhsType>
    class OpAddConcept
    {
        public:
            typedef typename BinaryExpressionTraits<LhsType, RhsType>::op_add_result_type result_type;
        
            void constraints()
            {
                LhsType* lhs_obj;
                RhsType* rhs_obj;
                result_type* result;

                // op_divide must be defined to use the / expression template.
                op_add(*lhs_obj, *rhs_obj, *result);
            }
    };

    template<typename LhsType, typename RhsType>
    class OpMultiplyConcept
    {
        public:
            typedef typename BinaryExpressionTraits<LhsType, RhsType>::op_multiply_result_type result_type;
        
            void constraints()
            {
                LhsType* lhs_obj;
                RhsType* rhs_obj;
                result_type* result;

                // op_divide must be defined to use the / expression template.
                op_multiply(*lhs_obj, *rhs_obj, *result);
            }
    };

    template<typename LhsType, typename RhsType>
    class OpDivideConcept
    {
        public:
            typedef typename BinaryExpressionTraits<LhsType, RhsType>::op_divide_result_type result_type;
        
            void constraints()
            {
                LhsType* lhs_obj;
                RhsType* rhs_obj;
                result_type* result;

                // op_divide must be defined to use the / expression template.
                op_divide(*lhs_obj, *rhs_obj, *result);
            }
    };

    template<typename LhsType, typename RhsType>
    class OpSubtractConcept
    {
        public:
            typedef typename BinaryExpressionTraits<LhsType, RhsType>::op_divide_result_type result_type;
        
            void constraints()
            {
                LhsType* lhs_obj;
                RhsType* rhs_obj;
                result_type* result;

                // op_divide must be defined to use the / expression template.
                op_subtract(*lhs_obj, *rhs_obj, *result);
            }
    };

}

#endif // NEKTAR_LIB_UTILITIES_EXPRT_CONCEPTS_H

/**
    $Log: ExpressionTemplateConcepts.hpp,v $
    Revision 1.1  2006/06/01 09:20:55  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:42  kirby
    *** empty log message ***

    Revision 1.2  2006/01/31 13:51:12  bnelson
    Updated for new configure.

    Revision 1.1  2006/01/10 14:50:31  bnelson
    Initial Revision.

**/

