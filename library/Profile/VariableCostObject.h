
#ifndef NEKTAR_PROFILE_VARIABLE_COST_OBJECT_H
#define NEKTAR_PROFILE_VARIABLE_COST_OBJECT_H

#include <boost/pool/pool.hpp>

#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>


class VariableCostObject
{
    public:
        VariableCostObject() : data(0)
        {
            if( !pool )
            {
                pool = new boost::pool<>(size*sizeof(double), 32);
            }
            data = (double*)pool->malloc();
            
            for(unsigned int j = 0; j < size; ++j)
            {
                data[j] = 0.0;
            }
        }

        VariableCostObject(const VariableCostObject& rhs) : data(0)
        {
            if( !pool )
            {
                pool = new boost::pool<>(size*sizeof(double), 32);
            }

            data = (double*)pool->malloc();
            for(unsigned int j = 0; j < size; ++j)
            {
                data[j] += rhs.data[j];
            }
        }

        VariableCostObject& operator=(const VariableCostObject& rhs)
        {
            if( &rhs == this )
            {
                return *this;
            }

            for(unsigned int j = 0; j < size; ++j)
            {
                data[j] += rhs.data[j];
            }

            return *this;
        }
        ~VariableCostObject()
        {
            pool->free(data);
            data = 0;
        }

        #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        template<typename PolicyType>
        VariableCostObject(const Nektar::Expression<PolicyType>& rhs)
        {
            if( !pool )
            {
                pool = new boost::pool<>(size*sizeof(double), 32);
            }

            data = (double*)pool->malloc();
            rhs.Evaluate(*this);
        }
        #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
        
        #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        template<typename PolicyType>
        VariableCostObject& operator=(const Nektar::Expression<PolicyType>& rhs)
        {
            rhs.Evaluate(*this);
            return *this;
        }
        #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
        

        static int slowdownFactor;
        static unsigned int size;

        double* data;
 
        
        static boost::pool<>* pool;

};

unsigned int VariableCostObject::size = 10;
int VariableCostObject::slowdownFactor = 1;
boost::pool<>* VariableCostObject::pool = 0;

void NekAdd(VariableCostObject& result, const VariableCostObject& lhs, const VariableCostObject& rhs)
{
    for(unsigned int i = 0; i < VariableCostObject::size; ++i)
    {
        result.data[i] = lhs.data[i] + rhs.data[i];
    }
}
    
void NekAddEqual(VariableCostObject& result, const VariableCostObject& rhs)
{
    for(unsigned int i = 0; i < VariableCostObject::size; ++i)
    {
        result.data[i] += rhs.data[i];
    }
}

VariableCostObject NekAdd(const VariableCostObject& lhs, const VariableCostObject& rhs)
{
    VariableCostObject result;
    NekAdd(result, lhs, rhs);
    return result;
}

GENERATE_ADDITION_OPERATOR(VariableCostObject, 0, VariableCostObject, 0);

#endif //NEKTAR_PROFILE_VARIABLE_COST_OBJECT_H

