
#ifndef NEKTAR_PROFILE_INT_WRAPPER_H
#define NEKTAR_PROFILE_INT_WRAPPER_H

#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>

class IntWrapper
{
    public:
        IntWrapper() : m_value(0) {}
        
        explicit IntWrapper(int v) : m_value(v) {}
                
        #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        template<typename PolicyType>
        IntWrapper(const Nektar::Expression<PolicyType>& rhs)
        {
            rhs.Evaluate(*this);
        }
        #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
        
        int GetValue() const { return m_value; }
        
        #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        template<typename PolicyType>
        IntWrapper& operator=(const Nektar::Expression<PolicyType>& rhs)
        {
            rhs.Evaluate(*this);
            return *this;
        }
        #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
        
        IntWrapper& operator+=(const IntWrapper& rhs)
        {
            m_value += rhs.m_value;
            return *this;
        }
        
        void SetValue(int v) { m_value = v; }
    private:
        int m_value;
};
      
void NekAdd(IntWrapper& result, const IntWrapper& lhs, const IntWrapper& rhs)
{
    result.SetValue(lhs.GetValue() + rhs.GetValue());
}
    
void NekAddEqual(IntWrapper& result, const IntWrapper& rhs)
{
    result += rhs;
}

IntWrapper NekAdd(const IntWrapper& lhs, const IntWrapper& rhs)
{
    return IntWrapper(lhs.GetValue() + rhs.GetValue());
}
    
GENERATE_ADDITION_OPERATOR(IntWrapper, 0, IntWrapper, 0);

#endif //NEKTAR_PROFILE_INT_WRAPPER_H
