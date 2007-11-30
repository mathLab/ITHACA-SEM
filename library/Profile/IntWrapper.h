
#ifndef NEKTAR_PROFILE_INT_WRAPPER_H
#define NEKTAR_PROFILE_INT_WRAPPER_H

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

class IntWrapper
{
    public:
        IntWrapper() : m_value(0) {}
        
        explicit IntWrapper(int v) : m_value(v) {}
                
        #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        template<typename PolicyType>
        explicit IntWrapper(const Nektar::Expression<PolicyType>& rhs)
        {
            rhs.Apply(*this);
        }
        #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
        
        int GetValue() const { return m_value; }
        
        #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        template<typename PolicyType>
        IntWrapper& operator=(const Nektar::Expression<PolicyType>& rhs)
        {
            rhs.Apply(*this);
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

void AddIntWrapper(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2);
void AddIntWrapper(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3);
void AddIntWrapper(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4);
void AddIntWrapper(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5);
void AddIntWrapper(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5,
                  const IntWrapper& r6);                  
void AddIntWrapper(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5,
                  const IntWrapper& r6,
                  const IntWrapper& r7);

void AddIntWrapperManualAccum(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2);
void AddIntWrapperManualAccum(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3);
void AddIntWrapperManualAccum(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4);
void AddIntWrapperManualAccum(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5);
void AddIntWrapperManualAccum(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5,
                  const IntWrapper& r6);                  
void AddIntWrapperManualAccum(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5,
                  const IntWrapper& r6,
                  const IntWrapper& r7);
                  
void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2);
void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3);
void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4);
void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5);
void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5,
                  const IntWrapper& r6);                  
void AddIntWrapperExprTemp(IntWrapper& result, const IntWrapper& r1, const IntWrapper& r2,
                  const IntWrapper& r3,
                  const IntWrapper& r4,
                  const IntWrapper& r5,
                  const IntWrapper& r6,
                  const IntWrapper& r7);  
                                
#endif //NEKTAR_PROFILE_INT_WRAPPER_H
