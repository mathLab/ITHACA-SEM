///////////////////////////////////////////////////////////////////////////////
//
// File: TestBinaryExpressionEvaluator_ConstantConstantSpecialization0.cpp
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

#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <UnitTests/CountedObject.h>

namespace Nektar
{
    namespace ConstantConstantSpecialization0
    {
        class A : public CountedObject<A>
        {
            public:
                A() : CountedObject<A>() {}
                A(int v) : CountedObject<A>(v)  {}
                A(const A& rhs) : CountedObject<A>(rhs) {}
                A& operator=(const A& rhs) { CountedObject<A>::operator=(rhs); return *this; }
        };

        class B : public CountedObject<B>
        {
            public:
                B() : CountedObject<B>() {}
                B(int v) : CountedObject<B>(v)  {}
                B(const B& rhs) : CountedObject<B>(rhs) {}
                B& operator=(const B& rhs) { CountedObject<B>::operator=(rhs); return *this; }
        };
        
        A NekAdd(const A& lhs, const B& rhs)
        {
            return A(lhs.value + rhs.value);
        }
        
        void NekAdd(A& result, const A& lhs, const B& rhs)
        {
            result.value = lhs.value + rhs.value;
        }
        
        GENERATE_ADDITION_OPERATOR(A, 0, B, 0);
        
        BOOST_AUTO_TEST_CASE(TestAllSameType)
        {
            typedef ConstantExpressionPolicy<A> LhsPolicyType;
            typedef ConstantExpressionPolicy<B> RhsPolicyType;
            
            A obj1(78);
            B obj2(65);
            
            Expression<LhsPolicyType> lhs = MakeExpr(obj1);
            Expression<RhsPolicyType> rhs = MakeExpr(obj2);
            
            BOOST_STATIC_ASSERT(( BinaryExpressionEvaluator<LhsPolicyType, RhsPolicyType, A, AddOp>::ClassId == 100 ));
            
            A result;
            Accumulator<A> accum(result);
            
            CountedObject<A>::ClearCounters();
            CountedObject<B>::ClearCounters();
            
            BinaryExpressionEvaluator<LhsPolicyType, RhsPolicyType, A, AddOp>::Eval(lhs, rhs, accum);
            
            BOOST_CHECK_EQUAL((*accum).value, 78+65);
            
            CountedObject<A>::Check(0, 0, 0, 0, 0, 0);
            CountedObject<B>::Check(0, 0, 0, 0, 0, 0);
            
        }
    }
}
