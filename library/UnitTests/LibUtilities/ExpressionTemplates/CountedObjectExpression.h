///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_UNIT_TESTS_COUNTED_OBJECT_EXPRESSION_EXPRESSION_TEMPLATE_UNIT_TESTS_COUNTED_OBJECT_EXPRESSION_H
#define NEKTAR_UNIT_TESTS_COUNTED_OBJECT_EXPRESSION_EXPRESSION_TEMPLATE_UNIT_TESTS_COUNTED_OBJECT_EXPRESSION_H

#include <UnitTests/CountedObject.h>
#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/Operators.hpp>

#include <boost/typeof/typeof.hpp>
#include  BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

BOOST_TYPEOF_REGISTER_TEMPLATE(Nektar::CountedObject, 1);

namespace Nektar
{
    template<typename T>
    expt::Node<expt::Node<CountedObject<T> >, expt::AddOp, expt::Node<CountedObject<T> > >
    operator+(const CountedObject<T>& lhs, const CountedObject<T>& rhs)
    {
        typedef expt::Node<CountedObject<T> > lhs_type;
        typedef expt::Node<CountedObject<T> > rhs_type;
        
        lhs_type lhs_node(lhs);
        rhs_type rhs_node(rhs);
        return expt::Node<lhs_type, expt::AddOp, rhs_type>(lhs_node, rhs_node);
    }

    template<typename T>
    CountedObject<T> Add(const CountedObject<T>& lhs, const CountedObject<T>& rhs)
    {
        return CountedObject<T>(lhs.value + rhs.value);
    }

    template<typename T>
    void AddEqual(CountedObject<T>& accumulator, const CountedObject<T>& rhs)
    {
        accumulator.value += rhs.value;
    }

    template<typename T>
    void Add(CountedObject<T>& accumulator, const CountedObject<T>& lhs, const CountedObject<T>& rhs)
    {
        accumulator.value = lhs.value + rhs.value;
    }
}

#endif //NEKTAR_UNIT_TESTS_COUNTED_OBJECT_EXPRESSION_EXPRESSION_TEMPLATE_UNIT_TESTS_COUNTED_OBJECT_EXPRESSION_H
