///////////////////////////////////////////////////////////////////////////////
//
// File: CommutativeTraits.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_COMMUTATIVE_TRAITS_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_COMMUTATIVE_TRAITS_HPP
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <boost/type_traits.hpp>

namespace Nektar
{
    template<typename FirstType,
                template <typename, typename> class OpType,
                typename SecondType>
    class IsCommutative : public boost::false_type
    {
        public:
            static const bool IsComm = false;
    };
    
    template<typename FirstType, typename SecondType>
    class IsCommutative<FirstType, AddOp, SecondType> : public boost::true_type
    {
        public:
            static const bool IsComm = true;
    };
    
    template<typename FirstType, typename SecondType>
    class IsCommutative<FirstType, MultiplyOp, SecondType> : public boost::true_type
    {
        public:
            static const bool IsComm = true;
    };
    
    template<typename LhsPolicy, template <typename, typename> class OpType, typename RhsPolicy>
    class IsCommutative<Expression<LhsPolicy>, OpType, Expression<RhsPolicy> > : public boost::mpl::if_
                                                                                <
                                                                                    boost::is_base_of
                                                                                    <
                                                                                        boost::true_type,
                                                                                        IsCommutative<typename Expression<LhsPolicy>::ResultType, OpType, typename RhsPolicy::ResultType>
                                                                                    >,
                                                                                    boost::true_type,
                                                                                    boost::false_type
                                                                                >
     {
      public:
            static const bool IsComm = true;
     };
                                                                
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_COMMUTATIVE_TRAITS_HPP
