///////////////////////////////////////////////////////////////////////////////
//
// File: BinaryOperators.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_OPERATORS_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_OPERATORS_HPP

#include <LibUtilities/ExpressionTemplates/ArithmeticTraits.hpp>
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
#include <LibUtilities/ExpressionTemplates/Accumulator.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionConcepts.hpp>

#include <boost/concept_check.hpp>
#include <boost/call_traits.hpp>
#include <boost/utility/result_of.hpp>
#include <string>

namespace Nektar
{
    // Forward declarations needed to declare inverse operators.
    template<typename LhsType, typename RhsType>
    class AddOp;
    
    template<typename LhsType, typename RhsType>
    class SubtractOp;
            
    template<typename LhsType, typename RhsType>
    class MultiplyOp;
            
    template<typename LhsType, typename RhsType>
    class DivideOp;
        
    /// \brief An expression to negate an object of ParameterType.
    /// Parameter type is the actual object, not the expression that may lead to it.
    template<typename LhsType, typename RhsType>
    class AddOp 
    {
        public:
            typedef AdditionTraits<LhsType, RhsType> TraitsType;
            typedef typename TraitsType::ResultType ResultType;

            template<typename L, typename R>
            class Rebind
            {
                public:
                    typedef AddOp<L,R> type;
            };

            static void Apply(Accumulator<ResultType>& result,
                                typename boost::call_traits<LhsType>::const_reference lhs, 
                                typename boost::call_traits<RhsType>::const_reference rhs)
            {
                TraitsType::Add(*result, lhs, rhs);
            }

            static void ApplyEqual(Accumulator<LhsType>& result,
                                    typename boost::call_traits<RhsType>::const_reference rhs)
            {
                //TraitsType::AddEqual(*result, rhs);
                if( result.IsInitialized() )
                {
                    *result += rhs;
                }
                else
                {
                    *result = rhs;
                }
            }
            
            static void ApplyLeftEqual(Accumulator<ResultType>& result,
                                        typename boost::call_traits<LhsType>::const_reference lhs)
            {
                if( result.IsInitialized() )
                {
                    TraitsType::AddLeftEqual(*result, lhs);
                }
                else
                {
                    *result = lhs;
                }
            }

            static const std::string& AsString()
            {
                return s_StringRep;
            }

            static const unsigned int Priority = 1;
        private:
            static std::string s_StringRep;
    };

    template<typename LhsType, typename RhsType>
    std::string AddOp<LhsType, RhsType>::s_StringRep("+");
            
    template<typename LhsType, typename RhsType>
    class MultiplyOp 
    {
        public:
            typedef MultiplicationTraits<LhsType, RhsType> TraitsType;
            typedef typename TraitsType::ResultType ResultType;
            
            template<typename L, typename R>
            class Rebind
            {
                public:
                    typedef MultiplyOp<L,R> type;
            };
            
            static void Apply(Accumulator<ResultType>& result,
                                typename boost::call_traits<LhsType>::const_reference lhs, 
                                typename boost::call_traits<RhsType>::const_reference rhs)
            {
                TraitsType::Multiply(*result, lhs, rhs);
            }
            
            static void ApplyEqual(typename boost::call_traits<LhsType>::const_reference lhs, 
                                    typename boost::call_traits<RhsType>::const_reference rhs,
                                    Accumulator<ResultType>& result)
            {
                TraitsType::MultiplyEqual(lhs, rhs, *result);
            }
            
            static void ApplyEqual(Accumulator<LhsType>& result,
                                    typename boost::call_traits<RhsType>::const_reference rhs)
            {
                if( result.IsInitialized() )
                {
                    *result *= rhs;
                }
                else
                {
                    *result = rhs;
                }
            }

            static const std::string& AsString()
            {
                return s_StringRep;
            }
            
            static const unsigned int Priority = 2;
        private:
            static std::string s_StringRep;
    };
    
    template<typename LhsType, typename RhsType>
    std::string MultiplyOp<LhsType, RhsType>::s_StringRep("*");

    template<typename LhsType, typename RhsType>
    class DivideOp 
    {
        public:
            typedef DivisionTraits<LhsType, RhsType> TraitsType;
            typedef typename TraitsType::ResultType ResultType;

            
            template<typename L, typename R>
            class Rebind
            {
                public:
                    typedef DivideOp<L,R> type;
            };

            
            static void Apply(Accumulator<ResultType>& result,
                                typename boost::call_traits<LhsType>::const_reference lhs, 
                                typename boost::call_traits<RhsType>::const_reference rhs)
            {
                TraitsType::Divide(*result, lhs, rhs);
            }
            
            static void ApplyEqual(typename boost::call_traits<LhsType>::const_reference lhs, 
                                    typename boost::call_traits<RhsType>::const_reference rhs,
                                    Accumulator<ResultType>& result)
            {
                TraitsType::DivideEqual(lhs, rhs, *result);
            }

            static void ApplyEqual(Accumulator<LhsType>& result,
                                    typename boost::call_traits<RhsType>::const_reference rhs)
            {
                //TraitsType::DivideEqual(*result, rhs);
                if( result.IsInitialized() )
                {
                    *result /= rhs;
                }
                else
                {
                    *result = rhs;
                }
            }
            
            static const std::string& AsString()
            {
                return s_StringRep;
            }
            
            static const unsigned int Priority = 2;
        private:
            static std::string s_StringRep;
    };

    template<typename LhsType, typename RhsType>
    std::string DivideOp<LhsType, RhsType>::s_StringRep("/");
            
    template<typename LhsType, typename RhsType>
    class SubtractOp 
    {
        public:
            typedef SubtractionTraits<LhsType, RhsType> TraitsType;
            typedef typename TraitsType::ResultType ResultType;
            
            
            template<typename L, typename R>
            class Rebind
            {
                public:
                    typedef SubtractOp<L,R> type;
            };
            
            
            static void Apply(Accumulator<ResultType>& result,
                                typename boost::call_traits<LhsType>::const_reference lhs, 
                                typename boost::call_traits<RhsType>::const_reference rhs)
            {
                TraitsType::Subtract(*result, lhs, rhs);
            }
            
            static void ApplyEqual(typename boost::call_traits<LhsType>::const_reference lhs, 
                                    typename boost::call_traits<RhsType>::const_reference rhs,
                                    Accumulator<ResultType>& result)
            {
                TraitsType::SubtractEqual(lhs, rhs, *result);
            }

            static void ApplyEqual(Accumulator<LhsType>& result,
                                    typename boost::call_traits<RhsType>::const_reference rhs)
            {
                if( result.IsInitialized() )
                {
                    *result -= rhs;
                }
                else
                {
                    *result = rhs;
                }
            }
            

            static const std::string& AsString()
            {
                return s_StringRep;
            }
            
            static const unsigned int Priority = 1;
        private:
            static std::string s_StringRep;

    };
    
    template<typename LhsType, typename RhsType>
    std::string SubtractOp<LhsType, RhsType>::s_StringRep("-");
}

#endif // NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_BINARY_EXPRESSION_OPERATORS_HPP


