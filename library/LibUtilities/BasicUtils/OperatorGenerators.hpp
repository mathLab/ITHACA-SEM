///////////////////////////////////////////////////////////////////////////////
//
// File: OperatorGenerators
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
// Description: Generates global operators for classes based upon traits.
//
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_OPERATOR_GENERATORS_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_OPERATOR_GENERATORS_HPP

#include <LibUtilities/BasicUtils/BinaryExpressionTraits.hpp>
#include <boost/concept_check.hpp>

namespace Nektar
{
    template<typename LhsType, typename RhsType, typename OpType>
    class OperatorGenerator;
    
//     template<typename LhsType, typename RhsType>
//     struct PlusEqualConcept
//     {
//         void constraints()
//         {
//             lhs += rhs;
//         }
//         
//         LhsType lhs;
//         RhsType rhs;
//     };
    
    template<typename LhsType, typename RhsType>
    class OperatorGenerator<LhsType, RhsType, AddOp>
    {
        public:
            typedef typename BinaryExpressionTraits<LhsType, RhsType, AddOp>::ResultType ResultType;
            
            friend ResultType operator+(typename boost::call_traits<LhsType>::const_reference lhs,
                                        typename boost::call_traits<RhsType>::const_reference rhs)
            {
                //boost::function_requires<PlusEqualConcept<LhsType, RhsType> >();
                ResultType result;
                NekAdd(result, lhs, rhs);
                return result;
            }
    };
    
    template<typename LhsType, typename RhsType>
    class OperatorGenerator<LhsType, RhsType, SubtractOp>
    {
        public:
            typedef typename BinaryExpressionTraits<LhsType, RhsType, SubtractOp>::ResultType ResultType;
            
            friend ResultType operator-(typename boost::call_traits<LhsType>::const_reference lhs,
                                        typename boost::call_traits<RhsType>::const_reference rhs)
            {
                ResultType result;
                NekSubtract(result, lhs, rhs);
                return result;
            }
    };
    
    template<typename LhsType, typename RhsType>
    class OperatorGenerator<LhsType, RhsType, DivideOp>
    {
        public:
            typedef typename BinaryExpressionTraits<LhsType, RhsType, DivideOp>::ResultType ResultType;
            
            friend ResultType operator/(typename boost::call_traits<LhsType>::const_reference lhs,
                                        typename boost::call_traits<RhsType>::const_reference rhs)
            {
                ResultType result;
                NekDivide(result, lhs, rhs);
                return result;
            }
    };
    
//     template<typename LhsType, typename RhsType>
//     class OperatorGenerator<LhsType, RhsType, MultiplyOp>
//     {
//         public:
//             typedef typename BinaryExpressionTraits<LhsType, RhsType, MultiplyOp>::ResultType ResultType;
//             
//             friend ResultType operator*(typename boost::call_traits<LhsType>::const_reference lhs,
//                                         typename boost::call_traits<RhsType>::const_reference rhs)
//             {
//                 ResultType result;
//                 NekMultiply(result, lhs, rhs);
//                 return result;
//             }
//     };
    
    #define GENERATE_MULTIPLICATION_OPERATOR(LhsType, RhsType) \
            BinaryExpressionTraits<LhsType, RhsType, MultiplyOp>::ResultType \
            operator*(const LhsType& lhs, const RhsType& rhs) \
            { \
                BinaryExpressionTraits<LhsType, RhsType, MultiplyOp>::ResultType result; \
                NekMultiply(result, lhs, rhs); \
                return result; \
            }
    
    #define GENERATE_MULTIPLICATION_OPERATOR_L3(LhsType, RhsType) \
            template<typename L1, typename L2, typename L3> \
            typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType, MultiplyOp>::ResultType \
            operator*(const LhsType<L1, L2, L3>& lhs, const RhsType& rhs) \
            { \
                typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType, MultiplyOp>::ResultType result; \
                NekMultiply(result, lhs, rhs); \
                return result; \
            }
    
    #define GENERATE_MULTIPLICATION_OPERATOR_L3_R0(LhsType) \
            template<typename L1, typename L2, typename L3, typename RhsType> \
            typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType, MultiplyOp>::ResultType \
            operator*(const LhsType<L1, L2, L3>& lhs, const RhsType& rhs) \
            { \
                typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType, MultiplyOp>::ResultType result; \
                NekMultiply(result, lhs, rhs); \
                return result; \
            }
    
    #define GENERATE_MULTIPLICATION_OPERATOR_L0_R3(RhsType) \
            template<typename R1, typename R2, typename R3, typename LhsType> \
            typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, MultiplyOp>::ResultType \
            operator*(const LhsType& lhs, const RhsType<R1, R2, R3>& rhs) \
            { \
                typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, MultiplyOp>::ResultType result; \
                NekMultiply(result, lhs, rhs); \
                return result; \
            }
    
//     template<typename LhsType, template<typename, typename, typename> class RhsType, typename OpType>
//     class OperatorGeneratorR3;
//     
//     template<typename LhsType, template<typename, typename, typename> class RhsType>
//     class OperatorGenerator<LhsType, RhsType, MultiplyOp>
//     {
//         public:
//             //typedef typename BinaryExpressionTraits<LhsType, RhsType, MultiplyOp>::ResultType ResultType;
//             
//             template<typename R1, typename R2, typename R3>
//             friend typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3> >::ResultType 
//             operator*(typename boost::call_traits<LhsType>::const_reference lhs,
//                       typename boost::call_traits<RhsType<R1, R2, R3> >::const_reference rhs)
//             {
//                 typedef typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3> >::ResultType ResultType;
//                 ResultType result;
//                 ResultType::NekMultiply(result, lhs, rhs);
//                 return result;
//             }
//     };
    
    template<template<typename, typename, typename> class LhsType, typename RhsType, typename OpType>
    class OperatorGeneratorL3;
    
    template<template<typename, typename, typename> class LhsType, typename RhsType>
    class OperatorGeneratorL3<LhsType, RhsType, MultiplyOp>
    {
        public:           
            template<typename L1, typename L2, typename L3>
            friend typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType, MultiplyOp>::ResultType 
            operator*(const LhsType<L1, L2, L3>& lhs,
                      const RhsType& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType, MultiplyOp>::ResultType ResultType;
                ResultType result;
                NekMultiply(result, lhs, rhs);
                return result;
            }
    };
    
    
    template<template<typename, typename, typename> class LhsType, 
             template<typename, typename, typename> class RhsType, 
             typename OpType>
    class OperatorGeneratorL3R3;
    
    template<template<typename, typename, typename> class LhsType, 
             template<typename, typename, typename> class RhsType>
    class OperatorGeneratorL3R3<LhsType, RhsType, MultiplyOp>
    {
        public:           
            template<typename L1, typename L2, typename L3,
                     typename R1, typename R2, typename R3>
            friend typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType<R1, R2, R3>, MultiplyOp>::ResultType 
            operator*(const LhsType<L1, L2, L3>& lhs,
                      const RhsType<R1, R2, R3>& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType<R1, R2, R3>, MultiplyOp>::ResultType ResultType;
                ResultType result;
                NekMultiply(result, lhs, rhs);
                return result;
            }
    };
    
    template<template<typename, typename, typename> class LhsType, 
             template<typename, typename, typename> class RhsType>
    class OperatorGeneratorL3R3<LhsType, RhsType, AddOp>
    {
        public:           
            template<typename L1, typename L2, typename L3,
                     typename R1, typename R2, typename R3>
            friend typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType<R1, R2, R3>, AddOp>::ResultType 
            operator+(const LhsType<L1, L2, L3>& lhs,
                      const RhsType<R1, R2, R3>& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType<R1, R2, R3>, AddOp>::ResultType ResultType;
                ResultType result;
                NekAdd(result, lhs, rhs);
                return result;
            }
    };
    
    template<template<typename, typename, typename> class LhsType, 
             template<typename, typename, typename> class RhsType>
    class OperatorGeneratorL3R3<LhsType, RhsType, SubtractOp>
    {
        public:           
            template<typename L1, typename L2, typename L3,
                     typename R1, typename R2, typename R3>
            friend typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType<R1, R2, R3>, SubtractOp>::ResultType 
            operator-(const LhsType<L1, L2, L3>& lhs,
                      const RhsType<R1, R2, R3>& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType<R1, R2, R3>, SubtractOp>::ResultType ResultType;
                ResultType result;
                NekSubtract(result, lhs, rhs);
                return result;
            }
    };
    
    template<template<typename, typename, typename> class LhsType, 
             template<typename, typename, typename> class RhsType>
    class OperatorGeneratorL3R3<LhsType, RhsType, DivideOp>
    {
        public:           
            template<typename L1, typename L2, typename L3,
                     typename R1, typename R2, typename R3>
            friend typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType<R1, R2, R3>, DivideOp>::ResultType 
            operator/(const LhsType<L1, L2, L3>& lhs,
                      const RhsType<R1, R2, R3>& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType<L1, L2, L3>, RhsType<R1, R2, R3>, DivideOp>::ResultType ResultType;
                ResultType result;
                NekDivide(result, lhs, rhs);
                return result;
            }
    };
    
    template<typename LhsType, 
             template<typename, typename, typename> class RhsType>
    class OperatorGeneratorR3
    {
        public:
            template<typename R1, typename R2, typename R3>
            friend typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, DivideOp>::ResultType 
            operator/(const LhsType& lhs,
                      const RhsType<R1, R2, R3>& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, DivideOp>::ResultType ResultType;
                ResultType result;
                NekDivide(result, lhs, rhs);
                return result;
            }
            
            template<typename R1, typename R2, typename R3>
            friend typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, MultiplyOp>::ResultType 
            operator*(const LhsType& lhs,
                      const RhsType<R1, R2, R3>& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, MultiplyOp>::ResultType ResultType;
                ResultType result;
                NekMultiply(result, lhs, rhs);
                return result;
            }
            
            template<typename R1, typename R2, typename R3>
            friend typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, AddOp>::ResultType 
            operator+(const LhsType& lhs,
                      const RhsType<R1, R2, R3>& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, AddOp>::ResultType ResultType;
                ResultType result;
                NekAdd(result, lhs, rhs);
                return result;
            }
            
            template<typename R1, typename R2, typename R3>
            friend typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, SubtractOp>::ResultType 
            operator-(const LhsType& lhs,
                      const RhsType<R1, R2, R3>& rhs)
            {
                typedef typename BinaryExpressionTraits<LhsType, RhsType<R1, R2, R3>, SubtractOp>::ResultType ResultType;
                ResultType result;
                NekSubtract(result, lhs, rhs);
                return result;
            }
    };
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_OPERATOR_GENERATORS_HPP

