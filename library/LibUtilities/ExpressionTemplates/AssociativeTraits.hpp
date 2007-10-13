///////////////////////////////////////////////////////////////////////////////
//
// File: AssociativeTraits.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpressionPolicyFwd.hpp>
#include <boost/utility/enable_if.hpp>

namespace Nektar
{

    template<typename FirstType,
                template <typename, typename> class FirstOp, 
                typename SecondType,
                template <typename, typename> class SecondOp,
                typename ThirdType>
    class AssociativeTraits
    {
        public:
            static const bool IsStrictlyAssociative = false;
            static const bool IsAssociativeWithOpChange = false;
            static const bool IsAssociative = IsStrictlyAssociative || IsAssociativeWithOpChange;
    };
    
    template<typename FirstType, typename SecondType, typename ThirdType>
    class AssociativeTraits<FirstType, AddOp, SecondType, AddOp, ThirdType>
    {
        public:
            static const bool IsStrictlyAssociative = true;
            static const bool IsAssociativeWithOpChange = false;
            static const bool IsAssociative = IsStrictlyAssociative || IsAssociativeWithOpChange;
    };
    
    template<typename FirstType, typename SecondType, typename ThirdType>
    class AssociativeTraits<FirstType, AddOp, SecondType, SubtractOp, ThirdType>
    {
        public:
            static const bool IsStrictlyAssociative = true;
            static const bool IsAssociativeWithOpChange = false;
            static const bool IsAssociative = IsStrictlyAssociative || IsAssociativeWithOpChange;
    };
    
    template<typename FirstType, typename SecondType, typename ThirdType>
    class AssociativeTraits<FirstType, SubtractOp, SecondType, AddOp, ThirdType>
    {
        public:
            static const bool IsStrictlyAssociative = false;
            static const bool IsAssociativeWithOpChange = true;
            static const bool IsAssociative = IsStrictlyAssociative || IsAssociativeWithOpChange;
            typedef SubtractOp<SecondType, ThirdType> OpChangeType;
    };

    template<typename FirstType, typename SecondType, typename ThirdType>
    class AssociativeTraits<FirstType, SubtractOp, SecondType, SubtractOp, ThirdType>
    {
        public:
            static const bool IsStrictlyAssociative = false;
            static const bool IsAssociativeWithOpChange = true;
            static const bool IsAssociative = IsStrictlyAssociative || IsAssociativeWithOpChange;
            typedef AddOp<SecondType, ThirdType> OpChangeType;
    };
    
    template<typename FirstType, typename SecondType, typename ThirdType>
    class AssociativeTraits<FirstType, MultiplyOp, SecondType, MultiplyOp, ThirdType>
    {
        public:
            static const bool IsStrictlyAssociative = true;
            static const bool IsAssociativeWithOpChange = false;
            static const bool IsAssociative = IsStrictlyAssociative || IsAssociativeWithOpChange;
    };
    
    template<typename FirstType, typename SecondType, typename ThirdType>
    class AssociativeTraits<FirstType, MultiplyOp, SecondType, DivideOp, ThirdType>
    {
        public:
            static const bool IsStrictlyAssociative = true;
            static const bool IsAssociativeWithOpChange = false;
            static const bool IsAssociative = IsStrictlyAssociative || IsAssociativeWithOpChange;
    };
    
    template<typename FirstType, typename SecondType, typename ThirdType>
    class AssociativeTraits<FirstType, DivideOp, SecondType, MultiplyOp, ThirdType>
    {
        public:
            static const bool IsStrictlyAssociative = false;
            static const bool IsAssociativeWithOpChange = true;
            static const bool IsAssociative = IsStrictlyAssociative || IsAssociativeWithOpChange;
            typedef DivideOp<SecondType, ThirdType> OpChangeType;
    };

    template<typename AssociativeTraitsType, typename DefaultType, typename enabled = void>
    class ChooseOpChangeType
    {
        public:
            typedef DefaultType Type;
    };

    template<typename AssociativeTraitsType, typename DefaultType>
    class ChooseOpChangeType<AssociativeTraitsType, DefaultType, 
        typename boost::disable_if<boost::is_same<typename AssociativeTraitsType::OpChangeType, DefaultType> >::type >
    {
        public:
            typedef typename AssociativeTraitsType::OpChangeType Type;
    };
}

#endif // NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRAITS_HPP
