///////////////////////////////////////////////////////////////////////////////
//
// File: NullOp.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_NULL_OP_HPP
#define NEKTAR_LIB_UTILITIES_NULL_OP_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <LibUtilities/ExpressionTemplates/UnaryExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>

#include <boost/call_traits.hpp>

namespace Nektar
{
    class NullOp
    {
        public:
    };
    
    class ConstantNullOp
    {
    };
    
    template<typename DataType>
    class UnaryNullOp
    {
        public:
    };
    
    template<typename L, typename R>
    class BinaryNullOpTraits
    {
        public:
            static const bool HasOpEqual = true;
            static const bool HasOpLeftEqual = false;
    };
    
    template<typename LhsType, typename RhsType>
    class BinaryNullOp
    {
        public:
            typedef BinaryNullOpTraits<LhsType, RhsType> TraitsType;
            
            static void ApplyEqual(Accumulator<LhsType>& result,
                                    typename boost::call_traits<RhsType>::const_reference rhs)
            {
                *result = rhs;
            }
            
            template<typename L, typename R>
            class Rebind
            {
                public:
                    typedef BinaryNullOp<L, R> type;
            };
    };

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif // NEKTAR_LIB_UTILITIES_NULL_OP_HPP

/**
    $Log: NullOp.hpp,v $
    Revision 1.6  2007/08/16 02:14:21  bnelson
    Moved expression templates to the Nektar namespace.

    Revision 1.5  2007/01/30 23:37:17  bnelson
    *** empty log message ***

    Revision 1.4  2007/01/29 01:34:05  bnelson
    Updates to compile in windows.

    Revision 1.3  2007/01/16 17:37:56  bnelson
    Wrapped everything with #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

    Revision 1.2  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.1  2006/11/06 17:07:19  bnelson
    Continued work on creating temporaries as needed when sub-expression types don't match the type of the accumulator.


 **/
