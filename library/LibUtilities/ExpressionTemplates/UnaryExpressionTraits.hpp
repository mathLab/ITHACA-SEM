///////////////////////////////////////////////////////////////////////////////
//
// File: UnaryExpressionTraits.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_UNARY_EXPRESSION_TRAITS_HPP
#define NEKTAR_LIB_UTILITIES_UNARY_EXPRESSION_TRAITS_HPP

namespace Nektar
{
    /// \brief Defines the abstract UnaryExpressionTraits class. 
    ///
    /// UnaryExpressionTraits are used when using unary expressions with 
    /// expression templates.  They specify the return type of a unary 
    /// expression.
    ///
    /// Every UnaryExpressionTemplate will be looking for a different typedef
    /// inside this class.  Consult the unary expression you are using for details.
    template<typename DataType>
    class UnaryExpressionTraits
    {
        public:
            typedef DataType NegationType;
    };
    
    template<typename DataType>
    class NegateTraits;
}

#endif // NEKTAR_LIB_UTILITIES_UNARY_EXPRESSION_TRAITS_HPP

/**
    $Log: UnaryExpressionTraits.hpp,v $
    Revision 1.3  2007/01/16 17:37:56  bnelson
    Wrapped everything with #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

    Revision 1.2  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.1  2006/08/25 01:33:48  bnelson
    no message

**/
