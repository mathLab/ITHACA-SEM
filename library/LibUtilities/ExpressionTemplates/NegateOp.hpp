///////////////////////////////////////////////////////////////////////////////
//
// File: NetageOP.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_NEGATE_EXPRESSION_HPP
#define NEKTAR_LIB_UTILITIES_NEGATE_EXPRESSION_HPP

#include <LibUtilities/ExpressionTemplates/UnaryExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>

#include <boost/call_traits.hpp>

namespace Nektar
{

    namespace expt
    {
        /// \brief An expression to negate an object of ParameterType.
        /// Parameter type is the actual object, not the expression that may lead to it.
        template<typename ParameterType>
        class NegateOp
        {
            public:
                typedef typename UnaryExpressionTraits<ParameterType>::NegationType ResultType;
                typedef ParameterType InputType;
                typedef typename ExpressionMetadataChooser<ParameterType>::MetadataType MetadataType;

                static void Apply(typename boost::call_traits<ParameterType>::reference result)
                {
                    negate(result);
                }

                static MetadataType CreateUnaryMetadata(const MetadataType& rhs)
                {
                    return MetadataType::CreateForNegation(rhs);
                }

            private:
        };
    }
}

#endif // NEKTAR_LIB_UTILITIES_NEGATE_EXPRESSION_HPP

/**
    $Log: NegateOp.hpp,v $
    Revision 1.5  2006/09/14 02:08:59  bnelson
    Fixed gcc compile errors

    Revision 1.4  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

    Revision 1.3  2006/08/28 02:39:53  bnelson
    *** empty log message ***

    Revision 1.2  2006/08/27 02:11:29  bnelson
    Added support for negating an expression.

    Revision 1.1  2006/08/25 01:33:48  bnelson
    no message

**/
