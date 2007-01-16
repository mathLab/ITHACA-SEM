///////////////////////////////////////////////////////////////////////////////
//
// File: ExpressionMetadata.hpp
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

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_METADATA_H
#define NEKTAR_LIB_UTILITIES_EXPRESSION_METADATA_H

#include <boost/utility.hpp>
#include <boost/type_traits.hpp>

#include <LibUtilities/ExpressionTemplates/ExpressionTraits.hpp>

namespace Nektar
{
    namespace expt
    {
        template<typename DataType>
        class DefaultExpressionMetadata
        {
            public:
                DefaultExpressionMetadata(typename boost::call_traits<DataType>::const_reference)
                {
                }

                template<typename T>
                DefaultExpressionMetadata(const DefaultExpressionMetadata<T>& rhs)
                {
                }
                
                template<typename L, typename R>
                DefaultExpressionMetadata(const DefaultExpressionMetadata<L>& lhs, 
                                          const DefaultExpressionMetadata<R>& rhs)
                {
                }
                
                DefaultExpressionMetadata(const DefaultExpressionMetadata<DataType>& rhs)
                {
                }
                
                DefaultExpressionMetadata()
                {
                }
                

            private:
        };

       
        template<typename TraitsType, typename enabled = void>
        class ExpressionMetadataChooser
        {
            public:
                typedef DefaultExpressionMetadata<typename TraitsType::result_type> MetadataType;
        };
        
        template<typename TraitsType>
        class ExpressionMetadataChooser<TraitsType,
                                        typename boost::enable_if<boost::is_same<typename TraitsType::MetadataType, typename TraitsType::MetadataType> >::type >
                                        //typename boost::enable_if<boost::is_same<typename TraitsType::MetadataType, typename TraitsType::MetadataType> >::type >
        {
            public:
                typedef typename TraitsType::MetadataType MetadataType;
        };
        
    }
}

#endif // NEKTAR_LIB_UTILITIES_EXPRESSION_METADATA_H
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
/**
    $Log: ExpressionMetadata.hpp,v $
    Revision 1.3  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.2  2006/11/12 17:58:52  bnelson
    *** empty log message ***

    Revision 1.1  2006/09/14 02:09:00  bnelson
    Fixed gcc compile errors

**/

