///////////////////////////////////////////////////////////////////////////////
//
// File: NekVectorTypeTraits.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_TYPE_TRAITS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_TYPE_TRAITS_HPP

#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>

#include <boost/type_traits.hpp>
#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace Nektar
{
    template<typename DataType, typename enabled = void>
    class NekVectorTypeTraits;

    template<typename DataType>
    class NekVectorTypeTraits<DataType, typename boost::enable_if<boost::is_floating_point<DataType> >::type >
    {
        public:
            static 
            typename boost::call_traits<DataType>::value_type 
            abs(typename boost::call_traits<DataType>::const_reference v)
            {
                return fabs(v);
            }
    };
    
    template<typename DataType>
    class NekVectorTypeTraits<DataType, typename boost::enable_if<boost::is_integral<DataType> >::type >
    {
        public:
            static 
            typename boost::call_traits<DataType>::value_type 
            abs(typename boost::call_traits<DataType>::const_reference v)
            {
                return abs(v);
            }
    };
    
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_TYPE_TRAITS_HPP
