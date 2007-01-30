///////////////////////////////////////////////////////////////////////////////
//
// File: Expression.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_HPP

#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <LibUtilities/ExpressionTemplates/Accumulator.hpp>

#include <algorithm>

namespace Nektar
{
    namespace expt
    {
//         template<typename ExpressionPolicy>
//         class Expression;
        template<typename PolicyType>
        class Expression
        {
            public:
                typedef typename PolicyType::ResultType ResultType;
                typedef typename PolicyType::MetadataType MetadataType;
                
            public:
                explicit Expression(typename boost::call_traits<typename PolicyType::DataType>::const_reference data) :
                    m_data(data),
                    m_metadata()
                {
                    PolicyType::InitializeMetadata(m_data, m_metadata);
                }

                Expression(const Expression<PolicyType>& rhs) :
                    m_data(rhs.m_data),
                    m_metadata(rhs.m_metadata)
                {
                }

                ~Expression() {}

                void Apply(typename boost::call_traits<ResultType>::reference result) const
                {
                    Accumulator<ResultType> accum(result);
                    PolicyType::Apply(accum, m_data);
                }

                template<template <typename, typename> class ParentOpType>
                void Apply(Accumulator<ResultType>& accum) const
                {
                    PolicyType::template Apply<ParentOpType>(accum, m_data);
                }
                
                template<template <typename> class ParentOpType>
                void Apply(Accumulator<ResultType>& accum) const
                {
                    PolicyType::template Apply<ParentOpType>(accum, m_data);
                }

//                 template<template <typename, typename> class ParentOpType>
//                 void ApplyEqual(Accumulator<ResultType>& accum) const
//                 {
//                     PolicyType::template ApplyEqual<ParentOpType>(accum, m_data);
//                 }
                
                const MetadataType& GetMetadata() const
                {
                    return m_metadata;
                }

                void Print(std::ostream& os) const
                {
                    PolicyType::Print(os, m_data);
                }

                typename boost::call_traits<typename PolicyType::DataType>::reference operator*()
                {
                    return m_data;
                }
                
                typename boost::call_traits<typename PolicyType::DataType>::const_reference operator*() const
                {
                    return m_data;
                }
                
            private:
                Expression<PolicyType>& operator=(const Expression<PolicyType>& rhs);

                typename PolicyType::DataType m_data;
                MetadataType m_metadata;
        };
    }

}

#endif // NEKTAR_LIB_UTILITIES_EXPRESSION_HPP
/**
    $Log: Expression.hpp,v $
    Revision 1.6  2007/01/16 17:37:55  bnelson
    Wrapped everything with #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

    Revision 1.5  2007/01/16 05:29:50  bnelson
    Major improvements for expression templates.

    Revision 1.4  2006/09/14 02:08:59  bnelson
    Fixed gcc compile errors

    Revision 1.3  2006/09/11 03:24:24  bnelson
    Updated expression templates so they are all specializations of an Expression object, using policies to differentiate.

    Revision 1.2  2006/08/28 02:39:53  bnelson
    *** empty log message ***

    Revision 1.1  2006/08/25 01:33:47  bnelson
    no message

**/
