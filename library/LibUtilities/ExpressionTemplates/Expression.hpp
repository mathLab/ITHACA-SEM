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

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/mpl/and.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/assert.hpp>

#include <LibUtilities/ExpressionTemplates/Accumulator.hpp>

#include <algorithm>

namespace Nektar
{
        /// \brief An expression is an arbitrary combination of operations acting on arbitrary amountsw of data to 
        /// produce a result.  The expressions stores the operations and data but does not evaluate it until requested.
        /// 
        /// \param PolicyType The type of expression.
        template<typename PolicyType>
        class Expression
        {
            public:
                // The type of the data which will be the result of evaluating the expression.
                typedef typename PolicyType::ResultType ResultType;

                // The type of the unevaluated subexpressions or data.
                typedef typename PolicyType::DataType DataType;

                typedef typename PolicyType::MetadataType MetadataType;
                
                typedef typename boost::remove_const<typename boost::remove_reference<DataType>::type >::type BaseDataType;
            public:
                /// This constructor allows temporaries.  I can't find a way to prevent this 
                /// from happening without also preventing constant parameters, so this will
                /// need to be a "don't do that" type warning.
                template<typename T>
                explicit Expression(const T& data) :
                    m_data(data),
                    m_metadata()
                {
                    PolicyType::InitializeMetadata(m_data, m_metadata);
                }

                explicit Expression(BaseDataType& data) :
                    m_data(data),
                    m_metadata()
                {
                    //using boost::mpl::and_;
                    //using boost::is_reference;
                    //using boost::is_same;
                    //using boost::is_const;

                    //BOOST_MPL_ASSERT( ( is_const<BaseDataType> ) );

                    PolicyType::InitializeMetadata(m_data, m_metadata);
                }

                //template<typename T>
                //explicit Expression(const T& data) :
                //    m_data(data),
                //    m_metadata()
                //{
                //    //BOOST_STATIC_ASSERT(0);
                //    using boost::mpl::and_;
                //    using boost::is_reference;
                //    using boost::is_same;

                //    BOOST_MPL_ASSERT_NOT( (and_<is_reference<DataType>, is_same<BaseDataType, T> >) );
                //    PolicyType::InitializeMetadata(m_data, m_metadata);
                //}

                Expression(const Expression<PolicyType>& rhs) :
                    m_data(rhs.m_data),
                    m_metadata(rhs.m_metadata)
                {
                }

                ~Expression() {}

                void Apply(typename boost::call_traits<ResultType>::reference result) const
                {
                   if( !PolicyType::ContainsReference(result, m_data) )
                   {
                        Accumulator<ResultType> accum(result);
                        PolicyType::Apply(accum, m_data);
                   }
                   else
                   {
                       ResultType temp;
                       Accumulator<ResultType> accum(result);
                       PolicyType::Apply(accum, m_data);
                       result = temp;
                   }
                        
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

                typename boost::call_traits<DataType>::reference operator*()
                {
                    return m_data;
                }
                
                typename boost::call_traits<DataType>::const_reference operator*() const
                {
                    return m_data;
                }
                
            private:
                Expression<PolicyType>& operator=(const Expression<PolicyType>& rhs);

                DataType m_data;
                MetadataType m_metadata;
        };
        
        template<typename DataType, typename PolicyType>
        void Assign(DataType& result, const Expression<PolicyType>& expr)
        {
            expr.Apply(result);
        }

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif // NEKTAR_LIB_UTILITIES_EXPRESSION_HPP
/**
    $Log: Expression.hpp,v $
    Revision 1.15  2007/12/19 05:09:21  bnelson
    First pass at detecting aliasing.  Still need to test performance implications.

    Revision 1.14  2007/11/13 18:07:25  bnelson
    Added Assign helper function for those classes the user can't modfiy.

    Revision 1.13  2007/10/04 03:48:54  bnelson
    *** empty log message ***

    Revision 1.12  2007/10/03 02:57:39  bnelson
    Removed the restriction on passing temporaries to expressions.

    Revision 1.11  2007/08/16 02:14:21  bnelson
    Moved expression templates to the Nektar namespace.

    Revision 1.10  2007/07/26 00:07:36  bnelson
    Fixed linux compiler errors.

    Revision 1.8  2007/07/20 00:36:05  bnelson
    *** empty log message ***

    Revision 1.7  2007/01/30 23:37:16  bnelson
    *** empty log message ***

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
