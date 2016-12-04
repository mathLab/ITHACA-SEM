///////////////////////////////////////////////////////////////////////////////
//
// File: NekPoint.hpp
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
// Description: Generic N-Dimensional Point.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_POINT_HPP
#define NEKTAR_LIB_UTILITIES_NEK_POINT_HPP

#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/Concepts.hpp>
#include <LibUtilities/LinearAlgebra/Space.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <boost/concept_check.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/call_traits.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/assert.hpp>

#include <math.h>
#include <functional>
#include <algorithm>

namespace Nektar
{
    template<typename data_type>
    class NekPoint
    {
        public:
            typedef data_type DataType;
            typedef ThreeD dim;

        public:
            NekPoint()
            {
                // This may be suboptimal if DataType isn't numeric.
                // If we use them then maybe we could look at an enable_if
                // template to choose a better constructor.
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    // If you get a compile error pointing you here then
                    // the DataType being stored in the point doesn't have an
                    // accessible operator= or default copy constructor.
                    m_data[i] = DataType();
                }
            }

            NekPoint(const std::string& pointValues)
            {
                bool result = fromString(pointValues, *this);
                ASSERTL0(result, "Unable to create a point from a string.");
            }

            NekPoint(typename boost::call_traits<DataType>::param_type x,
                        typename boost::call_traits<DataType>::param_type y,
                        typename boost::call_traits<DataType>::param_type z)
            {
                m_data[0] = x;
                m_data[1] = y;
                m_data[2] = z;
            }

            explicit NekPoint(typename boost::call_traits<DataType>::const_reference a)
            {
                boost::function_requires< Nektar::AssignableConcept<DataType> >();

                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    m_data[i] = a;
                }
            }


            NekPoint(const NekPoint<DataType>& rhs)
            {
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    m_data[i] = rhs.m_data[i];
                }
            }

            ~NekPoint()
            {
            }

 
            NekPoint<DataType>& operator=(const NekPoint<DataType>& rhs)
            {
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    m_data[i] = rhs.m_data[i];
                }
                return *this;
            }

            /// \brief Returns the number of dimensions for the point.
            static unsigned int dimension() { return dim::Value; }

            /// \brief Returns i^{th} element.
            /// \param i The element to return.
            /// \pre i < dim
            /// \return A reference to the i^{th} element.
            ///
            /// Retrieves the i^{th} element.  Since it returns a reference you may
            /// assign a new value (i.e., p(2) = 3.2;)
            ///
            /// This operator performs range checking.
            typename boost::call_traits<DataType>::reference operator()(unsigned int i)
            {
                ASSERTL0( (i>=0) && (i<dim::Value), "Invalid access to NekPoint data via parenthesis operator: index out of range");
                return m_data[i];
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int i) const
            {
                ASSERTL0( (i>=0) && (i<dim::Value), "Invalid access to NekPoint data via parenthesis operator: index out of range");
                return m_data[i];
            }

            typename boost::call_traits<DataType>::reference operator[](unsigned int i)
            {
                return m_data[i];
            }

            typename boost::call_traits<DataType>::const_reference operator[](unsigned int i) const
            {
                return m_data[i];
            }

            typename boost::call_traits<DataType>::const_reference x() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 1);
                return m_data[0];
            }

            typename boost::call_traits<DataType>::const_reference y() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 2);
                return (*this)[1];
            }

            typename boost::call_traits<DataType>::const_reference z() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 3);
                return (*this)[2];
            }

            typename boost::call_traits<DataType>::const_reference a() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 1);
                return m_data[0];
            }

            typename boost::call_traits<DataType>::const_reference b() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 2);
                return (*this)[1];
            }

            typename boost::call_traits<DataType>::const_reference c() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 3);
                return (*this)[2];
            }

            typename boost::call_traits<DataType>::const_reference r() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 1);
                return m_data[0];
            }

            typename boost::call_traits<DataType>::const_reference s() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 2);
                return (*this)[1];
            }

            typename boost::call_traits<DataType>::const_reference t() const
            {
                BOOST_STATIC_ASSERT(dim::Value >= 3);
                return (*this)[2];
            }

            void SetX(typename boost::call_traits<DataType>::const_reference val)
            {
                BOOST_STATIC_ASSERT(dim::Value >= 1);
                m_data[0] = val;
            }

            void SetY(typename boost::call_traits<DataType>::const_reference val)
            {
                BOOST_STATIC_ASSERT(dim::Value >= 2);
                m_data[1] = val;
            }

            void SetZ(typename boost::call_traits<DataType>::const_reference val)
            {
                BOOST_STATIC_ASSERT(dim::Value >= 2);
                m_data[2] = val;
            }

            typename boost::call_traits<DataType>::reference x()
            {
                BOOST_STATIC_ASSERT(dim::Value >= 1);
                return (*this)(0);
            }

            typename boost::call_traits<DataType>::reference y()
            {
                BOOST_STATIC_ASSERT(dim::Value >= 2);
                return (*this)(1);
            }

            typename boost::call_traits<DataType>::reference z()
            {
                BOOST_STATIC_ASSERT(dim::Value >= 3);
                return (*this)(2);
            }

            const DataType* GetPtr() const
            {
                return &m_data[0];
            }

            bool operator==(const NekPoint<DataType>& rhs) const
            {
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    // If you get a compile error here then you have to
                    // add a != operator to the DataType class.
                    if( m_data[i] != rhs.m_data[i] )
                    {
                        return false;
                    }
                }
                return true;
            }

            bool operator!=(const NekPoint<DataType>& rhs) const
            {
                return !(*this == rhs);
            }

            /// Arithmetic Routines

            // Unitary operators
            void negate()
            {
                for(int i=0; i < dim::Value; ++i)
                {
                    (*this)[i] = -(*this)[i];
                }
            }

            NekPoint<DataType> operator-() const
            {
                NekPoint<DataType> result(*this);
                result.negate();
                return result;
            }


            NekPoint<DataType>& operator+=(const NekPoint<DataType>& rhs)
            {
                for(unsigned int i=0; i < dim::Value; ++i)
                {
                    m_data[i] += rhs.m_data[i];
                }
                return *this;
            }

            NekPoint<DataType>& operator+=(typename boost::call_traits<DataType>::param_type rhs)
            {
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    m_data[i] += rhs;
                }
                return *this;
            }

            NekPoint<DataType>& operator-=(const NekPoint<DataType>& rhs)
            {
                for(unsigned int i=0; i < dim::Value; ++i)
                {
                    m_data[i] -= rhs.m_data[i];
                }
                return *this;
            }


            NekPoint<DataType>& operator-=(typename boost::call_traits<DataType>::param_type rhs)
            {
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    m_data[i] -= rhs;
                }
                return *this;
            }

            NekPoint<DataType>& operator*=(typename boost::call_traits<DataType>::param_type rhs)
            {
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    m_data[i] *= rhs;
                }
                return *this;
            }

            NekPoint<DataType>& operator/=(typename boost::call_traits<DataType>::param_type rhs)
            {
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    m_data[i] /= rhs;
                }
                return *this;
            }

            std::string AsString() const
            {
                std::string result = "(";
                for(unsigned int i = 0; i < dim::Value; ++i)
                {
                    result += boost::lexical_cast<std::string>(m_data[i]);
                    if( i < dim::Value-1 )
                    {
                        result += ", ";
                    }
                }
                result += ")";
                return result;
            }

        private:
            DataType m_data[dim::Value];
    };

    // Operators for expression templates
    template<typename DataType>
    void negate(NekPoint<DataType>& rhs)
    {
        rhs.negate();
    }
    
    template<typename DataType>
    NekPoint<DataType>
    operator+(const NekPoint<DataType>& lhs, const NekPoint<DataType>& rhs)
    {
        NekPoint<DataType> result(lhs);
        result += rhs;
        return result;
    }

    template<typename DataType>
    NekPoint<DataType>
    operator+(typename boost::call_traits<DataType>::const_reference lhs, const NekPoint<DataType>& rhs)
    {
        NekPoint<DataType> result(rhs);
        result += lhs;
        return result;
    }

    template<typename DataType>
    NekPoint<DataType>
    operator+(const NekPoint<DataType>& lhs, typename boost::call_traits<DataType>::const_reference rhs)
    {
        NekPoint<DataType> result(lhs);
        result += rhs;
        return result;
    }

    template<typename DataType>
    NekPoint<DataType>
    operator-(const NekPoint<DataType>& lhs, const NekPoint<DataType>& rhs)
    {
        NekPoint<DataType> result(lhs);
        result -= rhs;
        return result;
    }

    template<typename DataType>
    NekPoint<DataType>
    operator-(typename boost::call_traits<DataType>::const_reference lhs, const NekPoint<DataType>& rhs)
    {
        NekPoint<DataType> result(-rhs);
        result += lhs;
        return result;
    }

    template<typename DataType>
    NekPoint<DataType>
    operator-(const NekPoint<DataType>& lhs, typename boost::call_traits<DataType>::const_reference rhs)
    {
        NekPoint<DataType> result(lhs);
        result -= rhs;
        return result;
    }

    template<typename DataType, typename dim, typename space, typename ScalarType>
    NekPoint<DataType>
    operator*(const ScalarType& lhs, const NekPoint<DataType>& rhs)
    {
        NekPoint<DataType> result(rhs);
        result *= lhs;
        return result;
    }

    template<typename DataType, typename dim, typename space, typename ScalarType>
    NekPoint<DataType>
    operator*(const NekPoint<DataType>& lhs, const ScalarType& rhs)
    {
        NekPoint<DataType> result(lhs);
        result *= rhs;
        return result;
    }

    template<typename DataType>
    NekPoint<DataType>
    operator/(const NekPoint<DataType>& lhs, typename boost::call_traits<DataType>::param_type rhs)
    {
        NekPoint<DataType> result(lhs);
        result /= rhs;
        return result;
    }

    template<typename DataType>
    typename boost::call_traits<DataType>::value_type distanceBetween(const NekPoint<DataType>& lhs,
                            const NekPoint<DataType>& rhs)
    {
        DataType result = 0.0;
        for(unsigned int i = 0; i < 3; ++i)
        {
            DataType temp = lhs[i] - rhs[i];
            result += temp*temp;
        }
        return sqrt(result);
    }

    template<typename DataType>
    bool fromString(const std::string& str, NekPoint<DataType>& result)
    {
        try
        {
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            boost::char_separator<char> sep("(<,>) ");
            tokenizer tokens(str, sep);
            unsigned int i = 0;
            for(tokenizer::iterator iter = tokens.begin(); iter != tokens.end(); ++iter)
            {
                result[i] = boost::lexical_cast<DataType>(*iter);
                ++i;
            }

            return i == 3;
        }
        catch(boost::bad_lexical_cast&)
        {
            return false;
        }
    }

    template<typename DataType>
    std::ostream& operator<<(std::ostream& os, const NekPoint<DataType>& p)
    {
        os << p.AsString();
        return os;
    }

    //template<typename DataType>
    //NekPoint<DataType> operator+(const NekPoint<DataType>& lhs, const NekPoint<DataType>& rhs)
    //{
    //    NekPoint<DataType> result(lhs);
    //    result += rhs;
    //    return result;
    //}

    //template<typename DataType>
    //NekPoint<DataType,dim,space> operator+(const NekPoint<DataType,dim,space>& P, const NekVector<DataType>& V){
    //ASSERTL0(P.dimension() == V.size(),"NekPoint::operator", "Point and Vector dimensions do not match");

    //NekPoint<DataType,dim> temp(P);
    //for(int i=0;i<P.dimension();++i)
    //temp(i) = temp(i) + V(i);
    //return temp;
    //}
//
//   template<typename DataType, int dim>
//   NekPoint<DataType,dim> operator+(const NekVector<DataType>& V, const NekPoint<DataType,dim>& P){
//     ASSERTL0(dim == V.size(),"NekPoint::operator", "Point and Vector dimensions do not match");
//
//     NekPoint<DataType,dim> temp(P);
//     for(int i=0;i<dim;++i)
    //temp(i) += V(i);
//     return temp;
//   }
//
//   template<typename DataType, int dim>

//   NekPoint<DataType,dim> operator-(const NekPoint<DataType,dim>& P, const NekVector<DataType>& V){
//     ASSERTL0(dim == V.size(),"NekPoint::operator", "Point and Vector dimensions do not match");
//
//     NekPoint<DataType,dim> temp(P);
//     for(int i=0;i<dim;++i)
    //temp(i) -= V(i);
//     return temp;
//   }
}

#endif // NEKTAR_LIB_UTILITIES_NEK_POINT_HPP

/**
    $Log: NekPoint.hpp,v $
    Revision 1.18  2008/03/23 16:33:01  bnelson
    *** empty log message ***

    Revision 1.17  2008/03/03 02:28:39  bnelson
    Changed OneD, TwoD, and ThreeD to classes instead of enums to support type parameters in NekVector instead of unsigned int for the dimensions.

    Added a new NekVector<DataType> to allow wrapping of ConstArrays.

    Revision 1.16  2008/01/03 04:16:09  bnelson
    Changed method name in the expression library from Apply to Evaluate.

    Revision 1.15  2007/08/16 02:09:56  bnelson
    Moved expression templates to the Nektar namespace.

    Revision 1.14  2007/01/23 03:12:50  jfrazier
    Added more conditional compilation directives for expression templates.

    Revision 1.13  2006/09/30 15:18:37  bnelson
    no message

    Revision 1.12  2006/09/16 23:53:36  bnelson
    Modified the negation operation to reflect changes in the unary expression templates.

    Revision 1.11  2006/09/15 02:18:39  bnelson
    Fixed a problem with operator*

    Revision 1.10  2006/09/14 02:06:16  bnelson
    Fixed gcc compiler errors.

    Revision 1.9  2006/09/11 03:26:26  bnelson
    Updated to use new policy based expression templates.

    Revision 1.8  2006/09/10 20:40:24  bnelson
    Changed DataType to data_type

    Revision 1.7  2006/09/08 03:37:18  bnelson
    Fixed an ambiguous case with the templated copy constructor.

    Revision 1.6  2006/08/28 02:40:21  bnelson
    *** empty log message ***

    Revision 1.5  2006/08/25 01:27:04  bnelson
    Added construction and assignment from expressions.

    Added the unary negation expression template.

    Revision 1.4  2006/08/14 02:40:24  bnelson
    no message

    Revision 1.3  2006/08/14 02:29:49  bnelson
    Updated points, vectors, and matrix classes to work with ElVis.  Added a variety of methods to all of these classes.

    Revision 1.2  2006/06/01 13:44:28  kirby
    *** empty log message ***

    Revision 1.1  2006/06/01 09:12:42  kirby
    *** empty log message ***

    Revision 1.1  2006/05/04 18:57:44  kirby
    *** empty log message ***

    Revision 1.5  2006/04/11 02:01:32  bnelson
    Fleshed out the interface.

    Revision 1.4  2006/04/06 03:45:22  bnelson
    Added a third template parameter space.

    Revision 1.3  2006/03/25 00:52:43  jfrazier
    Minor formatting stuff to correct indenting.

    Revision 1.2  2006/01/31 13:51:13  bnelson
    Updated for new configure.

    Revision 1.1  2006/01/10 14:50:31  bnelson
    Initial Revision.

**/

