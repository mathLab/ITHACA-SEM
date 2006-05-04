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

#include <LibUtilities/ExpressionTemplates.hpp>
#include <LibUtilities/ArithmeticConcepts.hpp>
#include <LibUtilities/ErrorUtil.hpp>
#include <LibUtilities/Concepts.hpp>

#include <boost/concept_check.hpp>

#include <functional>
#include <algorithm>

namespace Nektar
{
    namespace LibUtilities
    {

        template<typename DataType, unsigned int dim, unsigned int space = 0>
        class NekPoint
        {
            public:
                NekPoint()
                {
                    // This may be suboptimal if DataType isn't numeric.
                    // If we use them then maybe we could look at an enable_if
                    // template to choose a better constructor.
                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        // If you get a compile error pointing you here then
                        // the DataType being stored in the point doesn't have an
                        // accessible operator= or default copy constructor.
                        m_data[i] = DataType();
                    }
                }

                explicit NekPoint(typename boost::call_traits<DataType>::param_type a)
                {
                    boost::function_requires< Nektar::AssignableConcept<DataType> >();

                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        m_data[i] = a;
                    }
                }

                NekPoint(const NekPoint<DataType, dim, space>& rhs)
                {
                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        m_data[i] = rhs.m_data[i];
                    }
                }

                ~NekPoint()
                {
                }

                NekPoint<DataType, dim, space>& operator=(const NekPoint<DataType, dim, space>& rhs)
                {
                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        m_data[i] = rhs.m_data[i];
                    }
                    return *this;
                }

                /// \brief Returns the number of dimensions for the point.
                static unsigned int dimension() { return dim; }

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
                    ASSERTL1( (i>=0) && (i<dim), "Invalid access to _data via parenthesis operator");
                    return m_data[i];
                }

                typename boost::call_traits<DataType>::const_reference operator()(unsigned int i) const
                {
                    ASSERTL1( (i>=0) && (i<dim), "Invalid access to _data via parenthesis operator");
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
                    BOOST_STATIC_ASSERT(dim >= 1);
                    return m_data[0];
                }

                typename boost::call_traits<DataType>::const_reference y() const
                {
                    BOOST_STATIC_ASSERT(dim >= 2);
                    return (*this)(1);
                }

                typename boost::call_traits<DataType>::const_reference z() const
                {
                    BOOST_STATIC_ASSERT(dim >= 3);
                    return (*this)(2);
                }

                typename boost::call_traits<DataType>::reference x()
                {
                    BOOST_STATIC_ASSERT(dim >= 1);
                    return (*this)(0);
                }

                typename boost::call_traits<DataType>::reference y()
                {
                    BOOST_STATIC_ASSERT(dim >= 2);
                    return (*this)(1);
                }

                typename boost::call_traits<DataType>::reference z()
                {
                    BOOST_STATIC_ASSERT(dim >= 3);
                    return (*this)(2);
                }

                const DataType* getPtr() const
                {
                    return &m_data[0];
                }

                bool operator==(const NekPoint<DataType, dim, space>& rhs)
                {
                    for(unsigned int i = 0; i < dim; ++i)
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

                /// Arithmetic Routines

                // Unitary operators
                NekPoint<DataType, dim, space> operator-() const
                {
                    NekPoint<DataType, dim, space> temp(*this);
                    for(int i=0; i < dim; ++i)
                    {
                        temp(i) = -temp(i);
                    }
                    return temp;
                }

                NekPoint<DataType, dim, space>& operator+=(const NekPoint<DataType, dim, space>& rhs)
                {
                    for(unsigned int i=0; i < dim; ++i)
                    {
                        m_data[i] += rhs.m_data[i];
                    }
                    return *this;
                }

                NekPoint<DataType, dim, space>& operator+=(typename boost::call_traits<DataType>::param_type rhs)
                {
                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        m_data[i] += rhs;
                    }
                    return *this;
                }

                NekPoint<DataType, dim, space>& operator-=(const NekPoint<DataType, dim, space>& rhs)
                {
                    for(unsigned int i=0; i < dim; ++i)
                    {
                        m_data[i] -= rhs.m_data[i];
                    }
                    return *this;
                }


                NekPoint<DataType, dim, space>& operator-=(typename boost::call_traits<DataType>::param_type rhs)
                {
                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        m_data[i] -= rhs;
                    }
                    return *this;
                }

                NekPoint<DataType, dim, space>& operator*=(typename boost::call_traits<DataType>::param_type rhs)
                {
                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        m_data[i] *= rhs;
                    }
                    return *this;
                }

                NekPoint<DataType, dim, space>& operator/=(typename boost::call_traits<DataType>::param_type rhs)
                {
                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        m_data[i] /= rhs;
                    }
                    return *this;
                }

            private:
                DataType m_data[dim];
        };

        template<typename DataType, unsigned int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator+(const NekPoint<DataType, dim, space>& lhs, const NekPoint<DataType, dim>& rhs)
        {
            NekPoint<DataType, dim, space> result(lhs);
            result += rhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator+(typename boost::call_traits<DataType>::param_type lhs, const NekPoint<DataType, dim, space>& rhs)
        {
            NekPoint<DataType, dim, space> result(rhs);
            result += lhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator+(const NekPoint<DataType, dim, space>& lhs, typename boost::call_traits<DataType>::param_type rhs)
        {
            NekPoint<DataType, dim, space> result(lhs);
            result += rhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator-(const NekPoint<DataType, dim, space>& lhs, const NekPoint<DataType, dim>& rhs)
        {
            NekPoint<DataType, dim, space> result(lhs);
            result -= rhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator-(typename boost::call_traits<DataType>::param_type lhs, const NekPoint<DataType, dim, space>& rhs)
        {
            NekPoint<DataType, dim, space> result(rhs);
            result -= lhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator-(const NekPoint<DataType, dim, space>& lhs, typename boost::call_traits<DataType>::param_type rhs)
        {
            NekPoint<DataType, dim, space> result(lhs);
            result -= rhs;
            return result;
        }

        template<typename DataType, int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator*(typename boost::call_traits<DataType>::param_type lhs, const NekPoint<DataType, dim, space>& rhs)
        {
            NekPoint<DataType, dim, space> result(rhs);
            result *= lhs;
            return result;
        }

        template<typename DataType, int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator*(const NekPoint<DataType, dim, space>& lhs, typename boost::call_traits<DataType>::param_type rhs)
        {
            NekPoint<DataType, dim, space> result(lhs);
            result *= rhs;
            return result;
        }

        template<typename DataType, int dim, unsigned int space>
        NekPoint<DataType, dim, space>
        operator/(const NekPoint<DataType, dim, space>& lhs, typename boost::call_traits<DataType>::param_type rhs)
        {
            NekPoint<DataType, dim, space> result(lhs);
            result /= rhs;
            return result;
        }

    //   template<typename DataType, int dim>
    //   NekPoint<DataType,dim> operator+(const NekPoint<DataType,dim>& P, const NekVector<DataType>& V){
    //     ASSERTL0(P.dimension() == V.size(),"NekPoint::operator", "Point and Vector dimensions do not match");
        //
    //     NekPoint<DataType,dim> temp(P);
    //     for(int i=0;i<P.dimension();++i)
    //      temp(i) = temp(i) + V(i);
    //     return temp;
    //   }
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
}

#endif // NEKTAR_LIB_UTILITIES_NEK_POINT_HPP

/**
    $Log: NekPoint.hpp,v $
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

