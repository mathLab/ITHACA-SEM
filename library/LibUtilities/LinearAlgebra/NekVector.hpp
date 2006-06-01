///////////////////////////////////////////////////////////////////////////////
//
// File: NekVector.hpp
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
// Description: Generic N-Dimensional Vector.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP
#define NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#include <functional>
#include <algorithm>
#include <math.h>

#include <boost/call_traits.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility/enable_if.hpp>

namespace Nektar
{
    namespace LibUtilities
    {

        // \param DataType The type of data held by each element of the vector.
        // \param dim The number of elements in the vector.  If set to 0, the vector
        //            will have a variable number of elements.
        // \param space The space of the vector.
        template<typename DataType, unsigned int dim = 0, unsigned int space = 0>
        class NekVector
        {
            public:
                NekVector() :
                    m_impl()
                {
                }

                // \brief Creates a vector with given size and initial value.
                //        This constructor is only valid for variable sized vectors.
                NekVector(unsigned int size, typename boost::call_traits<DataType>::param_type a) :
                    m_impl(size, a)
                {
                }

                explicit NekVector(typename boost::call_traits<DataType>::param_type a) :
                    m_impl(a)
                {
                }

                NekVector(const NekVector<DataType, 0, space>& rhs) :
                    m_impl(rhs.m_impl)
                {
                }

                NekVector(unsigned int size, const DataType* const ptr) :
                    m_impl(size, ptr)
                {
                }

                ~NekVector()
                {
                }

                NekVector<DataType, dim, space>& operator=(const NekVector<DataType, dim, space>& rhs)
                {
                    m_impl = rhs.m_impl;

                    return *this;
                }

                /// \brief Returns the number of dimensions for the point.
                unsigned int dimension() const
                {
                    return m_impl.dimension();
                }

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
                    ASSERTL1((i >= 0) && (i < dimension()), "Invalid access to m_data via parenthesis operator");
                    return m_impl.data[i];
                }

                typename boost::call_traits<DataType>::const_reference operator()(unsigned int i) const
                {
                    ASSERTL1(( i >= 0) && (i < dimension()), "Invalid access to m_data via parenthesis operator");
                    return m_impl.data[i];
                }

                typename boost::call_traits<DataType>::reference operator[](unsigned int i)
                {
                    return m_impl.data[i];
                }

                typename boost::call_traits<DataType>::const_reference operator[](unsigned int i) const
                {
                    return m_impl.data[i];
                }

                typename boost::call_traits<DataType>::const_reference x() const
                {
                    ASSERTL1( dimension() >= 1, "Invalid use of Vector::x");
                    return (*this)(0);
                }

                typename boost::call_traits<DataType>::const_reference y() const
                {
                    ASSERTL1( dimension() >= 2, "Invalid use of Vector::y");
                    return (*this)(1);
                }

                typename boost::call_traits<DataType>::const_reference z() const
                {
                    ASSERTL1( dimension() >= 3, "Invalid use of Vector::z");
                    return (*this)(2);
                }

                typename boost::call_traits<DataType>::reference x()
                {
                    ASSERTL1( dimension() >= 1, "Invalid use of Vector::x");
                    return (*this)(0);
                }

                typename boost::call_traits<DataType>::reference y()
                {
                    ASSERTL1( dimension() >= 2, "Invalid use of Vector::y");
                    return (*this)(1);
                }

                typename boost::call_traits<DataType>::reference z()
                {
                    ASSERTL1( dimension() >= 3, "Invalid use of Vector::z");
                    return (*this)(2);
                }

                bool operator==(const NekVector<DataType, dim, space>& rhs)
                {
                    if( dimension() != rhs.dimension() )
                    {
                        return false;
                    }

                    for(unsigned int i = 0; i < dimension(); ++i)
                    {
                        // If you get a compile error here then you have to
                        // add a != operator to the DataType class.
                        if( m_impl.data[i] != rhs[i] )
                        {
                            return false;
                        }
                    }
                    return true;
                }

                bool operator!=(const NekVector<DataType, dim, space>& rhs)
                {
                    return !(*this == rhs);
                }

                /// Arithmetic Routines

                // Unitary operators
                NekVector<DataType, dim, space> operator-() const
                {
                    NekVector<DataType, dim, space> temp(*this);
                    for(int i=0; i < dimension(); ++i)
                    {
                        temp(i) = -temp(i);
                    }
                    return temp;
                }


                NekVector<DataType, dim, space>& operator+=(const NekVector<DataType, dim, space>& rhs)
                {
                    ASSERTL1( dimension() == rhs.dimensin(), "NekVector::operator+=, dimension of the two operands must be identical.");
                    for(unsigned int i=0; i < dimension(); ++i)
                    {
                        (*this)[i] += rhs[i];
                    }
                    return *this;
                }

                NekVector<DataType, dim, space>& operator-=(const NekVector<DataType, dim, space>& rhs)
                {
                    ASSERTL1( dimension() == rhs.dimensin(), "NekVector::operator-=, dimension of the two operands must be identical.");
                    for(unsigned int i=0; i < dimension(); ++i)
                    {
                        (*this)[i] -= rhs[i];
                    }
                    return *this;
                }

                NekVector<DataType, dim, space>& operator*=(typename boost::call_traits<DataType>::param_type rhs)
                {
                    for(unsigned int i=0; i < dimension(); ++i)
                    {
                        (*this)[i] *= rhs;
                    }
                    return *this;
                }

                DataType magnitude() const
                {
                    DataType result = DataType(0);

                    for(unsigned int i = 0; i < dimension(); ++i)
                    {
                        result += (*this)[i]*(*this)[i];
                    }
                    return sqrt(result);
                }

                DataType dot(const NekVector<DataType, dim, space>& rhs) const
                {
                    ASSERTL1( dimension() == rhs.dimensin(), "NekVector::dot, dimension of the two operands must be identical.");

                    DataType result = DataType(0);
                    for(unsigned int i = 0; i < dimension(); ++i)
                    {
                        result += (*this)[i]*rhs[i];
                    }

                    return result;
                }

                void normalize()
                {
                    DataType m = magnitude();
                    if( m > DataType(0) )
                    {
                        for(unsigned int i = 0; i < dimension(); ++i)
                        {
                            (*this)[i] /= m;
                        }
                    }
                }

            private:
                template<typename ImplDataType, unsigned int ImplSize, unsigned int ImplSpace>
                class VectorImpl
                {
                    public:
                        VectorImpl() :
                            data()
                        {
                            for(unsigned int i = 0; i < ImplSize; ++i)
                            {
                                data[i] = ImplDataType(0);
                            }
                        }

                        explicit VectorImpl(typename boost::call_traits<ImplDataType>::param_type a) :
                            data()
                        {
                            for(unsigned int i = 0; i < ImplSize; ++i)
                            {
                                data[i] = a;
                            }
                        }

                        VectorImpl(unsigned int size, typename boost::call_traits<ImplDataType>::param_type a) :
                            data()
                        {
                            unsigned int end = std::min(size, ImplSize);
                            for(unsigned int i = 0; i < end; ++i)
                            {
                                data[i] = a;
                            }

                            for(unsigned int i = end; i < ImplSize; ++i)
                            {
                                data[i] = ImplDataType(0);
                            }
                        }

                        VectorImpl(const DataType* const ptr) :
                            data()
                        {
                            for(unsigned int i = 0; i < ImplSize; ++i)
                            {
                                data[i] = ptr[i];
                            }
                        }

                        VectorImpl(const VectorImpl<ImplDataType, ImplSize, ImplSpace>& rhs) :
                            data()
                        {
                            for(unsigned int i = 0; i < ImplSize; ++i)
                            {
                                data[i] = rhs.data[i];
                            }
                        }

                        VectorImpl<ImplDataType, ImplSize, ImplSpace>& operator=(const VectorImpl<ImplDataType, ImplSize, ImplSpace>& rhs)
                        {
                            for(unsigned int i = 0; i < ImplSize; ++i)
                            {
                                data[i] = rhs.data[i];
                            }
                            return *this;
                        }

                        ~VectorImpl() {}

                        unsigned int dimension() const { return ImplSize; }
                        ImplDataType data[ImplSize];
                };

                // \brief Specialization for variable sized vectors.
                template<typename ImplDataType, unsigned int ImplSpace>
                class VectorImpl<ImplDataType, 0, ImplSpace>
                {
                    public:
                        VectorImpl() :
                            data(new ImplDataType[1]),
                            size(1)
                        {
                            data[0] = ImplDataType(0);
                        }

                        VectorImpl(unsigned int s, typename boost::call_traits<ImplDataType>::param_type a) :
                            data(new DataType[s]),
                            size(s)
                        {
                            for(unsigned int i = 0; i < size; ++i)
                            {
                                data[i] = a;
                            }
                        }

                        explicit VectorImpl(typename boost::call_traits<ImplDataType>::param_type a) :
                            data(new DataType[1]),
                            size(1)
                        {
                            data[0] = a;
                        }

                        VectorImpl(unsigned int s, const DataType* const ptr) :
                            data(new DataType[s]),
                            size(s)
                        {
                            for(unsigned int i = 0; i < size; ++i)
                            {
                                data[i] = ptr[i];
                            }
                        }


                        VectorImpl(const VectorImpl<ImplDataType, 0, ImplSpace>& rhs) :
                                data(new ImplDataType[rhs.size]),
                            size(rhs.size)
                        {
                            for(unsigned int i = 0; i < size; ++i)
                            {
                                data[i] = rhs.data[i];
                            }
                        }

                        ~VectorImpl() {}

                        VectorImpl<ImplDataType, 0, ImplSpace>& operator=(const VectorImpl<ImplDataType, 0, ImplSpace>& rhs)
                        {
                            boost::shared_array<ImplDataType> temp = boost::shared_array<DataType>(new ImplDataType[rhs.size]);
                            for(unsigned int i = 0; i < rhs.size; ++i)
                            {
                                temp[i] = rhs.data[i];
                            }

                            data = temp;
                            size = rhs.size;
                            return *this;
                        }

                        unsigned int dimension() const { return size; }

                        boost::shared_array<ImplDataType> data;
                        unsigned int size;
                };

                VectorImpl<DataType, dim, space> m_impl;
        };

        template<typename DataType, unsigned int dim, unsigned int space>
        NekVector<DataType, dim, space>
        operator+(const NekVector<DataType, dim, space>& lhs, const NekVector<DataType, dim>& rhs)
        {
            NekVector<DataType, dim, space> result(lhs);
            result += rhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekVector<DataType, dim, space>
        operator-(const NekVector<DataType, dim, space>& lhs, const NekVector<DataType, dim>& rhs)
        {
            NekVector<DataType, dim, space> result(lhs);
            result -= rhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekVector<DataType, dim, space>
        operator*(const NekVector<DataType, dim, space>& lhs, typename boost::call_traits<DataType>::param_type rhs)
        {
            NekVector<DataType, dim, space> result(lhs);
            result *= rhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekVector<DataType, dim, space>
        operator*(typename boost::call_traits<DataType>::param_type lhs, const NekVector<DataType, dim, space>& rhs)
        {
            NekVector<DataType, dim, space> result(rhs);
            result *= lhs;
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        NekVector<DataType, dim, space>
        normalize(const NekVector<DataType, dim, space>& rhs)
        {
            NekVector<DataType, dim, space> result(rhs);
            result.normalize();
            return result;
        }

        template<typename DataType, unsigned int dim, unsigned int space>
        DataType dot(const NekVector<DataType, dim, space>& lhs, const NekVector<DataType, dim>& rhs)
        {
            return lhs.dot(rhs);
        }
    }
}

#endif // NEKTAR_LIB_UTILITIES_NEK_VECTOR_HPP

/**
    $Log: NekVector.hpp,v $
    Revision 1.1  2006/06/01 09:12:42  kirby
    *** empty log message ***

    Revision 1.2  2006/05/25 03:02:40  bnelson
    Added Matrix/Vector multiplication.

    Revision 1.1  2006/05/04 18:57:44  kirby
    *** empty log message ***

    Revision 1.1  2006/04/11 02:00:43  bnelson
    Initial Revision


**/

