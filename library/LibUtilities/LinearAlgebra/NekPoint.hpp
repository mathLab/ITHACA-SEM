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

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/call_traits.hpp>

#include <cmath>
#include <type_traits>
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
                ASSERTL0(i < dim::Value,
                        "Invalid access to NekPoint data via parenthesis "
                        "operator: index out of range");
                return m_data[i];
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int i) const
            {
                ASSERTL0(i < dim::Value,
                        "Invalid access to NekPoint data via parenthesis "
                        "operator: index out of range");
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
                static_assert(dim::Value >= 1, "invalid dimension");
                return m_data[0];
            }

            typename boost::call_traits<DataType>::const_reference y() const
            {
                static_assert(dim::Value >= 2, "invalid dimension");
                return (*this)[1];
            }

            typename boost::call_traits<DataType>::const_reference z() const
            {
                static_assert(dim::Value >= 3, "invalid dimension");
                return (*this)[2];
            }

            typename boost::call_traits<DataType>::const_reference a() const
            {
                static_assert(dim::Value >= 1, "invalid dimension");
                return m_data[0];
            }

            typename boost::call_traits<DataType>::const_reference b() const
            {
                static_assert(dim::Value >= 2, "invalid dimension");
                return (*this)[1];
            }

            typename boost::call_traits<DataType>::const_reference c() const
            {
                static_assert(dim::Value >= 3, "invalid dimension");
                return (*this)[2];
            }

            typename boost::call_traits<DataType>::const_reference r() const
            {
                static_assert(dim::Value >= 1, "invalid dimension");
                return m_data[0];
            }

            typename boost::call_traits<DataType>::const_reference s() const
            {
                static_assert(dim::Value >= 2, "invalid dimension");
                return (*this)[1];
            }

            typename boost::call_traits<DataType>::const_reference t() const
            {
                static_assert(dim::Value >= 3, "invalid dimension");
                return (*this)[2];
            }

            void SetX(typename boost::call_traits<DataType>::const_reference val)
            {
                static_assert(dim::Value >= 1, "invalid dimension");
                m_data[0] = val;
            }

            void SetY(typename boost::call_traits<DataType>::const_reference val)
            {
                static_assert(dim::Value >= 2, "invalid dimension");
                m_data[1] = val;
            }

            void SetZ(typename boost::call_traits<DataType>::const_reference val)
            {
                static_assert(dim::Value >= 2, "invalid dimension");
                m_data[2] = val;
            }

            typename boost::call_traits<DataType>::reference x()
            {
                static_assert(dim::Value >= 1, "invalid dimension");
                return (*this)(0);
            }

            typename boost::call_traits<DataType>::reference y()
            {
                static_assert(dim::Value >= 2, "invalid dimension");
                return (*this)(1);
            }

            typename boost::call_traits<DataType>::reference z()
            {
                static_assert(dim::Value >= 3, "invalid dimension");
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

}

#endif // NEKTAR_LIB_UTILITIES_NEK_POINT_HPP
