///////////////////////////////////////////////////////////////////////////////
//
// File: NekConstantSizedVector.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_CONSTANT_SIZED_VECTOR_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_CONSTANT_SIZED_VECTOR_HPP

#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorTypeTraits.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorCommon.hpp>

#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>


#include <algorithm>
#include <boost/static_assert.hpp>


namespace Nektar
{

    // \param DataType The type of data held by each element of the vector.
    // \param dim The number of elements in the vector.  If set to 0, the vector
    //            will have a variable number of elements.
    // \param space The space of the vector.
    template<typename DataType, unsigned int dim, unsigned int space>
    class NekVector
    {
        public:
            /// \brief Creates a constant sized vector with each element initialized to
            ///        the default value for DataType.
            NekVector();

            /// \brief Creates a constant sized vector with each element set to a.
            /// \param a The value to assign to each element of the vector.
            NekVector(typename boost::call_traits<DataType>::const_reference a);

            /// \brief Creates a vector from the elements in a delimited string.
            /// \param vectorValues A string the the vector values.
            /// 
            /// 
            explicit NekVector(const std::string& vectorValues);

            NekVector(typename boost::call_traits<DataType>::const_reference x,
                      typename boost::call_traits<DataType>::const_reference y,
                      typename boost::call_traits<DataType>::const_reference z);

            NekVector(typename boost::call_traits<DataType>::const_reference x,
                      typename boost::call_traits<DataType>::const_reference y,
                      typename boost::call_traits<DataType>::const_reference z,
                      typename boost::call_traits<DataType>::const_reference w);

            NekVector(const NekVector<DataType, dim, space>& rhs);

            explicit NekVector(const DataType* const ptr);

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            template<typename ExpressionPolicyType>
            NekVector(const Expression<ExpressionPolicyType>& rhs) :
                m_impl(rhs.GetMetadata().Rows, DataType())
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekVector<DataType, dim, space> > ));
                rhs.Evaluate(*this);
            }
#endif
            ~NekVector()
            {
#ifdef _DEBUG
                std::fill_n(m_impl, dim, DataType());
#endif //_DEBUG
            }

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            template<typename ExpressionPolicyType>
            NekVector<DataType, dim, space>& operator=(const Expression<ExpressionPolicyType>& rhs)
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekVector<DataType, dim, space> > ));
                rhs.Evaluate(*this);
                return *this;
            }
#endif

            NekVector<DataType, dim, space>& operator=(const NekVector<DataType, dim, space>& rhs)
            {
                std::copy(rhs.m_impl, rhs.m_impl+dim, m_impl);
                return *this;
            }

            /// \brief Returns the number of dimensions for the point.
            unsigned int GetDimension() const
            {
                return dim;
            }

            /// \brief Treating the vector as a column vector, how many rows it has.
            unsigned int GetRows() const
            {
                return dim;
            }

            DataType* GetRawPtr()
            {
                return &m_impl[0];
            }
            
            const DataType* GetRawPtr() const
            {
                return &m_impl[0];
            }
            
            typedef DataType* iterator;
            typedef const DataType* const_iterator;
            
            iterator begin() { return GetRawPtr(); }
            iterator end() { return GetRawPtr() + GetDimension(); }
            
            const_iterator begin() const { return GetRawPtr(); }
            const_iterator end() const { return GetRawPtr() + GetDimension(); }
            
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
                ASSERTL1((i >= 0) && (i < GetDimension()), "Invalid access to m_data via parenthesis operator");
                return m_impl[i];
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int i) const
            {
                ASSERTL1(( i >= 0) && (i < GetDimension()), "Invalid access to m_data via parenthesis operator");
                return m_impl[i];
            }

            typename boost::call_traits<DataType>::reference operator[](unsigned int i)
            {
                return m_impl[i];
            }

            typename boost::call_traits<DataType>::const_reference operator[](unsigned int i) const
            {
                return m_impl[i];
            }

            typename boost::call_traits<DataType>::const_reference x() const
            {
                BOOST_STATIC_ASSERT(dim >= 1);
                return m_impl[0];
            }

            typename boost::call_traits<DataType>::const_reference y() const
            {
                BOOST_STATIC_ASSERT(dim >= 2);
                return m_impl[1];
            }

            typename boost::call_traits<DataType>::const_reference z() const
            {
                BOOST_STATIC_ASSERT(dim >= 3);
                return m_impl[2];
            }

            typename boost::call_traits<DataType>::const_reference w() const
            {
                BOOST_STATIC_ASSERT(dim >= 4);
                return m_impl[3];
            }

            typename boost::call_traits<DataType>::reference x() 
            {
                BOOST_STATIC_ASSERT(dim >= 1);
                return m_impl[0];
            }

            typename boost::call_traits<DataType>::reference y() 
            {
                BOOST_STATIC_ASSERT(dim >= 2);
                return m_impl[1];
            }

            typename boost::call_traits<DataType>::reference z() 
            {
                BOOST_STATIC_ASSERT(dim >= 3);
                return m_impl[2];
            }

            typename boost::call_traits<DataType>::reference w() 
            {
                BOOST_STATIC_ASSERT(dim >= 4);
                return m_impl[3];
            }

            void SetX(typename boost::call_traits<DataType>::const_reference val)
            {
                BOOST_STATIC_ASSERT(dim >= 1);
                m_impl[0] = val;
            }

            void SetY(typename boost::call_traits<DataType>::const_reference val)
            {
                BOOST_STATIC_ASSERT(dim >= 2);
                m_impl[1] = val;
            }

            void SetZ(typename boost::call_traits<DataType>::const_reference val)
            {
                BOOST_STATIC_ASSERT(dim >= 3);
                m_impl[2] = val;
            }

            void SetW(typename boost::call_traits<DataType>::const_reference val)
            {
                BOOST_STATIC_ASSERT(dim >= 4);
                m_impl[3] = val;
            }

            /// Arithmetic Routines

            // Unitary operators
            NekVector<DataType, dim, space> operator-() const { return Negate(*this); }

            NekVector<DataType, dim, space>& operator+=(const NekVector<DataType, dim, space>& rhs)
            {
                PlusEqual(*this, rhs);
                return *this;
            }

            NekVector<DataType, dim, space>& operator-=(const NekVector<DataType, dim, space>& rhs)
            {
                MinusEqual(*this, rhs);
                return *this;
            }

            NekVector<DataType, dim, space>& operator*=(typename boost::call_traits<DataType>::const_reference rhs)
            {
                TimesEqual(*this, rhs);
                return *this;
            }
            
            NekVector<DataType, dim, space>& operator/=(typename boost::call_traits<DataType>::const_reference rhs)
            {
                DivideEqual(*this, rhs);
                return *this;
            }

            DataType Magnitude() const { return Nektar::Magnitude(*this); }
            
            DataType Dot(const NekVector<DataType, dim, space>& rhs) const { return Nektar::Dot(*this, rhs); }
            
            void Normalize() { return Nektar::Normalize(*this); }
            

            NekVector<DataType, dim, space> Cross(const NekVector<DataType, dim, space>& rhs) const
            {
                return Nektar::Cross(*this, rhs);
            }

            std::string AsString() const { return Nektar::AsString(*this); }

            // Norms
            DataType L1Norm() const { return Nektar::L1Norm(*this); }
            DataType L2Norm() const { return Nektar::L2Norm(*this); }
            DataType InfinityNorm() const { return Nektar::InfinityNorm(*this); }
            
        private:
            DataType m_impl[dim];

    };    

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>::NekVector() :
        m_impl()
    {
        std::fill_n(m_impl, dim, DataType());
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>::NekVector(typename boost::call_traits<DataType>::const_reference a) :
        m_impl()
    {
        std::fill_n(m_impl, dim, a);
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>::NekVector(const std::string& vectorValues) :
        m_impl()
    {
        std::vector<DataType> values = FromString<DataType>(vectorValues);

        ASSERTL0(values.size() == dim, "Error converting string values to vector");

        std::copy(values.begin(), values.end(), &m_impl[0]);
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>::NekVector(typename boost::call_traits<DataType>::const_reference x,
              typename boost::call_traits<DataType>::const_reference y,
              typename boost::call_traits<DataType>::const_reference z) :
        m_impl()
    {
        BOOST_STATIC_ASSERT(dim == 3);
        m_impl[0] = x;
        m_impl[1] = y;
        m_impl[2] = z;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>::NekVector(typename boost::call_traits<DataType>::const_reference x,
              typename boost::call_traits<DataType>::const_reference y,
              typename boost::call_traits<DataType>::const_reference z,
              typename boost::call_traits<DataType>::const_reference w) :
        m_impl()
    {
        BOOST_STATIC_ASSERT(dim == 4);
        m_impl[0] = x;
        m_impl[1] = y;
        m_impl[2] = z;
        m_impl[3] = w;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>::NekVector(const NekVector<DataType, dim, space>& rhs) :
        m_impl()
    {
        std::copy(rhs.m_impl, rhs.m_impl+dim, m_impl);
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space>::NekVector(const DataType* const ptr) :
        m_impl()
    {
        std::copy(ptr, ptr+dim, &m_impl[0]);
    }

}

#endif // NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_CONSTANT_SIZED_VECTOR_HPP
