///////////////////////////////////////////////////////////////////////////////
//
// File: NekDiagonalMatrix.hpp
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
// Description: Generic Matrix
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

#include <boost/shared_ptr.hpp>
#include <boost/call_traits.hpp>
#include <boost/concept_check.hpp>

#include <LibUtilities/NekMatrix.hpp>
#include <LibUtilities/NekMemoryManager.hpp>
#include <LibUtilities/ErrorUtil.hpp>
#include <LibUtilities/NekVector.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        template<typename DataType, unsigned int space>
        class NekMatrix<DataType, eDiagonal, space>
        {

            public:
                NekMatrix(unsigned int rows, unsigned int columns) :
                    m_rows(rows),
                    m_columns(columns),
                    m_data(new DataType[rows])
                {
                    ASSERTLO(rows == columns, "Diagonal matrices must be square.");
                }

                NekMatrix(unsigned int rows, unsigned int columns, const DataType* const ptr) :
                    m_rows(rows),
                    m_columns(columns),
                    m_data(new DataType[rows])
                {
                    ASSERTLO(rows == columns, "Diagonal matrices must be square.");
                    std::copy(ptr, ptr+rows, begin());
                }

                NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::param_type d) :
                    m_rows(rows),
                    m_columns(columns),
                    m_data(new DataType[rows])
                {
                    std::fill(begin(), end(), d);
                }

                unsigned int rows() const
                {
                    return m_rows;
                }

                unsigned int columns() const
                {
                    return m_columns;
                }

                inline typename boost::call_traits<DataType>::reference operator()
                    (unsigned int rowNumber, unsigned int colNumber)
                {
                    if( rowNumber == colNumber )
                    {
                        return m_data[rowNumber];
                    }
                    else
                    {
                        // TOOD - Find a good way to deal with this.  Probably
                        // need a proxy object with an assertion error.
                        static DataType result(0);
                        return result;
                    }
                }

                inline typename boost::call_traits<DataType>::const_reference operator()
                        (unsigned int rowNumber, unsigned int colNumber) const
                {
                    if( rowNumber == colNumber )
                    {
                        return m_data[rowNumber*m_columns + colNumber];
                    }
                    else
                    {
                        static DataType result(0);
                        return result;
                    }
                }


                typename boost::call_traits<DataType>::reference getValue
                        (unsigned int rowNumber, unsigned int colNumber)
                {
                    return (*this)(rowNumber, colNumber);
                }

                typename boost::call_traits<DataType>::const_reference getValue
                        (unsigned int rowNumber, unsigned int colNumber) const
                {
                    return (*this)(rowNumber, colNumber);
                }

                void setValue(unsigned int rowNumber, unsigned int colNumber, typename boost::call_traits<DataType>::param_type rhs)
                {
                    if( rowNumber == colNumber )
                    {
                        (*this)(rowNumber, colNumber) = rhs;
                    }
                    else
                    {
                        ASSERTL0(false, "Can't assign into non-diagonal element of a diagonal matrix.");
                    }
                }

                DataType* getPtr(unsigned int rowNumber, unsigned int colNumber)
                {
                    if( rowNumber == colNumber )
                    {
                        return &m_data[m_columns];
                    }
                    else
                    {
                        ASSERTL0(false, "Can't acces a non-diagonal element pointer of a diagonal matrix.");
                    }
                }

                typedef DataType* iterator;
                typedef const DataType* const_iterator;

                iterator begin()
                {
                    return &m_data[0];
                }

                iterator end()
                {
                    return &m_data[m_rows];
                }

                const_iterator begin() const
                {
                    return &m_data[0];
                }

                const_iterator end() const
                {
                    return &m_data[m_rows*m_columns];
                }

                NekMatrixForm GetForm() const
                {
                    return eDiagonal;
                }

                NekMatrix<DataType, eFull, space> operator+=(const NekMatrix<DataType, eFull, space>& rhs)
                {
                    ASSERTL0(rows() == rhs.rows() && columns() == rhs.columns(), "Matrix dimensions must agree in operator+");
                    DataType* lhs_data = begin();
                    const DataType* rhs_data = rhs.begin();

                    for(unsigned int i = 0; i < m_rows*m_columns; ++i)
                    {
                        lhs_data[i] += rhs_data[i];
                    }

                    return *this;
                }

                NekMatrix<DataType, eFull, space> operator+=(const NekMatrix<DataType, eDiagonal, space>& rhs)
                {
                    ASSERTL0(rows() == rhs.rows() && columns() == rhs.columns(), "Matrix dimensions must agree in operator+");
                    DataType* lhs_data = begin();
                    const DataType* rhs_data = rhs.begin();

                    return *this;
                }

            private:
                unsigned int m_rows;
                unsigned int m_columns;
                boost::shared_array<DataType> m_data;
        };
    }
}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

/**
    $Log: NekDiagonalMatrix.hpp,v $

 **/


