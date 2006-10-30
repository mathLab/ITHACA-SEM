///////////////////////////////////////////////////////////////////////////////
//
// File: NekFullMatrix.hpp
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
// Description: Full matrix
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_FULL_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_FULL_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <boost/call_traits.hpp>

#include <algorithm>

namespace Nektar
{
    template<typename DataType, unsigned int space>
    class NekMatrix<DataType, eFull, eNormal, space>
    {
        public:
            typedef NekMatrix<DataType, eFull, eNormal, space> ThisType;


            /// Necessary to allow construction of default constructed matrices.
            friend class MemoryManager;
            static const MemoryManager::MemoryPoolEnabler MemoryPoolEnabled = MemoryManager::eEnabled;

        public:
            /// \brief Create a matrix with the given size, initialized to the default value of DataType.
            ///
            /// It is not possible to create a block matrix with this constructor since the block size is not specified.
            NekMatrix(unsigned int rows, unsigned int columns) :
                m_rows(),
                m_columns(),
                m_data()
            {
                Initialize(rows, columns, DataType(0));
            }

            NekMatrix(unsigned int rows, unsigned int columns, DataType* ptr, MatrixDataHolderType t = eCopy) :
                m_rows(rows),
                m_columns(columns),
                m_data()
            {
                if( t == eCopy )
                {
                    m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
                    std::copy(ptr, ptr + m_rows*m_columns, begin());
                }
                else
                {
                    m_data = boost::shared_array<DataType>(ptr, DeleteNothing<DataType>());
                }
            }
             
            NekMatrix(unsigned int rows, unsigned int columns, const DataType* ptr) :
                m_rows(rows),
                m_columns(columns),
                m_data()
            {
                m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
                std::copy(ptr, ptr + m_rows*m_columns, begin());
            }
                
            /// It is not possible to create a block matrix with this constructor since the block size is not specified.
            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d) :
                m_rows(),
                m_columns(),
                m_data()
            {
                Initialize(rows, columns, d);
            }


            NekMatrix(const NekMatrix<DataType, eFull, eNormal, space>& rhs) :
                m_rows(rhs.m_rows),
                m_columns(rhs.m_columns),
                m_data(MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns))
            {
                std::copy(rhs.begin(), rhs.end(), begin());
            }

            template<typename ExpressionPolicyType>
            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                m_rows(rhs.GetMetadata().Rows),
                m_columns(rhs.GetMetadata().Columns),
                m_data(MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns))
        {
            BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eFull, eNormal, space> > ));
            rhs.Apply(*this);
        }

        NekMatrix<DataType, eFull, eNormal, space>& operator=(const NekMatrix<DataType, eFull, eNormal, space>& rhs)
        {
            NekMatrix<DataType, eFull, eNormal, space> temp(rhs);
            Swap(temp);
            return *this;
        }


        NekMatrix<DataType, eFull, eNormal, space>& operator=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
        {
            m_rows = rhs.GetRows();
            m_columns = rhs.GetColumns();
            m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
            std::fill(begin(), end(), DataType(0));
            
            for(unsigned int i = 0; i < m_rows; ++i)
            {
                (*this)(i,i) = rhs(i,i);    
            }
            return *this;
        }
        
        template<typename ExpressionPolicyType>
        NekMatrix<DataType, eFull, eNormal, space>& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
        {
            BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eFull, eNormal, space> > ));
            m_rows = rhs.GetMetadata().Rows;
            m_columns = rhs.GetMetadata().Columns;
            m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
            rhs.Apply(*this);
            return *this;
        }

        ~NekMatrix() {}

        unsigned int GetRows() const { return m_rows; }
        unsigned int GetColumns() const { return m_columns; }

        typename boost::call_traits<DataType>::reference operator()(unsigned int rowNumber, unsigned int colNumber)
        {
            ASSERTL2(rowNumber < m_rows, "Invalid row number to NekMatrix::operator()");
            ASSERTL2(colNumber < m_columns, "Invalid column number to NekMatrix::operator()");
            return m_data[rowNumber*m_columns + colNumber];
        }

        typename boost::call_traits<DataType>::const_reference operator()(unsigned int rowNumber, unsigned int colNumber) const
        {
            ASSERTL2(rowNumber < m_rows, "Invalid row number to NekMatrix::operator()");
            ASSERTL2(colNumber < m_columns, "Invalid column number to NekMatrix::operator()");
            return m_data[rowNumber*m_columns + colNumber];
        }

        DataType* GetPtr(unsigned int rowNumber, unsigned int colNumber)
        {
            return m_data.get() + rowNumber*m_columns + colNumber;
        }

        boost::shared_array<DataType> GetPtr()
        {
            return m_data;
        }

        boost::shared_array<const DataType> GetPtr() const
        {
            return m_data.get();
        }

        typedef DataType* iterator;
        typedef const DataType* const_iterator;

        iterator begin() { return m_data.get(); }

        iterator end() { return m_data.get() + m_columns*m_rows; }

        const_iterator begin() const { return m_data.get(); }

        const_iterator end() const { return m_data.get() + m_columns*m_rows; }

        void Negate()
        {
            for(iterator iter = begin(); iter != end(); ++iter)
            {
                *iter = -*iter;
            }
        }
        
        void Transpose()
        {
            for(unsigned int row = 0; row < m_rows; ++row)
            {
                for(unsigned int column = row+1; column < m_columns; ++column)
                {
                    std::swap( (*this)(row, column), (*this)(column, row));
                }
            }
        }

        expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eFull, eNormal, space> > >, expt::NegateOp> > operator-() const
        {
            return expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eFull, eNormal, space> > >, expt::NegateOp> >(
                    expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eFull, eNormal, space> > >(*this));
        }


        NekMatrix<DataType, eFull, eNormal, space> operator+=(const NekMatrix<DataType, eFull, eNormal, space>& rhs)
        {
            ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+");
            DataType* lhs_data = begin();
            const DataType* rhs_data = rhs.begin();

            for( ; lhs_data < end(); ++lhs_data, ++rhs_data )
            {
                *lhs_data += *rhs_data;
            }

            return *this;
        }

        // This is wrong as well.  What if this is diagonal?
        // enable if on the output.
        NekMatrix<DataType, eFull, eNormal, space> operator+=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
        {
            ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+");

            for(unsigned int i = 0; i < rhs.GetRows(); ++i)
            {
                (*this)(i,i) += rhs(i, i);
            }

            return *this;
        }


        NekMatrix<DataType, eFull, eNormal, space>& operator*=(const NekMatrix<DataType, eFull, eNormal, space>& rhs)
        {
            ASSERTL0(GetColumns() == rhs.GetRows(), "Invalid matrix dimensions in operator*");

            NekMatrix<DataType, eFull, eNormal, space> result(GetRows(), rhs.GetColumns());

            for(unsigned int i = 0; i < result.GetRows(); ++i)
            {
                for(unsigned int j = 0; j < result.GetColumns(); ++j)
                {
                    DataType t = DataType(0);

                    // Set the result(i,j) element.
                    for(unsigned int k = 0; k < GetColumns(); ++k)
                    {
                        t += (*this)(i,k)*rhs(k,j);
                    }
                    result(i,j) = t;
                }
            }

            Swap(result);

            return *this;
        }

        void Initialize(unsigned int rows, unsigned int columns)
        {
            m_rows = rows;
            m_columns = columns;
            m_data = MemoryManager::AllocateSharedArray<DataType>(m_rows*m_columns);
        }
        
        void Initialize(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d)
        {
            Initialize(rows, columns);
            std::fill(begin(), end(), d);
        }

    private:
        /// \brief Constructor used for block matrices.  Care must be used with this constructor
        ///         to initialize the matrix before use.
        NekMatrix() {}


        void Swap(NekMatrix<DataType, eFull, eNormal, space>& rhs)
        {
            std::swap(m_rows, rhs.m_rows);
            std::swap(m_columns, rhs.m_columns);
            std::swap(m_data, rhs.m_data);
        }

        unsigned int m_rows;
        unsigned int m_columns;
        boost::shared_array<DataType> m_data;
    };
    


}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_FULL_MATRIX_HPP

/**
    $Log:$
 **/

