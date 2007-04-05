///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrixObject.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_OBJECT_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_OBJECT_HPP

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekDiagonalMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekFullMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekBlockDiagonalMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekBlockFullMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>

namespace Nektar
{
    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space, typename enabled>
    class NekMatrix
    {
        public:
            typedef NekMatrix<DataType, form, BlockType, space, void> ThisType;
            typedef NekMatrixStoragePolicy<DataType, form> StoragePolicy;
            typedef NekMatrixArithmeticPolicy<DataType, form> ArtihmeticPolicy;
            typedef NekMatrixAssignmentPolicy<DataType, form> AssignmentPolicy;

        public:
            /// \brief Creates an empty, 0x0 matrix.  
            NekMatrix() :
                m_rows(0),
                m_columns(0),
                m_data()
            {
            }
    
            /// \brief Create a matrix with the given size, initialized to the default value of DataType.
            ///
            /// It is not possible to create a block matrix with this constructor since the block size is not specified.
            NekMatrix(unsigned int rows, unsigned int columns) :
                m_rows(rows),
                m_columns(columns),
                m_data()
            {
                // If you get a compiler error here, then you are trying to use this constructor with a block matrix.
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                DataType d(0);
                Initialize(rows, columns, d);
            }

            NekMatrix(unsigned int rows, unsigned int columns, DataType* ptr, PointerWrapper t = eCopy) :
                m_rows(rows),
                m_columns(columns),
                m_data()
            {
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                Initialize(rows, columns, ptr, t);
            }

            NekMatrix(unsigned int rows, unsigned int columns, const DataType* const ptr) :
                m_rows(rows),
                m_columns(columns),
                m_data()
            {
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                Initialize(rows, columns, ptr);
            }

            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d) :
                m_rows(),
                m_columns(),
                m_data()
            {
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                Initialize(rows, columns, d);
            }

            NekMatrix(const ThisType& rhs) :
                m_rows(rhs.m_rows),
                m_columns(rhs.m_columns),
                m_data()
            {
                BOOST_STATIC_ASSERT(BlockType == eNormal);
                Initialize(m_rows, m_columns, rhs.m_data.get());
            }

            ThisType& operator=(const ThisType& rhs)
            {
                ThisType temp(rhs);
                Swap(temp);
                return *this;
            }

            template<NekMatrixForm rhs_form>
            ThisType& operator=(const NekMatrix<DataType, rhs_form, BlockType, space>& rhs)
            {
                AssignmentPolicy::Assign(*this, rhs);
                return *this;
            }

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

            template<typename ExpressionPolicyType>
            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                m_rows(rhs.GetMetadata().Rows),
                m_columns(rhs.GetMetadata().Columns),
                m_data()
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, ThisType> ));
                Initialize(m_rows, m_columns);
                rhs.Apply(*this);
            }

            template<typename ExpressionPolicyType>
            NekMatrix& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, ThisType> ));
                m_rows = rhs.GetMetadata().Rows;
                m_columns = rhs.GetMetadata().Columns;
                Initialize(m_rows, m_columns);
                rhs.Apply(*this);
                return *this;
            }
#endif

            ~NekMatrix() 
            {
                m_rows = 0;
                m_columns = 0;
                m_data.reset();
            }


            unsigned int GetRows() const { return m_rows; }
            unsigned int GetColumns() const { return m_columns; }

            typename boost::call_traits<DataType>::reference operator()(unsigned int rowNumber, unsigned int colNumber)
            {
                ASSERTL2(rowNumber < m_rows, "Invalid row number to NekMatrix::operator()");
                ASSERTL2(colNumber < m_columns, "Invalid column number to NekMatrix::operator()");
                return StoragePolicy::GetData(rowNumber, colNumber, m_rows, m_columns, m_data);
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int rowNumber, unsigned int colNumber) const
            {
                ASSERTL2(rowNumber < m_rows, "Invalid row number to NekMatrix::operator()");
                ASSERTL2(colNumber < m_columns, "Invalid column number to NekMatrix::operator()");
                return StoragePolicy::GetConstData(rowNumber, colNumber, m_rows, m_columns, m_data);
            }

            //DataType* GetPtr(unsigned int rowNumber, unsigned int colNumber)
            //{
            //    return m_data.get() + rowNumber*m_columns + colNumber;
            //}

            SharedArray<DataType>& GetPtr()
            {
                return m_data;
            }

            const SharedArray<DataType>& GetPtr() const
            {
                return m_data;
            }

            typedef DataType* iterator;
            typedef const DataType* const_iterator;

            iterator begin() 
            { 
                return m_data.get();
            }

            iterator end() 
            { 
                return m_data.get() + StoragePolicy::NumStorageElements(m_rows, m_columns);
            }

            const_iterator begin() const 
            { 
                return m_data.get();
            }

            const_iterator end() const 
            { 
                return m_data.get() + StoragePolicy::NumStorageElements(m_rows, m_columns);
            }

            void Negate()
            {
                for(iterator iter = begin(); iter != end(); ++iter)
                {
                    *iter = -*iter;
                }
            }
            
            void Transpose()
            {
                StoragePolicy::Transpose(m_rows, m_columns, m_data);
            }

            void Invert()
            {
                StoragePolicy::Invert(m_rows, m_columns, m_data);
            }

            /// \brief Performs an LU factorization of the matrix.
            ///
            /// See getrf.
            //void Factorize()
            //{
            //    StoragePolicy::Factorize(m_data, m_rows, m_columns);
            //}

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<ThisType> >, expt::NegateOp> > operator-() const
            {
                return expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<ThisType> >, expt::NegateOp> >(
                        expt::Expression<expt::ConstantExpressionPolicy<ThisType> >(*this));
            }
#endif
            template<NekMatrixForm rhs_form>
            ThisType& operator+=(const NekMatrix<DataType, rhs_form, BlockType, space>& rhs)
            {
                ArtihmeticPolicy::PlusEqual(*this, rhs);
                return *this;
            }

            //// This is wrong as well.  What if this is diagonal?
            //// enable if on the output.
            //NekMatrix<DataType, eFull, eNormal, space> operator+=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            //{
            //    ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");

            //    for(unsigned int i = 0; i < rhs.GetRows(); ++i)
            //    {
            //        (*this)(i,i) += rhs(i, i);
            //    }

            //    return *this;
            //}

            ThisType& operator-=(const ThisType& rhs)
            {
                ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator-=");
                DataType* lhs_data = begin();
                const DataType* rhs_data = rhs.begin();

                for( ; lhs_data < end(); ++lhs_data, ++rhs_data )
                {
                    *lhs_data -= *rhs_data;
                }

                return *this;
            }

            //// This is wrong as well.  What if this is diagonal?
            //// enable if on the output.
            //NekMatrix<DataType, eFull, eNormal, space> operator-=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            //{
            //    ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator-=");

            //    for(unsigned int i = 0; i < rhs.GetRows(); ++i)
            //    {
            //        (*this)(i,i) -= rhs(i, i);
            //    }

            //    return *this;
            //}

            ThisType& operator*=(const ThisType& rhs)
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

            

        private:
            void Initialize(unsigned int rows, unsigned int columns)
            {
                m_rows = rows;
                m_columns = columns;
                unsigned int storageSize = StoragePolicy::NumStorageElements(m_rows, m_columns);
                m_data = MemoryManager::AllocateSharedArray<DataType>(storageSize);
            }
            
            void Initialize(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference d)
            {
                Initialize(rows, columns);
                std::fill(begin(), end(), d);
            }

            void Initialize(unsigned int rows, unsigned int columns, DataType* ptr, 
                                   PointerWrapper t)
            {
                m_rows = rows;
                m_columns = columns;
                if( t == eCopy )
                {
                    Initialize(rows, columns, ptr);
                }
                else
                {
                    m_data = SharedArray<DataType>(ptr, StoragePolicy::NumStorageElements(rows, columns),
                        DeleteNothing<DataType>());
                }
            }

            void Initialize(unsigned int rows, unsigned int columns, const DataType* const ptr)
            {
                Initialize(rows, columns);
                std::copy(ptr, ptr + StoragePolicy::NumStorageElements(rows, columns), begin());
            }

            void Swap(NekMatrix<DataType, eFull, eNormal, space>& rhs)
            {
                std::swap(m_rows, rhs.m_rows);
                std::swap(m_columns, rhs.m_columns);
                std::swap(m_data, rhs.m_data);
            }

            unsigned int m_rows;
            unsigned int m_columns;
            SharedArray<DataType> m_data;
    };
    
    

    template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
    std::ostream& operator<<(std::ostream& os, const NekMatrix<DataType, form, BlockType, space>& rhs)
    {
        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            os << "[";
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                os << rhs(i,j);
                if( j != rhs.GetColumns() - 1 )
                {
                    os << ", ";
                }
            }
            os << "]";
            if( i != rhs.GetRows()-1 )
            {
                os << std::endl;
            }
        }
        return os;
    }

}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_OBJECT_HPP
