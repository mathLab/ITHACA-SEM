///////////////////////////////////////////////////////////////////////////////
//
// File: NekBlockDiagonalMatrix.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_DIAGONAL_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_DIAGONAL_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekBlockMatrix.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/Memory/DeleteNothing.hpp>

#include <boost/call_traits.hpp>

#include <algorithm>
#include <vector>

namespace Nektar
{
    template<typename DataType, MatrixBlockType BlockType, unsigned int space>
    class NekMatrix<DataType, eDiagonal, BlockType, space, typename boost::disable_if_c<BlockType==eNormal, void>::type >
    {
        public:
            typedef NekMatrix<DataType, eDiagonal, BlockType, space, void> ThisType;
            typedef NekMatrix<DataType, eFull, eNormal, space, void> InnerMatrixType;
            typedef typename BlockMatrixDataType<BlockType, NekMatrix<DataType, eFull, eNormal, space, void> >::ResultType InnerDataType;
            
            
        public:
            NekMatrix(unsigned int numberOfElements, unsigned int blockRows, unsigned int blockColumns) :
                m_numberOfElements(numberOfElements),
                m_rows(numberOfElements*blockRows),
                m_columns(numberOfElements*blockColumns),
                m_blockColumns(m_rows, blockColumns),
                m_blockRows(m_columns, blockRows),
                m_data(MemoryManager::AllocateSharedArray<InnerDataType>(m_numberOfElements))
            {
                ASSERTL0(blockRows == blockColumns, "ERROR: Block diagonal matrices must consist of square blocks.");
                for(unsigned int i = 0; i < m_numberOfElements; ++i)
                {
                    GetObj(m_data[i])->Initialize(blockRows, blockColumns);
                }
            }
            
            NekMatrix(const NekMatrix<DataType, eDiagonal, BlockType, space>& rhs) :
                m_numberOfElements(rhs.m_numberOfElements),
                m_rows(rhs.m_rows),
                m_columns(rhs.m_columns),
                m_blockColumns(rhs.m_blockColumns),
                m_blockRows(rhs.m_blockRows),
                m_data(MemoryManager::AllocateSharedArray<InnerDataType>(m_numberOfElements))
            {
                for(unsigned int i = 0; i < m_numberOfElements; ++i)
                {
                    m_data[i] = rhs.m_data[i];
                }
            }
                
            NekMatrix<DataType, eDiagonal, BlockType, space>& operator=(const NekMatrix<DataType, eDiagonal, BlockType, space>& rhs)
            {
                NekMatrix<DataType, eDiagonal, BlockType, space> temp(rhs);
                Swap(temp);
                return *this;
            }
            
            template<typename ExpressionPolicyType>
            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                m_numberOfElements(rhs.GetMetadata().Rows),
                m_rows(rhs.GetMetadata().Rows),
                m_columns(rhs.GetMetadata().Columns),
                m_blockColumns(rhs.GetMetadata().BlockColumns),
                m_blockRows(rhs.GetMetadata().BlockRows),
                m_data(MemoryManager::AllocateSharedArray<InnerDataType>(m_numberOfElements))
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eDiagonal, BlockType, space> > ));
                rhs.Apply(*this);
            }
            
            template<typename ExpressionPolicyType>
            NekMatrix<DataType, eDiagonal, eNormal, space>& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eDiagonal, BlockType, space> > ));
                
                m_numberOfElements = rhs.GetMetadata().Rows;
                m_rows = rhs.GetMetadata().Rows;
                m_columns = rhs.GetMetadata().Columns;
                m_blockColumns = rhs.GetMetadata().BlockColumns;
                m_blockRows = rhs.GetMetadata().BlockRows;
                m_data = MemoryManager::AllocateSharedArray<InnerDataType>(m_numberOfElements);

                rhs.Apply(*this);
                return *this;
            }
                
            
            typename boost::call_traits<DataType>::reference operator()(unsigned int rowNumber, unsigned int colNumber)
            {
                unsigned int blockRow = rowNumber/m_blockRows[0];
                unsigned int blockColumn = colNumber/m_blockColumns[0];
                unsigned int innerBlockRow = rowNumber%m_blockRows[0];
                unsigned int innerBlockColumn = colNumber%m_blockColumns[0];
                
                if( blockRow == blockColumn )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via operator()");
                    return GetBlock(blockRow, blockColumn)(innerBlockRow, innerBlockColumn);
                }
                else
                {
                    static InnerDataType zeroElement(m_blockRows[0], m_blockColumns[0], DataType(0));
                    return zeroElement(innerBlockRow, innerBlockColumn);
                }
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int rowNumber, unsigned int colNumber) const
            {
                unsigned int blockRow = rowNumber/m_blockRows[0];
                unsigned int blockColumn = colNumber/m_blockColumns[0];
                unsigned int innerBlockRow = rowNumber%m_blockRows[0];
                unsigned int innerBlockColumn = colNumber%m_blockColumns[0];
                
                if( blockRow == blockColumn )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via operator()");
                    return GetBlock(blockRow, blockColumn)(innerBlockRow, innerBlockColumn);
                }
                else
                {
                    static InnerDataType zeroElement(m_blockRows[0], m_blockColumns[0], DataType(0));
                    return zeroElement(innerBlockRow, innerBlockColumn);
                }
            }

            typename boost::call_traits<InnerDataType>::reference GetBlock(unsigned int rowNumber, unsigned int colNumber)
            {
                if( rowNumber == colNumber )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via GetBlock");
                    return m_data[rowNumber];
                }
                else
                {
                    static InnerDataType zeroElement(m_blockRows[0], m_blockColumns[0], DataType(0));
                    return zeroElement;
                }
            }

            typename boost::call_traits<InnerDataType>::const_reference GetBlock(unsigned int rowNumber, unsigned int colNumber) const
            {
                if( rowNumber == colNumber )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via GetBlock");
                    return m_data[rowNumber];
                }
                else
                {
                    static InnerDataType zeroElement(m_blockRows[0], m_blockColumns[0], DataType(0));
                    return zeroElement;
                }
            }
            
            unsigned int GetRows() const { return m_rows; }
            unsigned int GetColumns() const { return m_columns; }
            const std::vector<unsigned int> GetBlockRows() const { return m_blockRows; }
            const std::vector<unsigned int> GetBlockColumns() const { return m_blockColumns; }

            /// \brief Return a full matrix version.
            boost::shared_array<DataType> GetPtr()
            {
                boost::shared_array<DataType> result(new DataType[GetRows()*GetColumns()]);
                //std::fill(result, result+GetRows()*GetColumns(), DataType(0));
                
                for(unsigned int i = 0; i < GetRows(); ++i)
                {
                    for(unsigned int j = 0; j < GetColumns(); ++j)
                    {
                        result[i*GetColumns() + j] = (*this)(i,j);
                    }
                }
                
                return result;
            }
            
            ThisType& operator+=(const ThisType& rhs)
            {
                for(unsigned int i = 0; i < m_numberOfElements; ++i)
                {
                    m_data[i] += rhs.m_data[i];
                }
                
                return *this;
            }
            
            ThisType& operator*=(const ThisType& rhs)
            {
                ASSERTL1(GetColumns() == rhs.GetRows(), "Invalid matrix dimensions in operator*");
                ASSERTL1(m_numberOfElements == rhs.m_numberOfElements, "Invalid matrix dimensions in block diagonal operator*=");
                
                for(unsigned int i = 0; i < m_numberOfElements; ++i)
                {
                    m_data[i] *= rhs.m_data[i];
                }
                
                return *this;
            }
            
        private:
            
            void Swap(NekMatrix<DataType, eDiagonal, BlockType, space>& rhs)
            {
                std::swap(m_numberOfElements, rhs.m_numberOfElements);
                std::swap(m_rows, rhs.m_rows);
                std::swap(m_columns, rhs.m_columns);
                std::swap(m_blockColumns, rhs.m_blockColumns);
                std::swap(m_blockRows, rhs.m_blockRows);
                std::swap(m_data, rhs.m_data);
            }
                
            unsigned int m_numberOfElements;
            unsigned int m_rows;
            unsigned int m_columns;
            
            // For column i, the number of columns in that block.
            std::vector<unsigned int> m_blockColumns;
            
            // For row i, the number of rows in that block.
            std::vector<unsigned int> m_blockRows;
            
            boost::shared_array<InnerDataType> m_data;
            
    };
/*                public:
                    typedef NekMatrix<DataType, eDiagonal, eNormal, space> ThisType;


            /// Necessary to allow construction of default constructed matrices.
                    friend class MemoryManager;
                    static const MemoryManager::MemoryPoolEnabler MemoryPoolEnabled = MemoryManager::eEnabled;

                public:

            /// \brief Create an identity matrix.
                    ///
            /// It is not possible to create a block matrix with this constructor since the block size is not specified.
            explicit NekMatrix(unsigned int numberOfElements) :
                m_numberOfElements(numberOfElements),
            m_data(MemoryManager::AllocateSharedArray<DataType>(numberOfElements))
            {
                std::fill(begin(), end(), DataType(1));
            }

            NekMatrix(unsigned int numberOfElements, const DataType* const ptr) :
                    m_numberOfElements(numberOfElements),
            m_data(MemoryManager::AllocateSharedArray<DataType>(numberOfElements))
            {
                std::copy(ptr, ptr+m_numberOfElements, begin());
            }

            NekMatrix(unsigned int numberOfElements, DataType* ptr, MatrixDataHolderType t = eCopy) :
                    m_numberOfElements(numberOfElements),
            m_data()
            {
                if( t == eCopy )
                {
                    m_data = MemoryManager::AllocateSharedArray<DataType>(numberOfElements);
                    std::copy(ptr, ptr+m_numberOfElements, begin());
                }
                else
                {
                    m_data = boost::shared_array<DataType>(ptr, DeleteNothing<DataType>());
                }
            }

            NekMatrix(unsigned int numberOfElements, typename boost::call_traits<DataType>::const_reference d) :
                    m_numberOfElements(numberOfElements),
            m_data(MemoryManager::AllocateSharedArray<DataType>(numberOfElements))
            {
                std::fill(begin(), end(), d);
            }


            NekMatrix(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs) :
                    m_numberOfElements(rhs.m_numberOfElements),
            m_data(MemoryManager::AllocateSharedArray<DataType>(m_numberOfElements))
            {
                std::copy(rhs.begin(), rhs.end(), begin());
            }


            template<typename ExpressionPolicyType>
            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                    m_numberOfElements(rhs.GetMetadata().Rows),
            m_data(MemoryManager::AllocateSharedArray<DataType>(m_numberOfElements))
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eDiagonal, eNormal, space> > ));
                rhs.Apply(*this);
            }

            NekMatrix<DataType, eDiagonal, eNormal, space>& operator=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            {
                NekMatrix<DataType, eDiagonal, eNormal, space> temp(rhs);
                Swap(temp);
                return *this;
            }

            template<typename ExpressionPolicyType>
                    NekMatrix<DataType, eDiagonal, eNormal, space>& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eDiagonal, eNormal, space> > ));
                
                m_numberOfElements = rhs.GetMetadata().Rows;
                m_data = MemoryManager::AllocateSharedArray<DataType>(m_numberOfElements);
                rhs.Apply(*this);
                return *this;
            }

            ~NekMatrix() {}

            unsigned int GetRows() const { return m_numberOfElements; }
            unsigned int GetColumns() const { return m_numberOfElements; }

            typename boost::call_traits<DataType>::reference operator()(unsigned int rowNumber, unsigned int colNumber)
            {
                if( rowNumber == colNumber )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via operator()");
                    return m_data[rowNumber];
                }
                else
                {
                    static DataType zeroElement(0);
                    return zeroElement;
                }
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int rowNumber, unsigned int colNumber) const
            {
                if( rowNumber == colNumber )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via operator()");
                    return m_data[rowNumber];
                }
                else
                {
                    static DataType zeroElement(0);
                    return zeroElement;
                }
            }

            DataType* GetPtr(unsigned int rowNumber, unsigned int colNumber)
            {
                ASSERTL1(rowNumber == colNumber, "NekMatrix::GetPtr requires equal row an column numbers for diagonal matrices.");
                return m_data.get() + rowNumber;
            }

            DataType* GetPtr()
            {
                return m_data.get();
            }

            const DataType* GetPtr() const
            {
                return m_data.get();
            }

            typedef DataType* iterator;
            typedef const DataType* const_iterator;

            iterator begin()
            {
                return m_data.get();
            }

            iterator end()
            {
                return m_data.get() + m_numberOfElements;
            }

            const_iterator begin() const
            {
                return m_data.get();
            }

            const_iterator end() const
            {
                return m_data.get() + m_numberOfElements;
            }

            void Negate()
            {
                for(unsigned int i = 0; i < m_numberOfElements; ++i)
                {
                    m_data[i] = -m_data[i];
                }
            }
            
            void Transpose()
            {
                // A diagonal matrix is its own transpose.
            }

            expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eDiagonal, eNormal, space> > >, expt::NegateOp> > operator-() const
            {
                return expt::Expression<expt::UnaryExpressionPolicy<expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eDiagonal, eNormal, space> > >, expt::NegateOp> >(
                        expt::Expression<expt::ConstantExpressionPolicy<NekMatrix<DataType, eDiagonal, eNormal, space> > >(*this));
            }


            // This is wrong as well.  What if this is diagonal?
            // enable if on the output.
            NekMatrix<DataType, eDiagonal, eNormal, space> operator+=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            {
                ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+");

                for(unsigned int i = 0; i < rhs.GetRows(); ++i)
                {
                    (*this)(i,i) += rhs(i,i);
                }

                return *this;
            }

            // Full *= full = full
            // Full *= diagonal = full

            // diag *= diag = diag
            // diag *= full = full


//             template<NekMatrixForm rhsForm>
//             NekMatrix<DataType, eFull, BlockType, space>& operator*=(const NekMatrix<DataType, rhsForm, BlockType, space>& rhs)
//             {
//                 ASSERTL0(GetColumns() == rhs.GetRows(), "Invalid matrix dimensions in operator*");
            // 
//                 NekMatrix<DataType, eFull, BlockType, space> result(GetRows(), rhs.GetColumns());
            // 
//                 for(unsigned int i = 0; i < result.GetRows(); ++i)
//                 {
//                     for(unsigned int j = 0; j < result.GetColumns(); ++j)
//                     {
//                         DataType t = DataType(0);
            // 
//                         // Set the result(i,j) element.
//                         for(unsigned int k = 0; k < GetColumns(); ++k)
//                         {
//                             t += (*this)(i,k)*rhs(k,j);
//                         }
//                         result(i,j) = t;
//                     }
//                 }
            // 
//                 Swap(result);
            // 
//                 return *this;
//             }

            NekMatrix<DataType, eDiagonal, eNormal, space>& operator*=(const NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
            {
                ASSERTL0(GetRows() == rhs.GetRows() && GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator*=");

                for(unsigned int i = 0; i < rhs.GetRows(); ++i)
                {
                    (*this)(i,i) *= rhs(i,i);
                }

                return *this;
            }
            
                private:
            /// \brief Constructor used for block matrices.  Care must be used with this constructor
            ///         to initialize the matrix before use.
            //NekMatrix();
            
            //void Initialize(unsigned int rows, unsigned int columns);
            
                    void Swap(NekMatrix<DataType, eDiagonal, eNormal, space>& rhs)
                    {
                        std::swap(m_numberOfElements, rhs.m_numberOfElements);
                        std::swap(m_data, rhs.m_data);
                    }

                    unsigned int m_numberOfElements;
                    boost::shared_array<DataType> m_data;

            };*/
            
};
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_DIAGONAL_MATRIX_HPP

/**
    $Log:$
 **/
 
