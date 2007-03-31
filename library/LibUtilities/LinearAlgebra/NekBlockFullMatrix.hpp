///////////////////////////////////////////////////////////////////////////////
//
// File: NekBlockFullMatrix.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_FULL_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_FULL_MATRIX_HPP

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
    class NekMatrix<DataType, eFull, BlockType, space, typename boost::disable_if_c<BlockType==eNormal, void>::type >
    {
        public:
            typedef NekMatrix<DataType, eFull, BlockType, space, void> ThisType;
            typedef NekMatrix<DataType, eFull, eNormal, space, void> InnerMatrixType;
            typedef typename BlockMatrixDataType<BlockType, NekMatrix<DataType, eFull, eNormal, space, void> >::ResultType InnerDataType;
            
            
        public:
            NekMatrix(unsigned int rows, unsigned int columns, unsigned int blockRows, unsigned int blockColumns) :
                m_rows(rows),
                m_columns(columns),
                m_blockColumns(m_rows, blockColumns),
                m_blockRows(m_columns, blockRows),
                m_data(MemoryManager::AllocateSharedArray<InnerDataType>(m_rows*m_columns))
            {
                ASSERTL0(blockRows%rows == 0, "ERROR: Blocks in a block matrix must fit evenly in the requested size.");
                ASSERTL0(blockColumns%columns == 0, "ERROR: Blocks in a block matrix must fit evenly in the requested size.");

                for(unsigned int i = 0; i < m_rows; ++i)
                {
                    for(unsigned int j = 0; j < m_rows; ++j)
                    {
                        GetObj(GetBlock(i,j))->Initialize(blockRows, blockColumns);
                    }
                }
            }
            
            NekMatrix(const NekMatrix<DataType, eFull, BlockType, space>& rhs) :
                m_rows(rhs.m_rows),
                m_columns(rhs.m_columns),
                m_blockColumns(rhs.m_blockColumns),
                m_blockRows(rhs.m_blockRows),
                m_data(MemoryManager::AllocateSharedArray<InnerDataType>(m_rows*m_columns))
            {
                for(unsigned int i = 0; i < m_rows; ++i)
                {
                    for(unsigned int j = 0; j < m_rows; ++j)
                    {
                        GetBlock(i,j) = rhs.GetBlock(i,j);
                    }
                }
           }
                
            NekMatrix<DataType, eFull, BlockType, space>& operator=(const NekMatrix<DataType, eFull, BlockType, space>& rhs)
            {
                NekMatrix<DataType, eFull, BlockType, space> temp(rhs);
                Swap(temp);
                return *this;
            }

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

            template<typename ExpressionPolicyType>
            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                m_rows(rhs.GetMetadata().Rows),
                m_columns(rhs.GetMetadata().Columns),
                m_blockColumns(rhs.GetMetadata().BlockColumns),
                m_blockRows(rhs.GetMetadata().BlockRows),
                m_data(MemoryManager::AllocateSharedArray<InnerDataType>(m_rows*m_columns))
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eFull, BlockType, space> > ));
                rhs.Apply(*this);
            }
            
            template<typename ExpressionPolicyType>
            NekMatrix<DataType, eFull, eNormal, space>& operator=(const expt::Expression<ExpressionPolicyType>& rhs)
            {
                BOOST_MPL_ASSERT(( boost::is_same<typename expt::Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, eFull, BlockType, space> > ));
                
                m_rows = rhs.GetMetadata().Rows;
                m_columns = rhs.GetMetadata().Columns;
                m_blockColumns = rhs.GetMetadata().BlockColumns;
                m_blockRows = rhs.GetMetadata().BlockRows;
                m_data = MemoryManager::AllocateSharedArray<InnerDataType>(m_rows*m_columns);

                rhs.Apply(*this);
                return *this;
            }
#endif              
            
            typename boost::call_traits<DataType>::reference operator()(unsigned int rowNumber, unsigned int colNumber)
            {
                unsigned int blockRow = rowNumber/m_blockRows[0];
                unsigned int blockColumn = colNumber/m_blockColumns[0];
                unsigned int innerBlockRow = rowNumber%m_blockRows[0];
                unsigned int innerBlockColumn = colNumber%m_blockColumns[0];
                

                return GetBlock(blockRow, blockColumn)(innerBlockRow, innerBlockColumn);
            }

            typename boost::call_traits<DataType>::const_reference operator()(unsigned int rowNumber, unsigned int colNumber) const
            {
                unsigned int blockRow = rowNumber/m_blockRows[0];
                unsigned int blockColumn = colNumber/m_blockColumns[0];
                unsigned int innerBlockRow = rowNumber%m_blockRows[0];
                unsigned int innerBlockColumn = colNumber%m_blockColumns[0];
                

                return GetBlock(blockRow, blockColumn)(innerBlockRow, innerBlockColumn);
            }

            typename boost::call_traits<InnerDataType>::reference GetBlock(unsigned int rowNumber, unsigned int colNumber)
            {
                return m_data[rowNumber*m_columns + colNumber];
            }

            typename boost::call_traits<InnerDataType>::const_reference GetBlock(unsigned int rowNumber, unsigned int colNumber) const
            {
                return m_data[rowNumber*m_columns + colNumber];
            }
            
            unsigned int GetRows() const { return m_rows; }
            unsigned int GetColumns() const { return m_columns; }
            const std::vector<unsigned int> GetBlockRows() const { return m_blockRows; }
            const std::vector<unsigned int> GetBlockColumns() const { return m_blockColumns; }

            /// \brief Return a full matrix version.
            SharedArray<DataType> GetPtr()
            {
                SharedArray<DataType> result(new DataType[GetRows()*GetColumns()]);
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
                for(unsigned int i = 0; i < m_rows; ++i)
                {
                    for(unsigned int j = 0; j < m_columns; ++j)
                    {
                        GetBlock(i,j) += rhs.GetBlock(i,j);
                    }
                }
                
                return *this;
            }

            
        private:
            
            void Swap(NekMatrix<DataType, eDiagonal, BlockType, space>& rhs)
            {
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
    
            SharedArray<InnerDataType> m_data;
            
    };
};
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_FULL_MATRIX_HPP

/**
    $Log: NekBlockFullMatrix.hpp,v $
    Revision 1.2  2007/01/23 03:12:49  jfrazier
    Added more conditional compilation directives for expression templates.

    Revision 1.1  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

 **/
 
