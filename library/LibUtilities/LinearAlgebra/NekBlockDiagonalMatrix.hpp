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
                m_data(m_numberOfElements)
            {
                ASSERTL0(blockRows == blockColumns, "ERROR: Block diagonal matrices must consist of square blocks.");
                for(unsigned int i = 0; i < m_numberOfElements; ++i)
                {
                    m_data[i]->Initialize(blockRows, blockColumns);
                }
            }
            
            NekMatrix(const NekMatrix<DataType, eDiagonal, BlockType, space>& rhs) :
                m_numberOfElements(rhs.m_numberOfElements),
                m_rows(rhs.m_rows),
                m_columns(rhs.m_columns),
                m_blockColumns(rhs.m_blockColumns),
                m_blockRows(rhs.m_blockRows),
                m_data(m_numberOfElements)
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

            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            template<typename ExpressionPolicyType>
            NekMatrix(const expt::Expression<ExpressionPolicyType>& rhs) :
                m_numberOfElements(rhs.GetMetadata().Rows),
                m_rows(rhs.GetMetadata().Rows),
                m_columns(rhs.GetMetadata().Columns),
                m_blockColumns(rhs.GetMetadata().BlockColumns),
                m_blockRows(rhs.GetMetadata().BlockRows),
                m_data(m_numberOfElements)
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
                m_data = Array<OneD, InnerDataType>(m_numberOfElements);

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
                
                if( blockRow == blockColumn )
                {
                    ASSERTL2(rowNumber < m_numberOfElements, "Illegal access to NekMatrix via operator()");
                    return GetBlock(blockRow, blockColumn)->operator()(innerBlockRow, innerBlockColumn);
                }
                else
                {
                    static InnerDataType zeroElement(m_blockRows[0], m_blockColumns[0], DataType(0));
                    return zeroElement->operator()(innerBlockRow, innerBlockColumn);
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
                    return GetBlock(blockRow, blockColumn)->operator()(innerBlockRow, innerBlockColumn);
                }
                else
                {
                    static InnerDataType zeroElement(m_blockRows[0], m_blockColumns[0], DataType(0));
                    return zeroElement->operator()(innerBlockRow, innerBlockColumn);
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
//             Array<OneD, DataType> GetPtr()
//             {
//                 SharedArray<DataType> result(new DataType[GetRows()*GetColumns()]);
//                 //std::fill(result, result+GetRows()*GetColumns(), DataType(0));
//                 
//                 for(unsigned int i = 0; i < GetRows(); ++i)
//                 {
//                     for(unsigned int j = 0; j < GetColumns(); ++j)
//                     {
//                         result[i*GetColumns() + j] = (*this)(i,j);
//                     }
//                 }
//                 
//                 return result;
//             }
            
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
            
            Array<OneD, InnerDataType> m_data;
            
    };
};
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_DIAGONAL_MATRIX_HPP

/**
    $Log: NekBlockDiagonalMatrix.hpp,v $
    Revision 1.4  2007/03/31 15:38:45  bnelson
    *** empty log message ***

    Revision 1.3  2007/01/23 03:12:49  jfrazier
    Added more conditional compilation directives for expression templates.

    Revision 1.2  2006/11/01 04:07:07  bnelson
    Changed block matrices to use the ConsistentObjectAccess object to store matrices or pointers to matrices so that the same pointer syntax works for both.

    Revision 1.1  2006/10/30 05:11:15  bnelson
    Added preliminary linear system and block matrix support.

 **/
 
