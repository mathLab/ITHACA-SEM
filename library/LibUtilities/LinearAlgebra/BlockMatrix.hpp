///////////////////////////////////////////////////////////////////////////////
//
// File: BlockMatrix.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/BasicUtils/BinaryExpressionTraits.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/FullMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/DiagonalMatrixStoragePolicy.hpp>

#include <boost/shared_ptr.hpp>

namespace Nektar
{
    template<typename DataType, typename InnerStorageType, typename InnerMatrixType, typename OuterStorageType>
    class NekMatrix<NekMatrix<DataType, InnerStorageType, InnerMatrixType>, OuterStorageType, BlockMatrixTag> : public Matrix<DataType>
    {
        public:
            typedef Matrix<DataType> BaseType;
            typedef NekMatrix<DataType, InnerStorageType, InnerMatrixType> InnerType;
            typedef NekMatrix<InnerType, OuterStorageType, BlockMatrixTag> ThisType;
            
            // Each inner matrix type can possible return references or value types from GetValue.
            // Query the type here to find out.
            typedef typename InnerType::GetValueType GetValueType;
            typedef typename InnerType::ConstGetValueType ConstGetValueType;
            
            typedef MatrixStoragePolicy<boost::shared_ptr<InnerType>, OuterStorageType> StoragePolicy;
                
        public:
            NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      unsigned int rowsPerBlock, unsigned int columnsPerBlock) :
                BaseType(numberOfBlockRows*rowsPerBlock, numberOfBlockColumns*columnsPerBlock),
                m_data(StoragePolicy::Initialize(numberOfBlockRows, numberOfBlockColumns, boost::shared_ptr<InnerType>())),
                m_rowSizes(numberOfBlockRows+1),
                m_columnSizes(numberOfBlockColumns+1),
                m_storageSize(this->GetRows()*this->GetColumns()),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns)
            {
                for(unsigned int i = 1; i <= numberOfBlockRows; ++i)
                {
                    m_rowSizes[i-1] = i*rowsPerBlock-1;
                }
                
                for(unsigned int i = 1; i <= numberOfBlockColumns; ++i)
                {
                    m_columnSizes[i-1] = i*columnsPerBlock-1;
                }
            }
            
            NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      unsigned int* rowsPerBlock, unsigned int* columnsPerBlock) :
                BaseType(std::accumulate(rowsPerBlock, rowsPerBlock + numberOfBlockRows, 0),
                         std::accumulate(columnsPerBlock, columnsPerBlock + numberOfBlockColumns, 0)),
                m_data(StoragePolicy::Initialize(numberOfBlockRows, numberOfBlockColumns, boost::shared_ptr<InnerType>())),
                m_rowSizes(numberOfBlockRows+1),
                m_columnSizes(numberOfBlockColumns+1),
                m_storageSize(this->GetRows()*this->GetColumns()),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns)
            {
                m_rowSizes[0] = rowsPerBlock[0] - 1;
                for(unsigned int i = 1; i < numberOfBlockRows; ++i)
                {
                    m_rowSizes[i] = rowsPerBlock[i] + m_rowSizes[i-1];
                }
                
                m_columnSizes[0] = columnsPerBlock[0] - 1;
                for(unsigned int i = 1; i < numberOfBlockColumns; ++i)
                {
                    m_columnSizes[i] = columnsPerBlock[i] + m_columnSizes[i-1];
                }
            }
                
            boost::shared_ptr<const InnerType> GetBlock(unsigned int row, unsigned int column) const        
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return StoragePolicy::GetValue(m_numberOfBlockRows, m_numberOfBlockColumns, row, column, m_data);
            }
            
            boost::shared_ptr<InnerType> GetBlock(unsigned int row, unsigned int column)       
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return StoragePolicy::GetValue(m_numberOfBlockRows, m_numberOfBlockColumns, row, column, m_data);
            }
            
            void SetBlock(unsigned int row, unsigned int column, boost::shared_ptr<InnerType>& m)
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                StoragePolicy::SetValue(m_numberOfBlockRows, m_numberOfBlockColumns, row, column, m_data, m);
            }
            
            
            
            
            ConstGetValueType operator()(unsigned int row, unsigned int col) const
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                unsigned int blockRow = std::lower_bound(m_rowSizes.begin(), m_rowSizes.end(), row) - m_rowSizes.begin();
                unsigned int blockColumn = std::lower_bound(m_columnSizes.begin(), m_columnSizes.end(), col) - m_columnSizes.begin();
                unsigned int actualRow = row-(m_rowSizes[blockRow]-1);
                unsigned int actualCol = col-(m_columnSizes[blockColumn]-1);
                
                ASSERTL2(GetBlock(blockRow, blockColumn), std::string("Attempting to access block (") +
                    boost::lexical_cast<std::string>(blockRow) + std::string(", ") + 
                    boost::lexical_cast<std::string>(blockColumn) + std::string(") of a block matrix but it is null."));
                    
                return GetBlock(blockRow, blockColumn)->operator()(actualRow, actualCol);
            }
            
            GetValueType operator()(unsigned int row, unsigned int col)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                unsigned int blockRow = std::lower_bound(m_rowSizes.begin(), m_rowSizes.end(), row) - m_rowSizes.begin();
                unsigned int blockColumn = std::lower_bound(m_columnSizes.begin(), m_columnSizes.end(), col) - m_columnSizes.begin();
                unsigned int actualRow = row-(m_rowSizes[blockRow]-1);
                unsigned int actualCol = col-(m_columnSizes[blockColumn]-1);
                
                ASSERTL2(GetBlock(blockRow, blockColumn), std::string("Attempting to access block (") +
                    boost::lexical_cast<std::string>(blockRow) + std::string(", ") + 
                    boost::lexical_cast<std::string>(blockColumn) + std::string(") of a block matrix but it is null."));
                    
                return GetBlock(blockRow, blockColumn)->operator()(actualRow, actualCol);
            }
            
            unsigned int GetStorageSize() const 
            {
                return m_storageSize;
            }
             
            MatrixStorage GetStorageType() const
            {
                return static_cast<MatrixStorage>(ConvertToMatrixStorageEnum<OuterStorageType>::Value);
            }            
            
            unsigned int GetNumberOfBlockRows() const { return m_numberOfBlockRows; }
            unsigned int GetNumberOfBlockColumns() const { return m_numberOfBlockColumns; }
            
        public:
        
        private:
            virtual typename boost::call_traits<DataType>::value_type v_GetValue(unsigned int row, unsigned int column) const 
            {
                return (*this)(row, column);
            }
            
            virtual unsigned int v_GetStorageSize() const 
            {
                return this->GetStorageSize();
            }
            
            virtual MatrixStorage v_GetStorageType() const
            {
                return this->GetStorageType();
            }
            
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                (*this)(row, column) = d;
            }
            
            Array<OneD, boost::shared_ptr<InnerType> > m_data;
            Array<OneD, unsigned int> m_rowSizes;
            Array<OneD, unsigned int> m_columnSizes;
            unsigned int m_storageSize;
            unsigned int m_numberOfBlockRows;
            unsigned int m_numberOfBlockColumns;            
    };
    
}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP
