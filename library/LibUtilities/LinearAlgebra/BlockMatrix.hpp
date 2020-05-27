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

#include <memory>
#include <tuple>

#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/MatrixFuncs.h>

namespace Nektar
{
    template<typename DataType, typename InnerMatrixType>
    class NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag> : public ConstMatrix<typename NekMatrix<DataType, InnerMatrixType>::NumberType>
    {
        public:
            typedef NekMatrix<DataType, InnerMatrixType> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> ThisType;
            typedef typename InnerType::NumberType NumberType;
            typedef ConstMatrix<NumberType> BaseType;

            // Each inner matrix type can possible return references or value types from GetValue.
            // Query the type here to find out.

            typedef typename InnerType::GetValueType GetValueType;
            typedef typename InnerType::ConstGetValueType ConstGetValueType;

        public:
            template<typename MatrixType>
            class iterator_base
            {
                public:
                    typedef typename MatrixType::InnerType IteratorInnerType;
                    typedef FullMatrixFuncs StoragePolicy;

                public:
                    iterator_base(MatrixType& m, unsigned int curRow, unsigned int curCol) :
                        m_matrix(m),
                        m_curRow(curRow),
                        m_curColumn(curCol)
                    {
                    }

                    iterator_base(MatrixType& m) :
                        m_matrix(m),
                        m_curRow(std::numeric_limits<unsigned int>::max()),
                        m_curColumn(std::numeric_limits<unsigned int>::max())
                    {
                    }

                    iterator_base(const iterator_base<MatrixType>& rhs) :
                        m_matrix(rhs.m_matrix),
                        m_curRow(rhs.m_curRow),
                        m_curColumn(rhs.m_curColumn)
                    {
                    }

                    void operator++()
                    {
                        if( m_curRow != std::numeric_limits<unsigned int>::max() )
                        {
                            std::tie(m_curRow, m_curColumn) = StoragePolicy::Advance(
                                m_matrix.GetRows(), m_matrix.GetColumns(), m_curRow, m_curColumn);
                        }
                    }

                    NumberType operator*()
                    {
                        return m_matrix(m_curRow, m_curColumn);
                    }

                    bool operator==(const iterator_base<MatrixType>& rhs)
                    {
                        return m_curRow == rhs.m_curRow && m_curColumn == rhs.m_curColumn;
                    }

                    bool operator!=(const iterator_base<MatrixType>& rhs)
                    {
                        return !(*this == rhs);
                    }

                private:
                    iterator_base<MatrixType>& operator=(const iterator_base<MatrixType>& rhs);

                    MatrixType& m_matrix;
                    unsigned int m_curRow;
                    unsigned int m_curColumn;

            };

            typedef iterator_base<ThisType> iterator;
            typedef iterator_base<const ThisType> const_iterator;

        public:
            LIB_UTILITIES_EXPORT explicit NekMatrix(MatrixStorage type = eFULL);

            LIB_UTILITIES_EXPORT NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      unsigned int rowsPerBlock, unsigned int columnsPerBlock,
                      MatrixStorage type = eFULL);

            LIB_UTILITIES_EXPORT NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      const unsigned int* rowsPerBlock, const unsigned int* columnsPerBlock,
                      MatrixStorage type = eFULL);

            LIB_UTILITIES_EXPORT NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      const Array<OneD, const unsigned int>& rowsPerBlock, const Array<OneD, const unsigned int>& columnsPerBlock,
                      MatrixStorage type = eFULL);

            LIB_UTILITIES_EXPORT NekMatrix(const Array<OneD, const unsigned int>& rowsPerBlock,
                      const Array<OneD, const unsigned int>& columnsPerBlock,
                      MatrixStorage type = eFULL);

            LIB_UTILITIES_EXPORT NekMatrix(const ThisType& rhs) ;

            LIB_UTILITIES_EXPORT unsigned int GetRequiredStorageSize() const;

            LIB_UTILITIES_EXPORT unsigned int CalculateBlockIndex(unsigned int row, unsigned int column) const;

            LIB_UTILITIES_EXPORT const InnerType* GetBlockPtr(unsigned int row, unsigned int column) const;

            LIB_UTILITIES_EXPORT std::shared_ptr<const InnerType> GetBlock(unsigned int row, unsigned int column) const;

            LIB_UTILITIES_EXPORT std::shared_ptr<InnerType>& GetBlock(unsigned int row, unsigned int column);

            LIB_UTILITIES_EXPORT void SetBlock(unsigned int row, unsigned int column, std::shared_ptr<InnerType>& m);

            LIB_UTILITIES_EXPORT ConstGetValueType operator()(unsigned int row, unsigned int col) const;

            LIB_UTILITIES_EXPORT unsigned int GetStorageSize() const;

            LIB_UTILITIES_EXPORT unsigned int GetNumberOfBlockRows() const;

            LIB_UTILITIES_EXPORT unsigned int GetNumberOfBlockColumns() const;

            LIB_UTILITIES_EXPORT unsigned int GetNumberOfRowsInBlockRow(unsigned int blockRow) const;

            LIB_UTILITIES_EXPORT unsigned int GetNumberOfColumnsInBlockColumn(unsigned int blockCol) const;

            LIB_UTILITIES_EXPORT void GetBlockSizes(Array<OneD, unsigned int>& rowSizes,
                                                    Array<OneD, unsigned int>& colSizes) const;

            LIB_UTILITIES_EXPORT iterator begin();
            LIB_UTILITIES_EXPORT iterator end();
            LIB_UTILITIES_EXPORT const_iterator begin() const;
            LIB_UTILITIES_EXPORT const_iterator end() const;

        public:
            LIB_UTILITIES_EXPORT static ThisType CreateWrapper(const ThisType& rhs);

            LIB_UTILITIES_EXPORT static std::shared_ptr<ThisType> CreateWrapper(const std::shared_ptr<ThisType>& rhs);

        private:

            LIB_UTILITIES_EXPORT static unsigned int GetNumberOfElementsInBlock(unsigned int block, unsigned int totalBlocks, const Array<OneD, unsigned int>& sizes);

            LIB_UTILITIES_EXPORT void Initialize(const unsigned int* rowsPerBlock, const unsigned int* columnsPerBlock);

            LIB_UTILITIES_EXPORT virtual typename boost::call_traits<NumberType>::value_type v_GetValue(unsigned int row, unsigned int column) const;

            LIB_UTILITIES_EXPORT virtual unsigned int v_GetStorageSize() const;

            LIB_UTILITIES_EXPORT virtual void v_Transpose();

            Array<OneD, std::shared_ptr<InnerType> > m_data;
            std::shared_ptr<InnerType> m_nullBlockPtr;
            Array<OneD, unsigned int> m_rowSizes;
            Array<OneD, unsigned int> m_columnSizes;
            unsigned int m_storageSize;
            unsigned int m_numberOfBlockRows;
            unsigned int m_numberOfBlockColumns;
            static NumberType m_zeroElement;
    };

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NumberType
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::m_zeroElement(0);

    template<typename InnerMatrixType>
    NekMatrix<InnerMatrixType, BlockMatrixTag>
    Transpose(NekMatrix<InnerMatrixType, BlockMatrixTag>& rhs)
    {
        NekMatrix<InnerMatrixType, BlockMatrixTag> result(rhs);
        result.Transpose();
        return result;
    }
}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP
