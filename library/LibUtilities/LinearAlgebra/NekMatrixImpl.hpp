///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrixImpl.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_MATRIX_IMPL_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_MATRIX_IMPL_HPP

namespace Nektar
{
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     class NekMatrix;
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     class NekMatrixImpl;
// 
//     template<typename DataType, MatrixBlockType BlockType, unsigned int space>
//     class NekMatrixImpl<DataType, FullMatrixTag, BlockType, space>
//     {
//         public:
//             typedef typename MatrixDataType<DataType, FullMatrixTag, BlockType, space>::Type BlockDataType;
//             
// 
//         public:
//             //static void CopyMatrixValues(DataType* data, const NekMatrix<DataType, DiagonalMatrixTag, space>& rhs, unsigned int rows, unsigned int columns)
//             static void CopyMatrixValues(BlockDataType* data, const NekMatrix<DataType, DiagonalMatrixTag, BlockType, space>& rhs, unsigned int rows, unsigned int columns)
//             {
//                 std::fill(data, end(data, rows, columns), BlockDataType(0));
//                 for(unsigned int i = 0; i < rows; ++i)
//                 {
//                     SetValue(data, rows, i, i, NekMatrixImpl<BlockDataType, DiagonalMatrixTag, BlockType, space>::GetValue(rhs.begin(), rows, i, i));
//                 }
//             }
// 
//             static inline BlockDataType* CreateMatrixStorage(unsigned int blockRows, unsigned int blockColumns, unsigned int* blockRowSizes, unsigned int* blockColumnSizes)
//             {
//                 BlockDataType* result = CreateMatrixStorage(blockRows, blockColumns);
//                 for(unsigned int i = 0; i < blockRows; ++i)
//                 {
//                     for(unsigned int j = 0; j < blockColumns; ++j)
//                     {
//                         GetValue(result, blockColumns, i, j).Initialize(blockRowSizes[i], blockColumnSizes[j]);
//                     }
//                 }
//                 return result;
//             }
//             
//             static void CopyMatrixValues(BlockDataType* data, const NekMatrix<DataType, FullMatrixTag, BlockType, space>& rhs, unsigned int rows, unsigned int columns)
//             {
//                 std::copy(rhs.begin(), rhs.end(), data);
//             }
// 
//             static void CopyMatrixValues(BlockDataType* dest, const BlockDataType* src, unsigned int rows, unsigned int columns)
//             {
//                 std::copy(src, src + rows*columns, dest);
//             }
//             
//             static inline BlockDataType* CreateMatrixStorage(unsigned int rows, unsigned int columns)
//             {
//                 return Nektar::MemoryManager::AllocateArray<BlockDataType>(rows*columns);
//             }
// 
//             static inline void DeallocateMatrixStorage(BlockDataType*& data, unsigned int rows, unsigned int columns)
//             {
//                 Nektar::MemoryManager::DeallocateArray<BlockDataType>(data, rows*columns);
//             }
// 
//             static inline typename NekMatrix<DataType, FullMatrixTag, BlockType, space>::iterator end(BlockDataType* data, unsigned int rows, unsigned int columns)
//             {
//                 return &data[rows*columns];
//             }
// 
//             static inline typename NekMatrix<DataType, FullMatrixTag, BlockType, space>::iterator begin(BlockDataType* data, unsigned int rows, unsigned int columns)
//             {
//                 return data;
//             }
// 
//             static inline typename boost::call_traits<BlockDataType>::reference GetValue(
//                     BlockDataType* data, unsigned int matrixColumns, unsigned int rowNumber, unsigned int colNumber)
//             {
//                 return data[rowNumber*matrixColumns + colNumber];
//             }
// 
//             static inline typename boost::call_traits<BlockDataType>::reference GetValue(
//                     const BlockDataType* data, unsigned int matrixColumns, unsigned int rowNumber, unsigned int colNumber)
//             {
//                 return data[rowNumber*matrixColumns + colNumber];
//             }
// 
//             static inline void SetValue(BlockDataType* data, unsigned int matrixRows, unsigned int rowNumber,
//                                         unsigned int colNumber, typename boost::call_traits<BlockDataType>::param_type rhs)
//             {
//                 NekMatrixImpl<DataType, FullMatrixTag, BlockType, space>::GetValue(data, matrixRows, rowNumber, colNumber) = rhs;
//             }
// 
//             static inline BlockDataType* GetPtr(BlockDataType* data, unsigned int matrixColumns, unsigned int rowNumber, unsigned int colNumber)
//             {
//                 return &data[rowNumber*matrixColumns + colNumber];
//             }
//             
//             static void Transpose(BlockDataType* data, unsigned int rows, unsigned int columns)
//             {
//                 for(unsigned int row = 0; row < rows; ++row)
//                 {
//                     for(unsigned int column = row+1; column < columns; ++column)
//                     {
//                         std::swap(GetValue(data, columns, row, column), GetValue(data, columns, column, row));
//                     }
//                 }
//             }
//     };
// 
//     template<typename DataType, MatrixBlockType BlockType, unsigned int space>
//     class NekMatrixImpl<DataType, DiagonalMatrixTag, BlockType, space>
//     {
//         public:
//             typedef typename MatrixDataType<DataType, DiagonalMatrixTag, BlockType, space>::Type BlockDataType;
//             
//         public:
//             static inline BlockDataType* CreateMatrixStorage(unsigned int rows, unsigned int columns)
//             {
//                 ASSERTL0(rows == columns, "Digaonal matrices must be square.");
//                 return Nektar::MemoryManager::AllocateArray<BlockDataType>(rows);
//             }
// 
//             static inline BlockDataType* CreateMatrixStorage(unsigned int blockRows, unsigned int blockColumns, unsigned int* blockRowSizes, unsigned int* blockColumnSizes)
//             {
//                 ASSERTL0(blockRows == blockColumns, "Digaonal matrices must be square.");
//                 BlockDataType* result = CreateMatrixStorage(blockRows, blockColumns);
//                 for(unsigned int i = 0; i < blockColumns; ++i)
//                 {
//                     result[i].Initialize(blockRowSizes[i], blockColumnSizes[i]);
//                 }
//                 return result;
//             }
//             
//             static void CopyMatrixValues(BlockDataType* dest, const BlockDataType* src, unsigned int rows, unsigned int columns)
//             {
//                 std::copy(src, src + rows, dest);
//             }
//             
//             static inline void DeallocateMatrixStorage(BlockDataType*& data, unsigned int rows, unsigned int /*columns*/)
//             {
//                 Nektar::MemoryManager::DeallocateArray<BlockDataType>(data, rows);
//             }
// 
//             static inline typename NekMatrix<DataType, DiagonalMatrixTag, BlockType, space>::iterator end(BlockDataType* data, unsigned int rows, unsigned int /*columns*/)
//             {
//                 return &data[rows];
//             }
// 
//             static inline typename boost::call_traits<BlockDataType>::reference GetValue(
//                     BlockDataType* data, unsigned int /*matrixColumns*/, unsigned int rowNumber, unsigned int /*columnNumber*/)
//             {
//                 return data[rowNumber];
//             }
// 
//             static inline typename boost::call_traits<BlockDataType>::const_reference GetValue(
//                     const BlockDataType* data, unsigned int /*matrixColumns*/, unsigned int rowNumber, unsigned int /*columnNumber*/)
//             {
//                 return data[rowNumber];
//             }
// 
//             static inline void SetValue(BlockDataType* data, unsigned int matrixRows, unsigned int rowNumber,
//                                         unsigned int colNumber, typename boost::call_traits<BlockDataType>::param_type rhs)
//             {
//                 ASSERTL0(rowNumber==colNumber, "Can only SetValue on the diagonal for diagonal matrices.") 
//                 NekMatrixImpl<DataType, DiagonalMatrixTag, BlockType, space>::GetValue(data, matrixRows, rowNumber, colNumber) = rhs;
//             }
// 
//             static inline BlockDataType* GetPtr(BlockDataType* data, unsigned int /*matrixColumns*/, unsigned int rowNumber, unsigned int /*columnNumber*/)
//             {
//                 return &data[rowNumber];
//             }
//     };
    
//     dtemplate<typename InnerDataType, typename InnerForm, typename InnerSpace, NekMatrixForm form, unsigned int space>
//     class NekMatrixImpl<NekMatrix<InnerDataType, InnerForm, InnerSpace>, form, space>
//     {
//         static inline DataType* CreateMatrixStorage(unsigned int rows, unsigned int columns)
//         {
//             
//             return Nektar::MemoryManager::AllocateArray<DataType>(rows);
//         }
// 
//         static inline void DeallocateMatrixStorage(DataType*& data, unsigned int rows, unsigned int /*columns*/)
//         {
//             Nektar::MemoryManager::DeallocateArray<DataType>(data, rows);
//         }
// 
//         static inline typename NekMatrix<DataType, DiagonalMatrixTag, space>::iterator end(DataType* data, unsigned int /*rows*/, unsigned int columns)
//         {
//             return &data[columns];
//         }
// 
//         static inline typename boost::call_traits<DataType>::reference GetValue(
//         DataType* data, unsigned int /*matrixColumns*/, unsigned int /*rowNumber*/, unsigned int columnNumber)
//         {
//             return data[columnNumber];
//         }
// 
//         static inline typename boost::call_traits<DataType>::const_reference GetValue(
//             const DataType* data, unsigned int /*matrixColumns*/, unsigned int /*rowNumber*/, unsigned int columnNumber)
//         {
//             return data[columnNumber];
//         }
// 
//         static inline void SetValue(DataType* data, unsigned int matrixRows, unsigned int rowNumber,
//             unsigned int colNumber, typename boost::call_traits<DataType>::param_type rhs)
//         {
//             ASSERTL0(rowNumber==colNumber, "Can only SetValue on the diagonal for diagonal matrices.") 
//                     NekMatrixImpl<DataType, DiagonalMatrixTag, space>::GetValue(data, matrixRows, rowNumber, colNumber) = rhs;
//         }
// 
//         static inline DataType* GetPtr(DataType* data, unsigned int /*matrixColumns*/, unsigned int /*rowNumber*/, unsigned int columnNumber)
//         {
//             return &data[columnNumber];
//         }
//     };

}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_MATRIX_IMPL_HPP


/**
    $Log: NekMatrixImpl.hpp,v $
    Revision 1.4  2007/07/27 03:15:30  bnelson
    Removed MatrixBlockType.h

    Revision 1.3  2007/06/10 23:42:16  bnelson
    Matrix updates.

    Revision 1.2  2007/03/29 18:59:05  bnelson
    Refactoring in preparation for scaled matrices.  Fixed transpose problem.

    Revision 1.1  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.


**/
