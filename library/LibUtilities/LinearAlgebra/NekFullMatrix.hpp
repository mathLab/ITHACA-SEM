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
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <boost/call_traits.hpp>

#include <algorithm>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif


namespace Nektar
{
//     template<typename DataType>
//     class NekMatrixStoragePolicy<DataType, FullMatrixTag> 
//     {
//         public: 
//             static void Transpose(unsigned int& rows, unsigned int& columns,
//                                   Array<OneD, DataType>& data)
//             {
//                 Array<OneD, DataType> temp(data.GetSize());
// 
//                 for(unsigned int row = 0; row < rows; ++row)
//                 {
//                     for(unsigned int column = 0; column < columns; ++column)
//                     {
//                         unsigned int firstIndex = CalculateIndex(row, column, rows, columns);
//                         unsigned int secondIndex = CalculateIndex(column, row, columns, rows);
// 
//                         temp[secondIndex] = data[firstIndex];
//                     }
//                 }
// 
//                 std::swap(rows, columns);
//                 std::swap(data, temp);
//             }
// 
// 
//             static typename boost::call_traits<DataType>::reference GetData(unsigned int row, unsigned int column, unsigned int matrixRows, unsigned int matrixColumns,
//                                                                             Array<OneD, DataType>& data)
//             {
//                 return data[CalculateIndex(row, column, matrixRows, matrixColumns)];
//             }
// 
//             static typename boost::call_traits<DataType>::const_reference GetConstData(unsigned int row, unsigned int column, unsigned int matrixRows, unsigned int matrixColumns,
//                                                                             const ConstArray<OneD, DataType>& data)
//             {
//                 return data[CalculateIndex(row, column, matrixRows, matrixColumns)];
//             }
// 
//             static unsigned int CalculateIndex(unsigned int row, unsigned int column, unsigned int matrixRows, unsigned int matrixColumns)
//             {
//                 return row*matrixColumns + column;
//             }
// 
//             static unsigned int NumStorageElements(unsigned int rows, unsigned int cols)
//             {
//                 return rows*cols;
//             }
// 
//             static void Invert(unsigned int rows, unsigned int columns, Array<OneD, DataType>& data)
//             {
// #ifdef NEKTAR_USING_LAPACK
//                 ASSERTL0(rows == columns, "Matrix Inversion only works for square arrays.");
// 
//                 /// Incoming data is row major, make it column major for lapack calls.
//                 Transpose(rows, columns, data);
// 
//                 int m = rows;
//                 int n = columns;
//                 int pivotSize = std::max(1, std::min(m, n));
// 
//                 Array<OneD, int> ipivot(pivotSize);
//                 int info = 0;
//                 Lapack::Dgetrf(m, n, data.get(), m, ipivot.get(), info);
// 
//                 if( info < 0 )
//                 {
//                     std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + 
//                             "th parameter had an illegal parameter for dgetrf";
//                     ASSERTL0(false, message.c_str());
//                 }
//                 else if( info > 0 )
//                 {
//                     std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
//                             boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
//                     ASSERTL0(false, message.c_str());
//                 }
// 
//                 unsigned int workSize = 64*n;
// 
//                 Array<OneD, NekDouble> work(workSize);
//                 Lapack::Dgetri(n, data.get(), n, ipivot.get(), work.get(), workSize, info);
// 
//                 if( info < 0 )
//                 {
//                     std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + 
//                             "th parameter had an illegal parameter for dgetri";
//                     ASSERTL0(false, message.c_str());
//                 }
//                 else if( info > 0 )
//                 {
//                     std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
//                             boost::lexical_cast<std::string>(info) + " is 0 from dgetri";
//                     ASSERTL0(false, message.c_str());
//                 }
// 
//                 // Put it back to row major form.
//                 Transpose(rows, columns, data);
// 
// #else
//                 // TODO
//                 BOOST_STATIC_ASSERT(sizeof(DataType) == 0);
// #endif //NEKTAR_USING_LAPACK
// 
//             }
// 
//             static void Factorize(unsigned int rows, unsigned int columns, Array<OneD, DataType>& data)
//             {
// 
//             }
// 
//         private:
//             
//     };
// 
//     template<typename DataType>
//     class NekMatrixArithmeticPolicy<DataType, FullMatrixTag>
//     {
//         public:
//             template<MatrixBlockType BlockType, unsigned int space>
//             static void PlusEqual(NekMatrix<DataType, FullMatrixTag, BlockType, space>& lhs, 
//                                   const NekMatrix<DataType, FullMatrixTag, BlockType, space>& rhs)
//             {
//                 ASSERTL0(lhs.GetRows() == rhs.GetRows() && lhs.GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");
//                 typename NekMatrix<DataType, FullMatrixTag, BlockType, space>::iterator lhs_data = lhs.begin();
//                 typename NekMatrix<DataType, FullMatrixTag, BlockType, space>::const_iterator rhs_data = rhs.begin();                
// 
//                 for( ; lhs_data < lhs.end(); ++lhs_data, ++rhs_data )
//                 {
//                     *lhs_data += *rhs_data;
//                 }
//             }
// 
//             template<MatrixBlockType BlockType, unsigned int space>
//             static void PlusEqual(NekMatrix<DataType, FullMatrixTag, BlockType, space>& lhs, 
//                                   const NekMatrix<DataType, DiagonalMatrixTag, BlockType, space>& rhs)
//             {
//                 ASSERTL0(lhs.GetRows() == rhs.GetRows() && lhs.GetColumns() == rhs.GetColumns(), "Matrix dimensions must agree in operator+=");
// 
//                 for(unsigned int i = 0; i < lhs.GetRows(); ++i)
//                 {
//                     lhs(i,i) += rhs(i,i);
//                 }
//             }
//     };
// 
//     template<typename DataType>
//     class NekMatrixAssignmentPolicy<DataType, FullMatrixTag>
//     {
//         public:
//             template<MatrixBlockType BlockType, unsigned int space>
//             static void Assign(NekMatrix<DataType, FullMatrixTag, BlockType, space>& lhs, 
//                                const NekMatrix<DataType, DiagonalMatrixTag, BlockType, space>& rhs)
//             {
//                 lhs = NekMatrix<DataType, FullMatrixTag, BlockType, space>(rhs.GetRows(), rhs.GetColumns(), DataType(0));
//                 for(unsigned int i = 0; i < rhs.GetRows(); ++i)
//                 {
//                     lhs(i,i) = rhs(i,i);
//                 }
//             }
//     };

}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_FULL_MATRIX_HPP

/**
    $Log: NekFullMatrix.hpp,v $
    Revision 1.15  2007/05/15 00:14:10  bnelson
    Updated to use the new Array object.

    Revision 1.14  2007/04/21 04:09:43  bnelson
    *** empty log message ***

    Revision 1.13  2007/04/10 16:43:49  bnelson
    *** empty log message ***

    Revision 1.12  2007/04/09 03:16:53  bnelson
    *** empty log message ***

    Revision 1.11  2007/04/05 05:12:44  bnelson
    *** empty log message ***

    Revision 1.10  2007/04/04 02:26:41  bnelson
    *** empty log message ***

    Revision 1.9  2007/04/04 02:11:07  bnelson
    Added inversion

    Revision 1.8  2007/03/29 18:59:05  bnelson
    Refactoring in preparation for scaled matrices.  Fixed transpose problem.

    Revision 1.7  2007/02/15 06:56:54  bnelson
    *** empty log message ***

    Revision 1.6  2007/02/04 04:27:43  bnelson
    Updated linear systems to work with normal objects instead of only shared pointers.

    Revision 1.5  2007/01/29 01:30:20  bnelson
    Removed memory manager requirements.

    Revision 1.4  2007/01/23 03:12:50  jfrazier
    Added more conditional compilation directives for expression templates.

    Revision 1.3  2006/11/08 04:16:14  bnelson
    Added subtraction operators.

    Revision 1.2  2006/11/01 04:07:08  bnelson
    Changed block matrices to use the ConsistentObjectAccess object to store matrices or pointers to matrices so that the same pointer syntax works for both.

    Revision 1.1  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

 **/

