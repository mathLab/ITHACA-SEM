///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixTraits.hpp
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
// Description: Interface classes for matrices
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_TRAITS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_TRAITS_HPP

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/BasicUtils/BinaryExpressionTraits.hpp>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace Nektar
{
    // Define the expression traits for a matrix before defining the matrix or the 
    // OperatorGenerator won't work.
    
    /// \brief Gives all expressions involving the same storage types the same result type.
//     template<typename DataType, typename LhsStorageType, typename RhsStorageType, typename OpType>
//     class BinaryExpressionTraits<NekMatrix<DataType, LhsStorageType, StandardMatrixTag>, NekMatrix<DataType, RhsStorageType, StandardMatrixTag>, OpType>
//     {
//         public:
//             typedef NekMatrix<DataType, FullMatrixTag, StandardMatrixTag> ResultType;
//     };
    
    // TODO - Temporary until all the traits are defined, then we can take this one out.
//     template<typename DataType, typename StorageType, typename OpType>
//     class BinaryExpressionTraits<NekMatrix<DataType, StorageType, StandardMatrixTag>, NekMatrix<DataType, StorageType, StandardMatrixTag>, OpType>
//     {
//         public:
//             typedef NekMatrix<DataType, StorageType, StandardMatrixTag> ResultType;
//     };

    template<typename LhsDataType, typename RhsDataType, typename LhsStorageType, typename LhsMatrixType>
    class BinaryExpressionTraits<NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>, RhsDataType, MultiplyOp>
    {
        public:
            typedef typename NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>::NumberType NumberType;
            typedef NekMatrix<NumberType, LhsStorageType, StandardMatrixTag> ResultType;
    };
                        
    template<typename LhsDataType, typename RhsDataType, typename RhsStorageType, typename RhsMatrixType>
    class BinaryExpressionTraits<LhsDataType, NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>, MultiplyOp>
    {
        public:
            typedef typename NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>::NumberType NumberType;
            typedef NekMatrix<NumberType, RhsStorageType, StandardMatrixTag> ResultType;
    };
    
    // All row references are for the spreadsheet "Matrix Operation and ResultTypes".
    
    // The following specializations are for rows 2-10.  These are for plus and minus on the left side of the expression.
    template<typename LhsDataType, typename RhsDataType, typename RhsStorageType, typename LhsMatrixType, typename RhsMatrixType>
    class BinaryExpressionTraits<NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>, NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>, AddOp>
    {
        public:
            typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
            typedef NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> ResultType;
    };
    
    template<typename LhsDataType, typename RhsDataType, typename RhsStorageType, typename LhsMatrixType, typename RhsMatrixType>
    class BinaryExpressionTraits<NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>, NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>, SubtractOp>
    {
        public:
            typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
            typedef NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> ResultType;
    };
    
    // The following are for rows
    // 11,20, 29, 38, 47, 56, 65, 74, and are for full matrices on the right side of a plus or minus.
    //
    // We have to use boost::disable_if to disable full matrices on the left hand side.  Otherwise there would 
    // be a conflict with the previous specializations.
    template<typename LhsDataType, typename RhsDataType, typename LhsStorageType, typename LhsMatrixType, typename RhsMatrixType>
    class BinaryExpressionTraits<NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>, 
                                 NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>, 
                                 AddOp,
                                 typename boost::disable_if<boost::is_same<LhsStorageType, FullMatrixTag> >::type >
    {
        public:
            typedef typename NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>::NumberType NumberType;
            typedef NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> ResultType;
    };
    
    template<typename LhsDataType, typename RhsDataType, typename LhsStorageType, typename LhsMatrixType, typename RhsMatrixType>
    class BinaryExpressionTraits<NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>, NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>, SubtractOp>
    {
        public:
            typedef typename NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>::NumberType NumberType;
            typedef NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> ResultType;
    };
    
    
    // All multiplications, rows 83 - 163.
    template<typename LhsDataType, typename RhsDataType, typename LhsStorageType, typename RhsStorageType, typename LhsMatrixType, typename RhsMatrixType>
    class BinaryExpressionTraits<NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>, NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>, MultiplyOp>
    {
        public:
            typedef typename NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>::NumberType NumberType;
            typedef NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> ResultType;
    };
}
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_TRAITS_HPP
