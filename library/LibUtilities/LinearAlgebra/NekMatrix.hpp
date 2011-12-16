///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrix.hpp
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
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/MatrixOperations.hpp>
#include <ExpressionTemplates/CommutativeTraits.hpp>

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
// Only need NekVector if we are using expression templates so we can define the 
// commutative traits.
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#endif
namespace Nektar
{        
    template<typename DataType, typename FormType>
    std::ostream& operator<<(std::ostream& os, const NekMatrix<DataType, FormType>& rhs)
    {
        int oswidth = 9;
        int osprecision = 6;

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            os << "[";
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                os.width(oswidth);
                os.precision(osprecision);
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

    template<typename DataType, typename FormType>
    std::ostream& operator>>(std::ostream& os, const NekMatrix<DataType, FormType>& rhs)
    {
        NekDouble tol = 1e-12;
        os <<  "[" << std::endl;

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                if((NekDouble)rhs(i,j) > tol)
                {
                    os << '+';
                }
                else if((NekDouble)rhs(i,j) < -tol)
                {
                    os << '*';
                }
                else
                {
                    os << '-';
                }
            }
            if( i != rhs.GetRows()-1 )
            {
                os << std::endl;
            }
        }
        os <<  "]" << std::endl;
        return os;
    }
}


    
    #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
namespace expt
{

        template<typename LhsDataType, typename LhsMatrixType,
            typename RhsDataType, typename RhsMatrixType>
        struct CommutativeTraits<NekMatrix<LhsDataType, LhsMatrixType> ,
            expt::MultiplyOp, NekMatrix<RhsDataType, RhsMatrixType> > : public boost::false_type
        {
        };

    template<typename LhsDataType, typename LhsMatrixType,
        typename RhsDataType>
    struct CommutativeTraits<NekMatrix<LhsDataType, LhsMatrixType> ,
        expt::MultiplyOp, NekVector<RhsDataType> > : public boost::false_type
    {
    };
}

    #endif

namespace Nektar
{
//     
//     // Now define general purpose operators.
//     
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     void negate(NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         rhs.Negate();
//     }
//         
//     
// 
// #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
// 
// #else
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     NekMatrix<DataType, form, BlockType, space> Invert(const NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         NekMatrix<DataType, form, BlockType, space> result(rhs);
//         result.Invert();
//         return result;
//     }
// #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
// 
//     /////////////////////////////////////////
//     // Multiplication
//     /////////////////////////////////////////
// //     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
// //     void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
// //                 const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
// //                 typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::MultiplicationResultType& result,
// //                 typename boost::enable_if<boost::is_same<DataType, double> >::type* p = NULL )
// //     {
// // #ifdef NEKTAR_USING_BLAS
// //         dgemm(lhs.GetRows(), lhs.GetColumns(), rhs.GetColumns(), lhs.begin(), rhs.begin(), result.begin());
// // #else
// //         result = lhs;
// //         result *= rhs;
// // #endif
// //     }
// 
// //     template<typename DataType, NekMatrixForm lhsForm, NekMatrixForm rhsForm, MatrixBlockType BlockType, unsigned int space>
// //     void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, 
// //                   const NekMatrix<DataType, rhsForm, BlockType, space>& rhs,
// //                   typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekMatrix<DataType, rhsForm, BlockType, space> >::MultiplicationResultType& result,
// //                 typename boost::disable_if<boost::is_same<DataType, double> >::type* p = NULL )
// //     {
// //         result = lhs;
// //         result *= rhs;
// //     }
// 
// //     template<typename DataType, NekMatrixForm lhsForm, unsigned int vectorDim, MatrixBlockType BlockType, unsigned int space>
// //     void multiply(const NekMatrix<DataType, lhsForm, BlockType, space>& lhs, const NekVector<DataType, vectorDim, space>& rhs,
// //                   typename expt::BinaryExpressionTraits<NekMatrix<DataType, lhsForm, BlockType, space>, NekVector<DataType, vectorDim, space> >::MultiplicationResultType& result)
// //     {
// //         ASSERTL0(lhs.GetColumns() == rhs.GetDimension(), "Invalid matrix dimensions in operator*");
// //         for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
// //         {
// //             DataType t = DataType(0);
// //             for(unsigned int j = 0; j < rhs.GetDimension(); ++j)
// //             {
// //                 t += lhs(i,j)*rhs(j);
// //             }
// //             result(i) = t;
// //         }
// //     }
// 
//                
// 
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     bool operator==(const NekMatrix<DataType, form, BlockType, space>& lhs,
//                     const NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         if( lhs.GetRows() != rhs.GetRows() )
//         {
//             return false;
//         }
// 
//         if( lhs.GetColumns() != rhs.GetColumns() )
//         {
//             return false;
//         }
// 
//         typename NekMatrix<DataType, form, BlockType, space>::const_iterator lhs_iter = lhs.begin();
//         typename NekMatrix<DataType, form, BlockType, space>::const_iterator rhs_iter = rhs.begin();
// 
//         for( ; lhs_iter != lhs.end(); ++lhs_iter, ++rhs_iter )
//         {
//             if( *lhs_iter != *rhs_iter )
//             {
//                 return false;
//             }
//         }
// 
//         return true;
//     }
//     
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     bool operator!=(const NekMatrix<DataType, form, BlockType, space>& lhs,
//                      const NekMatrix<DataType, form, BlockType, space>& rhs)
//     {
//         return !(lhs == rhs);
//     }

}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

