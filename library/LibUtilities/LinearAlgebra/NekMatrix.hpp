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

namespace Nektar
{        
    template<typename DataType, typename FormType>
    std::ostream& operator<<(std::ostream& os, const NekMatrix<DataType, FormType>& rhs)
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

    
    #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
        //template<typename T, typename R>
        //struct CommutativeTraits<T, MultiplyOp, R,
        //    typename boost::enable_if
        //    <
        //        boost::mpl::and_
        //        <
        //            IsMatrix<T>,
        //            IsMatrix<R>
        //        >
        //    >::type
        //>: public boost::false_type
        //{
        //};

        template<typename LhsDataType, typename LhsMatrixType,
            typename RhsDataType, typename RhsMatrixType>
        struct CommutativeTraitsSpecialization<NekMatrix<LhsDataType, LhsMatrixType> ,
            MultiplyOp, NekMatrix<RhsDataType, RhsMatrixType> > : public boost::false_type
        {
        };
    #endif
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

/**
    $Log: NekMatrix.hpp,v $
    Revision 1.37  2008/11/13 01:56:34  bnelson
    Started BLAS specializations using expression templates.

    Revision 1.36  2008/11/01 22:36:06  bnelson
    Removed uneeded files.

    Revision 1.35  2008/11/01 19:15:28  bnelson
    Updated matrices so the storage policy is no longer a template parameter.  Removed the template parameter from the LinearSystem class.

    Revision 1.34  2008/04/01 15:56:38  bnelson
    Added global transpose for scaled matrices.

    Revision 1.33  2008/01/20 20:09:20  bnelson
    *** empty log message ***

    Revision 1.32  2008/01/20 03:59:36  bnelson
    Expression template updates.

    Revision 1.31  2007/10/03 03:00:14  bnelson
    Added precompiled headers.

    Revision 1.30  2007/09/12 04:01:18  bnelson
    Updates for column major matrices.

    Revision 1.29  2007/08/27 04:41:29  bnelson
    Updated the existing storage policies in preparation for banded matrices.

    Revision 1.28  2007/08/16 02:12:33  bnelson
    *** empty log message ***

    Revision 1.27  2007/06/17 22:57:01  bnelson
    Added MatrixOperation files.

    Revision 1.26  2007/06/10 23:42:15  bnelson
    Matrix updates.

    Revision 1.25  2007/04/05 05:12:45  bnelson
    *** empty log message ***

    Revision 1.24  2007/04/04 02:11:08  bnelson
    Added inversion

    Revision 1.23  2007/04/03 04:00:27  bnelson
    Started adding inversion operations to NekMatrix.

    Revision 1.22  2007/03/29 18:59:05  bnelson
    Refactoring in preparation for scaled matrices.  Fixed transpose problem.

    Revision 1.21  2007/01/29 01:31:07  bnelson
    *** empty log message ***

    Revision 1.20  2007/01/23 03:12:50  jfrazier
    Added more conditional compilation directives for expression templates.

    Revision 1.19  2007/01/17 01:11:21  bnelson
    Removed old code.

    Revision 1.18  2007/01/16 05:30:33  bnelson
    Major improvements for expression templates.

    Revision 1.17  2006/11/08 04:16:14  bnelson
    Added subtraction operators.

    Revision 1.16  2006/11/06 17:09:10  bnelson
    *** empty log message ***

    Revision 1.15  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

    Revision 1.14  2006/10/04 03:07:59  bnelson
    Added a check around blas.h to not break existing code.

    Revision 1.13  2006/10/04 03:02:36  bnelson
    Fixed a conflict problem from the previous commit.

    Revision 1.12  2006/10/02 01:16:14  bnelson
    Started working on adding BLAS and LAPACK

    Revision 1.11  2006/09/30 15:18:37  bnelson
    no message

    Revision 1.10  2006/09/21 01:03:31  bnelson
    Added addition and subtraction expression templates.

    Revision 1.9  2006/09/16 23:53:35  bnelson
    Modified the negation operation to reflect changes in the unary expression templates.

    Revision 1.8  2006/09/14 02:06:16  bnelson
    Fixed gcc compiler errors.

    Revision 1.7  2006/09/11 03:26:26  bnelson
    Updated to use new policy based expression templates.

    Revision 1.6  2006/08/28 02:40:21  bnelson
    *** empty log message ***

    Revision 1.5  2006/08/25 03:05:16  bnelson
    Fixed gcc compile errors.

    Revision 1.4  2006/08/25 01:28:53  bnelson
    Changed the way specialized matrices are handled.

    Added support for a matrix wrapper around a pre-allocated array of data.

    Added the NekMemoryManager for allocations.

    Revision 1.3  2006/08/14 02:29:49  bnelson
    Updated points, vectors, and matrix classes to work with ElVis.  Added a variety of methods to all of these classes.

    Revision 1.2  2006/06/01 13:44:28  kirby
    *** empty log message ***

    Revision 1.1  2006/06/01 09:12:41  kirby
    *** empty log message ***

    Revision 1.11  2006/05/31 23:24:21  bnelson
    Moved the diagonal matrix to NekMatrix.hpp.

    Updated method names for the coding standard.

    Revision 1.10  2006/05/31 04:20:16  bnelson
    Changed matrix implementation so the form is a template parameter.

    Revision 1.9  2006/05/29 04:32:18  bnelson
    Changed the data holder to boost::shared_array.

    Revision 1.8  2006/05/29 03:45:04  bnelson
    Updated operator+= to be more efficient using iterators.

    Revision 1.7  2006/05/29 03:40:12  bnelson
    Changed the implementation from individual NekMatrixImpl objects for each type of matrix to an enumeration to store the type.

    Revision 1.6  2006/05/25 03:02:40  bnelson
    Added Matrix/Vector multiplication.

    Revision 1.5  2006/05/18 04:21:06  bnelson
    Removed the ability to specify the rows and columns in the template parameter list.  If this behavior is desired we'll need to create a fixed size array class.

    Added multiplication to the arrays.

    Revision 1.4  2006/05/15 05:06:55  bnelson
    Removed use of MemoryManager pending review of some memory problems.

    Revision 1.3  2006/05/15 04:13:36  bnelson
    no message

    Revision 1.2  2006/05/14 21:32:03  bnelson
 *** empty log message ***

    Revision 1.1  2006/05/04 18:57:43  kirby
 *** empty log message ***

    Revision 1.1  2006/04/11 02:00:43  bnelson
    Initial Revision


 **/


