///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrixFwd.hpp
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
// Description: Matrix Forward Declarations
//
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_MATRIX_FWD_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_MATRIX_FWD_HPP

#include <LibUtilities/LinearAlgebra/MatrixType.h>
#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/type_traits.hpp>

#include <boost/typeof/typeof.hpp>
#include  BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

namespace Nektar
{
    template<typename DataType>
    class ConstMatrix;
    
    template<typename DataType>
    class Matrix;
    
    template<typename DataType, typename MatType = StandardMatrixTag>
    class NekMatrix;

    template<typename DataType, typename InnerMatrixType>
    class NekMatrix<NekMatrix<DataType, InnerMatrixType>, ScaledMatrixTag>;
    
    template<typename DataType, typename InnerMatrixType>
    class NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>;
    
    template<typename DataType>
    class NekMatrix<DataType, StandardMatrixTag>;
    
    typedef boost::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > SharedNekMatrixPtr;
    typedef NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag> DNekScalMat;
    typedef boost::shared_ptr<DNekScalMat> DNekScalMatSharedPtr;
    
    // Type registration must occur for the expression template machinery to 
    // automatically detect the types of matrix operations.
    BOOST_TYPEOF_REGISTER_TEMPLATE(NekMatrix, 2);
    
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
    template<typename T>
    struct IsMatrix : public boost::false_type {};
    
    template<typename DataType, typename MatrixType>
    struct IsMatrix<NekMatrix<DataType, MatrixType> > : public boost::true_type {};
    
    template<typename DataType>
    struct IsMatrix<ConstMatrix<DataType> > : public boost::true_type {};
    
    template<typename DataType>
    struct IsMatrix<Matrix<DataType> > : public boost::true_type {};
#endif
    
};
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_MATRIX_FWD_HPP

/**
    $Log: NekMatrixFwd.hpp,v $
    Revision 1.15  2008/03/09 22:39:29  bnelson
    Added IsMatrix.

    Revision 1.14  2007/08/16 02:08:12  bnelson
    Removed typeof registration

    Revision 1.13  2007/07/25 23:47:45  bnelson
    *** empty log message ***

    Revision 1.12  2007/07/22 23:03:28  bnelson
    Backed out Nektar::ptr.

    Revision 1.11  2007/07/20 00:24:13  bnelson
    Replaced boost::shared_ptr with Nektar::ptr

    Revision 1.10  2007/07/12 04:04:14  bnelson
    *** empty log message ***

    Revision 1.9  2007/06/24 17:59:33  bnelson
    *** empty log message ***

    Revision 1.8  2007/06/10 23:42:15  bnelson
    Matrix updates.

    Revision 1.7  2007/03/29 18:59:05  bnelson
    Refactoring in preparation for scaled matrices.  Fixed transpose problem.

    Revision 1.6  2007/02/15 06:56:55  bnelson
    *** empty log message ***

    Revision 1.5  2006/12/17 22:36:35  bnelson
    Removed Macintosh line endings.

**/


