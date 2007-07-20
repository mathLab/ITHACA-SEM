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

#include <LibUtilities/BasicUtils/SharedPtr.hpp>

namespace Nektar
{
    template<typename DataType>
    class ConstMatrix;
    
    template<typename DataType>
    class Matrix;
    
    template<typename DataType, typename StorageType = FullMatrixTag, typename MatType = StandardMatrixTag>
    class NekMatrix;

    template<typename DataType, typename StorageType, typename OwnedMatrixType>
    class NekMatrix<NekMatrix<DataType, StorageType, OwnedMatrixType>, StorageType, ScaledMatrixTag>;
    
    template<typename DataType, typename StorageType, typename InnerMatrixType>
    class NekMatrix<NekMatrix<DataType, StorageType, InnerMatrixType>, StorageType, BlockMatrixTag>;
    
    template<typename DataType, typename StorageType>
    class NekMatrix<DataType, StorageType, StandardMatrixTag>;
    
    typedef ptr<NekMatrix<NekDouble> > SharedNekMatrixPtr;
    typedef NekMatrix<NekMatrix<NekDouble>, FullMatrixTag, ScaledMatrixTag> DNekScalMat;
    typedef ptr<DNekScalMat> DNekScalMatSharedPtr;
};
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_MATRIX_FWD_HPP

/**
    $Log: NekMatrixFwd.hpp,v $
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


