///////////////////////////////////////////////////////////////////////////////
//
// File: NekBlockMatrix.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/ConsistentObjectAccess.hpp>

#include <boost/shared_ptr.hpp>

namespace Nektar
{

//     template<MatrixBlockType BlockType, typename MatrixType>
//     class BlockMatrixDataType
//     {
//         public:
//             typedef ConsistentObjectAccess<MatrixType> ResultType;
//     };
// 
//     template<typename MatrixType>
//             class BlockMatrixDataType<ePointerBlock, MatrixType>
//     {
//         public:
//             typedef ConsistentObjectAccess<boost::shared_ptr<MatrixType> > ResultType;
//     };


};
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_BLOCK_MATRIX_HPP

/**
    $Log: NekBlockMatrix.hpp,v $
    Revision 1.5  2007/07/22 23:03:27  bnelson
    Backed out Nektar::ptr.

    Revision 1.4  2007/07/20 00:24:12  bnelson
    Replaced boost::shared_ptr with Nektar::ptr

    Revision 1.3  2007/06/10 23:42:15  bnelson
    Matrix updates.

    Revision 1.2  2006/11/01 04:07:08  bnelson
    Changed block matrices to use the ConsistentObjectAccess object to store matrices or pointers to matrices so that the same pointer syntax works for both.

    Revision 1.1  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

 **/
 
