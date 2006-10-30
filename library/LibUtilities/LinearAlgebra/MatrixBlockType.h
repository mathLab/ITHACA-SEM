///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixImplType.h
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

#ifndef NEKTAR_LIB_UTILIITIES_LINEAR_ALBEGRA_MATRIX_BLOCK_TYPE_H
#define NEKTAR_LIB_UTILIITIES_LINEAR_ALBEGRA_MATRIX_BLOCK_TYPE_H

#include <LibUtilities/LinearAlgebra/NekMatrixForm.h>

#include <boost/shared_ptr.hpp>

namespace Nektar
{
            
    enum MatrixBlockType 
    { 
        eNormal, 
        eBlock,
        ePointerBlock
    };
    
//     template<typename DataType, NekMatrixForm form, MatrixBlockType BlockType, unsigned int space>
//     class NekMatrix;
//             
//     template<typename DataType, NekMatrixForm form, MatrixBlockType type, unsigned int space>
//     class MatrixDataType
//     {
//         public:
//             typedef DataType Type;
//     };
// 
//     template<typename DataType, NekMatrixForm form, unsigned int space>
//     class MatrixDataType<DataType, form, eBlock, space>
//     {
//         public:
//             typedef NekMatrix<DataType, form, eNormal, space> Type;
//     };
// 
//     template<typename DataType, NekMatrixForm form, unsigned int space>
//     class MatrixDataType<DataType, form, ePointerBlock, space>
//     {
//         public:
//             typedef boost::shared_ptr<NekMatrix<DataType, form, eNormal, space> > Type;
//     };
}

#endif //NEKTAR_LIB_UTILIITIES_LINEAR_ALBEGRA_MATRIX_BLOCK_TYPE_H
