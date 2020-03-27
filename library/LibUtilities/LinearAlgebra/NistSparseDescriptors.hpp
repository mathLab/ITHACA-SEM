///////////////////////////////////////////////////////////////////////////////
//
// File: NistSparseDescriptors.hpp
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
// Description: NIST sparse BLAS matrix representation descriptors in the order
//              that matches MatrixStorage enum values (where applicable)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NIST_SPARSE_DESCRIPTORS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NIST_SPARSE_DESCRIPTORS_HPP

namespace Nektar
{

    // Converts Nektar++ dense matrix type enums to
    // NIST Sparse Blas descriptors
    static const int NistSpBlasDescra[5][9] =
    {
        // MatrixStorage::eFULL
        { 0,  // general matrix
          0,  // neither lower- nor upper-triangular
          0,  // non-unit diagonal
          0,  // array base is C++-compatible
          0, 0, 0, 0, 0 // not in use
        },
        // MatrixStorage::eDIAGONAL
        { 5,  // diagonal
          0,  // neither lower- nor upper-triangular
          0,  // non-unit diagonal
          0,  // array base is C++-compatible
          0, 0, 0, 0, 0 // not in use
        },
        // MatrixStorage::eUPPER_TRIANGULAR
        { 3,  // triangular
          2,  // upper-triangular
          0,  // non-unit diagonal
          0,  // array base is C++-compatible
          0, 0, 0, 0, 0 // not in use
        },
        // MatrixStorage::eLOWER_TRIANGULAR
        { 3,  // triangular
          1,  // lower-triangular
          0,  // non-unit diagonal
          0,  // array base is C++-compatible
          0, 0, 0, 0, 0 // not in use
        },
        // MatrixStorage::eSYMMETRIC
        { 1,  // symmetric
          2,  // upper-triangular part to be stored
          //1,  // lower-triangular part to be stored
          0,  // non-unit diagonal
          0,  // array base is C++-compatible
          0, 0, 0, 0, 0 // not in use
        }
    };

} // namespace

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NIST_SPARSE_DESCRIPTORS_HPP
