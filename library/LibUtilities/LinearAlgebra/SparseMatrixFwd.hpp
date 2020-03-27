///////////////////////////////////////////////////////////////////////////////
//
// File: SparseMatrixFwd.hpp
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
// Description: common typedefs and forward declarations for sparse matrices
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_MATRIX_FWD_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_MATRIX_FWD_HPP

#include <map>
#include <vector>
#include <utility>
#include <fstream>
#include <memory>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>


namespace Nektar
{
    typedef unsigned int    IndexType;

    typedef Array<OneD, IndexType>              IndexVector;

    // Elemental COO: each entry is a nonzero number
    typedef std::pair<IndexType, IndexType>     CoordType;
    typedef NekDouble                           COOEntryType;
    typedef std::map<CoordType, NekDouble>      COOMatType;
    typedef COOMatType::const_iterator          COOMatTypeConstIt;
    typedef std::shared_ptr<COOMatType>         COOMatTypeSharedPtr;
    typedef Array<OneD, COOMatType>             COOMatVector;

    // Block COO (BCO): each entry is a dense submatrix (of same size)
    typedef Array<OneD, NekDouble>              BCOEntryType;
    typedef std::map<CoordType, BCOEntryType >  BCOMatType;
    typedef BCOMatType::const_iterator          BCOMatTypeConstIt;
    typedef std::shared_ptr<BCOMatType>         BCOMatTypeSharedPtr;
    typedef Array<OneD, BCOMatType>             BCOMatVector;


    template<typename DataType> class StorageSmvBsr;

    template<typename SparseStorageType> class NekSparseMatrix;
    template<typename SparseStorageType> class NekSparseDiagBlkMatrix;

    

}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_MATRIX_FWD_HPP
