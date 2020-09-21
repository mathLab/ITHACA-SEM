///////////////////////////////////////////////////////////////////////////////
//
// File: NekTypeDefs.hpp
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
// Description: TypeDefs for Nek Matrices and vector.
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_TYPEDEFS_HPP
#define NEKTAR_LIB_UTILITIES_NEK_TYPEDEFS_HPP

#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>

#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

namespace Nektar
{
    typedef NekVector<NekDouble>                DNekVec;

    typedef NekMatrix<NekDouble, StandardMatrixTag>
            DNekMat;
    typedef NekMatrix<NekDouble, StandardMatrixTag> DenseMatrix;

    typedef NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>
            DNekScalMat;
    typedef NekMatrix<DenseMatrix, ScaledMatrixTag> ScaledMatrix;

    typedef NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, BlockMatrixTag>
            DNekBlkMat;
    typedef NekMatrix<DenseMatrix, BlockMatrixTag> BlockMatrix;

    typedef NekMatrix<NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, BlockMatrixTag>, BlockMatrixTag>
            BlkMatDNekBlkMat;
    typedef NekMatrix<NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>, BlockMatrixTag>
            DNekScalBlkMat;
    typedef NekMatrix<NekMatrix<NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>, BlockMatrixTag>, BlockMatrixTag>
            BlkMatDNekScalBlkMat;

    typedef std::shared_ptr<DNekMat>               DNekMatSharedPtr;
    typedef std::shared_ptr<DNekScalMat>           DNekScalMatSharedPtr;
    typedef std::shared_ptr<DNekBlkMat>            DNekBlkMatSharedPtr;
    typedef std::shared_ptr<BlkMatDNekBlkMat>      BlkMatDNekBlkMatSharedPtr;
    typedef std::shared_ptr<DNekScalBlkMat>        DNekScalBlkMatSharedPtr;
    typedef std::shared_ptr<BlkMatDNekScalBlkMat>  BlkMatDNekScalBlkMatSharedPtr;
    typedef std::shared_ptr<BlkMatDNekScalBlkMat>  BlkMatDNekScalBlkMatSharedPtr;


    static DNekMatSharedPtr NullDNekMatSharedPtr;
    static DNekScalMatSharedPtr NullDNekScalMatSharedPtr;
    static DNekScalBlkMatSharedPtr NullDNekScalBlkMatSharedPtr;

    typedef LinearSystem                      DNekLinSys;
    typedef std::shared_ptr<DNekLinSys>       DNekLinSysSharedPtr;

    typedef LinearSystem                      DNekScalLinSys;
    typedef std::shared_ptr<DNekScalLinSys>   DNekScalLinSysSharedPtr;

    static Array<OneD, DNekBlkMatSharedPtr>  NullArrayDNekBlkMatSharedPtr;

}

#endif //NEKTAR_LIB_UTILITIES_NEK_TYPEDEFS_HPP
