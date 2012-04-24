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
// Description: TypeDefs for Nek Matrices and vector.
// 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_TYPEDEFS_HPP
#define NEKTAR_LIB_UTILITIES_NEK_TYPEDEFS_HPP

#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <LibUtilities/LinearAlgebra/SparseStandardMatrix.hpp>
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
    typedef NekSparseMatrix<NekDouble>
            DNekSparseMat;

    typedef boost::shared_ptr<DNekMat>          DNekMatSharedPtr;
    typedef boost::shared_ptr<DNekScalMat>      DNekScalMatSharedPtr;
    typedef boost::shared_ptr<DNekBlkMat>    DNekBlkMatSharedPtr;
    typedef boost::shared_ptr<BlkMatDNekBlkMat>  BlkMatDNekBlkMatSharedPtr;
    typedef boost::shared_ptr<DNekScalBlkMat>      DNekScalBlkMatSharedPtr;
    typedef boost::shared_ptr<BlkMatDNekScalBlkMat>  BlkMatDNekScalBlkMatSharedPtr;
    typedef boost::shared_ptr<BlkMatDNekScalBlkMat>  BlkMatDNekScalBlkMatSharedPtr;
    typedef boost::shared_ptr<DNekSparseMat>     DNekSparseMatSharedPtr;


    static DNekMatSharedPtr NullDNekMatSharedPtr;
    static DNekScalMatSharedPtr NullDNekScalMatSharedPtr;
    static DNekScalBlkMatSharedPtr NullDNekScalBlkMatSharedPtr;

    typedef LinearSystem                        DNekLinSys;
    typedef boost::shared_ptr<DNekLinSys>       DNekLinSysSharedPtr;

    typedef LinearSystem                        DNekScalLinSys; 
    typedef boost::shared_ptr<DNekScalLinSys>   DNekScalLinSysSharedPtr;

}

#endif //NEKTAR_LIB_UTILITIES_NEK_TYPEDEFS_HPP

/**
    $Log: NekTypeDefs.hpp,v $
    Revision 1.11  2008/11/01 19:15:28  bnelson
    Updated matrices so the storage policy is no longer a template parameter.  Removed the template parameter from the LinearSystem class.

    Revision 1.10  2008/03/12 15:22:45  pvos
    Clean up of the code

    Revision 1.9  2007/10/03 11:37:50  sherwin
    Updates relating to static condensation implementation

    Revision 1.8  2007/07/26 08:40:49  sherwin
    Update to use generalised i/o hooks in Helmholtz1D

    Revision 1.7  2007/07/22 23:03:28  bnelson
    Backed out Nektar::ptr.

    Revision 1.6  2007/07/20 00:24:13  bnelson
    Replaced boost::shared_ptr with Nektar::ptr

    Revision 1.5  2007/07/13 09:02:21  sherwin
    Mods for Helmholtz solver

    Revision 1.4  2007/06/10 23:42:16  bnelson
    Matrix updates.

    Revision 1.3  2007/02/24 09:08:41  sherwin
    Updated to include definition of NekLinSysSharedPtr

    Revision 1.2  2007/02/17 01:11:27  bnelson
    *** empty log message ***

    Revision 1.1  2007/01/23 23:21:26  sherwin
    Modes so that we can use Lapack: definitions and added NekTypeDefs.hpp

 **/


