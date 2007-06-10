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
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
namespace Nektar
{
    typedef NekMatrix<NekDouble,FullMatrixTag>        DNekMat;
    typedef NekVector<NekDouble>              DNekVec;
    typedef LinearSystem <DNekMat>         DNekLinSys;

    typedef boost::shared_ptr<DNekMat>     DNekMatSharedPtr;
    typedef boost::shared_ptr<DNekLinSys>  DNekLinSysSharedPtr;
}

#endif //NEKTAR_LIB_UTILITIES_NEK_TYPEDEFS_HPP

/**
    $Log: NekTypeDefs.hpp,v $
    Revision 1.3  2007/02/24 09:08:41  sherwin
    Updated to include definition of NekLinSysSharedPtr

    Revision 1.2  2007/02/17 01:11:27  bnelson
    *** empty log message ***

    Revision 1.1  2007/01/23 23:21:26  sherwin
    Modes so that we can use Lapack: definitions and added NekTypeDefs.hpp

 **/


