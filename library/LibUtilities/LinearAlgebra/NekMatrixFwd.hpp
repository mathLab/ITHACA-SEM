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

#include <LibUtilities/LinearAlgebra/MatrixBlockType.h>
#include <LibUtilities/LinearAlgebra/PointerWrapper.h>

#include <boost/shared_ptr.hpp>

namespace Nektar
{
  
    /// \brief Enumeration used internally to determine if a pointer can be deleted or not.
    //enum MatrixDataDeltetableType { eDeletable, eNotDeletable };
    
    template<typename DataType, NekMatrixForm form = eFull, MatrixBlockType BlockType = eNormal, unsigned int space = 0, typename enabled=void>
    class NekMatrix;
    
    template<typename DataType, NekMatrixForm form>
    class NekMatrixStoragePolicy;

    template<typename DataType, NekMatrixForm form>
    class NekMatrixArithmeticPolicy;

    template<typename DataType, NekMatrixForm form>
    class NekMatrixAssignmentPolicy;

    typedef boost::shared_ptr<NekMatrix<double> > SharedNekMatrixPtr;
};
    
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_MATRIX_FWD_HPP

/**
    $Log: NekMatrixFwd.hpp,v $
    Revision 1.6  2007/02/15 06:56:55  bnelson
    *** empty log message ***

    Revision 1.5  2006/12/17 22:36:35  bnelson
    Removed Macintosh line endings.

**/


