///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrixMetadata.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_METADATA_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_METADATA_HPP

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace Nektar
{
   

}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_METADATA_HPP

/**
    $Log: NekMatrixMetadata.hpp,v $
    Revision 1.11  2008/11/11 03:35:07  bnelson
    Updated NekMatrix metadata to capture the largest matrix size during an expression.

    Revision 1.10  2008/01/20 03:59:36  bnelson
    Expression template updates.

    Revision 1.9  2007/10/03 03:00:14  bnelson
    Added precompiled headers.

    Revision 1.8  2007/08/16 02:11:57  bnelson
    *** empty log message ***

    Revision 1.7  2007/06/10 23:42:16  bnelson
    Matrix updates.

    Revision 1.6  2007/01/16 05:30:34  bnelson
    Major improvements for expression templates.

    Revision 1.5  2006/11/08 04:16:14  bnelson
    Added subtraction operators.

    Revision 1.4  2006/10/30 05:11:16  bnelson
    Added preliminary linear system and block matrix support.

    Revision 1.3  2006/09/30 15:18:37  bnelson
    no message

    Revision 1.2  2006/09/14 02:06:16  bnelson
    Fixed gcc compiler errors.

    Revision 1.1  2006/09/11 03:26:27  bnelson
    Updated to use new policy based expression templates.

 **/


