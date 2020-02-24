///////////////////////////////////////////////////////////////////////////////
//
// File  NekNonlinSys.h
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
// Description: NekNonlinSys header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_H

#include <boost/core/ignore_unused.hpp>

#include <boost/circular_buffer.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <iomanip>
#include <LibUtilities/LinearAlgebra/NekNonlinLinSys.h>

namespace Nektar
{
    namespace LibUtilities
    {
        class NekNonlinSys;

        typedef std::shared_ptr<NekNonlinSys> NekNonlinSysSharedPtr;

        class  NekNonlinSys : public NekNonlinLinSys 
        {
            public:
                LIB_UTILITIES_EXPORT NekNonlinSys(
                    const LibUtilities::SessionReaderSharedPtr  &pSession,
                    const LibUtilities::CommSharedPtr           &vComm,
                    const int                                   nDimen);
                LIB_UTILITIES_EXPORT ~NekNonlinSys();
                
            protected:
            
            private:
        };
    }
}
#endif
