///////////////////////////////////////////////////////////////////////////////
//
// File: testLinearSystem.h
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
// Description: Test code for NekVector
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_UNIT_TESTS_TEST_LINEAR_SYSTEM_H
#define NEKTAR_UNIT_TESTS_TEST_LINEAR_SYSTEM_H

namespace Nektar
{
    namespace LinearSystemUnitTests
    {
        void testDiagonalSystem();
        void testFullSystem();
        void testSolvingBlockDiagonalMatrices();
        void testMixedInputParameterTypes();
    }
}

#endif //NEKTAR_UNIT_TESTS_TEST_LINEAR_SYSTEM_H

/**
    $Log: testLinearSystem.h,v $
    Revision 1.5  2007/07/10 01:22:39  bnelson
    Unit tests for full linear system with both x and b as input parameters.

    Revision 1.4  2007/02/13 02:47:23  bnelson
    *** empty log message ***

    Revision 1.3  2006/10/30 05:08:13  bnelson
    Added preliminary linear system and block matrix support.

    Revision 1.2  2006/10/02 01:20:38  bnelson
    Started working on adding BLAS and LAPACK

    Revision 1.1  2006/09/30 15:38:29  bnelson
    no message

**/
