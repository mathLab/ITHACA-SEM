///////////////////////////////////////////////////////////////////////////////
//
// File: testNekVector.h
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


#ifndef NEKTAR_UNIT_TESTS_TEST_NEK_VECTOR_H
#define NEKTAR_UNIT_TESTS_TEST_NEK_VECTOR_H

namespace Nektar
{
    namespace UnitTests
    {
        void testNekVectorConstruction();
        void testNekVectorOperators();
        void testNekVectorArithmetic();
        void testNorms();
        void TestMatrixVectorMultiply();
        void TestVectorConstructorsWithSizeArguments();
    }
}

#endif // NEKTAR_UNIT_TESTS_TEST_NEK_POINT_H

/**
    $Log: testNekVector.h,v $
    Revision 1.2  2006/11/18 17:18:46  bnelson
    Added L1, L2, and Infinity norm tests.

    Revision 1.1  2006/05/04 18:59:56  kirby
    *** empty log message ***

    Revision 1.1  2006/04/11 02:02:13  bnelson
    Added more tests.

**/

