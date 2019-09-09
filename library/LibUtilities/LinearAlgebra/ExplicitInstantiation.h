///////////////////////////////////////////////////////////////////////////////
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_LIBUTILITIES_LINEARALGEBRA_EXPLICIT_INSTANTIATION_H
#define NEKTAR_LIBUTILITIES_LINEARALGEBRA_EXPLICIT_INSTANTIATION_H

#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/array/elem.hpp>
#include <boost/preprocessor/array/size.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/punctuation/paren.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/logical/bool.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/comparison/greater.hpp>
#include <boost/preprocessor/array/pop_front.hpp>
#include <boost/preprocessor/array/push_back.hpp>

#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

#define BOOST_PP_TUPLE_REM_0()

// Macros to make creating explicit instantiations of all possible matrix types for methods easier.
#define NEKTAR_ALL_MATRIX_TYPES (6, (const DNekMat&, const DNekScalMat&, const DNekBlkMat&, const BlkMatDNekBlkMat&, const DNekScalBlkMat&, const BlkMatDNekScalBlkMat&))
#define NEKTAR_BLOCK_MATRIX_TYPES (4, (const DNekBlkMat&, const BlkMatDNekBlkMat&, const DNekScalBlkMat&, const BlkMatDNekScalBlkMat&))
#define NEKTAR_STANDARD_AND_SCALED_MATRICES (2, (const DNekMat&, const DNekScalMat&))

#define NEKTAR_PRINT_ARRAY(z, n, data) \
    BOOST_PP_ARRAY_ELEM(n, data) \
    BOOST_PP_COMMA_IF(BOOST_PP_LESS(n, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE(data), 1)))

#define NEKTAR_CREATE_EXPLICIT_INTSTANTIATION(z, n, data) \
    template LIB_UTILITIES_EXPORT \
    BOOST_PP_ARRAY_ELEM(0, BOOST_PP_ARRAY_ELEM(2, data)) BOOST_PP_ARRAY_ELEM(0, data) BOOST_PP_LPAREN()  \
    BOOST_PP_REPEAT(BOOST_PP_ARRAY_SIZE(BOOST_PP_ARRAY_ELEM(3, data)), NEKTAR_PRINT_ARRAY, BOOST_PP_ARRAY_ELEM(3, data)) \
    BOOST_PP_COMMA_IF(BOOST_PP_GREATER(BOOST_PP_ARRAY_SIZE(BOOST_PP_ARRAY_ELEM(3, data)), 0)) \
    BOOST_PP_ARRAY_ELEM(n, BOOST_PP_ARRAY_ELEM(1, data))\
    BOOST_PP_COMMA_IF(BOOST_PP_GREATER(BOOST_PP_ARRAY_SIZE(BOOST_PP_ARRAY_ELEM(4, data)), 0)) \
    BOOST_PP_REPEAT(BOOST_PP_ARRAY_SIZE(BOOST_PP_ARRAY_ELEM(4, data)), NEKTAR_PRINT_ARRAY, BOOST_PP_ARRAY_ELEM(4, data)) \
    BOOST_PP_RPAREN() ;

#define NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(MethodName, MatrixTypes, ReturnType, BeforeArgs, AfterArgs) \
    BOOST_PP_REPEAT(BOOST_PP_ARRAY_SIZE(MatrixTypes), NEKTAR_CREATE_EXPLICIT_INTSTANTIATION, (5, (MethodName, MatrixTypes, ReturnType, BeforeArgs, AfterArgs)))


#define NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES_INNER(z, n, data) \
    BOOST_PP_REPEAT(BOOST_PP_ARRAY_SIZE(BOOST_PP_ARRAY_ELEM(2, data)),  NEKTAR_CREATE_EXPLICIT_INTSTANTIATION, (5, (BOOST_PP_ARRAY_ELEM(0, data), BOOST_PP_ARRAY_ELEM(2, data), BOOST_PP_ARRAY_ELEM(3, data), BOOST_PP_ARRAY_PUSH_BACK(BOOST_PP_ARRAY_ELEM(4, data), BOOST_PP_ARRAY_ELEM(n, BOOST_PP_ARRAY_ELEM(1, data))),  BOOST_PP_ARRAY_ELEM(5, data))))


// Assumes the matrices are adjacent in parameter list.
#define NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES(MethodName, FirstMatrixTypes, SecondMatrixTypes, ReturnType, BeforeArgs, AfterArgs) \
    BOOST_PP_REPEAT(BOOST_PP_ARRAY_SIZE(FirstMatrixTypes), NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_TWO_MATRICES_INNER, (6, (MethodName, FirstMatrixTypes, SecondMatrixTypes, ReturnType, BeforeArgs, AfterArgs)))


#endif //NEKTAR_LIBUTILITIES_LINEARALGEBRA_EXPLICIT_INSTANTIATION_H
