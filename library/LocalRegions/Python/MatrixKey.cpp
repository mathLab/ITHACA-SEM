////////////////////////////////////////////////////////////////////////////////
//
//  File: MatrixKey.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Python wrapper for MatrixKey.
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <LocalRegions/MatrixKey.h>

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

using namespace Nektar;
using namespace Nektar::LocalRegions;

MatrixKey *MatrixKey_Init(const StdRegions::MatrixType matType,
                          const LibUtilities::ShapeType shapeType,
                          const StdRegions::StdExpansionSharedPtr exp,
                          const py::object &constFactorMap,
                          const py::object &varCoeffMap)
{
    StdRegions::ConstFactorMap tmp = StdRegions::NullConstFactorMap;
    StdRegions::VarCoeffMap tmp2 = StdRegions::NullVarCoeffMap;

    if (!constFactorMap.is_none())
    {
        tmp = py::extract<StdRegions::ConstFactorMap>(constFactorMap);
    }

    if (!varCoeffMap.is_none())
    {
        tmp2 = py::extract<StdRegions::VarCoeffMap>(varCoeffMap);
    }

    return new MatrixKey(matType, shapeType, *exp, tmp, tmp2);
}

/**
 * @brief Export for MatrixKey enumeration.
 */
void export_MatrixKey()
{
    py::class_<MatrixKey, py::bases<StdRegions::StdMatrixKey>>(
        "MatrixKey", py::no_init)
        .def("__init__",
             py::make_constructor(
                 &MatrixKey_Init, py::default_call_policies(),
                 (py::arg("matType"), py::arg("shapeType"), py::arg("exp"),
                  py::arg("constFactorMap") = py::object(),
                  py::arg("varCoeffMap") = py::object())))
        ;
}
