////////////////////////////////////////////////////////////////////////////////
//
//  File: StdMatrixKey.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description: Python wrapper for StdMatrixKey.
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdMatrixKey.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

#include <boost/python/suite/indexing/map_indexing_suite.hpp>

using namespace Nektar;
using namespace Nektar::StdRegions;

StdMatrixKey *StdMatrixKey_Init(const MatrixType matType,
                                const LibUtilities::ShapeType shapeType,
                                const StdExpansionSharedPtr exp,
                                const py::object &constFactorMap)
{
    ConstFactorMap tmp = NullConstFactorMap;

    if (!constFactorMap.is_none())
    {
        tmp = py::extract<ConstFactorMap>(constFactorMap);
    }

    return new StdMatrixKey(matType, shapeType, *exp, tmp);
}

/**
 * @brief Export for StdMatrixKey enumeration.
 */
void export_StdMatrixKey()
{
    NEKPY_WRAP_ENUM(MatrixType, MatrixTypeMap);
    NEKPY_WRAP_ENUM(ConstFactorType, ConstFactorTypeMap);

    // Wrapper for constant factor map.
    py::class_<ConstFactorMap>("ConstFactorMap")
        .def(py::map_indexing_suite<ConstFactorMap, true>());

    py::class_<StdMatrixKey>("StdMatrixKey", py::no_init)
        .def("__init__", py::make_constructor(
                 &StdMatrixKey_Init, py::default_call_policies(),
                 (py::arg("matType"), py::arg("shapeType"), py::arg("exp"),
                  py::arg("constFactorMap") = py::object())))

        .def("GetMatrixType", &StdMatrixKey::GetMatrixType)
        .def("GetShapeType",  &StdMatrixKey::GetShapeType)
        .def("GetNcoeffs",    &StdMatrixKey::GetNcoeffs)
        .def("GetBasis",      &StdMatrixKey::GetBasis)
        ;
}
