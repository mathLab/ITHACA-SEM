///////////////////////////////////////////////////////////////////////////////
//
// File: Basis.cpp
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
// Description: Python wrapper for Basis.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar::LibUtilities;

BasisSharedPtr Basis_Create(const BasisKey &pts)
{
    return BasisManager()[pts];
}

py::tuple Basis_GetZW(BasisSharedPtr pts)
{
    return py::make_tuple(pts->GetZ(), pts->GetW());
}

/**
 * @brief Basis exports.
 */
void export_Basis()
{
    // Enumerator for basis type
    NEKPY_WRAP_ENUM(BasisType, BasisTypeMap);

    py::class_<BasisKey>("BasisKey", py::init<const BasisType&, const int,
                         const PointsKey&>())

        .def("GetNumModes", &BasisKey::GetNumModes)
        .def("GetTotNumModes", &BasisKey::GetTotNumModes)
        .def("GetNumPoints", &BasisKey::GetNumPoints)
        .def("GetTotNumPoints", &BasisKey::GetTotNumPoints)
        .def("GetBasisType", &BasisKey::GetBasisType)
        .def("GetPointsKey", &BasisKey::GetPointsKey)
        .def("GetPointsType", &BasisKey::GetPointsType)
        .def("Collocation", &BasisKey::Collocation)
        ;

    py::class_<Basis,
               std::shared_ptr<Basis> >(
               "Basis", py::no_init)

        .def("Create", &Basis_Create)
        .staticmethod("Create")

        .def("GetNumModes", &Basis::GetNumModes)
        .def("GetTotNumModes", &Basis::GetTotNumModes)
        .def("GetNumPoints", &Basis::GetNumPoints)
        .def("GetTotNumPoints", &Basis::GetTotNumPoints)
        .def("GetBasisType", &Basis::GetBasisType)
        .def("GetPointsKey", &Basis::GetPointsKey)
        .def("GetBasisKey", &Basis::GetBasisKey)
        .def("GetPointsType", &Basis::GetBasisType)
        .def("Initialize", &Basis::Initialize)

        .def("GetZ", &Basis::GetZ,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetW", &Basis::GetW,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetZW", &Basis_GetZW)

        .def("GetBdata", &Basis::GetBdata,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetDbdata", &Basis::GetDbdata,
             py::return_value_policy<py::copy_const_reference>())
        ;
}
