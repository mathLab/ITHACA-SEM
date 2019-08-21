////////////////////////////////////////////////////////////////////////////////
//
//  File: ShapeType.cpp
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
//  Description: Python wrapper for ShapeType.
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;

template<typename T>
T AllReduce(CommSharedPtr &comm, T toReduce, ReduceOperator oper)
{
    comm->AllReduce(toReduce, oper);
    return toReduce;
}

/**
 * @brief Export for Comm communicator.
 */
void export_Comm()
{
    // Export ReduceOperator enum
    NEKPY_WRAP_ENUM(ReduceOperator, ReduceOperatorMap);

    py::class_<Comm, std::shared_ptr<Comm>,
               boost::noncopyable>("Comm", py::no_init)
        .def("GetSize", &Comm::GetSize)
        .def("GetRank", &Comm::GetRank)
        .def("GetType", &Comm::GetType,
             py::return_value_policy<py::copy_const_reference>())
        .def("AllReduce", &AllReduce<double>)
        .def("AllReduce", &AllReduce<int>)
        .def("AllReduce", &AllReduce<long>)
        .def("AllReduce", &AllReduce<Array<OneD, NekDouble>>)
        ;
}
