///////////////////////////////////////////////////////////////////////////////
//
// File: Node.cpp
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
// Description: Python wrapper for Node.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <NekMeshUtils/MeshElements/Node.h>

using namespace Nektar;
using namespace Nektar::NekMeshUtils;

void KeyError()
{
    PyErr_SetString(PyExc_KeyError, "Key not found");
}

template<class T>
struct unordered_set_item
{
    using K = typename T::key_type;
    static bool contains(T const &x, K const &k)
    {
        return x.find(k) != x.end();
    }
    static void add(T &x, K const &k)
    {
        x.insert(k);
    }
};

void export_Node()
{
    py::class_<Node,
               std::shared_ptr<Node>,
               boost::noncopyable>(
                   "Node", py::init<int, NekDouble, NekDouble, NekDouble>())
        .def("GetID", &Node::GetID)
        .def("SetID", &Node::SetID)
        .def("Distance", &Node::Distance)
        .def("GetLoc", &Node::GetLoc)
        .def("abs2", &Node::abs2)
        .def_readwrite("x", &Node::m_x)
        .def_readwrite("y", &Node::m_y)
        .def_readwrite("z", &Node::m_z)
        .def_readwrite("id", &Node::m_id)
        ;

    // Create converter for NodeSet
    py::class_<NodeSet>("NodeSet")
        .def("__len__", &NodeSet::size)
        .def("clear", &NodeSet::clear)
        .def("__iter__", py::iterator<NodeSet>())
        .def("__contains__", &unordered_set_item<NodeSet>::contains)
        .def("add", &unordered_set_item<NodeSet>::add,
             py::with_custodian_and_ward<1,2>())
        ;
}
