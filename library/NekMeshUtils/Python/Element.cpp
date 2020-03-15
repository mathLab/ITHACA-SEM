///////////////////////////////////////////////////////////////////////////////
//
// File: Element.cpp
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
// Description: Python wrapper for Element.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <NekMeshUtils/MeshElements/Element.h>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace Nektar;
using namespace Nektar::NekMeshUtils;

ElmtConfig *ElmtConfig_Init(
    LibUtilities::ShapeType shapeType,
    unsigned int            order,
    bool                    faceNodes,
    bool                    volNodes,
    const py::object       &reorient,
    const py::object       &edgeNodeType,
    const py::object       &faceNodeType)
{
    // These are the three default arguments that are present in the ElmtConfig
    // constructor.
    bool orient = true;
    LibUtilities::PointsType eType = LibUtilities::ePolyEvenlySpaced;
    LibUtilities::PointsType fType = LibUtilities::ePolyEvenlySpaced;

    // See whether we were supplied any of the optional keyword arguments from
    // the Python side. If so then override the default arguments.
    if (!reorient.is_none())
    {
        orient = py::extract<bool>(reorient);
    }

    if (!edgeNodeType.is_none())
    {
        eType = py::extract<LibUtilities::PointsType>(edgeNodeType);
    }

    if (!faceNodeType.is_none())
    {
        fType = py::extract<LibUtilities::PointsType>(faceNodeType);
    }

    // Create a new ElmtConfig object. We return a pointer deliberately, since
    // boost::python has to store either a pointer or a shared_ptr to this
    // object.
    return new ElmtConfig(
        shapeType, order, faceNodes, volNodes, orient, eType, fType);
}

std::shared_ptr<Element> Element_Create(ElmtConfig *conf,
                                        py::list &nodes,
                                        py::list &tags)
{
    int numNodes = py::len(nodes), numTags = py::len(tags);

    std::vector<NodeSharedPtr> nodes_(numNodes);
    for (int i = 0; i < numNodes; ++i)
    {
        nodes_[i] = py::extract<NodeSharedPtr>(nodes[i]);
    }

    std::vector<int> tags_(numTags);
    for (int i = 0; i < numTags; ++i)
    {
        tags_[i] = py::extract<int>(tags[i]);
    }

    return GetElementFactory().CreateInstance(conf->m_e, *conf, nodes_, tags_);
}

struct ElementMapHelper
{
    using T = ElementMap;
    using V = std::vector<ElementSharedPtr>;

    static V& get(T &x, int i)
    {
        if (i < 0 || i > 3)
        {
            IndexError();
        }
        return x[i];
    }
    static int len(T const &x)
    {
        return 4;
    }
    static void IndexError()
    {
        PyErr_SetString(PyExc_IndexError, "Index out of range");
    }
};

void export_Element()
{
    // Export element configuration struct
    py::class_<ElmtConfig>("ElmtConfig", py::no_init)
        .def("__init__", py::make_constructor(
                 &ElmtConfig_Init, py::default_call_policies(),
                 (py::arg("shapeType"), py::arg("order"), py::arg("faceNodes"),
                  py::arg("volNodes"), py::arg("reorient") = py::object(),
                  py::arg("edgeNodeType") = py::object(),
                  py::arg("faceNodeType") = py::object())));

    // Export element base class
    py::class_<Element,
               std::shared_ptr<Element>,
               boost::noncopyable>(
                   "Element", py::no_init)
        .def("GetId", &Element::GetId)
        .def("GetDim", &Element::GetDim)
        .def("GetShapeType", &Element::GetShapeType)
        .def("GetTag", &Element::GetTag)

        // Factory methods
        .def("Create", &Element_Create)
        .staticmethod("Create")
        ;

    // Export handler for element mesh storage (ElementMap)
    py::class_<ElementMap>("ElementMap")
        .def("__len__", &ElementMapHelper::len)
        .def("__getitem__", &ElementMapHelper::get,
             py::return_internal_reference<>())
        ;

    // Export handler for list of elements (std::vector<ElementSharedPtr>)
    py::class_<std::vector<ElementSharedPtr>>("ElementVector")
        .def(py::vector_indexing_suite<std::vector<ElementSharedPtr>>())
    ;
}
