////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraph.cpp
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
//  Description: Python wrapper for MeshGraph.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/Python/NekPyConfig.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

/*
 * @brief Lightweight wrapper around MeshGraph::Read to avoid wrapping
 * DomainRange struct.
 */
MeshGraphSharedPtr MeshGraph_Read(
    const LibUtilities::SessionReaderSharedPtr &session)
{
    return MeshGraph::Read(session);
}

/**
 * @brief MeshGraph exports.
 */
void export_MeshGraph()
{
    py::class_<SegGeomMap>("SegGeomMap")
        .def(py::map_indexing_suite<SegGeomMap, true>());
    py::class_<QuadGeomMap>("QuadGeomMap")
        .def(py::map_indexing_suite<QuadGeomMap, true>());
    py::class_<TriGeomMap>("TriGeomMap")
        .def(py::map_indexing_suite<TriGeomMap, true>());
    py::class_<MeshGraph,
               std::shared_ptr<MeshGraph>,
               boost::noncopyable>(
                   "MeshGraph", py::no_init)

        .def("Read", MeshGraph_Read)
        .staticmethod("Read")

        .def("GetMeshDimension", &MeshGraph::GetMeshDimension)
        .def("GetAllSegGeoms", &MeshGraph::GetAllSegGeoms,
             py::return_internal_reference<>())
        .def("GetAllQuadGeoms", &MeshGraph::GetAllQuadGeoms,
             py::return_internal_reference<>())
        .def("GetAllTriGeoms",  &MeshGraph::GetAllTriGeoms,
             py::return_internal_reference<>())

        ;
}
