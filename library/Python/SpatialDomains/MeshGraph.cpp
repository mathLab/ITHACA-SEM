#include <SpatialDomains/MeshGraph.h>
#include <NekPyConfig.hpp>
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
