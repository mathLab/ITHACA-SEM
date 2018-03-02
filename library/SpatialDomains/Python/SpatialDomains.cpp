#include <NekPyConfig.hpp>

void export_Geometry();
void export_Geometry1D();
void export_Geometry2D();
void export_MeshGraph();
void export_SegGeom();
void export_QuadGeom();
void export_TriGeom();

BOOST_PYTHON_MODULE(_SpatialDomains)
{
    np::initialize();

    export_Geometry();
    export_Geometry1D();
    export_Geometry2D();
    export_MeshGraph();
    export_SegGeom();
    export_QuadGeom();
    export_TriGeom();
}
