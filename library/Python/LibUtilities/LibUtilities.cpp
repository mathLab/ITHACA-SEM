#include <NekPyConfig.hpp>

void export_Basis();
void export_Points();
void export_SessionReader();
void export_ShapeType();

template<typename T>
void export_SharedArray();

BOOST_PYTHON_MODULE(_LibUtilities)
{
    np::initialize();

    export_Basis();
    export_Points();
    export_SessionReader();
    export_ShapeType();
    export_SharedArray<double>();
}
