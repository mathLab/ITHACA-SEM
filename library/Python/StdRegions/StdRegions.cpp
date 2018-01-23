#include <NekPyConfig.hpp>

void export_StdExpansion();
void export_StdSegExp();
void export_StdTriExp();
void export_StdQuadExp();

BOOST_PYTHON_MODULE(_StdRegions)
{
    np::initialize();

    export_StdExpansion();
    export_StdSegExp();
    export_StdTriExp();
    export_StdQuadExp();
}
