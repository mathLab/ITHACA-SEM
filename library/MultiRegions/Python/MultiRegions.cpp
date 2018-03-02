#include <NekPyConfig.hpp>

void export_ExpList();
void export_ExpList2D();

BOOST_PYTHON_MODULE(_MultiRegions)
{
    np::initialize();

    export_ExpList();
    export_ExpList2D();
}
