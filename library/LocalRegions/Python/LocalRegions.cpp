#include <NekPyConfig.hpp>

void export_Expansion();
void export_SegExp();
void export_TriExp();
void export_QuadExp();

BOOST_PYTHON_MODULE(_LocalRegions)
{
    np::initialize();

    export_Expansion();
    export_SegExp();
    export_TriExp();
    export_QuadExp();
}
