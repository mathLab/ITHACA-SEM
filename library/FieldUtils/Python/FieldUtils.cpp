#include <LibUtilities/Python/NekPyConfig.hpp>

void export_Field();
void export_Module();

BOOST_PYTHON_MODULE(_FieldUtils)
{
    np::initialize();
    export_Field();
    export_Module();
}
