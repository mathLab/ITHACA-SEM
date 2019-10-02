#include <LibUtilities/Python/NekPyConfig.hpp>

void export_Field();
void export_Module();
void export_InputModules();
void export_ProcessModules();
void export_OutputModules();

BOOST_PYTHON_MODULE(_FieldUtils)
{
    np::initialize();
    export_Field();
    export_Module();
    export_InputModules();
    export_ProcessModules();
    export_OutputModules();
}
