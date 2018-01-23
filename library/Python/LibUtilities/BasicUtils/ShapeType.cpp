#include <vector>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <NekPyConfig.hpp>

using namespace Nektar::LibUtilities;

/**
 * @brief Export for ShapeType enumeration.
 */
void export_ShapeType()
{
    py::enum_<ShapeType> tmp("ShapeType");
    for (int i = 0; i < (int)SIZE_ShapeType; ++i)
    {
        tmp.value(ShapeTypeMap[i], (ShapeType) i);
    }
    tmp.export_values();
}
