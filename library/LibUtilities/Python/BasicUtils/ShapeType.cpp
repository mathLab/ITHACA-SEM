#include <vector>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <NekPyConfig.hpp>

using namespace Nektar::LibUtilities;

/**
 * @brief Export for ShapeType enumeration.
 */
void export_ShapeType()
{
    NEKPY_WRAP_ENUM(ShapeType, ShapeTypeMap)
}
