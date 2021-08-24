#include "Operator.hpp"

namespace Nektar {
namespace MatrixFree {

OperatorFactory &GetOperatorFactory()
{
    static OperatorFactory tmp;
    return tmp;
}

std::string GetOpstring(LibUtilities::ShapeType shape, bool deformed)
{

    std::string op_string = "_";

    if(shape == LibUtilities::eSegment)
    {
        op_string += "Seg";
    }
    else if(shape == LibUtilities::eTriangle)
    {
        op_string += "Tri";
    }
    else if(shape == LibUtilities::eQuadrilateral)
    {
        op_string += "Quad";
    }
    else if(shape == LibUtilities::eTetrahedron)
    {
        op_string += "Tet";
    }
    else if(shape == LibUtilities::ePyramid)
    {
        op_string += "Pyr";
    }
    else if(shape == LibUtilities::ePrism)
    {
        op_string += "Prism";
    }
    else if(shape == LibUtilities::eHexahedron)
    {
        op_string += "Hex";
    }

    if (deformed)
    {
        op_string += "_Deformed";
    }
    else
    {
        op_string += "_Regular";
    }

    return op_string;
}


}
}
