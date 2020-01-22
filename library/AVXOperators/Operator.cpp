#include "Operator.hpp"

namespace Nektar {
namespace AVX {

OperatorFactory &GetOperatorFactory()
{
    static OperatorFactory tmp;
    return tmp;
}

}
}