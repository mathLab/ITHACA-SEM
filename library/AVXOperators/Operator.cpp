#include "Operator.hpp"

OperatorFactory &GetOperatorFactory()
{
    static OperatorFactory tmp;
    return tmp;
}
