#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERINTERFACES_HPP
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERINTERFACES_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace SolverUtils
{

class FluidInterface
{
public:
    /// Extract array with velocity from physfield
    SOLVER_UTILS_EXPORT virtual void GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> >       &velocity) = 0;

    SOLVER_UTILS_EXPORT virtual bool HasConstantDensity() = 0;

    /// Extract array with density from physfield
    SOLVER_UTILS_EXPORT virtual void GetDensity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, NekDouble>                     &density) = 0;

    /// Extract array with pressure from physfield
    SOLVER_UTILS_EXPORT virtual void GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, NekDouble>                     &pressure) = 0;
};

}
}

#endif
