///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationScheme.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Header file of time integration scheme base class
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationTypes.h>

#define LUE LIB_UTILITIES_EXPORT

namespace Nektar
{
namespace LibUtilities
{

/// Datatype of the NekFactory used to instantiate classes derived from the
/// EquationSystem class.
typedef NekFactory<std::string,
                   TimeIntegrationScheme, std::string, unsigned int,
                   std::vector<NekDouble> > TimeIntegrationSchemeFactory;

// Allows a code to create a TimeIntegrator. Usually used like this:
//
//    LibUtilities::TimeIntegrationSchemeSharedPtr timeIntegrationScheme =
//      LibUtilities::GetTimeIntegrationSchemeFactory().CreateInstance(
//                  "IMEX", "dirk", 3, std::vector<unsigned int>{2,3} );
LUE TimeIntegrationSchemeFactory &GetTimeIntegrationSchemeFactory();

/**
 * @brief Base class for time integration schemes.
 */
class TimeIntegrationScheme
{
public:
    // Access methods
    LUE virtual std::string              GetFullName  () const;
    LUE virtual std::string              GetName      () const = 0;
    LUE virtual std::string              GetVariant   () const = 0;
    LUE virtual unsigned int             GetOrder     () const = 0;
    LUE virtual std::vector< NekDouble > GetFreeParams() const = 0;

    LUE virtual TimeIntegrationSchemeType GetIntegrationSchemeType() const = 0;

    LUE virtual NekDouble GetTimeStability() const = 0;

    LUE virtual unsigned int GetNumIntegrationPhases() const = 0;

    // Gets the solution Vector
    inline virtual const TripleArray &GetSolutionVector() const = 0;
    // Sets the solution Vector
    inline virtual void SetSolutionVector(const int Offset, const DoubleArray &y) = 0;

    // The worker methods
    /**
     * \brief Explicit integration of an ODE.
     *
     * This function explicitely perfroms a signle integration step of the ODE
     * system:
     * \f[
     * \frac{d\boldsymbol{y}}{dt}=\boldsymbol{f}(t,\boldsymbol{y})
     * \f]
     *
     * \param timestep The size of the timestep, i.e. \f$\Delta t\f$.
     * \param f an object of the class FuncType, where FuncType should have a
     * method FuncType::ODEforcing
     *       to evaluate the right hand side \f$f(t,\boldsymbol{y})\f$ of the
     * ODE.
     * \param y on input: the vectors \f$\boldsymbol{y}^{[n-1]}\f$ and
     * \f$t^{[n-1]}\f$ (which corresponds to the
     *    solution at the old time level)
     * \param y on output:  the vectors \f$\boldsymbol{y}^{[n]}\f$ and
     * \f$t^{[n]}\f$ (which corresponds to the
     *    solution at the old new level)
     * \return The actual solution \f$\boldsymbol{y}^{n}\f$ at the new time
     * level
     *    (which in fact is also embedded in the argument y).
     */
    LUE virtual void InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0,
        const NekDouble time, const TimeIntegrationSchemeOperators &op) = 0;

    LUE virtual ConstDoubleArray &TimeIntegrate(
        const int timestep, const NekDouble delta_t,
        const TimeIntegrationSchemeOperators &op) = 0;

    LUE virtual void print(std::ostream &os) const = 0;

    // Friend classes
    LUE friend std::ostream &operator<<(std::ostream &os,
        const TimeIntegrationScheme &rhs);
    LUE friend std::ostream &operator<<(std::ostream &os,
        const TimeIntegrationSchemeSharedPtr &rhs);

protected:
    // These methods should never be used directly, only used by child classes.
    LUE TimeIntegrationScheme(std::string variant, unsigned int order,
                              std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant, order, freeParams);
    }

    LUE TimeIntegrationScheme(const TimeIntegrationScheme &in)
    {
        boost::ignore_unused(in);

        NEKERROR(ErrorUtil::efatal, "Copy Constructor for the "
                                    "TimeIntegrationScheme class should not be "
                                    "called");
    }

    virtual ~TimeIntegrationScheme()
    {
    }

}; // end class TimeIntegrationScheme

LUE std::ostream &operator<<(std::ostream &os,
                             const TimeIntegrationScheme &rhs);

LUE std::ostream &operator<<(std::ostream &os,
                             const TimeIntegrationSchemeSharedPtr &rhs);

} // end of namespace LibUtilities
} // end of namespace Nektar
