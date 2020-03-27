///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationScheme.cpp
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
// Description: implementation of time integration key class
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>

namespace Nektar
{
namespace LibUtilities
{

TimeIntegrationSchemeFactory &GetTimeIntegrationSchemeFactory()
{
    static TimeIntegrationSchemeFactory instance;
    return instance;
}

std::ostream &operator<<(std::ostream &os,
                         const TimeIntegrationSchemeSharedPtr &rhs)
{
    os << *rhs.get();
    return os;
}

std::ostream &operator<<(std::ostream &os, const TimeIntegrationScheme &rhs)
{

    os << "Time Integration Scheme: " << rhs.GetFullName() << ".\n"
       << "        Has " << rhs.m_integration_phases.size() << " phases.\n";

    for (int i = 0; i < rhs.m_integration_phases.size(); i++)
    {
        os << "            - "
           << rhs.m_integration_phases[i]->m_parent->GetFullName()
           << "\n";
    }
    return os;
}

// Access Methods
std::string TimeIntegrationScheme::GetFullName () const
{
    return GetName() + GetVariant() + "Order" + std::to_string(GetOrder());
}

std::string TimeIntegrationScheme::GetVariant() const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    return m_integration_phases[m_integration_phases.size() - 1]->m_variant;
}

unsigned int TimeIntegrationScheme::GetOrder() const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    return m_integration_phases[m_integration_phases.size() - 1]->m_order;
}

std::vector< NekDouble > TimeIntegrationScheme::GetFreeParams() const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    return m_integration_phases[m_integration_phases.size() - 1]->m_freeParams;
}

TimeIntegrationSchemeType TimeIntegrationScheme::GetIntegrationSchemeType()
    const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    return m_integration_phases[m_integration_phases.size() - 1]->m_schemeType;
}

unsigned int TimeIntegrationScheme::GetNumIntegrationPhases() const
{
    return m_integration_phases.size();
}

// The worker methods
TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr TimeIntegrationScheme::
    InitializeScheme(const NekDouble deltaT,
                     TimeIntegrationScheme::ConstDoubleArray &y_0,
                     const NekDouble time,
                     const TimeIntegrationSchemeOperators &op)
{
    return m_integration_phases.back()->InitializeData(deltaT, y_0, time, op);
}

TimeIntegrationScheme::ConstDoubleArray &TimeIntegrationScheme::TimeIntegrate(
    const int timestep,
    const NekDouble delta_t,
    TimeIntegrationSolutionSharedPtr &solvector,
    const TimeIntegrationSchemeOperators &op)
{
    int nPhases = m_integration_phases.size();

    TimeIntegrationSchemeDataSharedPtr &data =
        m_integration_phases[std::min(timestep, nPhases - 1)];

    return data->TimeIntegrate(delta_t, solvector, op);
}

// Methods specific to exponential integration schemes.
void TimeIntegrationScheme::SetExponentialCoefficients(
    Array<OneD, std::complex<NekDouble>> &Lambda)
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    // Assumption: the one-dimensional Lambda matrix is a diagonal
    // matrix thus values are non zero if and only i=j. As such, the
    // diagonal Lambda values are stored an array of complex numbers.

    // Assume that each phase is an exponential integrator.
    for (int i = 0; i < m_integration_phases.size(); i++)
    {
        m_integration_phases[i]->m_L = Lambda;

        // Anytime the coefficents are updated reset the nVars to be
        // assured that the exponential matrices are recalculated
        // (e.g. the number of variables may remain the same but the
        // coefficients have changed).
        m_integration_phases[i]->m_lastNVars = 0;
    }
}

void TimeIntegrationScheme::SetupSchemeExponentialData(TimeIntegrationSchemeData *phase, NekDouble deltaT) const
{
    boost::ignore_unused(phase, deltaT);

    ASSERTL0(false, "No SetupSchemeExponentialData method for scheme " +
             GetFullName());
}

inline NekDouble TimeIntegrationScheme::factorial(unsigned int n) const
{
  return (n == 1 || n == 0) ? 1 : n * factorial(n - 1);
}

std::complex<NekDouble>
TimeIntegrationScheme::phi_function(const unsigned int order,
                                    const std::complex<NekDouble> z) const
{
    // Central to the implementation of exponential integrators is the
    // evaluation of exponential-like functions, commonly denoted by φ
    // functions. It is convenient to define φ0(z) = e^z, in which case the
    // functions obey the recurrence relation.

    // 0: exp(z);
    // 1: (exp(z)     - 1.0) / (z);
    // 2: (exp(z) - z - 1.0) / (z * z);

    if( z == 0.0 )
    {
        return 1.0 / factorial( order );
    }

    if( order == 0 )
    {
        return exp( z );
    }
    else
    {
        return (phi_function( order-1, z) -
                1.0 / factorial( order-1 ) ) / z;
    }


    return 0;
}

} // end namespace LibUtilities
} // end namespace NekTar
