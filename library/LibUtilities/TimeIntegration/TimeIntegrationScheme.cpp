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
#include <cmath>
#include <algorithm>
using namespace std;

#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>

namespace Nektar
{
namespace LibUtilities
{
// !!! Always make sure that this matches TimeIntegrationMethod enum... !!!
//
const char *const TimeIntegrationScheme::TimeIntegrationMethodMap[36] = {
    "NoTimeIntegrationMethod",
    "AdamsBashforthOrder1",
    "AdamsBashforthOrder2",
    "AdamsBashforthOrder3",
    "AdamsBashforthOrder4",
    "AdamsMoultonOrder1",
    "AdamsMoultonOrder2",
    "BDFImplicitOrder1",
    "BDFImplicitOrder2",
    "ClassicalRungeKutta4",
    "RungeKutta4",
    "RungeKutta5",
    "RungeKutta3_SSP",
    "RungeKutta2_ImprovedEuler",
    "RungeKutta2_SSP",
    "ForwardEuler",
    "BackwardEuler",
    "IMEXOrder1",
    "IMEXOrder2",
    "IMEXOrder3",
    "IMEXOrder4",
    "Midpoint",
    "RungeKutta2",
    "DIRKOrder2",
    "DIRKOrder3",
    "CNAB",
    "IMEXGear",
    "MCNAB",
    "IMEXdirk_1_1_1",
    "IMEXdirk_1_2_1",
    "IMEXdirk_1_2_2",
    "IMEXdirk_2_2_2",
    "IMEXdirk_2_3_2",
    "IMEXdirk_2_3_3",
    "IMEXdirk_3_4_3",
    "IMEXdirk_4_4_3"};

TimeIntegrationMethod TimeIntegrationScheme::methodFromName(const string &name)
{
    if (name == "AdamsBashforthOrder1")
    {
        return eAdamsBashforthOrder1;
    }
    else if (name == "AdamsBashforthOrder2")
    {
        return eAdamsBashforthOrder2;
    }
    else if (name == "AdamsBashforthOrder3")
    {
        return eAdamsBashforthOrder3;
    }
    else if (name == "AdamsBashforthOrder4")
    {
        return eAdamsBashforthOrder4;
    }
    else if (name == "AdamsMoultonOrder1")
    {
        return eAdamsMoultonOrder1;
    }
    else if (name == "AdamsMoultonOrder2")
    {
        return eAdamsMoultonOrder2;
    }
    else if (name == "BDFImplicitOrder1")
    {
        return eBDFImplicitOrder1;
    }
    else if (name == "BDFImplicitOrder2")
    {
        return eBDFImplicitOrder2;
    }
    else if (name == "ClassicalRungeKutta4")
    {
        return eClassicalRungeKutta4;
    }
    else if (name == "RungeKutta4")
    {
        return eRungeKutta4;
    }
    else if (name == "RungeKutta5")
    {
        return eRungeKutta5;
    }
    else if (name == "RungeKutta3_SSP")
    {
        return eRungeKutta3_SSP;
    }
    else if (name == "RungeKutta2_ImprovedEuler")
    {
        return eRungeKutta2_ImprovedEuler;
    }
    else if (name == "RungeKutta2_SSP")
    {
        return eRungeKutta2_SSP;
    }
    else if (name == "ForwardEuler")
    {
        return eForwardEuler;
    }
    else if (name == "BackwardEuler")
    {
        return eBackwardEuler;
    }
    else if (name == "IMEXOrder1")
    {
        return eIMEXOrder1;
    }
    else if (name == "IMEXOrder2")
    {
        return eIMEXOrder2;
    }
    else if (name == "IMEXOrder3")
    {
        return eIMEXOrder3;
    }
    else if (name == "IMEXOrder4")
    {
        return eIMEXOrder4;
    }
    else if (name == "Midpoint")
    {
        return eMidpoint;
    }
    else if (name == "RungeKutta2")
    {
        return eRungeKutta2;
    }
    else if (name == "DIRKOrder2")
    {
        return eDIRKOrder2;
    }
    else if (name == "DIRKOrder3")
    {
        return eDIRKOrder3;
    }
    else if (name == "CNAB")
    {
        return eCNAB;
    }
    else if (name == "IMEXGear")
    {
        return eIMEXGear;
    }
    else if (name == "MCNAB")
    {
        return eMCNAB;
    }
    else if (name == "IMEXdirk_1_1_1")
    {
        return eIMEXdirk_1_1_1;
    }
    else if (name == "IMEXdirk_1_2_1")
    {
        return eIMEXdirk_1_2_1;
    }
    else if (name == "IMEXdirk_1_2_2")
    {
        return eIMEXdirk_1_2_2;
    }
    else if (name == "IMEXdirk_2_2_2")
    {
        return eIMEXdirk_2_2_2;
    }
    else if (name == "IMEXdirk_2_3_2")
    {
        return eIMEXdirk_2_3_2;
    }
    else if (name == "IMEXdirk_2_3_3")
    {
        return eIMEXdirk_2_3_3;
    }
    else if (name == "IMEXdirk_3_4_3")
    {
        return eIMEXdirk_3_4_3;
    }
    else if (name == "IMEXdirk_4_4_3")
    {
        return eIMEXdirk_4_4_3;
    }
    else
    {
        string msg =
            "'" + name +
            "' is not a known TimeIntegrationMethod. (Check spelling?)";
        NEKERROR(ErrorUtil::efatal, msg);
        return eNoTimeIntegrationMethod;
    }
}

string TimeIntegrationScheme::nameFromMethod(const TimeIntegrationMethod method)
{
    return TimeIntegrationMethodMap[method];
}

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
    os << "Time Integration Scheme: "
       << TimeIntegrationScheme::nameFromMethod(rhs.GetIntegrationMethod())
       << ".\n";
    os << "        Has " << rhs.m_integration_phases.size() << " phases.\n";
    for (int i = 0; i < rhs.m_integration_phases.size(); i++)
    {
        os << "            - "
           << TimeIntegrationScheme::nameFromMethod(
                  rhs.m_integration_phases[i]->GetIntegrationMethod())
           << "\n";
    }
    return os;
}

TimeIntegrationScheme::ConstDoubleArray &TimeIntegrationScheme::TimeIntegrate(
    const int timestep, const NekDouble delta_t,
    TimeIntegrationSolutionSharedPtr &solvector,
    const TimeIntegrationSchemeOperators &op)
{
    int phases = GetNumIntegrationPhases();
    TimeIntegrationSchemeDataSharedPtr &data =
        m_integration_phases[std::min(timestep, phases - 1)];
    return data->TimeIntegrate(delta_t, solvector, op);
}

TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr TimeIntegrationScheme::
    InitializeScheme(const NekDouble deltaT,
                     TimeIntegrationScheme::ConstDoubleArray &y_0,
                     const NekDouble time,
                     const TimeIntegrationSchemeOperators &op)
{
    return m_integration_phases.back()->InitializeData(deltaT, y_0, time, op);
}

TimeIntegrationSchemeType TimeIntegrationScheme::GetIntegrationSchemeType()
    const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")
    return m_integration_phases[m_integration_phases.size() - 1]->m_schemeType;
}

void TimeIntegrationScheme::SetExponentialCoefficients(Array<TwoD, NekDouble> &Lambda)
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    // Assumption: the two-dimensional Lambda matrix is a diagonal
    // matrix thus values are non zero if and only i=j. As such, the
    // diagonal Lambda values are stored as two vectors so to
    // accomodate complex numbers lambda[0] real, lambda[1] imaginary.

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
    boost::ignore_unused(phase);
    boost::ignore_unused(deltaT);

    ASSERTL0(false, "No SetupSchemeExponentialData method for scheme " +
             GetName());
};

NekDouble TimeIntegrationScheme::exp_function(NekDouble deltaT,
                                              NekDouble L_Real,
                                              NekDouble L_Imaginary) const
{
    boost::ignore_unused( L_Imaginary );

    return exp( deltaT * L_Real );
}

NekDouble TimeIntegrationScheme::psi_function(unsigned int i,
                                              NekDouble deltaT,
                                              NekDouble L_Real,
                                              NekDouble L_Imaginary) const
{
    NekDouble z = deltaT * L_Real;

    NekDouble expZ = exp_function(deltaT, L_Real, L_Imaginary);

    switch( i )
    {
    case 1:
        return  (expZ - 1.0) / (z);
        break;
    case 2:
        return (expZ - z - 1.0) / (z * z);
        break;
    case 3:
        return (expZ - (z*z)/2.0 - z - 1.0) / (z * z * z);
        break;
    default:
        ASSERTL0(false,
                 "No phi function defined for index " + std::to_string(i) );
        break;
    }

    return 0;
}

} // end namespace LibUtilities
} // end namespace NekTar
