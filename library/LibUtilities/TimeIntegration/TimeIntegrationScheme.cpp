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

    os << "Time Integration Scheme: " << rhs.GetName() << ".\n"
       << "        Has " << rhs.m_integration_phases.size() << " phases.\n";

    for (int i = 0; i < rhs.m_integration_phases.size(); i++)
    {
        os << "            - "
           << rhs.m_integration_phases[i]->m_parent->GetName()
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

} // end namespace LibUtilities
} // end namespace NekTar
