///////////////////////////////////////////////////////////////////////////////
//
// File: IMEXGearTimeIntegrationScheme.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Header file of time integration scheme wrappers
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#pragma once

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirkTimeIntegrationSchemes.h>

#define LUE LIB_UTILITIES_EXPORT

namespace Nektar
{
namespace LibUtilities
{

class IMEXGearTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXGearTimeIntegrationScheme(std::string variant, unsigned int order,
                                  std::vector<NekDouble> freeParams) :
        TimeIntegrationScheme("", 2, freeParams)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(variant);

        m_integration_phases    = TimeIntegrationSchemeDataVector(2);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));
        m_integration_phases[1] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirkTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0], 2, std::vector<NekDouble> {2, 2});
        IMEXGearTimeIntegrationScheme::SetupSchemeData(m_integration_phases[1]);
    }

    virtual ~IMEXGearTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(variant);

        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXGearTimeIntegrationScheme>::AllocateSharedPtr("", 2, freeParams);

        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEX");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_variant = "Gear";
        phase->m_order = 2;
        phase->m_name = std::string("IMEXGearOrder" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 3;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        NekDouble twothirds = 2.0 / 3.0;

        phase->m_A[0] = Array<TwoD, NekDouble>(phase->m_numstages,
                                               phase->m_numstages, twothirds);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 0.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        phase->m_B[0][0][0] = twothirds;
        phase->m_B[1][2][0] = 1.0;
        phase->m_U[0][0]    = 2.0 * twothirds;
        phase->m_U[0][1]    = -0.5 * twothirds;
        phase->m_U[0][2]    = twothirds;

        phase->m_V[0][0] = 2.0 * twothirds;
        phase->m_V[0][1] = -0.5 * twothirds;
        phase->m_V[0][2] = twothirds;
        phase->m_V[1][0] = 1.0;

        phase->m_numMultiStepValues = 2;
        phase->m_numMultiStepDerivs = 1;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;
        phase->m_timeLevelOffset[1] = 1;
        phase->m_timeLevelOffset[2] = 0;

        phase->CheckAndVerify();
    }

}; // end class IMEXGearTimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
