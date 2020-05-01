///////////////////////////////////////////////////////////////////////////////
//
// File: EulerTimeIntegrationSchemes.h
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
// Description: Combined header file for all basic Euler based time integration
// schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_EULER_TIME_INTEGRATION_SCHEME
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_EULER_TIME_INTEGRATION_SCHEME

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationAlgorithmGLM.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGLM.h>

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// Euler

class EulerTimeIntegrationScheme : public TimeIntegrationSchemeGLM
{
public:
    EulerTimeIntegrationScheme(std::string variant, unsigned int order,
                                       std::vector<NekDouble> freeParams) :
        TimeIntegrationSchemeGLM(variant, 1, freeParams)
    {
        boost::ignore_unused(order);

        ASSERTL1(variant == "Backward" || variant == "Forward",
                 "Euler Time integration scheme unknown variant: "
                 + variant + ". Must be 'Backward' or 'Forward'");

        m_integration_phases    = TimeIntegrationAlgorithmGLMVector(1);
        m_integration_phases[0] = TimeIntegrationAlgorithmGLMSharedPtr(
            new TimeIntegrationAlgorithmGLM(this));

        EulerTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0], variant);
    }

    virtual ~EulerTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(order);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
          EulerTimeIntegrationScheme>::AllocateSharedPtr(variant, 1, freeParams);
        return p;
    }

    static std::string className;

    LUE virtual std::string GetFullName() const
    {
        return m_integration_phases[m_integration_phases.size()-1]->m_name;
    }

    LUE virtual std::string GetName() const
    {
        return std::string("Euler");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        if( GetVariant() == "Backward" )
        {
            return 1.0;
        }
        else
        {
            return 2.784;
        }
    }

    LUE static void SetupSchemeData(TimeIntegrationAlgorithmGLMSharedPtr &phase,
                                    std::string variant)
    {
        double A = 0;

        if( variant == "Backward" ) {
            phase->m_schemeType = eDiagonallyImplicit;

            A = 1;
        }
        else if( variant == "Forward" ) {
            phase->m_schemeType = eExplicit;

            A = 0;
        }

        phase->m_variant  = variant;
        phase->m_order = 1;
        phase->m_name = variant + std::string("EulerOrder") +
          std::to_string(phase->m_order);

        phase->m_numsteps  = 1;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, A  );
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 1.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->CheckAndVerify();
    }

}; // end class EulerTimeIntegrator

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class BackwardEulerTimeIntegrationScheme :
    public EulerTimeIntegrationScheme
{
public:
    BackwardEulerTimeIntegrationScheme(std::string variant, unsigned int order,
                                       std::vector<NekDouble> freeParams) :
      EulerTimeIntegrationScheme("Backward", 1, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            EulerTimeIntegrationScheme>::AllocateSharedPtr("Backward", 1, freeParams);
        return p;
    }

    static std::string className;

}; // end class BackwardEulerTimeIntegrationScheme

class ForwardEulerTimeIntegrationScheme :
    public EulerTimeIntegrationScheme
{
public:
    ForwardEulerTimeIntegrationScheme(std::string variant, unsigned int order,
                                       std::vector<NekDouble> freeParams) :
      EulerTimeIntegrationScheme("Forward", 1, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            EulerTimeIntegrationScheme>::AllocateSharedPtr("Forward", 1, freeParams);
        return p;
    }

    static std::string className;

}; // end class ForwardEulerTimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar

#endif
