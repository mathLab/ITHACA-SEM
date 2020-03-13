///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationWrapper.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
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
// Description: Time integration scheme wrappers; one class per time integration
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar {
namespace LibUtilities {

    TimeIntegrationWrapperFactory &GetTimeIntegrationWrapperFactory()
    {
        typedef Loki::SingletonHolder<TimeIntegrationWrapperFactory,
            Loki::CreateUsingNew,
            Loki::NoDestroy,
            Loki::SingleThreaded> Type;
        return Type::Instance();
    }

    TimeIntegrationWrapper::TimeIntegrationWrapper()
    {

    }

    void TimeIntegrationWrapper::v_InitObject()
    {

    }

    // --------------
    // IMEXOrder1
    // --------------
    string TimeIntegrationIMEXOrder1::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXOrder1", TimeIntegrationIMEXOrder1::create);
    void TimeIntegrationIMEXOrder1::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXOrder1);
        m_method       = eIMEXOrder1;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // IMEXOrder2
    // --------------
    string TimeIntegrationIMEXOrder2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXOrder2", TimeIntegrationIMEXOrder2::create);
    void TimeIntegrationIMEXOrder2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXOrder1);
        TimeIntegrationSchemeKey IntKey1(eIMEXOrder2);
        m_method       = eIMEXOrder2;
        m_intSteps     = 2;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // IMEXOrder3
    // --------------
    string TimeIntegrationIMEXOrder3::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXOrder3", TimeIntegrationIMEXOrder3::create);
    void TimeIntegrationIMEXOrder3::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_3_4_3);
        TimeIntegrationSchemeKey IntKey1(eIMEXdirk_3_4_3);
        TimeIntegrationSchemeKey IntKey2(eIMEXOrder3);
        m_method       = eIMEXOrder3;
        m_intSteps     = 3;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
        m_intScheme[2] = TimeIntegrationSchemeManager()[IntKey2];
    }

    // --------------
    // ForwardEuler
    // --------------
    string TimeIntegrationForwardEuler::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "ForwardEuler", TimeIntegrationForwardEuler::create);
    void TimeIntegrationForwardEuler::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eForwardEuler);
        m_method       = eForwardEuler;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // BackwardEuler
    // --------------
    string TimeIntegrationBackwardEuler::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "BackwardEuler", TimeIntegrationBackwardEuler::create);
    void TimeIntegrationBackwardEuler::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eBackwardEuler);
        m_method       = eBackwardEuler;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // BDFImplicitOrder1
    // --------------
    string TimeIntegrationBDFImplicitOrder1::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "BDFImplicitOrder1", TimeIntegrationBDFImplicitOrder1::create);
    void TimeIntegrationBDFImplicitOrder1::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eBDFImplicitOrder1);
        m_method       = eBDFImplicitOrder1;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // BDFImplicitOrder2
    // --------------
    string TimeIntegrationBDFImplicitOrder2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "BDFImplicitOrder2", TimeIntegrationBDFImplicitOrder2::create);
    void TimeIntegrationBDFImplicitOrder2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eBackwardEuler);
        TimeIntegrationSchemeKey IntKey1(eBDFImplicitOrder2);
        m_method       = eBDFImplicitOrder2;
        m_intSteps     = 2;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // IMEXdirk_1_1_1
    // --------------
    string TimeIntegrationIMEXdirk_1_1_1::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXdirk_1_1_1", TimeIntegrationIMEXdirk_1_1_1::create);
    void TimeIntegrationIMEXdirk_1_1_1::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_1_1_1);
        m_method       = eIMEXdirk_1_1_1;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // IMEXdirk_1_2_1
    // --------------
    string TimeIntegrationIMEXdirk_1_2_1::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXdirk_1_2_1", TimeIntegrationIMEXdirk_1_2_1::create);
    void TimeIntegrationIMEXdirk_1_2_1::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_1_2_1);
        m_method       = eIMEXdirk_1_2_1;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // IMEXdirk_1_2_2
    // --------------
    string TimeIntegrationIMEXdirk_1_2_2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXdirk_1_2_2", TimeIntegrationIMEXdirk_1_2_2::create);
    void TimeIntegrationIMEXdirk_1_2_2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_1_2_2);
        m_method       = eIMEXdirk_1_2_2;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // IMEXdirk_4_4_3
    // --------------
    string TimeIntegrationIMEXdirk_4_4_3::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXdirk_4_4_3", TimeIntegrationIMEXdirk_4_4_3::create);
    void TimeIntegrationIMEXdirk_4_4_3::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_4_4_3);
        m_method       = eIMEXdirk_4_4_3;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // IMEXdirk_2_2_2
    // --------------
    string TimeIntegrationIMEXdirk_2_2_2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXdirk_2_2_2", TimeIntegrationIMEXdirk_2_2_2::create);
    void TimeIntegrationIMEXdirk_2_2_2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_2_2_2);
        m_method       = eIMEXdirk_2_2_2;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // IMEXdirk_2_3_3
    // --------------
    string TimeIntegrationIMEXdirk_2_3_3::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXdirk_2_3_3", TimeIntegrationIMEXdirk_2_3_3::create);
    void TimeIntegrationIMEXdirk_2_3_3::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_2_3_3);
        m_method       = eIMEXdirk_2_3_3;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // IMEXdirk_2_3_2
    // --------------
    string TimeIntegrationIMEXdirk_2_3_2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXdirk_2_3_2", TimeIntegrationIMEXdirk_2_3_2::create);
    void TimeIntegrationIMEXdirk_2_3_2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_2_3_2);
        m_method       = eIMEXdirk_2_3_2;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // IMEXdirk_3_4_3
    // --------------
    string TimeIntegrationIMEXdirk_3_4_3::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXdirk_3_4_3", TimeIntegrationIMEXdirk_3_4_3::create);
    void TimeIntegrationIMEXdirk_3_4_3::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_3_4_3);
        m_method       = eIMEXdirk_3_4_3;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // DIRKOrder2
    // --------------
    string TimeIntegrationDIRKOrder2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "DIRKOrder2", TimeIntegrationDIRKOrder2::create);
    void TimeIntegrationDIRKOrder2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eDIRKOrder2);
        m_method       = eDIRKOrder2;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // DIRKOrder3
    // --------------
    string TimeIntegrationDIRKOrder3::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "DIRKOrder3", TimeIntegrationDIRKOrder3::create);
    void TimeIntegrationDIRKOrder3::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eDIRKOrder3);
        m_method       = eDIRKOrder3;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // Midpoint
    // --------------
    string TimeIntegrationMidpoint::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "Midpoint", TimeIntegrationMidpoint::create);
    void TimeIntegrationMidpoint::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eMidpoint);
        m_method       = eMidpoint;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // RungeKutta2
    // --------------
    string TimeIntegrationRungeKutta2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "RungeKutta2", TimeIntegrationRungeKutta2::create);
    void TimeIntegrationRungeKutta2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eRungeKutta2);
        m_method       = eRungeKutta2;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // RungeKutta2_ImprovedEuler
    // --------------
    string TimeIntegrationRungeKutta2_ImprovedEuler::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "RungeKutta2_ImprovedEuler",
            TimeIntegrationRungeKutta2_ImprovedEuler::create);
    void TimeIntegrationRungeKutta2_ImprovedEuler::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eRungeKutta2_ImprovedEuler);
        m_method       = eRungeKutta2_ImprovedEuler;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // RungeKutta2_SSP
    // --------------
    string TimeIntegrationRungeKutta2_SSP::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "RungeKutta2_SSP", TimeIntegrationRungeKutta2_SSP::create);
    void TimeIntegrationRungeKutta2_SSP::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eRungeKutta2_SSP);
        m_method       = eRungeKutta2_SSP;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // RungeKutta3_SSP
    // --------------
    string TimeIntegrationRungeKutta3_SSP::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "RungeKutta3_SSP",
            TimeIntegrationRungeKutta3_SSP::create);
    void TimeIntegrationRungeKutta3_SSP::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eRungeKutta3_SSP);
        m_method       = eRungeKutta3_SSP;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // ClassicalRungeKutta4
    // --------------
    string TimeIntegrationClassicalRungeKutta4::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "ClassicalRungeKutta4",
            TimeIntegrationClassicalRungeKutta4::create);
    void TimeIntegrationClassicalRungeKutta4::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eClassicalRungeKutta4);
        m_method       = eClassicalRungeKutta4;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // RungeKutta4
    // --------------
    string TimeIntegrationRungeKutta4::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "RungeKutta4",
            TimeIntegrationRungeKutta4::create);
    void TimeIntegrationRungeKutta4::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eRungeKutta4);
        m_method       = eRungeKutta4;
        m_intSteps     = 1;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // AdamsBashforthOrder2
    // --------------
    string TimeIntegrationAdamsBashforthOrder2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "AdamsBashforthOrder2",
            TimeIntegrationAdamsBashforthOrder2::create);
    void TimeIntegrationAdamsBashforthOrder2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eForwardEuler);
        TimeIntegrationSchemeKey IntKey1(eAdamsBashforthOrder2);
        m_method       = eAdamsBashforthOrder2;
        m_intSteps     = 2;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // AdamsBashforthOrder3
    // --------------
    string TimeIntegrationAdamsBashforthOrder3::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "AdamsBashforthOrder3",
            TimeIntegrationAdamsBashforthOrder3::create);
    void TimeIntegrationAdamsBashforthOrder3::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eForwardEuler);
        TimeIntegrationSchemeKey IntKey1(eAdamsBashforthOrder3);
        m_method       = eAdamsBashforthOrder3;
        m_intSteps     = 2;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // AdamsMoultonOrder2
    // --------------
    string TimeIntegrationAdamsMoultonOrder2::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "AdamsMoultonOrder2",
            TimeIntegrationAdamsMoultonOrder2::create);
    void TimeIntegrationAdamsMoultonOrder2::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXOrder1);
        TimeIntegrationSchemeKey IntKey1(eAdamsMoultonOrder2);
        m_method       = eAdamsMoultonOrder2;
        m_intSteps     = 2;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // IMEXGear
    // --------------
    string TimeIntegrationIMEXGear::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "IMEXGear", TimeIntegrationIMEXGear::create);
    void TimeIntegrationIMEXGear::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_2_2_2);
        TimeIntegrationSchemeKey IntKey1(eIMEXGear);
        m_method       = eIMEXGear;
        m_intSteps     = 2;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // CNAB
    // --------------
    string TimeIntegrationCNAB::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "CNAB", TimeIntegrationCNAB::create);
    void TimeIntegrationCNAB::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_3_4_3);
        TimeIntegrationSchemeKey IntKey1(eIMEXdirk_3_4_3);
        TimeIntegrationSchemeKey IntKey2(eCNAB);
        m_method       = eCNAB;
        m_intSteps     = 3;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
        m_intScheme[2] = TimeIntegrationSchemeManager()[IntKey2];
    }

    // --------------
    // MCNAB
    // --------------
    string TimeIntegrationMCNAB::className =
        GetTimeIntegrationWrapperFactory().RegisterCreatorFunction(
            "MCNAB", TimeIntegrationMCNAB::create);
    void TimeIntegrationMCNAB::v_InitObject()
    {
        TimeIntegrationSchemeKey IntKey0(eIMEXdirk_3_4_3);
        TimeIntegrationSchemeKey IntKey1(eIMEXdirk_3_4_3);
        TimeIntegrationSchemeKey IntKey2(eMCNAB);
        m_method       = eMCNAB;
        m_intSteps     = 3;
        m_intScheme    = vector<TimeIntegrationSchemeSharedPtr>(m_intSteps);
        m_intScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        m_intScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
        m_intScheme[2] = TimeIntegrationSchemeManager()[IntKey2];
    }
}
}
