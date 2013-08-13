/*
 * TimeIntegrationWrapper.cpp
 *
 *  Created on: 16 Jul 2013
 *      Author: cc
 */

#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar {
namespace LibUtilities {

    TimeIntegrationWrapperFactory& GetTimeIntegrationWrapperFactory()
    {
        typedef Loki::SingletonHolder<TimeIntegrationWrapperFactory,
            Loki::CreateUsingNew,
            Loki::NoDestroy > Type;
        return Type::Instance();
    }

    TimeIntegrationWrapper::TimeIntegrationWrapper()
    {

    }

    void TimeIntegrationWrapper::v_InitObject()
    {

    }


    // --------------
    // IMEX Order 1
    // --------------
    string TimeIntegrationIMEXOrder1::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXOrder1", TimeIntegrationIMEXOrder1::create);

    void TimeIntegrationIMEXOrder1::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXOrder1;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXOrder1);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEX Order 2
    // --------------
    string TimeIntegrationIMEXOrder2::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXOrder2", TimeIntegrationIMEXOrder2::create);

    void TimeIntegrationIMEXOrder2::v_InitObject()
    {
        m_intSteps = 2;
        m_method = eIMEXOrder2;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXOrder1);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eIMEXOrder2);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }


    // --------------
    // IMEX Order 3
    // --------------
    string TimeIntegrationIMEXOrder3::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXOrder3", TimeIntegrationIMEXOrder3::create);

    void TimeIntegrationIMEXOrder3::v_InitObject()
    {
        m_intSteps = 3;
        m_method = eIMEXOrder3;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_3_4_3);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eIMEXdirk_3_4_3);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
        TimeIntegrationSchemeKey       IntKey2(eIMEXOrder3);
        m_integrationScheme[2] = TimeIntegrationSchemeManager()[IntKey2];
    }


    // --------------
    // ForwardEuler
    // --------------
    string TimeIntegrationForwardEuler::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("ForwardEuler", TimeIntegrationForwardEuler::create);

    void TimeIntegrationForwardEuler::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eForwardEuler;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eForwardEuler);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // BackwardEuler
    // --------------
    string TimeIntegrationBackwardEuler::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("BackwardEuler", TimeIntegrationBackwardEuler::create);

    void TimeIntegrationBackwardEuler::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eBackwardEuler;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eBackwardEuler);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // BDFImplicitOrder1
    // --------------
    string TimeIntegrationBDFImplicitOrder1::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("BDFImplicitOrder1", TimeIntegrationBDFImplicitOrder1::create);

    void TimeIntegrationBDFImplicitOrder1::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eBDFImplicitOrder1;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eBDFImplicitOrder1);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // BDFImplicitOrder2
    // --------------
    string TimeIntegrationBDFImplicitOrder2::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("BDFImplicitOrder2", TimeIntegrationBDFImplicitOrder2::create);

    void TimeIntegrationBDFImplicitOrder2::v_InitObject()
    {
        m_intSteps = 2;
        m_method = eBDFImplicitOrder2;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eBackwardEuler);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eBDFImplicitOrder2);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }


    // --------------
    // BDFImplicitOrder2
    // --------------
    string TimeIntegrationRungeKutta2_ImprovedEuler::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("RungeKutta2_ImprovedEuler", TimeIntegrationRungeKutta2_ImprovedEuler::create);

    void TimeIntegrationRungeKutta2_ImprovedEuler::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eRungeKutta2_ImprovedEuler;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eRungeKutta2_ImprovedEuler);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEXdirk_1_1_1
    // --------------
    string TimeIntegrationIMEXdirk_1_1_1::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXdirk_1_1_1", TimeIntegrationIMEXdirk_1_1_1::create);

    void TimeIntegrationIMEXdirk_1_1_1::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXdirk_1_1_1;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_1_1_1);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEXdirk_1_2_1
    // --------------
    string TimeIntegrationIMEXdirk_1_2_1::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXdirk_1_2_1", TimeIntegrationIMEXdirk_1_2_1::create);

    void TimeIntegrationIMEXdirk_1_2_1::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXdirk_1_2_1;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_1_2_1);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEXdirk_1_2_2
    // --------------
    string TimeIntegrationIMEXdirk_1_2_2::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXdirk_1_2_2", TimeIntegrationIMEXdirk_1_2_2::create);

    void TimeIntegrationIMEXdirk_1_2_2::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXdirk_1_2_2;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_1_2_2);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEXdirk_4_4_3
    // --------------
    string TimeIntegrationIMEXdirk_4_4_3::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXdirk_4_4_3", TimeIntegrationIMEXdirk_4_4_3::create);

    void TimeIntegrationIMEXdirk_4_4_3::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXdirk_4_4_3;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_4_4_3);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEXdirk_2_2_2
    // --------------
    string TimeIntegrationIMEXdirk_2_2_2::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXdirk_2_2_2", TimeIntegrationIMEXdirk_2_2_2::create);

    void TimeIntegrationIMEXdirk_2_2_2::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXdirk_2_2_2;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_2_2_2);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEXdirk_2_3_3
    // --------------
    string TimeIntegrationIMEXdirk_2_3_3::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXdirk_2_3_3", TimeIntegrationIMEXdirk_2_3_3::create);

    void TimeIntegrationIMEXdirk_2_3_3::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXdirk_2_3_3;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_2_3_3);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEXdirk_2_3_2
    // --------------
    string TimeIntegrationIMEXdirk_2_3_2::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXdirk_2_3_2", TimeIntegrationIMEXdirk_2_3_2::create);

    void TimeIntegrationIMEXdirk_2_3_2::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXdirk_2_3_2;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_2_3_2);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // IMEXdirk_3_4_3
    // --------------
    string TimeIntegrationIMEXdirk_3_4_3::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXdirk_3_4_3", TimeIntegrationIMEXdirk_3_4_3::create);

    void TimeIntegrationIMEXdirk_3_4_3::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eIMEXdirk_3_4_3;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_3_4_3);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }


    // --------------
    // DIRKOrder2
    // --------------
    string TimeIntegrationDIRKOrder2::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("DIRKOrder2", TimeIntegrationDIRKOrder2::create);

    void TimeIntegrationDIRKOrder2::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eDIRKOrder2;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eDIRKOrder2);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // DIRKOrder3
    // --------------
    string TimeIntegrationDIRKOrder3::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("DIRKOrder3", TimeIntegrationDIRKOrder3::create);

    void TimeIntegrationDIRKOrder3::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eDIRKOrder3;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eDIRKOrder3);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // ClassicalRungeKutta4
    // --------------
    string TimeIntegrationClassicalRungeKutta4::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("ClassicalRungeKutta4", TimeIntegrationClassicalRungeKutta4::create);

    void TimeIntegrationClassicalRungeKutta4::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eClassicalRungeKutta4;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eClassicalRungeKutta4);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // Midpoint
    // --------------
    string TimeIntegrationMidpoint::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("Midpoint", TimeIntegrationMidpoint::create);

    void TimeIntegrationMidpoint::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eMidpoint;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eMidpoint);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // RungeKutta2_ModifiedEuler
    // --------------
    string TimeIntegrationRungeKutta2_ModifiedEuler::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("RungeKutta2_ModifiedEuler", TimeIntegrationRungeKutta2_ModifiedEuler::create);

    void TimeIntegrationRungeKutta2_ModifiedEuler::v_InitObject()
    {
        m_intSteps = 1;
        m_method = eRungeKutta2_ModifiedEuler;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eRungeKutta2_ModifiedEuler);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

    // --------------
    // AdamsBashforthOrder2
    // --------------
    string TimeIntegrationAdamsBashforthOrder2::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("AdamsBashforthOrder2", TimeIntegrationAdamsBashforthOrder2::create);

    void TimeIntegrationAdamsBashforthOrder2::v_InitObject()
    {
        m_intSteps = 2;
        m_method = eAdamsBashforthOrder2;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eForwardEuler);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eAdamsBashforthOrder2);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // AdamsBashforthOrder3
    // --------------
    string TimeIntegrationAdamsBashforthOrder3::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("AdamsBashforthOrder3", TimeIntegrationAdamsBashforthOrder3::create);

    void TimeIntegrationAdamsBashforthOrder3::v_InitObject()
    {
        m_intSteps = 2;
        m_method = eAdamsBashforthOrder3;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eForwardEuler);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eAdamsBashforthOrder3);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // AdamsMoultonOrder2
    // --------------
    string TimeIntegrationAdamsMoultonOrder2::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("AdamsMoultonOrder2", TimeIntegrationAdamsMoultonOrder2::create);

    void TimeIntegrationAdamsMoultonOrder2::v_InitObject()
    {
        m_intSteps = 2;
        m_method = eAdamsMoultonOrder2;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXOrder1);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eAdamsMoultonOrder2);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // IMEXGear
    // --------------
    string TimeIntegrationIMEXGear::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("IMEXGear", TimeIntegrationIMEXGear::create);

    void TimeIntegrationIMEXGear::v_InitObject()
    {
        m_intSteps = 2;
        m_method = eIMEXGear;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_2_2_2);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eIMEXGear);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
    }

    // --------------
    // CNAB
    // --------------
    string TimeIntegrationCNAB::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("CNAB", TimeIntegrationCNAB::create);

    void TimeIntegrationCNAB::v_InitObject()
    {
        m_intSteps = 3;
        m_method = eCNAB;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_3_4_3);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eIMEXdirk_3_4_3);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
        TimeIntegrationSchemeKey       IntKey2(eCNAB);
        m_integrationScheme[2] = TimeIntegrationSchemeManager()[IntKey2];
    }

    // --------------
    // MCNAB
    // --------------
    string TimeIntegrationMCNAB::className  = GetTimeIntegrationWrapperFactory().RegisterCreatorFunction("MCNAB", TimeIntegrationMCNAB::create);

    void TimeIntegrationMCNAB::v_InitObject()
    {
        m_intSteps = 3;
        m_method = eMCNAB;
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eIMEXdirk_3_4_3);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
        TimeIntegrationSchemeKey       IntKey1(eIMEXdirk_3_4_3);
        m_integrationScheme[1] = TimeIntegrationSchemeManager()[IntKey1];
        TimeIntegrationSchemeKey       IntKey2(eMCNAB);
        m_integrationScheme[2] = TimeIntegrationSchemeManager()[IntKey2];
    }
}
}
