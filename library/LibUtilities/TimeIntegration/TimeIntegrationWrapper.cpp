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
        m_integrationScheme = std::vector<TimeIntegrationSchemeSharedPtr> (m_intSteps);
        TimeIntegrationSchemeKey       IntKey0(eRungeKutta2_ImprovedEuler);
        m_integrationScheme[0] = TimeIntegrationSchemeManager()[IntKey0];
    }

}
}
