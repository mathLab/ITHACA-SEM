/*
 * TimeIntegrationWrapper.h
 *
 *  Created on: 16 Jul 2013
 *      Author: cc
 */

#ifndef TIMEINTEGRATIONWRAPPER_H_
#define TIMEINTEGRATIONWRAPPER_H_

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>

namespace Nektar {
    namespace LibUtilities {

        class TimeIntegrationWrapper;
        /// Datatype of the NekFactory used to instantiate classes derived from
        /// the EquationSystem class.
        typedef NekFactory<
                std::string, TimeIntegrationWrapper
            > TimeIntegrationWrapperFactory;
        TimeIntegrationWrapperFactory& GetTimeIntegrationWrapperFactory();

        typedef boost::shared_ptr<TimeIntegrationWrapper> TimeIntegrationWrapperSharedPtr;

        class TimeIntegrationWrapper
        {
            public:
                virtual ~TimeIntegrationWrapper() {}

                inline void InitObject()
                {
                    v_InitObject();
                }

                TimeIntegrationSolutionSharedPtr InitializeScheme(
                        const NekDouble                      timestep,
                              TimeIntegrationScheme::ConstDoubleArray &y_0,
                        const NekDouble                      time    ,
                        const TimeIntegrationSchemeOperators &op)
                {
                    return m_integrationScheme[m_intSteps-1]->InitializeScheme(timestep, y_0, time, op);
                }

                TimeIntegrationScheme::ConstDoubleArray&
                TimeIntegrate(
                        const int                               timestep,
                        const NekDouble                         delta_t,
                              TimeIntegrationSolutionSharedPtr &solvector,
                        const TimeIntegrationSchemeOperators   &op)
                {
                    return m_integrationScheme[min(timestep,m_intSteps-1)]->TimeIntegrate(delta_t, solvector, op);
                }

                TimeIntegrationMethod GetIntegrationMethod()
                {
                    return m_method;
                }

                unsigned int GetIntegrationSteps()
                {
                    return m_intSteps;
                }

            protected:
                TimeIntegrationMethod m_method;
                int m_intSteps;
                std::vector<TimeIntegrationSchemeSharedPtr> m_integrationScheme;

                /// Constructor
                TimeIntegrationWrapper();

                virtual void v_InitObject();

            private:

        };


        // --------
        // IMEX Order 1
        // --------
        class TimeIntegrationIMEXOrder1;
        typedef boost::shared_ptr<TimeIntegrationIMEXOrder1> TimeIntegrationIMEXOrder1SharedPtr;

        class TimeIntegrationIMEXOrder1 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXOrder1>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXOrder1>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXOrder1() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationIMEXOrder2 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXOrder2>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXOrder2>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXOrder2() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationIMEXOrder3 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXOrder3>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXOrder3>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXOrder3() {}

            protected:
                virtual void v_InitObject();

        };


        class TimeIntegrationForwardEuler : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationForwardEuler>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationForwardEuler>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationForwardEuler() {}

            protected:
                virtual void v_InitObject();

        };


        class TimeIntegrationBackwardEuler : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationBackwardEuler>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationBackwardEuler>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationBackwardEuler() {}

            protected:
                virtual void v_InitObject();

        };


        class TimeIntegrationBDFImplicitOrder1 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationBDFImplicitOrder1>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationBDFImplicitOrder1>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationBDFImplicitOrder1() {}

            protected:
                virtual void v_InitObject();

        };


        class TimeIntegrationBDFImplicitOrder2 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationBDFImplicitOrder2>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationBDFImplicitOrder2>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationBDFImplicitOrder2() {}

            protected:
                virtual void v_InitObject();
        };


        class TimeIntegrationRungeKutta2_ImprovedEuler : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationRungeKutta2_ImprovedEuler>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationRungeKutta2_ImprovedEuler>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationRungeKutta2_ImprovedEuler() {}

            protected:
                virtual void v_InitObject();
        };
    }
}
#endif /* TIMEINTEGRATIONWRAPPER_H_ */
