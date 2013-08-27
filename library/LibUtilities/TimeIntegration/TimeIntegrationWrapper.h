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
        LIB_UTILITIES_EXPORT TimeIntegrationWrapperFactory&
                                        GetTimeIntegrationWrapperFactory();

        typedef boost::shared_ptr<TimeIntegrationWrapper> TimeIntegrationWrapperSharedPtr;

        class TimeIntegrationWrapper
        {
            public:
                LIB_UTILITIES_EXPORT virtual ~TimeIntegrationWrapper() {}

                LIB_UTILITIES_EXPORT inline void InitObject()
                {
                    v_InitObject();
                }

                LIB_UTILITIES_EXPORT
                TimeIntegrationSolutionSharedPtr InitializeScheme(
                        const NekDouble                      timestep,
                              TimeIntegrationScheme::ConstDoubleArray &y_0,
                        const NekDouble                      time    ,
                        const TimeIntegrationSchemeOperators &op)
                {
                    return m_integrationScheme[m_intSteps-1]->InitializeScheme(timestep, y_0, time, op);
                }

                LIB_UTILITIES_EXPORT
                TimeIntegrationScheme::ConstDoubleArray&
                TimeIntegrate(
                        const int                               timestep,
                        const NekDouble                         delta_t,
                              TimeIntegrationSolutionSharedPtr &solvector,
                        const TimeIntegrationSchemeOperators   &op)
                {
                    return m_integrationScheme[min(timestep,m_intSteps-1)]->TimeIntegrate(delta_t, solvector, op);
                }

                LIB_UTILITIES_EXPORT
                TimeIntegrationMethod GetIntegrationMethod()
                {
                    return m_method;
                }

                LIB_UTILITIES_EXPORT
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

        class TimeIntegrationIMEXdirk_1_1_1 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXdirk_1_1_1>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXdirk_1_1_1>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXdirk_1_1_1() {}

            protected:
                virtual void v_InitObject();

        };


        class TimeIntegrationIMEXdirk_1_2_1 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXdirk_1_2_1>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXdirk_1_2_1>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXdirk_1_2_1() {}

            protected:
                virtual void v_InitObject();

        };

        
        class TimeIntegrationIMEXdirk_1_2_2 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXdirk_1_2_2>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXdirk_1_2_2>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXdirk_1_2_2() {}

            protected:
                virtual void v_InitObject();

        };

        
        class TimeIntegrationIMEXdirk_4_4_3 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXdirk_4_4_3>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXdirk_4_4_3>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXdirk_4_4_3() {}

            protected:
                virtual void v_InitObject();

        };

        
        class TimeIntegrationIMEXdirk_2_2_2 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXdirk_2_2_2>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXdirk_2_2_2>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXdirk_2_2_2() {}

            protected:
                virtual void v_InitObject();

        };

        
        class TimeIntegrationIMEXdirk_2_3_3 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXdirk_2_3_3>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXdirk_2_3_3>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXdirk_2_3_3() {}

            protected:
                virtual void v_InitObject();

        };

        
        class TimeIntegrationIMEXdirk_2_3_2 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXdirk_2_3_2>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXdirk_2_3_2>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXdirk_2_3_2() {}

            protected:
                virtual void v_InitObject();

        };

        
        class TimeIntegrationIMEXdirk_3_4_3 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXdirk_3_4_3>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXdirk_3_4_3>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXdirk_3_4_3() {}

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

        class TimeIntegrationRungeKutta2_ModifiedEuler : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationRungeKutta2_ModifiedEuler>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationRungeKutta2_ModifiedEuler>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationRungeKutta2_ModifiedEuler() {}

            protected:
                virtual void v_InitObject();
        };

        class TimeIntegrationClassicalRungeKutta4 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationClassicalRungeKutta4>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationClassicalRungeKutta4>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationClassicalRungeKutta4() {}

            protected:
                virtual void v_InitObject();
        };

        class TimeIntegrationDIRKOrder2 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationDIRKOrder2>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationDIRKOrder2>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationDIRKOrder2() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationDIRKOrder3 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationDIRKOrder3>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationDIRKOrder3>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationDIRKOrder3() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationMidpoint : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationMidpoint>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationMidpoint>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationMidpoint() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationAdamsBashforthOrder2 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationAdamsBashforthOrder2>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationAdamsBashforthOrder2>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationAdamsBashforthOrder2() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationAdamsBashforthOrder3 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationAdamsBashforthOrder3>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationAdamsBashforthOrder3>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationAdamsBashforthOrder3() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationAdamsMoultonOrder2 : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationAdamsMoultonOrder2>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationAdamsMoultonOrder2>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationAdamsMoultonOrder2() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationIMEXGear : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationIMEXGear>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationIMEXGear>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationIMEXGear() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationCNAB : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationCNAB>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationCNAB>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationCNAB() {}

            protected:
                virtual void v_InitObject();

        };

        class TimeIntegrationMCNAB : public TimeIntegrationWrapper
        {
            public:
                friend class MemoryManager<TimeIntegrationMCNAB>;

                /// Creates an instance of this class
                static TimeIntegrationWrapperSharedPtr create() {
                    TimeIntegrationWrapperSharedPtr p = MemoryManager<TimeIntegrationMCNAB>::AllocateSharedPtr();
                    p->InitObject();
                    return p;
                }
                /// Name of class
                static std::string className;

                virtual ~TimeIntegrationMCNAB() {}

            protected:
                virtual void v_InitObject();

        };
    }
}
#endif /* TIMEINTEGRATIONWRAPPER_H_ */
