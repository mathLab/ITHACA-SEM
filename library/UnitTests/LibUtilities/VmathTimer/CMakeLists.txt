SET(VmathTimerSources
    main.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(VmathTimer
    COMPONENT unit-test DEPENDS LibUtilities SOURCES ${VmathTimerSources})
