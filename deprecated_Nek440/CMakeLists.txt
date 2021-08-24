
MACRO(ADD_SOLVER_EXECUTABLE name component source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})
    
    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)

#    SET_COMMON_PROPERTIES(${name})


    TARGET_LINK_LIBRARIES(${name} SolverUtils)

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})
#    INSTALL(TARGETS ${name} RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component} OPTIONAL)
ENDMACRO(ADD_SOLVER_EXECUTABLE name component source)


cmake_minimum_required (VERSION 2.6)
project (ITHACASEM)



FIND_PACKAGE(Nektar++ REQUIRED)
ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})
INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS} ${NEKTAR++_TP_LIBRARY_DIRS})

SET(IncNavierStokesSolverSource    ./EquationSystems/CoupledLinearNS_trafoP.cpp   ./EquationSystems/CoupledLinearNS_TT.cpp    ./EquationSystems/CoupledLinearNS_ROM.cpp      ./EquationSystems/CoupledLinearNS.cpp       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp       ./EquationSystems/IncNavierStokes.cpp       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./EquationSystems/VelocityCorrectionSchemeWeakPressure.cpp       ./EquationSystems/VCSMapping.cpp       ./EquationSystems/Extrapolate.cpp       ./EquationSystems/StandardExtrapolate.cpp
       ./EquationSystems/MappingExtrapolate.cpp       ./EquationSystems/SubSteppingExtrapolate.cpp       ./EquationSystems/SubSteppingExtrapolateWeakPressure.cpp       ./EquationSystems/WeakPressureExtrapolate.cpp
       ./AdvectionTerms/AdjointAdvection.cpp       ./AdvectionTerms/LinearisedAdvection.cpp       ./AdvectionTerms/NavierStokesAdvection.cpp       ./AdvectionTerms/SkewSymmetricAdvection.cpp
       ./AdvectionTerms/AlternateSkewAdvection.cpp       ./AdvectionTerms/NoAdvection.cpp       ./Filters/FilterEnergy.cpp       ./Filters/FilterReynoldsStresses.cpp
       ./Filters/FilterMovingBody.cpp       ./Forcing/ForcingMovingBody.cpp	       ./Forcing/ForcingStabilityCoupledLNS.cpp	       ./myIncNavierStokesSolver.cpp       )


ADD_SOLVER_EXECUTABLE(ITHACASEM solvers ${IncNavierStokesSolverSource})

#add_executable(ITHACASEM myIncNavierStokesSolver.cpp)

TARGET_LINK_LIBRARIES(ITHACASEM ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})


SET(IncNavierStokesSolverSourceDeflation    ./EquationSystems/CoupledLinearNS_trafoP_Deflation.cpp   ./EquationSystems/CoupledLinearNS_TT_Deflation.cpp    ./EquationSystems/CoupledLinearNS_ROM.cpp      ./EquationSystems/CoupledLinearNS.cpp       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp       ./EquationSystems/IncNavierStokes.cpp       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./EquationSystems/VelocityCorrectionSchemeWeakPressure.cpp       ./EquationSystems/VCSMapping.cpp       ./EquationSystems/Extrapolate.cpp       ./EquationSystems/StandardExtrapolate.cpp
       ./EquationSystems/MappingExtrapolate.cpp       ./EquationSystems/SubSteppingExtrapolate.cpp       ./EquationSystems/SubSteppingExtrapolateWeakPressure.cpp       ./EquationSystems/WeakPressureExtrapolate.cpp
       ./AdvectionTerms/AdjointAdvection.cpp       ./AdvectionTerms/LinearisedAdvection.cpp       ./AdvectionTerms/NavierStokesAdvection.cpp       ./AdvectionTerms/SkewSymmetricAdvection.cpp
       ./AdvectionTerms/AlternateSkewAdvection.cpp       ./AdvectionTerms/NoAdvection.cpp       ./Filters/FilterEnergy.cpp       ./Filters/FilterReynoldsStresses.cpp
       ./Filters/FilterMovingBody.cpp       ./Forcing/ForcingMovingBody.cpp	       ./Forcing/ForcingStabilityCoupledLNS.cpp	       ./myIncNavierStokesSolver.cpp       )


ADD_SOLVER_EXECUTABLE(ITHACASEM_Deflation solvers ${IncNavierStokesSolverSourceDeflation})

#add_executable(ITHACASEM myIncNavierStokesSolver.cpp)

TARGET_LINK_LIBRARIES(ITHACASEM_Deflation ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})


