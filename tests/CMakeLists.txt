CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/tests/Tester.cpp.in
               ${CMAKE_BINARY_DIR}/tests/Tester.cpp
	       @ONLY)
SET(TESTER_SOURCES
  Metric.cpp
  MetricEigenvalue.cpp
  MetricFile.cpp
  MetricL2.cpp
  MetricLInf.cpp
  MetricPrecon.cpp
  MetricRegex.cpp
  TestData.cpp
  ${CMAKE_BINARY_DIR}/tests/Tester.cpp
  sha1.cpp
)

SET(TESTER_HEADERS
  Metric.h
  MetricEigenvalue.h
  MetricFile.h
  MetricL2.h
  MetricLInf.h
  MetricPrecon.h
  MetricRegex.h
  TestData.h
  Tester.h
  sha1.h
)

ADD_DEFINITIONS(-DBUILD_PATH="${CMAKE_BINARY_DIR}")

IF(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  ADD_DEFINITIONS(-DRELWITHDEBINFO)
ENDIF()


ADD_NEKTAR_EXECUTABLE(Tester COMPONENT tests SOURCES ${TESTER_SOURCES})

TARGET_LINK_LIBRARIES(Tester
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_REGEX_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  optimized ${TINYXML_LIBRARY} debug ${TINYXML_LIBRARY}
)
ADD_DEPENDENCIES(Tester boost tinyxml-2.6.2)

IF( ${CMAKE_SYSTEM} MATCHES "Linux.*")
  TARGET_LINK_LIBRARIES(Tester optimized rt debug rt)
  SET_PROPERTY(TARGET Tester APPEND_STRING PROPERTY COMPILE_FLAGS " -pthread")
  SET_PROPERTY(TARGET Tester APPEND_STRING PROPERTY LINK_FLAGS " -pthread")
ENDIF()

IF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
  SET_TARGET_PROPERTIES(Tester
    PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup -Wl,-rpath,${Boost_LIBRARY_DIRS}")
ENDIF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")

SET_TARGET_PROPERTIES(Tester
    PROPERTIES COMPILE_FLAGS "-UNEKTAR_USE_MPI")
