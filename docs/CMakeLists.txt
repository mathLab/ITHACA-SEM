ADD_SUBDIRECTORY(user-guide)

IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/developer-guide/CMakeLists.txt)
    ADD_SUBDIRECTORY(developer-guide)
ENDIF ()

IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tutorial/CMakeLists.txt)
    ADD_SUBDIRECTORY(tutorial)
ENDIF ()
