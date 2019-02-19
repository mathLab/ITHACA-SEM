ADD_SUBDIRECTORY(user-guide)
ADD_SUBDIRECTORY(developers-guide)

IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tutorial/CMakeLists.txt)
    ADD_SUBDIRECTORY(tutorial)
ENDIF ()
