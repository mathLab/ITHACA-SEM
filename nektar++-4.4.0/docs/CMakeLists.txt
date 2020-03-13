ADD_SUBDIRECTORY(user-guide)
ADD_SUBDIRECTORY(developer-guide)

IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tutorial/CMakeLists.txt)
    ADD_SUBDIRECTORY(tutorial)
ENDIF ()
