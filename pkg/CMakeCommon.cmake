macro (add_deb_package)
    set(options "")
    set(oneValueArgs NAME SUMMARY DESCRIPTION)
    set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
    cmake_parse_arguments(ARG "${options}"
            "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${ARG_NAME}-deb)
    set(PKG_PROJ ${ARG_NAME})
    set(PKG_NAME ${ARG_NAME})
    #    set(PKG_DESC_FILE
    #    "${CMAKE_CURRENT_SOURCE_DIR}/desc/libnektar++-utilities.txt")
    set(PKG_SUMM ${ARG_SUMMARY})
    set(PKG_DESC ${ARG_DESCRIPTION})
    configure_file(CMakeListsDpkg.txt.in
                ${BUILD_DIR}/CMakeLists.txt @ONLY)
    add_custom_target(
        pkg-deb-${ARG_NAME}
        rm -f ${BUILD_DIR}/CPackConfig.cmake
        COMMAND ${CMAKE_COMMAND} .
        COMMAND ${CMAKE_CPACK_COMMAND}
        WORKING_DIRECTORY ${BUILD_DIR}
    )
    add_dependencies(pkg-deb pkg-deb-${ARG_NAME})
endmacro (add_deb_package)

macro (add_rpm_package)
    set(options "")
    set(oneValueArgs NAME SUMMARY DESCRIPTION)
    set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
    cmake_parse_arguments(ARG "${options}"
            "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${ARG_NAME}-rpm)
    set(PKG_PROJ ${ARG_NAME})
    set(PKG_NAME ${ARG_NAME})
    #    set(PKG_DESC_FILE
    #    "${CMAKE_CURRENT_SOURCE_DIR}/desc/libnektar++-utilities.txt")
    set(PKG_SUMM ${ARG_SUMMARY})
    set(PKG_DESC ${ARG_DESCRIPTION})
    configure_file(CMakeListsRpm.txt.in
                ${BUILD_DIR}/CMakeLists.txt @ONLY)
    add_custom_target(
        pkg-rpm-${ARG_NAME}
        rm -f ${BUILD_DIR}/CPackConfig.cmake
        COMMAND ${CMAKE_COMMAND} .
        COMMAND ${CMAKE_CPACK_COMMAND}
        WORKING_DIRECTORY ${BUILD_DIR}
    )
    add_dependencies(pkg-rpm pkg-rpm-${ARG_NAME})
endmacro (add_rpm_package)


