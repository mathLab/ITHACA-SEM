macro (add_deb_package)
    set(options "")
    set(oneValueArgs NAME SUMMARY DESCRIPTION)
    set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
    cmake_parse_arguments(PKG "${options}"
            "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PKG_NAME}-deb)
    configure_file(CMakeListsDpkg.txt.in
                ${BUILD_DIR}/CMakeLists.txt @ONLY)
    add_custom_target(
        pkg-deb-${PKG_NAME}
        rm -f ${BUILD_DIR}/CPackConfig.cmake
        COMMAND ${CMAKE_COMMAND} .
        COMMAND ${CMAKE_CPACK_COMMAND}
        WORKING_DIRECTORY ${BUILD_DIR}
    )
    add_dependencies(pkg-deb pkg-deb-${PKG_NAME})
endmacro (add_deb_package)

macro (add_rpm_package)
    set(options "")
    set(oneValueArgs NAME SUMMARY DESCRIPTION)
    set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
    cmake_parse_arguments(PKG "${options}"
            "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PKG_NAME}-rpm)
    configure_file(CMakeListsRpm.txt.in
                ${BUILD_DIR}/CMakeLists.txt @ONLY)
    add_custom_target(
        pkg-rpm-${PKG_NAME}
        rm -f ${BUILD_DIR}/CPackConfig.cmake
        COMMAND ${CMAKE_COMMAND} .
        COMMAND ${CMAKE_CPACK_COMMAND}
        WORKING_DIRECTORY ${BUILD_DIR}
    )
    add_dependencies(pkg-rpm pkg-rpm-${PKG_NAME})
endmacro (add_rpm_package)

macro (add_tgz_package)
    set(options "")
    set(oneValueArgs NAME SUMMARY DESCRIPTION)
    set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
    cmake_parse_arguments(PKG "${options}"
            "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PKG_NAME}-tgz)
    configure_file(CMakeListsTgz.txt.in
                ${BUILD_DIR}/CMakeLists.txt @ONLY)
    add_custom_target(
        pkg-tgz-${PKG_NAME}
        rm -f ${BUILD_DIR}/CPackConfig.cmake
        COMMAND ${CMAKE_COMMAND} .
        COMMAND ${CMAKE_CPACK_COMMAND}
        WORKING_DIRECTORY ${BUILD_DIR}
    )
    add_dependencies(pkg-tgz pkg-tgz-${PKG_NAME})
endmacro (add_tgz_package)


