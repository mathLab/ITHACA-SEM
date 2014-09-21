function (find_lib_files PKG_INSTALL_LIBS PKG_INSTALL_LIBS_FILES)
    # Find library file and add the versioned form of each library
    set(PKG_INSTALL_LIBS_FILES)
    foreach(l ${PKG_INSTALL_LIBS})
        get_target_property(TARGET_LOCATION ${l} LOCATION)
        if (NOT TARGET_LOCATION)
            message(FATAL_ERROR "Target '${l}' could not be found.")
        endif ()
        list(APPEND PKG_INSTALL_LIBS_FILES ${TARGET_LOCATION})
        if (APPLE)
            list(APPEND PKG_INSTALL_LIBS_FILES 
                        ${TARGET_LOCATION}.${VERSION_MAJOR_MINOR})
        else ()
            list(APPEND PKG_INSTALL_LIBS_FILES 
                        ${TARGET_LOCATION}.${NEKTAR_VERSION})
        endif()
    endforeach()
    set(PKG_INSTALL_LIBS_FILES ${PKG_INSTALL_LIBS_FILES} PARENT_SCOPE)
endfunction ()

function (find_bin_files PKG_INSTALL_BINS PKG_INSTALL_BINS_FILES)
    # Find binary files
    set(PKG_INSTALL_BINS_FILES)
    foreach(b ${PKG_INSTALL_BINS})
        get_target_property(TARGET_LOCATION ${b} LOCATION)
        if (NOT TARGET_LOCATION)
            message(FATAL_ERROR "Target '${b}' could not be found.")
        endif ()
        list(APPEND PKG_INSTALL_BINS_FILES ${TARGET_LOCATION})
        list(APPEND PKG_INSTALL_BINS_FILES 
                        ${TARGET_LOCATION}-${NEKTAR_VERSION})
    endforeach()
    set(PKG_INSTALL_BINS_FILES ${PKG_INSTALL_BINS_FILES} PARENT_SCOPE)
endfunction ()

macro (add_deb_package)
    if (DPKG)
        set(options "")
        set(oneValueArgs NAME SUMMARY DESCRIPTION)
        set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
        cmake_parse_arguments(PKG "${options}"
            "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

        set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PKG_NAME}-deb)

        find_lib_files("${PKG_INSTALL_LIBS}" PKG_INSTALL_LIBS_FILES)
        find_bin_files("${PKG_INSTALL_BINS}" PKG_INSTALL_BINS_FILES)

        # Configure project for this package
        configure_file(CMakeListsDpkg.txt.in
                ${BUILD_DIR}/CMakeLists.txt @ONLY)

        add_custom_target(
            pkg-deb-${PKG_NAME}
            rm -f ${BUILD_DIR}/CPackConfig.cmake
            COMMAND ${CMAKE_COMMAND} .
            COMMAND ${CMAKE_CPACK_COMMAND} --config CPackConfig.cmake
            WORKING_DIRECTORY ${BUILD_DIR}
        )
        if (PKG_INSTALL_LIBS OR PKG_INSTALL_BINS)
            add_dependencies(pkg-deb-${PKG_NAME}
                ${PKG_INSTALL_LIBS} ${PKG_INSTALL_BINS})
        endif ()
        add_dependencies(pkg-deb pkg-deb-${PKG_NAME})
    endif ()
endmacro (add_deb_package)

macro (add_rpm_package)
    if (RPMBUILD)
        set(options "")
        set(oneValueArgs NAME SUMMARY DESCRIPTION)
        set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
        cmake_parse_arguments(PKG "${options}"
            "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

        set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PKG_NAME}-rpm)

        find_lib_files("${PKG_INSTALL_LIBS}" PKG_INSTALL_LIBS_FILES)
        find_bin_files("${PKG_INSTALL_BINS}" PKG_INSTALL_BINS_FILES)

        configure_file(CMakeListsRpm.txt.in
                ${BUILD_DIR}/CMakeLists.txt @ONLY)
        add_custom_target(
            pkg-rpm-${PKG_NAME}
            rm -f ${BUILD_DIR}/CPackConfig.cmake
            COMMAND ${CMAKE_COMMAND} .
            COMMAND ${CMAKE_CPACK_COMMAND} --config CPackConfig.cmake
            WORKING_DIRECTORY ${BUILD_DIR}
        )
        if (PKG_INSTALL_LIBS OR PKG_INSTALL_BINS)
            add_dependencies(pkg-rpm-${PKG_NAME}
                ${PKG_INSTALL_LIBS} ${PKG_INSTALL_BINS})
        endif ()
        add_dependencies(pkg-rpm pkg-rpm-${PKG_NAME})
    endif ()
endmacro (add_rpm_package)

macro (add_tgz_package)
    set(options "")
    set(oneValueArgs NAME SUMMARY DESCRIPTION)
    set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
    cmake_parse_arguments(PKG "${options}"
            "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PKG_NAME}-tgz)

    find_lib_files("${PKG_INSTALL_LIBS}" PKG_INSTALL_LIBS_FILES)
    find_bin_files("${PKG_INSTALL_BINS}" PKG_INSTALL_BINS_FILES)

    configure_file(CMakeListsTgz.txt.in
                ${BUILD_DIR}/CMakeLists.txt @ONLY)
    add_custom_target(
        pkg-tgz-${PKG_NAME}
        rm -f ${BUILD_DIR}/CPackConfig.cmake
        COMMAND ${CMAKE_COMMAND} .
        COMMAND ${CMAKE_CPACK_COMMAND} --config CPackConfig.cmake
        WORKING_DIRECTORY ${BUILD_DIR}
    )
    if (PKG_INSTALL_LIBS OR PKG_INSTALL_BINS)
        add_dependencies(pkg-tgz-${PKG_NAME}
            ${PKG_INSTALL_LIBS} ${PKG_INSTALL_BINS})
    endif()
    add_dependencies(pkg-tgz pkg-tgz-${PKG_NAME})
endmacro (add_tgz_package)

# Base packaging target
add_custom_target(pkg)

# Check if we can build DEB files
find_program(DPKG "dpkg")
mark_as_advanced(DPKG)
find_program(DPKGSHLIBDEPS "dpkg-shlibdeps")
mark_as_advanced(DPKGSHLIBDEPS)
if (DPKG)
    if (NOT DPKGSHLIBDEPS)
        MESSAGE(FATAL_ERROR "dpkg-shlibdeps program not found but is required.")
    endif ()
    add_custom_target(pkg-deb)
    add_dependencies(pkg pkg-deb)
endif (DPKG)

# Check if we can build RPM files
find_program(RPMBUILD "rpmbuild")
mark_as_advanced(RPMBUILD)
if (RPMBUILD)
    add_custom_target(pkg-rpm)
    add_dependencies(pkg pkg-rpm)
endif (RPMBUILD)

# Binary archive target
add_custom_target(pkg-tgz)
add_dependencies(pkg pkg-tgz)
