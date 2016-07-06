function (write_lib_files PKG_INSTALL_LIBS OUTPUT_FILE)
    # Find library file and add the versioned form of each library
    set(PKG_INSTALL_LIBS_FILES)
    foreach(l ${PKG_INSTALL_LIBS})
        IF(${CMAKE_MAJOR_VERSION} LESS 3)
            get_target_property(TARGET_LOCATION ${l} LOCATION)
        ELSE ()
            SET(TARGET_LOCATION $<TARGET_LINKER_FILE:${l}>)
        ENDIF()
        if (NOT TARGET_LOCATION)
            message(FATAL_ERROR "Target '${l}' could not be found.")
        endif ()
        list(APPEND PKG_INSTALL_LIBS_FILES ${TARGET_LOCATION})
        if (APPLE)
            list(APPEND PKG_INSTALL_LIBS_FILES $<TARGET_FILE:${l}>)
        else ()
            list(APPEND PKG_INSTALL_LIBS_FILES 
                        ${TARGET_LOCATION}.${NEKTAR_VERSION})
        endif()
    endforeach()

    # Output the list of files to be installed in the package
    IF(${CMAKE_MAJOR_VERSION} LESS 3)
        file(WRITE "${OUTPUT_FILE}" "${PKG_INSTALL_LIBS_FILES}")
    ELSE ()
        file(GENERATE OUTPUT "${OUTPUT_FILE}"
             CONTENT "${PKG_INSTALL_LIBS_FILES}")
    ENDIF ()
endfunction ()

function (write_bin_files PKG_INSTALL_BINS OUTPUT_FILE)
    # Find binary files
    set(PKG_INSTALL_BINS_FILES)
    foreach(b ${PKG_INSTALL_BINS})
        IF(EXISTS ${b})
            get_filename_component(TARGET_LOCATION ${b} REALPATH)
        ELSEIF(${CMAKE_MAJOR_VERSION} LESS 3)
            get_target_property(TARGET_LOCATION ${b} LOCATION)
        ELSE ()
            SET(TARGET_LOCATION $<TARGET_FILE:${b}>)
        ENDIF ()
        if (NOT TARGET_LOCATION)
            message(FATAL_ERROR "Target '${b}' could not be found.")
        endif ()
        list(APPEND PKG_INSTALL_BINS_FILES ${TARGET_LOCATION})
    endforeach()

    # Output the list of files to be installed in the package
    IF(${CMAKE_MAJOR_VERSION} LESS 3)
        file(WRITE "${OUTPUT_FILE}" "${PKG_INSTALL_BINS_FILES}")
    ELSE ()
        file(GENERATE OUTPUT "${OUTPUT_FILE}"
             CONTENT "${PKG_INSTALL_BINS_FILES}")
    ENDIF ()
endfunction ()

macro (add_nektar_package)
    set(options "")
    set(oneValueArgs NAME SUMMARY DESCRIPTION TYPE)
    set(multiValueArgs INSTALL_LIBS INSTALL_BINS BREAKS CONFLICTS DEPENDS)
    cmake_parse_arguments(PKG "${options}"
        "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PKG_NAME}-${PKG_TYPE})

    file(MAKE_DIRECTORY "${BUILD_DIR}/targets")
    write_lib_files("${PKG_INSTALL_LIBS}"
                    "${BUILD_DIR}/targets/install_libs.txt")
    write_bin_files("${PKG_INSTALL_BINS}"
                    "${BUILD_DIR}/targets/install_bins.txt")

    if(PKG_TYPE STREQUAL "deb")
        set(PKG_GENERATOR "DEB")
    elseif(PKG_TYPE STREQUAL "rpm")
        set(PKG_GENERATOR "RPM")
    elseif(PKG_TYPE STREQUAL "pkgmaker")
        set(PKG_GENERATOR "PackageMaker")
    elseif(PKG_TYPE STREQUAL "tgz")
        set(PKG_GENERATOR "TGZ")
    else()
        message(ERROR "Unknown package type: ${PKG_TYPE}")
    endif()

    # Configure project for this package
    configure_file(NektarPackage.cmake.in ${BUILD_DIR}/CMakeLists.txt @ONLY)

    add_custom_target(
        pkg-${PKG_TYPE}-${PKG_NAME}
        rm -f ${BUILD_DIR}/CPackConfig.cmake
        COMMAND ${CMAKE_COMMAND} .
        COMMAND ${CMAKE_CPACK_COMMAND} --config CPackConfig.cmake -G ${PKG_GENERATOR}
        WORKING_DIRECTORY ${BUILD_DIR}
        )
    if (PKG_INSTALL_LIBS OR PKG_INSTALL_BINS)
        add_dependencies(pkg-${PKG_TYPE}-${PKG_NAME}
            ${PKG_INSTALL_LIBS} ${PKG_INSTALL_BINS})
    endif ()
    add_dependencies(pkg-${PKG_TYPE} pkg-${PKG_TYPE}-${PKG_NAME})
endmacro()

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

# Check if we can build PackageMaker .pkg files. PackageMaker is technically
# deprecated but still works at least up to OS X 10.11. It could therefore live
# in a bunch of locations, so check common places.
find_program(PKGMAKER "PackageMaker"
    PATHS /Applications/Xcode.app/Contents/Applications/PackageMaker.app/Contents/MacOS
          /Applications/Utilities/PackageMaker.app/Contents/MacOS
          /Applications/PackageMaker.app/Contents/MacOS
          /Developer/Applications/Utilities/PackageMaker.app/Contents/MacOS
          /Developer/Applications/PackageMaker.app/Contents/MacOS)
if (PKGMAKER)
    add_custom_target(pkg-pkgmaker)
    add_dependencies(pkg pkg-pkgmaker)
endif()

# Binary archive target
add_custom_target(pkg-tgz)
add_dependencies(pkg pkg-tgz)
