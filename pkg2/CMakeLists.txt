set(CPACK_PACKAGE_VENDOR "Imperial College London")
set(CPACK_PACKAGE_CONTACT 
    "Nektar++ users mailing list <nektar-users@imperial.ac.uk>")

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
set(CPACK_PACKAGE_VERSION "${NEKTAR_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${NEKTAR_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${NEKTAR_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${NEKTAR_VERSION_PATCH}")

# Debian-specific options
set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_DEBIAN_PACKAGE_DEBUG ON)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER 
    "Chris Cantwell <c.cantwell@imperial.ac.uk>")
set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://www.nektar.info")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

# RPM-specific options
set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
set(CPACK_RPM_PACKAGE_LICENSE "MIT")
set(CPACK_RPM_PACKAGE_DEBUG 0)

# Set up generator-specific logic (i.e. that may change the above depending on
# generator) using additional CPack configuration file.
configure_file(${CMAKE_SOURCE_DIR}/cmake/NektarCPackConfig.cmake.in
               ${CMAKE_BINARY_DIR}/NektarCPackConfig.cmake @ONLY)
set(CPACK_PROJECT_CONFIG_FILE ${CMAKE_BINARY_DIR}/NektarCPackConfig.cmake)

# This lovely block of code is designed to patch up our executables ready for OS
# X deployment. The strategy is to create a piece of INSTALL code that will
# patch up each executable and copy any libraries that we depend on, so that
# things like ARPACK can be compiled using MacPorts/Homebrew and we don't have
# to bother messing around compiling them ourselves.
#
# The strategy is:
#
# - Examine the library dependencies using the get_prerequisites function from
#   the GetPrerequisites CMake module (that essentially runs otool recursively)
# - Copy each of the libraries into dist/lib
# - Run install_name_tool on each library to remove the absolute path to the
#   installed library.
# - This is then replaced by @rpath/../lib/libName.dylib so that @rpath can be
#   embedded into the executable.
# - Change the library install name so that it follows the same convention
# - Extract the rpaths from our executable using the get_item_rpaths function
#   from BundleUtilities
# - Finally, set the rpath of each executable to be @executable_path so that
#   everything becomes relative to the executable
#
# All of this makes an entirely self-contained version of Nektar++ that can be
# extracted anywhere on the system and still retain its dependencies.

# Finally, include the CPack module
include(CPack)

if(APPLE)
    foreach (component ${CPACK_COMPONENTS_ALL})
        install(CODE "
    include(GetPrerequisites)
    include(BundleUtilities)
    file(GLOB apps \"\${CMAKE_INSTALL_PREFIX}/bin/*\")
    message(STATUS \"APPLIST: \${CMAKE_INSTALL_PREFIX}/bin/*\")
    foreach (app \${apps})
        get_filename_component(exepath \${app} DIRECTORY)
        get_prerequisites(\"\${app}\" prereqs 1 1 \"\${exepath}\" \"\" \"\" 1)

        message(STATUS \"Fixing executable: \${app}\")

        # First pass -- copy libraries and figure out -change command for
        # install_name_tool. Note that since install_name_tool doesn't complain
        # if the thing you're changing isn't there, we just throw everything
        # at every library.
        unset(changes)
        unset(changes_lib)

        # Check if thirdparty component has been set up yet
        foreach(req \${prereqs})
            message(STATUS \"asdf: \${req}\")
            get_filename_component(reqname \${req} NAME)
            set(libdest \${CMAKE_INSTALL_PREFIX}/lib/nektar++-${NEKTAR_VERSION}/\${reqname})
            set(changes \${changes} \"-change\" \"\${req}\" \"@executable_path/../lib/nektar++-${NEKTAR_VERSION}/\${reqname}\")
            set(changes_lib \${changes_lib} \"-change\" \"\${req}\" \"@loader_path/../lib/nektar++-${NEKTAR_VERSION}/\${reqname}\")

            # Copy this library
            if (NOT EXISTS \${libdest})
                file(COPY \${req}
                     DESTINATION \${CMAKE_INSTALL_PREFIX}/lib/nektar++-${NEKTAR_VERSION}
                     FILE_PERMISSIONS OWNER_WRITE OWNER_READ)

                # If the library was symlinked, we follow the symlink and then 
                # copy that, too.
                if (IS_SYMLINK \${req})
                    # resolve symlink
                    get_filename_component(req_abs \${req} REALPATH)
                    file(COPY \${req_abs}
                         DESTINATION \${CMAKE_INSTALL_PREFIX}/lib/nektar++-${NEKTAR_VERSION}
                         FILE_PERMISSIONS OWNER_WRITE OWNER_READ)
                endif()
            endif()
        endforeach()

        # Second pass -- fix up library to use @loader_path/../lib/nektar++-${NEKTAR_VERSION}/libName.dylib
        foreach(req \${prereqs})
            set(searchname \${req})

            if (IS_SYMLINK \${req})
                get_filename_component(req_abs \${req} REALPATH)
                set(req \${req_abs})
            endif()

            get_filename_component(reqname \${req} NAME)
            set(libdest \${CMAKE_INSTALL_PREFIX}/lib/nektar++-${NEKTAR_VERSION}/\${reqname})
            set(cmd install_name_tool \${changes_lib} \"\${libdest}\")
            execute_process(COMMAND \${cmd} RESULT_VARIABLE install_name_tool_result)

            # change library install name ID
            execute_process(COMMAND install_name_tool -id @loader_path/../lib/nektar++-${NEKTAR_VERSION}/\${reqname} \${libdest})

            # change @loader_path (used by some boost libs through homebrew)
        endforeach()

        # Third pass -- fix executable library names
        foreach(req \${prereqs})
            set(cmd install_name_tool \${changes} \"\${app}\")
            execute_process(COMMAND \${cmd} RESULT_VARIABLE install_name_tool_result)
        endforeach()

        # Fix rpath -- unset variables that have been used on previous passes
        unset(app_rpaths)
        get_item_rpaths(\"\${app}\" app_rpaths)
        unset(app_rpaths_unique)
        unset(rpath_changes)
        foreach(rpath \${app_rpaths})
            gp_append_unique (app_rpaths_unique \${rpath})
        endforeach()
        foreach(rpath \${app_rpaths_unique})
            set(rpath_changes \${rpath_changes} -delete_rpath \${rpath})
        endforeach()
        set(rpath_changes \${rpath_changes} -add_rpath @executable_path)
        execute_process(COMMAND install_name_tool \${rpath_changes} \"\${app}\")
    endforeach()" COMPONENT ${component})
    endforeach()
endif()
