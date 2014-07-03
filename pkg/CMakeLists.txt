include (CMakeCommon.cmake)

add_custom_target(pkg)

# LibUtilities
get_target_property(TARGET_LOCATION LibUtilities LOCATION)
list(APPEND libnektar++-utilities_LIBS ${TARGET_LOCATION})
list(APPEND nektar++-LIBS ${TARGET_LOCATION})

# Check if we can build debian files
find_program(DPKG "dpkg")
if (DPKG)
    add_custom_target(pkg-deb)
    add_dependencies(pkg pkg-deb)

    add_deb_package(
        NAME libnektar++-utilities
        SUMMARY "Nektar++ library utilities"
        DESCRIPTION "This is the LibUtilities library."
        INSTALL_LIBS "${libnektar++-utilities_LIBS}")

endif (DPKG)


find_program(RPMBUILD "rpmbuild")
if (RPMBUILD)
    add_custom_target(pkg-rpm)
    add_dependencies(pkg pkg-rpm)

    add_rpm_package(
        NAME libnektar++-utilities
        SUMMARY "Nektar++ library utilities"
        DESCRIPTION "THis is the LibUtilities library."
        INSTALL_LIBS "${libnektar++-utilities_LIBS}")

endif (RPMBUILD)
