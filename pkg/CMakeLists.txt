include (CMakeCommon.cmake)

add_custom_target(pkg)

# LibUtilities
get_target_property(TARGET_LOCATION LibUtilities LOCATION)
list(APPEND libnektar++-utilities_LIBS ${TARGET_LOCATION})
list(APPEND nektar++_LIBS ${TARGET_LOCATION})

# StdRegions
get_target_property(TARGET_LOCATION StdRegions LOCATION)
list(APPEND libnektar++-stdregions_LIBS ${TARGET_LOCATION})
list(APPEND nektar++_LIBS ${TARGET_LOCATION})

# SpatialDomains
get_target_property(TARGET_LOCATION SpatialDomains LOCATION)
list(APPEND libnektar++-spatialdomains_LIBS ${TARGET_LOCATION})
list(APPEND nektar++_LIBS ${TARGET_LOCATION})

# LocalRegions
get_target_property(TARGET_LOCATION LocalRegions LOCATION)
list(APPEND libnektar++-localregions_LIBS ${TARGET_LOCATION})
list(APPEND nektar++_LIBS ${TARGET_LOCATION})

# MultiRegions
get_target_property(TARGET_LOCATION MultiRegions LOCATION)
list(APPEND libnektar++-multiregions_LIBS ${TARGET_LOCATION})
list(APPEND nektar++_LIBS ${TARGET_LOCATION})

# SolverUtils
get_target_property(TARGET_LOCATION SolverUtils LOCATION)
list(APPEND libnektar++-solverutils_LIBS ${TARGET_LOCATION})
list(APPEND nektar++_LIBS ${TARGET_LOCATION})

# Check if we can build debian files
find_program(DPKG "dpkg")
if (DPKG)
    add_custom_target(pkg-deb)
    add_dependencies(pkg pkg-deb)

    add_deb_package(
        NAME libnektar++-utilities
        SUMMARY "Nektar++ library utilities"
        DESCRIPTION 
        "This library provides core routines including linear algebra and
         integration with ThirdParty libraries."
        INSTALL_LIBS "${libnektar++-utilities_LIBS}")
    add_deb_package(
        NAME libnektar++-stdregions
        SUMMARY "Nektar++ StdRegions library"
        DESCRIPTION 
        "This library provides construction of the reference expansions for the
         various 1D, 2D and 3D regions."
        DEPENDS "libnektar++-utilities (= ${NEKTAR_VERSION})"
        INSTALL_LIBS "${libnektar++-stdregions_LIBS}")
    add_deb_package(
        NAME libnektar++-spatialdomains
        SUMMARY "Nektar++ SpatialDomains library"
        DESCRIPTION
        "This library provides the mappings between reference regions and
         physical regions in the domain."
        DEPENDS "libnektar++-stdregions (= ${NEKTAR_VERSION}), libnektar++-utilities (= ${NEKTAR_VERSION})"
        INSTALL_LIBS "${libnektar++-spatialdomains_LIBS}")
    add_deb_package(
        NAME libnektar++-localregions
        SUMMARY "Nektar++ LocalRegions library"
        DESCRIPTION 
        "This library provides physical space expansions on the various
         supported regions."
        DEPENDS "libnektar++-stdregions (= ${NEKTAR_VERSION}), libnektar++-utilities (= ${NEKTAR_VERSION}), libnektar++-spatialdomains (= ${NEKTAR_VERSION})"
        INSTALL_LIBS "${libnektar++-localregions_LIBS}")
    add_deb_package(
        NAME libnektar++-multiregions
        SUMMARY "Nektar++ MultiRegions library"
        DESCRIPTION
        "This library provides global expansions on multi-element domains."
        DEPENDS "libnektar++-stdregions (= ${NEKTAR_VERSION}), libnektar++-utilities (= ${NEKTAR_VERSION}), libnektar++-spatialdomains (= ${NEKTAR_VERSION}), libnektar++-localregions (= ${NEKTAR_VERSION})"
        INSTALL_LIBS "${libnektar++-multiregions_LIBS}")
    add_deb_package(
        NAME libnektar++-solverutils
        SUMMARY "Nektar++ SolverUtils library"
        DESCRIPTION 
        "This library provides support classes and routines for constructing
         complete spectral/hp element solvers."
        DEPENDS "libnektar++-stdregions (= ${NEKTAR_VERSION}), libnektar++-utilities (= ${NEKTAR_VERSION}), libnektar++-spatialdomains (= ${NEKTAR_VERSION}), libnektar++-localregions (= ${NEKTAR_VERSION}), libnektar++-multiregions (= ${NEKTAR_VERSION})"
        INSTALL_LIBS "${libnektar++-solverutils_LIBS}")

endif (DPKG)


find_program(RPMBUILD "rpmbuild")
if (RPMBUILD)
    add_custom_target(pkg-rpm)
    add_dependencies(pkg pkg-rpm)

    add_rpm_package(
        NAME libnektar++-utilities
        SUMMARY "Nektar++ library utilities"
        DESCRIPTION 
        "This library provides core routines including linear algebra and
         integration with ThirdParty libraries."
        INSTALL_LIBS "${libnektar++-utilities_LIBS}")
    add_rpm_package(
        NAME libnektar++-stdregions
        SUMMARY "Nektar++ StdRegions library"
        DESCRIPTION
        "This library provides construction of the reference expansions for the
         various 1D, 2D and 3D regions."
        DEPENDS "libnektar++-utilities (=${NEKTAR_VERSION})"
        INSTALL_LIBS "${libnektar++-stdregions_LIBS}")
    add_rpm_package(
        NAME libnektar++-spatialdomains
        SUMMARY "Nektar++ SpatialDomains library"
        DESCRIPTION 
        "This library provides the mappings between reference regions and
         physical regions in the domain."
        INSTALL_LIBS "${libnektar++-spatialdomains_LIBS}")
    add_rpm_package(
        NAME libnektar++-localregions
        SUMMARY "Nektar++ LocalRegions library"
        DESCRIPTION
        "This library provides physical space expansions on the various
         supported regions."
        INSTALL_LIBS "${libnektar++-localregions_LIBS}")
    add_rpm_package(
        NAME libnektar++-multiregions
        SUMMARY "Nektar++ MultiRegions library"
        DESCRIPTION
        "This library provides global expansions on multi-element domains."
        INSTALL_LIBS "${libnektar++-multiregions_LIBS}")
    add_rpm_package(
        NAME libnektar++-solverutils
        SUMMARY "Nektar++ SolverUtils library"
        DESCRIPTION 
        "This library provides support classes and routines for constructing
         complete spectral/hp element solvers."
        INSTALL_LIBS "${libnektar++-solverutils_LIBS}")


endif (RPMBUILD)

add_custom_target(pkg-tgz)
add_dependencies(pkg pkg-tgz)

add_tgz_package(
    NAME nektar++
    SUMMARY "Nektar++ libraries and solvers"
    DESCRIPTION "This is all of Nektar++"
    INSTALL_LIBS "${nektar++_LIBS}")
