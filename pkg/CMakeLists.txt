include (CMakeCommon.cmake)

add_custom_target(pkg)

# LibUtilities
list(APPEND libnektar++-utilities_LIBS LibUtilities)
list(APPEND nektar++_LIBS LibUtilities)

# StdRegions
list(APPEND libnektar++-stdregions_LIBS StdRegions)
list(APPEND nektar++_LIBS StdRegions)

# SpatialDomains
list(APPEND libnektar++-spatialdomains_LIBS SpatialDomains)
list(APPEND nektar++_LIBS SpatialDomains)

# LocalRegions
list(APPEND libnektar++-localregions_LIBS LocalRegions)
list(APPEND nektar++_LIBS LocalRegions)

# MultiRegions
list(APPEND libnektar++-multiregions_LIBS MultiRegions)
list(APPEND nektar++_LIBS MultiRegions)

# SolverUtils
list(APPEND libnektar++-solverutils_LIBS SolverUtils)
list(APPEND nektar++_LIBS SolverUtils)

# IncNavierStokesSolver
list(APPEND nektar++-incnavierstokes-solver_BINS IncNavierStokesSolver)
list(APPEND nektar++_BINS IncNavierStokesSolver)

# CardiacEPSolver
list(APPEND nektar++-cardiacep-solver_BINS CardiacEPSolver)
list(APPEND nektar++-BINS CardiacEPSolver)

# CompressibleFlowSolver
list(APPEND nektar++-compressibleflow-solver_BINS CompressibleFlowSolver)
list(APPEND nektar++-BINS CompressibleFlowSolver)

# Check if we can build debian files
find_program(DPKG "dpkg")
mark_as_advanced(DPKG)
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
    #add_deb_package(
    #    NAME nektar++-incnavierstokes-solver
    #    SUMMARY "Nektar++ incompressible Navier-Stokes flow solver"
    #    DESCRIPTION
    #    "Solves the incompressible Navier-Stokes equations using the Nektar++
    #     high-order spectral/hp element framework"
    #    DEPENDS "libnektar++-stdregions (= ${NEKTAR_VERSION}), libnektar++-     utilities (= ${NEKTAR_VERSION}), libnektar++-spatialdomains (=                  ${NEKTAR_VERSION}), libnektar++-localregions (= ${NEKTAR_VERSION}), libnektar++-multiregions (= ${NEKTAR_VERSION})"
    #    INSTALL_BINS "${nektar++-incnavierstokes-solver_BINS}")
endif (DPKG)


find_program(RPMBUILD "rpmbuild")
mark_as_advanced(RPMBUILD)
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
    INSTALL_LIBS "${nektar++_LIBS}"
)

# Configure source package
add_custom_target(pkg-src
    ${CMAKE_CPACK_COMMAND} --config CPackSourceConfig.cmake
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
add_dependencies(pkg pkg-src)

set(CPACK_SOURCE_TBZ2 "ON")
set(CPACK_SOURCE_TGZ "ON")
set(CPACK_SOURCE_TZ "OFF")
set(CPACK_SOURCE_ZIP "ON")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "nektar++-${NEKTAR_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES
    "/CVS/"
    "/\\\\.svn/"
    "/\\\\.bzr/"
    "/\\\\.hg/"
    "/\\\\.git/"
    "/build/"
    "/builds/"
    "/ThirdParty/"
    "/distribute\\\\.sh"
    "/dist-exclude"
    "/Testing/"
    "/solvers/ImageWarpingSolver/"
    "/solvers/PulseWaveSolver/"
    "/solvers/VortexWaveInteraction/"
    "/utilities/PreProcessing/Extras/"
    "/utilities/PostProcessing/Extras/"
)


include(CPack)
