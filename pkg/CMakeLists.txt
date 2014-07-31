include (CMakeCommon.cmake)

add_custom_target(pkg)

# LibUtilities
list(APPEND libnektar++-utilities_LIBS LibUtilities)

# StdRegions
list(APPEND libnektar++-stdregions_LIBS StdRegions)

# SpatialDomains
list(APPEND libnektar++-spatialdomains_LIBS SpatialDomains)

# LocalRegions
list(APPEND libnektar++-localregions_LIBS LocalRegions)

# MultiRegions
list(APPEND libnektar++-multiregions_LIBS MultiRegions)

# SolverUtils
list(APPEND libnektar++-solverutils_LIBS SolverUtils)

# IncNavierStokesSolver
list(APPEND nektar++-incnavierstokes-solver_BINS IncNavierStokesSolver)

# CardiacEPSolver
list(APPEND nektar++-cardiacep-solver_BINS CardiacEPSolver)

# CompressibleFlowSolver
list(APPEND nektar++-compressibleflow-solver_BINS CompressibleFlowSolver)

# ADRSolver
list(APPEND nektar++-adr-solver_BINS ADRSolver)

# DiffusionSolver
list(APPEND nektar++-diffusion-solver_BINS DiffusionSolver 
                                           DiffusionSolverTimeInt)

# ShallowWaterSolver
list(APPEND nektar++-shallowwater-solver_BINS ShallowWaterSolver)

# Utilities
list(APPEND nektar++-utilities_BINS MeshConvert FieldConvert FldToVtk
    FldToTecplot XmlToVtk XmlToTecplot ProbeFld)

# Nektar++ binary archive
list(APPEND nektar++_LIBS LibUtilities StdRegions SpatialDomains LocalRegions
                          MultiRegions SolverUtils)
list(APPEND nektar++_BINS ADRSolver IncNavierStokesSolver CardiacEPSolver
                          CompressibleFlowSolver DiffusionSolver 
                          DiffusionSolverTimeInt ShallowWaterSolver
                          ${nektar++-utilities_BINS})

# Check if we can build debian files
find_program(DPKG "dpkg")
mark_as_advanced(DPKG)
if (DPKG)
    add_custom_target(pkg-deb)
    add_dependencies(pkg pkg-deb)

    # add_deb_package macro used to create debian packaging project.
    # Note that the formatting of the DESCRIPTION field is VERY specific. It
    # must start with a new line, each line must start with a space, and there
    # must be NO blank line at the end of the description.
    add_deb_package(
        NAME libnektar++-utilities
        SUMMARY "Nektar++ library utilities"
        DESCRIPTION "
 This library provides core routines including linear algebra and integration
 with ThirdParty libraries."
        INSTALL_LIBS "${libnektar++-utilities_LIBS}")
    set(NEK_DEP "libnektar++-utilities (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME libnektar++-stdregions
        SUMMARY "Nektar++ StdRegions library"
        DESCRIPTION "
 This library provides construction of the reference expansions for the
 various 1D, 2D and 3D regions."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-stdregions_LIBS}")
    set(NEK_DEP "${NEK_DEP} libnektar++-stdregions (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME libnektar++-spatialdomains
        SUMMARY "Nektar++ SpatialDomains library"
        DESCRIPTION "
 This library provides the mappings between reference regions and physical
 regions in the domain."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-spatialdomains_LIBS}")
    set(NEK_DEP "${NEK_DEP} libnektar++-spatialdomains (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME libnektar++-localregions
        SUMMARY "Nektar++ LocalRegions library"
        DESCRIPTION "
 This library provides physical space expansions on the various supported
 regions."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-localregions_LIBS}")
    set(NEK_DEP "${NEK_DEP} libnektar++-localregions (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME libnektar++-multiregions
        SUMMARY "Nektar++ MultiRegions library"
        DESCRIPTION "
 This library provides global expansions on multi-element domains."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-multiregions_LIBS}")
    set(NEK_DEP "${NEK_DEP} libnektar++-multiregions (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME libnektar++-solverutils
        SUMMARY "Nektar++ SolverUtils library"
        DESCRIPTION "
 This library provides support classes and routines for constructing complete
 spectral/hp element solvers."
        DEPENDS ${NEK_DEP}
        INSTALL_LIBS "${libnektar++-solverutils_LIBS}")
    set(NEK_DEP "${NEK_DEP} libnektar++-solverutils (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME nektar++-incnavierstokes-solver
        SUMMARY "Nektar++ incompressible Navier-Stokes flow solver"
        DESCRIPTION "
 Solves the incompressible Navier-Stokes equations using the Nektar++ high-
 order spectral/hp element framework"
        DEPENDS "${NEK_DEP}"
        INSTALL_BINS "${nektar++-incnavierstokes-solver_BINS}")
    add_deb_package(
        NAME nektar++-adr-solver
        SUMMARY "Nektar++ advection-diffusion-reaction solver"
        DESCRIPTION
        "Solves a range of advection, diffusion, reaction PDE problems."
        DEPENDS "${NEK_DEP}"
        INSTALL_BINS "${nektar++-adr-solver_BINS}")
    add_deb_package(
        NAME nektar++-cardiacep-solver
        SUMMARY "Nektar++ cardiac electrophysiology solver"
        DESCRIPTION "
 Solves the cardiac electrophysiology monodomain equations."
        DEPENDS "${NEK_DEP}"
        INSTALL_BINS "${nektar++-cardiacep-solver_BINS}")
    add_deb_package(
        NAME nektar++-shallowwater-solver
        SUMMARY "Nektar++ shallow water equations solver"
        DESCRIPTION
        "Solves the shallow water equations."
        DEPENDS "${NEK_DEP}"
        INSTALL_BINS "${nektar++-shallowwater-solver_BINS}")
    add_deb_package(
        NAME nektar++-compressibleflow-solver
        SUMMARY "Nektar++ compressible Navier-Stokes flow solver"
        DESCRIPTION "
 Solves the compressible Euler and compressible Navier-Stokes equations."
        DEPENDS "${NEK_DEP}"
        INSTALL_BINS "${nektar++-compressibleflow-solver_BINS}")

    add_deb_package(
        NAME nektar++-utilities
        SUMMARY "Nektar++ pre- and post-processing utilities"
        DESCRIPTION "
 This package contains the mesh conversion utility and other utilities for
 converting output files into common visualisation formats"
        DEPENDS "${NEK_DEP}"
        INSTALL_BINS "${nektar++-utilities_BINS}")

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

    add_rpm_package(
        NAME nektar++-incnavierstokes-solver
        SUMMARY "Nektar++ incompressible Navier-Stokes flow solver"
        DESCRIPTION
        "Solves the incompressible Navier-Stokes equations using the Nektar++
         high-order spectral/hp element framework"
        INSTALL_BINS "${nektar++-incnavierstokes-solver_BINS}")
    add_rpm_package(
        NAME nektar++-adr-solver
        SUMMARY "Nektar++ advection-diffusion-reaction solver"
        DESCRIPTION
        "Solves a range of advection, diffusion, reaction PDE problems."
        INSTALL_BINS "${nektar++-adr-solver_BINS}")
    add_rpm_package(
        NAME nektar++-cardiacep-solver
        SUMMARY "Nektar++ cardiac electrophysiology solver"
        DESCRIPTION
        "Solves the cardiac electrophysiology monodomain equations."
        INSTALL_BINS "${nektar++-cardiacep-solver_BINS}")
    add_rpm_package(
        NAME nektar++-shallowwater-solver
        SUMMARY "Nektar++ shallow water equations solver"
        DESCRIPTION
        "Solves the shallow water equations."
        INSTALL_BINS "${nektar++-shallowwater-solver_BINS}")
    add_rpm_package(
        NAME nektar++-compressibleflow-solver
        SUMMARY "Nektar++ compressible Navier-Stokes flow solver"
        DESCRIPTION
        "Solves the compressible Euler and compressible Navier-Stokes 
         equations."
        INSTALL_BINS "${nektar++-compressibleflow-solver_BINS}")
    
    add_rpm_package(
        NAME nektar++-utilities
        SUMMARY "Nektar++ pre- and post-processing utilities"
        DESCRIPTION "
 This package contains the mesh conversion utility and other utilities for
 converting output files into common visualisation formats"
         INSTALL_BINS "${nektar++-utilities_BINS}")
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
