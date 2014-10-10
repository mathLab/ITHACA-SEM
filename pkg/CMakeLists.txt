include (CMakeCommon.cmake)

# Package content lists
list(APPEND libnektar++-utilities_LIBS LibUtilities)
list(APPEND libnektar++-stdregions_LIBS StdRegions)
list(APPEND libnektar++-spatialdomains_LIBS SpatialDomains)
list(APPEND libnektar++-localregions_LIBS LocalRegions)
list(APPEND libnektar++-multiregions_LIBS MultiRegions)
list(APPEND libnektar++-solverutils_LIBS SolverUtils)
list(APPEND nektar++-incnavierstokes-solver_BINS IncNavierStokesSolver)
list(APPEND nektar++-cardiacep-solver_BINS CardiacEPSolver)
list(APPEND nektar++-compressibleflow-solver_BINS CompressibleFlowSolver)
list(APPEND nektar++-adr-solver_BINS ADRSolver)
list(APPEND nektar++-diffusion-solver_BINS DiffusionSolver 
                                           DiffusionSolverTimeInt)
list(APPEND nektar++-shallowwater-solver_BINS ShallowWaterSolver)
list(APPEND nektar++-utilities_BINS MeshConvert FieldConvert FldToVtk
    FldToTecplot XmlToVtk XmlToTecplot ProbeFld)
list(APPEND nektar++-demos_BINS Helmholtz1D Helmholtz2D Helmholtz3D
    Helmholtz3DHomo1D Helmholtz3DHomo2D
    StdProject0D StdProject1D StdProject2D StdProject3D 
    StdProject_Diff1D StdProject_Diff2D StdProject_Diff3D
    SteadyAdvectionDiffusionReaction2D SteadyLinearAdvectionReaction2D
    LocProject1D LocProject2D LocProject3D LocProject_Diff1D
    LocProject_Diff2D LocProject_Diff3D
    HDGHelmholtz1D HDGHelmholtz2D HDGHelmholtz3D HDGHelmholtz3DHomo1D)

# Nektar++ binary archive
list(APPEND nektar++_LIBS LibUtilities StdRegions SpatialDomains LocalRegions
                          MultiRegions SolverUtils)
list(APPEND nektar++_BINS ADRSolver IncNavierStokesSolver CardiacEPSolver
                          CompressibleFlowSolver DiffusionSolver 
                          DiffusionSolverTimeInt ShallowWaterSolver
                          ${nektar++-utilities_BINS}
                          ${nektar++-demos_BINS})

# add_deb_package macro used to create debian packaging project.
# Note that the formatting of the DESCRIPTION field is VERY specific. It
# must start with a new line, each line must start with a space, and there
# must be NO blank line at the end of the description.
if (NEKTAR_BUILD_LIBRARY)
    ### DEBIAN PACKAGES ############################################
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
    set(NEK_DEP "${NEK_DEP}, libnektar++-stdregions (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME libnektar++-spatialdomains
        SUMMARY "Nektar++ SpatialDomains library"
        DESCRIPTION "
 This library provides the mappings between reference regions and physical
 regions in the domain."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-spatialdomains_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-spatialdomains (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME libnektar++-localregions
        SUMMARY "Nektar++ LocalRegions library"
        DESCRIPTION "
 This library provides physical space expansions on the various supported
 regions."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-localregions_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-localregions (= ${NEKTAR_VERSION})")

    add_deb_package(
        NAME libnektar++-multiregions
        SUMMARY "Nektar++ MultiRegions library"
        DESCRIPTION "
 This library provides global expansions on multi-element domains."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-multiregions_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-multiregions (= ${NEKTAR_VERSION})")

    if (NEKTAR_BUILD_DEMOS)
        add_deb_package(
            NAME nektar++-demos
            SUMMARY "Nektar++ demonstration applications"
            DESCRIPTION "
 This package contains a number of example programs for solving simple
 problems using the spectral/hp element method."
            DEPENDS ${NEK_DEP}
            INSTALL_BINS "${nektar++-demos_BINS}")
    endif (NEKTAR_BUILD_DEMOS)

    add_deb_package(
        NAME libnektar++-solverutils
        SUMMARY "Nektar++ SolverUtils library"
        DESCRIPTION "
 This library provides support classes and routines for constructing complete
 spectral/hp element solvers."
        DEPENDS ${NEK_DEP}
        INSTALL_LIBS "${libnektar++-solverutils_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-solverutils (= ${NEKTAR_VERSION})")

    if (NEKTAR_SOLVER_INCNAVIERSTOKES)
        add_deb_package(
            NAME nektar++-incnavierstokes-solver
            SUMMARY "Nektar++ incompressible Navier-Stokes flow solver"
            DESCRIPTION "
 Solves the incompressible Navier-Stokes equations using the Nektar++ high-
 order spectral/hp element framework"
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-incnavierstokes-solver_BINS}")
    endif (NEKTAR_SOLVER_INCNAVIERSTOKES)

    if (NEKTAR_SOLVER_ADR)
        add_deb_package(
            NAME nektar++-adr-solver
            SUMMARY "Nektar++ advection-diffusion-reaction solver"
            DESCRIPTION "
 Solves a range of advection, diffusion, reaction PDE problems."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-adr-solver_BINS}")
    endif (NEKTAR_SOLVER_ADR)

    if (NEKTAR_SOLVER_CARDIAC_EP)
        add_deb_package(
            NAME nektar++-cardiacep-solver
            SUMMARY "Nektar++ cardiac electrophysiology solver"
            DESCRIPTION "
 Solves the cardiac electrophysiology monodomain equations."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-cardiacep-solver_BINS}")
    endif (NEKTAR_SOLVER_CARDIAC_EP)

    if (NEKTAR_SOLVER_SHALLOW_WATER)
        add_deb_package(
            NAME nektar++-shallowwater-solver
            SUMMARY "Nektar++ shallow water equations solver"
            DESCRIPTION "
 Solves the shallow water equations."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-shallowwater-solver_BINS}")
    endif (NEKTAR_SOLVER_SHALLOW_WATER)

    if (NEKTAR_SOLVER_COMPRESSIBLE_FLOW)
        add_deb_package(
            NAME nektar++-compressibleflow-solver
            SUMMARY "Nektar++ compressible Navier-Stokes flow solver"
            DESCRIPTION "
 Solves the compressible Euler and compressible Navier-Stokes equations."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-compressibleflow-solver_BINS}")
    endif (NEKTAR_SOLVER_COMPRESSIBLE_FLOW)

    if (NEKTAR_BUILD_UTILITIES)
        add_deb_package(
            NAME nektar++-utilities
            SUMMARY "Nektar++ pre- and post-processing utilities"
            DESCRIPTION "
 This package contains the mesh conversion utility and other utilities for
 converting output files into common visualisation formats"
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-utilities_BINS}")
    endif (NEKTAR_BUILD_UTILITIES)

    # Top level meta package
    add_deb_package(
        NAME nektar++
        SUMMARY "Nektar++ complete meta package"
        DESCRIPTION "
 This meta-package depends on the entire Nektar++ suite."
        DEPENDS "${NEK_DEP}, nektar++-incnavierstokes-solver (= ${NEKTAR_VERSION}), nektar++-adr-solver (=${NEKTAR_VERSION}), nektar++-cardiacep-solver (=${NEKTAR_VERSION}), nektar++-shallowwater-solver (=${NEKTAR_VERSION}), nektar++-compressibleflow-solver (=${NEKTAR_VERSION}), nektar++-utilities (=${NEKTAR_VERSION})")

    ### RPM PACKAGES ###############################################
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

    if (NEKTAR_BUILD_DEMOS)
        add_rpm_package(
            NAME nektar++-demos
            SUMMARY "Nektar++ demonstration applications"
            DESCRIPTION 
            "This package contains a number of example programs for solving
             simple problems using the spectral/hp element method."
            INSTALL_BINS "${nektar++-demos_BINS}")
    endif (NEKTAR_BUILD_DEMOS)

    if (NEKTAR_SOLVER_INCNAVIERSTOKES)
        add_rpm_package(
            NAME nektar++-incnavierstokes-solver
            SUMMARY "Nektar++ incompressible Navier-Stokes flow solver"
            DESCRIPTION
            "Solves the incompressible Navier-Stokes equations using the    
             Nektar++ high-order spectral/hp element framework"
            INSTALL_BINS "${nektar++-incnavierstokes-solver_BINS}")
    endif (NEKTAR_SOLVER_INCNAVIERSTOKES)

    if (NEKTAR_SOLVER_ADR)
        add_rpm_package(
            NAME nektar++-adr-solver
            SUMMARY "Nektar++ advection-diffusion-reaction solver"
            DESCRIPTION
            "Solves a range of advection, diffusion, reaction PDE problems."
            INSTALL_BINS "${nektar++-adr-solver_BINS}")
    endif (NEKTAR_SOLVER_ADR)

    if (NEKTAR_SOLVER_CARDIAC_EP)
        add_rpm_package(
            NAME nektar++-cardiacep-solver
            SUMMARY "Nektar++ cardiac electrophysiology solver"
            DESCRIPTION
            "Solves the cardiac electrophysiology monodomain equations."
            INSTALL_BINS "${nektar++-cardiacep-solver_BINS}")
    endif (NEKTAR_SOLVER_CARDIAC_EP)

    if (NEKTAR_SOLVER_SHALLOW_WATER)
        add_rpm_package(
            NAME nektar++-shallowwater-solver
            SUMMARY "Nektar++ shallow water equations solver"
            DESCRIPTION
            "Solves the shallow water equations."
            INSTALL_BINS "${nektar++-shallowwater-solver_BINS}")
    endif (NEKTAR_SOLVER_SHALLOW_WATER)

    if (NEKTAR_SOLVER_COMPRESSIBLE_FLOW)
        add_rpm_package(
            NAME nektar++-compressibleflow-solver
            SUMMARY "Nektar++ compressible Navier-Stokes flow solver"
            DESCRIPTION
            "Solves the compressible Euler and compressible Navier-Stokes 
             equations."
            INSTALL_BINS "${nektar++-compressibleflow-solver_BINS}")
    endif (NEKTAR_SOLVER_COMPRESSIBLE_FLOW)

    if (NEKTAR_BUILD_UTILITIES)    
        add_rpm_package(
            NAME nektar++-utilities
            SUMMARY "Nektar++ pre- and post-processing utilities"
            DESCRIPTION "This package contains the mesh conversion utility and
            other utilities for converting output files into common
            visualisation formats"
            INSTALL_BINS "${nektar++-utilities_BINS}")
    endif (NEKTAR_BUILD_UTILITIES)

    add_rpm_package(
        NAME nektar++
        SUMMARY "Nektar++ complete meta package"
        DESCRIPTION "
 This meta-package depends on the entire Nektar++ suite."
        DEPENDS "nektar++-incnavierstokes-solver = ${NEKTAR_VERSION}, nektar++-adr-solver = ${NEKTAR_VERSION}, nektar++-cardiacep-solver = ${NEKTAR_VERSION}, nektar++-shallowwater-solver = ${NEKTAR_VERSION}, nektar++-compressibleflow-solver = ${NEKTAR_VERSION}, nektar++-utilities = ${NEKTAR_VERSION}")

    add_tgz_package(
        NAME nektar++
        SUMMARY "Nektar++ libraries and solvers"
        DESCRIPTION "This is all of Nektar++"
        INSTALL_LIBS "${nektar++_LIBS}"
    )   
endif (NEKTAR_BUILD_LIBRARY)

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
    "/library/Demos/MultiRegions/ExtraDemos/"
    "/library/Demos/MultiRegions/XmlFiles/"
    "/library/Demos/LocalRegions/XmlFiles/"
    "/solvers/ImageWarpingSolver/"
    "/solvers/VortexWaveInteraction/"
    "/utilities/PreProcessing/Extras/"
    "/utilities/PostProcessing/Extras/"
)


include(CPack)
