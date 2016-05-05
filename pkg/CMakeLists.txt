include (CMakeCommon.cmake)

# Package content lists

# Libraries
list(APPEND libnektar++-utilities_LIBS      LibUtilities)
list(APPEND libnektar++-stdregions_LIBS     StdRegions)
list(APPEND libnektar++-spatialdomains_LIBS SpatialDomains)
list(APPEND libnektar++-localregions_LIBS   LocalRegions)
list(APPEND libnektar++-nekmeshutils_LIBS   NekMeshUtils)
list(APPEND libnektar++-collections_LIBS    Collections)
list(APPEND libnektar++-multiregions_LIBS   MultiRegions)
list(APPEND libnektar++-globalmapping_LIBS  GlobalMapping)
list(APPEND libnektar++-solverutils_LIBS    SolverUtils)

# Solvers
list(APPEND nektar++-incnavierstokes-solver_BINS  IncNavierStokesSolver)
list(APPEND nektar++-cardiacep-solver_BINS        CardiacEPSolver)
list(APPEND nektar++-compressibleflow-solver_BINS CompressibleFlowSolver)
list(APPEND nektar++-adr-solver_BINS              ADRSolver)
list(APPEND nektar++-diffusion-solver_BINS        DiffusionSolver 
                                                  DiffusionSolverTimeInt)
list(APPEND nektar++-shallowwater-solver_BINS     ShallowWaterSolver)
list(APPEND nektar++-pulsewave-solver_BINS        PulseWaveSolver)

# Utilities
list(APPEND nektar++-utilities_BINS NekMesh FieldConvert)

# Demos
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
                          NekMeshUtils Collections MultiRegions GlobalMapping
                          SolverUtils)
list(APPEND nektar++_BINS ADRSolver IncNavierStokesSolver CardiacEPSolver
                          CompressibleFlowSolver DiffusionSolver
                          DiffusionSolverTimeInt ShallowWaterSolver
                          PulseWaveSolver)

if (NEKTAR_BUILD_DEMOS)
    list(APPEND nektar++_BINS ${nektar++-demos_BINS})
endif()

if (NEKTAR_BUILD_UTILITIES)
    list(APPEND nektar++_BINS ${nektar++-utilities_BINS})
endif()

if (NEKTAR_BUILD_LIBRARY AND DPKG)
    ### DEBIAN PACKAGES ############################################
    # Note that the formatting of the DESCRIPTION field is VERY specific. It
    # must start with a new line, each line must start with a space, and there
    # must be NO blank line at the end of the description.
    add_nektar_package(
        NAME libnektar++-utilities
        TYPE deb
        SUMMARY "Nektar++ library utilities"
        DESCRIPTION "
 This library provides core routines including linear algebra and integration
 with ThirdParty libraries."
        INSTALL_LIBS "${libnektar++-utilities_LIBS}")
    set(NEK_DEP "libnektar++-utilities (= ${NEKTAR_VERSION})")

    add_nektar_package(
        NAME libnektar++-stdregions
        TYPE deb
        SUMMARY "Nektar++ StdRegions library"
        DESCRIPTION "
 This library provides construction of the reference expansions for the
 various 1D, 2D and 3D regions."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-stdregions_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-stdregions (= ${NEKTAR_VERSION})")

    add_nektar_package(
        NAME libnektar++-spatialdomains
        TYPE deb
        SUMMARY "Nektar++ SpatialDomains library"
        DESCRIPTION "
 This library provides the mappings between reference regions and physical
 regions in the domain."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-spatialdomains_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-spatialdomains (= ${NEKTAR_VERSION})")

    add_nektar_package(
        NAME libnektar++-localregions
        TYPE deb
        SUMMARY "Nektar++ LocalRegions library"
        DESCRIPTION "
 This library provides physical space expansions on the various supported
 regions."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-localregions_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-localregions (= ${NEKTAR_VERSION})")

    add_nektar_package(
        NAME libnektar++-nekmeshutils
        TYPE deb
        SUMMARY "Nektar++ NekMeshUtils library"
        DESCRIPTION "
 This library provides mesh generation and routines for the NekMesh
 executable."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-nekmeshutils_LIBS}")

    add_nektar_package(
        NAME libnektar++-collections
        TYPE deb
        SUMMARY "Nektar++ Collections library"
        DESCRIPTION "
 This library provides efficient multi-expansion implementations of operators
 on the various supported regions."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-collections_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-collections (= ${NEKTAR_VERSION})")

    add_nektar_package(
        NAME libnektar++-multiregions
        TYPE deb
        SUMMARY "Nektar++ MultiRegions library"
        DESCRIPTION "
 This library provides global expansions on multi-element domains."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-multiregions_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-multiregions (= ${NEKTAR_VERSION})")

    add_nektar_package(
        NAME libnektar++-globalmapping
        TYPE deb
        SUMMARY "Nektar++ GlobalMapping library"
        DESCRIPTION "
 This library provides deformation routines for Nektar++ simulations."
        DEPENDS "${NEK_DEP}"
        INSTALL_LIBS "${libnektar++-globalmapping_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-globalmapping (= ${NEKTAR_VERSION})")

    if (NEKTAR_BUILD_DEMOS)
        add_nektar_package(
            NAME nektar++-demos
            TYPE deb
            SUMMARY "Nektar++ demonstration applications"
            DESCRIPTION "
 This package contains a number of example programs for solving simple
 problems using the spectral/hp element method."
            DEPENDS ${NEK_DEP}
            INSTALL_BINS "${nektar++-demos_BINS}")
    endif (NEKTAR_BUILD_DEMOS)

    add_nektar_package(
        NAME libnektar++-solverutils
        TYPE deb
        SUMMARY "Nektar++ SolverUtils library"
        DESCRIPTION "
 This library provides support classes and routines for constructing complete
 spectral/hp element solvers."
        DEPENDS ${NEK_DEP}
        INSTALL_LIBS "${libnektar++-solverutils_LIBS}")
    set(NEK_DEP "${NEK_DEP}, libnektar++-solverutils (= ${NEKTAR_VERSION})")

    if (NEKTAR_SOLVER_INCNAVIERSTOKES)
        add_nektar_package(
            NAME nektar++-incnavierstokes-solver
            TYPE deb
            SUMMARY "Nektar++ incompressible Navier-Stokes flow solver"
            DESCRIPTION "
 Solves the incompressible Navier-Stokes equations using the Nektar++ high-
 order spectral/hp element framework"
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-incnavierstokes-solver_BINS}")
    endif (NEKTAR_SOLVER_INCNAVIERSTOKES)

    if (NEKTAR_SOLVER_ADR)
        add_nektar_package(
            NAME nektar++-adr-solver
            TYPE deb
            SUMMARY "Nektar++ advection-diffusion-reaction solver"
            DESCRIPTION "
 Solves a range of advection, diffusion, reaction PDE problems."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-adr-solver_BINS}")
    endif (NEKTAR_SOLVER_ADR)

    if (NEKTAR_SOLVER_CARDIAC_EP)
        add_nektar_package(
            NAME nektar++-cardiacep-solver
            TYPE deb
            SUMMARY "Nektar++ cardiac electrophysiology solver"
            DESCRIPTION "
 Solves the cardiac electrophysiology monodomain equations."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-cardiacep-solver_BINS}")
    endif (NEKTAR_SOLVER_CARDIAC_EP)

    if (NEKTAR_SOLVER_SHALLOW_WATER)
        add_nektar_package(
            NAME nektar++-shallowwater-solver
            TYPE deb
            SUMMARY "Nektar++ shallow water equations solver"
            DESCRIPTION "
 Solves the shallow water equations."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-shallowwater-solver_BINS}")
    endif (NEKTAR_SOLVER_SHALLOW_WATER)

    if (NEKTAR_SOLVER_COMPRESSIBLE_FLOW)
        add_nektar_package(
            NAME nektar++-compressibleflow-solver
            TYPE deb
            SUMMARY "Nektar++ compressible Navier-Stokes flow solver"
            DESCRIPTION "
 Solves the compressible Euler and compressible Navier-Stokes equations."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-compressibleflow-solver_BINS}")
    endif (NEKTAR_SOLVER_COMPRESSIBLE_FLOW)

    if (NEKTAR_SOLVER_PULSEWAVE)
        add_nektar_package(
            NAME nektar++-pulsewave-solver
            TYPE deb
            SUMMARY "Nektar++ pulse wave solver"
            DESCRIPTION "
 Models pulse wave propagation through a 1-dimensional network of arteries."
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-pulsewave-solver_BINS}")
    endif (NEKTAR_SOLVER_PULSEWAVE)

    if (NEKTAR_BUILD_UTILITIES)
        add_nektar_package(
            NAME nektar++-utilities
            TYPE deb
            SUMMARY "Nektar++ pre- and post-processing utilities"
            DESCRIPTION "
 This package contains the mesh conversion utility and other utilities for
 converting output files into common visualisation formats"
            DEPENDS "${NEK_DEP}"
            INSTALL_BINS "${nektar++-utilities_BINS}")
    endif (NEKTAR_BUILD_UTILITIES)

    # Top level meta package
    add_nektar_package(
        NAME nektar++
        TYPE deb
        SUMMARY "Nektar++ complete meta package"
        DESCRIPTION "
 This meta-package depends on the entire Nektar++ suite."
        DEPENDS "${NEK_DEP}, nektar++-incnavierstokes-solver (= ${NEKTAR_VERSION}), nektar++-adr-solver (=${NEKTAR_VERSION}), nektar++-cardiacep-solver (=${NEKTAR_VERSION}), nektar++-shallowwater-solver (=${NEKTAR_VERSION}), nektar++-compressibleflow-solver (=${NEKTAR_VERSION}), nektar++-pulsewave-solver (=${NEKTAR_VERSION}),  nektar++-utilities (=${NEKTAR_VERSION})")
endif()

if (NEKTAR_BUILD_LIBRARY AND RPMBUILD)
    ### RPM PACKAGES ###############################################
    add_nektar_package(
        NAME libnektar++-utilities
        TYPE rpm
        SUMMARY "Nektar++ library utilities"
        DESCRIPTION
        "This library provides core routines including linear algebra and
         integration with ThirdParty libraries."
        INSTALL_LIBS "${libnektar++-utilities_LIBS}")
    add_nektar_package(
        NAME libnektar++-stdregions
        TYPE rpm
        SUMMARY "Nektar++ StdRegions library"
        DESCRIPTION
        "This library provides construction of the reference expansions for the
         various 1D, 2D and 3D regions."
        INSTALL_LIBS "${libnektar++-stdregions_LIBS}")
    add_nektar_package(
        NAME libnektar++-spatialdomains
        TYPE rpm
        SUMMARY "Nektar++ SpatialDomains library"
        DESCRIPTION
        "This library provides the mappings between reference regions and
         physical regions in the domain."
        INSTALL_LIBS "${libnektar++-spatialdomains_LIBS}")
    add_nektar_package(
        NAME libnektar++-localregions
        TYPE rpm
        SUMMARY "Nektar++ LocalRegions library"
        DESCRIPTION
        "This library provides physical space expansions on the various
         supported regions."
        INSTALL_LIBS "${libnektar++-localregions_LIBS}")
    add_nektar_package(
        NAME libnektar++-nekmeshutils
        TYPE rpm
        SUMMARY "Nektar++ NekMeshUtils library"
        DESCRIPTION
        "This library provides mesh generation and library routines for the 
         NekMesh executable."
        INSTALL_LIBS "${libnektar++-nekmeshutils_LIBS}")
    add_nektar_package(
        NAME libnektar++-collections
        TYPE rpm
        SUMMARY "Nektar++ Collections library"
        DESCRIPTION
        "This library provides multi-expansion implementations of operators
         on the various supported regions."
        INSTALL_LIBS "${libnektar++-collections_LIBS}")
    add_nektar_package(
        NAME libnektar++-multiregions
        TYPE rpm
        SUMMARY "Nektar++ MultiRegions library"
        DESCRIPTION
        "This library provides global expansions on multi-element domains."
        INSTALL_LIBS "${libnektar++-multiregions_LIBS}")
    add_nektar_package(
        NAME libnektar++-globalmapping
        TYPE rpm
        SUMMARY "Nektar++ GlobalMapping library"
        DESCRIPTION 
        "This library provides deformation routines for Nektar++ simulations."
        INSTALL_LIBS "${libnektar++-globalmapping_LIBS}")
    add_nektar_package(
        NAME libnektar++-solverutils
        TYPE rpm
        SUMMARY "Nektar++ SolverUtils library"
        DESCRIPTION
        "This library provides support classes and routines for constructing
         complete spectral/hp element solvers."
        INSTALL_LIBS "${libnektar++-solverutils_LIBS}")

    if (NEKTAR_BUILD_DEMOS)
        add_nektar_package(
            NAME nektar++-demos
            TYPE rpm
            SUMMARY "Nektar++ demonstration applications"
            DESCRIPTION
            "This package contains a number of example programs for solving
             simple problems using the spectral/hp element method."
            INSTALL_BINS "${nektar++-demos_BINS}")
    endif (NEKTAR_BUILD_DEMOS)

    if (NEKTAR_SOLVER_INCNAVIERSTOKES)
        add_nektar_package(
            NAME nektar++-incnavierstokes-solver
            TYPE rpm
            SUMMARY "Nektar++ incompressible Navier-Stokes flow solver"
            DESCRIPTION
            "Solves the incompressible Navier-Stokes equations using the
             Nektar++ high-order spectral/hp element framework"
            INSTALL_BINS "${nektar++-incnavierstokes-solver_BINS}")
    endif (NEKTAR_SOLVER_INCNAVIERSTOKES)

    if (NEKTAR_SOLVER_ADR)
        add_nektar_package(
            NAME nektar++-adr-solver
            TYPE rpm
            SUMMARY "Nektar++ advection-diffusion-reaction solver"
            DESCRIPTION
            "Solves a range of advection, diffusion, reaction PDE problems."
            INSTALL_BINS "${nektar++-adr-solver_BINS}")
    endif (NEKTAR_SOLVER_ADR)

    if (NEKTAR_SOLVER_CARDIAC_EP)
        add_nektar_package(
            NAME nektar++-cardiacep-solver
            TYPE rpm
            SUMMARY "Nektar++ cardiac electrophysiology solver"
            DESCRIPTION
            "Solves the cardiac electrophysiology monodomain equations."
            INSTALL_BINS "${nektar++-cardiacep-solver_BINS}")
    endif (NEKTAR_SOLVER_CARDIAC_EP)

    if (NEKTAR_SOLVER_SHALLOW_WATER)
        add_nektar_package(
            NAME nektar++-shallowwater-solver
            TYPE rpm
            SUMMARY "Nektar++ shallow water equations solver"
            DESCRIPTION
            "Solves the shallow water equations."
            INSTALL_BINS "${nektar++-shallowwater-solver_BINS}")
    endif (NEKTAR_SOLVER_SHALLOW_WATER)

    if (NEKTAR_SOLVER_COMPRESSIBLE_FLOW)
        add_nektar_package(
            NAME nektar++-compressibleflow-solver
            TYPE rpm
            SUMMARY "Nektar++ compressible Navier-Stokes flow solver"
            DESCRIPTION
            "Solves the compressible Euler and compressible Navier-Stokes
             equations."
            INSTALL_BINS "${nektar++-compressibleflow-solver_BINS}")
    endif (NEKTAR_SOLVER_COMPRESSIBLE_FLOW)

    if (NEKTAR_SOLVER_PULSEWAVE)
        add_nektar_package(
            NAME nektar++-pulsewave-solver
            TYPE rpm
            SUMMARY "Nektar++ pulse wave solver"
            DESCRIPTION
            "Models 1-dimensional pulse wave propagation through an arterial network."
            INSTALL_BINS "${nektar++-pulsewave-solver_BINS}")
    endif (NEKTAR_SOLVER_PULSEWAVE)

    if (NEKTAR_BUILD_UTILITIES)    
        add_nektar_package(
            NAME nektar++-utilities
            TYPE rpm
            SUMMARY "Nektar++ pre- and post-processing utilities"
            DESCRIPTION "This package contains the mesh conversion utility and
            other utilities for converting output files into common
            visualisation formats"
            INSTALL_BINS "${nektar++-utilities_BINS}")
    endif (NEKTAR_BUILD_UTILITIES)

    add_nektar_package(
        NAME nektar++
        TYPE rpm
        SUMMARY "Nektar++ complete meta package"
        DESCRIPTION "
 This meta-package depends on the entire Nektar++ suite."
        DEPENDS "nektar++-incnavierstokes-solver = ${NEKTAR_VERSION}, nektar++-adr-solver = ${NEKTAR_VERSION}, nektar++-cardiacep-solver = ${NEKTAR_VERSION}, nektar++-shallowwater-solver = ${NEKTAR_VERSION}, nektar++-compressibleflow-solver = ${NEKTAR_VERSION}, nektar++-pulsewave-solver = ${NEKTAR_VERSION}, nektar++-utilities = ${NEKTAR_VERSION}")
endif()

if (NEKTAR_BUILD_LIBRARY AND PKGMAKER)
    if (NEKTAR_USE_MPI)
        MESSAGE(WARNING "Shipping OS X with MPI is unlikely to work!")
        FIND_PROGRAM(HAVE_MPIRUN mpirun)
        IF(HAVE_MPIRUN)
            LIST(APPEND nektar++_BINS ${HAVE_MPIRUN})
        ENDIF()
	UNSET(HAVE_MPIRUN CACHE)
    endif ()
    if (NEKTAR_USE_PETSC)
        MESSAGE(WARNING "Shipping OS X with PETSc is unlikely to work!")
    endif ()
    add_nektar_package(
        NAME nektar++-serial
        TYPE pkgmaker
        SUMMARY "Nektar++ libraries and solvers"
        DESCRIPTION "This is all of Nektar++"
        INSTALL_LIBS "${nektar++_LIBS}"
        INSTALL_BINS "${nektar++_BINS}")
ENDIF()

# TGZ package
add_nektar_package(
    NAME nektar++
    TYPE tgz
    SUMMARY "Nektar++ libraries and solvers"
    DESCRIPTION "This is all of Nektar++"
    INSTALL_LIBS "${nektar++_LIBS}"
    INSTALL_BINS "${nektar++_BINS}")

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
    "/\\\\.clang-format"
    "/\\\\.gitignore"
    "/\\\\.gitmodules"
    "/build/"
    "/builds/"
    "/docs/tutorial"
    "/ThirdParty/"
    "/Testing/"
    "/library/Demos/MultiRegions/ExtraDemos/"
    "/library/Demos/MultiRegions/XmlFiles/"
    "/library/Demos/LocalRegions/XmlFiles/"
    "/solvers/ImageWarpingSolver/"
    "/solvers/VortexWaveInteraction/"
    "/utilities/Extras/"
)


include(CPack)
