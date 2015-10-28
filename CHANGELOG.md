Changelog
=========

v4.3.0
------


v4.2.0
------

**Library:**
- Add Runge-Kutta SSP schemes for 2nd/3rd order using keys `RungeKutta2_SSP` and
  `RungeKutta3_SSP`. `ClassicalRungeKutta4` is now called `RungeKutta4`. (!481)
- Add rudimentary support for 3D CAD models using OpenCascade - work in progress
  (!486)
- Allow filters to evaluate expressions in their parameter definitions (!489)
- Fix block preconditioner to work with periodic boundary conditions (!420)
- Dump a backtrace when crash occurs and Nektar++ is compiled in FullDebug mode
  (!495)
- Stop the execution of a time-dependent solver if NaN is detected in the
  solution field (!496)
- Fixes to improve robustness of interpolation routines (!499)
- Allow solvers to use multi-level static condensation with Xxt, most useful
  when running a 3DH1D simulation (!502)

**IncNavierStokesSolver:**
- Various fixes for the coupled stability solver (!508)

**MeshConvert:**
- Add module to extract prismatic boundary layer elements from mixed prism-tet
  mesh (!493).

**FieldConvert:**
- Add a processing module to calculate height of an element connected to a
  surface, allowing for calculation of y plus values (!488)
- Fixes for equispaced output (!510)


v4.1.0
------

**Library:**
- Add support for interpolating point data from .pts files (!433)
- Fixes for curvilinear element normals (!443)
- Fix consistency issues between FFT and MVM approaches for homogeneous
  expansions (!444)
- Fix a bug in Tecplot output (!445)
- Fix a bug with PETSc and MPI_Finalize (!456)
- Fix bugs with mesh partitioning (!449, !480)
- Fix a bug with non-symmetric SVV parameters for curvilinear elements (!451)
- Fix detection of Intel MKL 2013/2015 (453)
- Fix linearised stability solver in parallel (!454)
- Add a filter for 1D energy spectra (!457)
- Add an incomplete developer guide containing most information from the wiki
  (!459)
- Change user defined boundary conditions to remove dependency on enumerator
  inside SpatialDomains (!460)
- Add a new collections library for optimised evaluation of operators (!461)
- Change minimum version of boost to 1.52.
- Add initial multithreading support (!463)
- Fix third-party boost compilation on OS X (!467)
- Disable some regression tests on 32-bit systems (!468)
- Fix memory issues inside collections (!473)
- Fix collections autotuning (!476)
- Fix VtkToPng utility (!477)
- Add PulseWaveSolver to packaging (!478)
- Fix bug in iterative static condensation solver (!483)
- Fix zlib install path on OS X (!484)
- Fix documentation HTML styling for user and developer guide (!485)
- Add fixes to support native Nektar++ extension in VisIt visulisation software
  (!490)
- Fix warnings on OS X (!491)

**CardiacEPSolver:**
- Fixes for stimuli (!442, !446), conductivity (!441), cell restarts (!458)
- Add a new filter for outputting cell states at specific points over time (!465)

**Linear elastic solver (new):**
- Add solver for linear elasticity equations (!400)

**IncNavierStokesSolver:**
- Add support for moving bodies (!344, !448)
- Fixes for modal energy filter (!427)
- Fix import of mesh file in the Adaptive SFD driver (!440) and other general
  fixes (!452)
- Documentation for high order pressure and outflow boundary conditions (!447)
- Update examples to use correct forcing terms (!470)
- Fixes for half-mode stability (!471)
- Fix static initialisation problem in extrapolation classes (!492)

**CompressibleFlowSolver:**
- Add support for sponge region (!396)
- Add support for adiabiatic walls (!430)
- Add utility to generate boundary layer from similarity solution (!438)

**ShallowWaterSolver:**
- Added a DG solver for the Boussinesq equations of Peregrine (!431)

**APESolver:**
- Add support for variable speed of sound (!438)

**MeshConvert:**
- Fix Star file input for highly stretched elements (!455)
- Add Star input from binary format (!474)
- Tidy up files to align with FieldConvert (!479)

**FieldConvert:**
- Major re-organisation of modules, most post-processing utilities now available
  within FieldConvert (!475)

v4.0.1
------



v4.0.0
------



v3.4.0
------

**Library:**
- New parallel output format. Parallel files are now stored in directories which
  contain partition information. (!100, !102, !236, !242, !249, !256).
- gzip-compressed XML mesh files are now supported with extension .xml.gz (!116,
  !140, !186).
- HDG solvers now run in parallel and have post-processing utilities (!188,
  !230).
- Partitioning can be done only on root process if shared filesystem is
  present with use of `--shared-filesystem` command line option (!220, !250).
- A variety of preconditioners are now supported, including linear space and
  low-energy preconditioning (!148).
- Many changes to geometric factors storage and interpolation (!99, !197).
- Improvements to identification of invalid elements (!208, !227).
- Removed elemental storage to reduce memory consumption by 30-50% for large
  problems (!240).
- Various performance and design improvements for discontinuous formulation (!134).
- Periodic boundary conditions are supported in 3D for both continuous and
  discontinuous formulations (!139, !150, !152, !155, !196).
- Utilities added to mesh converter to help identify pairs of periodic faces
  (!214).
- Preconditioner support for periodic boundary conditions (!231, !239).
- New radiation boundary condition type (!74).
- Some solvers (compressible flow solver, advection-diffusion-reaction solver)
  now support dealiasing options (!78, !146, !167).
- BLAS and vectorisation performance improvements for static-condensed iterative
  solver (!86, !109).
- New driver to improve steady state convergence and add parallel support (!91,
  !235).
- Updated to METIS v5.1.0 (!97, !142, !189).
- Iterative solvers now use previous timestep (when available) to improve
  convergence speed (!106).
- Added CPU timing for timestep loop (!156).
- Added provenance information (date, time, code version, git revision, etc) to
  field file output (!179).
- Disabled long-running regression tests by default (!183).
- Support for command line arguments without parameters (!187).
- Added support for reading boundary conditions from files, and appropriate
  utilities in MeshConvert to extract surfaces (!226).
- Updated XXt and Gs libraries to latest version (!232).
- Fix singularity check for Poisson equations (!74, !154).
- Fixes for 2D Gauss points (!73, !149, !157).
- Fixes to parallel I/O (!77, !218, !264).
- Fixes for parallel implementation (!93, !107, !121, !169, !217, !245, !246).
- Fixes for normal calculation (!94, !135).
- Improved compilation techniques, particularly when compiler includes MPI
  automatically (!80, !82, !84, !85, !113, !114, !131, !141, !166, !210, !241).
- Updated zlib to v1.2.7 (!115).
- Fix for boost 1.5.3 compilation (!120).
- Most compiler warnings silenced with clang/gcc (!81, !92, !103, !123, !201,
  !243).
- Attempts to improve mesh partitioning/load balancing (!160, !170, !175).
- Fixes for Newton iteration to interpolate inside deformed elements (!216,
  !251).
- Fixed curved tetrahedron and hexahedron issue (!219, !248).
- Fixed reading of field files for tetrahedron (!228).
- Fixed uninitialised variable inside SessionrReader (!233).
- Various improvements to support use of Nektar++ externally (!111, !260, !261).
- Fixed base flow reading (!112).

**CardiacEPSolver:**
- Cardiac electrophysiology solver improvements (!87, !95, !96, !108, !119,
  !165, !173, !174, !199, !222).

**CompressibleFlowSolver:**
- Compressible Navier-Stokes equations are now available for both DG and FR
  discretisations (!110, !125, !128).
- Meshes with spatially varying p in both 2D and 3D are now supported (!158).
- Homogeneous Fourier extension is now supported (!180).
- Various fixes (!90, !98, !147, !172).

**DiffusionSolver (new):**
- Added small solver to demonstrate usage of higher library levels outside of
  EquationSystem (!225).

**IncNavierStokesSolver:**
- Major refactoring of time-integration classes (!181, !184).
- Summary information now generated via callbacks (!182).
- Implemented new generic forcing function classes (!194).
- Current time now written out in field files (!198).
- Major refactoring of incompressible Navier-Stokes solver to improve
  readability and performance (212, !213).
- Spectral vanishing viscosity for stabilisation (!101, !104, !211, !263).
- Added filter to compute aerodynamic forces on surfaces (!168, !203, !204).
- Added filter to compute kinetic energy and enstrophy (!207, !257).

**ShallowWaterSolver:**
- Various improvements/modernisations to shallow water solver (!190).

**Utilities:**
- VTK to PNG converter (!122)
- Added scalar gradient utility (!129, !252).
- Added utility to calculate Q-criterion field (!153).
- Added support to XmlToVtk to write Jacobian field (!223).
- Added utility to calculate wall shear stress (!224).
- Fixed vorticity calculator (!138).

**MeshConvert:**
- Added face-interior quadrature and 2D/3D manifold support to spherigon code
  (!130).
- Fixes for boundary layer refinement and prism-to-tetrahedron splitting (!137,
  !202, !206, !244).

**FieldConvert (new):**
- Added new FieldConvert utility which will eventually encompass most existing
  utilities (!255).
