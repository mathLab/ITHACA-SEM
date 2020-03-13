Changelog
=========

v4.4.0
------
**Library**:
- Add support for variable polynomial order for 3D simulations with continuous
  Galerkin discretisation (!604)
- Bump version of gsmpi to suppress autotuning output unless `--verbose` is
  specified (!652)
- Add support for variable polynomial order with periodic boundary conditions
  (!658)
- Statistics are now printed for lowest level of multi-level static condensation
  (!656)
- Sped up interpolataion from pts files and fixed parallel pts import (!584)
- Increased required boost version to 1.56.0 (!584)
- New FieldUtils library allows support for most `FieldConvert` post-processing
  operations during simulation using a new filter (!589)
- Adjust CMake dependencies to reduce compile time (!671)
- Homogeneous1D dealiasing improvements (!622)
- Add support for HDF5 as an alternative output to XML-based output, including
  refactoring of FieldIO, improvements to MPI interface and added communicators
  to boundary conditions (!615)
- Allow expansions to be loaded directly from field file (!617)
- New options for load balancing (DOF or BOUNDARY) in mesh partitioner (!617)
- Rework nodal utilities to support nodal prismatic elements (!660)
- Update Body/Field forces at each timestep (!665)
- Update nodalutil to include quad and hex elements and introduce SPI nodal
  points (!696)
- Add ability to restart time-averaging and Reynolds stresses from checkpoint
  file (!678)
- Extend ExtractDataToCoeffs to support interpolation between basis types for
  quads and hexahedra (!682)
- Enabled MUMPS support in PETSc if a Fortran compiler was found and added 3D
  support to the Helmholtz smoother used e.g. in FieldConverts C0Projection
  module (!714)
- Fix bug in `Vmath::FillWhiteNoise` which caused `ForcingNoise` to have
  a repeated pattern (!718)
- Fix bug in the calculation of the RHS magnitude in CG solver (!721)
- Fix bug in MPI detection for recent CMake on OS X (!725)
- Fix bug in CMake Homebrew and MacPorts detection for OS X (!729)
- Fix bug in FieldUtils when using half mode expansions (!734)
- Do not read the same fld/pts files again for every variable (!670)
- Fix bug in CMake PETSc detection for Ubuntu 16.04/Debian 9 (!735)

**ADRSolver:**
- Add a projection equation system for C^0 projections (!675)

**APESolver:**
- Use a continuous basefield projection and revert to constant c formulation (!664)
- Added ability to compute CFL number (!664)
- Output Sourceterm (!664)
- Use the Forcing framework to define source terms (!665)

**IncNavierStokesSolver:**
- Add ability to simulate additional scalar fields (!624)
- Improve performance when using homogeneous dealiasing (!622)
- Fix linearised advection for full 3D cases (!708)
- Added a weak pressure formulation following Guermond & Shen (!713)
- Added a convective like outflow boundary condition from Dong (!713)
- Added the ability to specifiy Womersley boundary conditions for pulsatile flow (!472)

**CardiacEPSolver:**
- Added a Python translator utility to generate cell models from CellML (!723)

**FieldConvert:**
- Allow equi-spaced output for 1D and 2DH1D fields (!613)
- Update quality metric to include scaled Jacobian output (!695)
- Allow multiple XML files to be specified in InterpField module (!705)
- Fix issues with isocontour module (!719)
- Fix issue with interpolator routine (!746)

**NekMesh:**
- Modify curve module to allow for spline input (!628)
- Add STL surface writer module (!668)
- New module for inserting an alternate high-order surface into the working
  mesh (!669)
- Add curve projection routines to CAD system (!697)
- Extensive clean-up of NekMeshUtils/MeshElements and extension of makeorder to
  consider CAD information (!698)
- Improvements to mesh linearisation module (!659)
- Add support for Gmsh high-order output (!679)
- Move CAD classes to factory format (!676)
- Add module to check topology of the mesh along with boundary connectivity
  to detect problems such as hanging nodes (!691)
- Add option to `linearise` module to linearise only prisms (!688)
- Add reader for Nek5000 mesh files (!680)
- Add option to `linearise` to use element quality (!690)
- Add flag to `insertsurface` process for non-conforming geometries (!700)
- Bug fix to get two meshgen regression tests working (!700)
- Remove libANN in deference to boost::geometry (!703)
- Refactor library to use NekMesh modules for CAD generation (!704)
- Add `varopti` process module to optimise meshes (!711)
- Add a mesh extract option to the linearise module to visualise the result
  (!712)
- 2D to 3D mesh extrusion module (!715)
- Add new two-dimensional mesher from NACA code or step file (!720)
- Fix inverted boundary layer in 2D (!736)
- More sensible element sizing with boundary layers in 2D (!736)
- Change variable names in mcf file to make more sense (!736)
- Fix issues in varopti module so that in can be compiled without meshgen on
  (!736)
- Replace LAPACK Eigenvalue calculation with handwritten function in 
  varopti (!738)
- Improved node-colouring algorithm for better load-balancing 
  in varopti (!738)
- Simplified calculation of the energy functional in varopti for improved
  performance (!738)

**FieldConvert:**
- Move all modules to a new library, FieldUtils, to support post-processing
  during simulations (!589)
- Add module to stretch homogeneous direction (!609)
- Add module to add composite ID of elements as a field (!674)
- Add reader for Nek5000 field files (!680)

**Tester:**
- Fix output not displayed on segfault or system error (!745)

v4.3.5
------
**Library:**
- Fix bug in DG with hybrid meshes (!694)
- Fix issue with parallel output (!699)
- Fix performance issue with iterative full solver (!693)
- Enforced precision on history point output (!706)

**Documentation**
- Update build instructions in user guide for Windows (!692)

**Tester**
- Fix bug in tester when no parameters specified for test executable (!701)

v4.3.4
------
**Library:**
- Fix performance issue with `v_ExtractDataToCoeffs` for post-processing of
  large simulations (!672)
- Added additional assertions to ensure homogeneous simulations have an even
  number of planes per process (!666)
- Fix compilation with NEKTAR_USE_MESHGEN option
- Fix IterativeFull solver in parallel (!685)
- Fix error message for missing fld file (!689)

**IncNavierStokesSolver:**
- Fix 2nd order time-integration for VCSMapping (!687)

v4.3.4
------
**Library:**
- Fix performance issue with `v_ExtractDataToCoeffs` for post-processing of large
  simulations (!672)

v4.3.3
------
**Library**:
- Auto-detect a shared filesystem and removed --shared-filesystem option (!654)
- Fix filters when using adaptive driver to avoid output being overwritten after
  each adaptive update (!588)
- Minor fix to suppress Xxt output unless `--verbose` is specified (!642)
- Fix of DirectFull solver in case where only Neumann boundary conditions
  are imposed. (!655)

**FieldConvert**:
- Fix to avoid repeated import of field file (!649)
- Fix issue with C^0 projection (!644)
- Fix verbose output when using --procid (!648)

**NekMesh:**
- Fix namespace issue in Star-CCM+ input header in NekMesh (!661)

**CompressibleFlowSolver**:
- Fix issue with residual output (!647)
- Issues with 1D Euler solver fixed (!565)
- Fix deadlocking issue with boundary conditions (!657)

**Packaging**:
- Fix NekMesh dependencies for DEB package (!650)
- Fix PETSc build on newer linux distributions (!646)

v4.3.2
------
**Library**:
- Add small optimisation for DriverAdaptive (!618)
- Updated FFTW build to use the compiler used for building Nektar++ (!629)
- Fix numbering bug in periodic boundary conditions (!631)
- Print error message for invalid equation also in release version (!634)
- HistoryPoints filter now uses closest plane to requested z-coordinate and
  output is produced in physical space (!621).
- Fix minor performance issue with time integration schemes (!632)
- Fix FilterCheckpoint filter to be consistent with `IO_CheckSteps` (!633)
- Fix CMake configuration for building on Windows 10 with VS 2015 (!641)
- Fix `IO_CheckSteps` to avoid missing first checkpoint (!639)
- Fix bug in iterative solver where only root process would ASSERT when
  exceeding the maximum number of iterations (!636)

**FieldConvert**:
- Fix appearence of duplicate messages when running in parallel (!626)
- Fix issue with efficiency when using large number of 3DH1D planes (!627)
- Add module for combining average fields (!620)
- Fix wall shear stress processing module for parallel execution (!635)

**Packaging**:
- Fixes for DEB package dependencies (!630)

v4.3.1
------
**Library**:
- Add `THIRDPARTY_USE_SSL` option to disable use of SSL on systems where CMake
  is not compiled with SSL support. (!602)
- Fixed a number of documentation issues (!586, !593, !596)
- Fix Homogeneous transform when unshuffling is not used. (!599)
- Fix namespace pollution in library header files. (!601)
- Fix issue with METIS compilation on clang 7.3 (!603)
- Fix issue with heterogeneous quadrilaterals (!607)
- Fix bug in modified Arnoldi algorithm causing convergence to be reported when
  number of vectors is less than `nvec` (!608)
- Fix uninitialised array bug in AssemblyMap (!598)
- Fix issue with LAPACK call in eigenvalue calculation (!610)
- Fix FieldConvert processing of partitions in serial (!612)
- Fix use of multi-level static condensation in parallel with periodic
  boundary conditions (!614)
- Fix NaN detection to work in parallel (!605)
- Add additional constructor to ContField3DHomogeneous1D for FieldConvert
  extract module. (!590)

**NekMesh**:
- Fix incorrect link directory on CCMIO library.

**FieldConvert**:
- Fix to FLD input to update the field definitions always, not just when a range
  is specified. (!611)

**Tester**:
- Remove requirement for executable to be specified in .tst file if it is
  overridden on the command-line (!595)

**Packaging**:
- Fix dependency resolution on generation of DEB packages. (!616)

v4.3.0
------
**Library:**
- Changed default XML format to compress mesh data (!533, !547)
- Various fixes for 3D homogeneous post-processing (!531, !529, !528, !526, !521)
- Fix boundary condition imposition for 3D homogeneous 2D HelmSolve (!545)
- Fix range with variable p option (!522)
- Fix bug with hexahedra of heterogeneous order (!520) and reading files (!522)
- Fix history point output formatting (!518)
- Fix for OS X 10.11 (!512)
- Fix `HexGeom::v_GetDir` to support heterogeneous basis functions (!520)
- Added new `NekMeshUtils` library to support new `NekMesh` executable and
  associated CAD routines. Old CAD wrappers in LibUtilities now moved to
  `NekMeshUtils` (!527)
- Fix initialisation bug in ExpList2DH1D and ExpListHomogeneous2D (!528, !529)
- Fix bug in ExpList1D which may lead to invalid .vtu files (!531)
- Make `GetBoundaryToElmtMap` consistent for 3DH1D (!526)
- Add support for PETSc matrix shell to use Nektar++ operations/preconditioners
  (!537)
- Fix bug with initial conditions of CG simulations using variable P (!543)
- Fix bug in 3DH2D with non-zero Dirichlet boundary conditions (!545)
- Added in a method to convert equispaced interpolated points back to
  coefficients which requires the introduction of a new StdRegions matrix.(!561)
- Empty XML tags which would override non-empty XML tags are now ignored (!581)
- Add contribution guide (!551)
- Add a filter to calculate exponential moving averages (!566)

**APESolver:**
- Fix restarting from checkpoint file (!517)

**IncNavierStokesSolver**
- Fix floquet stability analysis for HalfMode case (!536)
- Add a filter to calculate Reynolds stresses (!566)

**FieldConvert:**
- Extended surface distance module to support hexahedra and quads (!524)
- Small fixes in interpolation routine (!515)
- Add support for surface extraction in 3DH1D case (!521)
- Add support for isocontour extraction for 3DH1D (!525)
- Add process module to calculate high-order mesh quality metric (!527).
- Add module to extract one of the planes of 3DH1D (!542)
- Add module to enable mean mode of 3DH1D to be extracted (!530)
- Fix bug in C^0 projection (!541))
- Add command line option to set number of homogeneous planes (!540)
- Add module to project set of points to a fld file(!561)
- Add support for interpolating to a box of points and fix ability to run
  interppointstofld module in parallel when using a plane or box option (!561)
- Add option to output equi-spaced points in VTU format (!550)
- Add module innerproduct (!568)
- Add command line option of `--part-only` and `--part-only-overlapping` (!569)

**NekMesh:**
- `MeshConvert` is now renamed to `NekMesh` to reflect new mesh generation
  functionality (!527).
- Enable face curvature inside core MeshConvert objects (!511)
- Add linearise processing module to remove all curvature from high order
  elements (!509)

**Documentation:**
- Added git submodule for including Nektar++ tutorials in the source tree (!507)

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
- A range of fixes for the coupled stability solver, which now works in parallel
  (!508)

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

**Library:**
- Change hybrid parallelisation to use command line options (!368)
- Add support for multi-variable functions in expression evaluator: new
  functions include rad and ang for polar coordinates (!375)
- Add more documentation (!376, !383)
- Various OS X (!377, !378, !382, !425), compiler warning (!432), documentation
  (!434) Windows 7 (!391, !407), CMake (!392, !415), packaging (!435, !436) and
  Intel compiler (!414, !416) fixes
- Refactor of CG and DG assembly maps (!380)
- Fixes for PETSc running in serial (!381, !420)
- Fixes for running Arnoldi solver in parallel (!384)
- Enable MPI tests on Cray machines such as ARCHER (!386)
- Fix issues with extracting face physical values (!393)
- Fix threshold filter (!395)
- HDG can now use block preconditioner (!397)
- Fix issue with singular vertices in parallel (!398)
- Timing executables now use `Timer` class from LibUtilities (!402)
- Fix manifold history points again (!410)
- Fix time output inside energy filter (!412)
- Fix GetExpIndex function (!417)
- Fixes to external project compilation (!419)
- Fixes from CPC paper review (!422)
- Fixes for scotch partitioner tests (!423)
- Fixes for ACML BLAS libraries (!424)
- Allow prepartitioned meshes to be used (!426)
- Enable variable names to be remapped inside files to different names in XML
  functions (!428)

**APESolver:**
- Fixes for tests (!404)
- Add support for advection classes (!408)

**CardiacEPSolver:**
- Add benchmark (!411)
- Fix cardiac exmplaes (!418)

**CompressibleFlowSolver:**
- Add filter for kinetic energy/enstrophy calculation (!388)

**FieldConvert:**
- Support equi-spaced output for simplex elements to reduce storage (!421)

**IncNavierStokesSolver:**
- Unify advection classes with those in `SolverUtils` (!403, !408)

**MeshConvert:**
- Boundary layer refinement now supports hexahedra (!390)
- Improve support for Gmsh high order elements (!401)
- Many fixes for face-interior curvature (!401)
- Add rudimentary test suite (!401)
- New module for imposing curvature based on a scalar function (!401)

v4.0.0
------

**Library:**
- Update boost to 1.55 (!289)
- Fix parallel history points on manifold (!298)
- Add support for scotch partitioner (!311)
- Fixes for thirdparty builds (!319, !330, !353)
- Fix CMake >= 3.0.0 warnings (!320)
- Add support for PETSc library and tidy up global system classes (!322)
- Fixes for 1D Helmholtz solver (!326)
- Fixes for history points (!327) and solver output (!331)
- Fix issue with mesh IDs that do not start from zero (!354)

**CardiacEPSolver:**
- Simplify support for global conductiity (!295)

**FieldConvert:**
- Fixes for parallel operation and interpolation of points (!351)

**IncNavierStokesSolver:**
- Fixes for sponge layer (!272)
- Fix setting of initial conditions (!298)

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
