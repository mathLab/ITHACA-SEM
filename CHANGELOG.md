Changelog
=========

v5.1.0
------
**Library**
- Restructure library to use local coefficient storage down to the GlobalLinSys
  level. Removed GlobalCeoffs functionality (!963)
- Add interior penalty method to DG framework (!1101)
- Add an error filter for the time-evolution of the L2 and Linf errors (!1147)

**FieldConvert**
- Refactored time integration code using factory pattern (!1034)
- Fix to preprocessor logic for boost with Visual Studio >= 2015 (!1115)
- Fix type consistency and real comparison in SharedArray.hpp, replaced
  num_elements with size() (!1127)
- Use base MPI functions instead of the GS library in the trace exchange
  for parallel DG simulations (!1112)
  num_elements with size() (!1127, !1137, !1141)

**CardiacEPSolver**
- Added additional parameter sets to Fenton-Karma model (!1119)

**NekMesh**
- Improved boundary layer splitting and output to CADfix (!938)
- Improve .geo reader and support 3D geometries with voids (!1031)
- Added r-adaptation code (!1109)

**BuildSystem**
- Toggle build type (!1135)

v5.0.1
------
**Library**
- Fix incorrect coordinate dimension used in history point filter (!1118)
- Fix compile errors with GCC 9.x (!1108)
- Correct the Energy/Enstropy integral for the 3DH1 flow (!1132)
- Added IsRealEqual method to compare real numbers with relative tolerance.
  Started using it in SharedArray and in NekMesh to fix peralign-extrude tool
  chain (!1134)

**IncNavierStokesSolver**
- Change the baseflow time in the Adjoint advection (!1133)

**FieldConvert**
- Fix OutputTecplot skipping final plane in 3DH1D (!1016)
- Fix Interppoints in 3DH1D (!1140)

**NekMesh**
- Fix compile errors when using intel cc (!1114)

**Documentation**
- Fix error in compilation of developer guide (!1136)

**CI**
- Added checked conversion from double to int in SessionReader (!1113)
- Switched to Gitlab CI (!1120, !1120, !1128, !1129, !1131, !1141)

v5.0.0
------
**Library**
- Added in sum factorisation version for pyramid expansions and orthogonal
  expansion in pyramids (!750)
- Added detection of 'abort' file to cleanly terminate simulation early (!772)
- Significant overhaul of CMake infrastructure (!770, !804)
- Fix ThridpartyCCM options (!802)
- Fix Windows CRLF tokens in GEO reader and improve comment handling (!805)
- Use chrono in Timer (!807)
- Fix caching of FUNCTION tags that read from file and provide the same
  functionality in FUNCTIONs defined for forcings (!759)
- Transition to C++11 (!795, !847)
- Add patch to tinyxml to fix size_t vs int bug (!820, !1006)
- Add ARPACK thirdparty build capabilities (!828)
- Added native support for csv files in addititon to pts (!760, !835, !906)
- Utilize LAPACK_DIR env variable to find the native blas/lapack install (!827)
- Extend AeroForces filter to compressible flows (!815)
- Remove StdExpansion use from MultiRegion (use Expansions instead). (!831)
- Move steady state check and CFL output from solvers to SolverUtils (!832)
- Remove DG advection implementation from EquationSystem (!832)
- Simplify RawType typedefs (!840)
- Remove unused files from BasicUtils (!841)
- Remove checks for old boost versions which are no longer supported (!841)
- Refactor ParseUtils to be more consistent (!843, !896, !908)
- Added support for using the distance to a specific region (e.g. outlet) in the
  function definitions for the Absorption Forcing (!769)
- Improve performance of DisContField2D::v_ExtractTracePhys (!824)
- Fix small bug in Jacobian Energy (!857)
- fix variable name overriding in file functions (!870)
- Adds CFI CAD engine back-end (!864)
- Adds CFI Mesh IO support (!864)
- Cleanup of CAD system data structures (!864)
- Fix mac OSX on buildbots (!876)
- Fix error from (!826) (!876)
- Fix minor bug in ARPACK thirdparty build cmake (!874)
- Added in sum factorisation version for pyramid expnasions and orthogonal
  expansion in pyramids (!750)
- Adjust boost third-party compilation to account for different toolset
  choices (!886)
- Switch MeshGraph to use factory pattern and add HDF5 geometry support (!900,
  !904, !941)
- Restructure the low energy preconditioner to handle pyramidic and variable
  p expansions (!920)
- Remove requirement for modmetis, switch to SCOTCH by default (!899)
- Switch MeshGraph to use factory pattern and add HDF5 geometry support
  (!900, !904)
- Fix bug in MeshPartition.cpp which caused incorrect array access when
  WeightPartitions was used in parallel (!923)
- Removed instance count from beginning of Array storage to improve memory
  alignment (!921)
- Fix naming issue of duplicate Unit tests (!924)
- Fix warnings about missing virtual destructors in abstract classes (!932)
- Fix ability to have periodic boundary conditions that are aligned by a
  rotation rather than just a translation (!933)
- Added a coupling interface to exchange data between solvers at run time
  and a DummySolver to test the implementations (!853, !931, !950, !973, !1017)
- Fix compilation issue with newer Boost versions and clang (!940)
- If only `NEKTAR_BUILD_LIBRARY` is enabled, only libraries up to and including
  `MultiRegions` will be built by default (!945)
- Dont add doxygen documentation to the all target (!834)
- Fix missing metadata import from Hdf5 files (!971)
- Fix missing flags for periodic BC in DiffusionLDG (!985)
- Add the moving reference frame as a forcing (!987)
- Added rtree for element bounding box lookup to accelerate interpolation (!996,
  !1066)
- Fix integration weights on prisms and pyramids if not using the default
  integration rule (!998)
- Fix missing ContainsPoint in Pyramid expansion (!1000)
- Added path prefixes to find packaged Scotch (!979, !1008)
- Add HDF5 geometry format (!977)
- Combine and generalise demo code in StdRegions and LocalRegions (!993)
- Fix for error output to allow for custom error streams (!944)
- Fixed bug in ReOrientQuadFacePhysMap (!1003)
- Add NekPy Python interface (!962, !990, !989, !1004, !1014, !1061, !1070)
- Fix edge case for ThirdPartyScotch and FindScoth (!1009)
- Fix to populate m_elmtToExpId map if not already set up in GetExpIndex (!1019)
- Added flag to skip periodic BCs while filling Dirichlet BCs in
  ContField3D.cpp (!1018)
- Fix bounding box for interpolation (!1033)
- Added IMEXOrder4, RK5 and AB4 time integration schemes (!1037)
- Fix TriExp.cpp orientation bug (!1048)
- Fix XML attributes in conditions.cpp to be unordered (!1015)
- Fix issue with HDF5 mesh input in serial (!1049)
- Add estimate of filters CPU time (!1044)
- Update CompressibleFlowSolver/Examples/Test_IsentropicVortex1.xml example (!1045)
- Add error if HDG used with periodic BCs (!1071)
- Fix issues related to leading factors, arithmetic order and associativity of
  exponential operator in expression evaluator (!1066)
- Remove use of `using namespace std` in header files (!1066)
- Add error messages for use of ARPACK in serial (!1079)
- Generalise ContainsPoint routine (!1078)
- Homogenized fallthrough to fix issues with gcc 7.4.0 (!1084)

**NekMesh**:
- Add feature to read basic 2D geo files as CAD (!731)
- Add periodic boundary condition meshing in 2D (!733)
- Adjust boundary layer thickness in corners in 2D (!739)
- Add non-O BL meshing in 2D (!757)
- Add ability to compile CCIO library but tar file is not yet openly
  available whist we seek permission from Simens (!799)
- Fix issue with reading CCM files due to definition of default arrays
  rather than a vector (!797)
- Fix inverted triangles and small memory issue in surface meshing (!798)
- Update for the CAD system, more advance self-healing and analysis (!822)
- Additional curve types in GEO reader: BSpline, Circle, Ellipse (!800)
- Fix default command line argument value (!823)
- Add projection meshing module which can curve linear meshes with CAD (!826)
- XML meshes now write with provenance information, including information about
  their source, for debugging purposes (!872)
- Force 3-node loops to avoid degenerate 1-triangle faces (!875)
- Smooth BL normals in 2D when normals intersect or cause invalid macro BL
  elements (!877)
- Revert triangle code to ThirdParty library (!883)
- Fix coinciding nodes issue with very fine meshes (!883)
- Skip CFI groups of bodies and non-numbered nodes (!891)
- Add ability to space out 2D BL nodes to better fit local target Delta (!890)
- Fix automatic peralign call in 2D periodic meshing (!888)
- Fix BL splitting call from MCF (!910)
- Support CFI combined lines (!917)
- Order nodes in Gmsh output (!912)
- Fix manifold face curvature nodes (!913)
- Fix writing 1D surfaces (!930)
- Fix surface string parsin in BL splitting (!937)
- Enable use of distributed packages for triangle and TetGen (!953)
- Fix issue with MLSC after Scotch conversion (!943)
- Add support for Gmsh 4.0 mesh file format (!964)
- Fix issue with extracting 1D curved surface from 2D file (!984)
- Fix surface extraction, added regression test (!994)
- Fix 2D meshing running out of memory due to missing else (!1012)
- Add support for .msh v4.1 file input (!1054)
- Added penalty term to LDG and LDGNS, slight generalization of LDG (!1080)

**FieldConvert**:
- Add input module for Semtex field files (!777)
- Fixed interppoints module (!760)
- Fix OutputTecplot in 2DH1D (!818)
- Move StreamFunction utility to a FieldConvert module (!809)
- Allow using expansion from session file with new `--useSessionExpansion`
  command line option (!842)
- Extend wss module to compressible flows (!810)
- Allow explicitly setting bool options of FieldConvert modules as false (!811)
- Enable output to multiple files (!844)
- Allow using xml file without expansion tag in FieldConvert (!849)
- Add Lambda 2 vortex detection criteria (!882)
- Add module for modifying/adding fields from expressions (!889, !903)
- Add module for evaluating the mean of variables on the domain (!894)
- Add module for counting the total number of DOF (!948)
- Fixed wss module for compressible flows (!958)
- Made Sutherland's law non-dimensional (!972)
- Add module for removing fields from .fld files (!978)
- Fixed nparts option in FieldConvert and automated Info.xml generation (!995)
- Added if statement to fix case of 1D/2D manifold interpolation in 1D/2D space,
  added check on dimensions for interpolation, fixed seg interp (!999)
- Fixed scaling for compressed xml, fixed error printout for mesh only (!1040)
- Add field conversion from Halfmode to SingleMode (!1032)
- Fix double precision output in .dat format (!1059)
- Add phase sampling feature in FilterFieldConvert (!1068)

**IncNavierStokesSolver**
- Replace steady-state check based on difference of norms by check based on
  norm of the difference, to be consistent with the compressible solver (!832)
- Updated SVV to allow for the DGKernel extension (!851)
- Pre-calculate Time invariant portion of Womersley Solution (!814)
- Fix for independent setting of SVV in Homogeneous direction (!936)
- Write flow field based on CFL threshold (!1025)
- Fix unsteady Stokes solver (!1074)

**CompressibleFlowSolver**
- Add 3D regression tests (!567)
- Introduce forcing for quasi-1D Euler simulations (!771)
- Allow performing axi-symmetric Euler and NS simulations (!771, !866)
- Add ability to use an exponential filtering for stabilization with
  seg, quad and hex elements (!771, !862)
- Fix compressible solver with NUMMODES=1 (!868)
- Introduce equations of state to account for real gas effects (!880)
- Made Sutherland's law non-dimensional (!972)
- Modified pressure outlet BCs to allow for the reference static pressure to be
  set from the VALUE fields (!981)
- hp scaling for Laplacian AV (!1013)
- Removed smooth AV (!1072)

**AcousticSolver:**
- Added two new boundary conditions to the APE system: RiemannInvariantBC
  and WhiteNoise (!782)
- Store base flow fields in a discontinuous projection (!918)
- Enabled 1D cases (!918)
- The APE system now uses u_i, c^2 and rho as base flow fields (!918)
- Added the Linearized Euler Equations (LEE) (!918)

**ADRSolver:**
- Fix forcing from file for Poisson solver (!1029)

**APESolver:**
- APESolver was replaced with AcousticSolver (!918)

**PulseWaveSolver**
- Added two new boundary conditions: AInflow and UInflow

**CardiacEPSolver**
- Converted FentonKarma model to dimensional form and added variants (!1011)

**Documentation**:
- Added an initial developer's guide (!1001)
- Updated user guide to reflect current implementation (!1051)
- Added manpages for key solvers and utilities (!1051)

**Tester**
- Fix build with boost 1.67 (!947)
- Various change to tests to decrease test time (!1053)
- Extend to support MPI tests with multiple executables (!1085)

**Packaging:**
- Add Dockerfiles and gitlab CI configuration for automatic builds (!1021,
  !1092, !1098)

v4.4.2
------
**Library**
- Fix evaluation of points (e.g. HistoryPoints, Interpolation to pts) close to
  the interface of two elements (!836)
- Fix deadlock in Hdf5 with homogeneous expansions (!858)
- Fix a few memory leaks in polylib (!863)
- Fix a crash when Interpolator is called on an empty field (!869)
- Fix petsc compile without MPI (!873)
- Fix calculation of BLPoints (!892)
- Fix deadlock in DiffusionLDG (!885)
- Fix uninitialised coefficients in DirectFull solver (!898)
- Updated PETSc to 3.7.7 (!916)
- Fix typecast to an integer which set Lz < 1 to zero when postprocess hdf5 output (!922)
- Fix program options errors on Windows in debug mode (!986)
- Fix potential clobbered output of ModArnoldi EVs when run in parallel (!983)

**IncNavierStokesSolver**
- Add a test for imaginary shift to be only used with Homogenous and SingleMode on. (!928)

**NekMesh**
- Fix missing periodic boundary meshing and boundary layer mesh adjustment
  configurations in 2D (!859)
- Fix 2D BL splitting where out-of-plane nodes would be created (!887)

**Documentation**:
- Fix sign of the viscous term in the velocity correction scheme equations in
  the user guide (!856)
- Fixed anonymous clone URL (!909)
- Add information on the limitations of Imaginary Shift for stability. (!928)

**FieldConvert**
- Allow passing input name with trailing separator (!879)
- Fix the interpcoord option  of the interppointdatatofld module (!952)

**Utilities**
- Fix VtkToPng to account for deprecated VTK API for VTK version > 8.1 (!925)

v4.4.1
------
**Library**
- Remove m_offset_elmt_id and GetOffsetElmtId which fixed problems in 2D when
  quad elements are listed before tri elements (!758)
- Remove the duplicate output of errorutil (!756)
- Fix BLAS CMake dependencies (!763)
- Fix interpolation issue with Lagrange basis functions (!768)
- Fix issue with average fields not working with different polynomial order
  fields (!776)
- Fix rounding of integer parameters (!774)
- Fix Hdf5 output in FilterFieldConvert (!781)
- Fixed extreme memory consumption of Interpolator when interpolating from pts
  to fld or between different meshes (!783)
- Fix deadlock with HDF5 input (!786)
- Fix missing entriess in LibUtilities::kPointsTypeStr (!792)
- Fix compiler warnings with CommDataType (!793)
- Fix ability to set default implementation in Collections and added an option
  to set eNoCollections in FieldConvert as default (!789)
- Fix performance issue in ProcessIsoContour in relation to memory consumption
  (!821)
- Fix performance issue with ExtractPhysToBndElmt (!796)
- Fix available classes being listed multiple times (!817)
- Fix Intel compiler warnings (!837)
- Fix overwriting and backup of chk/fld files on slow file systes (!741)
- Fix DriverAdaptive with second order IMEX (!850)
- Fixed typo in eIMEXGear part (!854)
- Added regression tests for IMEXOrder1, IMEXOrder2, IMEXOrder3, MCNAB,
  IMEXGear, CNAB, 2nd order IMEX-DIRK, 3rd order IMEX-DIRK (!854)
- Fix bug due to subtractive cancellation in polylib routines (!778)


**FieldConvert:**
- Fix issue with field ordering in the interppointdatatofld module (!754)
- Fix issue with FieldConvert when range flag used (!761)
- Fix issue when using output-points combined with noequispaced (!775)
- Fix equispacedoutput for 3DH1D with triangles (!787)

**NekMesh**:
- Fix memory consumption issue with Gmsh output (!747, !762)
- Rework meshing control so that if possible viewable meshes will be dumped
  when some part of the system fails (!756)
- Add manifold meshing option (!756)
- Fix issue with older rea input files (!765)
- Fix memory leak in variational optimiser, add small optimisations (!785)
- Check the dimensionality of the CAD system before running the 2D generator (!780)
- Fix uninitialised memory bug in Nek5000 input module (!801)

**IncNavierStokesSolver**
- Fix an initialisation issue when using an additional advective field (!779)
- Fix MovingBody boundary condition (!852)

**Utilities**
- Fix vtkToFld missing dependency which prevented compiling with VTK 7.1 (!808)

**Documentation**
- Added missing details on artificial viscosity and dealising to compressible
  flow solver user guide (!846)

**Packaging**
- Added missing package for FieldUtils library (!755)

**ADRSolver:**
- Fix UnsteadyAdvectionDiffusion with DG (!855)

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
- Fix warnings with Intel compiler (!742)

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
- Add basic gmsh cad (.geo) reader to the meshing system (!731)
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
