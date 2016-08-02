# Nektar++ will, by default, download and build its ThirdParty
# dependencies. This download strategy has been added so that Homebrew can
# download resources and check hashes independently.
class NektarThirdPartyDownloadStrategy < CurlDownloadStrategy
  def stage
    cp cached_location, basename_without_params
  end
end

class Nektar < Formula
  desc "Nektar++ spectral/hp element framework"
  homepage "https://www.nektar.info/"
  url "http://www.nektar.info/downloads/file/nektar-4.3.2.tar.gz"
  sha256 "afa20b74ed19165c9356f5ec5e33dac7de9bf48e7f2076db6488c7c03d331dd0"
  head "https://gitlab.nektar.info/nektar/nektar.git"

  option "with-demos", "Compile Nektar++ demo executables"

  depends_on "cmake"  => :run
  depends_on "boost"
  depends_on "tinyxml"
  depends_on "homebrew/dupes/zlib"
  depends_on "vtk"    => :recommended
  depends_on :mpi     => :recommended
  depends_on "arpack" => :recommended
  depends_on "fftw"   => :recommended
  depends_on "petsc"  => :recommended
  depends_on "scotch" => :optional

  # Modified version of the METIS library
  resource "modmetis" do
    url "https://www.nektar.info/thirdparty/modmetis-5.1.0_2.tar.bz2",
      :using => NektarThirdPartyDownloadStrategy
    sha256 "c6496cb32c892b5ab1f78484040c1f8d53b48da232c7a7e085239ea5bd431392"
  end

  # Loki header library
  resource "loki" do
    url "https://www.nektar.info/thirdparty/loki-0.1.3.tar.bz2",
      :using => NektarThirdPartyDownloadStrategy
    sha256 "1f7aed37eec4afb113f60507955e9621808d4e34b0cb9a3c89c793e57888b65e"
  end

  # GsMpi library extracted from Nek5000 code
  resource "gsmpi" do
    url "https://www.nektar.info/thirdparty/gsmpi-1.2.tar.bz2",
      :using => NektarThirdPartyDownloadStrategy
    sha256 "ffca1d418cb7e4353de89ee11fa9fdaa878d0de5110b95b9e13e43761098ec8e"
  end

  def install
    # Copy third-party archives to ThirdParty directory. Nektar++ uses CMake
    # ExternalProject to unpack them.
    tp = File.join(buildpath, "ThirdParty")
    Dir.mkdir(tp)
    resources.each do |r|
      r.stage do
        cp File.basename(r.url)[/[^?]+/], tp
      end
    end

    args = std_cmake_args + ["-DNEKTAR_BUILD_TESTS=OFF",
                             "-DNEKTAR_BUILD_UNIT_TESTS=OFF",
                             "-DZLIB_ROOT=#{Formula["zlib"].opt_prefix}"]

    args << "-DNEKTAR_BUILD_DEMOS=ON"  if build.with?    "demos"
    args << "-DNEKTAR_BUILD_DEMOS=OFF" if build.without? "demos"
    args << "-DNEKTAR_USE_MPI=ON"      if build.with?    :mpi
    args << "-DNEKTAR_USE_ARPACK=ON"   if build.with?    "arpack"
    args << "-DNEKTAR_USE_FFTW=ON"     if build.with?    "fftw"
    args << "-DNEKTAR_USE_VTK=ON"      if build.with?    "vtk"
    args << "-DNEKTAR_USE_SCOTCH=ON"   if build.with?    "scotch"

    if build.with? "petsc"
      petscdir = File.join(Formula["petsc"].opt_prefix, "real")
      args << "-DNEKTAR_USE_PETSC=ON"
      args << "-DPETSC_DIR=#{petscdir}"
    end

    mkdir "build" do
      system "cmake", "..", *args
      system "make"
      components = %w[ThirdParty lib solvers util dev]
      components.each { |c| system "cmake", "-DCOMPONENT=#{c}", "-P", "cmake_install.cmake" }
    end
  end

  test do
    (testpath/"helm.cpp").write <<-EOS
#include <iostream>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph2D.h>

using namespace Nektar;
using namespace std;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);
    SpatialDomains::MeshGraphSharedPtr graph =
        SpatialDomains::MeshGraph::Read(vSession);
    MultiRegions::ContField2DSharedPtr fld =
        MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(
            vSession, graph, vSession->GetVariable(0));

    // Set up forcing
    const int nq = fld->GetNpoints(), nm = fld->GetNcoeffs();
    Array<OneD, NekDouble> xc0(nq), xc1(nq), xc2(nq);
    fld->GetCoords(xc0, xc1, xc2);
    vSession->GetFunction("Forcing", 0)->Evaluate(xc0, xc1, xc2, fld->UpdatePhys());

    // Solve Helmholtz
    StdRegions::ConstFactorMap factors;
    FlagList flags;
    factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
    Vmath::Zero(nm, fld->UpdateCoeffs(), 1);
    fld->HelmSolve(fld->GetPhys(), fld->UpdateCoeffs(), flags, factors);

    // Evaluate L^2 error
    fld->BwdTrans(fld->GetCoeffs(), fld->UpdatePhys());
    Array<OneD, NekDouble> sol(nq);
    vSession->GetFunction("ExactSolution", 0)->Evaluate(xc0, xc1, xc2, sol);
    cout << fld->L2(fld->GetPhys(), sol) << endl;
    return 0;
}
EOS

    (testpath/"CMakeLists.txt").write <<-EOS
find_package(Nektar++ REQUIRED NO_MODULE NO_DEFAULT_PATH NO_CMAKE_BUILDS_PATH NO_CMAKE_PACKAGE_REGISTRY)
include_directories(${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
add_definitions(${NEKTAR++_DEFINITIONS})
link_directories(${NEKTAR++_LIBRARY_DIRS} ${NEKTAR++_TP_LIBRARY_DIRS})
add_executable(helm helm.cpp)
target_link_libraries(helm ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
EOS

    (testpath/"input.xml").write <<-EOS
<NEKTAR>
    <GEOMETRY DIM="2" SPACE="2">
        <VERTEX COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAGD/aj0gwMjAzYAEKeCas8AjDj0AcDLNjNt4exWLG7Dy7PhlUfwh52HObCAAd2/XB1AAerEJkA</VERTEX>
        <EDGE COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAFjDhoJhw0Mw4aBljQ1MP4rDj4bGh8mHnsaO6BqeNA48PUcaLxYfZzoYnD9HOj8WHuAgBAqACe</EDGE>
        <ELEMENT>
            <Q COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAFjFCaCUoz4xBngdKsUJoNTZ4dSnNAaU40c5jRxLmgNDea+QAUwABZ</Q>
        </ELEMENT>
        <COMPOSITE>
            <C ID="0"> Q[0-3] </C>
            <C ID="1"> E[0,3,5,6,7,8,10,11] </C>
        </COMPOSITE>
        <DOMAIN> C[0] </DOMAIN>
    </GEOMETRY>
    <EXPANSIONS>
        <E COMPOSITE="C[0]" NUMMODES="9" TYPE="MODIFIED" FIELDS="u" />
    </EXPANSIONS>
    <CONDITIONS>
        <PARAMETERS> <P> Lambda = 1 </P> </PARAMETERS>
        <VARIABLES> <V ID="0"> u </V> </VARIABLES>
        <BOUNDARYREGIONS> <B ID="0"> C[1] </B> </BOUNDARYREGIONS>
        <BOUNDARYCONDITIONS>
            <REGION REF="0">
                <D VAR="u" VALUE="sin(PI*x)*sin(PI*y)" />
            </REGION>
        </BOUNDARYCONDITIONS>
        <FUNCTION NAME="Forcing">
            <E VAR="u" VALUE="-(Lambda + 2*PI*PI)*sin(PI*x)*sin(PI*y)" />
        </FUNCTION>
        <FUNCTION NAME="ExactSolution">
            <E VAR="u" VALUE="sin(PI*x)*sin(PI*y)" />
        </FUNCTION>
    </CONDITIONS>
</NEKTAR>
EOS
    nekpp_version = `#{bin}/ADRSolver --version | awk '{print $3}'`.strip
    system "cmake", "-DNektar++_DIR=#{lib}/nektar++-#{nekpp_version}/cmake/", "."
    system "make"
    assert (`./helm input.xml`.to_f < 1.0e-8)
  end
end
