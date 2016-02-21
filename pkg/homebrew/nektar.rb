class Nektar < Formula
  class NektarThirdPartyDownloadStrategy < CurlDownloadStrategy
    def stage
      cp cached_location, basename_without_params
    end
  end

  desc "Nektar++ spectral/hp element framework"
  homepage "http://www.nektar.info/"
  url "http://ae-nektar.ae.ic.ac.uk/~dmoxey/nektar-4.2.0.tar.gz"
  version "4.2.0"
  sha256 "5ae78f8fae4f0f1bfab9fe94fdb5c1b356f9a8acf8e2bca1680e3c1f04529873"

  depends_on "cmake" => :build
  depends_on "boost"
  depends_on "tinyxml"
  depends_on "zlib"
  depends_on "vtk"
  depends_on :mpi => :recommended
  depends_on "arpack" => :recommended
  depends_on "fftw" => :recommended
  depends_on "petsc" => :recommended
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
    tp = File.join(buildpath, 'ThirdParty')
    Dir.mkdir(tp)
    resources.each do |r|
      r.stage do
        cp File.basename(r.url)[/[^?]+/], tp
      end
    end

    args = std_cmake_args + [ "-DNEKTAR_BUILD_DEMOS=OFF",
                              "-DNEKTAR_BUILD_TESTS=OFF",
                              "-DNEKTAR_BUILD_UNIT_TESTS=OFF",
                              "-DZLIB_ROOT=/usr/local/opt/zlib" ]

    if build.with? :mpi
      args << "-DNEKTAR_USE_MPI=ON"
    end

    if build.with? "arpack"
      args << "-DNEKTAR_USE_ARPACK=ON"
    end

    if build.with? "fftw"
      args << "-DNEKTAR_USE_FFTW=ON"
    end

    if build.with? "vtk"
      args << "-DNEKTAR_USE_VTK=ON"
    end

    if build.with? "petsc"
      petscdir = File.join(Formula["petsc"].opt_prefix, "real")
      args << "-DNEKTAR_USE_PETSC=ON"
      args << "-DPETSC_DIR=#{petscdir}"
    end

    if build.with? "scotch"
      args << "-DNEKTAR_USE_SCOTCH=ON"
    end

    mkdir "build" do
      system "cmake", *args, ".."
      system "make"
      components = %w(ThirdParty lib solvers util dev)
      components.each do |c| system "cmake", "-DCOMPONENT=#{c}", "-P", "cmake_install.cmake" end
    end
  end

  test do
    system "false"
  end
end
