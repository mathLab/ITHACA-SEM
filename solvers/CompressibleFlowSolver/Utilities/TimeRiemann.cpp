#include "../RiemannSolvers/RoeSolver.h"

#include <AVXOperators/AVXUtil.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/Likwid.hpp>

#include "../RiemannSolvers/RoeSolver.h"
#include "../RiemannSolvers/RoeSolverSIMD.h"

using namespace Nektar;

int main(int argc, char const *argv[])
{

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("Riemann");
    LIKWID_MARKER_REGISTER("RotationTo");
    LIKWID_MARKER_REGISTER("v_Solve");
    LIKWID_MARKER_REGISTER("RotationFrom");

    using vec_t = AVX::VecData<double, AVX::SIMD_WIDTH_SIZE>;

    size_t nEle;
    if (argc < 2)
    {
        nEle = 10;
    }
    else
    {
        nEle = std::stoi(argv[1]);
    }

    std::cout << "number of faces\t" << nEle
        << "\t(assuming 4*4 nodes per face)\n";

    // auto riemannSolver = RoeSolver();
    auto riemannSolver = RoeSolverSIMD();
    // Setting up parameters for Riemann solver
    NekDouble gamma = 1.4;
    riemannSolver.SetParam("gamma", [&gamma]()
        -> NekDouble& {return gamma;});

    size_t spaceDim = 3;
    size_t npts = 4*4*nEle; // so that avx spillover loop is engaged

    // Set up locations of velocity vector.
    Array<OneD, Array<OneD, NekDouble>> vecLocs(1);
    vecLocs[0] = Array<OneD, NekDouble>(spaceDim);
    for (size_t i = 0; i < spaceDim; ++i)
    {
        vecLocs[0][i] = 1 + i;
    }
    riemannSolver.SetAuxVec("vecLocs", [&vecLocs]()
        -> const Array<OneD, const Array<OneD, NekDouble>>&
        {return vecLocs;});

    // setup normals
    Array<OneD, Array<OneD, NekDouble>> normals(spaceDim);
    for (size_t i = 0; i < spaceDim; ++i)
    {
       normals[i] = Array<OneD, NekDouble>(npts);
    }
    riemannSolver.SetVector("N", [&normals]()
        -> const Array<OneD, const Array<OneD, NekDouble>>&
        {return normals;});

    size_t nFields = spaceDim + 2;
    Array<OneD, Array<OneD, NekDouble>> fwd(nFields), bwd(nFields),
        flx(nFields), flxRef(nFields);
    for (size_t i = 0; i < nFields; ++i)
    {
        fwd[i] = Array<OneD, NekDouble>(npts);
        bwd[i] = Array<OneD, NekDouble>(npts);
        flx[i] = Array<OneD, NekDouble>(npts);
        flxRef[i] = Array<OneD, NekDouble>(npts);
    }

    // up to now it is "boiler plate" code to set up the test
    // below are the conditions for the test

    for (size_t i = 0; i < npts; ++i)
    {
        // density
        NekDouble rho = 0.9;
        fwd[0][i] = rho;
        bwd[0][i] = rho;
        // x-momentum
        NekDouble rhou = rho*1.0;
        fwd[1][i] = rhou;
        bwd[1][i] = rhou;
        // y-momentum
        NekDouble rhov = rho*2.0;
        fwd[2][i] = rhov;
        bwd[2][i] = rhov;
        // z-momentum
        NekDouble rhow = rho*3.0;
        fwd[3][i] = rhow;
        bwd[3][i] = rhow;
        // energy
        NekDouble p = 1.0;
        NekDouble rhoe = p / (gamma - 1.0);
        NekDouble E = rhoe + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow) / rho;
        fwd[nFields-1][i] = E;
        bwd[nFields-1][i] = E;
        // set face normal along x
        normals[0][i] = 1.0;
        // Ref solution
        flxRef[0][i] = rhou;
        flxRef[1][i] = rhou*rhou/rho + p;
        flxRef[2][i] = rhou*rhov/rho;
        flxRef[3][i] = rhou*rhow/rho;
        flxRef[nFields-1][i] = (E+p)*rhou/rho;
    }

    // number of experiments
    constexpr size_t experiments = 1 << 12;

    LIKWID_MARKER_START("Riemann");
    for (size_t j = 0; j < experiments; ++j)
    {
        // time
        riemannSolver.Solve(spaceDim, fwd, bwd, flx);
    }
    LIKWID_MARKER_STOP("Riemann");
    // get likwid events
    constexpr short CPU_CLK_UNHALTED_REF_id = 2;
    int nevents{20};
    std::vector<double> events(nevents);
    double time;
    int count;
    //
    LIKWID_MARKER_GET("Riemann", &nevents, events.data(), &time, &count);
    // print out CPE
    double cpeRiemann = events[CPU_CLK_UNHALTED_REF_id]/npts/experiments;
    std::cout << "Riemann likwid CPE\t" << cpeRiemann << '\n';
    //
    LIKWID_MARKER_GET("RotationTo", &nevents, events.data(), &time, &count);
    // print out CPE
    double cpeRotationTo = events[CPU_CLK_UNHALTED_REF_id]/npts/experiments;
    std::cout << "RotationTo likwid CPE\t" << cpeRotationTo << '\n';
    //
    LIKWID_MARKER_GET("v_Solve", &nevents, events.data(), &time, &count);
    // print out CPE
    double cpev_Solve = events[CPU_CLK_UNHALTED_REF_id]/npts/experiments;
    std::cout << "v_Solve likwid CPE\t" << cpev_Solve << '\n';
    //
    LIKWID_MARKER_GET("RotationFrom", &nevents, events.data(), &time, &count);
    // print out CPE
    double cpeRotationFrom = events[CPU_CLK_UNHALTED_REF_id]/npts/experiments;
    std::cout << "RotationFrom likwid CPE\t" << cpeRotationFrom << '\n';
    // avoid opt out
    std::cout << flx[0][0] << std::endl;



LIKWID_MARKER_CLOSE;

}