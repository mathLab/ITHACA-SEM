#include "../RiemannSolvers/RoeSolver.h"

#include <LibUtilities/SimdLib/tinysimd.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/Likwid.hpp>

using namespace Nektar;

int main(int argc, char const *argv[])
{

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("scalar");
    LIKWID_MARKER_REGISTER("vector");
    LIKWID_MARKER_REGISTER("vectorOfVector");

    using namespace tinysimd;
    using vec_t = simd<NekDouble>;

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

    size_t sizeScalar = 4*4 * nEle;
    size_t nVars = 5;
    size_t spaceDim = nVars - 2;
    size_t sizeVec = sizeScalar / vec_t::width;

    double gamma = 1.4;

    std::vector<std::vector<double>>
        Fwd(nVars),
        Bwd(nVars),
        Flux(nVars);

    for (size_t i = 0; i < nVars; ++i)
    {
        Fwd[i] = std::vector<double>(sizeScalar);
        Bwd[i] = std::vector<double>(sizeScalar);
        Flux[i] = std::vector<double>(sizeScalar);

        for (size_t j = 0; j < sizeScalar; ++j)
        {
            Fwd[i][j] = 1.;
            Bwd[i][j] = 1.;
            // fix energy to avoid negative pressure
            if (i == nVars - 1)
            {
                Fwd[i][j] = 10.;
                Bwd[i][j] = 10.;
            }
            Flux[i][j] = 0.;
        }
    }

    std::vector<std::vector<vec_t, allocator<vec_t>>>
        alignedFwd(nVars),
        alignedBwd(nVars),
        alignedFlux(nVars);

    for (size_t i = 0; i < nVars; ++i)
    {
        alignedFwd[i] = std::vector<vec_t, allocator<vec_t>>(sizeVec);
        alignedBwd[i] = std::vector<vec_t, allocator<vec_t>>(sizeVec);
        alignedFlux[i] = std::vector<vec_t, allocator<vec_t>>(sizeVec);

        for (size_t j = 0; j < sizeVec; ++j)
        {
            alignedFwd[i][j] = 1.;
            alignedBwd[i][j] = 1.;
            // fix energy to avoid negative pressure
            if (i == nVars - 1)
            {
                alignedFwd[i][j] = 10.;
                alignedBwd[i][j] = 10.;
            }
            alignedFlux[i][j] = 0.;
        }
    }

    // number of experiments
    constexpr size_t experiments = 1 << 18;

    LIKWID_MARKER_START("scalar");
    for (size_t j = 0; j < experiments; ++j)
    {
        // time scalar
        // loop
        for (size_t i = 0; i < sizeScalar; ++i)
        {
            double rhoL{}, rhouL{}, rhovL{}, rhowL{}, EL{};
            double rhoR{}, rhouR{}, rhovR{}, rhowR{}, ER{};

            // load
            rhoL  = Fwd[0][i];
            rhouL = Fwd[1][i];
            EL    = Fwd[spaceDim+1][i];
            rhoR  = Bwd[0][i];
            rhouR = Bwd[1][i];
            ER    = Bwd[spaceDim+1][i];

            if (spaceDim == 2)
            {
                rhovL = Fwd[2][i];
                rhovR = Bwd[2][i];
            }
            else if (spaceDim == 3)
            {
                rhovL = Fwd[2][i];
                rhowL = Fwd[3][i];
                rhovR = Bwd[2][i];
                rhowR = Bwd[3][i];
            }

            double rhof{}, rhouf{}, rhovf{}, rhowf{}, Ef{};

            RoeKernel(
                rhoL, rhouL, rhovL, rhowL, EL,
                rhoR, rhouR, rhovR, rhowR, ER,
                rhof, rhouf, rhovf, rhowf, Ef,
                gamma);

            // store
            Flux[0][i] = rhof;
            Flux[1][i] = rhouf;
            Flux[nVars-1][i] = Ef;
            if (spaceDim == 2)
            {
                Flux[2][i] = rhovf;
            }
            else if (spaceDim == 3)
            {
                Flux[2][i] = rhovf;
                Flux[3][i] = rhowf;
            }

        } // loop
    }
    LIKWID_MARKER_STOP("scalar");
    // get likwid events
    constexpr short CPU_CLK_UNHALTED_REF_id = 2;
    int nevents{20};
    std::vector<double> events(nevents);
    double time;
    int count;

    boost::ignore_unused(time,count);
    LIKWID_MARKER_GET("scalar", &nevents, events.data(), &time, &count);
    // print out CPE
    double cpeScalar = events[CPU_CLK_UNHALTED_REF_id]/sizeScalar/experiments;
    std::cout << "scalar likwid CPE\t" << cpeScalar << '\n';
    // avoid opt out
    std::cout << Flux[0][0] << std::endl;


    LIKWID_MARKER_START("vector");
    // time SIMD
    for (size_t j = 0; j < experiments; ++j)
    {
        // loop
        for (size_t i = 0; i < sizeScalar; i+=vec_t::width)
        {
            vec_t rhoL{}, rhouL{}, rhovL{}, rhowL{}, EL{};
            vec_t rhoR{}, rhouR{}, rhovR{}, rhowR{}, ER{};

            // load
            rhoL  .load(&(Fwd[0][i]), is_not_aligned);
            rhouL .load(&(Fwd[1][i]), is_not_aligned);
            EL    .load(&(Fwd[spaceDim+1][i]), is_not_aligned);
            rhoR  .load(&(Bwd[0][i]), is_not_aligned);
            rhouR .load(&(Bwd[1][i]), is_not_aligned);
            ER    .load(&(Bwd[spaceDim+1][i]), is_not_aligned);

            if (spaceDim == 2)
            {
                rhovL.load(&(Fwd[2][i]), is_not_aligned);
                rhovR.load(&(Bwd[2][i]), is_not_aligned);
            }
            else if (spaceDim == 3)
            {
                rhovL.load(&(Fwd[2][i]), is_not_aligned);
                rhowL.load(&(Fwd[3][i]), is_not_aligned);
                rhovR.load(&(Bwd[2][i]), is_not_aligned);
                rhowR.load(&(Bwd[3][i]), is_not_aligned);
            }

            vec_t rhof{}, rhouf{}, rhovf{}, rhowf{}, Ef{};

            RoeKernel(
                rhoL, rhouL, rhovL, rhowL, EL,
                rhoR, rhouR, rhovR, rhowR, ER,
                rhof, rhouf, rhovf, rhowf, Ef,
                gamma);

            // store
            rhof.store(&(Flux[0][i]), is_not_aligned);
            rhouf.store(&(Flux[1][i]), is_not_aligned);
            Ef.store(&(Flux[nVars-1][i]), is_not_aligned);
            if (spaceDim == 2)
            {
                rhovf.store(&(Flux[2][i]), is_not_aligned);
            }
            else if (spaceDim == 3)
            {
                rhovf.store(&(Flux[2][i]), is_not_aligned);
                rhowf.store(&(Flux[3][i]), is_not_aligned);
            }

        } // loop
    }
    LIKWID_MARKER_STOP("vector");
    // get likwid events
    LIKWID_MARKER_GET("vector", &nevents, events.data(), &time, &count);
    // print out CPE
    double cpeVector = events[CPU_CLK_UNHALTED_REF_id]/sizeScalar/experiments;
    std::cout << "vector likwid CPE\t" << cpeVector << '\t'
        << cpeScalar/cpeVector << " %\n";
    // avoid opt out
    std::cout << Flux[0][0] << std::endl;

    LIKWID_MARKER_START("vectorOfVector");
    // time SIMD
    for (size_t j = 0; j < experiments; ++j)
    {
        // loop
        for (size_t i = 0; i < sizeVec; ++i)
        {
            vec_t rhoL{}, rhouL{}, rhovL{}, rhowL{}, EL{};
            vec_t rhoR{}, rhouR{}, rhovR{}, rhowR{}, ER{};

            // load
            rhoL  = alignedFwd[0][i];
            rhouL = alignedFwd[1][i];
            EL    = alignedFwd[spaceDim+1][i];
            rhoR  = alignedBwd[0][i];
            rhouR = alignedBwd[1][i];
            ER    = alignedBwd[spaceDim+1][i];

            if (spaceDim == 2)
            {
                rhovL = alignedFwd[2][i];
                rhovR = alignedBwd[2][i];
            }
            else if (spaceDim == 3)
            {
                rhovL = alignedFwd[2][i];
                rhowL = alignedFwd[3][i];
                rhovR = alignedBwd[2][i];
                rhowR = alignedBwd[3][i];
            }

            vec_t rhof{}, rhouf{}, rhovf{}, rhowf{}, Ef{};

            RoeKernel(
                rhoL, rhouL, rhovL, rhowL, EL,
                rhoR, rhouR, rhovR, rhowR, ER,
                rhof, rhouf, rhovf, rhowf, Ef,
                gamma);

            // store
            alignedFlux[0][i] = rhof;
            alignedFlux[1][i] = rhouf;
            alignedFlux[nVars-1][i] = Ef;
            if (spaceDim == 2)
            {
                alignedFlux[2][i] = rhovf;
            }
            else if (spaceDim == 3)
            {
                alignedFlux[2][i] = rhovf;
                alignedFlux[3][i] = rhowf;
            }

        } // loop
    }
    LIKWID_MARKER_STOP("vectorOfVector");
    // get likwid events
    LIKWID_MARKER_GET("vectorOfVector", &nevents, events.data(), &time, &count);
    // print out CPE
    double cpevectorOfVector = events[CPU_CLK_UNHALTED_REF_id]/sizeScalar/experiments;
    std::cout << "vectorOfVector likwid CPE\t" << cpevectorOfVector << '\t'
        << cpeScalar/cpevectorOfVector << " %\n";
    // avoid opt out
    std::cout << alignedFlux[0][0][0] << std::endl;

LIKWID_MARKER_CLOSE;

}
