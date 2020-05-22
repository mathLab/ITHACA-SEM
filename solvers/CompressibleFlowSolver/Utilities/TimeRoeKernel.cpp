#include "../RiemannSolvers/RoeSolver.h"

#include <AVXOperators/AVXUtil.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/Likwid.hpp>

using namespace Nektar;

int main(int argc, char const *argv[])
{

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("scalar");
    LIKWID_MARKER_REGISTER("vector");

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

    size_t sizeScalar = 4*4 * nEle;
    size_t nVars = 5;
    size_t spaceDim = nVars - 2;
    size_t sizeVec = sizeScalar / AVX::SIMD_WIDTH_SIZE;

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

    std::vector<AVX::AlignedVector<vec_t>>
        alignedFwd(nVars),
        alignedBwd(nVars),
        alignedFlux(nVars);

    for (size_t i = 0; i < nVars; ++i)
    {
        alignedFwd[i] = AVX::AlignedVector<vec_t>(sizeVec);
        alignedBwd[i] = AVX::AlignedVector<vec_t>(sizeVec);
        alignedFlux[i] = AVX::AlignedVector<vec_t>(sizeVec);

        for (size_t j = 0; j < sizeVec; ++j)
        {
            for (size_t k = 0; k < AVX::SIMD_WIDTH_SIZE; ++k)
            {
                alignedFwd[i][j].m_data[k] = 1.;
                alignedBwd[i][j].m_data[k] = 1.;
                // fix energy to avoid negative pressure
                if (i == nVars - 1)
                {
                    alignedFwd[i][j].m_data[k] = 10.;
                    alignedBwd[i][j].m_data[k] = 10.;
                }
                alignedFlux[i][j].m_data[k] = 0.;
            }
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
    LIKWID_MARKER_GET("scalar", &nevents, events.data(), &time, &count);
    // print out CPE
    std::cout << "likwid CPE\t"
        << events[CPU_CLK_UNHALTED_REF_id]/sizeScalar/experiments << '\n';
    // avoid opt out
    std::cout << Flux[0][0] << std::endl;


    LIKWID_MARKER_START("vector");
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
    LIKWID_MARKER_STOP("vector");
    // get likwid events
    LIKWID_MARKER_GET("vector", &nevents, events.data(), &time, &count);
    // print out CPE
    std::cout << "likwid CPE\t"
        << events[CPU_CLK_UNHALTED_REF_id]/sizeScalar/experiments << '\n';
    // avoid opt out
    std::cout << alignedFlux[0][0].m_data[0] << std::endl;

LIKWID_MARKER_CLOSE;

}