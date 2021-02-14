#include <LibUtilities/SimdLib/tinysimd.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/Likwid.hpp>

using namespace Nektar;

int main(int argc, char const *argv[])
{

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("GathrArray");
    LIKWID_MARKER_REGISTER("GathrScalar");
    LIKWID_MARKER_REGISTER("GathrSimd");


    size_t nPts;
    if (argc < 2)
    {
        nPts = 100;
    }
    else
    {
        nPts = std::stoi(argv[1]);
    }

    std::cout << "number of points\t" << nPts << '\n';


    // number of experiments
    constexpr size_t experiments = 1 << 18;
    {
        // data should be randomized
        Array<OneD, NekDouble> data{nPts*10, 0.0};
        Array<OneD, NekDouble> dataTrace{nPts, 0.0};
        Array<OneD, size_t> indexTrace{nPts, 0.0};

        // inizialize index (should be randomized)
        for (size_t i = 0; i < nPts; ++i)
        {
            indexTrace[i] = i;
        }

        LIKWID_MARKER_START("GathrArray");
        for (size_t j = 0; j < experiments; ++j)
        {
            // time
            Vmath::Gathr(nPts, data, indexTrace, dataTrace);
        }
        LIKWID_MARKER_STOP("GathrArray");
        // get likwid events
        constexpr short CPU_CLK_UNHALTED_REF_id = 2;
        int nevents{20};
        std::vector<double> events(nevents);
        //
        LIKWID_MARKER_GET("GathrArray", &nevents, events.data(), &time, &count);
        // print out CPE
        double cpeGathrArray = events[CPU_CLK_UNHALTED_REF_id]/nPts/experiments;
        std::cout << "GathrArray likwid CPE\t" << cpeGathrArray << '\n';
        std::cout << dataTrace[0] << std::endl;
    }

    {
        // data should be randomized
        Array<OneD, NekDouble> data{nPts*10, 0.0};
        Array<OneD, NekDouble> dataTrace{nPts, 0.0};
        Array<OneD, size_t> indexTrace{nPts, 0.0};

        // inizialize index (should be randomized)
        for (size_t i = 0; i < nPts; ++i)
        {
            indexTrace[i] = i;
        }

        LIKWID_MARKER_START("GathrSimd");
        for (size_t j = 0; j < experiments; ++j)
        {
            // time
            Vmath::SIMD::Gathr(nPts, data.data(), indexTrace.data(), dataTrace.data());
        }
        LIKWID_MARKER_STOP("GathrSimd");
        // get likwid events
        constexpr short CPU_CLK_UNHALTED_REF_id = 2;
        int nevents{20};
        std::vector<double> events(nevents);
        //
        LIKWID_MARKER_GET("GathrSimd", &nevents, events.data(), &time, &count);
        // print out CPE
        double cpeGathrSimd = events[CPU_CLK_UNHALTED_REF_id]/nPts/experiments;
        std::cout << "GathrSimd likwid CPE\t" << cpeGathrSimd << '\n';
        std::cout << dataTrace[0] << std::endl;
    }

    {
        // data should be randomized
        Array<OneD, NekDouble> data{nPts*10, 0.0};
        Array<OneD, NekDouble> dataTrace{nPts, 0.0};
        Array<OneD, size_t> indexTrace{nPts, 0.0};

        // inizialize index (should be randomized)
        for (size_t i = 0; i < nPts; ++i)
        {
            indexTrace[i] = i;
        }

        LIKWID_MARKER_START("GathrScalar");
        for (size_t j = 0; j < experiments; ++j)
        {
            // time
            Vmath::Gathr(nPts, data.data(), indexTrace.data(), dataTrace.data());
        }
        LIKWID_MARKER_STOP("GathrScalar");
        // get likwid events
        constexpr short CPU_CLK_UNHALTED_REF_id = 2;
        int nevents{20};
        std::vector<double> events(nevents);
        //
        LIKWID_MARKER_GET("GathrScalar", &nevents, events.data(), &time, &count);
        // print out CPE
        double cpeGathrScalar = events[CPU_CLK_UNHALTED_REF_id]/nPts/experiments;
        std::cout << "GathrScalar likwid CPE\t" << cpeGathrScalar << '\n';
        std::cout << dataTrace[0] << std::endl;
    }


LIKWID_MARKER_CLOSE;

}
