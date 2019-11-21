#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <Collections/Collection.h>
#include <SpatialDomains/MeshGraph.h>

// Add likwid support
#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

using namespace std;
using namespace Nektar;

MultiRegions::ExpListSharedPtr SetupExpList(
    int                                  N,
    LibUtilities::SessionReaderSharedPtr session,
    SpatialDomains::MeshGraphSharedPtr   graph,
    Collections::ImplementationType      impType,
    bool                                 useOrtho)
{
    // Accommodate use of orthogonal basis without having to specify in the mesh
    // file.
    graph->SetExpansionsToPolyOrder(N);

    if (useOrtho)
    {
        SpatialDomains::ExpansionMap expMap = graph->GetExpansions();
        for (auto &expIt : expMap)
        {
            auto exp = expIt.second;
            auto bkeyVec = exp->m_basisKeyVector;

            std::vector<LibUtilities::BasisType> newbtype;

            switch(exp->m_geomShPtr->GetShapeType())
            {
                case LibUtilities::eQuadrilateral:
                    newbtype = { LibUtilities::eOrtho_A, LibUtilities::eOrtho_A };
                    break;
                case LibUtilities::eTriangle:
                    newbtype = { LibUtilities::eOrtho_A, LibUtilities::eOrtho_B };
                    break;
                case LibUtilities::eHexahedron:
                    newbtype = { LibUtilities::eOrtho_A, LibUtilities::eOrtho_A, LibUtilities::eOrtho_A };
                    break;
                case LibUtilities::ePrism:
                    newbtype = { LibUtilities::eOrtho_A, LibUtilities::eOrtho_A, LibUtilities::eOrtho_B };
                    break;
                case LibUtilities::eTetrahedron:
                    newbtype = { LibUtilities::eOrtho_A, LibUtilities::eOrtho_B, LibUtilities::eOrtho_C };
                    break;
            }

            for (int i = 0; i < exp->m_geomShPtr->GetShapeDim(); ++i)
            {
                exp->m_basisKeyVector[i] = LibUtilities::BasisKey(
                    newbtype[i], bkeyVec[i].GetNumModes(),
                    bkeyVec[i].GetPointsKey());
            }
        }
    }

    MultiRegions::ExpListSharedPtr expList;

    switch(graph->GetMeshDimension())
    {
        case 1:
            expList = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(
                session, graph);
            break;
        case 2:
            expList = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(
                session, graph);
            break;
        case 3:
            expList = MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(
                session, graph);
            break;
        default:
            ASSERTL0(false, "Mesh must be 2D or 3D.");
    }

    expList->CreateCollections(impType);
    return expList;
}

int main(int argc, char *argv[])
{
    LibUtilities::SessionReader::RegisterCmdLineFlag(
        "data", "d", "Print in data format");
    LibUtilities::SessionReaderSharedPtr session
        = LibUtilities::SessionReader::CreateInstance(argc, argv);
    LibUtilities::CommSharedPtr comm = session->GetComm();

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("Compute");

    // Read in mesh
    SpatialDomains::MeshGraphSharedPtr graph =
        SpatialDomains::MeshGraph::Read(session);

    int Ntest, order, verbose_p, deformed_p, ortho_p;
    session->LoadParameter("Ntest", Ntest, 1);
    session->LoadParameter("order", order, 3);
    session->LoadParameter("verbose", verbose_p, 1);
    session->LoadParameter("deformed", deformed_p, 0);
    session->LoadParameter("ortho", ortho_p, 0);

    bool verbose = verbose_p != 0;
    bool deformed = deformed_p != 0;
    bool ortho = ortho_p != 0;

    bool fmt = session->DefinesCmdLineArgument("data");
    string sl = fmt ? "# " : "";

    MultiRegions::ExpListSharedPtr expList =
        SetupExpList(order + 1, session, graph, Collections::eIterPerExp, ortho);

    const int nq = expList->GetNpoints();
    const int nm = expList->GetNcoeffs();
    const int nElmt = expList->GetExpSize();
    const int nmElmt = nm / nElmt;
    const int nmMatElmt = nmElmt * nmElmt;
    const int dim = expList->GetExp(0)->GetShapeDimension();

    // Scan through list of elements and check if any are . If they are,
    // we assume they _all_ have curvature for the purposes of testing.

    std::vector<double> matData(nmMatElmt * nElmt);
    for (int i = 0; i < nElmt; ++i)
    {
        auto elmt = expList->GetExp(i);
        LocalRegions::MatrixKey matKey(
            StdRegions::eLaplacian, elmt->DetShapeType(), *elmt);

        if (verbose && comm->GetRank() == 0)
        {
            LibUtilities::PrintProgressbar(i*100/nElmt, 100, "Generating matrices");
        }

        const double *matrixRaw = elmt->GetLocStaticCondMatrix(matKey)->GetBlock(0,0)->GetRawPtr();
        double *matDataElmt = & matData[i * nmMatElmt];
        std::memcpy(matDataElmt, matrixRaw, nmMatElmt);
    }

    if (verbose && comm->GetRank() == 0)
    {
        std::cout << std::endl;
    }

    Array<OneD, NekDouble> in(nm), out(nm);

    LibUtilities::Timer t;
    t.Start();
    for (int n = 0; n < Ntest; ++n)
    {
        NekDouble *inptr = &in[0], *outptr = &out[0];
        for (int i = 0; i < nElmt; ++i)
        {
            Blas::Dgemv('N', nmElmt, nmElmt, 1.0, &matData[i*nmMatElmt], nmElmt, inptr, 1, 0.0, outptr, 1);
            inptr += nmElmt;
            outptr += nmElmt;
        }
    }
    t.Stop();

    // Elapsed time
    NekDouble elapsed = t.TimePerTest(Ntest);
    comm->AllReduce(elapsed, LibUtilities::ReduceSum);
    elapsed /= comm->GetSize();

    double dofs = nm;
    comm->AllReduce(dofs, LibUtilities::ReduceSum);
    dofs /= elapsed;

    if (comm->GetRank() == 0)
    {
        if(verbose)
        {
            std::cout << "throughput = " << dofs << " dof/s" << std::endl;
        }
        else
        {
            std::cout << dofs << std::endl;
        }

    }

    LIKWID_MARKER_CLOSE;

    comm->Finalise();

    return 0;
}
