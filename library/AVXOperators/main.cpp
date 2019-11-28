#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include <libxsmm.h>

#include "Operator.hpp"
#include "AVXUtil.hpp"

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
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

std::pair<MultiRegions::ExpListSharedPtr, MultiRegions::AssemblyMapSharedPtr> SetupExpList(
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
        SpatialDomains::ExpansionInfoMap expMap = graph->GetExpansions();
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
    MultiRegions::AssemblyMapSharedPtr asmMap;

    switch(graph->GetMeshDimension())
    {
        case 1:
        {
            auto contField = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(
                session, graph);
            asmMap = contField->GetLocalToGlobalMap();
            expList = contField;
            break;
        }
        case 2:
        {
            auto contField = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(
                session, graph);
            asmMap = contField->GetLocalToGlobalMap();
            expList = contField;
            break;
        }
        case 3:
        {
            auto contField = MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(
                session, graph);
            asmMap = contField->GetLocalToGlobalMap();
            expList = contField;
            break;
        }
        default:
            ASSERTL0(false, "Mesh must be 2D or 3D.");
    }

    expList->CreateCollections(impType);
    return std::make_pair(expList, asmMap);
}

NekDouble Linf(const Array<OneD, const NekDouble> &input,
               const Array<OneD, const NekDouble> &exact)
{
    ASSERTL0(input.num_elements() == exact.num_elements(), "wut");
    NekDouble maxdiff = 0.0;
    for (int i = 0; i < exact.num_elements(); ++i)
    {
        if(std::abs(input[i] - exact[i]) > 1e-12)
        {
            //std::cout << i / (64) << "| " << i % (64) << ": " << input[i] << ", " << exact[i] << std::endl;
        }
        maxdiff = std::max(maxdiff, std::abs(input[i] - exact[i]));
    }
    return maxdiff;
}

template<class OpType, class... Args>
std::tuple<double, double, double> RunTest(
    std::shared_ptr<OpType> op, LibUtilities::CommSharedPtr comm, int Ntest, Args... args)
{
    bool rank = comm->GetRank() == 0;

    // warmup
    for (int i = 0; i < 10; ++i)
    {
        (*op)(args...);
    }
    comm->Block();

    LibUtilities::Timer timer;
    timer.Start();
    LIKWID_MARKER_START("Compute");
    for (int i = 0; i < Ntest; ++i)
    {
        (*op)(args...);
        comm->Block();
    }
    LIKWID_MARKER_STOP("Compute");
    timer.Stop();

    // Elapsed time
    NekDouble elapsed = timer.TimePerTest(Ntest);

    // Flops
    NekDouble gflops = op->GFlops();
    NekDouble dofs = op->Ndof();

    comm->AllReduce(gflops, LibUtilities::ReduceSum);
    comm->AllReduce(dofs, LibUtilities::ReduceSum);

    dofs /= elapsed;
    gflops /= elapsed;

    return std::make_tuple(elapsed, gflops, dofs);
}

std::string get_opstring(LibUtilities::ShapeType shape, int test, bool deformed=false)
{
    std::string op_string;
    if(test == 0){
        op_string = "IProduct";
    }
    else if(test == 1){
        op_string = "BwdTrans";
    }
    else if(test == 2){
        op_string = "PhysDeriv";
    }
    else if (test == 3){
        op_string = "Helmholtz";
    }
    else if (test == 4){
        op_string = "HelmholtzGlobal";
    }
    else{
        ASSERTL0(false, "Invalid operator");
    }

    op_string += "_";

    if(shape == LibUtilities::eTriangle){
        op_string += "Tri";
    }
    else if(shape == LibUtilities::eQuadrilateral){
        op_string += "Quad";
    }
    else if(shape == LibUtilities::eTetrahedron){
        op_string += "Tet";
    }
    else if(shape == LibUtilities::ePrism){
        op_string += "Prism";
    }
    else if(shape == LibUtilities::eHexahedron){
        op_string += "Hex";
    }

    if (deformed && test != 1) //BwdTrans can use the regualr implementation
    {
        op_string += "_Deformed";
    }
    else
    {
        op_string += "_Regular";
    }

    return op_string + "_AVX";
}

void print_result_testing(std::tuple<double,double,double> &result, int test, double error, double error2=0.0, double error3=0.0)
{
    std::cout << "time = " << std::get<0>(result) << ", ";
    std::cout << "gflops = " << std::get<1>(result) << ", ";
    std::cout << "dof/s= " << std::get<2>(result) << std::endl;

    std::cout << "error= " << error;

    if(test == 2){
        std::cout << " ," << "error2= " << error2;
        if(error3 != 0){
            std::cout << " ," << "error3= " << error3;
        }
    }
    std::cout << std::endl;
}

//Easily transformed into a csv file.
void print_result_auto(std::tuple<double,double,double> &result){
    std::cout << std::get<1>(result) << ",";
    std::cout << std::get<2>(result) << std::endl;
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

    libxsmm_init();

    // MultiRegions::AVX_VECTOR_SIZE = 4;

    // Read in mesh
    SpatialDomains::MeshGraphSharedPtr graph =
        SpatialDomains::MeshGraph::Read(session);

    int Ntest, test, order, report_p, deformed_p, ortho_p;
    session->LoadParameter("Ntest", Ntest, 1);
    session->LoadParameter("order", order, 3);
    session->LoadParameter("test", test, 0);
    session->LoadParameter("report", report_p, 1);
    session->LoadParameter("deformed", deformed_p, 0);
    session->LoadParameter("ortho", ortho_p, 0);

    bool report = report_p != 0;
    bool deformed = deformed_p != 0;
    bool ortho = ortho_p != 0;

    bool fmt = session->DefinesCmdLineArgument("data");
    string sl = fmt ? "# " : "";

    auto pair = SetupExpList(order + 1, session, graph, Collections::eIterPerExp, ortho);

    MultiRegions::ExpListSharedPtr expList = pair.first;
    MultiRegions::AssemblyMapSharedPtr asmMap = pair.second;

    const int nElmtOrig = expList->GetExpSize();

    // Pad appropriately.
    //const int nElmt = 4 * (nElmtOrig / 4);
    const int nElmt = nElmtOrig + 4 - (nElmtOrig % 4);
    const int nq    = expList->GetExp(0)->GetTotPoints() * nElmt;
    const int nm    = expList->GetExp(0)->GetNcoeffs() * nElmt;
    const int dim   = expList->GetExp(0)->GetShapeDimension();

    // Scan through list of elements and check if any are . If they are,
    // we assume they _all_ have curvature for the purposes of testing.
    int nqElmt = expList->GetExp(0)->GetTotPoints();
    int nmElmt = expList->GetExp(0)->GetNcoeffs();
    auto shapeType = expList->GetExp(0)->DetShapeType();

    bool deformed_fill = deformed;

    for (int i = 0; i < nElmtOrig; ++i)
    {
        if (expList->GetExp(i)->GetMetricInfo()->GetGtype() ==
            SpatialDomains::eDeformed)
        {
            deformed = true;
            deformed_fill = false; //For testing without an actual deformed mesh
            break;
        }

        ASSERTL0(expList->GetExp(i)->DetShapeType() == shapeType,
                 "Shape type should be uniform within the mesh.");
        ASSERTL0(expList->GetExp(0)->GetTotPoints() == nqElmt,
                 "Number of quadrature points should be uniform within the mesh");
        ASSERTL0(expList->GetExp(0)->GetNcoeffs() == nmElmt,
                 "Number of modes should be uniform within the mesh");
    }

    // Assemble Jacobian list and derivative factors.
    Array<OneD, NekDouble> jac;
    Array<TwoD, NekDouble> df;

    if (deformed)
    {
        jac = Array<OneD, NekDouble>(nElmt * nqElmt, 0.0);
        df = Array<TwoD, NekDouble>(dim * dim, nElmt * nqElmt, 0.0);

        for (int i = 0; i < nElmtOrig; ++i)
        {
            auto jacElmt = expList->GetExp(i)->GetMetricInfo()->GetJac(
                expList->GetExp(i)->GetPointsKeys());
            auto dfElmt = expList->GetExp(i)->GetMetricInfo()->GetDerivFactors(
                expList->GetExp(i)->GetPointsKeys());

            if (deformed_fill)
            {
                Vmath::Fill(nqElmt, jacElmt[0], &jac[i*nqElmt], 1);
                for (int j = 0; j < dim * dim; ++j)
                {
                    Vmath::Fill(nqElmt, dfElmt[j][0],&df[j][i*nqElmt], 1);
                }
            }
            else
            {
                ASSERTL0(jacElmt.num_elements() == nqElmt,
                         "Incompatible number of elements for Jacobian");
                Vmath::Vcopy(nqElmt, &jacElmt[0], 1, &jac[i*nqElmt], 1);
                for (int j = 0; j < dim * dim; ++j)
                {
                    Vmath::Vcopy(nqElmt, &dfElmt[j][0], 1, &df[j][i*nqElmt], 1);
                }
            }
        }
    }
    else
    {
        jac = Array<OneD, NekDouble>(nElmt, 0.0);
        df = Array<TwoD, NekDouble>(dim * dim, nElmt, 0.0);

        for (int i = 0; i < nElmtOrig; ++i)
        {
            auto jacElmt = expList->GetExp(i)->GetMetricInfo()->GetJac(
                expList->GetExp(i)->GetPointsKeys());
            auto dfElmt = expList->GetExp(i)->GetMetricInfo()->GetDerivFactors(
                expList->GetExp(i)->GetPointsKeys());

            jac[i] = jacElmt[0];

            for(int j = 0; j < dim * dim; ++j)
            {
                df[j][i] = dfElmt[j][0];
            }
        }
    }

    // Basis vector.
    std::vector<LibUtilities::BasisSharedPtr> basis(dim);
    for(int i = 0; i < dim; i++)
    {
        basis[i] = expList->GetExp(0)->GetBasis(i);
    }

    // Generate operator string and create operator.
    auto op_string = get_opstring(shapeType, test, deformed);
    auto oper = GetOperatorFactory().CreateInstance(op_string, basis, nElmt);

    // If the operator needs the Jacobian, provide it here
    if (oper->NeedsJac())
    {
        oper->SetJac(jac);
    }

    // If the operator needs the derivative factors, provide it here
    if (oper->NeedsDF())
    {
        oper->SetDF(df);
    }

    if (test == 0) //IProduct
    {
        Array<OneD, NekDouble> in(nq, 0.0), out(nm, 0.0), outRef(nm, 0.0);

        ////////////// Create and run new implementation
        std::shared_ptr<IProduct> op = std::dynamic_pointer_cast<IProduct>(oper);
        ASSERTL0(op, "Failed to cast pointer.");

        op->Ref(expList, in, outRef);
        TransposeData(nElmt, op->VectorWidth(), in);

        auto result = RunTest(op, comm, Ntest, in, out);

        InvTransposeData(nElmt, op->VectorWidth(), out);

        if (comm->GetRank() == 0)
        {
            if(report){
                print_result_testing(result, test, Linf(out, outRef));
            }
            else{
                print_result_auto(result);
            }
        }
    }
    else if(test == 1) //BwdTrans
    {
        Array<OneD, NekDouble> in(nm, 0.0), out(nq, 0.0), outRef(nq, 0.0);

        ////////////// Create and run new implementation
        std::shared_ptr<BwdTrans> op = std::dynamic_pointer_cast<BwdTrans>(oper);
        ASSERTL0(op, "Failed to cast pointer.");

        op->Ref(expList, in, outRef);
        TransposeData(nElmt, op->VectorWidth(), in);

        auto result = RunTest(op, comm, Ntest, in, out);

        InvTransposeData(nElmt, op->VectorWidth(), out);

        if (comm->GetRank() == 0)
        {
            if(report){
                print_result_testing(result, test, Linf(out, outRef));
            }
            else{
                print_result_auto(result);
            }
        }
    }
    else if(test == 2) //PhysDeriv
    {
        Array<OneD, NekDouble> in(nq, 0.0), d0(nq, 0.0), d1(nq, 0.0), d2(nq, 0.0);
        Array<OneD, NekDouble> ref_d0(nq, 0.0), ref_d1(nq, 0.0), ref_d2(nq, 0.0);

        std::shared_ptr<PhysDeriv> op = std::dynamic_pointer_cast<PhysDeriv>(oper);
        ASSERTL0(op, "Failed to cast pointer.");

        if(dim == 2)
        {
            op->Ref(expList, in, ref_d0, ref_d1);
        }
        else
        {
            op->Ref(expList, in, ref_d0, ref_d1, ref_d2);
        }

        TransposeData(nElmt, op->VectorWidth(), in);

        std::tuple<double, double, double> result;

        if(dim == 2)
        {
            result = RunTest(op, comm, Ntest, in, d0, d1);
        }
        else
        {
            result = RunTest(op, comm, Ntest, in, d0, d1, d2);
        }

        InvTransposeData(nElmt, op->VectorWidth(), d0);
        InvTransposeData(nElmt, op->VectorWidth(), d1);
        InvTransposeData(nElmt, op->VectorWidth(), d2);

        if (comm->GetRank() == 0)
        {
            if(report){
                if(dim == 2)
                {
                    print_result_testing(result, test, Linf(d0, ref_d0), Linf(d1, ref_d1));
                }
                else{
                    print_result_testing(result, test, Linf(d0, ref_d0), Linf(d1, ref_d1), Linf(d2, ref_d2));
                }
            }
            else{
                print_result_auto(result);
            }
        }

    }
    else if(test == 3) // Helmholtz
    {
        Array<OneD, NekDouble> in(nm, 0.0), out(nm, 0.0), outRef(nm, 0.0);

        ////////////// Create and run new implementation
        std::shared_ptr<Helmholtz> op = std::dynamic_pointer_cast<Helmholtz>(oper);
        ASSERTL0(op, "Failed to cast pointer.");

        op->Ref(expList, in, outRef);
        TransposeData(nElmt, op->VectorWidth(), in);

        auto result = RunTest(op, comm, Ntest, in, out);

        InvTransposeData(nElmt, op->VectorWidth(), out);

        if (comm->GetRank() == 0)
        {
            if (report)
            {
                print_result_testing(result, test, Linf(out, outRef));
            }
            else
            {
                print_result_auto(result);
            }
        }
    }
    else if(test == 4)
    {
        oper->set_asmMap(asmMap);

        int nGlobal = asmMap->GetNumGlobalCoeffs() + 1; // potential for padding

        Array<OneD, NekDouble> in(nGlobal, 0.0), out(nGlobal, 0.0), outRef(nGlobal, 0.0);
        AVXAssembly<4> l2g(asmMap, expList, nElmt);

        ////////////// Create and run new implementation
        std::shared_ptr<HelmholtzGlobal<4>> op = std::dynamic_pointer_cast<HelmholtzGlobal<4>>(oper);
        ASSERTL0(op, "Failed to cast pointer.");

        op->Ref2(expList, in, outRef, asmMap);

        //TransposeData(nElmt, op->VectorWidth(), inLocal);

        auto result = RunTest(op, comm, Ntest, in, out, l2g);

        //InvTransposeData(nElmt, op->VectorWidth(), outLocal);

        if (comm->GetRank() == 0)
        {
            if(report){
                print_result_testing(result, test, Linf(out, outRef));
            }
            else{
                print_result_auto(result);
            }
        }
    }

    LIKWID_MARKER_CLOSE;

    comm->Finalise();

    return 0;
}
