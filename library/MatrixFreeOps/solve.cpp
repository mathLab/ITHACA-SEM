#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include <libxsmm.h>

#include "Operator.hpp"
#include "AVXAssembly.h"
#include "AVXUtil.hpp"

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/Timer.h>
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

bool solver_verbose = true;

void ConjugateGradient(
    const int nElmt,
    const Array<OneD, const NekDouble> &rhs,
    Array<OneD, NekDouble> &sol,
    LibUtilities::CommSharedPtr comm,
    std::shared_ptr<Helmholtz> oper,
    MultiRegions::AssemblyMapSharedPtr asmMap)
{
    // Conjugate gradient loop. Can use original unique map since this is
    // unchanged from AVX operations.
    auto &uniqueMap = asmMap->GetGlobalToUniversalMapUnique();

    const int nGlobal = asmMap->GetNumGlobalCoeffs();
    const int nBlocks = nElmt / oper->VectorWidth();
    const int nmElmt  = asmMap->GetNumLocalCoeffs() / nElmt;
    const int nm      = asmMap->GetNumLocalCoeffs();
    Array<OneD, NekDouble> tmp2(nm), tmp3(nm);

    // Zero initial solution field.
    Vmath::Zero(nGlobal, &sol[0], 1);

    // Firstly build a Jacobi preconditioner
    Array<OneD, NekDouble> tmpin(nm, 0.0), tmpout(nm), jacobi_local(nm), jacobi(nGlobal, 0.0);

    for (int i = 0; i < nmElmt; ++i)
    {
        int offset = i * oper->VectorWidth();
        for (int j = 0; j < nBlocks; ++j, offset += nmElmt * oper->VectorWidth())
        {
            for (int k = 0; k < oper->VectorWidth(); ++k)
            {
                tmpin[offset+k] = 1.0;
            }
        }

        if (solver_verbose && comm->GetRank() == 0)
            std::cout << "precon step " << i << "/" << nmElmt << std::endl;

        (*oper)(tmpin, tmpout);

        offset = i * oper->VectorWidth();
        for (int j = 0; j < nBlocks; ++j, offset += nmElmt * oper->VectorWidth())
        {
            for (int k = 0; k < oper->VectorWidth(); ++k)
            {
                jacobi_local[offset+k] = tmpout[offset+k];
                tmpin[offset+k] = 0.0;
            }
        }
    }

    // Assemble local contributions & calculate inverse
    //AVXAssembly<4> l2g(asmMap, nElmt);
    AVXAssemblyOld<4,16> l2g(asmMap, nElmt);

    l2g.Assemble(jacobi_local, jacobi);

    for (int i = 0; i < nGlobal; ++i)
    {
        jacobi[i] = 1.0 / jacobi[i];
    }

    Array<OneD, NekDouble> r(nGlobal), rnew(nGlobal), z(nGlobal), znew(nGlobal), p(nGlobal);
    for (int i = 0; i < nGlobal; ++i)
    {
        r[i] = -rhs[i];
        z[i] = jacobi[i] * r[i];
        p[i] = z[i];
    }

    bool root = comm->GetRank() == 0;

    int k = 0;
    const int maxiter = 1000;
    LibUtilities::Timer t;
    t.Start();

    while (k < maxiter)
    {
        if (solver_verbose && root) std::cout << k << std::endl;
        // Evaluate operator
        l2g.Scatter(p, tmp2);
        (*oper)(tmp2, tmp3);
        Vmath::Zero(nGlobal, &tmp2[0], 1);
        l2g.Assemble(tmp3, tmp2);

        std::vector<double> alpha(2, 0.0), beta(2, 0.0);
        for (int i = 0; i < nGlobal; ++i)
        {
            if (uniqueMap[i] == 0)
            {
                continue;
            }

            alpha[0] += r[i] * z[i];
            alpha[1] += p[i] * tmp2[i];
        }

        comm->AllReduce(alpha, LibUtilities::ReduceSum);
        NekDouble alph = alpha[0] / alpha[1];
        NekDouble res = 0.0;

        for (int i = 0; i < nGlobal; ++i)
        {
            sol[i] += alph * p[i];
            rnew[i] = r[i] - alph * tmp2[i];
            if (uniqueMap[i] != 0)
            {
                res += rnew[i] * rnew[i];
            }
        }

        comm->AllReduce(res, LibUtilities::ReduceSum);
        if (solver_verbose && comm->GetRank() == 0)
        {
            std::cout << "k = " << k << " res = " << res << std::endl;
        }

        if (res < 1e-10)
        {
            break;
        }

        for (int i = 0; i < nGlobal; ++i)
        {
            znew[i] = jacobi[i] * rnew[i];

            if (uniqueMap[i] == 0)
            {
                continue;
            }

            beta[0] += znew[i] * rnew[i];
            beta[1] += z[i] * r[i];
        }

        comm->AllReduce(beta, LibUtilities::ReduceSum);
        NekDouble bet = beta[0] / beta[1];

        for (int i = 0; i < nGlobal; ++i)
        {
            p[i] = znew[i] + bet * p[i];
            z[i] = znew[i];
            r[i] = rnew[i];
        }

        ++k;
    }
    t.Stop();

    // Calculate average elapsed time
    NekDouble elapsed = t.TimePerTest(maxiter);
    comm->AllReduce(elapsed, LibUtilities::ReduceSum);
    elapsed /= comm->GetSize();

    NekDouble dof = oper->Ndof();
    comm->AllReduce(dof, LibUtilities::ReduceSum);

    if (comm->GetRank() == 0)
    {
        if(solver_verbose)
        {
            std::cout << "time taken = " << t.TimePerTest(maxiter) << std::endl;
            std::cout << "total dof = " << dof << std::endl;
            std::cout << "dof/s = " << dof / elapsed << std::endl;

            if (k == maxiter)
            {
                std::cout << "failed to converge" << std::endl;
            }
        }
        else
        {
            std::cout << dof / elapsed << std::endl;
        }
    }

}

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
    MultiRegions::AssemblyMapSharedPtr asmMap;

    switch(graph->GetMeshDimension())
    {
        case 2:
        {
            auto tmp = MemoryManager<MultiRegions::ContField2D>
                ::AllocateSharedPtr(session, graph);
            asmMap = tmp->GetLocalToGlobalMap();
            expList = tmp;
            break;
        }
        case 3:
        {
            auto tmp = MemoryManager<MultiRegions::ContField3D>
                ::AllocateSharedPtr(session, graph);
            asmMap = tmp->GetLocalToGlobalMap();
            expList = tmp;
            break;
        }
        default:
            ASSERTL0(false, "Mesh must be 2D or 3D.");
            break;
    }

    //expList->CreateCollections(impType);
    return std::make_pair(expList, asmMap);
}

std::string get_opstring(LibUtilities::ShapeType shape, bool deformed=false)
{
    std::string op_string = "Helmholtz";

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

    if (deformed)
    {
        op_string += "_Deformed";
    }
    else
    {
        op_string += "_Regular";
    }

    return op_string + "_AVX";
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

    int Ntest, order, verbose_p, deformed_p, ortho_p;
    session->LoadParameter("Ntest", Ntest, 1);
    session->LoadParameter("order", order, 3);
    session->LoadParameter("verbose", verbose_p, 1);
    session->LoadParameter("deformed", deformed_p, 0);
    session->LoadParameter("ortho", ortho_p, 0);

    solver_verbose = verbose_p != 0;
    bool deformed = deformed_p != 0;
    bool ortho = ortho_p != 0;

    bool fmt = session->DefinesCmdLineArgument("data");
    string sl = fmt ? "# " : "";

    auto tmppair =
        SetupExpList(order + 1, session, graph, Collections::eIterPerExp, ortho);

    auto expList = tmppair.first;
    auto asmMap = tmppair.second;

    const int nq = expList->GetNpoints();
    const int nm = expList->GetNcoeffs();
    const int nElmt = expList->GetExpSize();
    const int dim = expList->GetExp(0)->GetShapeDimension();

    // Scan through list of elements and check if any are . If they are,
    // we assume they _all_ have curvature for the purposes of testing.

    int nqElmt = expList->GetExp(0)->GetTotPoints();
    int nmElmt = expList->GetExp(0)->GetNcoeffs();
    auto shapeType = expList->GetExp(0)->DetShapeType();

    bool deformed_fill = deformed;

    for (int i = 0; i < nElmt; ++i)
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
        jac = Array<OneD, NekDouble>(nElmt * nqElmt);
        df = Array<TwoD, NekDouble>(dim * dim, nElmt * nqElmt);

        for (int i = 0; i < nElmt; ++i)
        {
            auto jacElmt = expList->GetExp(i)->GetMetricInfo()->GetJac(
                expList->GetExp(i)->GetPointsKeys());
            auto dfElmt = expList->GetExp(i)->GetMetricInfo()->GetDerivFactors(
                expList->GetExp(i)->GetPointsKeys());

            if(deformed_fill){
                Vmath::Fill(nqElmt, jacElmt[0], &jac[i*nqElmt], 1);
                for (int j = 0; j < dim * dim; ++j)
                {
                    Vmath::Fill(nqElmt, dfElmt[j][0],&df[j][i*nqElmt], 1);
                }
            }
            else{
                ASSERTL0(jacElmt.num_elements() == nqElmt, "Incompatible number of elements for Jacobian");
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
        jac = Array<OneD, NekDouble>(nElmt);
        df = Array<TwoD, NekDouble>(dim * dim, nElmt);

        for (int i = 0; i < nElmt; ++i)
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
    auto op_string = get_opstring(shapeType, deformed);
    auto oper = std::dynamic_pointer_cast<Helmholtz>(
        GetOperatorFactory().CreateInstance(op_string, basis, nElmt));

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

    ////////////////////// Jacobi preconditioner
    const int nGlobal = asmMap->GetNumGlobalCoeffs();

    ////////////////////// Set RHS forcing function
    Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq, 0.0), tmp(nq), tmp2(nm), tmp3(nm), rhs(nm, 0.0);
    expList->GetCoords(xc, yc, zc);

    for (int i = 0; i < nq; ++i)
    {
        tmp[i] = sin(xc[i]) * cos(yc[i]) * cos(zc[i]);
    }

    // Inner product and assemble data for global RHS
    expList->IProductWRTBase(tmp, tmp2);
    asmMap->Assemble(tmp2, rhs);

    Array<OneD, NekDouble> sol(nGlobal, 0.0);

    ConjugateGradient(nElmt, rhs, sol, comm, oper, asmMap);

    LIKWID_MARKER_CLOSE;
    comm->Finalise();

    return 0;
}
