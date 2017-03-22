#include <cstdio>
#include <cstdlib>
#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>

#include <StdRegions/StdRegions.hpp>

using namespace Nektar;
using namespace StdRegions;
using namespace std;

NekDouble Tri_sol(NekDouble x, NekDouble y, int order1, int order2);

int main(int argc, char *argv[])
{
    int i;
    int order, nq1, nq2;
    LibUtilities::PointsType Qtype1, Qtype2;
    LibUtilities::BasisType btype1, btype2;
    StdExpansion *E;

    if (argc != 4)
    {
        fprintf(stderr, "Usage: EquiToCoeff2D   order  nq1 nq2  \n");
        exit(1);
    }

    btype1 = (LibUtilities::BasisType)LibUtilities::eModified_A;
    btype2 = (LibUtilities::BasisType)LibUtilities::eModified_B;

    order = atoi(argv[1]);
    nq1   = atoi(argv[2]);
    nq2   = atoi(argv[3]);

    Qtype1 = LibUtilities::eGaussLobattoLegendre;
    Qtype2 = LibUtilities::eGaussRadauMAlpha1Beta0;

    //-----------------------------------------------
    // Define a segment expansion based on basis definition
    const LibUtilities::PointsKey Pkey1(nq1, Qtype1);
    const LibUtilities::PointsKey Pkey2(nq2, Qtype2);
    const LibUtilities::BasisKey b1(btype1, order, Pkey1);
    const LibUtilities::BasisKey b2(btype2, order, Pkey2);
    E = new StdTriExp(b1, b2);

    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(nq1 * nq2);
    Array<OneD, NekDouble> y = Array<OneD, NekDouble>(nq1 * nq2);
    E->GetCoords(x, y);

    Array<OneD, NekDouble> sol = Array<OneD, NekDouble>(nq1 * nq2);
    //----------------------------------------------
    // Define solution to be projected
    for (i = 0; i < nq1 * nq2; ++i)
    {
        sol[i] = Tri_sol(x[i], y[i], order, order);
    }
    //----------------------------------------------

    Array<OneD, NekDouble> esol(order * (order + 1) / 2);
    Array<OneD, NekDouble> coeffs(order * (order + 1) / 2);

    //---------------------------------------------
    // Interpolate to equispaced points
    E->PhysInterpToSimplexEquiSpaced(sol, esol, order);
    //---------------------------------------------

    //--------------------------------------------
    // Project from equispaced points back to coeffs
    E->EquiSpacedToCoeffs(esol, coeffs);
    //--------------------------------------------

    Array<OneD, NekDouble> phys = Array<OneD, NekDouble>(nq1 * nq2);
    //--------------------------------------------
    // Project from equispaced points back to coeffs
    E->BwdTrans(coeffs, phys);
    //--------------------------------------------

    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << E->Linf(phys, sol) << endl;
    cout << "L 2 error:        " << E->L2(phys, sol) << endl;
    //--------------------------------------------
}

NekDouble Tri_sol(NekDouble x, NekDouble y, int order1, int order2)
{
    int l, k;
    NekDouble sol = 0;

    for (k = 0; k < order1; ++k)
    {
        for (l = 0; l < order2 - k; ++l)
        {
            sol += pow(x, k) * pow(y, l);
        }
    }

    return sol;
}
