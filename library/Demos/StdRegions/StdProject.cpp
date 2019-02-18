#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <StdRegions/StdPointExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdNodalTetExp.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

namespace po = boost::program_options;

NekDouble Point_sol(NekDouble z, int order1, BasisType btype1);

NekDouble Seg_sol(NekDouble z, int order1, BasisType btype1);

NekDouble Tri_sol(NekDouble x, NekDouble y, int order1, int order2);

NekDouble Quad_sol(NekDouble x, NekDouble y, int order1, int order2,
                   BasisType btype1, BasisType btype2);

NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2,
                  int order3);

NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z, int order1,
                    int order2, int order3);

NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z, int order1,
                  int order2, int order3);

NekDouble Seg_Dsol(NekDouble z, int order1, BasisType btype1);

NekDouble Tri_Dsol(NekDouble x, NekDouble y, int order1, int order2);

NekDouble Quad_Dsol(NekDouble x, NekDouble y, int order1, int order2,
                    BasisType btype1, BasisType btype2);

NekDouble Tet_Dsol(NekDouble x, NekDouble y, NekDouble z,
                   int order1, int order2, int order3);

NekDouble Prism_Dsol(NekDouble x, NekDouble y, NekDouble z,
                     int order1, int order2, int order3);

NekDouble Hex_Dsol(NekDouble x, NekDouble y, NekDouble z, int order1,
                   int order2, int order3);

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
static double pow_loc(const double val, const int i)
{
    return (i < 0) ? 1.0 : pow(val, i);
}

int main(int argc, char *argv[])
{
    string shape, ntype;
    vector<string> basis(3, "NoBasisType");
    vector<int> order(3, 0), points(3, 0);
    po::options_description desc("Available options");
    desc.add_options()("help,h",
                       "Produce this help message and list "
                       "basis and shape types.")
            ("nodal,n", po::value<string>(&ntype),
             "Optional nodal type, autofills shape and basis choices.")
            ("shape,s", po::value<string>(&shape),
             "Region shape to project function on.")
            ("basis,b", po::value<vector<string>>(&basis)->multitoken(),
             "Basis type, separate by spaces for higher dimensions.")
            ("order,o", po::value<vector<int>>(&order)->multitoken()
                     ->required(),
             "Order of basis sets, separate by spaces for higher dimensions.")
            ("points,p", po::value<vector<int>>(&points)->multitoken()
                     ->required(),
             "Number of quadrature points, separate by spaces for "
             "higher dimensions.")
            ("diff,d",
             "Project derivative.");

    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help"))
        {
            cout << desc << endl << "All valid selections for nodal "
                                    "type, -n [ --nodal-type ], are:" << endl;
            for (int i = 22; i < SIZE_PointsType; ++i) //starts at nodal points
            {
                cout << kPointsTypeStr[i] << endl;
            };
            cout << endl << "All valid selections for shape "
                            "type, -s [ --shape ], are:" << endl;
            for (int i = 1; i < SIZE_ShapeType; ++i)
            {
                cout << ShapeTypeMap[i] << endl;
            };
            cout << endl << "All valid selections for basis type, -b "
                            "[ --basis ], are:" << endl;
            for (int i = 1; i < SIZE_BasisType; ++i)
            {
                cout << BasisTypeMap[i] << endl;
            };
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        cerr << desc;
        return 0;
    }

    //Convert string input argument to nodal type
    PointsType nodaltype = eNoPointsType;
    ShapeType stype = eNoShapeType;
    vector<BasisType> btype(3, eNoBasisType);
    if (vm.count("nodal"))
    {
        for (int i = 22; i < SIZE_PointsType; ++i) // starts at nodal points
        {
            if (boost::iequals(kPointsTypeStr[i], ntype))
            {
                nodaltype = (PointsType) i;
                break;
            }
            ASSERTL0(i != SIZE_PointsType - 1, ("The nodal type '" + ntype +
                                                "' does not exist"))
        }
        switch (nodaltype)
        {
            case eNodalTriElec:
            case eNodalTriFekete:
            case eNodalTriSPI:
            case eNodalTriEvenlySpaced:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_B;
                stype = eTriangle;
                break;
            case eNodalQuadElec:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_B;
                stype = eQuadrilateral;
                break;
            case eNodalTetElec:
            case eNodalTetSPI:
            case eNodalTetEvenlySpaced:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_B;
                btype[2] = eOrtho_C;
                stype = eTetrahedron;
                break;
            case eNodalPrismElec:
            case eNodalPrismSPI:
            case eNodalPrismEvenlySpaced:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_A;
                btype[2] = eOrtho_B;
                stype = ePrism;
                break;
            case eNodalHexElec:
                btype[0] = eOrtho_A;
                btype[1] = eOrtho_A;
                btype[2] = eOrtho_A;
                stype = eHexahedron;
                break;
            default:
                ASSERTL0(!nodaltype, ("The nodal type '" + ntype +
                                      "' is invalid."));
                break;
        }
    }

    //Convert string input argument to shape type
    if (stype == eNoShapeType)
    {
        for (int i = 1; i < SIZE_ShapeType; ++i)
        {
            if (boost::iequals(ShapeTypeMap[i], shape))
            {
                stype = (ShapeType) i;
                break;
            }
            ASSERTL0(i != SIZE_ShapeType - 1, ("The shape type '" + shape +
                                               "' does not exist"))
        }
    }

    if (stype == ePoint && vm.count("diff"))
    {
        NEKERROR(ErrorUtil::efatal,
                 "It is not possible to run the diff version for shape: point")
    }

    //Check arguments supplied equals dimension
    int dimension = ShapeTypeDimMap[stype];
    ASSERTL0(order.size() == dimension,
             "Number of orders supplied should match shape dimension");
    ASSERTL0(points.size() == dimension,
             "Number of points supplied should match shape dimension");

    if (!vm.count("nodal"))
    {
        ASSERTL0(basis.size() == dimension,
                 "Number of bases supplied should match shape dimension");
        //Convert string input argument to basis types
        for (int i = 0; i < dimension; ++i)
        {
            for (int j = 1; j < SIZE_BasisType; ++j)
            {
                if (boost::iequals(BasisTypeMap[j], basis[i]))
                {
                    btype[i] = (BasisType) j;
                    break;
                }
                ASSERTL0(j != SIZE_BasisType - 1, ("The basis type '" + basis[i]
                                                   + "' does not exist"))
            }
        }
    }

    //check basis selection is permitted for chosen shape
    map<ShapeType, vector<vector<BasisType>>> allowableBasis;
    allowableBasis[ePoint] = {
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
                    eFourierHalfModeRe, eFourierHalfModeIm}
    };
    allowableBasis[eSegment] = {
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
                    eFourierHalfModeRe, eFourierHalfModeIm}
    };
    allowableBasis[eTriangle] = {
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[eQuadrilateral] = {
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
                    eFourierHalfModeRe, eFourierHalfModeIm},
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
                    eFourierHalfModeRe, eFourierHalfModeIm}
    };
    allowableBasis[eTetrahedron] = {
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_C, eModified_C, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[ePyramid] = {
            {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
            {eOrthoPyr_C, eModifiedPyr_C, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[ePrism] = {
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[eHexahedron] = {
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial},
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial},
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
                    eLegendre, eChebyshev, eMonomial}
    };

    for (int i = 0; i < dimension; ++i)
    {
        for (int j = 0; j < (allowableBasis[stype][i].size()); ++j)
        {
            if (allowableBasis[stype][i][j] == btype[i])
            {
                break;
            }
            ASSERTL0(j != allowableBasis[stype][i].size() - 1,
                     ("The basis type '" + (string) BasisTypeMap[btype[i]] +
                      "' is invalid for basis argument " + to_string(i + 1) +
                      " for shape '" + ShapeTypeMap[stype] + "'."))
        }
    }

    //Declaration of other variables needed
    vector<PointsType> pointstype(3, eNoPointsType);
    StdExpansion *E = nullptr;

    //Assign points type according to basis type selection
    for (int i = 0; i < dimension; ++i)
    {
        if (btype[i] == eFourier)
        {
            pointstype[i] = eFourierEvenlySpaced;
        }
        else if (btype[i] == eFourierSingleMode ||
                 btype[i] == eFourierHalfModeRe ||
                 btype[i] == eFourierHalfModeIm)
        {
            pointstype[i] = eFourierSingleModeSpaced;
        }
        else
        {
            if (i == 1 && (stype == eTriangle || stype == eTetrahedron))
            {
                pointstype[i] = eGaussRadauMAlpha1Beta0;
            }
            else if (i == 2 && (stype == eTetrahedron || stype == ePyramid))
            {
                pointstype[i] = eGaussRadauMAlpha2Beta0;
            }
            else if (i == 2 && stype == ePrism)
            {
                pointstype[i] = eGaussRadauMAlpha1Beta0;
            }
            else
            {
                pointstype[i] = eGaussLobattoLegendre;
            }
        }
    }

    vector<PointsKey> pkey;
    vector<BasisKey> bkey;
    for (int i = 0; i < dimension; ++i)
    {
        pkey.push_back(PointsKey(points[i], pointstype[i]));
        bkey.push_back(BasisKey(btype[i], order[i], pkey[i]));
    }

    switch (stype)
    {
        case ePoint:
        {
            E = new StdPointExp(bkey[0]);
            break;
        }
        case eSegment:
        {
            E = new StdSegExp(bkey[0]);
            break;
        }
        case eTriangle:
        {
            E = nodaltype != eNoPointsType ? new StdNodalTriExp(bkey[0],
                                                                bkey[1],
                                                                nodaltype)
                                           : new StdTriExp(bkey[0], bkey[1]);
            break;
        }
        case eQuadrilateral:
        {
            E = new StdQuadExp(bkey[0], bkey[1]);
            break;
        }
        case eTetrahedron:
        {
            E = nodaltype != eNoPointsType ? new StdNodalTetExp(bkey[0],
                                                                bkey[1],
                                                                bkey[2],
                                                                nodaltype)
                                           : new StdTetExp(bkey[0], bkey[1],
                                                           bkey[2]);
            break;
        }
        case ePyramid:
        {
            E = new StdPyrExp(bkey[0], bkey[1], bkey[2]);
            break;
        }
        case ePrism:
        {
            E = nodaltype != eNoPointsType ? new StdNodalPrismExp(bkey[0],
                                                                  bkey[1],
                                                                  bkey[2],
                                                                  nodaltype)
                                           : new StdPrismExp(bkey[0], bkey[1],
                                                             bkey[2]);
            break;
        }
        case eHexahedron:
        {
            E = new StdHexExp(bkey[0], bkey[1], bkey[2]);
            break;
        }
        default:
            break;
    }

    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> y = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> z = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> dx = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> dy = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> dz = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> sol = Array<OneD, NekDouble>(
            (unsigned) E->GetTotPoints());

    switch (dimension)
    {
        case 0:
        case 1:
        {
            E->GetCoords(z);
            break;
        }

        case 2:
        {
            E->GetCoords(x, y);
            break;
        }

        case 3:
        {
            E->GetCoords(x, y, z);
            break;
        }
        default:
            break;
    }

    switch (stype)
    {
        case ePoint:
        {
            for (int i = 0; i < E->GetTotPoints(); ++i)
            {
                sol[i] = Point_sol(z[i], order[0], btype[0]);
            }
            break;
        }

        case eSegment:
        {
            for (int i = 0; i < E->GetTotPoints(); ++i)
            {
                sol[i] = Seg_sol(z[i], order[0], btype[0]);
            }
            break;
        }

        case eTriangle:
        {
            for (int i = 0; i < E->GetTotPoints(); ++i)
            {
                sol[i] = Tri_sol(x[i], y[i], order[0], order[1]);
            }
            break;
        }

        case eQuadrilateral:
        {
            for (int i = 0; i < E->GetTotPoints(); ++i)
            {
                sol[i] = Quad_sol(x[i], y[i], order[0], order[1], btype[0],
                                  btype[1]);
            }
            break;
        }

        case eTetrahedron:
        case ePyramid:
        {
            for (int i = 0; i < E->GetTotPoints(); ++i)
            {
                sol[i] = Tet_sol(x[i], y[i], z[i], order[0], order[1],
                                 order[2]);
            }
            break;
        }

        case ePrism:
        {
            for (int i = 0; i < E->GetTotPoints(); ++i)
            {
                sol[i] = Prism_sol(x[i], y[i], z[i], order[0], order[1],
                                   order[2]);
            }
            break;
        }

        case eHexahedron:
        {
            for (int i = 0; i < E->GetTotPoints(); ++i)
            {
                sol[i] = Hex_sol(x[i], y[i], z[i], order[0], order[1],
                                 order[2]);
            }
            break;
        }
        default:
            break;
    }

    Array<OneD, NekDouble> phys((unsigned) E->GetTotPoints());
    Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());

    if (vm.count("diff"))
    {
        if (dimension == 1)
        {
            E->PhysDeriv(sol, sol);
        }
        else if (dimension == 2)
        {
            E->PhysDeriv(sol, dx, dy);
            Vmath::Vadd(E->GetTotPoints(), dx, 1, dy, 1, sol, 1);
        }
        else if (dimension == 3)
        {
            E->PhysDeriv(sol, dx, dy, dz);
            Vmath::Vadd(E->GetTotPoints(), dx, 1, dy, 1, sol, 1);
            Vmath::Vadd(E->GetTotPoints(), dz, 1, sol, 1, sol, 1);
        }
    }

    //Project onto expansion
    E->FwdTrans(sol, coeffs);

    //Backward transform solution to get projected values
    E->BwdTrans(coeffs, phys);

    if (vm.count("diff"))
    {
        switch (stype)
        {
            case eSegment:
            {
                for (int i = 0; i < E->GetTotPoints(); ++i)
                {
                    sol[i] = Seg_Dsol(z[i], order[0], btype[0]);
                }
                break;
            }
            case eTriangle:
            {
                for (int i = 0; i < E->GetTotPoints(); ++i)
                {
                    sol[i] = Tri_Dsol(x[i], y[i], order[0], order[1]);
                }
                break;
            }
            case eQuadrilateral:
            {
                for (int i = 0; i < E->GetTotPoints(); ++i)
                {
                    sol[i] = Quad_Dsol(x[i], y[i], order[0], order[1], btype[0],
                                       btype[1]);
                }
                break;
            }
            case ePyramid:
            case eTetrahedron:
            {
                for (int i = 0; i < E->GetTotPoints(); ++i)
                {
                    sol[i] = Tet_Dsol(x[i], y[i], z[i], order[0], order[1],
                                      order[2]);
                }
                break;
            }

            case ePrism:
            {
                for (int i = 0; i < E->GetTotPoints(); ++i)
                {
                    sol[i] = Prism_Dsol(x[i], y[i], z[i], order[0], order[1],
                                        order[2]);
                }
                break;
            }

            case eHexahedron:
            {
                for (int i = 0; i < E->GetTotPoints(); ++i)
                {
                    sol[i] = Hex_Dsol(x[i], y[i], z[i], order[0], order[1],
                                      order[2]);
                }
                break;
            }
            default:
                break;
        }
    }


    //Print selected settings
    cout << (vm.count("diff") ? "Standard project of differential for "
                              : "Standard project for ");
    cout << dimension << "D " << ShapeTypeMap[stype] << " with ";
    for (int i = 0; i < dimension; ++i)
    {
        cout << BasisTypeMap[btype[i]] << " ";
    }
    cout << "of order ";
    for (int i = 0; i < dimension; ++i)
    {
        cout << order[i] << " ";
    }
    cout << "with # points ";
    for (int i = 0; i < dimension; ++i)
    {
        cout << points[i] << " ";
    }
    cout << "\b." << endl;

    //Calculate L_inf & L_2 error
    cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
    if (stype != ePoint)
    {
        cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }

    if (!vm.count("diff"))
    {
        //Evaluate solution at x = y = 0 and print error
        Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
        t[0] = -0.5;
        t[1] = -0.25;
        t[2] = -0.3;

        switch (stype)
        {
            case ePoint:
                sol[0] = Point_sol(t[0], order[0], btype[0]);
                break;
            case eSegment:
                sol[0] = Seg_sol(t[0], order[0], btype[0]);
                break;
            case eTriangle:
                sol[0] = Tri_sol(t[0], t[1], order[0], order[1]);
                break;
            case eQuadrilateral:
                sol[0] = Quad_sol(t[0], t[1], order[0], order[1], btype[0],
                                  btype[1]);
                break;
            case eTetrahedron:
            case ePyramid:
                sol[0] = Tet_sol(t[0], t[1], t[2], order[0], order[1],
                                 order[2]);
                break;
            case ePrism:
                sol[0] = Prism_sol(t[0], t[1], t[2], order[0], order[1],
                                   order[2]);
                break;
            case eHexahedron:
                sol[0] = Hex_sol(t[0], t[1], t[2], order[0], order[1],
                                 order[2]);
                break;
            default:
                break;
        }

        NekDouble nsol = E->PhysEvaluate(t, phys);

        cout << "Error at x = (";
        for (int i = 0; i < dimension; ++i)
        {
            cout << t[i] << ", ";
        }
        cout << "\b\b): " << nsol - sol[0] << endl;
    }

    return 0;
}

NekDouble Point_sol(NekDouble z, int order1, BasisType btype1)
{
    NekDouble sol = 0.0;
    if (btype1 == eFourier)
    {
        for (int l = 0; l < order1 / 2 - 1; ++l)
        {
            sol += sin(l * M_PI * z) + cos(l * M_PI * z);
        }
    }
    else if (btype1 == eFourierSingleMode)
    {
        sol += 0.45 * sin(M_PI * z) + 0.25 * cos(M_PI * z);
    }
    else if (btype1 == eFourierHalfModeRe)
    {
        sol += 0.45 * cos(M_PI * z);
    }
    else if (btype1 == eFourierHalfModeIm)
    {
        sol += 0.25 * sin(M_PI * z);
    }
    else
    {
        for (int l = 0; l < order1; ++l)
        {
            sol += pow_loc(z, l);
        }
    }

    return sol;
}

NekDouble Seg_sol(NekDouble z, int order1, BasisType btype1)
{
    NekDouble sol = 0.0;
    if (btype1 == eFourier)
    {
        for (int l = 0; l < order1 / 2 - 1; ++l)
        {
            sol += sin(l * M_PI * z) + cos(l * M_PI * z);
        }
    }
    else if (btype1 == eFourierSingleMode)
    {
        sol += 0.25 * sin(M_PI * z) + 0.25 * cos(M_PI * z);
    }
    else
    {
        for (int l = 0; l < order1; ++l)
        {
            sol += pow_loc(z, l);
        }
    }

    return sol;
}

NekDouble Tri_sol(NekDouble x, NekDouble y, int order1, int order2)
{
    NekDouble sol = 0.0;
    for (int k = 0; k < order1; ++k)
    {
        for (int l = 0; l < order2 - k; ++l)
        {
            sol += pow_loc(x, k) * pow_loc(y, l);
        }
    }

    return sol;
}

NekDouble
Quad_sol(NekDouble x, NekDouble y, int order1, int order2, BasisType btype1,
         BasisType btype2)
{
    NekDouble sol = 0.0;
    if (btype1 == eFourier)
    {
        if (btype2 == eFourier)
        {
            for (int k = 0; k < order1 / 2; ++k)
            {
                for (int l = 0; l < order2 / 2; ++l)
                {
                    sol += sin(k * M_PI * x) * sin(l * M_PI * y) +
                           sin(k * M_PI * x) * cos(l * M_PI * y) +
                           cos(k * M_PI * x) * sin(l * M_PI * y) +
                           cos(k * M_PI * x) * cos(l * M_PI * y);
                }
            }
        }
        else if (btype2 == eFourierSingleMode)
        {
            for (int k = 0; k < order1 / 2; ++k)
            {
                sol += sin(k * M_PI * x) * sin(M_PI * y) +
                       sin(k * M_PI * x) * cos(M_PI * y) +
                       cos(k * M_PI * x) * sin(M_PI * y) +
                       cos(k * M_PI * x) * cos(M_PI * y);
            }
        }
        else
        {
            for (int k = 0; k < order1 / 2; ++k)
            {
                for (int l = 0; l < order2; ++l)
                {
                    sol += sin(k * M_PI * x) * pow_loc(y, l) +
                           cos(k * M_PI * x) * pow_loc(y, l);
                }
            }
        }
    }
    else if (btype1 == eFourierSingleMode)
    {
        if (btype2 == eFourier)
        {
            for (int l = 0; l < order2 / 2; ++l)
            {
                sol += sin(M_PI * x) * sin(l * M_PI * y) +
                       sin(M_PI * x) * cos(l * M_PI * y) +
                       cos(M_PI * x) * sin(l * M_PI * y) +
                       cos(M_PI * x) * cos(l * M_PI * y);
            }

        }
        else if (btype2 == eFourierSingleMode)
        {
            sol += sin(M_PI * x) * sin(M_PI * y) +
                   sin(M_PI * x) * cos(M_PI * y) +
                   cos(M_PI * x) * sin(M_PI * y) +
                   cos(M_PI * x) * cos(M_PI * y);
        }
        else
        {
            for (int l = 0; l < order2; ++l)
            {
                sol += sin(M_PI * x) * pow_loc(y, l) +
                       cos(M_PI * x) * pow_loc(y, l);
            }
        }
    }
    else
    {
        if (btype2 == eFourier)
        {
            for (int k = 0; k < order1; ++k)
            {
                for (int l = 0; l < order2 / 2; ++l)
                {
                    sol += sin(l * M_PI * y) * pow_loc(x, k) +
                           cos(l * M_PI * y) * pow_loc(x, k);
                }
            }
        }
        else if (btype2 == eFourierSingleMode)
        {
            for (int k = 0; k < order1; ++k)
            {
                sol += sin(M_PI * y) * pow_loc(x, k) +
                       cos(M_PI * y) * pow_loc(x, k);
            }
        }
        else
        {
            for (int k = 0; k < order1; ++k)
            {
                for (int l = 0; l < order2; ++l)
                {
                    sol += pow_loc(x, k) * pow_loc(y, l);
                }
            }
        }
    }
    return sol;
}

NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2,
                  int order3)
{
    NekDouble sol = 0.0;
    for (int k = 0; k < order1; ++k)
    {
        for (int l = 0; l < order2 - k; ++l)
        {
            for (int m = 0; m < order3 - k - l; ++m)
            {
                sol += pow_loc(x, k) * pow_loc(y, l) * pow_loc(z, m);
            }
        }
    }

    return sol;
}

NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3)
{
    NekDouble sol = 0;
    for (int k = 0; k < order1; ++k)
    {
        for (int l = 0; l < order2; ++l)
        {
            for (int m = 0; m < order3 - k; ++m)
            {
                sol += pow_loc(x, k) * pow_loc(y, l) * pow_loc(z, m);
            }
        }
    }

    return sol;
}

NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3)
{
    NekDouble sol = 0.0;
    for (int k = 0; k < order1; ++k)
    {
        for (int l = 0; l < order2; ++l)
        {
            for (int m = 0; m < order3; ++m)
            {
                sol += pow_loc(x, k) * pow_loc(y, l) * pow_loc(z, m);
            }
        }
    }

    return sol;
}

NekDouble Seg_Dsol(NekDouble z, int order1, BasisType btype1)
{
    NekDouble sol = 0.0;
    if (btype1 != eFourier)
    {
        for (int k = 0; k < order1; ++k)
        {
            sol += k * pow_loc(z, k - 1);
        }
    }
    else
    {
        for (int k = 0; k < order1 / 2 - 1; ++k)
        {
            sol += k * M_PI * (cos(k * M_PI * z) - sin(k * M_PI * z));
        }
    }
    return sol;
}

NekDouble Tri_Dsol(NekDouble x, NekDouble y, int order1, int order2)
{
    int l, k;
    NekDouble sol = 0;

    for (k = 0; k < order1; ++k)
    {
        for (l = 0; l < order2 - k; ++l)
        {
            sol += k * pow_loc(x, k - 1) * pow_loc(y, l) +
                   l * pow_loc(x, k) * pow_loc(y, l - 1);
        }
    }

    return sol;
}

NekDouble Quad_Dsol(NekDouble x, NekDouble y, int order1, int order2,
                    LibUtilities::BasisType btype1,
                    LibUtilities::BasisType btype2)
{

    int k, l;
    NekDouble sol = 0;

    if (btype1 != LibUtilities::eFourier)
    {
        if (btype2 != LibUtilities::eFourier)
        {
            for (k = 0; k < order1; ++k)
            {
                for (l = 0; l < order2; ++l)
                {
                    sol += k * pow_loc(x, k - 1) * pow_loc(y, l)
                           + l * pow_loc(x, k) * pow_loc(y, l - 1);
                }
            }
        }
        else
        {
            for (k = 0; k < order1; ++k)
            {
                for (l = 0; l < order2 / 2; ++l)
                {
                    sol += k * pow_loc(x, k - 1) * sin(M_PI * l * y)
                           + M_PI * l * pow_loc(x, k) * cos(M_PI * l * y) +
                           +k * pow_loc(x, k - 1) * cos(M_PI * l * y)
                           - M_PI * l * pow_loc(x, k) * sin(M_PI * l * y);
                }
            }
        }
    }
    else
    {
        if (btype2 != LibUtilities::eFourier)
        {
            for (k = 0; k < order1 / 2; ++k)
            {
                for (l = 0; l < order2; ++l)
                {
                    sol += M_PI * k * cos(M_PI * k * x) * pow_loc(y, l)
                           + l * sin(M_PI * k * x) * pow_loc(y, l - 1) +
                           -M_PI * k * sin(M_PI * k * x) * pow_loc(y, l)
                           + l * sin(M_PI * k * x) * pow_loc(y, l - 1);
                }
            }
        }
        else
        {
            for (k = 0; k < order1 / 2; ++k)
            {
                for (l = 0; l < order2 / 2; ++l)
                {
                    sol += M_PI * k * cos(M_PI * k * x) * sin(M_PI * l * y)
                           + M_PI * l * sin(M_PI * k * x) * cos(M_PI * l * y)
                           + M_PI * k * cos(M_PI * k * x) * cos(M_PI * l * y)
                           - M_PI * l * sin(M_PI * k * x) * sin(M_PI * l * y)
                           - M_PI * k * sin(M_PI * k * x) * sin(M_PI * l * y)
                           + M_PI * l * cos(M_PI * k * x) * cos(M_PI * l * y)
                           - M_PI * k * sin(M_PI * k * x) * cos(M_PI * l * y)
                           - M_PI * l * cos(M_PI * k * x) * sin(M_PI * l * y);
                }
            }
        }
    }

    return sol;
}

NekDouble Tet_Dsol(NekDouble x, NekDouble y, NekDouble z,
                   int order1, int order2, int order3)
{
    NekDouble sol = 0;
    for (int k = 0; k < order1; ++k)
    {
        for (int l = 0; l < order2 - k; ++l)
        {
            for (int m = 0; m < order3 - k - l; ++m)
            {
                sol += k * pow_loc(x, k - 1) * pow_loc(y, l) * pow_loc(z, m)
                       + pow_loc(x, k) * l * pow_loc(y, l - 1) * pow_loc(z, m)
                       + pow_loc(x, k) * pow_loc(y, l) * m * pow_loc(z, m - 1);
            }
        }
    }
    return sol;
}

NekDouble Prism_Dsol(NekDouble x, NekDouble y, NekDouble z,
                     int order1, int order2, int order3)
{
    NekDouble sol = 0;
    for (int k = 0; k < order1; ++k)
    {
        for (int l = 0; l < order2; ++l)
        {
            for (int m = 0; m < order3 - k; ++m)
            {
                sol += k * pow_loc(x, k - 1) * pow_loc(y, l) * pow_loc(z, m)
                       + pow_loc(x, k) * l * pow_loc(y, l - 1) * pow_loc(z, m)
                       + pow_loc(x, k) * pow_loc(y, l) * m * pow_loc(z, m - 1);
            }
        }
    }
    return sol;
}

NekDouble Hex_Dsol(NekDouble x, NekDouble y, NekDouble z,
                   int order1, int order2, int order3)
{
    NekDouble sol = 0.0;
    NekDouble a;

    for (int i = 0; i < order1; ++i)
    {
        for (int j = 0; j < order2; ++j)
        {
            for (int k = 0; k < order3; ++k)
            {
                a = i * pow_loc(x, i - 1) * pow_loc(y, j) * pow_loc(z, k);
                sol += a;
                a = j * pow_loc(x, i) * pow_loc(y, j - 1) * pow_loc(z, k);
                sol += a;
                a = k * pow_loc(x, i) * pow_loc(y, j) * pow_loc(z, k - 1);
                sol += a;
            }
        }
    }
    return sol;
}