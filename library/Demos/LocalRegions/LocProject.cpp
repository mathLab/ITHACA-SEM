#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/TetExp.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;
using namespace Nektar::SpatialDomains;
using namespace Nektar::LocalRegions;

namespace po = boost::program_options;

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

GeometrySharedPtr CreateGeom(vector<double> coords, ShapeType stype);

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
static double pow_loc(const double val, const int i)
{
    return (i < 0) ? 1.0 : pow(val, i);
}

int main(int argc, char *argv[])
{
    string shape, ntype;
    vector<string> basis(3, "NoBasisType");
    vector<int> order, points;
    vector<double> coords;
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
            ("coords,c", po::value<vector<double>>(&coords)->multitoken()
                     ->required(),
             "Coordinates, separate by spaces for higher dimensions with "
             "grouping by vertex i.e. 'x1 y1 z1 x2 y2 z2 x3 y3 z3 ...")
            ("diff,d",
             "Project derivative.");

    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help"))
        {
            cout << desc;
            cout << endl << "All nodal types, -n [ --nodal ], are:" << endl;
            for (int i = 22; i < SIZE_PointsType; ++i)
            {
                cout << kPointsTypeStr[i] << endl;
            };
            cout << endl << "All shape types, -s [ --shape ], are:" << endl;
            for (int i = 2; i < SIZE_ShapeType; ++i)
            {
                cout << ShapeTypeMap[i] << endl;
            };
            cout << endl << "All basis types, -b [ --basis ], are:" << endl;
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
        for (int i = 1; i < SIZE_PointsType; ++i) // starts at nodal points
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
            default:
                ASSERTL0(!nodaltype, ("The nodal type '" + ntype +
                                      "' is invalid for LocProject."));
                break;
        }
    }

    //Convert string input argument to shape type
    if (stype == eNoShapeType)
    {
        for (int i = 2; i < SIZE_ShapeType; ++i)
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

    GeometrySharedPtr geomFull = CreateGeom(coords, stype);

    switch (stype)
    {
        case eSegment:
        {
            SegGeomSharedPtr geom = dynamic_pointer_cast<SegGeom>(geomFull);
            E = new SegExp(bkey[0], geom);
            break;
        }
        case eTriangle:
        {
            TriGeomSharedPtr geom = dynamic_pointer_cast<TriGeom>(geomFull);
            if (nodaltype == eNoPointsType)
            {
                E = new TriExp(bkey[0], bkey[1], geom);
            }
            else
            {
                E = new NodalTriExp(bkey[0], bkey[1], nodaltype, geom);
            }
            break;
        }
        case eQuadrilateral:
        {
            QuadGeomSharedPtr geom = dynamic_pointer_cast<QuadGeom>(geomFull);
            E = new QuadExp(bkey[0], bkey[1], geom);
            break;
        }
        case eTetrahedron:
        {
            TetGeomSharedPtr geom = dynamic_pointer_cast<TetGeom>(geomFull);
            E = new TetExp(bkey[0], bkey[1], bkey[2], geom);
            break;
        }
        case ePyramid:
        {
            PyrGeomSharedPtr geom = dynamic_pointer_cast<PyrGeom>(geomFull);
            E = new PyrExp(bkey[0], bkey[1], bkey[2], geom);
            break;
        }
        case ePrism:
        {
            PrismGeomSharedPtr geom = dynamic_pointer_cast<PrismGeom>(geomFull);
            E = new PrismExp(bkey[0], bkey[1], bkey[2], geom);
            break;
        }
        case eHexahedron:
        {
            HexGeomSharedPtr geom = dynamic_pointer_cast<HexGeom>(geomFull);
            E = new HexExp(bkey[0], bkey[1], bkey[2], geom);
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

    //Calculate L_inf & L_2 error
    cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
    cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;

    return 0;
}

GeometrySharedPtr CreateGeom(vector<double> coords,
                             LibUtilities::ShapeType shapeType)
{
    map<ShapeType, vector<vector<int>>> edgeDef, faceDef, volDef;

    edgeDef[eSegment] = {{0, 1}};
    edgeDef[eTriangle] = {{0, 1},
                          {1, 2},
                          {2, 0}};
    edgeDef[eQuadrilateral] = {{0, 1},
                               {1, 2},
                               {2, 3},
                               {3, 0}};
    edgeDef[eHexahedron] = {{0, 1},
                            {1, 2},
                            {3, 2},
                            {0, 3},
                            {0, 4},
                            {1, 5},
                            {2, 6},
                            {3, 7},
                            {4, 5},
                            {5, 6},
                            {7, 6},
                            {4, 7}};
    edgeDef[ePrism] = {{0, 1},
                       {1, 2},
                       {3, 2},
                       {0, 3},
                       {0, 4},
                       {1, 4},
                       {2, 5},
                       {3, 5},
                       {4, 5}};
    edgeDef[ePyramid] = {{0, 1},
                         {1, 2},
                         {3, 2},
                         {0, 3},
                         {0, 4},
                         {1, 4},
                         {2, 4},
                         {3, 4}};
    edgeDef[eTetrahedron] = {{0, 1},
                             {1, 2},
                             {0, 2},
                             {0, 3},
                             {1, 3},
                             {2, 3}};

    faceDef[eTriangle] = {{0, 1, 2}};
    faceDef[eQuadrilateral] = {{0, 1, 2, 3}};
    faceDef[eHexahedron] = {{0, 1, 2,  3},
                            {0, 5, 8,  4},
                            {1, 6, 9,  5},
                            {2, 7, 10, 6},
                            {3, 7, 11, 4},
                            {8, 9, 10, 11}};
    faceDef[ePrism] = {{0, 1, 2, 3},
                       {0, 5, 4},
                       {1, 6, 8, 5},
                       {2, 6, 7},
                       {3, 7, 8, 4}};
    faceDef[ePyramid] = {{0, 1, 2, 3},
                         {0, 5, 4},
                         {1, 6, 5},
                         {2, 6, 7},
                         {3, 7, 4}};
    faceDef[eTetrahedron] = {{0, 1, 2},
                             {0, 4, 3},
                             {1, 5, 4},
                             {2, 5, 3}};

    volDef[eHexahedron] = {{0, 1, 2, 3, 4, 5}};
    volDef[ePrism] = {{0, 1, 2, 3, 4}};
    volDef[ePyramid] = {{0, 1, 2, 3, 4}};
    volDef[eTetrahedron] = {{0, 1, 2, 3}};

    map<ShapeType, int> numVerts;
    numVerts[eSegment] = 2;
    numVerts[eTriangle] = 3;
    numVerts[eQuadrilateral] = 4;
    numVerts[eTetrahedron] = 4;
    numVerts[ePyramid] = 5;
    numVerts[ePrism] = 6;
    numVerts[eHexahedron] = 8;

    int dimension = ShapeTypeDimMap[shapeType];
    ASSERTL0(coords.size() == dimension * numVerts[shapeType],
             ("The number of coordinates supplied should match the shape type, "
              "you supplied " + to_string(coords.size() / dimension) + " and "
                                                                       "shape type " +
              ShapeTypeMap[shapeType] + " requires " +
              to_string(numVerts[shapeType])));


    PointGeomSharedPtr verts[numVerts[shapeType]];
    for (int i = 0; i < numVerts[shapeType]; ++i)
    {
        verts[i] = MemoryManager<PointGeom>::
        AllocateSharedPtr(dimension, i, coords[dimension * i],
                          coords[dimension * i + 1],
                          coords[dimension * i + 2]);
    }

    vector<SegGeomSharedPtr> edges;
    for (int i = 0; i < edgeDef[shapeType].size(); ++i)
    {
        vector<PointGeomSharedPtr> tmp;
        for (int j = 0; j < 2; ++j)
        {
            tmp.push_back(verts[edgeDef[shapeType][i][j]]);
        }
        edges.push_back(MemoryManager<SegGeom>::AllocateSharedPtr(i, dimension,
                                                                  &tmp[0]));
    }

    if (dimension == 1)
    {
        return edges[0];
    }

    vector<Geometry2DSharedPtr> faces;
    for (int i = 0; i < faceDef[shapeType].size(); ++i)
    {
        vector<SegGeomSharedPtr> tmp;
        for (int j = 0; j < faceDef[shapeType][i].size(); ++j)
        {
            tmp.push_back(edges[faceDef[shapeType][i][j]]);
        }
        if (faceDef[shapeType][i].size() == 3)
        {
            faces.push_back(
                    MemoryManager<TriGeom>::AllocateSharedPtr(i, &tmp[0]));
        }
        else
        {
            faces.push_back(
                    MemoryManager<QuadGeom>::AllocateSharedPtr(i, &tmp[0]));
        }

    }

    if (dimension == 2)
    {
        return faces[0];
    }

    vector<Geometry3DSharedPtr> volumes;
    for (int i = 0; i < volDef[shapeType].size(); ++i)
    {
        vector<Geometry2DSharedPtr> tmp;
        for (int j = 0; j < volDef[shapeType][i].size(); ++j)
        {
            tmp.push_back(faces[volDef[shapeType][i][j]]);
        }

        if (shapeType == eTetrahedron)
        {
            vector<TriGeomSharedPtr> tmp2;
            for (int j = 0; j < tmp.size(); ++j)
            {
                tmp2.push_back(dynamic_pointer_cast<TriGeom>(tmp[j]));
            }
            volumes.push_back(
                    MemoryManager<TetGeom>::AllocateSharedPtr(0, &tmp2[0]));
        }
        else if (shapeType == ePyramid)
        {
            volumes.push_back(
                    MemoryManager<PyrGeom>::AllocateSharedPtr(0, &tmp[0]));
        }
        else if (shapeType == ePrism)
        {
            volumes.push_back(
                    MemoryManager<PrismGeom>::AllocateSharedPtr(0, &tmp[0]));
        }
        else if (shapeType == eHexahedron)
        {
            vector<QuadGeomSharedPtr> tmp2;
            for (int j = 0; j < tmp.size(); ++j)
            {
                tmp2.push_back(dynamic_pointer_cast<QuadGeom>(tmp[j]));
            }
            volumes.push_back(
                    MemoryManager<HexGeom>::AllocateSharedPtr(0, &tmp2[0]));
        }
    }

    return volumes[0];
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