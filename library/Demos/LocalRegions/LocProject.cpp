///////////////////////////////////////////////////////////////////////////////
//
// File: NodalDemo.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Demo for testing functionality of LocProject
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <boost/core/ignore_unused.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <LibUtilities/Foundations/PointsType.h>
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

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff);

GeometrySharedPtr CreateGeom(vector<double> coords, ShapeType shapeType);

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
static double pow_loc(const double val, const int i)
{
    return (i < 0) ? 1.0 : pow(val, i);
}

int main(int argc, char *argv[])
{
    string shape, ntype;
    vector<string> basis(3, "NoBasisType"), pointstype(3, "NoPointsType");
    vector<int> order, points;
    vector<double> coords;
    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",
         "Produce this help message and list basis and shape types.")
        ("nodal,n",
         po::value<string>(&ntype),
         "Optional nodal type, autofills shape and basis choices.")
        ("shape,s",
         po::value<string>(&shape),
         "Region shape to project function on.")
        ("basis,b",
         po::value<vector<string>>(&basis)->multitoken(),
         "Basis type, separate by spaces for higher dimensions.")
        ("order,o",
         po::value<vector<int>>(&order)->multitoken()->required(),
         "Order of basis sets, separate by spaces for higher dimensions.")
        ("points,p",
         po::value<vector<int>>(&points)->multitoken()->required(),
         "Number of quadrature points, separate by spaces for "
         "higher dimensions.")
        ("pointstype,P",
         po::value<vector<string>>(&pointstype)->multitoken(),
         "Optional points type, separate by spaces for higher dimensions.")
        ("coords,c",
         po::value<vector<NekDouble>>(&coords)->multitoken()->required(),
         "Coordinates, separate by spaces for higher dimensions with "
         "grouping by vertex i.e. 'x1 y1 z1  x2 y2 z2  x3 y3 z3  ...")
        ("diff,d", "Project derivative.");

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
            cout << endl << "All points types, -P [ --pointstype ], are:"
                 << endl;
            for (int i = 1; i < SIZE_PointsType; ++i)
            {
                cout << kPointsTypeStr[i] << endl;
            };
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl << desc;
        return 0;
    }

    vector<PointsType> ptype;
    if (vm.count("pointstype"))
    {
        for (auto &p : pointstype)
        {
            PointsType tmp = eNoPointsType;
            for (int i = 1; i < SIZE_PointsType; ++i) // starts at nodal points
            {
                if (boost::iequals(kPointsTypeStr[i], p))
                {
                    tmp = static_cast<PointsType>(i);
                    break;
                }
                ASSERTL0(i != SIZE_PointsType - 1,
                         "The points type '" + p + "' does not exist");
            }

            ptype.push_back(tmp);
        }
    }

    // Convert string input argument to nodal type
    PointsType nodaltype = eNoPointsType;
    ShapeType stype = eNoShapeType;
    vector<BasisType> btype(3, eNoBasisType);
    if (vm.count("nodal"))
    {
        for (int i = 1; i < SIZE_PointsType; ++i) // starts at nodal points
        {
            if (boost::iequals(kPointsTypeStr[i], ntype))
            {
                nodaltype = static_cast<PointsType>(i);
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
                stype = static_cast<ShapeType>(i);
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

    // Check arguments supplied equals dimension
    const int dimension = ShapeTypeDimMap[stype];
    ASSERTL0(order.size() == dimension,
             "Number of orders supplied should match shape dimension");
    ASSERTL0(points.size() == dimension,
             "Number of points supplied should match shape dimension");
    ASSERTL0(ptype.size() == dimension || ptype.size() == 0,
             "Number of points types should match shape dimension if "
             "supplied.");

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
                    btype[i] = static_cast<BasisType>(j);
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
    allowableBasis[eSegment] = { allowableBasis[ePoint][0] };
    allowableBasis[eTriangle] = {
        {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
        {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
    };
    allowableBasis[eQuadrilateral] = {
        allowableBasis[eSegment][0], allowableBasis[eSegment][0]
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
        const unsigned int basisListLength = allowableBasis[stype][i].size();
        for (int j = 0; j < basisListLength; ++j)
        {
            if (allowableBasis[stype][i][j] == btype[i])
            {
                break;
            }
            ASSERTL0(j != basisListLength - 1,
                     ("The basis type '" +
                     static_cast<string>(BasisTypeMap[btype[i]]) +
                      "' is invalid for basis argument " + to_string(i + 1) +
                      " for shape '" + ShapeTypeMap[stype] + "'."))
        }
    }

    // Declaration of other variables needed
    StdExpansion *E = nullptr;

    // Assign points type according to basis type selection, if not already
    // assigned.
    if (ptype.size() == 0)
    {
        ptype.resize(dimension);
        for (int i = 0; i < dimension; ++i)
        {
            if (btype[i] == eFourier)
            {
                ptype[i] = eFourierEvenlySpaced;
            }
            else if (btype[i] == eFourierSingleMode ||
                     btype[i] == eFourierHalfModeRe ||
                     btype[i] == eFourierHalfModeIm)
            {
                ptype[i] = eFourierSingleModeSpaced;
            }
            else
            {
                if (i == 1 && (stype == eTriangle || stype == eTetrahedron))
                {
                    ptype[i] = eGaussRadauMAlpha1Beta0;
                }
                else if (i == 2 && (stype == eTetrahedron || stype == ePyramid))
                {
                    ptype[i] = eGaussRadauMAlpha2Beta0;
                }
                else if (i == 2 && stype == ePrism)
                {
                    ptype[i] = eGaussRadauMAlpha1Beta0;
                }
                else
                {
                    ptype[i] = eGaussLobattoLegendre;
                }
            }
        }
    }

    vector<PointsKey> pkey;
    vector<BasisKey> bkey;
    for (int i = 0; i < dimension; ++i)
    {
        pkey.emplace_back(PointsKey(points[i], ptype[i]));
        bkey.emplace_back(BasisKey(btype[i], order[i], pkey[i]));
    }

    GeometrySharedPtr geomFull = CreateGeom(coords, stype);

    switch (stype)
    {
        case eSegment:
        {
            SegGeomSharedPtr geom = static_pointer_cast<SegGeom>(geomFull);
            E = new SegExp(bkey[0], geom);
            break;
        }
        case eTriangle:
        {
            TriGeomSharedPtr geom = static_pointer_cast<TriGeom>(geomFull);
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
            QuadGeomSharedPtr geom = static_pointer_cast<QuadGeom>(geomFull);
            E = new QuadExp(bkey[0], bkey[1], geom);
            break;
        }
        case eTetrahedron:
        {
            TetGeomSharedPtr geom = static_pointer_cast<TetGeom>(geomFull);
            E = new TetExp(bkey[0], bkey[1], bkey[2], geom);
            break;
        }
        case ePyramid:
        {
            PyrGeomSharedPtr geom = static_pointer_cast<PyrGeom>(geomFull);
            E = new PyrExp(bkey[0], bkey[1], bkey[2], geom);
            break;
        }
        case ePrism:
        {
            PrismGeomSharedPtr geom = static_pointer_cast<PrismGeom>(geomFull);
            E = new PrismExp(bkey[0], bkey[1], bkey[2], geom);
            break;
        }
        case eHexahedron:
        {
            HexGeomSharedPtr geom = static_pointer_cast<HexGeom>(geomFull);
            E = new HexExp(bkey[0], bkey[1], bkey[2], geom);
            break;
        }
        default:
            break;
    }

    const auto totPoints = (unsigned) E->GetTotPoints();
    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> y = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> z = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dx = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dy = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dz = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> sol = Array<OneD, NekDouble>(totPoints);

    switch (dimension)
    {
        case 1:
        {
            E->GetCoords(x);
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
    }

    //get solution array
    for (int i = 0; i < totPoints; ++i)
    {
        sol[i] = Shape_sol(x[i], y[i], z[i], order, btype, stype, false);
    }

    Array<OneD, NekDouble> phys(totPoints);
    Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());

    if (vm.count("diff"))
    {
        switch (dimension)
        {
            case 1:
            {
                E->PhysDeriv(sol, sol);
                break;
            }
            case 2:
            {
                E->PhysDeriv(sol, dx, dy);
                Vmath::Vadd(totPoints, dx, 1, dy, 1, sol, 1);
                break;
            }
            case 3:
            {
                E->PhysDeriv(sol, dx, dy, dz);
                Vmath::Vadd(totPoints, dx, 1, dy, 1, sol, 1);
                Vmath::Vadd(totPoints, dz, 1, sol, 1, sol, 1);
                break;
            }
        }
    }

    //Project onto expansion
    E->FwdTrans(sol, coeffs);

    //Backward transform solution to get projected values
    E->BwdTrans(coeffs, phys);

    if (vm.count("diff"))
    {
        for (int i = 0; i < totPoints; ++i)
        {
            sol[i] = Shape_sol(x[i], y[i], z[i], order, btype, stype, true);
        }
    }

    //Calculate L_inf & L_2 error
    cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
    cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;

    return 0;
}

GeometrySharedPtr CreateGeom(vector<NekDouble> coords, ShapeType shapeType)
{
    map<ShapeType, vector<vector<unsigned int>>> edgeDef, faceDef;
    edgeDef[eSegment]       = {{0, 1}};
    edgeDef[eTriangle]      = {{0, 1}, {1, 2}, {2, 0}};
    edgeDef[eQuadrilateral] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    edgeDef[eHexahedron]    = {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {0, 4}, {1, 5},
                               {2, 6}, {3, 7}, {4, 5}, {5, 6}, {7, 6}, {4, 7}};
    edgeDef[ePrism]         = {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {0, 4}, {1, 4},
                               {2, 5}, {3, 5}, {4, 5}};
    edgeDef[ePyramid]       = {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {0, 4}, {1, 4},
                               {2, 4}, {3, 4}};
    edgeDef[eTetrahedron]   = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};

    faceDef[eTriangle]      = {{0, 1, 2}};
    faceDef[eQuadrilateral] = {{0, 1, 2, 3}};
    faceDef[eHexahedron]    = {{0, 1, 2,  3}, {0, 5, 8,  4}, {1, 6, 9,  5},
                               {2, 7, 10, 6}, {3, 7, 11, 4}, {8, 9, 10, 11}};
    faceDef[ePrism]         = {{0, 1, 2, 3}, {0, 5, 4}, {1, 6, 8, 5}, {2, 6, 7},
                               {3, 7, 8, 4}};
    faceDef[ePyramid]       = {{0, 1, 2, 3}, {0, 5, 4}, {1, 6, 5}, {2, 6, 7},
                               {3, 7, 4}};
    faceDef[eTetrahedron]   = {{0, 1, 2}, {0, 4, 3}, {1, 5, 4}, {2, 5, 3}};

    map<ShapeType, vector<unsigned int>> volDef;
    volDef[eHexahedron] = {0, 1, 2, 3, 4, 5};
    volDef[ePrism] = {0, 1, 2, 3, 4};
    volDef[ePyramid] = {0, 1, 2, 3, 4};
    volDef[eTetrahedron] = {0, 1, 2, 3};

    map<ShapeType, unsigned int> numVerts;
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
              "shape type " + ShapeTypeMap[shapeType] + " requires " +
              to_string(numVerts[shapeType])));

    //Set up vector of vertices
    vector<PointGeomSharedPtr> verts(numVerts[shapeType]);
    for (int i = 0; i < numVerts[shapeType]; ++i)
    {
        if (dimension == 1)
        {
            verts[i] = MemoryManager<PointGeom>::AllocateSharedPtr(
                    dimension, i, coords[dimension * i], 0, 0);
        }
        if (dimension == 2)
        {
            verts[i] = MemoryManager<PointGeom>::AllocateSharedPtr(
                    dimension, i, coords[dimension * i],
                    coords[dimension * i + 1], 0);
        }
        if (dimension == 3)
        {
            verts[i] = MemoryManager<PointGeom>::AllocateSharedPtr(
                    dimension, i, coords[dimension * i],
                    coords[dimension * i + 1],
                    coords[dimension * i + 2]);
        }
    }

    //Set up vector of edges
    const unsigned int numEdges = edgeDef[shapeType].size();
    vector<SegGeomSharedPtr> edges(numEdges);
    for (int i = 0; i < numEdges; ++i)
    {
        vector<PointGeomSharedPtr> tmp(2);
        for (int j = 0; j < 2; ++j)
        {
            tmp[j] = verts[edgeDef[shapeType][i][j]];
        }
        edges[i] = MemoryManager<SegGeom>::AllocateSharedPtr(
                i, dimension, &tmp[0]);
    }

    if (dimension == 1)
    {
        return edges[0];
    }

    //Set up vector of faces
    const unsigned int numFaces = faceDef[shapeType].size();
    vector<Geometry2DSharedPtr> faces(numFaces);
    for (int i = 0; i < numFaces; ++i)
    {
        const unsigned int numEdgesFace = faceDef[shapeType][i].size();
        vector<SegGeomSharedPtr> tmp(numEdgesFace);
        for (int j = 0; j < numEdgesFace; ++j)
        {
            tmp[j] = edges[faceDef[shapeType][i][j]];
        }
        if (numEdgesFace == 3)
        {
            faces[i] = MemoryManager<TriGeom>::AllocateSharedPtr(i, &tmp[0]);
        }
        else
        {
            faces[i] = MemoryManager<QuadGeom>::AllocateSharedPtr(i, &tmp[0]);
        }

    }

    if (dimension == 2)
    {
        return faces[0];
    }

    //Set up volume
    Geometry3DSharedPtr volume;
    vector<Geometry2DSharedPtr> tmp;
    for (auto volEdge : volDef[shapeType])
    {
        tmp.emplace_back(faces[volEdge]);
    }
    switch (shapeType)
    {
        case eTetrahedron:
        {
            const unsigned int numVol = tmp.size();
            vector<TriGeomSharedPtr> tmp2(numVol);
            for (int i = 0; i < numVol; ++i)
            {
                tmp2[i] = static_pointer_cast<TriGeom>(tmp[i]);
            }
            volume = MemoryManager<TetGeom>::AllocateSharedPtr(0, &tmp2[0]);
            break;
        }
        case ePyramid:
        {
            volume = MemoryManager<PyrGeom>::AllocateSharedPtr(0, &tmp[0]);
            break;
        }
        case ePrism:
        {
            volume = MemoryManager<PrismGeom>::AllocateSharedPtr(0, &tmp[0]);
            break;
        }
        case eHexahedron:
        {
            const unsigned int numVol = tmp.size();
            vector<QuadGeomSharedPtr> tmp2(numVol);
            for (int i = 0; i < numVol; ++i)
            {
                tmp2[i] = static_pointer_cast<QuadGeom>(tmp[i]);
            }
            volume = MemoryManager<HexGeom>::AllocateSharedPtr(0, &tmp2[0]);
            break;
        }
        default:
            break;
    }

    return volume;
}

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff)
{
    map<ShapeType, function<int(int, const vector<int> &)>> shapeConstraint2;
    shapeConstraint2[eSegment] =
            [](int,   const vector<int> &     ) { return 1; };
    shapeConstraint2[eTriangle] =
            [](int k, const vector<int> &order) { return order[1] - k; };
    shapeConstraint2[eQuadrilateral] =
            [](int,   const vector<int> &order) { return order[1]; };
    shapeConstraint2[eTetrahedron] =
            [](int k, const vector<int> &order) { return order[1] - k; };
    shapeConstraint2[ePyramid] =
            [](int k, const vector<int> &order) { return order[1] - k; };
    shapeConstraint2[ePrism] =
            [](int,   const vector<int> &order) { return order[1]; };
    shapeConstraint2[eHexahedron] =
            [](int,   const vector<int> &order) { return order[1]; };

    map<ShapeType, function<int(int, int, const vector<int> &order)>>
                            shapeConstraint3;
    shapeConstraint3[eSegment] =
            [](int,   int,   const vector<int> &     ) { return 1; };
    shapeConstraint3[eTriangle] =
            [](int,   int,   const vector<int> &     ) { return 1; };
    shapeConstraint3[eQuadrilateral] =
            [](int,   int,   const vector<int> &     ) { return 1; };
    shapeConstraint3[eTetrahedron] =
            [](int k, int l, const vector<int> &order) { return order[2] - k - l; };
    shapeConstraint3[ePyramid] =
            [](int k, int l, const vector<int> &order) { return order[2] - k - l; };
    shapeConstraint3[ePrism] =
            [](int k, int,   const vector<int> &order) { return order[2] - k; };
    shapeConstraint3[eHexahedron] =
            [](int,   int,   const vector<int> &order) { return order[2]; };

    NekDouble sol = 0.0;
    if (!diff)
    {
        if (btype[0] == eFourier && stype == eSegment)
        {
            for (int k = 0; k < order[0] / 2 - 1; ++k)
            {
                sol += sin(k * M_PI * x) + cos(k * M_PI * x);
            }
        }
        else if (btype[0] == eFourierSingleMode && stype == eSegment)
        {
            sol += 0.25 * sin(M_PI * x) + 0.25 * cos(M_PI * x);
        }
        else if (btype[0] == eFourier && stype == eQuadrilateral)
        {
            if (btype[1] == eFourier)
            {
                for (int k = 0; k < order[0] / 2; ++k)
                {
                    for (int l = 0; l < order[1] / 2; ++l)
                    {
                        sol += sin(k * M_PI * x) * sin(l * M_PI * y) +
                               sin(k * M_PI * x) * cos(l * M_PI * y) +
                               cos(k * M_PI * x) * sin(l * M_PI * y) +
                               cos(k * M_PI * x) * cos(l * M_PI * y);
                    }
                }
            }
            else if (btype[1] == eFourierSingleMode)
            {
                for (int k = 0; k < order[0] / 2; ++k)
                {
                    sol += sin(k * M_PI * x) * sin(M_PI * y) +
                           sin(k * M_PI * x) * cos(M_PI * y) +
                           cos(k * M_PI * x) * sin(M_PI * y) +
                           cos(k * M_PI * x) * cos(M_PI * y);
                }
            }
            else
            {
                for (int k = 0; k < order[0] / 2; ++k)
                {
                    for (int l = 0; l < order[1]; ++l)
                    {
                        sol += sin(k * M_PI * x) * pow_loc(y, l) +
                               cos(k * M_PI * x) * pow_loc(y, l);
                    }
                }
            }
        }
        else if (btype[0] == eFourierSingleMode && stype == eQuadrilateral)
        {
            if (btype[1] == eFourier)
            {
                for (int l = 0; l < order[1] / 2; ++l)
                {
                    sol += sin(M_PI * x) * sin(l * M_PI * y) +
                           sin(M_PI * x) * cos(l * M_PI * y) +
                           cos(M_PI * x) * sin(l * M_PI * y) +
                           cos(M_PI * x) * cos(l * M_PI * y);
                }

            }
            else if (btype[1] == eFourierSingleMode)
            {
                sol += sin(M_PI * x) * sin(M_PI * y) +
                       sin(M_PI * x) * cos(M_PI * y) +
                       cos(M_PI * x) * sin(M_PI * y) +
                       cos(M_PI * x) * cos(M_PI * y);
            }
            else
            {
                for (int l = 0; l < order[1]; ++l)
                {
                    sol += sin(M_PI * x) * pow_loc(y, l) +
                           cos(M_PI * x) * pow_loc(y, l);
                }
            }
        }
        else if (btype[1] == eFourier && stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0]; ++k)
            {
                for (int l = 0; l < order[1] / 2; ++l)
                {
                    sol += sin(l * M_PI * y) * pow_loc(x, k) +
                           cos(l * M_PI * y) * pow_loc(x, k);
                }
            }
        }
        else if (btype[1] == eFourierSingleMode && stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0]; ++k)
            {
                sol += sin(M_PI * y) * pow_loc(x, k) +
                       cos(M_PI * y) * pow_loc(x, k);
            }
        }
        else
        {
            for (int k = 0;
                 k < order[0]; ++k) //ShapeConstraint 1 is always < order1
            {
                for (int l = 0; l < shapeConstraint2[stype](k, order); ++l)
                {
                    for (int m = 0;
                         m < shapeConstraint3[stype](k, l, order); ++m)
                    {
                        sol += pow_loc(x, k) * pow_loc(y, l) * pow_loc(z, m);
                    }
                }
            }
        }
    }
    else if (diff)
    {
        if (btype[0] == eFourier && stype == eSegment)
        {
            for (int k = 0; k < order[0] / 2 - 1; ++k)
            {
                sol += k * M_PI * (cos(k * M_PI * z) - sin(k * M_PI * z));
            }
        }
        else if (btype[0] != eFourier && btype[1] == eFourier &&
                 stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0]; ++k)
            {
                for (int l = 0; l < order[1] / 2; ++l)
                {
                    sol += k * pow_loc(x, k - 1) * sin(M_PI * l * y)
                           + M_PI * l * pow_loc(x, k) * cos(M_PI * l * y) +
                           +k * pow_loc(x, k - 1) * cos(M_PI * l * y)
                           - M_PI * l * pow_loc(x, k) * sin(M_PI * l * y);
                }
            }
        }
        else if (btype[0] == eFourier && btype[1] != eFourier &&
                 stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0] / 2; ++k)
            {
                for (int l = 0; l < order[1]; ++l)
                {
                    sol += M_PI * k * cos(M_PI * k * x) * pow_loc(y, l)
                           + l * sin(M_PI * k * x) * pow_loc(y, l - 1) +
                           -M_PI * k * sin(M_PI * k * x) * pow_loc(y, l)
                           + l * sin(M_PI * k * x) * pow_loc(y, l - 1);
                }
            }
        }
        else if (btype[0] == eFourier && btype[1] == eFourier &&
                 stype == eQuadrilateral)
        {
            for (int k = 0; k < order[0] / 2; ++k)
            {
                for (int l = 0; l < order[1] / 2; ++l)
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
        else
        {
            NekDouble a;
            for (int k = 0;
                 k < order[0]; ++k) //ShapeConstraint 1 is always < order1
            {
                for (int l = 0; l < shapeConstraint2[stype](k, order); ++l)
                {
                    for (int m = 0;
                         m < shapeConstraint3[stype](k, l, order); ++m)
                    {
                        a = k * pow_loc(x, k - 1) * pow_loc(y, l) *
                            pow_loc(z, m);
                        sol += a;
                        a = l * pow_loc(x, k) * pow_loc(y, l - 1) *
                            pow_loc(z, m);
                        sol += a;
                        a = m * pow_loc(x, k) * pow_loc(y, l) *
                            pow_loc(z, m - 1);
                        sol += a;
                    }
                }
            }
        }
    }

    return sol;
}
