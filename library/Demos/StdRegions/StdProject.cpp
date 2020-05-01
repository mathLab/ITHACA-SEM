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
// Description: Demo for testing functionality of StdProject
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "StdDemoSupport.hpp"

namespace po = boost::program_options;

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff);

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
static double pow_loc(const double val, const int i)
{
    return (i < 0) ? 1.0 : pow(val, i);
}

int main(int argc, char *argv[])
{
    DemoSupport demo;

    demo.GetOptions().add_options()("diff,d", "Project derivative.");
    demo.ParseArguments(argc, argv);

    po::variables_map vm = demo.GetVariableMap();

    StdExpansion *E = demo.CreateStdExpansion();
    int dimension = E->GetShapeDimension();

    if (E == nullptr)
    {
        return 1;
    }

    std::vector<int> order;
    std::vector<BasisType> btype;
    LibUtilities::ShapeType stype = E->DetShapeType();

    for (int i = 0; i < E->GetShapeDimension(); ++i)
    {
        btype.push_back(E->GetBasisType(i));
        order.push_back(E->GetBasisNumModes(i));
    }

    if (stype == ePoint && vm.count("diff"))
    {
        NEKERROR(ErrorUtil::efatal,
                 "It is not possible to run the diff version for shape: point");
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
        default:
            break;
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
    if (stype != ePoint)
    {
        cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }

    if (!vm.count("diff") && stype != ePoint)
    {
        //Evaluate solution at x = y = 0 and print error
        Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
        t[0] = -0.5;
        t[1] = -0.25;
        t[2] = -0.3;
        sol[0] = Shape_sol(t[0], t[1], t[2], order, btype, stype, false);

        NekDouble nsol = E->PhysEvaluate(t, phys);

        cout << "Error at x = (";
        for (int i = 0; i < dimension; ++i)
        {
            cout << t[i] << ", ";
        }
        cout << "\b\b): " << nsol - sol[0] << endl;
    }

    // Calculate integral of known function to test different quadrature
    // distributions on each element.
    for (int i = 0; i < totPoints; ++i)
    {
        sol[i] = dimension == 1 ? exp(x[i]) : dimension == 2 ?
            exp(x[i]) * sin(y[i]) : exp(x[i] + y[i] + z[i]);
    }

    NekDouble exact = 0.0;
    switch(stype)
    {
        case eSegment:
            exact = M_E - 1.0 / M_E;
            break;
        case eTriangle:
            exact = -0.5 * (sin(1.0) + cos(1.0) + M_E * M_E *
                            (sin(1.0) - cos(1.0))) / M_E;
            break;
        case eQuadrilateral:
            exact = 2.0 * (M_E - 1.0 / M_E) * sin(1.0);
            break;
        case eTetrahedron:
            exact = 1.0 / M_E - 1.0 / M_E / M_E / M_E;
            break;
        case ePrism:
            exact = M_E - 1.0 / M_E / M_E / M_E;
            break;
        case ePyramid:
            exact = - 1.0 / M_E / M_E / M_E - 4.0 / M_E + M_E;
            break;
        case eHexahedron:
            exact = pow((M_E * M_E - 1.0) / M_E, 3.0);
            break;
        default:
            ASSERTL0(false, "Exact solution not known.");
            break;
    }
    std::cout << "Integral error: " << fabs(exact - E->Integral(sol))
              << std::endl;

    return 0;
}

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff)
{
    map<ShapeType, function<int(int, const vector<int> &)>> shapeConstraint2;
    shapeConstraint2[ePoint] =
            [](int,   const vector<int> &     ) { return 1; };
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
    shapeConstraint3[ePoint] =
            [](int,   int,   const vector<int> &     ) { return 1; };
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
