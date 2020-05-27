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
// Description: Demo for testing functionality of NodalUtil classes
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.h>

#include <LibUtilities/BasicUtils/ShapeType.hpp>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",         "Produce this help message.")
        ("integral,i",     "Evaluate the integral of a known function and "
                           "return error.")
        ("order,o",        po::value<int>(),
                           "Order to evaluate at.")
        ("type,t",         po::value<int>(),
                           "Points type of the nodal element.")
        ("deriv,d",        "Evaluate the derivative of a known function and "
                           "return error.")
        ("interp,p",       po::value<string>(),
                           "Interpolate function at a known point of the form"
                           "x,y,z and return error.");

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        cerr << e.what() << endl;
        cerr << desc;
        return 1;
    }

    if (!vm.count("order"))
    {
        cerr << "Must supply an element order." << endl;
        return 1;
    }
    int order = vm["order"].as<int>();

    map<PointsType, ShapeType> nodalTypes;
    map<PointsType, ShapeType>::iterator nIt;

    nodalTypes[eNodalTriElec]           = eTriangle;
    nodalTypes[eNodalTriFekete]         = eTriangle;
    nodalTypes[eNodalTriEvenlySpaced]   = eTriangle;
    nodalTypes[eNodalTetEvenlySpaced]   = eTetrahedron;
    nodalTypes[eNodalTetElec]           = eTetrahedron;
    nodalTypes[eNodalPrismEvenlySpaced] = ePrism;

    if (!vm.count("type"))
    {
        cerr << "Must supply a points type. Valid points types are:" << endl;

        for (nIt = nodalTypes.begin(); nIt != nodalTypes.end(); ++nIt)
        {
            cerr << "  - " << (int)nIt->first << " ("
                 << kPointsTypeStr[nIt->first] << ")" << endl;
        }

        return 1;
    }

    PointsType pointsType = (PointsType)vm["type"].as<int>();

    if (nodalTypes.find(pointsType) == nodalTypes.end())
    {
        cerr << "Invalid points type. Valid points types are:" << endl;

        for (nIt = nodalTypes.begin(); nIt != nodalTypes.end(); ++nIt)
        {
            cerr << "  - " << (int)nIt->first << " ("
                 << kPointsTypeStr[nIt->first] << ")" << endl;
        }

        return 1;
    }

    // Generate nodal points.
    PointsKey pointsKey(order + 1, pointsType);
    PointsSharedPtr points = PointsManager()[pointsKey];

    Array<OneD, NekDouble> r, s, t;
    points->GetPoints(r, s, t);

    // Generate nodal utility.
    ShapeType shape = nodalTypes[pointsType];
    NodalUtil *util = NULL;

    if (shape == eTriangle)
    {
        util = new NodalUtilTriangle(order, r, s);
    }
    else if (shape == eTetrahedron)
    {
        util = new NodalUtilTetrahedron(order, r, s, t);
    }
    else if (shape == ePrism)
    {
        util = new NodalUtilPrism(order, r, s, t);
    }
    else if(shape == eQuadrilateral)
    {
        util = new NodalUtilQuad(order, r, s);
    }

    ASSERTL1(util, "Unknown shape type!");
    const int nPoints = r.size();
    const int dim = (shape == eTriangle || shape == eQuadrilateral) ? 2 : 3;

    if (vm.count("integral"))
    {
        NekVector<NekDouble> weights = util->GetWeights();
        NekDouble integral = 0.0;

        for (int i = 0; i < nPoints; ++i)
        {
            NekDouble integrand = dim == 2 ?
                exp(r[i]) * sin(s[i]) : exp(r[i] + s[i] + t[i]);
            integral += weights[i] * integrand;
        }

        NekDouble exact = 0.0;
        switch(shape)
        {
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
            default:
                exact = 0.0;
                break;
        }

        cout << "L inf error : " << fabs(exact - integral) << endl;
    }
    else if (vm.count("deriv"))
    {
        Array<OneD, NekVector<NekDouble> > exact(dim);
        NekVector<NekDouble> input(nPoints), output(nPoints);

        for (int i = 0; i < dim; ++i)
        {
            exact[i] = NekVector<NekDouble>(nPoints);
        }

        if (dim == 2)
        {
            // Exact solution: u(x,y) = sin(x) * cos(y)
            for (int i = 0; i < nPoints; ++i)
            {
                input[i]    = sin(r[i]) * cos(s[i]);
                exact[0][i] = cos(r[i]) * cos(s[i]);
                exact[1][i] = -sin(r[i]) * sin(s[i]);
            }
        }
        else
        {
            // Exact solution: u(x,y) = sin(x) * cos(y) * sin(z)
            for (int i = 0; i < nPoints; ++i)
            {
                input[i]    = sin(r[i]) * cos(s[i]) * sin(t[i]);
                exact[0][i] = cos(r[i]) * cos(s[i]) * sin(t[i]);
                exact[1][i] = -sin(r[i]) * sin(s[i]) * sin(t[i]);
                exact[2][i] = sin(r[i]) * cos(s[i]) * cos(t[i]);
            }
        }

        string vars = "uvw";

        for (int i = 0; i < dim; ++i)
        {
            std::shared_ptr<NekMatrix<NekDouble> > deriv =
                util->GetDerivMatrix(i);
            output = *deriv * input;
            NekVector<NekDouble> tmp = output - exact[i];

            cout << "L 2 error (variable " << vars[i] << ") : " << tmp.L2Norm()
                 << endl;
        }
    }
    else if (vm.count("interp"))
    {
        string pointStr = vm["interp"].as<string>();
        vector<string> point;
        boost::split(point, pointStr, boost::is_any_of(","));

        if (point.size() != dim)
        {
            cerr << "Point " << pointStr << " does not have correct dimensions"
                 << endl;
            return 1;
        }

        Array<OneD, Array<OneD, NekDouble> > tmp(dim);
        for (int i = 0; i < dim; ++i)
        {
            tmp[i] = Array<OneD, NekDouble>(1);
            try
            {
                tmp[i][0] = boost::lexical_cast<NekDouble>(point[i]);
            }
            catch(boost::bad_lexical_cast &)
            {
                cerr << "Could not convert " << point[i] << " to a coordinate"
                     << endl;
                return 1;
            }
        }

        NekVector<NekDouble> input(nPoints);
        for (int i = 0; i < nPoints; ++i)
        {
            input[i] = dim == 2 ? exp(r[i]) * exp(s[i]) :
                exp(r[i]) * exp(s[i]) * exp(t[i]);
        }

        std::shared_ptr<NekMatrix<NekDouble> > interp =
            util->GetInterpolationMatrix(tmp);

        NekVector<NekDouble> output = *interp * input;
        NekDouble exact = dim == 2 ? exp(tmp[0][0]) * exp(tmp[1][0]) :
            exp(tmp[0][0]) * exp(tmp[1][0]) * exp(tmp[2][0]);

        cout << "L inf error : " << fabs(exact - output[0]) << endl;
    }
    else
    {
        cerr << "You must specify one of --integral, --deriv or --interp"
             << endl;
        cerr << desc;
        return 1;
    }

    return 0;
}
