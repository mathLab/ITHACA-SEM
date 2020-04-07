////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////
#include "CADSurf.h"
#include "CADCurve.h"

#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

namespace bg = boost::geometry;
typedef bg::model::d2::point_xy<double> point_xy;

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

void CADSurf::OrientateEdges(CADSurfSharedPtr surf,
                             vector<EdgeLoopSharedPtr> &ein)
{
    // this piece of code orientates the surface,
    // it used to be face mesh but its easier to have it here
    int np = 20;
    vector<vector<Array<OneD, NekDouble>>> loopt;
    for (int i = 0; i < ein.size(); i++)
    {
        vector<Array<OneD, NekDouble>> loop;
        for (int j = 0; j < ein[i]->edges.size(); j++)
        {
            Array<OneD, NekDouble> bnds = ein[i]->edges[j]->GetBounds();
            NekDouble dt = (bnds[1] - bnds[0]) / (np - 1);
            if (ein[i]->edgeo[j] == CADOrientation::eForwards)
            {
                for (int k = 0; k < np - 1; k++)
                {
                    NekDouble t = bnds[0] + dt * k;
                    Array<OneD, NekDouble> l  = ein[i]->edges[j]->P(t);
                    Array<OneD, NekDouble> uv = surf->locuv(l);
                    loop.push_back(uv);
                }
            }
            else
            {
                for (int k = np - 1; k > 0; k--)
                {
                    NekDouble t = bnds[0] + dt * k;
                    Array<OneD, NekDouble> l  = ein[i]->edges[j]->P(t);
                    Array<OneD, NekDouble> uv = surf->locuv(l);
                    loop.push_back(uv);
                }
            }
        }
        loopt.push_back(loop);
    }

    vector<bg::model::polygon<point_xy, false, true>> polygons;

    for (int i = 0; i < loopt.size(); i++)
    {
        bg::model::polygon<point_xy, false, true> polygon;
        vector<point_xy> points;
        for (int j = 0; j < loopt[i].size(); j++)
        {
            points.push_back(point_xy(loopt[i][j][0], loopt[i][j][1]));
        }
        // boost requires for closed polygons (last point == first point)
        points.push_back(point_xy(loopt[i][0][0], loopt[i][0][1]));

        bg::assign_points(polygon, points);

        NekDouble area = bg::area(polygon);

        ein[i]->area = area;

        point_xy cen(0.0, 0.0);
        bg::centroid(polygon, cen);

        ein[i]->center    = Array<OneD, NekDouble>(2);
        ein[i]->center[0] = cen.x();
        ein[i]->center[1] = cen.y();

        polygons.push_back(polygon);
    }

    // order by absoulte area
    int ct = 0;
    do
    {
        ct = 0;
        for (int i = 0; i < ein.size() - 1; i++)
        {
            if (fabs(ein[i]->area) < fabs(ein[i + 1]->area))
            {
                // swap
                swap(ein[i], ein[i + 1]);
                swap(loopt[i], loopt[i + 1]);
                swap(polygons[i], polygons[i + 1]);
                ct += 1;
            }
        }

    } while (ct > 0);

    // only need center points for inner loops
    for (int i = 1; i < ein.size(); i++)
    {
        point_xy p(ein[i]->center[0], ein[i]->center[1]);

        if (!bg::within(p, polygons[i]))
        {
            Array<OneD, NekDouble> n1 = loopt[i][0];
            Array<OneD, NekDouble> n2 = loopt[i][1];

            Array<OneD, NekDouble> N(2);
            NekDouble mag = sqrt((n1[0] - n2[0]) * (n1[0] - n2[0]) +
                                 (n1[1] - n2[1]) * (n1[1] - n2[1]));

            N[0] = (n2[1] - n1[1]) / mag;
            N[1] = -1.0 * (n2[0] - n1[0]) / mag;

            Array<OneD, NekDouble> P(2);
            P[0] = n1[0] + N[0];
            P[1] = n1[1] + N[1];

            ein[i]->center = P;

            p = point_xy(P[0], P[1]);

            ASSERTL0(boost::geometry::within(p, polygons[i]),
                     "point is not side loop");
        }
    }
}
}
}
