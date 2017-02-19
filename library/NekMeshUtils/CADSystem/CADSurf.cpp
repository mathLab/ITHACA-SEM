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
//  License for the specific language governing rights and limitations under
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
    vector<vector<Array<OneD, NekDouble> > > loopt;
    for (int i = 0; i < ein.size(); i++)
    {
        vector<Array<OneD, NekDouble> > loop;
        for (int j = 0; j < ein[i]->edges.size(); j++)
        {
            Array<OneD, NekDouble> bnds = ein[i]->edges[j]->GetBounds();
            NekDouble dt = (bnds[1] - bnds[0]) / (20 - 1);
            if (ein[i]->edgeo[j] == eForwards)
            {
                for (int k = 0; k < 20; k++)
                {
                    NekDouble t = bnds[0] + dt * k;
                    Array<OneD, NekDouble> l  = ein[i]->edges[j]->P(t);
                    Array<OneD, NekDouble> uv = surf->locuv(l);
                    loop.push_back(uv);
                }
            }
            else
            {
                for (int k = 19; k >= 0; k--)
                {
                    NekDouble t = bnds[0] + dt * k;
                    Array<OneD, NekDouble> l  = ein[i]->edges[j]->P(t);
                    Array<OneD, NekDouble> uv = surf->locuv(l);
                    loop.push_back(uv);
                }
            }
            loopt.push_back(loop);
        }
    }

    for (int i = 0; i < loopt.size(); i++)
    {
        NekDouble area = 0.0;
        NekDouble mn   = numeric_limits<double>::max();

        for (int j = 0; j < loopt[i].size(); j++)
        {
            mn = min(loopt[i][j][1], mn);
        }

        for (int j = 0; j < loopt[i].size(); j++)
        {
            loopt[i][j][1] -= mn;
        }

        for (int j = 0; j < loopt[i].size() - 1; j++)
        {
            area += (loopt[i][j + 1][0] - loopt[i][j][0]) *
                    (loopt[i][j][1] + loopt[i][j + 1][1]) / 2.0;
        }
        area += (loopt[i][0][0] - loopt[i][loopt[i].size() - 1][0]) *
                (loopt[i][loopt[i].size() - 1][1] + loopt[i][0][1]) / 2.0;

        ein[i]->area = area;
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
                ct += 1;
            }
        }

    } while (ct > 0);

    for (int i = 0; i < ein.size(); i++)
    {
        Array<OneD, NekDouble> n1info, n2info;
        n1info = loopt[i][0];
        n2info = loopt[i][1];

        Array<OneD, NekDouble> N(2);
        NekDouble mag = sqrt((n1info[0] - n2info[0]) * (n1info[0] - n2info[0]) +
                             (n1info[1] - n2info[1]) * (n1info[1] - n2info[1]));
        ASSERTL0(mag > 1e-30, "infinity");
        N[0] = -1.0 * (n2info[1] - n1info[1]) / mag;
        N[1] = (n2info[0] - n1info[0]) / mag;

        Array<OneD, NekDouble> P(2);
        P[0] = (n1info[0] + n2info[0]) / 2.0 + 1e-8 * N[0];
        P[1] = (n1info[1] + n2info[1]) / 2.0 + 1e-8 * N[1];

        // now test to see if p is inside or outside the shape
        // vector to the right
        int intercepts = 0;
        for (int j = 0; j < loopt[i].size() - 1; j++)
        {
            Array<OneD, NekDouble> nt1, nt2;
            nt1 = loopt[i][j];
            nt2 = loopt[i][j + 1];

            if (fabs(nt2[1] - nt1[1]) < 1e-30)
            {
                continue;
            }

            NekDouble lam = (P[1] - nt1[1]) / (nt2[1] - nt1[1]);
            NekDouble S   = nt1[0] - P[0] + (nt2[0] - nt1[0]) * lam;

            if (!(lam < 0) && !(lam > 1) && S > 0)
            {
                intercepts++;
            }
        }
        {
            Array<OneD, NekDouble> nt1, nt2;
            nt1 = loopt[i].back();
            nt2 = loopt[i][0];

            if (fabs(nt2[1] - nt1[1]) < 1e-30)
            {
                continue;
            }

            NekDouble lam = (P[1] - nt1[1]) / (nt2[1] - nt1[1]);
            NekDouble S   = nt1[0] - P[0] + (nt2[0] - nt1[0]) * lam;

            if (!(lam < 0) && !(lam > 1) && S > 0)
            {
                intercepts++;
            }
        }
        if (intercepts % 2 == 0)
        {
            P[0]       = (n1info[0] + n2info[0]) / 2.0 - 1e-6 * N[0];
            P[1]       = (n1info[1] + n2info[1]) / 2.0 - 1e-6 * N[1];
            intercepts = 0;
            for (int j = 0; j < loopt[i].size() - 1; j++)
            {
                Array<OneD, NekDouble> nt1, nt2;
                nt1 = loopt[i][j];
                nt2 = loopt[i][j + 1];

                if (fabs(nt2[1] - nt1[1]) < 1e-30)
                {
                    continue;
                }

                NekDouble lam = (P[1] - nt1[1]) / (nt2[1] - nt1[1]);
                NekDouble S   = nt1[0] - P[0] + (nt2[0] - nt1[0]) * lam;

                if (!(lam < 0) && !(lam > 1) && S > 0)
                {
                    intercepts++;
                }
            }
            {
                Array<OneD, NekDouble> nt1, nt2;
                nt1 = loopt[i].back();
                nt2 = loopt[i][0];

                if (fabs(nt2[1] - nt1[1]) < 1e-30)
                {
                    continue;
                }

                NekDouble lam = (P[1] - nt1[1]) / (nt2[1] - nt1[1]);
                NekDouble S   = nt1[0] - P[0] + (nt2[0] - nt1[0]) * lam;

                if (!(lam < 0) && !(lam > 1) && S > 0)
                {
                    intercepts++;
                }
            }
            if (intercepts % 2 == 0)
            {
                cerr << "still failed to find point inside loop" << endl;
            }
        }

        ein[i]->center = P;
    }

    if (ein[0]->area < 0) // reverse the first uvLoop
    {
        reverse(ein[0]->edgeo.begin(), ein[0]->edgeo.end());
        reverse(ein[0]->edges.begin(), ein[0]->edges.end());
        // need to flip edgeo
        for (int i = 0; i < ein[0]->edgeo.size(); i++)
        {
            if (ein[0]->edgeo[i] == eForwards)
            {
                ein[0]->edgeo[i] = eBackwards;
            }
            else
            {
                ein[0]->edgeo[i] = eForwards;
            }
        }
    }

    for (int i = 1; i < ein.size(); i++)
    {
        if (ein[i]->area > 0) // reverse the loop
        {
            ein[i]->area *= -1.0;
            reverse(ein[i]->edgeo.begin(), ein[i]->edgeo.end());
            reverse(ein[i]->edges.begin(), ein[i]->edges.end());

            // need to flip edgeo
            for (int j = 0; j < ein[i]->edgeo.size(); j++)
            {
                if (ein[i]->edgeo[j] == eForwards)
                {
                    ein[i]->edgeo[j] = eBackwards;
                }
                else
                {
                    ein[i]->edgeo[j] = eForwards;
                }
            }
        }
    }
}

}
}
