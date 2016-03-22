////////////////////////////////////////////////////////////////////////////////
//
//  File: OptimseFunctions.cpp
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
//  Description: functions and gradients for optimisation
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/SurfaceMeshing/OptimiseFunctions.h>
#include <NekMeshUtils/CADSystem/CADCurve.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

NekDouble Dot(Array<OneD, NekDouble> a, Array<OneD, NekDouble> b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Array<OneD, NekDouble> Take(Array<OneD, NekDouble> a, Array<OneD, NekDouble> b)
{
    Array<OneD, NekDouble> ret(3);
    ret[0] = a[0] - b[0];
    ret[1] = a[1] - b[1];
    ret[2] = a[2] - b[2];
    return ret;
}
Array<OneD, NekDouble> Times(NekDouble t, Array<OneD, NekDouble> a)
{
    Array<OneD, NekDouble> ret(3);
    ret[0] = a[0] * t;
    ret[1] = a[1] * t;
    ret[2] = a[2] * t;
    return ret;
}
Array<OneD, NekDouble> Add(Array<OneD, NekDouble> a, Array<OneD, NekDouble> b)
{
    Array<OneD, NekDouble> ret(3);
    ret[0] = a[0] + b[0];
    ret[1] = a[1] + b[1];
    ret[2] = a[2] + b[2];
    return ret;
}

Array<OneD, NekDouble> OptiEdge::Getxi()
{
    Array<OneD, NekDouble> xi;

    switch (o->GetType())
    {
        case curve:
            xi = Array<OneD, NekDouble>(all.num_elements() - 2);
            for (int i = 1; i < all.num_elements() - 1; i++)
            {
                xi[i - 1] = all[i];
            }
            break;

        case surf:
            xi = Array<OneD, NekDouble>(all.num_elements() - 4);
            for (int i = 2; i < all.num_elements() - 2; i++)
            {
                xi[i - 2] = all[i];
            }
            break;

        case vert:
            ASSERTL0(false, "Should not be able to pass vert");
    }
    return xi;
}

Array<OneD, NekDouble> OptiEdge::Getli()
{
    Array<OneD, NekDouble> li;
    Array<OneD, NekDouble> bnds;
    switch (o->GetType())
    {
        case curve:
            li   = Array<OneD, NekDouble>(all.num_elements() - 2);
            bnds = boost::dynamic_pointer_cast<CADCurve>(o)->Bounds();
            for (int i = 1; i < all.num_elements() - 1; i++)
            {
                li[i - 1] = bnds[0];
            }
            break;

        case surf:
            li   = Array<OneD, NekDouble>(all.num_elements() - 4);
            bnds = boost::dynamic_pointer_cast<CADSurf>(o)->GetBounds();
            for (int i = 2; i < all.num_elements() - 2; i++)
            {
                if (i % 2 == 0)
                {
                    li[i - 2] = bnds[0];
                }
                else
                {
                    li[i - 2] = bnds[2];
                }
            }
            break;

        case vert:
            ASSERTL0(false, "Should not be able to pass vert");
    }
    return li;
}

Array<OneD, NekDouble> OptiEdge::Getui()
{
    Array<OneD, NekDouble> ui;
    Array<OneD, NekDouble> bnds;
    switch (o->GetType())
    {
        case curve:
            ui   = Array<OneD, NekDouble>(all.num_elements() - 2);
            bnds = boost::dynamic_pointer_cast<CADCurve>(o)->Bounds();
            for (int i = 1; i < all.num_elements() - 1; i++)
            {
                ui[i - 1] = bnds[1];
            }
            break;

        case surf:
            ui   = Array<OneD, NekDouble>(all.num_elements() - 4);
            bnds = boost::dynamic_pointer_cast<CADSurf>(o)->GetBounds();
            for (int i = 2; i < all.num_elements() - 2; i++)
            {
                if (i % 2 == 0)
                {
                    ui[i - 2] = bnds[1];
                }
                else
                {
                    ui[i - 2] = bnds[3];
                }
            }
            break;

        case vert:
            ASSERTL0(false, "Should not be able to pass vert");
    }
    return ui;
}

NekDouble OptiEdge::F(Array<OneD, NekDouble> xitst)
{
    Array<OneD, NekDouble> val(all.num_elements());

    if (o->GetType() == curve)
    {
        val[0] = all[0];
        for (int i = 0; i < xitst.num_elements(); i++)
        {
            val[i + 1] = xitst[i];
        }
        val[all.num_elements() - 1] = all[all.num_elements() - 1];
    }
    else if (o->GetType() == surf)
    {
        val[0] = all[0];
        val[1] = all[1];
        for (int i = 0; i < xitst.num_elements(); i++)
        {
            val[i + 2] = xitst[i];
        }
        val[all.num_elements() - 2] = all[all.num_elements() - 2];
        val[all.num_elements() - 1] = all[all.num_elements() - 1];
    }

    NekDouble ret = 0.0;
    if (o->GetType() == curve)
    {
        CADCurveSharedPtr c = boost::dynamic_pointer_cast<CADCurve>(o);

        for (int i = 0; i < all.num_elements() - 1; i++)
        {
            Array<OneD, NekDouble> dis = Take(c->P(val[i + 1]), c->P(val[i]));
            NekDouble norm =
                dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2];
            ret += norm / (z[i + 1] - z[i]);
        }
    }
    else if (o->GetType() == surf)
    {
        CADSurfSharedPtr s = boost::dynamic_pointer_cast<CADSurf>(o);
        // need to organise the val array
        Array<OneD, Array<OneD, NekDouble> > uv(val.num_elements() / 2);
        for (int i = 0; i < val.num_elements() / 2; i++)
        {
            uv[i]    = Array<OneD, NekDouble>(2);
            uv[i][0] = val[i * 2 + 0];
            uv[i][1] = val[i * 2 + 1];
        }
        for (int i = 0; i < uv.num_elements() - 1; i++)
        {
            Array<OneD, NekDouble> dis = Take(s->P(uv[i + 1]), s->P(uv[i]));
            NekDouble norm =
                dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2];
            ret += norm / (z[i + 1] - z[i]);
        }
    }
    return ret;
}

DNekMat OptiEdge::dF(Array<OneD, NekDouble> xitst)
{
    Array<OneD, NekDouble> val(all.num_elements());

    if (o->GetType() == curve)
    {
        val[0] = all[0];
        for (int i = 0; i < xitst.num_elements(); i++)
        {
            val[i + 1] = xitst[i];
        }
        val[all.num_elements() - 1] = all[all.num_elements() - 1];
    }
    else if (o->GetType() == surf)
    {
        val[0] = all[0];
        val[1] = all[1];
        for (int i = 0; i < xitst.num_elements(); i++)
        {
            val[i + 2] = xitst[i];
        }
        val[all.num_elements() - 2] = all[all.num_elements() - 2];
        val[all.num_elements() - 1] = all[all.num_elements() - 1];
    }

    DNekMat ret;

    if (o->GetType() == curve)
    {
        CADCurveSharedPtr c = boost::dynamic_pointer_cast<CADCurve>(o);
        vector<Array<OneD, NekDouble> > r;
        vector<Array<OneD, NekDouble> > dr;

        for (int i = 0; i < all.num_elements(); i++)
        {
            Array<OneD, NekDouble> ri(3), dri(3);
            Array<OneD, NekDouble> d2 = c->D2(val[i]);
            for (int j = 0; j < 3; j++)
            {
                ri[j]  = d2[j];
                dri[j] = d2[j + 3];
            }
            r.push_back(ri);
            dr.push_back(dri);
        }

        DNekMat J(all.num_elements() - 2, 1, 0.0);
        for (int i = 0; i < all.num_elements() - 2; i++)
        {
            J(i, 0) =
                2.0 / (z[i + 1] - z[i]) * Dot(dr[i + 1], Take(r[i + 1], r[i])) -
                2.0 / (z[i + 2] - z[i + 1]) *
                    Dot(dr[i + 1], Take(r[i + 2], r[i + 1]));
        }

        ret = J;
    }
    else if (o->GetType() == surf)
    {
        CADSurfSharedPtr s = boost::dynamic_pointer_cast<CADSurf>(o);
        // need to organise the all array
        Array<OneD, Array<OneD, NekDouble> > uv(val.num_elements() / 2);
        for (int i = 0; i < val.num_elements() / 2; i++)
        {
            uv[i]    = Array<OneD, NekDouble>(2);
            uv[i][0] = val[i * 2 + 0];
            uv[i][1] = val[i * 2 + 1];
        }

        vector<Array<OneD, NekDouble> > r;
        vector<Array<OneD, NekDouble> > dru, drv;
        for (int i = 0; i < uv.num_elements(); i++)
        {
            Array<OneD, NekDouble> ri(3), drui(3), drvi(3);
            Array<OneD, NekDouble> d2 = s->D1(uv[i]);
            for (int j = 0; j < 3; j++)
            {
                ri[j]   = d2[j];
                drui[j] = d2[j + 3];
                drvi[j] = d2[j + 6];
            }
            r.push_back(ri);
            dru.push_back(drui);
            drv.push_back(drvi);
        }

        DNekMat J(2 * (uv.num_elements() - 2), 1, 0.0);
        for (int i = 0; i < uv.num_elements() - 2; i++)
        {
            J(2 * i + 0, 0) = 2.0 / (z[i + 1] - z[i]) *
                                  Dot(dru[i + 1], Take(r[i + 1], r[i])) +
                              2.0 / (z[i + 2] - z[i + 1]) *
                                  Dot(dru[i + 1], Take(r[i + 1], r[i + 2]));

            J(2 * i + 1, 0) = 2.0 / (z[i + 1] - z[i]) *
                                  Dot(drv[i + 1], Take(r[i + 1], r[i])) +
                              2.0 / (z[i + 2] - z[i + 1]) *
                                  Dot(drv[i + 1], Take(r[i + 1], r[i + 2]));
        }

        ret = J;
    }

    return ret;
}

void OptiEdge::Update(Array<OneD, NekDouble> xinew)
{
    if (o->GetType() == curve)
    {
        for (int i = 0; i < xinew.num_elements(); i++)
        {
            all[i + 1] = xinew[i];
        }
    }
    else if (o->GetType() == surf)
    {
        for (int i = 0; i < xinew.num_elements(); i++)
        {
            all[i + 2] = xinew[i];
        }
    }
}

Array<OneD, NekDouble> OptiFace::Getxi()
{
    Array<OneD, NekDouble> ret(ni * 2);
    for (int i = np - ni; i < np; i++)
    {
        ret[2 * (i - np + ni) + 0] = uv[i][0];
        ret[2 * (i - np + ni) + 1] = uv[i][1];
    }
    return ret;
}

Array<OneD, NekDouble> OptiFace::Getli()
{
    Array<OneD, NekDouble> li(ni * 2);
    Array<OneD, NekDouble> bnds = s->GetBounds();
    for (int i = np - ni; i < np; i++)
    {
        li[2 * (i - np + ni) + 0] = bnds[0];
        li[2 * (i - np + ni) + 1] = bnds[2];
    }
    return li;
}

Array<OneD, NekDouble> OptiFace::Getui()
{
    Array<OneD, NekDouble> ui(ni * 2);
    Array<OneD, NekDouble> bnds = s->GetBounds();
    for (int i = np - ni; i < np; i++)
    {
        ui[2 * (i - np + ni) + 0] = bnds[1];
        ui[2 * (i - np + ni) + 1] = bnds[3];
    }
    return ui;
}

NekDouble OptiFace::F(Array<OneD, NekDouble> xitst)
{
    Array<OneD, Array<OneD, NekDouble> > val = uv;
    for (int i = np - ni; i < np; i++)
    {
        val[i][0] = xitst[(i - np + ni) * 2 + 0];
        val[i][1] = xitst[(i - np + ni) * 2 + 1];
    }

    NekDouble ret = 0.0;
    set<pair<int, int> >::iterator it;
    for (it = spring.begin(); it != spring.end(); it++)
    {
        Array<OneD, NekDouble> dis =
            Take(s->P(uv[(*it).first]), s->P(uv[(*it).second]));
        NekDouble norm = dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2];
        ret += norm / z[(*it)];
    }

    return ret;
}

DNekMat OptiFace::dF(Array<OneD, NekDouble> xitst)
{
    Array<OneD, Array<OneD, NekDouble> > val = uv;
    for (int i = np - ni; i < np; i++)
    {
        val[i][0] = xitst[(i - np + ni) * 2 + 0];
        val[i][1] = xitst[(i - np + ni) * 2 + 1];
    }

    DNekMat ret(ni * 2, 1, 0.0);

    vector<Array<OneD, NekDouble> > r, dru, drv;

    for (int i = 0; i < val.num_elements(); i++)
    {
        Array<OneD, NekDouble> ri(3), du(3), dv(3);
        Array<OneD, NekDouble> d1 = s->D1(val[i]);
        for (int i = 0; i < 3; i++)
        {
            ri[i] = d1[i];
            du[i] = d1[i + 3];
            dv[i] = d1[i + 6];
        }
        r.push_back(ri);
        dru.push_back(du);
        drv.push_back(dv);
    }

    for (int i = 0; i < ni * 2; i++) // for each of the varibles
    {
        int var = floor(i / 2) + np - ni;
        int tp  = i % 2; // 0 is a u 1 is a v

        set<pair<int, int> >::iterator it;
        for (it = spring.begin(); it != spring.end();
             it++) // for each of the springs
        {
            Array<OneD, NekDouble> dr1, dr2;

            if ((*it).first == var)
            {
                if (tp == 0)
                {
                    dr1 = dru[(*it).first];
                }
                else
                {
                    dr1 = drv[(*it).first];
                }
            }
            else
            {
                dr1 = Array<OneD, NekDouble>(3, 0.0);
            }

            if ((*it).second == var)
            {
                if (tp == 0)
                {
                    dr2 = dru[(*it).second];
                }
                else
                {
                    dr2 = drv[(*it).second];
                }
            }
            else
            {
                dr2 = Array<OneD, NekDouble>(3, 0.0);
            }

            ret(i, 0) +=
                2.0 / z[(*it)] *
                Dot(Take(r[(*it).first], r[(*it).second]), Take(dr1, dr2));
        }
    }
    return ret;
}

void OptiFace::Update(Array<OneD, NekDouble> xinew)
{
    for (int i = np - ni; i < np; i++)
    {
        uv[i][0] = xinew[2 * (i - np + ni) + 0];
        uv[i][1] = xinew[2 * (i - np + ni) + 1];
    }
}
}
}
