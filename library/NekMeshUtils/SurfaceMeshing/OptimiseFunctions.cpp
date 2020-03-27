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
        case CADType::eCurve:
            xi = Array<OneD, NekDouble>(all.size() - 2);
            for (int i = 1; i < all.size() - 1; i++)
            {
                xi[i - 1] = all[i];
            }
            break;

        case CADType::eSurf:
            xi = Array<OneD, NekDouble>(all.size() - 4);
            for (int i = 2; i < all.size() - 2; i++)
            {
                xi[i - 2] = all[i];
            }
            break;

        case CADType::eVert:
        case CADType::eOther:
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
        case CADType::eCurve:
            li   = Array<OneD, NekDouble>(all.size() - 2);
            bnds = std::dynamic_pointer_cast<CADCurve>(o)->GetBounds();
            for (int i = 1; i < all.size() - 1; i++)
            {
                li[i - 1] = bnds[0];
            }
            break;

        case CADType::eSurf:
            li   = Array<OneD, NekDouble>(all.size() - 4);
            bnds = std::dynamic_pointer_cast<CADSurf>(o)->GetBounds();
            for (int i = 2; i < all.size() - 2; i++)
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

        case CADType::eVert:
        case CADType::eOther:
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
        case CADType::eCurve:
            ui   = Array<OneD, NekDouble>(all.size() - 2);
            bnds = std::dynamic_pointer_cast<CADCurve>(o)->GetBounds();
            for (int i = 1; i < all.size() - 1; i++)
            {
                ui[i - 1] = bnds[1];
            }
            break;

        case CADType::eSurf:
            ui   = Array<OneD, NekDouble>(all.size() - 4);
            bnds = std::dynamic_pointer_cast<CADSurf>(o)->GetBounds();
            for (int i = 2; i < all.size() - 2; i++)
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

        case CADType::eVert:
        case CADType::eOther:
            ASSERTL0(false, "Should not be able to pass vert");
    }
    return ui;
}

NekDouble OptiEdge::F(Array<OneD, NekDouble> xitst)
{
    Array<OneD, NekDouble> val(all.size());

    if (o->GetType() == CADType::eCurve)
    {
        val[0] = all[0];
        for (int i = 0; i < xitst.size(); i++)
        {
            val[i + 1] = xitst[i];
        }
        val[all.size() - 1] = all[all.size() - 1];
    }
    else if (o->GetType() == CADType::eSurf)
    {
        val[0] = all[0];
        val[1] = all[1];
        for (int i = 0; i < xitst.size(); i++)
        {
            val[i + 2] = xitst[i];
        }
        val[all.size() - 2] = all[all.size() - 2];
        val[all.size() - 1] = all[all.size() - 1];
    }

    NekDouble ret = 0.0;
    if (o->GetType() == CADType::eCurve)
    {
        CADCurveSharedPtr c = std::dynamic_pointer_cast<CADCurve>(o);

        for (int i = 0; i < all.size() - 1; i++)
        {
            Array<OneD, NekDouble> dis = Take(c->P(val[i + 1]), c->P(val[i]));
            NekDouble norm =
                dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2];
            ret += norm / (z[i + 1] - z[i]);
        }
    }
    else if (o->GetType() == CADType::eSurf)
    {
        CADSurfSharedPtr s = std::dynamic_pointer_cast<CADSurf>(o);
        // need to organise the val array
        Array<OneD, Array<OneD, NekDouble> > uv(val.size() / 2);
        for (int i = 0; i < val.size() / 2; i++)
        {
            uv[i]    = Array<OneD, NekDouble>(2);
            uv[i][0] = val[i * 2 + 0];
            uv[i][1] = val[i * 2 + 1];
        }
        for (int i = 0; i < uv.size() - 1; i++)
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
    Array<OneD, NekDouble> val(all.size());

    if (o->GetType() == CADType::eCurve)
    {
        val[0] = all[0];
        for (int i = 0; i < xitst.size(); i++)
        {
            val[i + 1] = xitst[i];
        }
        val[all.size() - 1] = all[all.size() - 1];
    }
    else if (o->GetType() == CADType::eSurf)
    {
        val[0] = all[0];
        val[1] = all[1];
        for (int i = 0; i < xitst.size(); i++)
        {
            val[i + 2] = xitst[i];
        }
        val[all.size() - 2] = all[all.size() - 2];
        val[all.size() - 1] = all[all.size() - 1];
    }

    DNekMat ret;

    if (o->GetType() == CADType::eCurve)
    {
        CADCurveSharedPtr c = std::dynamic_pointer_cast<CADCurve>(o);
        vector<Array<OneD, NekDouble> > r;
        vector<Array<OneD, NekDouble> > dr;

        for (int i = 0; i < all.size(); i++)
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

        DNekMat J(all.size() - 2, 1, 0.0);
        for (int i = 0; i < all.size() - 2; i++)
        {
            J(i, 0) =
                2.0 / (z[i + 1] - z[i]) * Dot(dr[i + 1], Take(r[i + 1], r[i])) -
                2.0 / (z[i + 2] - z[i + 1]) *
                    Dot(dr[i + 1], Take(r[i + 2], r[i + 1]));
        }

        ret = J;
    }
    else if (o->GetType() == CADType::eSurf)
    {
        CADSurfSharedPtr s = std::dynamic_pointer_cast<CADSurf>(o);
        // need to organise the all array
        Array<OneD, Array<OneD, NekDouble> > uv(val.size() / 2);
        for (int i = 0; i < val.size() / 2; i++)
        {
            uv[i]    = Array<OneD, NekDouble>(2);
            uv[i][0] = val[i * 2 + 0];
            uv[i][1] = val[i * 2 + 1];
        }

        vector<Array<OneD, NekDouble> > r;
        vector<Array<OneD, NekDouble> > dru, drv;
        for (int i = 0; i < uv.size(); i++)
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

        DNekMat J(2 * (uv.size() - 2), 1, 0.0);
        for (int i = 0; i < uv.size() - 2; i++)
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
    if (o->GetType() == CADType::eCurve)
    {
        for (int i = 0; i < xinew.size(); i++)
        {
            all[i + 1] = xinew[i];
        }
    }
    else if (o->GetType() == CADType::eSurf)
    {
        for (int i = 0; i < xinew.size(); i++)
        {
            all[i + 2] = xinew[i];
        }
    }
}
}
}
