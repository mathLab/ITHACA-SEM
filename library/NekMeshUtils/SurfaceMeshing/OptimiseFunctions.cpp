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

#include <NekMeshUtils/SurfaceMeshing/SurfaceMesh.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

NekDouble Dot(Array<OneD, NekDouble> a, Array<OneD, NekDouble> b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
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
    ret[0] = a[0] *t;
    ret[1] = a[1] *t;
    ret[2] = a[2] *t;
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

NekDouble SurfaceMesh::EdgeF(Array<OneD, NekDouble> all, Array<OneD, NekDouble> z,
                             CADObjSharedPtr o)
{
    NekDouble ret = 0.0;
    if(o->GetType() == curve)
    {
        CADCurveSharedPtr c = boost::dynamic_pointer_cast<CADCurve>(o);

        for(int i = 0; i < all.num_elements() - 1; i++)
        {
            Array<OneD, NekDouble> dis = Take(c->P(all[i+1]), c->P(all[i]));
            NekDouble norm = dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2];
            ret += norm/(z[i+1] - z[i]);
        }
    }
    else if(o->GetType() == surf)
    {
        CADSurfSharedPtr s = boost::dynamic_pointer_cast<CADSurf>(o);
        //need to organise the all array
        Array<OneD, Array<OneD, NekDouble> > uv(all.num_elements()/2);
        for(int i = 0; i < all.num_elements()/2; i++)
        {
            uv[i] = Array<OneD, NekDouble>(2);
            uv[i][0] = all[i*2+0];
            uv[i][1] = all[i*2+1];
        }
        for(int i = 0; i < uv.num_elements() - 1; i++)
        {
            Array<OneD, NekDouble> dis = Take(s->P(uv[i+1]), s->P(uv[i]));
            NekDouble norm = dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2];
            ret += norm/(z[i+1] - z[i]);
        }
    }
    return ret;
}

DNekMat SurfaceMesh::EdgeGrad(Array<OneD, NekDouble> all, Array<OneD, NekDouble> z,
                     CADObjSharedPtr o)
{
    if(o->GetType() == curve)
    {
        CADCurveSharedPtr c = boost::dynamic_pointer_cast<CADCurve>(o);
        vector<Array<OneD, NekDouble> > r;
        vector<Array<OneD, NekDouble> > dr;

        for(int i = 0; i < all.num_elements(); i++)
        {
            Array<OneD, NekDouble> ri(3), dri(3);
            Array<OneD, NekDouble> d2 = c->D2(all[i]);
            for(int j = 0; j < 3; j++)
            {
                ri[j] = d2[j];
                dri[j] = d2[j+3];
            }
            r.push_back(ri);
            dr.push_back(dri);
        }

        DNekMat J(all.num_elements() - 2 , 1, 0.0);
        for(int i = 0; i < all.num_elements() - 2; i++)
        {
            J(i,0) = 2.0/(z[i+1] - z[i]) * Dot(dr[i+1],Take(r[i+1],r[i])) -
                     2.0/(z[i+2] - z[i+1]) * Dot(dr[i+1],Take(r[i+2],r[i+1]));
        }

        return J;
    }
    else if(o->GetType() == surf)
    {
        CADSurfSharedPtr s = boost::dynamic_pointer_cast<CADSurf>(o);
        //need to organise the all array
        Array<OneD, Array<OneD, NekDouble> > uv(all.num_elements()/2);
        for(int i = 0; i < all.num_elements()/2; i++)
        {
            uv[i] = Array<OneD, NekDouble>(2);
            uv[i][0] = all[i*2+0];
            uv[i][1] = all[i*2+1];
        }

        vector<Array<OneD, NekDouble> > r;
        vector<Array<OneD, NekDouble> > dru, drv;
        for(int i = 0; i < uv.num_elements(); i++)
        {
            Array<OneD, NekDouble> ri(3), drui(3), drvi(3);
            Array<OneD, NekDouble> d2 = s->D1(uv[i]);
            for(int j = 0; j < 3; j++)
            {
                ri[j] = d2[j];
                drui[j] = d2[j+3];
                drvi[j] = d2[j+6];
            }
            r.push_back(ri);
            dru.push_back(drui);
            drv.push_back(drvi);
        }

        DNekMat J(2*(uv.num_elements() - 2), 1, 0.0);
        for(int i = 0; i < uv.num_elements() - 2; i++)
        {
            J(2*i+0,0) = 2.0/(z[i+1]-z[i]) * Dot(dru[i+1], Take(r[i+1], r[i])) +
                         2.0/(z[i+2]-z[i+1]) * Dot(dru[i+1], Take(r[i+1], r[i+2]));

            J(2*i+1,0) = 2.0/(z[i+1]-z[i]) * Dot(drv[i+1], Take(r[i+1], r[i])) +
                         2.0/(z[i+2]-z[i+1]) * Dot(drv[i+1], Take(r[i+1], r[i+2]));
        }

        return J;
    }
}

NekDouble SurfaceMesh::FaceEdgeF(Array<OneD, Array<OneD, NekDouble> > uv, Array<OneD, NekDouble> z,
                                 CADSurfSharedPtr s)
{
    NekDouble ret = 0.0;
    for(int i = 0; i < uv.num_elements() - 1; i++)
    {
        Array<OneD, NekDouble> dis = Take(s->P(uv[i+1]), s->P(uv[i]));
        NekDouble d = dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2];
        ret += d/(z[i+1] - z[i]);
    }

    return ret;
}

void SurfaceMesh::FaceEdgeJac(Array<OneD, Array<OneD, NekDouble> > uv,
                              Array<OneD, NekDouble> z, CADSurfSharedPtr s,
                              DNekMat &Jac)
{
    vector<Array<OneD, NekDouble> > r;
    vector<Array<OneD, NekDouble> > dru, drv;
    for(int i = 0; i < uv.num_elements(); i++)
    {
        Array<OneD, NekDouble> ri(3), drui(3), drvi(3);
        Array<OneD, NekDouble> d2 = s->D1(uv[i]);
        for(int j = 0; j < 3; j++)
        {
            ri[j] = d2[j];
            drui[j] = d2[j+3];
            drvi[j] = d2[j+6];
        }
        r.push_back(ri);
        dru.push_back(drui);
        drv.push_back(drvi);
    }

    DNekMat J(2*(uv.num_elements() - 2), 1, 0.0);
    for(int i = 0; i < uv.num_elements() - 2; i++)
    {
        J(2*i+0,0) = 2.0/(z[i+1]-z[i]) * Dot(dru[i+1], Take(r[i+1], r[i])) +
                     2.0/(z[i+2]-z[i+1]) * Dot(dru[i+1], Take(r[i+1], r[i+2]));

        J(2*i+1,0) = 2.0/(z[i+1]-z[i]) * Dot(drv[i+1], Take(r[i+1], r[i])) +
                     2.0/(z[i+2]-z[i+1]) * Dot(drv[i+1], Take(r[i+1], r[i+2]));
    }

    Jac = J;
}

/*void SurfaceMesh::EdgeOnFaceUpdate(Array<OneD, Array<OneD, NekDouble> > uv,
                                                       Array<OneD, NekDouble> x,
                                                       Array<OneD, NekDouble> gll,
                                                       CADSurfSharedPtr s,
                                                       DNekMat &B, DNekMat &J)
{
    for(int i = 0; i < x.num_elements(); i++)
    {
        J(i,0) *= -1.0; //minimise not maximise
    }
    DNekMat P = B*J;

    NekDouble pnorm = 0;
    for(int i = 0; i < x.num_elements(); i++)
    {
        pnorm += P(i,0)*P(i,0);
    }
    pnorm = sqrt(pnorm);

    if(pnorm < 1E-6)
    {
        B = DNekMat(x.num_elements(), x.num_elements(), 0.0);
        for(int k = 0; k < x.num_elements(); k++)
        {
            B(k,k) = 1.0;
        }
        return;
    }

    cout << endl << endl << P << endl;

    Array<OneD, NekDouble> bnds = s->GetBounds();

    NekDouble F = FaceEdgeF(uv, gll, s);
    NekDouble f = F*2.0;

    NekDouble alpha = 2.0;
    int ct = 0;
    //perform line search to obtain alpha
    while(f > F)
    {
        if(ct > 15)
            break;
        ct++;

        alpha*=0.5;
        Array<OneD, Array<OneD, NekDouble> > tst(uv.num_elements());
        int i;
        tst[0] = uv[0];
        for(i = 0; i < tst.num_elements() -2; i++)
        {
            tst[i+1] = Array<OneD, NekDouble>(2);
            if(uv[i+1][0] + alpha*P(i*2+0,0) < bnds[0])
            {
                tst[i+1][0] = bnds[0];
            }
            else if(uv[i+1][0] + alpha*P(i*2+0,0) > bnds[1])
            {
                tst[i+1][0] = bnds[1];
            }
            else
            {
                tst[i+1][0] = uv[i+1][0] + alpha*P(i*2+1,0);
            }

            if(uv[i+1][1] + alpha*P(i*2+1,0) < bnds[2])
            {
                tst[i+1][1] = bnds[2];
            }
            else if(uv[i+1][1] + alpha*P(i*2+1,0) > bnds[3])
            {
                tst[i+1][1] = bnds[3];
            }
            else
            {
                tst[i+1][1] = uv[i+1][1] + alpha*P(i*2+1,0);
            }
        }
        tst[i+1] = uv[i+1];

        f = FaceEdgeF(tst, gll, s);
    }
    cout << alpha << endl;

    //for(int i = 1; i < ti.num_elements() - 1; i++)
    //{
    //    ti[i] += alpha*P(i-1,0);
    //}

    //now need to update the inverse hessian B
    //update J
    DNekMat Jn;
    CurveEdgeJac(ti, gll, c, Jn);

    DNekMat y(x.num_elements(),1);
    DNekMat s(x.num_elements(),1);
    DNekMat yt(x.num_elements(),1);
    DNekMat st(x.num_elements(),1);

    for(int i = 0; i < x.num_elements(); i++)
    {
        s(i,0) = alpha*P(i,0);
        y(i,0) = Jn(i,0) - J(i,0);
        st(i,0) = alpha*P(i,0);
        yt(i,0) = Jn(i,0) - J(i,0);
    }

    st.Transpose();
    yt.Transpose();
    //get scalars

    DNekMat s1 = yt * B * y;
    DNekMat s2 = st * y;

    DNekMat sst = s * st;

    NekDouble mul = (s1(0,0) + s2(0,0))/s2(0,0)/s2(0,0);
    for(int i = 0; i < sst.GetRows(); i++)
    {
        for(int j = 0; j < sst.GetColumns(); j++)
        {
            sst(i,j) *= mul;
        }
    }

    DNekMat r1 = B * y * st;
    DNekMat r2 = s * yt * B;

    for(int i = 0; i < r1.GetRows(); i++)
    {
        for(int j = 0; j < r1.GetColumns(); j++)
        {
            r1(i,j) = (r1(i,j) + r2(i,j)) / s2(0,0);
        }
    }

    for(int i = 0; i < r1.GetRows(); i++)
    {
        for(int j = 0; j < r1.GetColumns(); j++)
        {
            B(i,j) = B(i,j) + sst(i,j) - r1(i,j);
        }
    }

    J = Jn;

}*/

void SurfaceMesh::FaceFaceJac(int p, Array<OneD, Array<OneD, NekDouble> > uv,
                              map<int, Array<OneD, NekDouble> > z,
                              map<int, vector<int> > n,
                              CADSurfSharedPtr s,
                              DNekMat &Jac, DNekMat &Hes)
{
    DNekMat J(2,1,0.0);
    DNekMat H(2,2,0.0);

    vector<Array<OneD, NekDouble> > r;
    Array<OneD, NekDouble> d2;

    for(int i = 0; i < uv.num_elements(); i++)
    {
        Array<OneD, NekDouble> d2 = s->P(uv[i]);
        r.push_back(d2);
    }
    d2 = s->D2(uv[p]);

    Array<OneD, NekDouble> du(3), dv(3), d2u(3), d2v(3), d2uv(3);
    for(int i = 0; i < 3; i++)
    {
        du[i] = d2[i+3];
        dv[i] = d2[i+6];
        d2u[i] = d2[i+9];
        d2v[i] = d2[i+12];
        d2uv[i] = d2[i+15];
    }

    Array<OneD, NekDouble> zero(3,0.0);

    for(int k = 0; k < 6; k++) //over number of directions of springs
    {

        J(0,0) += 2.0/ z[p][k] * Dot(Take(zero, du), Take(r[n[p][k]], r[p]));
        J(1,0) += 2.0/ z[p][k] * Dot(Take(zero, dv), Take(r[n[p][k]], r[p]));
    }


    for(int k = 0; k < 6; k++)
    {
        H(0,0) += 2.0 / z[p][k] * (Dot(Take(zero,du), Take(zero,du)) +
                                   Dot(Take(r[n[p][k]], r[p]), Take(zero, d2u)));

        H(1,1) += 2.0 / z[p][k] * (Dot(Take(zero,dv), Take(zero,dv)) +
                                   Dot(Take(r[n[p][k]], r[p]), Take(zero, d2v)));

        H(1,0) += 2.0 / z[p][k] * (Dot(Take(zero,dv), Take(zero,du)) +
                                   Dot(Take(r[n[p][k]], r[p]), Take(zero, d2uv)));

        H(0,1) += H(1,0);
    }


    Jac = J;
    Hes = H;
}

}
}
