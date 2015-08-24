////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
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

#include <MeshUtils/SurfaceMeshing/SurfaceMeshing.h>

using namespace std;
namespace Nektar{
namespace MeshUtils{

NekDouble SurfaceMeshing::BrentOpti(NekDouble ax, NekDouble bx,
                                    NekDouble cx, NekDouble &fx,
                                    NekDouble tol, int surf,
                                    Array<OneD, NekDouble> uvi,
                                    Array<OneD, NekDouble> df,
                                    Array<OneD, NekDouble> bounds,
                                    vector<Array<OneD,NekDouble> > bcs,
                                    NekDouble (SurfaceMeshing::*GetF)(
                                    NekDouble, NekDouble,
                                    vector<Array<OneD,NekDouble> >, int))
{
    //enter brent algoithm from book
    NekDouble a,b;
    NekDouble e1 = 0.0;
    NekDouble xmin=0.0;
    a = ax < cx ? ax : cx; b = ax > cx ? ax : cx;
    NekDouble x,w,v;
    x = w = v = bx;
    NekDouble fw,fv;
    fw = fv = fx;

    NekDouble xm,tol2,tol1,r,q,p,etemp,d,fu,u;
    NekDouble zeps = 1E-10;
    for(int it = 0; it < 100; it++)
    {
        xm = 0.5*(a+b);
        tol2 = 2.0*(tol1=tol*fabs(x)+zeps);
        if(fabs(x-xm) <= (tol2-0.5*(b-a)))
        {
            break;
        }
        if(fabs(e1) > tol1)
        {
            r = (x-w)*(fx-fv);
            q = (x-v)*(fx-fw);
            p = (x-v)*q-(x-w)*r;
            q = 2.0*(q-r);
            if(q>0.0) p = -p;
            q = fabs(q);
            etemp = e1;
            e1 = d;
            if(fabs(p) > fabs(0.5*q*etemp) ||
               p <= q*(a-x) || p >= q*(b-x))
                    d=0.3819660*(e1=(x>= xm ? a-x : b-x));
            else
            {
                d=p/q;
                u=x+d;
                if(u-a < tol2 || b-u < tol2)
                    d = tol1*(xm-x)/fabs(xm-x);
            }
        }
        else
        {
            d=0.3819660*(e1=(x>= xm ? a-x : b-x));
        }

        u = fabs(d) >= tol1 ? x+d : x+tol1*d/fabs(d);
        if(uvi[0]+df[0]*u < bounds[0] ||
           uvi[0]+df[0]*u > bounds[1] ||
           uvi[1]+df[1]*u < bounds[2] ||
           uvi[1]+df[1]*u > bounds[3])
        {
            break;
        }
        fu = (this->*GetF)(uvi[0]+df[0]*u,uvi[1]+df[1]*u,bcs,surf);
        if(fu <= fx)
        {
            if(u>=x) a=x; else b=x;
            v=w;
            w=x;
            x=u;
            fv=fw;
            fw=fx;
            fx=fu;
        }
        else
        {
            if(u < x) a=u; else b=u;
            if(fu<=fw || w ==x)
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else if(fu <= fv || v==x || v==w)
            {
                v=u;
                fv=fu;
            }
        }
    }
    xmin = x;
    return xmin;
}

void SurfaceMeshing::Find1DBounds(NekDouble &a, NekDouble &b,
                                  Array<OneD, NekDouble> uvi,
                                  Array<OneD, NekDouble> df,
                                  Array<OneD, NekDouble> bounds)
{
    bool aset = false; bool bset = false;
    NekDouble K, L, inter;
    //want a to be negative, against the gradient, but properly bounded!!

    //check edges of bounding box one by one for intersect, because paralell lines some cases can be ingnored
    //line 1 left edge;
    if(!(fabs(df[0]) < 1E-10)) //wouldnt exist on this edge
    {
        K = (bounds[0] - uvi[0]) / df[0];
        L = df[1] * K + uvi[1] - bounds[2];
        inter = bounds[2] + L;
        if(!(inter < bounds[2] || inter > bounds[3]))
        {
            //hit
            if(K < 0)
            {
                if(aset)
                {
                    cout << "error, should not be set" << endl;
                    cout << a << " " << b << " " << df[0] << " " << df[1] << endl;
                    exit(-1);
                }
                else
                {
                    a = K;
                    aset = true;
                }
            }
            else
            {
                if(bset)
                {
                    cout << "error, should not be set" << endl;
                    cout << a << " " << b << " " << df[0] << " " << df[1] << endl;
                    exit(-1);
                }
                else
                {
                    b = K;
                    bset = true;
                }
            }
        }
    }
    //line 2 right edge;
    if(!(fabs(df[0]) < 1E-10)) //wouldnt exist on this edge
    {
        K = (bounds[1] - uvi[0]) / df[0];
        L = df[1] * K + uvi[1] - bounds[2];
        inter = bounds[2] + L;
        if(!(inter < bounds[2] || inter > bounds[3]))
        {
            //hit
            if(K < 0)
            {
                if(aset)
                {
                    cout << "error, should not be set" << endl;
                    cout << a << " " << b << " " << df[0] << " " << df[1] << endl;
                    exit(-1);
                }
                else
                {
                    a = K;
                    aset = true;
                }
            }
            else
            {
                if(bset)
                {
                    cout << "error, should not be set" << endl;
                    cout << a << " " << b << " " << df[0] << " " << df[1] << endl;
                    exit(-1);
                }
                else
                {
                    b = K;
                    bset = true;
                }
            }
        }
    }
    //line 3 bottom edge;
    if(!(fabs(df[1]) < 1E-10)) //wouldnt exist on this edge
    {
        K = (bounds[2] - uvi[1]) / df[1];
        L = df[0] * K + uvi[0] - bounds[0];
        inter = bounds[0] + L;
        if(!(inter < bounds[0] || inter > bounds[1]))
        {
            //hit
            if(K < 0)
            {
                if(aset)
                {
                    cout << "error, should not be set" << endl;
                    cout << a << " " << b << " " << df[0] << " " << df[1] << endl;
                    exit(-1);
                }
                else
                {
                    a = K;
                    aset = true;
                }
            }
            else
            {
                if(bset)
                {
                    cout << "error, should not be set" << endl;
                    cout << a << " " << b << " " << df[0] << " " << df[1] << endl;
                    exit(-1);
                }
                else
                {
                    b = K;
                    bset = true;
                }
            }
        }
    }
    //line 4 top edge;
    if(!(fabs(df[1]) < 1E-10)) //wouldnt exist on this edge
    {
        K = (bounds[3] - uvi[1]) / df[1];
        L = df[0] * K + uvi[0] - bounds[0];
        inter = bounds[0] + L;
        if(!(inter < bounds[0] || inter > bounds[1]))
        {
            //hit
            if(K < 0)
            {
                if(aset)
                {
                    cout << "error, should not be set" << endl;
                    cout << a << " " << b << " " << df[0] << " " << df[1] << endl;
                    exit(-1);
                }
                else
                {
                    a = K;
                    aset = true;
                }
            }
            else
            {
                if(bset)
                {
                    cout << "error, should not be set" << endl;
                    cout << a << " " << b << " " << df[0] << " " << df[1] << endl;
                    exit(-1);
                }
                else
                {
                    b = K;
                    bset = true;
                }
            }
        }
    }
}

}
}
