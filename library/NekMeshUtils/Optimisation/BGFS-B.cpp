////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.cpp
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
//  Description: surfacemeshing object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/Optimisation/BGFS-B.h>

#include <limits>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

    void BGFSUpdate(function<NekDouble(Array<OneD, NekDouble>, Array<OneD, NekDouble>, CADObjSharedPtr)> F,
                    function<DNekMat(Array<OneD, NekDouble>, Array<OneD, NekDouble>, CADObjSharedPtr)> Jac,
                    Array<OneD, NekDouble> &all, Array<OneD, NekDouble> z,
                    CADObjSharedPtr o,
                    DNekMat &J, DNekMat &B, DNekMat &H)
    {

        Array<OneD, NekDouble> xi;
        Array<OneD, NekDouble> gi;
        Array<OneD, NekDouble> ui;
        Array<OneD, NekDouble> li;
        Array<OneD, NekDouble> bnds;
        //reduce all data down to xi vector
        switch (o->GetType())
        {
            case curve:
                xi = Array<OneD, NekDouble>(all.num_elements()-2);
                gi = Array<OneD, NekDouble>(all.num_elements()-2);
                li = Array<OneD, NekDouble>(all.num_elements()-2);
                ui = Array<OneD, NekDouble>(all.num_elements()-2);
                bnds = boost::dynamic_pointer_cast<CADCurve>(o)->Bounds();
                for(int i = 1; i < all.num_elements() - 1; i++)
                {
                    li[i-1] = bnds[0];
                    ui[i-1] = bnds[1];
                    xi[i-1] = all[i];
                }
                break;

            case surf:
                xi = Array<OneD, NekDouble>(all.num_elements()-4);
                gi = Array<OneD, NekDouble>(all.num_elements()-4);
                li = Array<OneD, NekDouble>(all.num_elements()-4);
                ui = Array<OneD, NekDouble>(all.num_elements()-4);
                bnds = boost::dynamic_pointer_cast<CADSurf>(o)->GetBounds();
                for(int i = 2; i < all.num_elements() - 2; i++)
                {
                    if(i % 2 == 0)
                    {
                        li[i-2] = bnds[0];
                        ui[i-2] = bnds[1];
                    }
                    else
                    {
                        li[i-2] = bnds[2];
                        ui[i-2] = bnds[3];
                    }
                    xi[i-2] = all[i];
                }
                break;

            case vert:
                ASSERTL0(false,"Should not be able to pass vert");
        }

        set<int> Fset;
        Array<OneD, NekDouble> ti(xi.num_elements());
        for(int i = 0; i < ti.num_elements(); i++)
        {
            if(J(i,0) < 0)
            {
                ti[i] = (xi[i] - ui[i]) / J(i,0);
            }
            else if(J(i,0) > 0)
            {
                ti[i] = (xi[i] - li[i]) / J(i,0);
            }
            else
            {
                ti[i] = numeric_limits<double>::max();
            }
            if(ti[i] > 0)
            {
                Fset.insert(i);
            }
        }

        //intitialise d
        DNekMat d(xi.num_elements(),1);
        for(int i = 0; i < xi.num_elements(); i++)
        {
            if(fabs(ti[i]) < 1E-10)
            {
                d(i,0) = 0;
            }
            else
            {
                d(i,0) = -1.0 * J(i,0);
            }
        }

        Array<OneD, NekDouble> xci(xi.num_elements());

        for(int i = 0; i < xci.num_elements(); i++)
        {
            if(xi[i] + d(i,0) < li[i])
            {
                xci[i] = li[i];
                Fset.erase(i);
                cout << "hit bounded" << endl;
            }
            else
            {
                xci[i] = xi[i] + d(i,0);
            }

            if(xi[i] + d(i,0) > ui[i])
            {
                xci[i] = ui[i];
                Fset.erase(i);
                cout << "hit bounded" << endl;
            }
            else
            {
                xci[i] = xi[i] + d(i,0);
            }
        }

        DNekMat Z(xci.num_elements(), xci.num_elements() , 0.0);

        set<int>::iterator it;
        for(int i = 0; i < xci.num_elements(); i++)
        {
            it = Fset.find(i);
            if(it != Fset.end())
            {
                Z(i,i) = 1.0;
            }
        }

        DNekMat dx(xci.num_elements(), 1, 0.0);
        for(int i = 0; i < xci.num_elements(); i++)
        {
            dx(i,0) = xci[i] - xi[i];
        }

        DNekMat rg = Z * (J + B * dx);

        cout << endl <<rg << endl << endl;

        DNekMat du = -1.0 * H * rg;

        NekDouble alpha = 1.0;
        for(it = Fset.begin(); it != Fset.end(); it++)
        {
            int i = (*it);
            if(li[i] - xci[i] > alpha * du(i,0))
            {
                alpha = min(alpha, (li[i] - xci[i])/du(i,0));
            }
            else if(ui[i] - xci[i] < alpha * du(i,0))
            {
                alpha = min(alpha, (ui[i] - xci[i])/du(i,0));
            }
        }

        DNekMat grad = alpha * du;

        Array<OneD, NekDouble> dk(xci.num_elements()), xibar(xci.num_elements());
        for(int i = 0; i < xci.num_elements(); i++)
        {
            set<int>::iterator f = Fset.find(i);
            if(f != Fset.end())
            {
                xibar[i] = xci[i] + grad(i,0);
            }
            else
            {
                xibar[i] = xci[i];
            }
        }

        Vmath::Vsub(xci.num_elements(),&xibar[0],1,&xi[0],1,&dk[0],1);

        NekDouble c = 0.0;
        for(int i = 0; i < dk.num_elements(); i++)
        {
            c += 1E-4 * J(i,0) * dk[i];
        }

        cout << H << endl << endl;
        for(int i = 0; i < dk.num_elements(); i++)
        {
            cout << dk[i] << endl;
        }
        cout << endl;

        //this section needs a case evaluation for edges on faces
        NekDouble lam = 2.0;
        int iterct = 0;
        NekDouble fo = F(all,z,o);
        NekDouble fn;
        Array<OneD, NekDouble> tst(all.num_elements());
        cout << "begining line search: " << fo << endl;
        do
        {
            if(iterct > 20)
            {
                cout << "failed line search" << endl;
                abort();
            }
            iterct++;

            lam*=0.5;

            switch (o->GetType())
            {
                case curve:
                    tst[0] = all[0];
                    for(int i = 0; i < xi.num_elements(); i++)
                    {
                        tst[i+1] = xi[i] + lam * dk[i];
                    }
                    tst[tst.num_elements()-1] = all[tst.num_elements()-1];
                    break;

                case surf:
                    tst[0] = all[0];
                    tst[1] = all[1];
                    for(int i = 0; i < xi.num_elements(); i++)
                    {
                        tst[i+2] = xi[i] + lam * dk[i];
                    }
                    tst[tst.num_elements()-2] = all[tst.num_elements()-2];
                    tst[tst.num_elements()-1] = all[tst.num_elements()-1];
                    break;

                case vert:
                    ASSERTL0(false, "what");
            }

            fn = F(tst,z,o);
            cout << fn << endl;
        }while(fn > fo + c);
        cout << "lam " << lam << endl;

        //tst at this point is the new all vector
        //now need to update hessians
        DNekMat Jn = Jac(tst,z,o);
        DNekMat y = Jn - J;
        DNekMat yT = Jn - J;
        yT.Transpose();
        DNekMat s(dk.num_elements(),1,0.0);
        for(int i = 0; i < dk.num_elements(); i++)
        {
            s(i,0) = lam * dk[i];
        }
        DNekMat sT = s;
        sT.Transpose();

        DNekMat d1 = yT * s;
        DNekMat d2 = sT * B * s;
        DNekMat d3 = sT * y;
        DNekMat n1 = yT * H * y;

        NekDouble ynorm = 0.0;
        for(int i = 0; i < dk.num_elements(); i++)
        {
            ynorm += y(i,0)*y(i,0);
        }

        if(d3(0,0) > 2.2E-16 * ynorm)
        {
            B = B + y * yT * (1.0 / d1(0,0)) - B * s * sT * B * (1.0 / d2(0,0));
            H = H + (d3(0,0) + n1(0,0)) / d3(0,0) / d3(0,0) * s * sT - 1.0/d3(0,0) * (H * y * sT + s * yT * H);
        }

        J = Jn;
        all = tst;

    }

}
}
