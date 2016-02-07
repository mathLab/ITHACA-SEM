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
                    Array<OneD, NekDouble> all, Array<OneD, NekDouble> z,
                    CADObjSharedPtr o,
                    DNekMat J, DNekMat &B)
    {
        cout << endl;

        Array<OneD, NekDouble> xi;
        Array<OneD, NekDouble> gi;
        Array<OneD, NekDouble> ui;
        Array<OneD, NekDouble> li;
        //reduce all data down to xi vector
        switch (o->GetType())
        {
            case curve:
                xi = Array<OneD, NekDouble>(all.num_elements()-2);
                gi = Array<OneD, NekDouble>(all.num_elements()-2);
                li = Array<OneD, NekDouble>(all.num_elements()-2);
                ui = Array<OneD, NekDouble>(all.num_elements()-2);
                Array<OneD, NekDouble> bnds = boost::dynamic_pointer_cast<CADCurve>(o)->Bounds();
                for(int i = 1; i < all.num_elements() - 1; i++)
                {
                    li[i-1] = bnds[0];
                    ui[i-1] = bnds[1];
                    xi[i-1] = all[i];
                }
                break;
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


        cout << endl;
        for(int i = 0; i < xci.num_elements(); i++)
        {
            cout << xci[i] << endl;
        }

        DNekMat ZT(Fset.size(), xci.num_elements(), 0.0);
        DNekMat Z(Fset.size(), xci.num_elements(), 0.0);

        set<int>::iterator it;
        int ct = 0;
        for(it = Fset.begin(); it != Fset.end(); it++)
        {
            ZT(ct,(*it)) = 1.0;
            Z(ct,(*it)) = 1.0;
            ct++;
        }
        ZT.Transpose();

        cout << endl << Z << endl << endl << ZT << endl << endl;

        DNekMat rB = ZT * B * Z;

        DNekMat dx(xci.num_elements(), 1, 0.0);
        for(int i = 0; i < xci.num_elements(); i++)
        {
            dx(i,0) = xci[i] - xi[i];
        }

        DNekMat rg = ZT * (J + B * dx);

        exit(-1);


    }

}
}
