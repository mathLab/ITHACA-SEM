///////////////////////////////////////////////////////////////////////////////
//
// File ComputeCriticalLayer.cpp
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
// Description: Compute location of critical layer from streak file
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>

#include <MultiRegions/ExpList.h>

using namespace std;
using namespace Nektar;

void Computestreakpositions(MultiRegions::ExpListSharedPtr &streak,
                            Array<OneD, NekDouble> &xc,
                            Array<OneD, NekDouble> &yc,
                            NekDouble cr,
                            NekDouble trans)
{
    int i;
    int npts = xc.size();

    int nq = streak->GetTotPoints();
    Array<OneD, NekDouble> derstreak(nq);
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    streak->GetCoords(x,y);

    streak->BwdTrans(streak->GetCoeffs(),streak->UpdatePhys());
    streak->PhysDeriv(MultiRegions::eY, streak->GetPhys(), derstreak);

    // set intiial xc to be equispaced over mesh and yc to be zero
    NekDouble x_max = Vmath::Vmax(nq,x,1);
    NekDouble x_min = Vmath::Vmin(nq,x,1);

    for(i = 0; i < npts; ++i)
    {
        xc[i]  = x_min + (x_max - x_min)*i/((NekDouble)(npts-1));
        yc[i] = 0.0;
    }


    int elmtid, offset,cnt;
    NekDouble U,dU;
    NekDouble F;
    NekDouble ConvTol  = 1e-9;
    NekDouble CoordTol = 1e-5;
    int maxiter = 100;
    Array<OneD, NekDouble> coord(2);

    // Do Newton iteration on y direction
    cerr << "[";
    for(int e=0; e<npts; e++)
    {
        coord[0] = xc[e];
        coord[1] = yc[e];

        if(!(e%10))
        {
            cerr << ".";
        }

        F = 1000;
        cnt = 0;
        while((abs(F)> ConvTol)&&(cnt < maxiter))
        {
            elmtid = streak->GetExpIndex(coord,CoordTol);
            offset = streak->GetPhys_Offset(elmtid);

            U  = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);
            dU = streak->GetExp(elmtid)->PhysEvaluate(coord, derstreak + offset);

            coord[1] = coord[1] - (U-cr)/dU;

            F = U-cr;
            cnt++;
        }
        ASSERTL0(cnt < maxiter, "Failed to converge Newton iteration");

        yc[e] = coord[1];
    }
    cerr << "]" << endl;

    if(trans != NekConstants::kNekUnsetDouble)
    {
        // output to interface file
        FILE *fp = fopen("interfacedat.geo","w");

        NekDouble y_max = Vmath::Vmax(nq,y,1);
        NekDouble y_min = Vmath::Vmin(nq,y,1);

        cnt = 1;
        fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0}; \n",
                cnt++,x_min,y_min);
        fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0}; \n",
                cnt++,x_max,y_min);
        fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0}; \n",
                cnt++,x_max,y_max);
        fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0}; \n",
                cnt++,x_min,y_max);

        for(i = 0; i < npts; ++i)
        {
            fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0};  \n",
                    cnt++,xc[i],yc[i]);
        }

        fclose(fp);


        // output to interface_up file as bend of vertical shift and 45 degrees shift
        fp = fopen("interfacedat_up.geo","w");


        NekDouble nx,ny,norm;

        fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0};  \n",cnt++,xc[0],yc[0]+trans);

        for(i = 1; i < npts-1; ++i)
        {
            norm = sqrt((xc[i+1]-xc[i-1])*(xc[i+1]-xc[i-1])+(yc[i+1]-yc[i-1])*(yc[i+1]-yc[i-1]));
            nx = (yc[i-1]-yc[i+1])/norm;
            ny = (xc[i+1]-xc[i-1])/norm;

            fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0};  \n",
                    cnt++,xc[i]+nx*trans,yc[i]+ny*trans);
        }

        fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0};  \n",cnt++,xc[npts-1],yc[npts-1]+trans);


        // output to interface_up file as bend of vertical shift and 45 degrees shift
        fp = fopen("interfacedat_dn.geo","w");

        trans = -trans;

        fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0};  \n",cnt++,xc[0],yc[0]+trans);

        for(i = 1; i < npts-1; ++i)
        {
            norm = sqrt((xc[i+1]-xc[i-1])*(xc[i+1]-xc[i-1])+(yc[i+1]-yc[i-1])*(yc[i+1]-yc[i-1]));
            nx = (yc[i-1]-yc[i+1])/norm;
            ny = (xc[i+1]-xc[i-1])/norm;

            fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0};  \n",
                    cnt++,xc[i]+nx*trans,yc[i]+ny*trans);
        }

        fprintf(fp,"Point(%d)={%12.10lf,%12.10lf,0,1.0};  \n",cnt++,xc[npts-1],yc[npts-1]+trans);
    }

}
