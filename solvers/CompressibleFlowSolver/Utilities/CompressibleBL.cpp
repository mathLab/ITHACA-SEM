///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleBL.cpp
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
// Description: Generate the compressible boundary layer similarity solution
// on a 2D/3D mesh.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/DisContField2D.h>
#include <MultiRegions/DisContField3D.h>

#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion.h>

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <SpatialDomains/MeshGraph.h>

#include <SolverUtils/SolverUtilsDeclspec.h>

using namespace std;
using namespace Nektar;

NekDouble m_Re;
NekDouble m_Mach;
NekDouble L;
NekDouble m_Tinf;
NekDouble m_Suth;
NekDouble m_Tw;
NekDouble m_Twall;
NekDouble m_Gamma;
NekDouble m_Pr;
NekDouble m_long;
NekDouble m_uInf;
NekDouble m_rhoInf;
NekDouble m_R;
NekDouble m_vInf;
NekDouble m_mu;
NekDouble m_To = 273.11;

const int m_xpoints = 1000001;

const NekDouble Nvisc = 1;
const NekDouble Omega = 1;
const NekDouble etamax = 10.0;
const NekDouble errtol = 1e-5;


/**
 * Calculate the compressible boundary layer using the similarity solution
 */
void COMPBL(Array<OneD, NekDouble> v,
            Array<OneD, NekDouble> dv)
{
    NekDouble c, dcdg, cp;

    if (Nvisc == 1)
    {
        c = sqrt(v[3]) * (1.0 + m_Suth) / (v[3] + m_Suth);
        dcdg = 1.0 / (2. * sqrt(v[3])) - sqrt(v[3]) / (v[3]+m_Suth);
        dcdg = dcdg * (1.0 + m_Suth) / (v[3] + m_Suth);
        cp = dcdg * v[4];
    }
    if (Nvisc == 2)
    {
        c    = pow(v[3], (Omega-1.0));
        dcdg = (Omega - 1.0) * pow(v[3], (Omega - 2.0));
        cp   = dcdg * v[4];
    }
    if (Nvisc == 3)
    {
        c  = sqrt(m_Twall) * (1.0 + m_Suth) / (m_Suth + m_Twall);
        cp = 0.0;
    }

    dv[0] = v[1];
    dv[1] = v[2];
    dv[2] = - v[2] * (cp + v[0]) / c;
    dv[3] = v[4];
    dv[4] = - v[4] * (cp + m_Pr * v[0]) / c -
        m_Pr * (m_Gamma - 1.0) * pow(m_Mach, 2.0) *
        pow(v[2], 2);
}

/**
 * Perform the RK4 integration
 */
void RK4(Array<OneD, NekDouble> y,
         Array<OneD, NekDouble> dydx,
         int                    n,
         NekDouble              x,
         NekDouble              h,
         Array<OneD, NekDouble> yout)
{
    boost::ignore_unused(x);

    int nmax = 5;

    Array<OneD, NekDouble> yt (nmax, 0.0);
    Array<OneD, NekDouble> dyt(nmax, 0.0);
    Array<OneD, NekDouble> dym(nmax, 0.0);
    NekDouble hh = h * 0.5;
    NekDouble h6 = h / 6;

    for (int i = 0; i < n ; i++)
    {
        yt[i] = y[i] + hh * dydx[i];
    }

    COMPBL(yt, dyt);

    for (int i = 0; i < n; i++)
    {
        yt[i] = y[i] + hh * dyt[i];
    }

    COMPBL(yt, dym);

    for (int i = 0; i < n; i++)
    {
        yt[i]  = y[i] + h * dym[i];
        dym[i] = dyt[i] + dym[i];
    }

    COMPBL(yt, dyt);

    for (int i = 0; i < n; i++)
    {
        yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2 * dym[i]);
    }
}


/**
 * Calculate initial guess for RK4
 */
void RKDUMB(Array<OneD, NekDouble>               vstart,
            int                                  nvar,
            NekDouble                            x1,
            NekDouble                            x2,
            int                                  m_xpoints,
            Array<OneD, NekDouble>               xx,
            Array<OneD, Array<OneD, NekDouble> > y)
{
    int nmax = 5;
    NekDouble x, h;
    Array<OneD, NekDouble> v (nmax, 0.0);
    Array<OneD, NekDouble> dv(nmax, 0.0);

    for (int i = 0; i < nvar; i++)
    {
        v[i]    = vstart[i];
        y[i][0] = v[i];
    }

    xx[0] = x1;
    x     = x1;
    h     = (x2-x1) / m_xpoints;

    for (int k = 0; k < m_xpoints; k++)
    {
        COMPBL(v, dv);
        RK4   (v, dv, nvar, x, h, v);

        if (x + h == x)
        {
            cout << "bug" << endl;
        }

        x = x + h;
        xx[k+1] = x;

        for (int i = 0; i < nvar; i++)
        {
            y[i][k+1] = v[i];
        }
    }
}

/**
 * Create the output file
 */
void OUTPUT(int                                  m_xpoints,
            Array <OneD, NekDouble >             xx,
            Array<OneD, Array<OneD, NekDouble> > ff,
            int                                  nQuadraturePts,
            Array <OneD, NekDouble >             x_QuadraturePts,
            Array <OneD, NekDouble >             y_QuadraturePts,
            Array <OneD, NekDouble >             u_QuadraturePts,
            Array <OneD, NekDouble >             v_QuadraturePts,
            Array <OneD, NekDouble >             rho_QuadraturePts,
            Array <OneD, NekDouble >             T_QuadraturePts)
{
    Array <OneD, NekDouble > z       (m_xpoints, 0.0);
    Array <OneD, NekDouble > v       (m_xpoints, 0.0);
    Array <OneD, NekDouble > dv      (m_xpoints, 0.0);
    Array <OneD, NekDouble > u       (m_xpoints, 0.0);
    Array <OneD, NekDouble > t       (m_xpoints, 0.0);
    Array <OneD, NekDouble > rho     (m_xpoints, 0.0);
    Array <OneD, NekDouble > mu      (m_xpoints, 0.0);
    Array <OneD, NekDouble > vv      (m_xpoints, 0.0);
    Array <OneD, NekDouble > velocity(m_xpoints, 0.0);
    Array <OneD, NekDouble > test    (m_xpoints, 0.0);


    NekDouble dd, dm, scale;
    NekDouble xcher, ycher;
    int index = -1;

    z[0]           = 0.0;
    NekDouble sumd = 0.0;

    for (int i = 1; i < m_xpoints ; i++)
    {
        z[i] = z[i-1] + 0.5 * (xx[i] - xx[i-1]) * (ff[3][i] + ff[3][i-1]);
        dm   = ff[3][i-1] - ff[1][i-1];
        dd   = ff[3][i] - ff[1][i];
        sumd = sumd + 0.5 * (xx[i] - xx[i-1]) * (dd + dm);
    }

    scale = sumd;

    ofstream file3;
    file3.open("physical_data.dat");

    NekDouble xin, rex, delsx, delta;

    for (int i = 0; i < m_xpoints; i++)
    {
        for (int k = 0; k < 5; k++)
        {
            v[k] = ff[k][i];
        }
        COMPBL(v, dv);
        u[i]        = ff[1][i];
        t[i]        = ff[3][i];
        rho[i]      = (1.0 / ff[3][i]);
        vv[i]       = -ff[0][i]/sqrt(m_uInf);
        mu[i]       = pow(t[i], 1.5) * (1 + m_Suth) / (t[i] + m_Suth) / (m_Re);
        velocity[i] = ff[0][i] ;
    }

    NekDouble scale2, coeff;

    for (int i = 0; i < nQuadraturePts; i++)
    {
        if (i%100000 == 0)
        {
            cout << "i" << "  " << i << "/" << nQuadraturePts << endl;
        }

        xcher  = x_QuadraturePts[i];
        ycher  = y_QuadraturePts[i];

        scale  = sumd;
        xin    = xcher;
        rex    = 0.5 * pow(((m_Re) / scale), 2) + (m_Re) * xin;
        delsx  = sqrt(2.0 / rex) * scale * (xin)* m_Pr;
        scale  = scale / delsx;
        delta  = 4.91 * sqrt((xin * m_mu) / (m_rhoInf * m_uInf));
        scale2 = ycher * (scale * delta) / sqrt(etamax) ;
        coeff  = 0.5 * sqrt( 2 / (xcher*m_Re)) ;

        if (scale2 > z[m_xpoints-3])
        {
            u_QuadraturePts[i]   = 1;
            rho_QuadraturePts[i] = 1;
            T_QuadraturePts[i]   = 1.0 / rho_QuadraturePts[i];
            v_QuadraturePts[i]   =  coeff * (z[m_xpoints-3] -
                                             velocity[m_xpoints-3]);

            file3 << xcher                 << "    "
                  << ycher                 << "    "
                  << velocity[m_xpoints-3] << "    "
                  << z[m_xpoints-3]        << "    "
                  << u[m_xpoints-3]
                  << endl;
        }
        else
        {
            for (int j = 0 ; j< m_xpoints-1; j++)
            {
                if ((z[j] <= scale2) && (z[j+1] > scale2))
                {
                    index = j;
                    break;
                }
            }
            if (index == -1)
            {
                ASSERTL0(false, "Could not determine index in CompressibleBL");
            }

            u_QuadraturePts[i]   = u[index];
            rho_QuadraturePts[i] = rho[index];
            T_QuadraturePts[i]   = 1.0/rho_QuadraturePts[i];
            v_QuadraturePts[i]   = coeff * (u[index]*scale2 - velocity[index]);
        }
    }
}


/**
 * Calculate the compressible boundary layer solution for a given 3D mesh
 * and dump the solution into a .rst file.
 */
int main(int argc, char *argv[])
{
    // Variable initialisation
    int nmax  = 5;
    int maxit = 10;

    string opt;

    int  i, j, k, numModes;

    Array<OneD, NekDouble>               xx(m_xpoints, 0.0);
    Array<OneD, Array<OneD, NekDouble> > ff(nmax);
    Array<OneD, NekDouble>               parameter(9, 0.0);

    for (i = 0; i < nmax; i++)
    {
        ff[i] = Array<OneD, NekDouble> (m_xpoints);
    }

    Array<OneD, NekDouble > vstart(nmax, 0.0);
    Array<OneD, NekDouble > v (2);
    Array<OneD, NekDouble > dv(2);
    Array<OneD, NekDouble > f (2);
    Array<OneD, NekDouble > f1(2);
    Array<OneD, NekDouble > f2(2);

    NekDouble al11, al21, al12, al22, det;

    // Reading the session file
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);

    // Read in mesh from input file and create an object of class MeshGraph
    SpatialDomains::MeshGraphSharedPtr graphShPt
        = SpatialDomains::MeshGraph::Read(vSession);

    int expdim = graphShPt->GetMeshDimension();

    int  nElements, nQuadraturePts = 0;
    Array<OneD, NekDouble> x_QuadraturePts;
    Array<OneD, NekDouble> y_QuadraturePts;
    Array<OneD, NekDouble> z_QuadraturePts;

    if (expdim == 2)
    {
        MultiRegions::ContField2DSharedPtr Domain;
        Domain = MemoryManager<MultiRegions::ContField2D>
            ::AllocateSharedPtr(vSession, graphShPt,
                                vSession->GetVariable(0));

        // Get the total number of elements
        nElements = Domain->GetExpSize();
        std::cout << "Number of elements                 = "
                  << nElements << std::endl;

        // Get the total number of quadrature points (depends on n. modes)
        nQuadraturePts = Domain->GetTotPoints();
        std::cout << "Number of quadrature points        = "
                  << nQuadraturePts << std::endl;

        // Coordinates of the quadrature points
        x_QuadraturePts = Array<OneD, NekDouble>(nQuadraturePts);
        y_QuadraturePts = Array<OneD, NekDouble>(nQuadraturePts);
        z_QuadraturePts = Array<OneD, NekDouble>(nQuadraturePts);
        Domain->GetCoords(x_QuadraturePts, y_QuadraturePts, z_QuadraturePts);
    }
    else if (expdim == 3)
    {
        MultiRegions::ContField3DSharedPtr Domain;
        Domain = MemoryManager<MultiRegions::ContField3D>
            ::AllocateSharedPtr(vSession, graphShPt, vSession->GetVariable(0));

        // Get the total number of elements
        nElements = Domain->GetExpSize();
        std::cout << "Number of elements                 = "
                  << nElements << std::endl;

        // Get the total number of quadrature points (depends on n. modes)
        nQuadraturePts = Domain->GetTotPoints();
        std::cout << "Number of quadrature points        = "
                  << nQuadraturePts << std::endl;

        // Coordinates of the quadrature points
        x_QuadraturePts = Array<OneD, NekDouble>(nQuadraturePts);
        y_QuadraturePts = Array<OneD, NekDouble>(nQuadraturePts);
        z_QuadraturePts = Array<OneD, NekDouble>(nQuadraturePts);
        Domain->GetCoords(x_QuadraturePts, y_QuadraturePts, z_QuadraturePts);
    }
    else
    {
        ASSERTL0(false, "Routine available for 2D and 3D problem only.")
    }

    // Loading parameters from session file
    vSession->LoadParameter("Re",          m_Re,     1.0);
    vSession->LoadParameter("Mach",        m_Mach,   1.0);
    vSession->LoadParameter("TInf",        m_Tinf,   1.0);
    vSession->LoadParameter("Twall",       m_Twall,  1.0);
    vSession->LoadParameter("Gamma",       m_Gamma,  1.0);
    vSession->LoadParameter("Pr",          m_Pr,     1.0);
    vSession->LoadParameter("L",           m_long,   1.0);
    vSession->LoadParameter("rhoInf",      m_rhoInf, 1.0);
    vSession->LoadParameter("uInf",        m_uInf,   1.0);
    vSession->LoadParameter("GasConstant", m_R,      1.0);
    vSession->LoadParameter("vInf",        m_vInf,   1.0);
    vSession->LoadParameter("mu",          m_mu,     1.0);

    // Rescaling factors
    m_Suth = 110.4 / m_Tinf;
    m_Tw   = m_Twall / m_Tinf;
    m_Re   = m_Re / m_long;

    cout << "Number of points" << "   " << m_xpoints << endl;

    // Defining the solution arrays
    Array<OneD, NekDouble> u_QuadraturePts  (nQuadraturePts, 0.0);
    Array<OneD, NekDouble> v_QuadraturePts  (nQuadraturePts, 0.0);
    Array<OneD, NekDouble> rho_QuadraturePts(nQuadraturePts, 0.0);
    Array<OneD, NekDouble> T_QuadraturePts  (nQuadraturePts, 0.0);

    // Calculation of the similarity variables
    if (m_Tw > 0)
    {
        vstart[3] = m_Tw;
    }
    if (m_Tw < 0.0)
    {
        v[1] = 1.0 + 0.5 * 0.84 * (m_Gamma - 1) * (m_Mach * m_Mach);
        v[0] = 0.47 * pow(v[1], 0.21);
    }
    else
    {
        v[1] = 0.062 * pow(m_Mach, 2) - 0.1 * (m_Tw - 1.0) *
        (10 + m_Mach) / (0.2 + m_Mach);
        v[0] = 0.45 - 0.01 * m_Mach + (m_Tw - 1.0) * 0.06;
        m_Twall = m_Tw;
    }

    dv[0] = v[0] * 0.01;

    if (m_Tw < 0.0)
    {
        dv[1] = v[1] * 0.01;
    }
    else
    {
        dv[1] = 0.1;
    }

    vstart[2] = v[0];

    if (m_Tw < 0)
    {
        vstart[3] = v[1];
        m_Twall = vstart[3];
    }
    else
    {
        vstart[4] = v[1];
    }

    RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);

    for (k = 0; k < maxit; k++)
    {
        vstart[2] = v[0];

        if (m_Tw < 0)
        {
            vstart[3] = v[1];
            m_Twall   = vstart[3];
        }
        else
        {
            vstart[4] = v[1];
        }

        RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);

        NekDouble err = fabs(ff[1][m_xpoints] - 1) +
        fabs(ff[3][m_xpoints] - 1);

        cout << "err" << scientific << setprecision(9) << "   " << err << endl;

        if (expdim == 2)
        {
            if (err < errtol)
            {
                cout << "Calculating" << endl;
                OUTPUT(m_xpoints, xx, ff, nQuadraturePts, x_QuadraturePts,
                       y_QuadraturePts, u_QuadraturePts, v_QuadraturePts,
                       rho_QuadraturePts, T_QuadraturePts);
                break;
            }
            else
            {
                f[0]      = ff[1][m_xpoints] - 1;
                f[1]      = ff[3][m_xpoints] - 1;
                vstart[2] = v[0] + dv[0];

                if (m_Tw < 0)
                {
                    vstart[3] = v[1];
                    m_Twall   = vstart[3];
                }
                else
                {
                    vstart[4] = v[1];
                }

                RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);

                f1[0] = ff[1][m_xpoints] - 1;
                f1[1] = ff[3][m_xpoints] - 1;

                vstart[2] = v[0];

                if (m_Tw < 0)
                {
                    vstart[3] = v[1] + dv[1];
                    m_Twall   = vstart[3];
                }
                else
                {
                    vstart[4] = v[1] + dv[1];
                }

                RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);

                f2[0] = ff[1][m_xpoints] - 1;
                f2[1] = ff[3][m_xpoints] - 1;

                al11 = (f1[0] - f[0]) / dv[0];
                al21 = (f1[1] - f[1]) / dv[0];
                al12 = (f2[0] - f[0]) / dv[1];
                al22 = (f2[1] - f[1]) / dv[1];
                det = al11 * al22 - al21 * al12;

                dv[0] = ( - al22 * f[0] + al12 * f[1]) / det;
                dv[1] = (al21 * f[0] - al11 * f[1]) / det;
                v[0]  = v[0] + dv[0];
                v[1]  = v[1] + dv[1];
            }
        }
        else if (expdim == 3)
        {
            {
                if (err < errtol)
                {
                    cout << "Calculating" << endl;
                    OUTPUT(m_xpoints, xx, ff, nQuadraturePts, x_QuadraturePts,
                           z_QuadraturePts, u_QuadraturePts, v_QuadraturePts,
                           rho_QuadraturePts, T_QuadraturePts);
                    break;
                }
                else
                {
                    f[0]      = ff[1][m_xpoints] - 1;
                    f[1]      = ff[3][m_xpoints] - 1;
                    vstart[2] = v[0] + dv[0];

                    if (m_Tw < 0)
                    {
                        vstart[3] = v[1];
                        m_Twall   = vstart[3];
                    }
                    else
                    {
                        vstart[4] = v[1];
                    }

                    RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);

                    f1[0] = ff[1][m_xpoints] - 1;
                    f1[1] = ff[3][m_xpoints] - 1;

                    vstart[2] = v[0];

                    if (m_Tw < 0)
                    {
                        vstart[3] = v[1] + dv[1];
                        m_Twall   = vstart[3];
                    }
                    else
                    {
                        vstart[4] = v[1] + dv[1];
                    }

                    RKDUMB(vstart, 5, 0.0, etamax, m_xpoints, xx, ff);

                    f2[0] = ff[1][m_xpoints] - 1;
                    f2[1] = ff[3][m_xpoints] - 1;

                    al11 = (f1[0] - f[0]) / dv[0];
                    al21 = (f1[1] - f[1]) / dv[0];
                    al12 = (f2[0] - f[0]) / dv[1];
                    al22 = (f2[1] - f[1]) / dv[1];
                    det = al11 * al22 - al21 * al12;

                    dv[0] = ( - al22 * f[0] + al12 * f[1]) / det;
                    dv[1] = (al21 * f[0] - al11 * f[1]) / det;
                    v[0]  = v[0] + dv[0];
                    v[1]  = v[1] + dv[1];
                }
            }
        }
    }

    // Verification of the compressible similarity solution
    ofstream verif;
    verif.open("similarity_solution.dat");
    for (i=0; i< nQuadraturePts; i++)
    {
        verif << scientific << setprecision(9) << x_QuadraturePts[i]
              << "  \t  " << y_QuadraturePts[i] << "  \t " ;
        verif << scientific << setprecision(9) << u_QuadraturePts[i]
              << "  \t  " << v_QuadraturePts[i] << "  \t " ;
        verif << scientific << setprecision(9) << rho_QuadraturePts[i]
              << "  \t  " << T_QuadraturePts[i] << endl;
    }
    verif.close();

    // Calculation of the physical variables
    for (i = 0; i < nQuadraturePts; i++)
    {
        rho_QuadraturePts[i] = rho_QuadraturePts[i] * m_rhoInf;
        u_QuadraturePts[i]   = u_QuadraturePts[i] * m_uInf;
        v_QuadraturePts[i]   = v_QuadraturePts[i] * m_uInf;
        T_QuadraturePts[i]   = T_QuadraturePts[i] * m_Tinf;

        T_QuadraturePts[i] = T_QuadraturePts[i] * rho_QuadraturePts[i] * m_R;
        T_QuadraturePts[i] = T_QuadraturePts[i] / (m_Gamma-1);
        T_QuadraturePts[i] = T_QuadraturePts[i] + 0.5 * rho_QuadraturePts[i] * (
            pow(u_QuadraturePts[i], 2.0) + pow(v_QuadraturePts[i], 2.0));

        u_QuadraturePts[i] = u_QuadraturePts[i] * rho_QuadraturePts[i];
        v_QuadraturePts[i] = v_QuadraturePts[i] * rho_QuadraturePts[i];
    }
    string file_name;
    if (expdim == 2)
    {
        MultiRegions::ContField2DSharedPtr Domain;
        Domain = MemoryManager<MultiRegions::ContField2D>
            ::AllocateSharedPtr(vSession, graphShPt, vSession->GetVariable(0));

        Array<OneD, MultiRegions::ExpListSharedPtr> Exp(4);
        MultiRegions::ExpList2DSharedPtr Exp2D_uk;
        Exp2D_uk = MemoryManager<MultiRegions::ExpList2D>
            ::AllocateSharedPtr(vSession, graphShPt);

        MultiRegions::ExpList2DSharedPtr Exp2D_vk;
        Exp2D_vk = MemoryManager<MultiRegions::ExpList2D>
            ::AllocateSharedPtr(vSession, graphShPt);

        MultiRegions::ExpList2DSharedPtr Exp2D_rhok;
        Exp2D_rhok = MemoryManager<MultiRegions::ExpList2D>
            ::AllocateSharedPtr(vSession, graphShPt);

        MultiRegions::ExpList2DSharedPtr Exp2D_Tk;
        Exp2D_Tk = MemoryManager<MultiRegions::ExpList2D>
            ::AllocateSharedPtr(vSession, graphShPt);

        // Filling the 2D expansion using a recursive algorithm based on the
        // mesh ordering
        LibUtilities::BasisSharedPtr Basis;
        Basis    = Domain->GetExp(0)->GetBasis(0);
        numModes = Basis->GetNumModes();

        std::cout << "Number of modes = " << numModes << std::endl;

        // Copying the ukGlobal vector in m_phys (with the same pattern of
        // m_phys)
        Vmath::Vcopy(nQuadraturePts, u_QuadraturePts, 1,
                     Exp2D_uk->UpdatePhys(), 1);
        Vmath::Vcopy(nQuadraturePts, v_QuadraturePts, 1,
                     Exp2D_vk->UpdatePhys(), 1);
        Vmath::Vcopy(nQuadraturePts, rho_QuadraturePts, 1,
                     Exp2D_rhok->UpdatePhys(), 1);
        Vmath::Vcopy(nQuadraturePts, T_QuadraturePts, 1,
                     Exp2D_Tk->UpdatePhys(), 1);

        // Initialisation of the ExpList Exp
        Exp[0] = Exp2D_rhok;
        Exp[1] = Exp2D_uk;
        Exp[2] = Exp2D_vk;
        Exp[3] = Exp2D_Tk;

        // Expansion coefficient extraction (necessary to write the .fld file)
        Exp[0]->FwdTrans(Exp2D_rhok->GetPhys(), Exp[0]->UpdateCoeffs());
        Exp[1]->FwdTrans(Exp2D_uk->GetPhys(),   Exp[1]->UpdateCoeffs());
        Exp[2]->FwdTrans(Exp2D_vk->GetPhys(),   Exp[2]->UpdateCoeffs());
        Exp[3]->FwdTrans(Exp2D_Tk->GetPhys(),   Exp[3]->UpdateCoeffs());

        // Definition of the name of the .fld file
        cout << argv[1] << endl;
        string tmp = argv[1];
        int len = tmp.size();
        for (i = 0; i < len-4; ++i)
        {
            file_name += argv[1][i];
        }
        file_name = file_name+".rst";

        // Definition of the Field
        std::vector<LibUtilities::FieldDefinitionsSharedPtr>
            FieldDef = Exp[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for (j = 0; j < 4; j++)
        {
            for (i = 0; i < FieldDef.size(); i++)
            {
                if (j == 0)
                {
                    FieldDef[i]->m_fields.push_back("rho");
                }
                else if (j == 1)
                {
                    FieldDef[i]->m_fields.push_back("rhou");
                }
                else if (j == 2 )
                {
                    FieldDef[i]->m_fields.push_back("rhov");
                }
                else if (j == 3 )
                {
                    FieldDef[i]->m_fields.push_back("E");
                }
                Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
            }
        }

        LibUtilities::Write(file_name, FieldDef, FieldData);
    }
    else if (expdim == 3)
    {
        MultiRegions::ContField3DSharedPtr Domain;
        Domain = MemoryManager<MultiRegions::ContField3D>
            ::AllocateSharedPtr(vSession, graphShPt, vSession->GetVariable(0));

        Array<OneD,NekDouble> w_QuadraturePts;
        w_QuadraturePts = Array<OneD,NekDouble>(nQuadraturePts, 0.0);
        Array<OneD, MultiRegions::ExpListSharedPtr> Exp(5);

        MultiRegions::ExpList3DSharedPtr Exp3D_uk;
        Exp3D_uk = MemoryManager<MultiRegions::ExpList3D>
            ::AllocateSharedPtr(vSession, graphShPt);

        MultiRegions::ExpList3DSharedPtr Exp3D_vk;
        Exp3D_vk = MemoryManager<MultiRegions::ExpList3D>
            ::AllocateSharedPtr(vSession, graphShPt);

        MultiRegions::ExpList3DSharedPtr Exp3D_wk;
        Exp3D_wk = MemoryManager<MultiRegions::ExpList3D>
            ::AllocateSharedPtr(vSession, graphShPt);

        MultiRegions::ExpList3DSharedPtr Exp3D_rhok;
        Exp3D_rhok = MemoryManager<MultiRegions::ExpList3D>
            ::AllocateSharedPtr(vSession, graphShPt);

        MultiRegions::ExpList3DSharedPtr Exp3D_Tk;
        Exp3D_Tk = MemoryManager<MultiRegions::ExpList3D>
            ::AllocateSharedPtr(vSession, graphShPt);

        // Filling the 3D expansion using a recursive algorithm based
        // on the mesh ordering
        LibUtilities::BasisSharedPtr Basis;
        Basis    = Domain->GetExp(0)->GetBasis(0);
        numModes = Basis->GetNumModes();

        std::cout<< "Number of modes = " << numModes << std::endl;

        // Copying the ukGlobal vector in m_phys (with the same pattern
        // of m_phys)
        Vmath::Vcopy(nQuadraturePts, rho_QuadraturePts, 1,
                     Exp3D_rhok->UpdatePhys(), 1);
        Vmath::Vcopy(nQuadraturePts, u_QuadraturePts, 1,
                     Exp3D_uk->UpdatePhys(), 1);
        Vmath::Vcopy(nQuadraturePts, w_QuadraturePts, 1,
                     Exp3D_vk->UpdatePhys(), 1);
        Vmath::Vcopy(nQuadraturePts, v_QuadraturePts, 1,
                     Exp3D_wk->UpdatePhys(), 1);
        Vmath::Vcopy(nQuadraturePts, T_QuadraturePts, 1,
                     Exp3D_Tk->UpdatePhys(), 1);

        // Initialisation of the ExpList Exp
        Exp[0] = Exp3D_rhok;
        Exp[1] = Exp3D_uk;
        Exp[2] = Exp3D_vk;
        Exp[3] = Exp3D_wk;
        Exp[4] = Exp3D_Tk;

        // Expansion coefficient extraction (necessary to write the .fld file)
        Exp[0]->FwdTrans(Exp3D_rhok->GetPhys(), Exp[0]->UpdateCoeffs());
        Exp[1]->FwdTrans(Exp3D_uk->GetPhys(), Exp[1]->UpdateCoeffs());
        Exp[2]->FwdTrans(Exp3D_vk->GetPhys(), Exp[2]->UpdateCoeffs());
        Exp[3]->FwdTrans(Exp3D_wk->GetPhys(), Exp[3]->UpdateCoeffs());
        Exp[4]->FwdTrans(Exp3D_Tk->GetPhys(), Exp[4]->UpdateCoeffs());

        // Definition of the name of the .fld file
        cout << argv[1] << endl;
        string tmp = argv[1];
        int len = tmp.size();
        for (i = 0; i < len-4; ++i)
        {
            file_name += argv[1][i];
        }
        file_name = file_name+".rst";

        // Definition of the Field
        std::vector<LibUtilities::FieldDefinitionsSharedPtr>
            FieldDef = Exp[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for (j = 0; j < 5; j++)
        {
            for (i = 0; i < FieldDef.size(); i++)
            {
                if (j == 0)
                {
                    FieldDef[i]->m_fields.push_back("rho");
                }
                else if (j == 1)
                {
                    FieldDef[i]->m_fields.push_back("rhou");
                }
                else if (j == 2 )
                {
                    FieldDef[i]->m_fields.push_back("rhov");
                }
                else if (j == 3 )
                {
                    FieldDef[i]->m_fields.push_back("rhow");
                }
                else if (j == 4 )
                {
                    FieldDef[i]->m_fields.push_back("E");
                }
                Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
            }
        }

        LibUtilities::Write(file_name, FieldDef, FieldData);
    }

    std::cout <<"----------------------------------------------------\n";
    std::cout <<"\n=================================================\n";
    std::cout <<"Similarity solution \n";
    std::cout <<"===================================================\n";
    std::cout <<"***************************************************\n";
    std::cout <<"DATA FROM THE SESSION FILE:\n";
    std::cout << "Reynolds number                  = " << m_Re
              << "\t[-]"   << std::endl;
    std::cout << "Mach number                      = " << m_Mach
              << "\t[-]"   << std::endl;
    std::cout << "Characteristic length            = " << m_long
              << "\t[m]" << std::endl;
    std::cout << "U_infinity                       = " << m_uInf
              << "\t[m/s]" << std::endl;
    std::cout <<"***************************************************\n";
    std::cout <<"---------------------------------------------------\n";
    std::cout <<"MESH and EXPANSION DATA:\n";
    std::cout << "Done." << std::endl;

    return 0;
}
