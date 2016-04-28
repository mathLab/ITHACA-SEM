////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.cpp
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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessVarOpti.h"

#include <boost/thread/mutex.hpp>

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx;

inline NekDouble GetElFunctional(ElDataSharedPtr d, optimiser opti, int dim)
{
    SpatialDomains::GeometrySharedPtr    geom = d->el->GetGeom(dim);
    StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();
    LibUtilities::PointsKeyVector        p    = chi->GetPointsKeys();
    mtx.lock();
    SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

    mtx.unlock();
    SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);

    const int pts = deriv[0][0].num_elements();

    ASSERTL0(pts == d->maps.size(), "what");

    Array<OneD, NekDouble> dW(pts);

    switch (opti)
    {
        case eLinEl:
        {
            for(int k = 0; k < pts; k++)
            {
                Array<OneD, NekDouble> jac(9,0.0);
                for (int i = 0; i < dim; ++i)
                {
                    for (int j = 0; j < dim; ++j)
                    {
                        jac[j+i*3] = deriv[i][j][k];
                    }
                }

                Array<OneD, NekDouble> jacIdeal(9,0.0);
                jacIdeal[0] = jac[0]*d->maps[k][0] + jac[3]*d->maps[k][1] + jac[6]*d->maps[k][2];
                jacIdeal[1] = jac[1]*d->maps[k][0] + jac[4]*d->maps[k][1] + jac[7]*d->maps[k][2];
                jacIdeal[2] = jac[2]*d->maps[k][0] + jac[5]*d->maps[k][1] + jac[8]*d->maps[k][2];

                jacIdeal[3] = jac[0]*d->maps[k][3] + jac[3]*d->maps[k][4] + jac[6]*d->maps[k][5];
                jacIdeal[4] = jac[1]*d->maps[k][3] + jac[4]*d->maps[k][4] + jac[7]*d->maps[k][5];
                jacIdeal[5] = jac[2]*d->maps[k][3] + jac[5]*d->maps[k][4] + jac[8]*d->maps[k][5];

                jacIdeal[6] = jac[0]*d->maps[k][6] + jac[3]*d->maps[k][7] + jac[6]*d->maps[k][8];
                jacIdeal[7] = jac[1]*d->maps[k][6] + jac[4]*d->maps[k][7] + jac[7]*d->maps[k][8];
                jacIdeal[8] = jac[2]*d->maps[k][6] + jac[5]*d->maps[k][7] + jac[8]*d->maps[k][8];

                NekDouble a = (jacIdeal[0]*jacIdeal[0]+jacIdeal[1]*jacIdeal[1]+jacIdeal[2]*jacIdeal[2]-1.0)*
                              (jacIdeal[0]*jacIdeal[0]+jacIdeal[1]*jacIdeal[1]+jacIdeal[2]*jacIdeal[2]-1.0);
                NekDouble b = (jacIdeal[3]*jacIdeal[3]+jacIdeal[4]*jacIdeal[4]+jacIdeal[5]*jacIdeal[5]-1.0)*
                              (jacIdeal[3]*jacIdeal[3]+jacIdeal[4]*jacIdeal[4]+jacIdeal[5]*jacIdeal[5]-1.0);
                NekDouble c = (jacIdeal[6]*jacIdeal[6]+jacIdeal[7]*jacIdeal[7]+jacIdeal[8]*jacIdeal[8]-1.0)*
                              (jacIdeal[6]*jacIdeal[6]+jacIdeal[7]*jacIdeal[7]+jacIdeal[8]*jacIdeal[8]-1.0);
                NekDouble d = (jacIdeal[0]*jacIdeal[6]+jacIdeal[1]*jacIdeal[7]+jacIdeal[2]*jacIdeal[8])*
                              (jacIdeal[0]*jacIdeal[6]+jacIdeal[1]*jacIdeal[7]+jacIdeal[2]*jacIdeal[8]);
                NekDouble e = (jacIdeal[3]*jacIdeal[6]+jacIdeal[4]*jacIdeal[7]+jacIdeal[5]*jacIdeal[8])*
                              (jacIdeal[3]*jacIdeal[6]+jacIdeal[4]*jacIdeal[7]+jacIdeal[5]*jacIdeal[8]);
                NekDouble f = (jacIdeal[0]*jacIdeal[3]+jacIdeal[1]*jacIdeal[4]+jacIdeal[3]*jacIdeal[5])*
                              (jacIdeal[0]*jacIdeal[3]+jacIdeal[1]*jacIdeal[4]+jacIdeal[3]*jacIdeal[5]);

                NekDouble trEtE = 0.25 * (a+b+c) + 0.5 * (d+e+f);
                NekDouble ljacDet;
                if(dim == 2)
                {
                    ljacDet = log(jacIdeal[0] * jacIdeal[4] - jacIdeal[3]*jacIdeal[1]);
                }
                else
                {
                    ljacDet = log(jacIdeal[0]*(jacIdeal[4]*jacIdeal[8]-jacIdeal[5]*jacIdeal[7])
                                 -jacIdeal[3]*(jacIdeal[1]*jacIdeal[8]-jacIdeal[2]*jacIdeal[7])
                                 +jacIdeal[6]*(jacIdeal[1]*jacIdeal[5]-jacIdeal[2]*jacIdeal[4]));
                }

                NekDouble nu = 0.45;
                NekDouble mu = 1.0/2.0/(1.0+nu);
                NekDouble K = 1.0 / 3.0 / (1.0 - 2.0 * nu);
                dW[k] = K *0.5 * ljacDet * ljacDet + mu * trEtE;

            }
            break;
        }
        case eHypEl:
        {
            for(int k = 0; k < pts; k++)
            {
                Array<OneD, NekDouble> jac(9,0.0);
                for (int i = 0; i < dim; ++i)
                {
                    for (int j = 0; j < dim; ++j)
                    {
                        jac[j+i*3] = deriv[i][j][k];
                    }
                }

                Array<OneD, NekDouble> jacIdeal(9,0.0);
                jacIdeal[0] = jac[0]*d->maps[k][0] + jac[3]*d->maps[k][1] + jac[6]*d->maps[k][2];
                jacIdeal[1] = jac[1]*d->maps[k][0] + jac[4]*d->maps[k][1] + jac[7]*d->maps[k][2];
                jacIdeal[2] = jac[2]*d->maps[k][0] + jac[5]*d->maps[k][1] + jac[8]*d->maps[k][2];

                jacIdeal[3] = jac[0]*d->maps[k][3] + jac[3]*d->maps[k][4] + jac[6]*d->maps[k][5];
                jacIdeal[4] = jac[1]*d->maps[k][3] + jac[4]*d->maps[k][4] + jac[7]*d->maps[k][5];
                jacIdeal[5] = jac[2]*d->maps[k][3] + jac[5]*d->maps[k][4] + jac[8]*d->maps[k][5];

                jacIdeal[6] = jac[0]*d->maps[k][6] + jac[3]*d->maps[k][7] + jac[6]*d->maps[k][8];
                jacIdeal[7] = jac[1]*d->maps[k][6] + jac[4]*d->maps[k][7] + jac[7]*d->maps[k][8];
                jacIdeal[8] = jac[2]*d->maps[k][6] + jac[5]*d->maps[k][7] + jac[8]*d->maps[k][8];

                NekDouble jacDet;
                if(dim == 2)
                {
                    jacDet = jacIdeal[0] * jacIdeal[4] - jacIdeal[3]*jacIdeal[1];
                }
                else
                {
                    jacDet = jacIdeal[0]*(jacIdeal[4]*jacIdeal[8]-jacIdeal[5]*jacIdeal[7])
                            -jacIdeal[3]*(jacIdeal[1]*jacIdeal[8]-jacIdeal[2]*jacIdeal[7])
                            +jacIdeal[6]*(jacIdeal[1]*jacIdeal[5]-jacIdeal[2]*jacIdeal[4]);
                }

                NekDouble ljacDet = log(jacDet);
                NekDouble nu = 0.45;
                NekDouble mu = 1.0/2.0/(1.0+nu);
                NekDouble lam = nu / (1.0 - 2.0 * nu) / (1.0+nu);
                dW[k] = 0.5 * mu * (jacDet * jacDet - 3.0) - mu * ljacDet + 0.5 * lam * ljacDet * ljacDet;

            }
            break;
        }
        case eRoca:
        {
            for(int k = 0; k < pts; k++)
            {
                Array<OneD, NekDouble> jac(9,0.0);
                for (int i = 0; i < dim; ++i)
                {
                    for (int j = 0; j < dim; ++j)
                    {
                        jac[j+i*3] = deriv[i][j][k];
                    }
                }

                Array<OneD, NekDouble> jacIdeal(9,0.0);
                jacIdeal[0] = jac[0]*d->maps[k][0] + jac[3]*d->maps[k][1] + jac[6]*d->maps[k][2];
                jacIdeal[1] = jac[1]*d->maps[k][0] + jac[4]*d->maps[k][1] + jac[7]*d->maps[k][2];
                jacIdeal[2] = jac[2]*d->maps[k][0] + jac[5]*d->maps[k][1] + jac[8]*d->maps[k][2];

                jacIdeal[3] = jac[0]*d->maps[k][3] + jac[3]*d->maps[k][4] + jac[6]*d->maps[k][5];
                jacIdeal[4] = jac[1]*d->maps[k][3] + jac[4]*d->maps[k][4] + jac[7]*d->maps[k][5];
                jacIdeal[5] = jac[2]*d->maps[k][3] + jac[5]*d->maps[k][4] + jac[8]*d->maps[k][5];

                jacIdeal[6] = jac[0]*d->maps[k][6] + jac[3]*d->maps[k][7] + jac[6]*d->maps[k][8];
                jacIdeal[7] = jac[1]*d->maps[k][6] + jac[4]*d->maps[k][7] + jac[7]*d->maps[k][8];
                jacIdeal[8] = jac[2]*d->maps[k][6] + jac[5]*d->maps[k][7] + jac[8]*d->maps[k][8];

                NekDouble jacDet, i2rmDet;
                if(dim == 2)
                {
                    jacDet = jacIdeal[0] * jacIdeal[4] - jacIdeal[3]*jacIdeal[1];
                    i2rmDet = 1.0 / (d->maps[k][0]*d->maps[k][4] - d->maps[k][1]*d->maps[k][3]);
                }
                else
                {
                    jacDet = jacIdeal[0]*(jacIdeal[4]*jacIdeal[8]-jacIdeal[5]*jacIdeal[7])
                            -jacIdeal[3]*(jacIdeal[1]*jacIdeal[8]-jacIdeal[2]*jacIdeal[7])
                            +jacIdeal[6]*(jacIdeal[1]*jacIdeal[5]-jacIdeal[2]*jacIdeal[4]);
                    ASSERTL0(false,"need to extend to do i2rmDet in 3d");
                }

                NekDouble frob = 0.0;

                frob += jacIdeal[0] * jacIdeal[0];
                frob += jacIdeal[1] * jacIdeal[1];
                frob += jacIdeal[2] * jacIdeal[2];
                frob += jacIdeal[3] * jacIdeal[3];
                frob += jacIdeal[4] * jacIdeal[4];
                frob += jacIdeal[5] * jacIdeal[5];
                frob += jacIdeal[6] * jacIdeal[6];
                frob += jacIdeal[7] * jacIdeal[7];
                frob += jacIdeal[8] * jacIdeal[8];

                NekDouble de = fabs(i2rmDet) * sqrt(1E-6 + 1E-3);

                NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                dW[k] = frob / dim / pow(sigma, 2.0/dim);
            }
            break;
        }
        case eWins:
        {
            for(int k = 0; k < pts; k++)
            {
                Array<OneD, NekDouble> jac(9,0.0);
                for (int i = 0; i < dim; ++i)
                {
                    for (int j = 0; j < dim; ++j)
                    {
                        jac[j+i*3] = deriv[i][j][k];
                    }
                }

                Array<OneD, NekDouble> jacIdeal(9,0.0);
                jacIdeal[0] = jac[0]*d->maps[k][0] + jac[3]*d->maps[k][1] + jac[6]*d->maps[k][2];
                jacIdeal[1] = jac[1]*d->maps[k][0] + jac[4]*d->maps[k][1] + jac[7]*d->maps[k][2];
                jacIdeal[2] = jac[2]*d->maps[k][0] + jac[5]*d->maps[k][1] + jac[8]*d->maps[k][2];

                jacIdeal[3] = jac[0]*d->maps[k][3] + jac[3]*d->maps[k][4] + jac[6]*d->maps[k][5];
                jacIdeal[4] = jac[1]*d->maps[k][3] + jac[4]*d->maps[k][4] + jac[7]*d->maps[k][5];
                jacIdeal[5] = jac[2]*d->maps[k][3] + jac[5]*d->maps[k][4] + jac[8]*d->maps[k][5];

                jacIdeal[6] = jac[0]*d->maps[k][6] + jac[3]*d->maps[k][7] + jac[6]*d->maps[k][8];
                jacIdeal[7] = jac[1]*d->maps[k][6] + jac[4]*d->maps[k][7] + jac[7]*d->maps[k][8];
                jacIdeal[8] = jac[2]*d->maps[k][6] + jac[5]*d->maps[k][7] + jac[8]*d->maps[k][8];

                NekDouble jacDet;
                if(dim == 2)
                {
                    jacDet = jacIdeal[0] * jacIdeal[4] - jacIdeal[3]*jacIdeal[1];
                }
                else
                {
                    jacDet = jacIdeal[0]*(jacIdeal[4]*jacIdeal[8]-jacIdeal[5]*jacIdeal[7])
                            -jacIdeal[3]*(jacIdeal[1]*jacIdeal[8]-jacIdeal[2]*jacIdeal[7])
                            +jacIdeal[6]*(jacIdeal[1]*jacIdeal[5]-jacIdeal[2]*jacIdeal[4]);
                    ASSERTL0(false,"need to extend to do i2rmDet in 3d");
                }

                if(jacDet < 1E-10)
                {
                    jacDet = 1E-10;
                }

                NekDouble frob = 0.0;

                frob += jacIdeal[0] * jacIdeal[0];
                frob += jacIdeal[1] * jacIdeal[1];
                frob += jacIdeal[2] * jacIdeal[2];
                frob += jacIdeal[3] * jacIdeal[3];
                frob += jacIdeal[4] * jacIdeal[4];
                frob += jacIdeal[5] * jacIdeal[5];
                frob += jacIdeal[6] * jacIdeal[6];
                frob += jacIdeal[7] * jacIdeal[7];
                frob += jacIdeal[8] * jacIdeal[8];

                dW[k] = frob / jacDet;
            }
            break;
        }
    }

    mtx.lock();
    NekDouble inte = chi->Integral(dW);
    mtx.unlock();
    return inte;
}

inline vector<Array<OneD, NekDouble> > MappingIdealToRef(SpatialDomains::GeometrySharedPtr geom,
                                 StdRegions::StdExpansionSharedPtr chi)
{
    int dim = geom->GetShapeDim();
    vector<Array<OneD, NekDouble> > ret;

    if(geom->GetShapeType() == LibUtilities::eQuadrilateral)
    {
        vector<Array<OneD, NekDouble> > xy;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(2);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xy.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u = b[0]->GetZ();
        Array<OneD, NekDouble> v = b[1]->GetZ();

        for(int j = 0; j < b[1]->GetNumPoints(); j++)
        {
            for(int i = 0; i < b[0]->GetNumPoints(); i++)
            {
                NekDouble a1 = 0.5*(1.0-u[i]), a2 = 0.5*(1.0+u[i]);
                NekDouble b1 = 0.5*(1.0-v[j]), b2 = 0.5*(1.0+v[j]);
                DNekMat dxdz(2,2,1.0,eFULL);

                dxdz(0,0) = 0.5*(-b1*xy[0][0] + b1*xy[1][0] + b2*xy[2][0] - b2*xy[3][0]);
                dxdz(1,0) = 0.5*(-b1*xy[0][1] + b1*xy[1][1] + b2*xy[2][1] - b2*xy[3][1]);

                dxdz(0,1) = 0.5*(-a1*xy[0][0] - a2*xy[1][0] + a2*xy[2][0] + a1*xy[3][0]);
                dxdz(1,1) = 0.5*(-a1*xy[0][1] - a2*xy[1][1] + a2*xy[2][1] + a1*xy[3][1]);

                NekDouble det = 1.0/(dxdz(0,0)*dxdz(1,1) - dxdz(1,0)*dxdz(0,1));

                dxdz.Invert();
                Array<OneD, NekDouble> r(9,0.0);
                r[0] = dxdz(0,0);
                r[1] = dxdz(1,0);
                r[3] = dxdz(0,1);
                r[4] = dxdz(1,1);
                ret.push_back(r);
            }
        }
    }
    else if(geom->GetShapeType() == LibUtilities::eTriangle)
    {
        vector<Array<OneD, NekDouble> > xy;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(2);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xy.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u = b[0]->GetZ();
        Array<OneD, NekDouble> v = b[1]->GetZ();

        for(int i = 0; i < b[0]->GetNumPoints(); i++)
        {
            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                DNekMat dxdz(2,2,1.0,eFULL);
                dxdz(0,0) = -xy[0][0]/2.0 + xy[1][0]/2.0;

                dxdz(0,1) = -xy[0][0]/2.0 + xy[2][0]/2.0;

                dxdz(1,0) = -xy[0][1]/2.0 + xy[1][1]/2.0;

                dxdz(1,1) = -xy[0][1]/2.0 + xy[2][1]/2.0;

                dxdz.Invert();
                Array<OneD, NekDouble> r(9,0.0);
                r[0] = dxdz(0,0);
                r[1] = dxdz(1,0);
                r[3] = dxdz(0,1);
                r[4] = dxdz(1,1);
                ret.push_back(r);
            }
        }
    }
    else if(geom->GetShapeType() == LibUtilities::eTetrahedron)
    {
        vector<Array<OneD, NekDouble> > xyz;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u = b[0]->GetZ();
        Array<OneD, NekDouble> v = b[1]->GetZ();
        Array<OneD, NekDouble> z = b[2]->GetZ();

        for(int i = 0; i < b[0]->GetNumPoints(); i++)
        {
            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for(int k = 0; k < b[2]->GetNumPoints(); k++)
                {
                    DNekMat dxdz(3,3,1.0,eFULL);
                    dxdz(0,0) = -xyz[0][0]/2.0 + xyz[1][0]/2.0;

                    dxdz(0,1) = -xyz[0][0]/2.0 + xyz[2][0]/2.0;

                    dxdz(0,2) = -xyz[0][0]/2.0 + xyz[3][0]/2.0;


                    dxdz(1,0) = -xyz[0][1]/2.0 + xyz[1][1]/2.0;

                    dxdz(1,1) = -xyz[0][1]/2.0 + xyz[2][1]/2.0;

                    dxdz(1,2) = -xyz[0][1]/2.0 + xyz[3][1]/2.0;


                    dxdz(2,0) = -xyz[0][2]/2.0 + xyz[1][2]/2.0;

                    dxdz(2,1) = -xyz[0][2]/2.0 + xyz[2][2]/2.0;

                    dxdz(2,2) = -xyz[0][2]/2.0 + xyz[3][2]/2.0;

                    dxdz.Invert();
                    Array<OneD, NekDouble> r(9,0.0);
                    r[0] = dxdz(0,0);
                    r[1] = dxdz(1,0);
                    r[3] = dxdz(0,1);
                    r[4] = dxdz(1,1);
                    r[2] = dxdz(2,0);
                    r[5] = dxdz(2,1);
                    r[6] = dxdz(0,2);
                    r[7] = dxdz(1,2);
                    r[8] = dxdz(2,2);
                    ret.push_back(r);
                }
            }
        }
    }
    else if(geom->GetShapeType() == LibUtilities::ePrism)
    {
        vector<Array<OneD, NekDouble> > xyz;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> eta1 = b[0]->GetZ();
        Array<OneD, NekDouble> eta2 = b[1]->GetZ();
        Array<OneD, NekDouble> eta3 = b[2]->GetZ();

        for(int k = 0; k < b[2]->GetNumPoints(); k++)
        {
            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for(int i = 0; i < b[0]->GetNumPoints(); i++)
                {
                    NekDouble xi1 = 0.5*(1+eta1[i])*(1-eta3[k])-1.0;
                    NekDouble a1 = 0.5*(1-xi1),     a2 = 0.5*(1+xi1);
                    NekDouble b1 = 0.5*(1-eta2[j]), b2 = 0.5*(1+eta2[j]);
                    NekDouble c1 = 0.5*(1-eta3[k]), c2 = 0.5*(1+eta3[k]);

                    DNekMat dxdz(3,3,1.0,eFULL);

                    dxdz(0,0) = 0.5*(-b1*xyz[0][0] + b1*xyz[1][0] + b2*xyz[2][0] - b2*xyz[3][0]);
                    dxdz(1,0) = 0.5*(-b1*xyz[0][1] + b1*xyz[1][1] + b2*xyz[2][1] - b2*xyz[3][1]);
                    dxdz(2,0) = 0.5*(-b1*xyz[0][2] + b1*xyz[1][2] + b2*xyz[2][2] - b2*xyz[3][2]);

                    dxdz(0,1) = 0.5*((a2-c1)*xyz[0][0] - a2*xyz[1][0] + a2*xyz[2][0] + (c1-a2)*xyz[3][0] - c2*xyz[4][0] + c2*xyz[5][0]);
                    dxdz(1,1) = 0.5*((a2-c1)*xyz[0][1] - a2*xyz[1][1] + a2*xyz[2][1] + (c1-a2)*xyz[3][1] - c2*xyz[4][1] + c2*xyz[5][1]);
                    dxdz(2,1) = 0.5*((a2-c1)*xyz[0][2] - a2*xyz[1][2] + a2*xyz[2][2] + (c1-a2)*xyz[3][2] - c2*xyz[4][2] + c2*xyz[5][2]);

                    dxdz(0,2) = 0.5*(-b1*xyz[0][0] - b2*xyz[3][0] + b1*xyz[4][0] + b2*xyz[5][0]);
                    dxdz(1,2) = 0.5*(-b1*xyz[0][1] - b2*xyz[3][1] + b1*xyz[4][1] + b2*xyz[5][1]);
                    dxdz(2,2) = 0.5*(-b1*xyz[0][2] - b2*xyz[3][2] + b1*xyz[4][2] + b2*xyz[5][2]);

                    dxdz.Invert();
                    Array<OneD, NekDouble> r(9,0.0);
                    r[0] = dxdz(0,0);
                    r[1] = dxdz(1,0);
                    r[3] = dxdz(0,1);
                    r[4] = dxdz(1,1);
                    r[2] = dxdz(2,0);
                    r[5] = dxdz(2,1);
                    r[6] = dxdz(0,2);
                    r[7] = dxdz(1,2);
                    r[8] = dxdz(2,2);
                    ret.push_back(r);
                }
            }
        }
    }
    else
    {
        ASSERTL0(false,"not coded");
    }

    return ret;
}

ModuleKey ProcessVarOpti::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "varopti"),
    ProcessVarOpti::create,
    "Optimise mesh locations.");

ProcessVarOpti::ProcessVarOpti(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["linearelastic"] =
        ConfigOption(true, "", "Optimise for linear elasticity");
    m_config["winslow"] =
        ConfigOption(true, "", "Optimise for winslow");
    m_config["roca"] =
        ConfigOption(true, "", "Optimise for roca method");
    m_config["hyperelastic"] =
        ConfigOption(true, "", "Optimise for hyper elasticity");
    m_config["numthreads"] =
        ConfigOption(false, "1", "Number of threads");
    m_config["nq"] =
        ConfigOption(false, "0", "Number of quad points");
}

ProcessVarOpti::~ProcessVarOpti()
{
}

void ProcessVarOpti::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessVarOpti: Optimising... " << endl;
    }

    if(m_config["linearelastic"].beenSet)
    {
        opti = eLinEl;
    }
    else if(m_config["winslow"].beenSet)
    {
        opti = eWins;
    }
    else if(m_config["roca"].beenSet)
    {
        opti = eRoca;
    }
    else if(m_config["hyperelastic"].beenSet)
    {
        opti = eHypEl;
    }
    else
    {
        ASSERTL0(false,"not opti type set");
    }

    m_mesh->m_nummode = m_config["nq"].as<int>();

    ASSERTL0(m_mesh->m_nummode > 2,"not specified high-order");

    if(m_mesh->m_expDim == 2 && m_mesh->m_spaceDim == 3)
    {
        ASSERTL0(false,"cannot deal with manifolds");
    }

    FillQuadPoints();

    GetElementMap();

    vector<vector<NodeSharedPtr> > freenodes = GetColouredNodes();

    ResidualSharedPtr res = boost::shared_ptr<Residual>(new Residual);
    res->val = 1.0;

    vector<vector<NodeOpti> > optiNodes;
    for(int i = 0; i < freenodes.size(); i++)
    {
        vector<NodeOpti> ns;
        for(int j = 0; j < freenodes[i].size(); j++)
        {
            NodeElMap::iterator it = nodeElMap.find(freenodes[i][j]->m_id);
            ASSERTL0(it != nodeElMap.end(),"could not find");
            ns.push_back(NodeOpti(freenodes[i][j],it->second,opti,res,m_mesh->m_spaceDim));
        }
        optiNodes.push_back(ns);
    }

    NekDouble functionalStart = 0.0;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        functionalStart += GetElFunctional(dataSet[i], opti, m_mesh->m_spaceDim);
    }

    cout << scientific << endl;
    cout << "starting energy: " << functionalStart << endl;

    int nThreads = m_config["numthreads"].as<int>();

    int ctr = 0;
    Thread::ThreadMaster tms;
    tms.SetThreadingType("ThreadManagerBoost");
    Thread::ThreadManagerSharedPtr tm = tms.CreateInstance(Thread::ThreadMaster::SessionJob, nThreads);

    while (res->val > 1e-3)
    {
        ctr++;
        res->val = 0.0;
        for(int i = 0; i < optiNodes.size(); i++)
        {
            vector<Thread::ThreadJob*> jobs;
            for(int j = 0; j < optiNodes[i].size(); j++)
            {
                jobs.push_back(optiNodes[i][j].GetJob());
            }
            //cout << " -- inner loop " << i+1 << "/" << optiNodes.size()
            //     << " of size: " << jobs.size() << endl;

            tm->SetNumWorkers(0);
            tm->QueueJobs(jobs);
            tm->SetNumWorkers(nThreads);
            tm->Wait();

            /*for(int j = 0; j < optiNodes[i].size(); j++)
            {
                optiNodes[i][j].Run();
            }*/
        }

        cout << ctr <<  "\tResidual: " << res->val << endl;

    }

    functionalStart = 0.0;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        functionalStart += GetElFunctional(dataSet[i], opti, m_mesh->m_spaceDim);
    }
    cout << "end energy: " << functionalStart << endl;
}

void ProcessVarOpti::NodeOpti::Optimise()
{
    //it doesnt matter at this point what the dim is so long that
    //in the 2d case z is left as zero
    Array<OneD, NekDouble> G = GetGrad();

    if(sqrt(G[0]*G[0] + G[1]*G[1] + G[2]*G[2]) > 1e-6)
    {
        //needs to optimise
        NekDouble currentW = GetFunctional();
        NekDouble xc = node->m_x;
        NekDouble yc = node->m_y;
        NekDouble zc = node->m_z;
        NekDouble alpha = 1.0;
        bool found = false;
        while(alpha > 1e-6)
        {
            node->m_x = xc - alpha * G[0];
            node->m_y = yc - alpha * G[1];
            node->m_z = zc - alpha * G[2];

            if(GetFunctional() < currentW)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }
        if(!found)
        {
            //reset the node
            node->m_x = xc;
            node->m_y = yc;
            node->m_z = zc;
            cout << "warning: had to reset node" << endl;
        }
        mtx.lock();
        res->val += sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                         (node->m_z-zc)*(node->m_z-zc));
        mtx.unlock();
    }
}

NekDouble dir[6][3] = {{1.0,0.0,0.0},{-1.0,0.0,0.0},
                       {0.0,1.0,0.0},{0.0,-1.0,0.0},
                       {0.0,0.0,1.0},{0.0,0.0,-1.0}};

Array<OneD, NekDouble> ProcessVarOpti::NodeOpti::GetGrad()
{
    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    NekDouble zc = node->m_z;
    NekDouble dx = 1e-6;

    vector<NekDouble> w;

    for(int i = 0; i < dim * 2; i++)
    {
        node->m_x += dir[i][0] * dx;
        node->m_y += dir[i][1] * dx;
        node->m_z += dir[i][2] * dx;
        w.push_back(GetFunctional());
        node->m_x = xc;
        node->m_y = yc;
        node->m_z = zc;
    }

    Array<OneD, NekDouble> ret(3,0.0);

    ret[0] = (w[0] - w[1]) / 2.0 / dx;
    ret[1] = (w[2] - w[3]) / 2.0 / dx;
    if(dim == 3)
    {
        ret[3] = (w[4] - w[5]) / 2.0 / dx;
    }

    return ret;
}

NekDouble ProcessVarOpti::NodeOpti::GetFunctional()
{
    NekDouble r = 0.0;
    for(int i = 0; i < data.size(); i++)
    {
        r += GetElFunctional(data[i], opti,dim);
    }
    return r;
}

vector<vector<NodeSharedPtr> > ProcessVarOpti::GetColouredNodes()
{
    //this figures out the dirclet nodes and colors the others into paralell sets
    NodeSet boundaryNodes;

    vector<NodeSharedPtr> remain;

    switch (m_mesh->m_spaceDim)
    {
        case 2:
        {
            EdgeSet::iterator it;
            for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
            {
                if((*it)->m_elLink.size() == 2)
                {
                    continue;
                }

                boundaryNodes.insert((*it)->m_n1);
                boundaryNodes.insert((*it)->m_n2);
                for(int i = 0; i < (*it)->m_edgeNodes.size(); i++)
                {
                    boundaryNodes.insert((*it)->m_edgeNodes[i]);
                }
            }

            for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
            {
                vector<NodeSharedPtr> n = (*it)->m_edgeNodes;
                n.push_back((*it)->m_n1);
                n.push_back((*it)->m_n2);
                for(int j = 0; j < n.size(); j++)
                {
                    NodeSet::iterator nit = boundaryNodes.find(n[j]);
                    if(nit == boundaryNodes.end())
                    {
                        remain.push_back(n[j]);
                    }
                }
            }
            break;
        }
        case 3:
        {
            FaceSet::iterator it;
            for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end(); it++)
            {
                if((*it)->m_elLink.size() == 2)
                {
                    continue;
                }

                vector<NodeSharedPtr> vs = (*it)->m_vertexList;
                for(int j = 0; j < vs.size(); j++)
                {
                    boundaryNodes.insert(vs[j]);
                }

                vector<EdgeSharedPtr> es = (*it)->m_edgeList;
                for(int j = 0; j < es.size(); j++)
                {
                    for(int k = 0; k < es[j]->m_edgeNodes.size(); k++)
                    {
                        boundaryNodes.insert(es[j]->m_edgeNodes[k]);
                    }
                }
            }

            for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end(); it++)
            {
                vector<NodeSharedPtr> n = (*it)->m_vertexList;

                vector<EdgeSharedPtr> es = (*it)->m_edgeList;
                for(int j = 0; j < es.size(); j++)
                {
                    for(int k = 0; k < es[j]->m_edgeNodes.size(); k++)
                    {
                        n.push_back(es[j]->m_edgeNodes[k]);
                    }
                }
                for(int j = 0; j < n.size(); j++)
                {
                    NodeSet::iterator nit = boundaryNodes.find(n[j]);
                    if(nit == boundaryNodes.end())
                    {
                        remain.push_back(n[j]);
                    }
                }
            }
            break;
        }
        default:
            ASSERTL0(false,"space dim issue");
    }

    vector<vector<NodeSharedPtr> > ret;

    while (remain.size() > 0)
    {
        vector<NodeSharedPtr> layer;
        set<int> locked;
        set<int> completed;
        for(int i = 0; i < remain.size(); i++)
        {
            NodeElMap::iterator it = nodeElMap.find(remain[i]->m_id);
            ASSERTL0(it != nodeElMap.end(),"could not find");
            bool islocked = false;
            for(int j = 0; j < it->second.size(); j++)
            {
                set<int>::iterator sit = locked.find(it->second[j]->el->GetId());
                if(sit != locked.end())
                {
                    islocked = true;
                    break;
                }
            }
            if(!islocked)
            {
                layer.push_back(remain[i]);
                completed.insert(remain[i]->m_id);
                for(int j = 0; j < it->second.size(); j++)
                {
                    locked.insert(it->second[j]->el->GetId());
                }
            }
        }

        vector<NodeSharedPtr> tmp = remain;
        remain.clear();
        for(int i = 0; i < tmp.size(); i++)
        {
            set<int>::iterator sit = completed.find(tmp[i]->m_id);
            if(sit == completed.end())
            {
                remain.push_back(tmp[i]);
            }
        }
        ret.push_back(layer);
    }
    return ret;
}

void ProcessVarOpti::GetElementMap()
{
    //build ideal maps and structs;
    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        ElDataSharedPtr d = boost::shared_ptr<ElData>(new ElData);
        d->el = el;

        SpatialDomains::GeometrySharedPtr    geom = el->GetGeom(m_mesh->m_spaceDim);
        StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();

        vector<LibUtilities::BasisKey> basisKeys;

        for (int i = 0; i < m_mesh->m_expDim; ++i)
        {
            basisKeys.push_back(chi->GetBasis(i)->GetBasisKey());
        }

        StdRegions::StdExpansionSharedPtr chiMod;
        switch(chi->DetShapeType())
        {
            case LibUtilities::eTriangle:
                chiMod = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(
                    basisKeys[0], basisKeys[1]);
                break;
            case LibUtilities::eQuadrilateral:
                chiMod = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(
                    basisKeys[0], basisKeys[1]);
                break;
            case LibUtilities::eTetrahedron:
                chiMod = MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(
                    basisKeys[0], basisKeys[1], basisKeys[2]);
                break;
            case LibUtilities::ePrism:
                chiMod = MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(
                    basisKeys[0], basisKeys[1], basisKeys[2]);
                break;
            default:
                ASSERTL0(false, "nope");
        }

        d->maps = MappingIdealToRef(geom, chiMod);

        dataSet.push_back(d);
    }

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

        vector<NodeSharedPtr> vs = el->GetVertexList();
        for(int j = 0; j < vs.size(); j++)
        {
            nodeElMap[vs[j]->m_id].push_back(dataSet[i]);
        }

        vector<EdgeSharedPtr> es = el->GetEdgeList();
        for(int j = 0; j < es.size(); j++)
        {
            for(int k = 0; k < es[j]->m_edgeNodes.size(); k++)
            {
                nodeElMap[es[j]->m_edgeNodes[k]->m_id].push_back(dataSet[i]);
            }
        }

        if(m_mesh->m_spaceDim == 3)
        {
            vector<FaceSharedPtr> fs = el->GetFaceList();
            for(int j = 0; j < fs.size(); j++)
            {
                for(int k = 0; k < fs[j]->m_faceNodes.size(); k++)
                {
                    nodeElMap[fs[j]->m_faceNodes[k]->m_id].push_back(dataSet[i]);
                }
            }
        }
    }
}

void ProcessVarOpti::FillQuadPoints()
{

    //this routine is going to take linear or curved egdes and then project
    //a set of GLL points onto the system (either way) which is of order of the
    //input variable nq
    int nq = m_mesh->m_nummode;
    LibUtilities::PointsKey ekey(nq,LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;
    LibUtilities::PointsManager()[ekey]->GetPoints(gll);
    EdgeSet::iterator eit;
    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        SpatialDomains::Geometry1DSharedPtr geom = (*eit)->GetGeom(m_mesh->m_spaceDim);
        geom->FillGeom();
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

        vector<NodeSharedPtr> hons;
        (*eit)->m_edgeNodes.clear();

        switch (m_mesh->m_spaceDim)
        {
            case 2:
            {
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

                Array<OneD, NekDouble> xc(nq+1);
                Array<OneD, NekDouble> yc(nq+1);

                xmap->BwdTrans(coeffs0,xc);
                xmap->BwdTrans(coeffs1,yc);

                for(int i = 1; i < nq - 1; i++)
                {
                    Array<OneD, NekDouble> xp(1);
                    xp[0] = gll[i];

                    hons.push_back(boost::shared_ptr<Node>(new Node(
                            0,xmap->PhysEvaluate(xp,xc),xmap->PhysEvaluate(xp,yc),0.0)));
                }
                break;
            }
            case 3:
            {
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
                Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

                Array<OneD, NekDouble> xc(nq);
                Array<OneD, NekDouble> yc(nq);
                Array<OneD, NekDouble> zc(nq);

                xmap->BwdTrans(coeffs0,xc);
                xmap->BwdTrans(coeffs1,yc);
                xmap->BwdTrans(coeffs2,zc);

                for(int i = 1; i < nq - 1; i++)
                {
                    Array<OneD, NekDouble> xp(1);
                    xp[0] = gll[i];

                    hons.push_back(boost::shared_ptr<Node>(new Node(
                            0,xmap->PhysEvaluate(xp,xc),xmap->PhysEvaluate(xp,yc),
                              xmap->PhysEvaluate(xp,zc))));
                }
                break;
            }
        }

        (*eit)->m_edgeNodes = hons;
        (*eit)->m_curveType = LibUtilities::eGaussLobattoLegendre;
    }

    if(m_mesh->m_spaceDim == 3)
    {
        //need to do the same for faces
        int np = nq * (nq + 1) / 2;
        int ni = (nq-2)*(nq-3) / 2;
        FaceSet::iterator fit;

        LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                     LibUtilities::eNodalTriElec);
        Array<OneD, NekDouble> u, v;
        LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

        vector<NodeSharedPtr> hons;
        (*fit)->m_faceNodes.clear();

        for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
        {
            SpatialDomains::Geometry2DSharedPtr geom = (*fit)->GetGeom(3);
            geom->FillGeom();
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();
            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
            Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

            Array<OneD, NekDouble> xc(nq*nq);
            Array<OneD, NekDouble> yc(nq*nq);
            Array<OneD, NekDouble> zc(nq*nq);

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);
            xmap->BwdTrans(coeffs2,zc);

            for(int j = np-ni; j < np; j++)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = u[j];
                xp[1] = v[j];

                Array<OneD, NekDouble> xyz(3);
                xyz[0] = xmap->PhysEvaluate(xp, xc);
                xyz[1] = xmap->PhysEvaluate(xp, yc);
                xyz[2] = xmap->PhysEvaluate(xp, zc);

                hons.push_back(boost::shared_ptr<Node>(new Node(0,
                    xyz[0],xyz[1],xyz[2])));
            }

            (*fit)->m_faceNodes = hons;
            (*fit)->m_curveType = LibUtilities::eNodalTriElec;
        }
    }

    int id = m_mesh->m_vertexSet.size();
    //enumerate the curved nodes for mapping purposes
    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        for(int j = 0; j < (*eit)->m_edgeNodes.size(); j++)
        {
            (*eit)->m_edgeNodes[j]->m_id = id++;
        }
    }
    if(m_mesh->m_spaceDim == 3)
    {
        FaceSet::iterator fit;
        for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
        {
            for(int j = 0; j < (*fit)->m_faceNodes.size(); j++)
            {
                (*fit)->m_faceNodes[j]->m_id = id++;
            }
        }
    }

}

}
}
