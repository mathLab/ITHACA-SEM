////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include "ElUtil.h"
#include "ProcessVarOpti.h"

#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx2;

ElUtil::ElUtil(ElementSharedPtr e, DerivUtilSharedPtr d, PtsHelperSharedPtr p,
               ResidualSharedPtr r, int n)
{
    m_el = e;
    derivUtil = d;
    ptsHelp = p;
    res = r;
    m_mode = n;
    m_dim = m_el->GetDim();
    vector<NodeSharedPtr> ns;
    m_el->GetCurvedNodes(ns);
    nodes.resize(ns.size());
    for (int i = 0; i < ns.size(); ++i)
    {
        nodes[i].resize(m_dim);
        nodes[i][0] = &ns[i]->m_x;

        if (m_dim >= 2)
        {
            nodes[i][1] = &ns[i]->m_y;
        }

        if (m_dim >= 3)
        {
            nodes[i][2] = &ns[i]->m_z;
        }
    }
}

vector<Array<OneD, NekDouble> > ElUtil::MappingIdealToRef()
{
    //need to make ideal element out of old element
    ElmtConfig ec = m_el->GetConf();
    ec.m_order  = 1;
    ec.m_faceNodes = false;
    ec.m_volumeNodes = false;
    ec.m_reorient = false;

    ElementSharedPtr E = GetElementFactory().CreateInstance(
                            ec.m_e, ec, m_el->GetVertexList(),
                            m_el->GetTagList());

    SpatialDomains::GeometrySharedPtr    geom = E->GetGeom(m_dim);
    geom->FillGeom();
    StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();

    vector<Array<OneD, NekDouble> > ret;

    if(geom->GetShapeType() == LibUtilities::eQuadrilateral)
    {
        ASSERTL0(false,"Not coded");
        /*vector<Array<OneD, NekDouble> > xy;
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
        }*/
    }
    else if(geom->GetShapeType() == LibUtilities::eTriangle)
    {
        LibUtilities::PointsKey pkey(m_mode,
                                     LibUtilities::eNodalTriElec);
        Array<OneD, NekDouble> u, v;
        LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

        Array<OneD, NekDouble> xc(chi->GetTotPoints());
        Array<OneD, NekDouble> yc(chi->GetTotPoints());

        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

        chi->BwdTrans(coeffs0,xc);
        chi->BwdTrans(coeffs1,yc);

        NekVector<NekDouble> X(ptsHelp->ptsLow),Y(ptsHelp->ptsLow);
        for(int j = 0; j < u.num_elements(); j++)
        {
            Array<OneD, NekDouble> xp(2);
            xp[0] = u[j];
            xp[1] = v[j];

            X(j) = chi->PhysEvaluate(xp, xc);
            Y(j) = chi->PhysEvaluate(xp, yc);
        }

        NekVector<NekDouble> x1i(ptsHelp->ptsHigh),y1i(ptsHelp->ptsHigh),
                             x2i(ptsHelp->ptsHigh),y2i(ptsHelp->ptsHigh);

        x1i = derivUtil->VdmD[0]*X;
        y1i = derivUtil->VdmD[0]*Y;
        x2i = derivUtil->VdmD[1]*X;
        y2i = derivUtil->VdmD[1]*Y;

        for(int i = 0 ; i < ptsHelp->ptsHigh; i++)
        {
            DNekMat dxdz(2,2,1.0,eFULL);
            dxdz(0,0) = x1i(i);
            dxdz(0,1) = x2i(i);
            dxdz(1,0) = y1i(i);
            dxdz(1,1) = y2i(i);

            Array<OneD, NekDouble> r(10,0.0);
            r[9] = dxdz(0,0)*dxdz(1,1)-dxdz(1,0)*dxdz(0,1);

            dxdz.Invert();

            r[0] = dxdz(0,0);
            r[1] = dxdz(1,0);
            r[3] = dxdz(0,1);
            r[4] = dxdz(1,1);
            ret.push_back(r);
        }
    }
    else if(geom->GetShapeType() == LibUtilities::eTetrahedron)
    {
        LibUtilities::PointsKey pkey(m_mode,
                                     LibUtilities::eNodalTetElec);
        Array<OneD, NekDouble> u, v, w;
        LibUtilities::PointsManager()[pkey]->GetPoints(u, v, w);

        Array<OneD, NekDouble> xc(chi->GetTotPoints());
        Array<OneD, NekDouble> yc(chi->GetTotPoints());
        Array<OneD, NekDouble> zc(chi->GetTotPoints());

        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

        chi->BwdTrans(coeffs0,xc);
        chi->BwdTrans(coeffs1,yc);
        chi->BwdTrans(coeffs2,zc);

        NekVector<NekDouble> X(ptsHelp->ptsLow),Y(ptsHelp->ptsLow),Z(ptsHelp->ptsLow);
        for(int j = 0; j < u.num_elements(); j++)
        {
            Array<OneD, NekDouble> xp(3);
            xp[0] = u[j];
            xp[1] = v[j];
            xp[2] = w[j];

            X(j) = chi->PhysEvaluate(xp, xc);
            Y(j) = chi->PhysEvaluate(xp, yc);
            Z(j) = chi->PhysEvaluate(xp, zc);
        }

        NekVector<NekDouble> x1i(ptsHelp->ptsHigh),y1i(ptsHelp->ptsHigh),z1i(ptsHelp->ptsHigh),
                             x2i(ptsHelp->ptsHigh),y2i(ptsHelp->ptsHigh),z2i(ptsHelp->ptsHigh),
                             x3i(ptsHelp->ptsHigh),y3i(ptsHelp->ptsHigh),z3i(ptsHelp->ptsHigh);

        x1i = derivUtil->VdmD[0]*X;
        y1i = derivUtil->VdmD[0]*Y;
        z1i = derivUtil->VdmD[0]*Z;
        x2i = derivUtil->VdmD[1]*X;
        y2i = derivUtil->VdmD[1]*Y;
        z2i = derivUtil->VdmD[1]*Z;
        x3i = derivUtil->VdmD[2]*X;
        y3i = derivUtil->VdmD[2]*Y;
        z3i = derivUtil->VdmD[2]*Z;

        for(int i = 0 ; i < ptsHelp->ptsHigh; i++)
        {
            DNekMat dxdz(3,3,1.0,eFULL);
            dxdz(0,0) = x1i(i);
            dxdz(0,1) = x2i(i);
            dxdz(0,2) = x3i(i);
            dxdz(1,0) = y1i(i);
            dxdz(1,1) = y2i(i);
            dxdz(1,2) = y3i(i);
            dxdz(2,0) = z1i(i);
            dxdz(2,1) = z2i(i);
            dxdz(2,2) = z3i(i);

            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                  -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                  +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));

            dxdz.Invert();

            r[0] = dxdz(0,0);
            r[1] = dxdz(1,0);
            r[2] = dxdz(2,0);
            r[3] = dxdz(0,1);
            r[4] = dxdz(1,1);
            r[5] = dxdz(2,1);
            r[6] = dxdz(0,2);
            r[7] = dxdz(1,2);
            r[8] = dxdz(2,2);
            ret.push_back(r);
        }
    }
    else if(geom->GetShapeType() == LibUtilities::ePrism)
    {
        ASSERTL0(false, "not coded");
        /*vector<Array<OneD, NekDouble> > xyz;
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
        }*/
    }
    else
    {
        ASSERTL0(false,"not coded");
    }

    return ret;
}

void ElUtil::Evaluate()
{
    NekDouble mx = -1.0 * numeric_limits<double>::max();
    NekDouble mn =  numeric_limits<double>::max();

    if(m_dim == 2)
    {
        NekVector<NekDouble> X(nodes.size()),Y(nodes.size());
        for(int j = 0; j < nodes.size(); j++)
        {
            X(j) = *nodes[j][0];
            Y(j) = *nodes[j][1];
        }

        NekVector<NekDouble> x1i(nodes.size()),y1i(nodes.size()),
                             x2i(nodes.size()),y2i(nodes.size());

        x1i = derivUtil->VdmDL[0]*X;
        y1i = derivUtil->VdmDL[0]*Y;
        x2i = derivUtil->VdmDL[1]*X;
        y2i = derivUtil->VdmDL[1]*Y;

        for(int j = 0; j < nodes.size(); j++)
        {
            NekDouble jacDet = x1i(j) * y2i(j) - x2i(j)*y1i(j);
            mx = max(mx,jacDet);
            mn = min(mn,jacDet);
        }
    }
    else if(m_dim == 3)
    {
        NekVector<NekDouble> X(nodes.size()),Y(nodes.size()),Z(nodes.size());
        for(int j = 0; j < nodes.size(); j++)
        {
            X(j) = *nodes[j][0];
            Y(j) = *nodes[j][1];
            Z(j) = *nodes[j][2];
        }

        NekVector<NekDouble> x1i(nodes.size()),y1i(nodes.size()),z1i(nodes.size()),
                             x2i(nodes.size()),y2i(nodes.size()),z2i(nodes.size()),
                             x3i(nodes.size()),y3i(nodes.size()),z3i(nodes.size());

        x1i = derivUtil->VdmDL[0]*X;
        y1i = derivUtil->VdmDL[0]*Y;
        z1i = derivUtil->VdmDL[0]*Z;
        x2i = derivUtil->VdmDL[1]*X;
        y2i = derivUtil->VdmDL[1]*Y;
        z2i = derivUtil->VdmDL[1]*Z;
        x3i = derivUtil->VdmDL[2]*X;
        y3i = derivUtil->VdmDL[2]*Y;
        z3i = derivUtil->VdmDL[2]*Z;

        for(int j = 0; j < nodes.size(); j++)
        {
            DNekMat dxdz(3,3,1.0,eFULL);
            dxdz(0,0) = x1i(j);
            dxdz(0,1) = x2i(j);
            dxdz(0,2) = x3i(j);
            dxdz(1,0) = y1i(j);
            dxdz(1,1) = y2i(j);
            dxdz(1,2) = y3i(j);
            dxdz(2,0) = z1i(j);
            dxdz(2,1) = z2i(j);
            dxdz(2,2) = z3i(j);

            NekDouble jacDet = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                   -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                   +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));

            mx = max(mx,jacDet);
            mn = min(mn,jacDet);
        }
    }

    mtx2.lock();
    if(mn < 0)
    {
        res->startInv++;
    }
    res->worstJac = min(res->worstJac,mn/mx);
    mtx2.unlock();

    mtx2.lock();
    maps = MappingIdealToRef();
    mtx2.unlock();
}

ElUtilJob* ElUtil::GetJob()
{
    return new ElUtilJob(this);
}

}
}
