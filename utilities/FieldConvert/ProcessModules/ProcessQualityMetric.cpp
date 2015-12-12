////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJacobianEnergy.cpp
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
//  Description: Compute energy of Jacobian.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessQualityMetric.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/Foundations/Interp.h>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessQualityMetric::className =
        GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "qualitymetric"),
                ProcessQualityMetric::create,
                "add quality metric to field.");

ProcessQualityMetric::ProcessQualityMetric(FieldSharedPtr f) :
    ProcessModule(f)
{

}

ProcessQualityMetric::~ProcessQualityMetric()
{
}

void ProcessQualityMetric::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessQualityMetric: Process Jacobian fld" << endl;
    }

    Array<OneD, NekDouble> &phys   = m_f->m_exp[0]->UpdatePhys();
    Array<OneD, NekDouble> &coeffs = m_f->m_exp[0]->UpdateCoeffs();
    int expdim = m_f->m_graph->GetMeshDimension();

    for(int i =0; i < m_f->m_exp[0]->GetExpSize(); ++i)
    {
        // copy Jacobian into field
        LocalRegions::ExpansionSharedPtr Elmt = m_f->m_exp[0]->GetExp(i);
        int offset = m_f->m_exp[0]->GetPhys_Offset(i);
        Array<OneD, NekDouble> q = GetQ(Elmt);
        Array<OneD, NekDouble> out = phys + offset;

        // Project onto output stuff
        LibUtilities::PointsKeyVector pFrom = Elmt->GetGeom()->GetPointsKeys();
        LibUtilities::PointsKeyVector pTo   = Elmt->GetPointsKeys();

        if(expdim == 2)
            LibUtilities::Interp2D(pFrom[0], pFrom[1], q, pTo[0], pTo[1], out);
        else if(expdim == 3)
            LibUtilities::Interp3D(pFrom[0], pFrom[1], pFrom[2], q,
                                   pTo[0], pTo[1], pTo[2], out);
        else
            ASSERTL0(false,"mesh dim makes no sense");
    }

    m_f->m_exp[0]->FwdTrans_IterPerExp(phys, coeffs);

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (int i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("QualityMetric");
        m_f->m_exp[0]->AppendFieldData(FieldDef[i], FieldData[i]);
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

inline vector<DNekMat> MappingIdealToRef(SpatialDomains::GeometrySharedPtr geom,
                                 StdRegions::StdExpansionSharedPtr chi)
{
    int dim = geom->GetShapeDim();
    vector<DNekMat> ret;

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

        for(int i = 0; i < b[0]->GetNumPoints(); i++)
        {
            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                DNekMat dxdz(2,2,1.0,eFULL);
                dxdz(0,0) = -xy[0][0]/4.0 + xy[0][0]/4.0*v[j]
                            +xy[1][0]/4.0 - xy[1][0]/4.0*v[j]
                            +xy[2][0]/4.0 + xy[2][0]/4.0*v[j]
                            -xy[3][0]/4.0 - xy[3][0]/4.0*v[j];

                dxdz(0,1) = +xy[0][0]/4.0*u[i] - xy[0][0]/4.0
                            -xy[1][0]/4.0*u[i] - xy[1][0]/4.0
                            +xy[2][0]/4.0*u[i] + xy[2][0]/4.0
                            -xy[3][0]/4.0*u[i] + xy[3][0]/4.0;

                dxdz(1,0) = -xy[0][1]/4.0 + xy[0][1]/4.0*v[j]
                            +xy[1][1]/4.0 - xy[1][1]/4.0*v[j]
                            +xy[2][1]/4.0 + xy[2][1]/4.0*v[j]
                            -xy[3][1]/4.0 - xy[3][1]/4.0*v[j];

                dxdz(1,1) = +xy[0][1]/4.0*u[i] - xy[0][1]/4.0
                            -xy[1][1]/4.0*u[i] - xy[1][1]/4.0
                            +xy[2][1]/4.0*u[i] + xy[2][1]/4.0
                            -xy[3][1]/4.0*u[i] + xy[3][1]/4.0;

                dxdz.Invert();
                ret.push_back(dxdz);
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
                dxdz(0,0) = -xy[0][0]/4.0 + xy[0][0]/4.0*v[j]
                            +xy[1][0]/4.0 - xy[1][0]/4.0*v[j]
                            -xy[2][0]/4.0 - xy[2][0]/4.0*v[j];

                dxdz(0,1) = +xy[0][0]/4.0*u[i] - xy[0][0]/4.0
                            -xy[1][0]/4.0*u[i] - xy[1][0]/4.0
                            -xy[2][0]/4.0*u[i] + xy[2][0]/4.0;

                dxdz(1,0) = -xy[0][1]/4.0 + xy[0][1]/4.0*v[j]
                            +xy[1][1]/4.0 - xy[1][1]/4.0*v[j]
                            -xy[2][1]/4.0 - xy[2][1]/4.0*v[j];

                dxdz(1,1) = +xy[0][1]/4.0*u[i] - xy[0][1]/4.0
                            -xy[1][1]/4.0*u[i] - xy[1][1]/4.0
                            -xy[2][1]/4.0*u[i] + xy[2][1]/4.0;

                dxdz.Invert();
                ret.push_back(dxdz);
            }
        }
    }
    else
    {
        ASSERTL0(false,"not coded");
    }



    return ret;
}

Array<OneD, NekDouble> ProcessQualityMetric::GetQ(LocalRegions::ExpansionSharedPtr e)
{
    SpatialDomains::GeometrySharedPtr    geom = e->GetGeom();
    StdRegions::StdExpansionSharedPtr    chi  = e->GetGeom()->GetXmap();
    LibUtilities::PointsKeyVector        p    = chi->GetPointsKeys();
    SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

    vector<DNekMat> i2rm = MappingIdealToRef(geom,chi);

    SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);

    const int pts = deriv[0][0].num_elements();
    const int nq  = chi->GetTotPoints();
    const int expDim = chi->GetNumBases();

    Array<OneD, NekDouble> eta(nq);

    for (int k = 0; k < pts; ++k)
    {
        DNekMat jac     (expDim, expDim, 0.0, eFULL);
        DNekMat jacIdeal(expDim, expDim, 0.0, eFULL);

        for (int i = 0; i < expDim; ++i)
        {
            for (int j = 0; j < expDim; ++j)
            {
                jac(j,i) = deriv[i][j][k];
            }
        }

        jacIdeal = jac * i2rm[k];

        NekDouble jacDet = jacIdeal(0,0) * jacIdeal(1,1) - jacIdeal(0,1)*jacIdeal(1,0);
        NekDouble frob = 0.0;

        for (int i = 0; i < expDim; ++i)
        {
            for (int j = 0; j < expDim; ++j)
            {
                frob += jacIdeal(i,j) * jacIdeal(i,j);
            }
        }

        NekDouble delta = 0.0; // fixme
        NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4*delta*delta));

        eta[k] = expDim * pow(sigma,2.0/expDim) / frob;
    }

    if (pts == 1)
    {
        Vmath::Fill(nq-1, eta[0], &eta[1], 1);
    }

    return eta;
}

}
}
