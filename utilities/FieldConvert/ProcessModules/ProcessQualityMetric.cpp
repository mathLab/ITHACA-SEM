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
        LibUtilities::Interp2D(pFrom[0], pFrom[1], q, pTo[0], pTo[1], out);
        
        //Vmath::Vcopy(q.num_elements(), &q[0], 1, &phys[offset], 1);
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

/*
 * @brief Generate mapping between a simplex and the reference element.
 *
 * Mapping that requires the coordinates of the three vertices (x,y) of the
 * triangular element and outputs the function coefficients that defines the
 * mapping to the reference element for each coordinate (xfunc and yfunc)
 *
 * @param x   co-ordinates of triangle/tetrahedron
 * @param xf  components of mapping
 */
inline DNekMat MappingIdealToRef(SpatialDomains::GeometrySharedPtr geom)
{
    int n = geom->GetNumVerts(), i, j;

    DNekMat map   (n, n, 1.0, eFULL);
    DNekMat mapref(n, n, 1.0, eFULL);

    // Extract coordinate information.
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n-1; ++j)
        {
            map(j,i) = (*geom->GetVertex(i))[j];
        }
    }

    // Set up reference triangle or tetrahedron mapping.
    if (n == 3)
    {
        mapref(0,0) = -1.0; mapref(1,0) = -1.0;
        mapref(0,1) =  1.0; mapref(1,1) = -1.0;
        mapref(0,2) = -1.0; mapref(1,2) =  1.0;
    }
    else if (n == 4)
    {
        mapref(0,0) = -1.0; mapref(1,0) = -1.0; mapref(2,0) = -1.0;
        mapref(0,1) =  1.0; mapref(1,1) = -1.0; mapref(2,1) = -1.0;
        mapref(0,2) = -1.0; mapref(1,2) =  1.0; mapref(2,2) = -1.0;
        mapref(0,3) = -1.0; mapref(1,3) = -1.0; mapref(2,3) =  1.0;
    }

    map.Invert();

    DNekMat newmap = mapref * map;
    DNekMat mapred(n-1, n-1, 1.0, eFULL);

    for (i = 0; i < n-1; ++i)
    {
        for (j = 0; j < n-1; ++j)
        {
            mapred(i,j) = newmap(i,j);
        }
    }

    return mapred;
}

Array<OneD, NekDouble> ProcessQualityMetric::GetQ(LocalRegions::ExpansionSharedPtr e)
{
    SpatialDomains::GeometrySharedPtr    geom = e->GetGeom();
    StdRegions::StdExpansionSharedPtr    chi  = e->GetGeom()->GetXmap();
    LibUtilities::PointsKeyVector        p    = chi->GetPointsKeys();
    SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

    DNekMat i2rm = MappingIdealToRef(geom);

    SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);

    const int pts = deriv[0][0].num_elements();
    const int nq  = chi->GetTotPoints();
    const int expDim = chi->GetNumBases();

    Array<OneD, NekDouble> jac(pts, 0.0);
    Array<OneD, NekDouble> eta(nq);

    switch (expDim)
    {
        case 2:
        {
            Vmath::Vvtvvtm(pts, &deriv[0][0][0], 1, &deriv[1][1][0], 1,
                           &deriv[1][0][0], 1, &deriv[0][1][0], 1,
                           &jac[0],         1);
            break;
        }
        case 3:
        {
            Array<OneD, NekDouble> tmp(pts, 0.0);
            Vmath::Vvtvvtm(pts, &deriv[1][1][0], 1, &deriv[2][2][0], 1,
                           &deriv[2][1][0], 1, &deriv[1][2][0], 1,
                           &tmp[0],         1);
            Vmath::Vvtvp  (pts, &deriv[0][0][0], 1, &tmp[0],         1,
                           &jac[0],         1, &jac[0],         1);
            Vmath::Vvtvvtm(pts, &deriv[2][1][0], 1, &deriv[0][2][0], 1,
                           &deriv[0][1][0], 1, &deriv[2][2][0], 1,
                           &tmp[0],         1);
            Vmath::Vvtvp  (pts, &deriv[1][0][0], 1, &tmp[0],         1,
                           &jac[0],         1, &jac[0],         1);
            Vmath::Vvtvvtm(pts, &deriv[0][1][0], 1, &deriv[1][2][0], 1,
                           &deriv[1][1][0], 1, &deriv[0][2][0], 1,
                           &tmp[0],         1);
            Vmath::Vvtvp  (pts, &deriv[2][0][0], 1, &tmp[0],         1,
                           &jac[0],         1, &jac[0],         1);
            break;
        }
    }

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

        jacIdeal = jac * i2rm;

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
        if(eta[k] > 1 || eta[k] < 0)
            cout << eta[k] << endl;
    }

    if (pts == 1)
    {
        Vmath::Fill(nq-1, eta[0], &eta[1], 1);
    }

    return eta;
}

}
}
