////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessQualityMetric.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description: Compute quality metric of Roca et al.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Interp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdTriExp.h>

#include "ProcessQualityMetric.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessQualityMetric::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "qualitymetric"),
        ProcessQualityMetric::create,
        "add quality metric to field.");

ProcessQualityMetric::ProcessQualityMetric(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["scaled"] =
        ConfigOption(true, "0", "use scaled jacobian instead");
}

ProcessQualityMetric::~ProcessQualityMetric()
{
}

void ProcessQualityMetric::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int nfields           = m_f->m_variables.size();
    m_f->m_variables.push_back("qualitymetric");
    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int NumHomogeneousDir = m_f->m_numHomogeneousDir;
    MultiRegions::ExpListSharedPtr exp;

    if (nfields)
    {
        m_f->m_exp.resize(nfields + 1);
        exp = m_f->AppendExpList(NumHomogeneousDir);

        m_f->m_exp[nfields] = exp;
    }
    else
    {
        exp = m_f->m_exp[0];
    }

    Array<OneD, NekDouble> &phys   = exp->UpdatePhys();
    Array<OneD, NekDouble> &coeffs = exp->UpdateCoeffs();

    for (int i = 0; i < exp->GetExpSize(); ++i)
    {
        // copy Jacobian into field
        LocalRegions::ExpansionSharedPtr Elmt = exp->GetExp(i);
        int offset = exp->GetPhys_Offset(i);
        Array<OneD, NekDouble> q = GetQ(Elmt,m_config["scaled"].as<bool>());
        Array<OneD, NekDouble> out = phys + offset;

        ASSERTL0(q.size() == Elmt->GetTotPoints(),
                 "number of points mismatch");
        Vmath::Vcopy(q.size(), q, 1, out, 1);
    }

    exp->FwdTrans_IterPerExp(phys, coeffs);
}

inline vector<DNekMat> MappingIdealToRef(SpatialDomains::GeometrySharedPtr geom,
                                         StdRegions::StdExpansionSharedPtr chi)
{
    vector<DNekMat> ret;

    if (geom->GetShapeType() == LibUtilities::eQuadrilateral)
    {
        vector<Array<OneD, NekDouble> > xy;
        for (int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(2);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xy.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u                          = b[0]->GetZ();
        Array<OneD, NekDouble> v                          = b[1]->GetZ();

        for (int j = 0; j < b[1]->GetNumPoints(); j++)
        {
            for (int i = 0; i < b[0]->GetNumPoints(); i++)
            {
                NekDouble a1 = 0.5 * (1.0 - u[i]), a2 = 0.5 * (1.0 + u[i]);
                NekDouble b1 = 0.5 * (1.0 - v[j]), b2 = 0.5 * (1.0 + v[j]);
                DNekMat dxdz(2, 2, 1.0, eFULL);

                dxdz(0, 0) = 0.5 * (-b1 * xy[0][0] + b1 * xy[1][0] +
                                    b2 * xy[2][0] - b2 * xy[3][0]);
                dxdz(1, 0) = 0.5 * (-b1 * xy[0][1] + b1 * xy[1][1] +
                                    b2 * xy[2][1] - b2 * xy[3][1]);

                dxdz(0, 1) = 0.5 * (-a1 * xy[0][0] - a2 * xy[1][0] +
                                    a2 * xy[2][0] + a1 * xy[3][0]);
                dxdz(1, 1) = 0.5 * (-a1 * xy[0][1] - a2 * xy[1][1] +
                                    a2 * xy[2][1] + a1 * xy[3][1]);

                dxdz.Invert();
                ret.push_back(dxdz);
            }
        }
    }
    else if (geom->GetShapeType() == LibUtilities::eTriangle)
    {
        vector<Array<OneD, NekDouble> > xy;
        for (int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(2);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xy.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u                          = b[0]->GetZ();
        Array<OneD, NekDouble> v                          = b[1]->GetZ();

        for (int i = 0; i < b[0]->GetNumPoints(); i++)
        {
            for (int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                DNekMat dxdz(2, 2, 1.0, eFULL);
                dxdz(0, 0) = -xy[0][0] / 2.0 + xy[1][0] / 2.0;

                dxdz(0, 1) = -xy[0][0] / 2.0 + xy[2][0] / 2.0;

                dxdz(1, 0) = -xy[0][1] / 2.0 + xy[1][1] / 2.0;

                dxdz(1, 1) = -xy[0][1] / 2.0 + xy[2][1] / 2.0;

                dxdz.Invert();
                ret.push_back(dxdz);
            }
        }
    }
    else if (geom->GetShapeType() == LibUtilities::eTetrahedron)
    {
        vector<Array<OneD, NekDouble> > xyz;
        for (int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u                          = b[0]->GetZ();
        Array<OneD, NekDouble> v                          = b[1]->GetZ();
        Array<OneD, NekDouble> z                          = b[2]->GetZ();

        for (int i = 0; i < b[0]->GetNumPoints(); i++)
        {
            for (int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for (int k = 0; k < b[2]->GetNumPoints(); k++)
                {
                    DNekMat dxdz(3, 3, 1.0, eFULL);
                    dxdz(0, 0) = -xyz[0][0] / 2.0 + xyz[1][0] / 2.0;

                    dxdz(0, 1) = -xyz[0][0] / 2.0 + xyz[2][0] / 2.0;

                    dxdz(0, 2) = -xyz[0][0] / 2.0 + xyz[3][0] / 2.0;

                    dxdz(1, 0) = -xyz[0][1] / 2.0 + xyz[1][1] / 2.0;

                    dxdz(1, 1) = -xyz[0][1] / 2.0 + xyz[2][1] / 2.0;

                    dxdz(1, 2) = -xyz[0][1] / 2.0 + xyz[3][1] / 2.0;

                    dxdz(2, 0) = -xyz[0][2] / 2.0 + xyz[1][2] / 2.0;

                    dxdz(2, 1) = -xyz[0][2] / 2.0 + xyz[2][2] / 2.0;

                    dxdz(2, 2) = -xyz[0][2] / 2.0 + xyz[3][2] / 2.0;

                    dxdz.Invert();
                    ret.push_back(dxdz);
                }
            }
        }
    }
    else if (geom->GetShapeType() == LibUtilities::ePrism)
    {
        vector<Array<OneD, NekDouble> > xyz;
        for (int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> eta1                       = b[0]->GetZ();
        Array<OneD, NekDouble> eta2                       = b[1]->GetZ();
        Array<OneD, NekDouble> eta3                       = b[2]->GetZ();

        for (int k = 0; k < b[2]->GetNumPoints(); k++)
        {
            for (int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for (int i = 0; i < b[0]->GetNumPoints(); i++)
                {
                    NekDouble xi1 = 0.5 * (1 + eta1[i]) * (1 - eta3[k]) - 1.0;
                    NekDouble a2  = 0.5 * (1 + xi1);
                    NekDouble b1  = 0.5 * (1 - eta2[j]),
                              b2  = 0.5 * (1 + eta2[j]);
                    NekDouble c1  = 0.5 * (1 - eta3[k]),
                              c2  = 0.5 * (1 + eta3[k]);

                    DNekMat dxdz(3, 3, 1.0, eFULL);

                    dxdz(0, 0) = 0.5 * (-b1 * xyz[0][0] + b1 * xyz[1][0] +
                                        b2 * xyz[2][0] - b2 * xyz[3][0]);
                    dxdz(1, 0) = 0.5 * (-b1 * xyz[0][1] + b1 * xyz[1][1] +
                                        b2 * xyz[2][1] - b2 * xyz[3][1]);
                    dxdz(2, 0) = 0.5 * (-b1 * xyz[0][2] + b1 * xyz[1][2] +
                                        b2 * xyz[2][2] - b2 * xyz[3][2]);

                    dxdz(0, 1) = 0.5 * ((a2 - c1) * xyz[0][0] - a2 * xyz[1][0] +
                                        a2 * xyz[2][0] + (c1 - a2) * xyz[3][0] -
                                        c2 * xyz[4][0] + c2 * xyz[5][0]);
                    dxdz(1, 1) = 0.5 * ((a2 - c1) * xyz[0][1] - a2 * xyz[1][1] +
                                        a2 * xyz[2][1] + (c1 - a2) * xyz[3][1] -
                                        c2 * xyz[4][1] + c2 * xyz[5][1]);
                    dxdz(2, 1) = 0.5 * ((a2 - c1) * xyz[0][2] - a2 * xyz[1][2] +
                                        a2 * xyz[2][2] + (c1 - a2) * xyz[3][2] -
                                        c2 * xyz[4][2] + c2 * xyz[5][2]);

                    dxdz(0, 2) = 0.5 * (-b1 * xyz[0][0] - b2 * xyz[3][0] +
                                        b1 * xyz[4][0] + b2 * xyz[5][0]);
                    dxdz(1, 2) = 0.5 * (-b1 * xyz[0][1] - b2 * xyz[3][1] +
                                        b1 * xyz[4][1] + b2 * xyz[5][1]);
                    dxdz(2, 2) = 0.5 * (-b1 * xyz[0][2] - b2 * xyz[3][2] +
                                        b1 * xyz[4][2] + b2 * xyz[5][2]);

                    dxdz.Invert();
                    ret.push_back(dxdz);
                }
            }
        }
    }
    else if (geom->GetShapeType() == LibUtilities::eHexahedron)
    {
        vector<Array<OneD, NekDouble> > xyz;
        for (int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> eta1                       = b[0]->GetZ();
        Array<OneD, NekDouble> eta2                       = b[1]->GetZ();
        Array<OneD, NekDouble> eta3                       = b[2]->GetZ();

        for (int k = 0; k < b[2]->GetNumPoints(); k++)
        {
            for (int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for (int i = 0; i < b[0]->GetNumPoints(); i++)
                {
                    NekDouble a1 = 0.5 * (1 - eta1[i]);
                    NekDouble a2 = 0.5 * (1 + eta1[i]);
                    NekDouble b1 = 0.5 * (1 - eta2[j]),
                              b2 = 0.5 * (1 + eta2[j]);
                    NekDouble c1 = 0.5 * (1 - eta3[k]),
                              c2 = 0.5 * (1 + eta3[k]);

                    DNekMat dxdz(3, 3, 1.0, eFULL);

                    dxdz(0, 0) =
                        -0.5 * b1 * c1 * xyz[0][0] + 0.5 * b1 * c1 * xyz[1][0] +
                        0.5 * b2 * c1 * xyz[2][0] - 0.5 * b2 * c1 * xyz[3][0] -
                        0.5 * b1 * c2 * xyz[5][0] + 0.5 * b1 * c2 * xyz[5][0] +
                        0.5 * b2 * c2 * xyz[6][0] - 0.5 * b2 * c2 * xyz[7][0];
                    dxdz(1, 0) =
                        -0.5 * b1 * c1 * xyz[0][1] + 0.5 * b1 * c1 * xyz[1][1] +
                        0.5 * b2 * c1 * xyz[2][1] - 0.5 * b2 * c1 * xyz[3][1] -
                        0.5 * b1 * c2 * xyz[5][1] + 0.5 * b1 * c2 * xyz[5][1] +
                        0.5 * b2 * c2 * xyz[6][1] - 0.5 * b2 * c2 * xyz[7][1];
                    dxdz(2, 0) =
                        -0.5 * b1 * c1 * xyz[0][2] + 0.5 * b1 * c1 * xyz[1][2] +
                        0.5 * b2 * c1 * xyz[2][2] - 0.5 * b2 * c1 * xyz[3][2] -
                        0.5 * b1 * c2 * xyz[5][2] + 0.5 * b1 * c2 * xyz[5][2] +
                        0.5 * b2 * c2 * xyz[6][2] - 0.5 * b2 * c2 * xyz[7][2];

                    dxdz(0, 1) =
                        -0.5 * a1 * c1 * xyz[0][0] - 0.5 * a2 * c1 * xyz[1][0] +
                        0.5 * a2 * c1 * xyz[2][0] + 0.5 * a1 * c1 * xyz[3][0] -
                        0.5 * a1 * c2 * xyz[5][0] - 0.5 * a2 * c2 * xyz[5][0] +
                        0.5 * a2 * c2 * xyz[6][0] + 0.5 * a1 * c2 * xyz[7][0];
                    dxdz(1, 1) =
                        -0.5 * a1 * c1 * xyz[0][1] - 0.5 * a2 * c1 * xyz[1][1] +
                        0.5 * a2 * c1 * xyz[2][1] + 0.5 * a1 * c1 * xyz[3][1] -
                        0.5 * a1 * c2 * xyz[5][1] - 0.5 * a2 * c2 * xyz[5][1] +
                        0.5 * a2 * c2 * xyz[6][1] + 0.5 * a1 * c2 * xyz[7][1];
                    dxdz(2, 1) =
                        -0.5 * a1 * c1 * xyz[0][2] - 0.5 * a2 * c1 * xyz[1][2] +
                        0.5 * a2 * c1 * xyz[2][2] + 0.5 * a1 * c1 * xyz[3][2] -
                        0.5 * a1 * c2 * xyz[5][2] - 0.5 * a2 * c2 * xyz[5][2] +
                        0.5 * a2 * c2 * xyz[6][2] + 0.5 * a1 * c2 * xyz[7][2];

                    dxdz(0, 0) =
                        -0.5 * b1 * a1 * xyz[0][0] - 0.5 * b1 * a2 * xyz[1][0] -
                        0.5 * b2 * a2 * xyz[2][0] - 0.5 * b2 * a1 * xyz[3][0] +
                        0.5 * b1 * a1 * xyz[5][0] + 0.5 * b1 * a2 * xyz[5][0] +
                        0.5 * b2 * a2 * xyz[6][0] + 0.5 * b2 * a1 * xyz[7][0];
                    dxdz(1, 0) =
                        -0.5 * b1 * a1 * xyz[0][1] - 0.5 * b1 * a2 * xyz[1][1] -
                        0.5 * b2 * a2 * xyz[2][1] - 0.5 * b2 * a1 * xyz[3][1] +
                        0.5 * b1 * a1 * xyz[5][1] + 0.5 * b1 * a2 * xyz[5][1] +
                        0.5 * b2 * a2 * xyz[6][1] + 0.5 * b2 * a1 * xyz[7][1];
                    dxdz(2, 0) =
                        -0.5 * b1 * a1 * xyz[0][2] - 0.5 * b1 * a2 * xyz[1][2] -
                        0.5 * b2 * a2 * xyz[2][2] - 0.5 * b2 * a1 * xyz[3][2] +
                        0.5 * b1 * a1 * xyz[5][2] + 0.5 * b1 * a2 * xyz[5][2] +
                        0.5 * b2 * a2 * xyz[6][2] + 0.5 * b2 * a1 * xyz[7][2];

                    dxdz.Invert();
                    ret.push_back(dxdz);
                }
            }
        }
    }
    else
    {
        ASSERTL0(false, "not coded");
    }

    return ret;
}

Array<OneD, NekDouble> ProcessQualityMetric::GetQ(
    LocalRegions::ExpansionSharedPtr e,
    bool                             s)
{
    SpatialDomains::GeometrySharedPtr geom    = e->GetGeom();
    StdRegions::StdExpansionSharedPtr chi     = e->GetGeom()->GetXmap();
    LibUtilities::PointsKeyVector p           = chi->GetPointsKeys();
    LibUtilities::PointsKeyVector pElem       = e->GetPointsKeys();
    SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();
    const int expDim                          = chi->GetNumBases();
    int nElemPts                              = 1;

    vector<LibUtilities::BasisKey> basisKeys;
    bool needsInterp = false;

    for (int i = 0; i < expDim; ++i)
    {
        nElemPts *= pElem[i].GetNumPoints();
        needsInterp =
            needsInterp || pElem[i].GetNumPoints() < p[i].GetNumPoints() - 1;
    }

    if (needsInterp)
    {
        stringstream err;
        err << "Interpolating from higher order geometry to lower order in "
            << "element " << geom->GetGlobalID();
        NEKERROR(ErrorUtil::ewarning, err.str());
    }

    for (int i = 0; i < expDim; ++i)
    {
        basisKeys.push_back(
            needsInterp
                ? chi->GetBasis(i)->GetBasisKey()
                : LibUtilities::BasisKey(chi->GetBasisType(i),
                                         chi->GetBasisNumModes(i), pElem[i]));
    }

    StdRegions::StdExpansionSharedPtr chiMod;
    switch (chi->DetShapeType())
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
        case LibUtilities::eHexahedron:
            chiMod = MemoryManager<StdRegions::StdHexExp>::AllocateSharedPtr(
                basisKeys[0], basisKeys[1], basisKeys[2]);
            break;
        default:
            ASSERTL0(false, "nope");
    }

    SpatialDomains::DerivStorage deriv = gfac->GetDeriv(pElem);

    const int pts = deriv[0][0].size();
    const int nq  = chiMod->GetTotPoints();

    ASSERTL0(pts == nq, "what");

    vector<DNekMat> i2rm = MappingIdealToRef(geom, chiMod);
    Array<OneD, NekDouble> eta(nq);

    for (int k = 0; k < pts; ++k)
    {
        DNekMat jac(expDim, expDim, 0.0, eFULL);
        DNekMat jacIdeal(expDim, expDim, 0.0, eFULL);

        for (int i = 0; i < expDim; ++i)
        {
            for (int j = 0; j < expDim; ++j)
            {
                jac(j, i) = deriv[i][j][k];
            }
        }

        jacIdeal = jac * i2rm[k];
        NekDouble jacDet = 1.0;

        if (expDim == 2)
        {
            jacDet = jacIdeal(0, 0) * jacIdeal(1, 1) -
                     jacIdeal(0, 1) * jacIdeal(1, 0);
        }
        else if (expDim == 3)
        {
            jacDet = jacIdeal(0, 0) * (jacIdeal(1, 1) * jacIdeal(2, 2) -
                                       jacIdeal(2, 1) * jacIdeal(1, 2)) -
                     jacIdeal(0, 1) * (jacIdeal(1, 0) * jacIdeal(2, 2) -
                                       jacIdeal(2, 0) * jacIdeal(1, 2)) +
                     jacIdeal(0, 2) * (jacIdeal(1, 0) * jacIdeal(2, 1) -
                                       jacIdeal(2, 0) * jacIdeal(1, 1));
        }
        else
        {
            NEKERROR(ErrorUtil::efatal, "invalid expansion dimension.");
        }

        if(s)
        {
            eta[k] = jacDet;
        }
        else
        {
            NekDouble frob = 0.0;

            for (int i = 0; i < expDim; ++i)
            {
                for (int j = 0; j < expDim; ++j)
                {
                    frob += jacIdeal(i,j) * jacIdeal(i,j);
                }
            }

            NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet));
            eta[k] = expDim * pow(sigma, 2.0/expDim) / frob;
        }
    }

    if(s)
    {
        NekDouble mx = -1.0 * numeric_limits<double>::max();
        NekDouble mn = numeric_limits<double>::max();
        for(int k = 0; k < pts; k++)
        {
            mx = max(mx,eta[k]);
            mn = min(mn,eta[k]);
        }
        for(int k = 0; k < pts; k++)
        {
            eta[k] = mn/mx;
        }
    }

    // Project onto output stuff
    if (needsInterp && pts != 1)
    {
        Array<OneD, NekDouble> tmp(nElemPts);

        if (expDim == 2)
        {
            LibUtilities::Interp2D(p[0], p[1], eta, pElem[0], pElem[1], tmp);
        }
        else if (expDim == 3)
        {
            LibUtilities::Interp3D(p[0], p[1], p[2], eta, pElem[0], pElem[1],
                                   pElem[2], tmp);
        }
        else
        {
            ASSERTL0(false, "mesh dim makes no sense");
        }

        eta = tmp;
    }

    if (pts == 1)
    {
        Vmath::Fill(nq - 1, eta[0], &eta[1], 1);
    }

    return eta;
}
}
}
