////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessSurfDistance.cpp
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
//  Description: Computes height of elements connected to a surface.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessSurfDistance.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessSurfDistance::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "surfdistance"),
        ProcessSurfDistance::create,
        "Computes height of element connected to a surface.");

ProcessSurfDistance::ProcessSurfDistance(FieldSharedPtr f)
    : ProcessBoundaryExtract(f)
{
}

ProcessSurfDistance::~ProcessSurfDistance()
{
}

void ProcessSurfDistance::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);
    ASSERTL0( !boost::iequals(m_config["bnd"].as<string>(), "All"),
        "ProcessSurfDistance needs bnd parameter with a single id.");

    int i, j, k, cnt;
    int surf   = m_config["bnd"].as<int>();
    int expdim = m_f->m_graph->GetMeshDimension();

    ASSERTL0(surf >= 0, "Invalid surface " + boost::lexical_cast<string>(surf));

    int nfields           = m_f->m_variables.size();
    m_f->m_variables.push_back("dist");

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

    // Grab boundary expansions.
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp =
        exp->GetBndCondExpansions();

    // Get map that takes us from boundary element to element.
    Array<OneD, int> BoundarytoElmtID, BoundarytoTraceID;
    exp->GetBoundaryToElmtMap(BoundarytoElmtID, BoundarytoTraceID);

    ASSERTL0(!(m_f->m_numHomogeneousDir),
            "Homogeneous expansions not supported");

    for (i = cnt = 0; i < BndExp.size(); ++i)
    {
        if (i != surf)
        {
            cnt += BndExp[i]->GetExpSize();
            continue;
        }

        for (j = 0; j < BndExp[i]->GetExpSize(); ++j, ++cnt)
        {
            int elmtNum     = BoundarytoElmtID[cnt];
            int facetNum    = BoundarytoTraceID[cnt];
            int oppositeNum = 0;

            // Get boundary and element expansions.
            LocalRegions::ExpansionSharedPtr bndElmt = BndExp[i]->GetExp(j);
            LocalRegions::ExpansionSharedPtr elmt =
                exp->GetExp(elmtNum);

            // Determine which face is opposite to the surface
            switch (elmt->DetShapeType())
            {
                case LibUtilities::eQuadrilateral:
                {
                    oppositeNum = (facetNum + 2) % 4;
                }
                break;

                case LibUtilities::ePrism:
                {
                    switch (facetNum)
                    {
                        case 1:
                            oppositeNum = 3;
                            break;
                        case 3:
                            oppositeNum = 1;
                            break;
                        default:
                            ASSERTL0(false, "Surface must be on a triangular "
                                            "face of the prism.");
                    }
                }
                break;

                case LibUtilities::eHexahedron:
                {
                    switch (facetNum)
                    {
                        case 0:
                            oppositeNum = 5;
                            break;
                        case 1:
                            oppositeNum = 3;
                            break;
                        case 2:
                            oppositeNum = 4;
                            break;
                        case 3:
                            oppositeNum = 1;
                            break;
                        case 4:
                            oppositeNum = 2;
                            break;
                        case 5:
                            oppositeNum = 0;
                            break;
                        default:
                            ASSERTL0(false, "Face out of bound.");
                    }
                }
                break;

                default:
                    ASSERTL0(false, "Element not supported");
            }

            int nq    = elmt->GetTotPoints();
            int nqBnd = bndElmt->GetTotPoints();

            Array<OneD, Array<OneD, NekDouble> > x(3);
            x[0] = Array<OneD, NekDouble>(nq);
            x[1] = Array<OneD, NekDouble>(nq);
            x[2] = Array<OneD, NekDouble>(nq);
            elmt->GetCoords(x[0], x[1], x[2]);

            Array<OneD, NekDouble> face1(nqBnd), face2(nqBnd);
            Array<OneD, NekDouble> dist =
                BndExp[i]->UpdatePhys() + BndExp[i]->GetPhys_Offset(j);

            // Zero existing value.
            Vmath::Zero(nqBnd, dist, 1);

            // Calculate distance between two faces of the element
            for (k = 0; k < expdim; ++k)
            {
                switch (expdim)
                {
                    case 2:
                    {
                        elmt->GetEdgePhysVals(facetNum, bndElmt, x[k], face1);
                        elmt->GetEdgePhysVals(oppositeNum, bndElmt, x[k],
                                              face2);
                        // Consider edge orientation
                        if (elmt->GetEorient(facetNum) ==
                            elmt->GetEorient(oppositeNum))
                        {
                            Vmath::Reverse(nqBnd, face2, 1, face2, 1);
                        }
                    }
                    break;
                    case 3:
                    {
                        // Use orientation from the surface for both faces
                        StdRegions::Orientation orientation =
                            elmt->GetForient(facetNum);
                        elmt->GetFacePhysVals(facetNum, bndElmt, x[k], face1,
                                              orientation);
                        elmt->GetFacePhysVals(oppositeNum, bndElmt, x[k], face2,
                                              orientation);
                    }
                    break;
                    default:
                        ASSERTL0(false, "Expansion not supported");
                }
                Vmath::Vsub(nqBnd, face1, 1, face2, 1, face1, 1);
                Vmath::Vvtvp(nqBnd, face1, 1, face1, 1, dist, 1, dist, 1);
            }
            Vmath::Vsqrt(nqBnd, dist, 1, dist, 1);
        }

        BndExp[i]->FwdTrans_IterPerExp(BndExp[i]->GetPhys(),
                                       BndExp[i]->UpdateCoeffs());
    }
}
}
}
