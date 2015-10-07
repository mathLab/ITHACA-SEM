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
//  Description: Computes height of elements connected to a surface.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessSurfDistance.h"

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessSurfDistance::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "surfdistance"),
        ProcessSurfDistance::create,
        "Computes height of element connected to a surface.");

ProcessSurfDistance::ProcessSurfDistance(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["bnd"] = ConfigOption(false,"-1","Boundary region to calculate height");
    f->m_writeBndFld = true;
    f->m_declareExpansionAsContField = true;
    m_f->m_fldToBnd = false;
}

ProcessSurfDistance::~ProcessSurfDistance()
{
}

void ProcessSurfDistance::Process(po::variables_map &vm)
{
    int i, j, k, cnt;
    int surf = m_config["bnd"].as<int>();

    ASSERTL0(surf >= 0, "Invalid surface "+boost::lexical_cast<string>(surf));

    // Add this boundary region to the list that we will output.
    m_f->m_bndRegionsToWrite.push_back(surf);

    // Remove existing fields.
    m_f->m_exp.resize(1);

    // Grab boundary expansions.
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp =
        m_f->m_exp[0]->GetBndCondExpansions();

    // Get map that takes us from boundary element to element.
    Array<OneD, int> BoundarytoElmtID, BoundarytoTraceID;
    m_f->m_exp[0]->GetBoundaryToElmtMap(BoundarytoElmtID, BoundarytoTraceID);

    if (m_f->m_fielddef.size() == 0)
    {
        m_f->m_fielddef = m_f->m_exp[0]->GetFieldDefinitions();
        m_f->m_fielddef[0]->m_fields.push_back("dist");
    }
    else
    {
        // Override field variable
        m_f->m_fielddef[0]->m_fields[0] = "dist";
    }

    for (i = cnt = 0; i < BndExp.num_elements(); ++i)
    {
        if (i != surf)
        {
            cnt += BndExp[i]->GetExpSize();
            continue;
        }

        for (j = 0; j < BndExp[i]->GetExpSize(); ++j, ++cnt)
        {
            int elmtNum  = BoundarytoElmtID [cnt];
            int facetNum = BoundarytoTraceID[cnt];

            // Get boundary and element expansions.
            LocalRegions::ExpansionSharedPtr bndElmt = BndExp[i]->GetExp(j);
            LocalRegions::ExpansionSharedPtr elmt =
                m_f->m_exp[0]->GetExp(elmtNum);

            ASSERTL0(elmt->DetShapeType() == LibUtilities::ePrism,
                     "Only prisms supported for now!");

            ASSERTL0(facetNum == 1 || facetNum == 3,
                     "Surface must be on a triangular face of the prism.");

            int nq    = elmt   ->GetTotPoints();
            int nqBnd = bndElmt->GetTotPoints();

            Array<OneD, Array<OneD, NekDouble> > x(3);
            x[0] = Array<OneD, NekDouble>(nq);
            x[1] = Array<OneD, NekDouble>(nq);
            x[2] = Array<OneD, NekDouble>(nq);
            elmt->GetCoords(x[0], x[1], x[2]);

            Array<OneD, NekDouble> face1(nqBnd), face3(nqBnd);
            Array<OneD, NekDouble> dist =
                BndExp[i]->UpdatePhys() + BndExp[i]->GetPhys_Offset(j);

            // Zero existing value.
            Vmath::Zero(nqBnd, dist, 1);

            // Calculate distance between two faces of prism.
            for (k = 0; k < 3; ++k)
            {
                elmt->GetFacePhysVals(1, bndElmt, x[k], face1);
                elmt->GetFacePhysVals(3, bndElmt, x[k], face3);
                Vmath::Vsub (nqBnd, face1, 1, face3, 1, face1, 1);
                Vmath::Vvtvp(nqBnd, face1, 1, face1, 1, dist, 1, dist, 1);
            }
        }

        BndExp[i]->FwdTrans(BndExp[i]->GetPhys(), BndExp[i]->UpdateCoeffs());
    }
}

}
}
