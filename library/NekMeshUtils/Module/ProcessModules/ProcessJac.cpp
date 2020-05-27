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

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessJac.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessJac::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "jac"),
    ProcessJac::create,
    "Process elements based on values of Jacobian.");

ProcessJac::ProcessJac(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["extract"] =
        ConfigOption(false, "0.0", "Extract non-valid elements from mesh.");
    m_config["list"] = ConfigOption(
        true, "0", "Print list of elements having negative Jacobian.");
}

ProcessJac::~ProcessJac()
{
}

void ProcessJac::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessJac: Calculating Jacobians... " << endl;
    }

    bool extract   = m_config["extract"].beenSet;
    bool printList = m_config["list"].as<bool>();
    NekDouble thres = m_config["extract"].as<NekDouble>();

    vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];

    if (extract)
    {
        m_mesh->m_element[m_mesh->m_expDim].clear();
    }

    if (printList)
    {
        cout << "Elements with negative Jacobian:" << endl;
    }

    int nNeg = 0;

    Array<OneD, int> bin(20, 0);

    // Iterate over list of elements of expansion dimension.
    for (int i = 0; i < el.size(); ++i)
    {
        // Create elemental geometry.
        SpatialDomains::GeometrySharedPtr geom =
            el[i]->GetGeom(m_mesh->m_spaceDim);

        // Generate geometric factors.
        SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

        LibUtilities::PointsKeyVector p = geom->GetXmap()->GetPointsKeys();
        SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);
        const int pts = deriv[0][0].size();
        Array<OneD,NekDouble> jc(pts);
        for (int k = 0; k < pts; ++k)
        {
            DNekMat jac(m_mesh->m_expDim, m_mesh->m_expDim, 0.0, eFULL);

            for (int l = 0; l < m_mesh->m_expDim; ++l)
            {
                for (int j = 0; j < m_mesh->m_expDim; ++j)
                {
                    jac(j,l) = deriv[l][j][k];
                }
            }

            if(m_mesh->m_expDim == 2)
            {
                jc[k] = jac(0,0) * jac(1,1) - jac(0,1)*jac(1,0);
            }
            else if(m_mesh->m_expDim == 3)
            {
                jc[k] =  jac(0,0) * (jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2)) -
                         jac(0,1) * (jac(1,0)*jac(2,2) - jac(2,0)*jac(1,2)) +
                         jac(0,2) * (jac(1,0)*jac(2,1) - jac(2,0)*jac(1,1));
            }
        }

        NekDouble scaledJac = Vmath::Vmin(jc.size(),jc,1) /
                              Vmath::Vmax(jc.size(),jc,1);

        bool valid = gfac->IsValid();

        if (extract && (scaledJac < thres || !valid))
        {
            m_mesh->m_element[m_mesh->m_expDim].push_back(el[i]);
        }

        // Get the Jacobian and, if it is negative, print a warning
        // message.
        if (!valid)
        {
            nNeg++;

            if (printList)
            {
                cout << "  - " << el[i]->GetId() << " ("
                     << LibUtilities::ShapeTypeMap[el[i]->GetConf().m_e] << ")"
                     << "  " << scaledJac
                     << endl;
            }
        }
    }

    if (extract)
    {
        m_mesh->m_element[m_mesh->m_expDim - 1].clear();
        ProcessVertices();
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
    }

    if (printList || m_mesh->m_verbose)
    {
        cout << "Total negative Jacobians: " << nNeg << endl;
    }
    else if (nNeg > 0)
    {
        cout << "WARNING: Detected " << nNeg << " element"
             << (nNeg == 1 ? "" : "s") << " with negative Jacobian." << endl;
    }
}
}
}
