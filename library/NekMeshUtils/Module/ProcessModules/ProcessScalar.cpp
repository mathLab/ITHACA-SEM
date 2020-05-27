////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessScalar.cpp
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
//  Description: Add scalar function curvature to a given surface.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/Interpreter/Interpreter.h>
#include <NekMeshUtils/MeshElements/Element.h>

#include "ProcessScalar.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessScalar::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "scalar"),
    ProcessScalar::create,
    "Impose a scalar function z=f(x,y) on a surface.");

ProcessScalar::ProcessScalar(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["surf"] = ConfigOption(
        false, "-1", "Tag identifying surface/composite to process.");
    m_config["nq"] =
        ConfigOption(false, "-1", "Number of quadrature points to generate.");
    m_config["scalar"] = ConfigOption(false, "", "Expression to evaluate.");
}

ProcessScalar::~ProcessScalar()
{
}

void ProcessScalar::Process()
{
    int i, j, k;
    string surf = m_config["surf"].as<string>();

    // Obtain vector of surface IDs from string.
    vector<unsigned int> surfs;
    ParseUtils::GenerateSeqVector(surf, surfs);
    sort(surfs.begin(), surfs.end());

    // If we're running in verbose mode print out a list of surfaces.
    if (m_mesh->m_verbose)
    {
        cout << "ProcessScalar: extracting surface"
             << (surfs.size() > 1 ? "s" : "") << " " << surf << endl;
    }

    const int nq = m_config["nq"].as<int>();
    string expr  = m_config["scalar"].as<string>();

    LibUtilities::Interpreter rEval;
    int rExprId = rEval.DefineFunction("x y z", expr);

    // Make a copy of all existing elements of one dimension lower.
    vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim - 1];

    // Iterate over list of surface elements.
    for (i = 0; i < el.size(); ++i)
    {
        // Work out whether this lies on our surface of interest.
        vector<int> inter, tags = el[i]->GetTagList();

        sort(tags.begin(), tags.end());
        set_intersection(surfs.begin(),
                         surfs.end(),
                         tags.begin(),
                         tags.end(),
                         back_inserter(inter));

        // It doesn't continue to next element.
        if (inter.size() != 1)
        {
            continue;
        }

        // Grab face link.
        FaceSharedPtr f = el[i]->GetFaceLink();

        // Update vertices
        for (j = 0; j < 4; ++j)
        {
            NodeSharedPtr n = f->m_vertexList[j];
            n->m_z          = rEval.Evaluate(rExprId, n->m_x, n->m_y, 0.0, 0.0);

            if (n->m_z < 1e-32)
            {
                n->m_z = 0;
            }
        }

        // Put curvature into edges
        for (j = 0; j < f->m_edgeList.size(); ++j)
        {
            NodeSharedPtr n1 = f->m_edgeList[j]->m_n1;
            NodeSharedPtr n2 = f->m_edgeList[j]->m_n2;
            Node disp        = *n2 - *n1;

            f->m_edgeList[j]->m_edgeNodes.clear();

            for (k = 1; k < nq - 1; ++k)
            {
                Node n = *n1 + disp * k / (nq - 1.0);
                n.m_z = rEval.Evaluate(rExprId, n.m_x, n.m_y, 0.0, 0.0);
                if (n.m_z < 1e-32)
                {
                    n.m_z = 0;
                }

                f->m_edgeList[j]->m_edgeNodes.push_back(
                    NodeSharedPtr(new Node(n)));
            }
        }
    }
}
}
}
