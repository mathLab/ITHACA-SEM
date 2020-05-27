////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCurvedEdges.cpp
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
//  Description: Abstract base class for creating curved edges on boundaries.
//
////////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <NekMeshUtils/MeshElements/Element.h>

#include "ProcessCurvedEdges.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{
/**
 * @brief Default constructor.
 */
ProcessCurvedEdges::ProcessCurvedEdges(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["surf"] =
        ConfigOption(false, "-1", "Tag identifying surface to process.");
    m_config["N"] = ConfigOption(false, "7", "Number of points along edge.");
}

/**
 * @brief Destructor.
 */
ProcessCurvedEdges::~ProcessCurvedEdges()
{
}

void ProcessCurvedEdges::Process()
{
    int surfTag         = m_config["surf"].as<int>();
    int prismedge[2][3] = {{0, 5, 4}, {2, 6, 7}};
    int dim = m_mesh->m_expDim;

    for (int i = 0; i < m_mesh->m_element[dim].size(); ++i)
    {
        ElementSharedPtr el = m_mesh->m_element[dim][i];
        int nSurf           = dim == 3 ? el->GetFaceCount() : el->GetEdgeCount();

        for (int j = 0; j < nSurf; ++j)
        {
            int bl = el->GetBoundaryLink(j);
            if (bl == -1)
            {
                continue;
            }

            ElementSharedPtr bEl = m_mesh->m_element[dim - 1][bl];
            vector<int> tags     = bEl->GetTagList();

            if (find(tags.begin(), tags.end(), surfTag) == tags.end())
            {
                continue;
            }

            switch (dim)
            {
                case 2:
                {
                    EdgeSharedPtr edge = el->GetEdge(j);
                    GenerateEdgeNodes(edge);
                }
                break;

                case 3:
                {
                    ASSERTL0(j == 1 || j == 3,
                            "Curved edge needs to be on prism triangular face");
                    // Check all edge interior points.
                    for (int k = 0; k < 3; ++k)
                    {
                        EdgeSharedPtr edge =
                                el->GetEdge(prismedge[(j - 1) / 2][k]);
                        GenerateEdgeNodes(edge);
                    }
                }
                break;

                default:
                    ASSERTL0(0,"Dimension not supported");
                break;
            }
        }
    }
}
}
}
