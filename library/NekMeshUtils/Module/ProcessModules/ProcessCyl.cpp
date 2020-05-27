////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCyl.cpp
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
//  Description: create cylinder curved edges
//
////////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <NekMeshUtils/MeshElements/Element.h>

#include "ProcessCyl.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessCyl::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "cyl"), ProcessCyl::create);

/**
 * @brief Default constructor.
 */
ProcessCyl::ProcessCyl(MeshSharedPtr m) : ProcessCurvedEdges(m)
{
    m_config["r"] = ConfigOption(false, "0.0", "Radius of cylinder.");
    m_config["xc"] = ConfigOption(false, "0.0", "Radius of cylinder.");
    m_config["yc"] = ConfigOption(false, "0.0", "Radius of cylinder.");
}

/**
 * @brief Destructor.
 */
ProcessCyl::~ProcessCyl()
{
}

void ProcessCyl::v_GenerateEdgeNodes(EdgeSharedPtr edge)
{
    NodeSharedPtr n1 = edge->m_n1;
    NodeSharedPtr n2 = edge->m_n2;

    int nq    = m_config["N"].as<int>();
    double r  = m_config["r"].as<double>();
    double xc  = m_config["xc"].as<double>();
    double yc  = m_config["yc"].as<double>();
    double t1 = atan2(n1->m_y - yc, n1->m_x - xc);
    double t2 = atan2(n2->m_y - yc, n2->m_x - xc);
    double dt;
    double dz;

    if (t1 < -M_PI / 2.0 && t2 > 0.0)
    {
        t1 += 2 * M_PI;
    }
    if (t2 < -M_PI / 2.0 && t1 > 0.0)
    {
        t2 += 2 * M_PI;
    }

    dt = (t2 - t1) / (nq - 1);
    dz = (n2->m_z - n1->m_z) / (nq - 1);

    edge->m_edgeNodes.resize(nq - 2);
    for (int i = 1; i < nq - 1; ++i)
    {
        edge->m_edgeNodes[i - 1] = NodeSharedPtr(new Node(
            0, xc + r * cos(t1 + i * dt),
               yc + r * sin(t1 + i * dt), n1->m_z + i * dz));
    }
    edge->m_curveType = LibUtilities::ePolyEvenlySpaced;
}
}
}
