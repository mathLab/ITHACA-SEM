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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <MeshUtils/MeshElements/MeshElements.h>

#include <SpatialDomains/MeshGraph.h>

#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

#include "ProcessOptiExtract.h"

using namespace std;
using namespace Nektar::MeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessOptiExtract::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "opti"), ProcessOptiExtract::create,
        "Pulls out blobs for linear elastic solver.");

ProcessOptiExtract::ProcessOptiExtract(MeshSharedPtr m) : ProcessModule(m)
{

}

ProcessOptiExtract::~ProcessOptiExtract()
{

}

void ProcessOptiExtract::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessOptiExtract: ... " << endl;
    }

    vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];

    m_mesh->m_element[m_mesh->m_expDim].clear();
    m_mesh->m_element[m_mesh->m_expDim-1].clear();

    vector<ElementSharedPtr> invalid;

    // get invalid elements
    for (int i = 0; i < el.size(); ++i)
    {
        // Create elemental geometry.
        SpatialDomains::GeometrySharedPtr geom =
            el[i]->GetGeom(m_mesh->m_spaceDim);

        // Generate geometric factors.
        SpatialDomains::GeomFactorsSharedPtr gfac =
            geom->GetGeomFactors();

        // Get the Jacobian and, if it is negative, print a warning
        // message.
        if (!gfac->IsValid())
        {
            invalid.push_back(el[i]);
        }
    }

    boost::unordered_set<int> inmesh;
    pair<boost::unordered_set<int>::iterator, bool> t;
    vector<ElementSharedPtr> totest;

    for(int i = 0; i < invalid.size(); i++)
    {
        t = inmesh.insert(invalid[i]->GetId());
        if(t.second)
            m_mesh->m_element[m_mesh->m_expDim].push_back(invalid[i]);

        vector<FaceSharedPtr> f = invalid[i]->GetFaceList();
        for(int j = 0; j < f.size(); j++)
        {
            for(int k = 0; k < f[j]->m_elLink.size(); k++)
            {
                if(f[j]->m_elLink[k].first->GetId() == invalid[i]->GetId())
                    continue;

                t = inmesh.insert(f[j]->m_elLink[k].first->GetId());
                if(t.second)
                {
                    m_mesh->m_element[m_mesh->m_expDim].push_back(f[j]->m_elLink[k].first);
                    totest.push_back(f[j]->m_elLink[k].first);
                }
            }
        }
    }

    for(int i = 0; i < 10; i++)
    {
        vector<ElementSharedPtr> tmp = totest;
        totest.clear();
        for(int j = 0; j < tmp.size(); j++)
        {
            vector<FaceSharedPtr> f = tmp[j]->GetFaceList();
            for(int k = 0; k < f.size(); k++)
            {
                for(int l = 0; l < f[k]->m_elLink.size(); l++)
                {
                    if(f[k]->m_elLink[l].first->GetId() == tmp[j]->GetId())
                        continue;

                    t = inmesh.insert(f[k]->m_elLink[l].first->GetId());
                    if(t.second)
                    {
                        m_mesh->m_element[m_mesh->m_expDim].push_back(f[k]->m_elLink[l].first);
                        totest.push_back(f[k]->m_elLink[l].first);
                    }
                }
            }
        }
    }

    ClearElementLinks();
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    el = m_mesh->m_element[m_mesh->m_expDim];
    for(int i = 0; i < el.size(); i++)
    {
        vector<FaceSharedPtr> f = el[i]->GetFaceList();
        for(int j = 0; j < f.size(); j++)
        {
            if(f[j]->m_elLink.size() == 1) //boundary element make new composite
            {
                ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);

                vector<int> tags;
                tags.push_back(1);
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                            LibUtilities::eTriangle, conf, f[j]->m_vertexList, tags);
                m_mesh->m_element[m_mesh->m_expDim-1].push_back(E);
            }
        }
    }

    if(m_mesh->m_verbose)
        cout << el.size() << " elements in blobs" << endl;

    ClearElementLinks();
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

}
}
