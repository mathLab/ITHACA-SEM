////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLinear.cpp
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
//  Description: linearises mesh.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessLinear.h"

using namespace std;

namespace Nektar
{
namespace Utilities
{
ModuleKey ProcessLinear::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "linearise"),
    ProcessLinear::create,
    "Linearises mesh.");

ProcessLinear::ProcessLinear(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["all"] =
        ConfigOption(true, "0", "remove curve nodes for all elements.");
    m_config["invalid"] =
        ConfigOption(true, "0", "remove curve nodes if element is invalid.");
}

ProcessLinear::~ProcessLinear()
{
}

void ProcessLinear::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessLinear: Linearising mesh... " << endl;
    }

    bool all     = m_config["all"].as<bool>();
    bool invalid = m_config["invalid"].as<bool>();

    ASSERTL0(all || invalid, "must specify option all or invalid");

    if (all)
    {
        EdgeSet::iterator eit;
        for (eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end();
             eit++)
        {
            (*eit)->m_edgeNodes.clear();
        }

        FaceSet::iterator fit;
        for (fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end();
             fit++)
        {
            (*fit)->m_faceNodes.clear();
        }

        for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            vector<NodeSharedPtr> empty;
            m_mesh->m_element[m_mesh->m_expDim][i]->SetVolumeNodes(empty);
        }
    }
    else if (invalid)
    {
        vector<NodeSharedPtr> zeroNodes;

        map<int,vector<FaceSharedPtr> > eidToFace;
        map<int,vector<ElementSharedPtr> > eidToElm;

        vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];

        for(int i = 0; i < el.size(); i++)
        {
            vector<EdgeSharedPtr> e = el[i]->GetEdgeList();
            for(int j = 0; j < e.size(); j++)
            {
                eidToElm[e[j]->m_id].push_back(el[i]);
            }
        }

        if(m_mesh->m_expDim > 2)
        {
            FaceSet::iterator it;
            for(it = m_mesh->m_faceSet.begin();
                it != m_mesh->m_faceSet.end(); it++)
            {
                vector<EdgeSharedPtr> es = (*it)->m_edgeList;
                for(int i = 0; i < es.size(); i++)
                {
                    eidToFace[es[i]->m_id].push_back((*it));
                }
            }
        }

        set<int> neigh;

        // Iterate over list of elements of expansion dimension.
        while(el.size() > 0)
        {
            for (int i = 0; i < el.size(); ++i)
            {
                // Create elemental geometry.
                SpatialDomains::GeometrySharedPtr geom =
                    el[i]->GetGeom(m_mesh->m_spaceDim);

                // Generate geometric factors.
                SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

                if (!gfac->IsValid())
                {
                    el[i]->SetVolumeNodes(zeroNodes);

                    vector<FaceSharedPtr> f = el[i]->GetFaceList();
                    for (int j = 0; j < f.size(); j++)
                    {
                        f[j]->m_faceNodes = zeroNodes;
                    }
                    vector<EdgeSharedPtr> e = el[i]->GetEdgeList();
                    for(int j = 0; j < e.size(); j++)
                    {
                        e[j]->m_edgeNodes = zeroNodes;
                    }
                    for(int j = 0; j < e.size(); j++)
                    {
                        map<int,vector<FaceSharedPtr> >::iterator it =
                                            eidToFace.find(e[j]->m_id);
                        for(int k = 0; k < it->second.size(); k++)
                        {
                            it->second[k]->m_faceNodes = zeroNodes;
                        }
                    }
                    for(int j = 0; j < e.size(); j++)
                    {
                        map<int,vector<ElementSharedPtr> >::iterator it =
                                            eidToElm.find(e[j]->m_id);
                        for(int k = 0; k < it->second.size(); k++)
                        {
                            neigh.insert(it->second[k]->GetId());
                        }
                    }
                }
            }

            vector<ElementSharedPtr> tmp = el;
            el.clear();
            set<int>::iterator it;
            for(int i = 0; i < tmp.size(); i++)
            {
                it = neigh.find(tmp[i]->GetId());
                if(it != neigh.end())
                {
                    el.push_back(tmp[i]);
                }
            }
            neigh.clear();
        }
    }
}
}
}
