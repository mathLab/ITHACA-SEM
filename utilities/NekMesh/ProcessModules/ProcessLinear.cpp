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
        ConfigOption(false, "0", "remove curve nodes if element is invalid.");
    m_config["prismonly"] =
        ConfigOption(false, "", "only acts on prims");
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
    bool invalid = m_config["invalid"].beenSet;
    NekDouble thr = m_config["invalid"].as<NekDouble>();

    ASSERTL0(all || invalid,
             "Must specify an option: all (to remove all curvature) or invalid "
             "(to remove curvature that makes elements invalid)");

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

        if (m_mesh->m_verbose)
        {
            cerr << "Removed all element curvature" << endl;
        }
    }
    else if (invalid)
    {
        map<int,vector<FaceSharedPtr> > eidToFace;
        map<int,vector<ElementSharedPtr> > eidToElm;

        vector<ElementSharedPtr> els = m_mesh->m_element[m_mesh->m_expDim];
        vector<ElementSharedPtr> el = els;

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
        vector<NodeSharedPtr> zeroNodes;
        boost::unordered_set<int> clearedEdges, clearedFaces, clearedElmts;

        // Iterate over list of elements of expansion dimension.
        while(el.size() > 0)
        {
            for (int i = 0; i < el.size(); ++i)
            {
                if(m_config["prismonly"].beenSet)
                {
                    if(el[i]->GetConf().m_e != LibUtilities::ePrism)
                    {
                        continue;
                    }
                }

                if (Invalid(el[i],thr))//(!gfac->IsValid())
                {
                    clearedElmts.insert(el[i]->GetId());;
                    el[i]->SetVolumeNodes(zeroNodes);

                    vector<FaceSharedPtr> f = el[i]->GetFaceList();
                    for (int j = 0; j < f.size(); j++)
                    {
                        f[j]->m_faceNodes.clear();
                        clearedFaces.insert(f[j]->m_id);
                    }
                    vector<EdgeSharedPtr> e = el[i]->GetEdgeList();
                    for(int j = 0; j < e.size(); j++)
                    {
                        e[j]->m_edgeNodes.clear();
                        clearedEdges.insert(e[j]->m_id);
                    }

                    if(m_mesh->m_expDim > 2)
                    {
                        for(int j = 0; j < e.size(); j++)
                        {
                            map<int,vector<FaceSharedPtr> >::iterator it =
                                eidToFace.find(e[j]->m_id);
                            for(int k = 0; k < it->second.size(); k++)
                            {
                                clearedEdges.insert(it->second[k]->m_id);
                                it->second[k]->m_faceNodes.clear();
                            }
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

            el.clear();
            set<int>::iterator it1;
            boost::unordered_set<int>::iterator it2;
            for(int i = 0; i < els.size(); i++)
            {
                it1 = neigh.find(els[i]->GetId());
                it2 = clearedElmts.find(els[i]->GetId());
                if(it1 != neigh.end() && it2 == clearedElmts.end())
                {
                    el.push_back(els[i]);
                }
            }
            neigh.clear();
        }

        if (m_mesh->m_verbose)
        {
            cerr << "Removed curvature from " << clearedElmts.size()
                 << " elements (" << clearedEdges.size() << " edges, "
                 << clearedFaces.size() << " faces)" << endl;
        }
    }
}

bool ProcessLinear::Invalid(ElementSharedPtr el, NekDouble thr)
{
    // Create elemental geometry.
    SpatialDomains::GeometrySharedPtr geom =
        el->GetGeom(m_mesh->m_spaceDim);

    // Generate geometric factors.
    SpatialDomains::GeomFactorsSharedPtr gfac =
        geom->GetGeomFactors();

    if(!gfac->IsValid())
    {
        return true;
    }

    LibUtilities::PointsKeyVector p = geom->GetPointsKeys();
    SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);
    const int pts = deriv[0][0].num_elements();
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

    NekDouble scaledJac = Vmath::Vmin(jc.num_elements(),jc,1) /
                          Vmath::Vmax(jc.num_elements(),jc,1);

    if(scaledJac < thr)
    {
        return true;
    }

    return false;
}
}
}
