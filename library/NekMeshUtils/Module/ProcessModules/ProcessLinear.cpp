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

using namespace Nektar::NekMeshUtils;

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
        ConfigOption(true, "0", "only acts on prims");
    m_config["extract"] =
        ConfigOption(false, "", "dump a mesh of the extracted elements");
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
        for (auto &edge : m_mesh->m_edgeSet)
        {
            edge->m_edgeNodes.clear();
        }

        for (auto &face : m_mesh->m_faceSet)
        {
            face->m_faceNodes.clear();
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
            for(auto &face : m_mesh->m_faceSet)
            {
                vector<EdgeSharedPtr> es = face->m_edgeList;
                for(int i = 0; i < es.size(); i++)
                {
                    eidToFace[es[i]->m_id].push_back(face);
                }
            }
        }

        set<int> neigh;
        vector<NodeSharedPtr> zeroNodes;
        std::unordered_set<int> clearedEdges, clearedFaces, clearedElmts;

        vector<ElementSharedPtr> dumpEls;

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
                    dumpEls.push_back(el[i]);
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
                            auto it = eidToFace.find(e[j]->m_id);
                            for(int k = 0; k < it->second.size(); k++)
                            {
                                clearedEdges.insert(it->second[k]->m_id);
                                it->second[k]->m_faceNodes.clear();
                            }
                        }
                    }

                    for(int j = 0; j < e.size(); j++)
                    {
                        auto it = eidToElm.find(e[j]->m_id);
                        for(int k = 0; k < it->second.size(); k++)
                        {
                            neigh.insert(it->second[k]->GetId());
                        }
                    }
                }
            }

            el.clear();
            for(int i = 0; i < els.size(); i++)
            {
                auto it1 = neigh.find(els[i]->GetId());
                auto it2 = clearedElmts.find(els[i]->GetId());
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

        if(m_config["extract"].beenSet)
        {
            MeshSharedPtr dmp = std::shared_ptr<Mesh>(new Mesh());
            dmp->m_expDim     = 3;
            dmp->m_spaceDim   = 3;
            dmp->m_nummode    = 2;

            dmp->m_element[3] = dumpEls;

            ModuleSharedPtr mod = GetModuleFactory().CreateInstance(
                ModuleKey(eOutputModule, "xml"), dmp);
            mod->RegisterConfig("outfile", m_config["extract"].as<string>().c_str());
            mod->ProcessVertices();
            mod->ProcessEdges();
            mod->ProcessFaces();
            mod->ProcessElements();
            mod->ProcessComposites();
            mod->Process();
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

    vector<NodeSharedPtr> ns = el->GetVertexList();
    ElmtConfig c = el->GetConf();
    c.m_order = 1;
    c.m_faceNodes = false;
    c.m_volumeNodes = false;
    c.m_reorient = false;

    ElementSharedPtr elL = GetElementFactory().CreateInstance(
                c.m_e, c, ns, el->GetTagList());
    SpatialDomains::GeometrySharedPtr geomL = elL->GetGeom(m_mesh->m_spaceDim);
    SpatialDomains::GeomFactorsSharedPtr gfacL = geomL->GetGeomFactors();

    LibUtilities::PointsKeyVector p = geom->GetXmap()->GetPointsKeys();
    SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);
    SpatialDomains::DerivStorage derivL = gfacL->GetDeriv(p);
    const int pts = deriv[0][0].size();
    Array<OneD,NekDouble> jc(pts);
    Array<OneD,NekDouble> jcL(pts);
    for (int k = 0; k < pts; ++k)
    {
        DNekMat jac(m_mesh->m_expDim, m_mesh->m_expDim, 0.0, eFULL);
        DNekMat jacL(m_mesh->m_expDim, m_mesh->m_expDim, 0.0, eFULL);

        for (int l = 0; l < m_mesh->m_expDim; ++l)
        {
            for (int j = 0; j < m_mesh->m_expDim; ++j)
            {
                jac(j,l) = deriv[l][j][k];
                jacL(j,l) = derivL[l][j][k];
            }
        }

        if(m_mesh->m_expDim == 2)
        {
            jc[k] = jac(0,0) * jac(1,1) - jac(0,1)*jac(1,0);
            jcL[k] = jacL(0,0) * jacL(1,1) - jacL(0,1)*jacL(1,0);
        }
        else if(m_mesh->m_expDim == 3)
        {
            jc[k] =  jac(0,0) * (jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2)) -
                     jac(0,1) * (jac(1,0)*jac(2,2) - jac(2,0)*jac(1,2)) +
                     jac(0,2) * (jac(1,0)*jac(2,1) - jac(2,0)*jac(1,1));
            jcL[k] =  jacL(0,0) * (jacL(1,1)*jacL(2,2) - jacL(2,1)*jacL(1,2)) -
                      jacL(0,1) * (jacL(1,0)*jacL(2,2) - jacL(2,0)*jacL(1,2)) +
                      jacL(0,2) * (jacL(1,0)*jacL(2,1) - jacL(2,0)*jacL(1,1));
        }
    }

    Array<OneD, NekDouble> j(pts);
    Vmath::Vdiv(jc.size(),jc,1,jcL,1,j,1);

    NekDouble scaledJac = Vmath::Vmin(j.size(),j,1) /
                          Vmath::Vmax(j.size(),j,1);

    //cout << scaledJac << endl;

    if(scaledJac < thr)
    {
        return true;
    }

    return false;
}
}
}
