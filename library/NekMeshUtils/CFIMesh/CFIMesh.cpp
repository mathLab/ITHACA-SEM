////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.cpp
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
//  Description: surfacemeshing object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include "CFIMesh.h"

#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/CADSystem/CFI/CADSystemCFI.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey CFIMesh::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "cfimesh"), CFIMesh::create,
    "Extracts mesh from cfi");

CFIMesh::CFIMesh(MeshSharedPtr m) : ProcessModule(m)
{
}

CFIMesh::~CFIMesh()
{
}

void CFIMesh::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << endl << "Loading mesh from CFI" << endl;
    }

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    /*CADSystemCFISharedPtr cad =
        boost::dynamic_pointer_cast<CADSystemCFI>(m_mesh->m_cad);
    map<string, int> nameToCurveId                  = cad->GetCFICurveId();
    map<string, int> nameToFaceId                   = cad->GetCFIFaceId();
    map<string, vector<string> > nameVertToListEdge = cad->GetVertId();
    cfi::Model *model = cad->GetCFIModel();
    NekDouble scal    = cad->GetScaling();

    map<int, NodeSharedPtr> nodes;
    vector<cfi::NodeDefinition> *cfinodes = model->getFenodes();

    cout << "Nodes " << cfinodes->size() << endl;

    // filter all mesh nodes into a indexed map and project to CAD
    for (vector<cfi::NodeDefinition>::iterator it = cfinodes->begin();
         it != cfinodes->end(); it++)
    {
        Array<OneD, NekDouble> xyz(3);
        cfi::Position ps = (*it).node->getXYZ();
        xyz[0]           = ps.x * scal;
        xyz[1]           = ps.y * scal;
        xyz[2]           = ps.z * scal;
        int id           = (*it).node->number;

        NodeSharedPtr n =
            boost::shared_ptr<Node>(new Node(id, xyz[0], xyz[1], xyz[2]));
        nodes.insert(pair<int, NodeSharedPtr>(id, n));

        // point built now add cad if needed

        cfi::MeshableEntity *p = (*it).parent;
        if (p->type == cfi::TYPE_LINE)
        {
            int cid             = nameToCurveId[p->getName()];
            CADCurveSharedPtr c = m_mesh->m_cad->GetCurve(cid);
            NekDouble t         = c->loct(xyz);
            n->SetCADCurve(cid, c, t);
            vector<pair<CADSurfSharedPtr, CADSystem::Orientation> > ss =
                c->GetAdjSurf();
            for (int j = 0; j < ss.size(); j++)
            {
                Array<OneD, NekDouble> uv = ss[j]->locuv(xyz);
                n->SetCADSurf(ss[j].first->GetId(), ss[j].first, uv);
            }
        }
        else if (p->type == cfi::TYPE_FACE)
        {
            int sid            = nameToFaceId[p->getName()];
            CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(sid);
            Array<OneD, NekDouble> uv = s->locuv(xyz);
            n->SetCADSurf(sid, s, uv);
        }
        else if (p->type == cfi::TYPE_POINT)
        {
            vector<string> cv = nameVertToListEdge[p->getName()];
            for (int i = 0; i < cv.size(); i++)
            {
                int cid             = nameToCurveId[cv[i]];
                CADCurveSharedPtr c = m_mesh->m_cad->GetCurve(cid);
                NekDouble t         = c->loct(xyz);
                n->SetCADCurve(cid, c, t);
                vector<pair<CADSurfSharedPtr, CADSystem::Orientation> > ss =
                    c->GetAdjSurf();
                for (int j = 0; j < ss.size(); j++)
                {
                    Array<OneD, NekDouble> uv = ss[j]->locuv(xyz);
                    n->SetCADSurf(ss[j].first->GetId(), ss[j].first, uv);
                }
            }
        }
        else
        {
            ASSERTL0(p->type == cfi::TYPE_BODY, "unsure on point type");
        }
    }

    int prefix = m_mesh->m_cad->GetNumSurf() > 100 ? 1000 : 100;

    vector<cfi::ElementDefinition>::iterator it;

    vector<cfi::ElementDefinition> *tets =
        model->getElements(cfi::SUBTYPE_TE4, 4);
    cout << "tets " << tets->size() << endl;
    int i = 0;

    for (it = tets->begin(); it != tets->end(); it++, i++)
    {
        vector<NodeSharedPtr> n;
        vector<cfi::Node *> ns = (*it).nodes;

        for (int j = 0; j < ns.size(); j++)
        {
            n.push_back(nodes[ns[j]->number]);
        }

        vector<int> tags;
        tags.push_back(prefix);
        ElmtConfig conf(LibUtilities::eTetrahedron, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTetrahedron, conf, n, tags);

        m_mesh->m_element[3].push_back(E);
    }

    // there must be a surface mesh as well or the high-order doesnt work

    vector<cfi::ElementDefinition> *tris =
        model->getElements(cfi::SUBTYPE_TR3, 3);
    ASSERTL0(tris->size() > 0, "no TR3 information found");
    cout << "tris " << tris->size() << endl;
    for (it = tris->begin(); it != tris->end(); it++)
    {
        cfi::ElementDefinition el = *it;
        vector<NodeSharedPtr> n;
        vector<cfi::Node *> ns = el.nodes;
        for (int j = 0; j < ns.size(); j++)
        {
            n.push_back(nodes[ns[j]->number]);
        }

        vector<int> tags;

        cfi::MeshableEntity *p = el.parent;
        tags.push_back(nameToFaceId[p->getName()]);

        ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTriangle, conf, n, tags);
        E->m_parentCAD = m_mesh->m_cad->GetSurf(nameToFaceId[p->getName()]);
        m_mesh->m_element[2].push_back(E);
    }

    // need to dumy process the surface edges
    EdgeSet localEdges;
    for (int i = 0; i < m_mesh->m_element[2].size(); ++i)
    {
        for (int j = 0; j < m_mesh->m_element[2][i]->GetEdgeCount(); ++j)
        {
            pair<EdgeSet::iterator, bool> testIns;
            EdgeSharedPtr ed = m_mesh->m_element[2][i]->GetEdge(j);

            testIns = localEdges.insert(ed);

            if (!testIns.second)
            {
                EdgeSharedPtr e2 = *(testIns.first);
                m_mesh->m_element[2][i]->SetEdge(j, e2);
            }
        }
    }

    // first look at the beams and add cad info, then if cad info is missing
    // edge belongs to same surface as the host triangle
    vector<cfi::ElementDefinition> *es =
        model->getElements(cfi::SUBTYPE_BE2, 2);
    ASSERTL0(es->size() > 0, "no BE2 information found");
    cout << "beams " << es->size() << endl;
    for (it = es->begin(); it != es->end(); it++)
    {
        cfi::ElementDefinition el = *it;
        vector<NodeSharedPtr> n;
        vector<cfi::Node *> ns = el.nodes;
        for (int j = 0; j < ns.size(); j++)
        {
            n.push_back(nodes[ns[j]->number]);
        }

        EdgeSharedPtr e = boost::shared_ptr<Edge>(new Edge(n[0], n[1]));

        EdgeSet::iterator f = localEdges.find(e);
        ASSERTL0(f != localEdges.end(), "edge not found");
        cfi::MeshableEntity *p = el.parent;
        (*f)->m_parentCAD =
            m_mesh->m_cad->GetCurve(nameToCurveId[p->getName()]);
    }

    for (int i = 0; i < m_mesh->m_element[2].size(); ++i)
    {
        for (int j = 0; j < m_mesh->m_element[2][i]->GetEdgeCount(); ++j)
        {
            EdgeSharedPtr ed = m_mesh->m_element[2][i]->GetEdge(j);

            if (!ed->m_parentCAD)
            {
                ed->m_parentCAD = m_mesh->m_element[2][i]->m_parentCAD;
            }
        }
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();*/
}
}
}
