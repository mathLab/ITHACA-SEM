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

    CADSystemCFISharedPtr cad =
        boost::dynamic_pointer_cast<CADSystemCFI>(m_mesh->m_cad);
    map<string, int> nameToCurveId = cad->GetCFICurveId();
    map<string, int> nameToFaceId  = cad->GetCFIFaceId();
    map<string, int> nameToVertId  = cad->GetCFIVertId();
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

        map<string, int>::iterator f;

        cfi::MeshableEntity *p = (*it).parent;
        if (p->type == cfi::TYPE_LINE)
        {
            vector<int> cvs;
            f = nameToCurveId.find(p->getName());
            if (f == nameToCurveId.end())
            {
                // if the imediate parent is not in the main CAD body
                // it is medial object, scan for better CAD entity
                int na = p->getAssignmentTotal();
                for (int i = 1; i <= na; i++)
                {
                    string nm         = p->getAssignmentName(i);
                    cfi::Assignment a = p->getAssignment(nm);
                    if (a.getType() == cfi::DATA_TYPE_ENTITY)
                    {
                        cfi::Entity *en = a.getEntityValue();
                        f               = nameToCurveId.find(en->getName());
                        if (f != nameToCurveId.end())
                        {
                            cvs.push_back(f->second);
                        }
                    }
                }
            }
            else
            {
                cvs.push_back(f->second);
            }

            for (int i = 0; i < cvs.size(); i++)
            {
                // found a number of possible cadidates
                // in theory should only be 1 but just incase will loop
                CADCurveSharedPtr c = m_mesh->m_cad->GetCurve(cvs[i]);

                NekDouble dis = c->DistanceTo(xyz);
                if (dis > 1e-6)
                {
                    // node is acutally not on curve but is only
                    // notionally connected to it through medial object
                    continue;
                }

                NekDouble t = c->loct(xyz);
                n->SetCADCurve(c, t);
                vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> >
                    ss = c->GetAdjSurf();
                for (int j = 0; j < ss.size(); j++)
                {
                    Array<OneD, NekDouble> uv = ss[j].first->locuv(xyz);
                    n->SetCADSurf(ss[j].first, uv);
                }
            }
        }
        else if (p->type == cfi::TYPE_FACE)
        {
            vector<int> sfs;
            f = nameToFaceId.find(p->getName());
            if (f == nameToFaceId.end())
            {
                int na = p->getAssignmentTotal();
                for (int i = 1; i <= na; i++)
                {
                    string nm         = p->getAssignmentName(i);
                    cfi::Assignment a = p->getAssignment(nm);
                    if (a.getType() == cfi::DATA_TYPE_ENTITY)
                    {
                        cfi::Entity *en = a.getEntityValue();
                        f               = nameToFaceId.find(en->getName());
                        if (f != nameToFaceId.end())
                        {
                            sfs.push_back(f->second);
                        }
                    }
                }
            }
            else
            {
                sfs.push_back(f->second);
            }

            for (int i = 0; i < sfs.size(); i++)
            {
                CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(sfs[i]);
                NekDouble dis      = s->DistanceTo(xyz);
                if (dis > 1e-6)
                {
                    continue;
                }
                Array<OneD, NekDouble> uv = s->locuv(xyz);
                n->SetCADSurf(s, uv);
            }
        }
        else if (p->type == cfi::TYPE_POINT)
        {
            vector<int> vts;
            f = nameToVertId.find(p->getName());
            if (f == nameToVertId.end())
            {
                int na = p->getAssignmentTotal();
                for (int i = 1; i <= na; i++)
                {
                    string nm         = p->getAssignmentName(i);
                    cfi::Assignment a = p->getAssignment(nm);
                    if (a.getType() == cfi::DATA_TYPE_ENTITY)
                    {
                        cfi::Entity *en = a.getEntityValue();
                        f               = nameToVertId.find(en->getName());
                        if (f != nameToVertId.end())
                        {
                            vts.push_back(f->second);
                        }
                    }
                }
            }
            else
            {
                vts.push_back(f->second);
            }

            for (int i = 0; i < vts.size(); i++)
            {
                CADVertSharedPtr v = m_mesh->m_cad->GetVert(vts[i]);
                NekDouble dis      = v->DistanceTo(xyz);
                if (dis > 1e-6)
                {
                    continue;
                }
                vector<CADCurveSharedPtr> cs = v->GetAdjCurves();
                for (int j = 0; j < cs.size(); j++)
                {
                    NekDouble t = cs[j]->loct(xyz);
                    n->SetCADCurve(cs[j], t);
                    vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> >
                        ss = cs[j]->GetAdjSurf();
                    for (int k = 0; k < ss.size(); k++)
                    {
                        Array<OneD, NekDouble> uv = ss[k].first->locuv(xyz);
                        n->SetCADSurf(ss[k].first, uv);
                    }
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

    vector<cfi::ElementDefinition> *prisms =
        model->getElements(cfi::SUBTYPE_PE6, 6);
    cout << "prisms " << prisms->size() << endl;
    i         = 0;
    int nm[6] = {3, 2, 5, 0, 1, 4};
    for (it = prisms->begin(); it != prisms->end(); it++, i++)
    {
        vector<NodeSharedPtr> n(6);
        vector<cfi::Node *> ns = (*it).nodes;

        for (int j = 0; j < ns.size(); j++)
        {
            n[nm[j]] = nodes[ns[j]->number];
        }

        vector<int> tags;
        tags.push_back(prefix + 1);
        ElmtConfig conf(LibUtilities::ePrism, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::ePrism, conf, n, tags);

        m_mesh->m_element[3].push_back(E);
    }

    vector<cfi::ElementDefinition> *hexs =
        model->getElements(cfi::SUBTYPE_HE8, 8);
    cout << "hexes " << hexs->size() << endl;
    i = 0;
    for (it = hexs->begin(); it != hexs->end(); it++, i++)
    {
        vector<NodeSharedPtr> n;
        vector<cfi::Node *> ns = (*it).nodes;

        for (int j = 0; j < ns.size(); j++)
        {
            n.push_back(nodes[ns[j]->number]);
        }

        vector<int> tags;
        tags.push_back(prefix + 2);
        ElmtConfig conf(LibUtilities::eHexahedron, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eHexahedron, conf, n, tags);

        m_mesh->m_element[3].push_back(E);
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    // ProcessComposites();

    vector<cfi::ElementDefinition> *tris =
        model->getElements(cfi::SUBTYPE_TR3, 3);
    cout << "tris " << tris->size() << endl;
    i = 0;
    for (it = tris->begin(); it != tris->end(); it++, i++)
    {
        vector<NodeSharedPtr> n;
        vector<cfi::Node *> ns = (*it).nodes;

        for (int j = 0; j < ns.size(); j++)
        {
            n.push_back(nodes[ns[j]->number]);
        }

        vector<int> tags;
        tags.push_back(0); // dummy for tempory
        ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTriangle, conf, n, tags);

        // going to create a face from the cfi element and find its counterpart
        // in the m_mesh faceset
        // if its a boundary element we can continue

        FaceSharedPtr fc = FaceSharedPtr(
            new Face(E->GetVertexList(), vector<NodeSharedPtr>(),
                     E->GetEdgeList(), LibUtilities::ePolyEvenlySpaced));

        FaceSet::iterator find = m_mesh->m_faceSet.find(fc);
        ASSERTL0(find != m_mesh->m_faceSet.end(), "surface element not found in mesh");

        FaceSharedPtr mf = *find;
        if(mf->m_elLink.size() == 1)
        {
            //boundary element, we want to use it
            map<string, int>::iterator f;
            cfi::MeshableEntity* p = (*it).parent;
            vector<int> sfs;
            f = nameToFaceId.find(p->getName());
            if (f == nameToFaceId.end())
            {
                int na = p->getAssignmentTotal();
                for (int i = 1; i <= na; i++)
                {
                    string nm         = p->getAssignmentName(i);
                    cfi::Assignment a = p->getAssignment(nm);
                    if (a.getType() == cfi::DATA_TYPE_ENTITY)
                    {
                        cfi::Entity *en = a.getEntityValue();
                        f               = nameToFaceId.find(en->getName());
                        if (f != nameToFaceId.end())
                        {
                            sfs.push_back(f->second);
                        }
                    }
                }
            }
            else
            {
                sfs.push_back(f->second);
            }

            //this might work
            ASSERTL0(sfs.size() == 1,"weirdness");
            E->m_parentCAD = m_mesh->m_cad->GetSurf(sfs[0]);
            tags.clear();
            tags.push_back(sfs[0]);
            E->SetTagList(tags);
            m_mesh->m_element[2].push_back(E);
        }
    }

    vector<cfi::ElementDefinition> *quads =
        model->getElements(cfi::SUBTYPE_QU4, 3);
    cout << "quads " << quads->size() << endl;
    i = 0;
    for (it = quads->begin(); it != quads->end(); it++, i++)
    {
        vector<NodeSharedPtr> n;
        vector<cfi::Node *> ns = (*it).nodes;

        for (int j = 0; j < ns.size(); j++)
        {
            n.push_back(nodes[ns[j]->number]);
        }

        vector<int> tags;
        tags.push_back(0); // dummy for tempory
        ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eQuadrilateral, conf, n, tags);

        // going to create a face from the cfi element and find its counterpart
        // in the m_mesh faceset
        // if its a boundary element we can continue

        FaceSharedPtr fc = FaceSharedPtr(
            new Face(E->GetVertexList(), vector<NodeSharedPtr>(),
                     E->GetEdgeList(), LibUtilities::ePolyEvenlySpaced));

        FaceSet::iterator find = m_mesh->m_faceSet.find(fc);
        ASSERTL0(find != m_mesh->m_faceSet.end(), "surface element not found in mesh");

        FaceSharedPtr mf = *find;
        if(mf->m_elLink.size() == 1)
        {
            //boundary element, we want to use it
            map<string, int>::iterator f;
            cfi::MeshableEntity* p = (*it).parent;
            vector<int> sfs;
            f = nameToFaceId.find(p->getName());
            if (f == nameToFaceId.end())
            {
                int na = p->getAssignmentTotal();
                for (int i = 1; i <= na; i++)
                {
                    string nm         = p->getAssignmentName(i);
                    cfi::Assignment a = p->getAssignment(nm);
                    if (a.getType() == cfi::DATA_TYPE_ENTITY)
                    {
                        cfi::Entity *en = a.getEntityValue();
                        f               = nameToFaceId.find(en->getName());
                        if (f != nameToFaceId.end())
                        {
                            sfs.push_back(f->second);
                        }
                    }
                }
            }
            else
            {
                sfs.push_back(f->second);
            }

            //this might work
            ASSERTL0(sfs.size() == 1,"weirdness");
            E->m_parentCAD = m_mesh->m_cad->GetSurf(sfs[0]);
            tags.clear();
            tags.push_back(sfs[0]);
            E->SetTagList(tags);
            m_mesh->m_element[2].push_back(E);
        }
    }

}
}
}
