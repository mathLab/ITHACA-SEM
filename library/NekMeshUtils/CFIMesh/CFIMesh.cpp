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

cfi::Entity *CFIMesh::FigureOutCADParent(cfi::NodeDefinition node,
                                         Array<OneD, NekDouble> xyz)
{
    cfi::MeshableEntity *p = node.parent;

    map<string, int>::iterator f;

    if (p->type == cfi::TYPE_LINE)
    {
        vector<string> cvs;
        vector<string> sfs;
        vector<string> cvs_possible;
        f = m_nameToCurveId.find(p->getName());
        if (f == m_nameToCurveId.end())
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
                    f               = m_nameToCurveId.find(en->getName());
                    if (f != m_nameToCurveId.end())
                    {
                        // have a curve candiate, check distance to curve
                        CADCurveSharedPtr c =
                            m_mesh->m_cad->GetCurve(f->second);

                        cvs_possible.push_back(f->first);

                        NekDouble dis = c->DistanceTo(xyz);
                        if (dis < 1e-6)
                        {
                            cvs.push_back(f->first);
                        }
                    }
                }
            }
        }
        else
        {
            cvs.push_back(f->first);
        }

        ASSERTL0(cvs.size() <= 1, "should not be more than 1 possible curve");

        if (cvs.size() == 0)
        {
            // need to be careful if the curve is not immediately found
            // the node could acutally be on a surface (think blocking sym
            // plane)
            // look at the adjacent surfaces in cvs_possible
            for (int i = 0; i < cvs_possible.size(); i++)
            {
                f = m_nameToCurveId.find(cvs_possible[i]);
                vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> >
                    ss = m_mesh->m_cad->GetCurve(f->second)->GetAdjSurf();
                for (int j = 0; j < ss.size(); j++)
                {
                    NekDouble dis = ss[j].first->DistanceTo(xyz);
                    if (dis < 1e-6)
                    {
                        sfs.push_back(ss[j].first->GetName());
                    }
                }
            }
        }

        if (cvs.size() == 1)
        {
            return m_model->getEntity(cvs[0]);
        }
        else if (cvs.size() == 0 && sfs.size() == 1)
        {
            return m_model->getEntity(sfs[0]);
        }
    }
    else if (p->type == cfi::TYPE_FACE)
    {
        // going to assume that if the parent is a face then scaning the
        // faces is sufficent
        vector<string> sfs;
        f = m_nameToFaceId.find(p->getName());
        if (f == m_nameToFaceId.end())
        {
            int na = p->getAssignmentTotal();
            for (int i = 1; i <= na; i++)
            {
                string nm         = p->getAssignmentName(i);
                cfi::Assignment a = p->getAssignment(nm);
                if (a.getType() == cfi::DATA_TYPE_ENTITY)
                {
                    cfi::Entity *en = a.getEntityValue();
                    f               = m_nameToFaceId.find(en->getName());
                    if (f != m_nameToFaceId.end())
                    {
                        CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(f->second);
                        NekDouble dis      = s->DistanceTo(xyz);
                        if (dis < 1e-6)
                        {
                            sfs.push_back(f->first);
                        }
                    }
                }
            }
        }
        else
        {
            sfs.push_back(f->first);
        }

        if (sfs.size() == 0)
        {
            return NULL;
        }
        ASSERTL0(sfs.size() == 1, "not sure");
        return m_model->getEntity(sfs[0]);
    }
    else if (p->type == cfi::TYPE_POINT)
    {
        vector<string> vts;
        vector<string> vts_possible;
        set<string> sfs; // need a set this time (duplicates)
        f = m_nameToVertId.find(p->getName());
        if (f == m_nameToVertId.end())
        {
            int na = p->getAssignmentTotal();
            for (int i = 1; i <= na; i++)
            {
                string nm         = p->getAssignmentName(i);
                cfi::Assignment a = p->getAssignment(nm);
                if (a.getType() == cfi::DATA_TYPE_ENTITY)
                {
                    cfi::Entity *en = a.getEntityValue();
                    f               = m_nameToVertId.find(en->getName());
                    if (f != m_nameToVertId.end())
                    {
                        CADVertSharedPtr v = m_mesh->m_cad->GetVert(f->second);
                        NekDouble dis      = v->DistanceTo(xyz);
                        if (dis < 1e-6)
                        {
                            vts.push_back(f->first);
                        }
                        // while it doesnt exist on this vert
                        // it may be on a surface which is adjacent to this
                        // vert
                        vts_possible.push_back(f->first);
                    }
                }
            }
        }
        else
        {
            vts.push_back(f->first);
        }

        if (vts.size() == 0)
        {
            for (int i = 0; i < vts_possible.size(); i++)
            {
                f = m_nameToVertId.find(vts_possible[i]);
                vector<CADCurveSharedPtr> cs =
                    m_mesh->m_cad->GetVert(f->second)->GetAdjCurves();
                for (int j = 0; j < cs.size(); j++)
                {
                    vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> >
                        ss = cs[j]->GetAdjSurf();

                    for (int k = 0; k < ss.size(); k++)
                    {
                        NekDouble dis = ss[k].first->DistanceTo(xyz);
                        if (dis < 1e-6)
                        {
                            sfs.insert(ss[k].first->GetName());
                        }
                    }
                }
            }
        }

        if (vts.size() == 1)
        {
            return m_model->getEntity(vts[0]);
        }
        else if (vts.size() == 0 && sfs.size() == 1)
        {
            set<string>::iterator ii = sfs.begin();
            return m_model->getEntity(*ii);
        }
        else if (vts.size() == 0 && sfs.size() == 2)
        {
            // strange time when node is acutally on a curve.
            // need to find the common curve between the two surfaces
            set<string>::iterator i1 = sfs.begin();
            map<string, int>::iterator si1 = m_nameToFaceId.find(*i1);
            i1++;
            map<string, int>::iterator si2 = m_nameToFaceId.find(*i1);
            CADSurfSharedPtr s1 = m_mesh->m_cad->GetSurf(si1->second);
            CADSurfSharedPtr s2 = m_mesh->m_cad->GetSurf(si2->second);
            vector<CADSystem::EdgeLoopSharedPtr> e1 = s1->GetEdges();
            vector<CADSystem::EdgeLoopSharedPtr> e2 = s2->GetEdges();
            for (int i = 0; i < e1.size(); i++)
            {
                for (int j = 0; j < e1[i]->edges.size(); j++)
                {
                    for (int k = 0; k < e2.size(); k++)
                    {
                        for (int l = 0; l < e2[k]->edges.size(); l++)
                        {
                            if (e1[i]->edges[j]->GetId() ==
                                e2[k]->edges[l]->GetId())
                            {
                                return m_model->getEntity(
                                    e1[i]->edges[j]->GetName());
                            }
                        }
                    }
                }
            }

            ASSERTL0(false, "failed");
        }
        else if (vts.size() == 0 && sfs.size() == 0)
        {
            return NULL;
        }
        else
        {
            ASSERTL0(false, "not sure");
        }
    }
    else
    {
        ASSERTL0(p->type == cfi::TYPE_BODY, "unsure on point type");
    }
    return NULL;
}

void CFIMesh::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << endl << "Loading mesh from CFI" << endl;
    }

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    m_cad           = boost::dynamic_pointer_cast<CADSystemCFI>(m_mesh->m_cad);
    m_nameToCurveId = m_cad->GetCFICurveId();
    m_nameToFaceId  = m_cad->GetCFIFaceId();
    m_nameToVertId  = m_cad->GetCFIVertId();
    m_model         = m_cad->GetCFIModel();
    NekDouble scal  = m_cad->GetScaling();

    map<int, NodeSharedPtr> nodes;
    vector<cfi::NodeDefinition> *cfinodes = m_model->getFenodes();

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

        cfi::Entity *p = FigureOutCADParent(*it, xyz);

        if (p == NULL)
        {
            continue;
        }

        if (p->type == cfi::TYPE_LINE)
        {
            f                   = m_nameToCurveId.find(p->getName());
            CADCurveSharedPtr c = m_mesh->m_cad->GetCurve(f->second);
            n->SetCADCurve(c, c->loct(xyz));
            vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> > ss =
                c->GetAdjSurf();
            for (int i = 0; i < ss.size(); i++)
            {
                n->SetCADSurf(ss[i].first, ss[i].first->locuv(xyz));
            }
        }
        else if (p->type == cfi::TYPE_FACE)
        {
            f                  = m_nameToFaceId.find(p->getName());
            CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(f->second);
            n->SetCADSurf(s, s->locuv(xyz));
        }
        else if (p->type == cfi::TYPE_POINT)
        {
            f                            = m_nameToVertId.find(p->getName());
            CADVertSharedPtr v           = m_mesh->m_cad->GetVert(f->second);
            vector<CADCurveSharedPtr> cs = v->GetAdjCurves();
            for (int i = 0; i < cs.size(); i++)
            {
                n->SetCADCurve(cs[i], cs[i]->loct(xyz));
                vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> >
                    ss = cs[i]->GetAdjSurf();
                for (int j = 0; j < ss.size(); j++)
                {
                    n->SetCADSurf(ss[j].first, ss[j].first->locuv(xyz));
                }
            }
        }
    }

    int prefix = m_mesh->m_cad->GetNumSurf() > 100 ? 1000 : 100;

    vector<cfi::ElementDefinition>::iterator it;

    vector<cfi::ElementDefinition> *tets =
        m_model->getElements(cfi::SUBTYPE_TE4, 4);
    cout << "tets " << tets->size() << endl;
    for (it = tets->begin(); it != tets->end(); it++)
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
        m_model->getElements(cfi::SUBTYPE_PE6, 6);
    cout << "prisms " << prisms->size() << endl;
    int nm[6] = {3, 2, 5, 0, 1, 4};
    for (it = prisms->begin(); it != prisms->end(); it++)
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
        m_model->getElements(cfi::SUBTYPE_HE8, 8);
    cout << "hexes " << hexs->size() << endl;
    for (it = hexs->begin(); it != hexs->end(); it++)
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
        m_model->getElements(cfi::SUBTYPE_TR3, 3);
    cout << "tris " << tris->size() << endl;
    for (it = tris->begin(); it != tris->end(); it++)
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
        ASSERTL0(find != m_mesh->m_faceSet.end(),
                 "surface element not found in mesh");

        FaceSharedPtr mf = *find;
        if (mf->m_elLink.size() == 1)
        {
            // boundary element, we want to use it
            map<string, int>::iterator f;
            cfi::MeshableEntity *p = (*it).parent;
            vector<int> sfs;
            f = m_nameToFaceId.find(p->getName());
            if (f == m_nameToFaceId.end())
            {
                int na = p->getAssignmentTotal();
                for (int i = 1; i <= na; i++)
                {
                    string nm         = p->getAssignmentName(i);
                    cfi::Assignment a = p->getAssignment(nm);
                    if (a.getType() == cfi::DATA_TYPE_ENTITY)
                    {
                        cfi::Entity *en = a.getEntityValue();
                        f               = m_nameToFaceId.find(en->getName());
                        if (f != m_nameToFaceId.end())
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

            // this might work
            ASSERTL0(sfs.size() == 1, "weirdness");

            // check each of the nodes know their on this surface
            for (int j = 0; j < n.size(); j++)
            {
                try
                {
                    n[j]->GetCADSurfInfo(sfs[0]);
                }
                catch (runtime_error &e)
                {
                    cout << "error" << endl;
                }
            }

            E->m_parentCAD = m_mesh->m_cad->GetSurf(sfs[0]);
            tags.clear();
            tags.push_back(sfs[0]);
            E->SetTagList(tags);
            m_mesh->m_element[2].push_back(E);
        }
    }

    vector<cfi::ElementDefinition> *quads =
        m_model->getElements(cfi::SUBTYPE_QU4, 4);
    cout << "quads " << quads->size() << endl;
    for (it = quads->begin(); it != quads->end(); it++)
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
        ASSERTL0(find != m_mesh->m_faceSet.end(),
                 "surface element not found in mesh");

        FaceSharedPtr mf = *find;
        if (mf->m_elLink.size() == 1)
        {
            // quads are wierd beucase of the CAD linkage, need to reconstruct
            // information for these (not ideal)
            vector<int> sfs;
            map<int, vector<int> > sufs;

            for (int j = 0; j < n.size(); j++)
            {
                vector<CADSurfSharedPtr> ss = n[j]->GetCADSurfs();
                for (int k = 0; k < ss.size(); k++)
                {
                    sufs[ss[k]->GetId()].push_back(0);
                }
            }

            map<int, vector<int> >::iterator ii;
            for (ii = sufs.begin(); ii != sufs.end(); ii++)
            {
                if (ii->second.size() == 4)
                {
                    sfs.push_back(ii->first);
                }
            }

            // this might work
            ASSERTL0(sfs.size() == 1, "weirdness");

            E->m_parentCAD = m_mesh->m_cad->GetSurf(sfs[0]);
            tags.clear();
            tags.push_back(sfs[0] + prefix * 2);
            E->SetTagList(tags);
            m_mesh->m_element[2].push_back(E);
        }
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    // surface edges are different to mesh edges
    // first of all need to process them so they are unique
    // and then find CAD for them from beams

    EdgeSet surfaceEdges;
    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        vector<EdgeSharedPtr> es = m_mesh->m_element[2][i]->GetEdgeList();
        for (int j = 0; j < es.size(); j++)
        {
            surfaceEdges.insert(es[j]);
        }
    }

    vector<cfi::ElementDefinition> *beams =
        m_model->getElements(cfi::SUBTYPE_BE2, 2);
    cout << "beams " << beams->size() << endl;
    int i = 0;
    for (it = beams->begin(); it != beams->end(); it++, i++)
    {
        vector<NodeSharedPtr> n;
        vector<cfi::Node *> ns = (*it).nodes;

        for (int j = 0; j < ns.size(); j++)
        {
            n.push_back(nodes[ns[j]->number]);
        }

        // going to create a edge from the cfi element and find its counterpart
        // in the m_mesh edgeset

        EdgeSharedPtr ec = EdgeSharedPtr(new Edge(n[0], n[1]));

        EdgeSet::iterator find = surfaceEdges.find(ec);
        if (find == surfaceEdges.end())
        {
            continue;
        }

        EdgeSharedPtr me = *find;
        vector<int> cvs;
        map<int, vector<int> > sufs;

        for (int j = 0; j < n.size(); j++)
        {
            vector<CADCurveSharedPtr> ss = n[j]->GetCADCurves();
            for (int k = 0; k < ss.size(); k++)
            {
                sufs[ss[k]->GetId()].push_back(0);
            }
        }

        map<int, vector<int> >::iterator ii;
        for (ii = sufs.begin(); ii != sufs.end(); ii++)
        {
            if (ii->second.size() == 2)
            {
                cvs.push_back(ii->first);
            }
        }

        if (cvs.size() > 0)
        {
            // some could be interior and therfore have no link
            ASSERTL0(cvs.size() == 1, "weirdness");
            me->m_parentCAD = m_mesh->m_cad->GetCurve(cvs[0]);
        }
    }
}
}
}
