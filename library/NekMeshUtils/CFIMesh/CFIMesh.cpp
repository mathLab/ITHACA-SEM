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

ofstream file;

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

    file.open("pts.3D");
    file << "x y z value" << endl;

    // filter all mesh nodes into a indexed map and project to CAD
    int k = 0;
    for (vector<cfi::NodeDefinition>::iterator it = cfinodes->begin();
         it != cfinodes->end(); it++, k++)
    {
        Array<OneD, NekDouble> xyz(3);
        cfi::Position ps = (*it).node->getXYZ();
        xyz[0]           = ps.x * scal;
        xyz[1]           = ps.y * scal;
        xyz[2]           = ps.z * scal;
        int id           = (*it).node->number;

        NodeSharedPtr n =
            boost::shared_ptr<Node>(new Node(k, xyz[0], xyz[1], xyz[2]));
        nodes.insert(pair<int, NodeSharedPtr>(id, n));

        // point built now add cad if needed

        cfi::MeshableEntity *p = (*it).parent;

        if (p->type == cfi::TYPE_LINE)
        {
            auto f = m_nameToCurveId.find(p->getName());
            if(f != m_nameToCurveId.end())
            {
                CADCurveSharedPtr c = m_mesh->m_cad->GetCurve(f->second);

                NekDouble t = c->loct(xyz);

                n->SetCADCurve(c, t);

                vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> >
                    ss = c->GetAdjSurf();
                for (int j = 0; j < ss.size(); j++)
                {
                    Array<OneD, NekDouble> uv(2);
                    uv = ss[j].first->locuv(xyz);
                    n->SetCADSurf(ss[j].first, uv);
                }
            }
        }
        else if (p->type == cfi::TYPE_FACE)
        {
            auto f = m_nameToFaceId.find(p->getName());
            if(f != m_nameToFaceId.end())
            {
                CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(f->second);

                Array<OneD, NekDouble> uv(2);
                uv = s->locuv(xyz);

                n->SetCADSurf(s, uv);
            }
        }
        else if (p->type == cfi::TYPE_POINT)
        {
            auto f = m_nameToVertId.find(p->getName());
            if(f != m_nameToVertId.end())
            {
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
    }

    ////
    //Really important fact. Nodes must be renumbered as they are read by the elements
    //such that vertical edges on the prism are sequential
    //In doing so we ensure the orienation will work
    //we dont want to renumber nodes that have already been numbered, hence the set
    //the set will be tacked by the cfiID as that is a constant

    map<int, set<LibUtilities::ShapeType> > cfiIdToTypes;

    vector<cfi::ElementDefinition> *prisms =
        m_model->getElements(cfi::SUBTYPE_PE6, 6);
    for(auto &it : *prisms)
    {
        vector<cfi::Node *> ns = it.nodes;
        for(int i = 0; i < ns.size(); i++)
        {
            cfiIdToTypes[ns[i]->number].insert(LibUtilities::ePrism);
        }
    }
    vector<cfi::ElementDefinition> *hexs =
        m_model->getElements(cfi::SUBTYPE_HE8, 8);
    for(auto &it : *hexs)
    {
        vector<cfi::Node *> ns = it.nodes;
        for(int i = 0; i < ns.size(); i++)
        {
            cfiIdToTypes[ns[i]->number].insert(LibUtilities::eHexahedron);
        }
    }
    vector<cfi::ElementDefinition> *tets =
        m_model->getElements(cfi::SUBTYPE_TE4, 4);
    for(auto &it : *tets)
    {
        vector<cfi::Node *> ns = it.nodes;
        for(int i = 0; i < ns.size(); i++)
        {
            cfiIdToTypes[ns[i]->number].insert(LibUtilities::eTetrahedron);
        }
    }

    ASSERTL0(nodes.size() == cfiIdToTypes.size(), "not all nodes marked");

    int id = 0;

    for(int i = 3; i > 0; i--)
    {
        for(auto &it : cfiIdToTypes)
        {
            if(it.second.size() == i)
            {
                NodeSharedPtr n  = nodes[it.first];
                if(id  <= 230)
                {
                    file << n->m_x << " " << n->m_y << " " << n->m_z << " " << id << endl;
                }

                nodes[it.first]->m_id = id++;
            }
        }
    }


    ASSERTL0(id == nodes.size(), "not all nodes numbered");


    int prefix = m_mesh->m_cad->GetNumSurf() > 100 ? 1000 : 100;

    vector<cfi::ElementDefinition>::iterator it;

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
    ProcessComposites();

    return;

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
        ElmtConfig conf(LibUtilities::eTriangle, 1, false, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTriangle, conf, n, tags);

        cfi::MeshableEntity *p = (*it).parent;

        auto f = m_nameToFaceId.find(p->getName());

        if(f != m_nameToFaceId.end())
        {
            E->m_parentCAD = m_mesh->m_cad->GetSurf(f->second);
            m_mesh->m_element[2].push_back(E);
        }

        // going to create a face from the cfi element and find its counterpart
        // in the m_mesh faceset
        // if its a boundary element we can continue

        /*FaceSharedPtr fc = FaceSharedPtr(
            new Face(E->GetVertexList(), vector<NodeSharedPtr>(),
                     E->GetEdgeList(), LibUtilities::ePolyEvenlySpaced));

        FaceSet::iterator fnd = m_mesh->m_faceSet.find(fc);
        ASSERTL0(fnd != m_mesh->m_faceSet.end(),
                 "surface element not found in mesh");

        FaceSharedPtr mf = *fnd;*/

        /*if (mf->m_elLink.size() == 1)
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
            ASSERTL0(sfs.size() == 1, "weirdness " + boost::lexical_cast<string>(sfs.size()));

            int error = 0;
            // check each of the nodes know their on this surface
            for (int j = 0; j < n.size(); j++)
            {
                vector<int> ids = n[j]->GetCADSurfsIds();
                vector<int>::iterator a = find(ids.begin(), ids.end(), sfs[0]);
                if(a == ids.end())
                {
                    error++;
                }
            }

            if(error !=0)
            {
                cout << "unable to make triangle high-order because of node on " << sfs[0] << endl;
            }
            else
            {
                E->m_parentCAD = m_mesh->m_cad->GetSurf(sfs[0]);
            }

            //if(!error)
            //{
            //    continue;
            //}

            tags.clear();
            tags.push_back(sfs[0]);
            E->SetTagList(tags);
            m_mesh->m_element[2].push_back(E);
        }*/
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
        tags.push_back(1); // dummy for tempory
        ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eQuadrilateral, conf, n, tags);

        //m_mesh->m_element[2].push_back(E);

        // going to create a face from the cfi element and find its counterpart
        // in the m_mesh faceset
        // if its a boundary element we can continue

        /*FaceSharedPtr fc = FaceSharedPtr(
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
            tags.clear();
            if(sfs.size() == 1)
            {
                E->m_parentCAD = m_mesh->m_cad->GetSurf(sfs[0]);
                tags.push_back(sfs[0] + prefix * 2);
            }
            else
            {
                cout << "cannot make quad high-order " << endl;
                tags.push_back(5000);
            }

            E->SetTagList(tags);
            m_mesh->m_element[2].push_back(E);
        }*/
    }

    //m_mesh->m_element[3].clear();
    m_mesh->m_expDim=2;

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    //return;

    // surface edges are different to mesh edges
    // first of all need to process them so they are unique
    // and then find CAD for them from beams

    /*EdgeSet surfaceEdges;
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
    }*/

}
}
}
